#include "ElectronRepulsionPrimRecSHSI.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_shsi(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_shsi,
                                  size_t                idx_eri_0_sfsi,
                                  size_t                idx_eri_1_sfsi,
                                  size_t                idx_eri_1_sgsh,
                                  size_t                idx_eri_0_sgsi,
                                  size_t                idx_eri_1_sgsi,
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

    /// Set up components of auxilary buffer : SFSI

    auto g_0_xxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi);

    auto g_0_xxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 1);

    auto g_0_xxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 2);

    auto g_0_xxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 3);

    auto g_0_xxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 4);

    auto g_0_xxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 5);

    auto g_0_xxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 6);

    auto g_0_xxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 7);

    auto g_0_xxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 8);

    auto g_0_xxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 9);

    auto g_0_xxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 10);

    auto g_0_xxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 11);

    auto g_0_xxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 12);

    auto g_0_xxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 13);

    auto g_0_xxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 14);

    auto g_0_xxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 15);

    auto g_0_xxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 16);

    auto g_0_xxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 17);

    auto g_0_xxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 18);

    auto g_0_xxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 19);

    auto g_0_xxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 20);

    auto g_0_xxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 21);

    auto g_0_xxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 22);

    auto g_0_xxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 23);

    auto g_0_xxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 24);

    auto g_0_xxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 25);

    auto g_0_xxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 26);

    auto g_0_xxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 27);

    auto g_0_xxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 28);

    auto g_0_xxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 30);

    auto g_0_xxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 33);

    auto g_0_xxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 37);

    auto g_0_xxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 42);

    auto g_0_xxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 48);

    auto g_0_xxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 56);

    auto g_0_xxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 57);

    auto g_0_xxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 59);

    auto g_0_xxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 62);

    auto g_0_xxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 66);

    auto g_0_xxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 71);

    auto g_0_xyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 85);

    auto g_0_xyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 87);

    auto g_0_xyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 88);

    auto g_0_xyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 90);

    auto g_0_xyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 91);

    auto g_0_xyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 92);

    auto g_0_xyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 94);

    auto g_0_xyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 95);

    auto g_0_xyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 96);

    auto g_0_xyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 97);

    auto g_0_xyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 99);

    auto g_0_xyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 100);

    auto g_0_xyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 101);

    auto g_0_xyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 102);

    auto g_0_xyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 103);

    auto g_0_xyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 105);

    auto g_0_xyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 106);

    auto g_0_xyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 107);

    auto g_0_xyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 108);

    auto g_0_xyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 109);

    auto g_0_xyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 110);

    auto g_0_xyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 111);

    auto g_0_xzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 142);

    auto g_0_xzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 144);

    auto g_0_xzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 145);

    auto g_0_xzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 147);

    auto g_0_xzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 148);

    auto g_0_xzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 149);

    auto g_0_xzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 151);

    auto g_0_xzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 152);

    auto g_0_xzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 153);

    auto g_0_xzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 154);

    auto g_0_xzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 156);

    auto g_0_xzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 157);

    auto g_0_xzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 158);

    auto g_0_xzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 159);

    auto g_0_xzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 160);

    auto g_0_xzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 161);

    auto g_0_xzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 162);

    auto g_0_xzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 163);

    auto g_0_xzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 164);

    auto g_0_xzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 165);

    auto g_0_xzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 166);

    auto g_0_xzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 167);

    auto g_0_yyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 168);

    auto g_0_yyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 169);

    auto g_0_yyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 170);

    auto g_0_yyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 171);

    auto g_0_yyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 172);

    auto g_0_yyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 173);

    auto g_0_yyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 174);

    auto g_0_yyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 175);

    auto g_0_yyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 176);

    auto g_0_yyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 177);

    auto g_0_yyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 178);

    auto g_0_yyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 179);

    auto g_0_yyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 180);

    auto g_0_yyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 181);

    auto g_0_yyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 182);

    auto g_0_yyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 183);

    auto g_0_yyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 184);

    auto g_0_yyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 185);

    auto g_0_yyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 186);

    auto g_0_yyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 187);

    auto g_0_yyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 188);

    auto g_0_yyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 189);

    auto g_0_yyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 190);

    auto g_0_yyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 191);

    auto g_0_yyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 192);

    auto g_0_yyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 193);

    auto g_0_yyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 194);

    auto g_0_yyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 195);

    auto g_0_yyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 197);

    auto g_0_yyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 199);

    auto g_0_yyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 202);

    auto g_0_yyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 206);

    auto g_0_yyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 211);

    auto g_0_yyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 217);

    auto g_0_yzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 224);

    auto g_0_yzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 226);

    auto g_0_yzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 228);

    auto g_0_yzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 229);

    auto g_0_yzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 231);

    auto g_0_yzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 232);

    auto g_0_yzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 233);

    auto g_0_yzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 235);

    auto g_0_yzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 236);

    auto g_0_yzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 237);

    auto g_0_yzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 238);

    auto g_0_yzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 240);

    auto g_0_yzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 241);

    auto g_0_yzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 242);

    auto g_0_yzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 243);

    auto g_0_yzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 244);

    auto g_0_yzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 246);

    auto g_0_yzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 247);

    auto g_0_yzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 248);

    auto g_0_yzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 249);

    auto g_0_yzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 250);

    auto g_0_yzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 251);

    auto g_0_zzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 252);

    auto g_0_zzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 253);

    auto g_0_zzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 254);

    auto g_0_zzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 255);

    auto g_0_zzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 256);

    auto g_0_zzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 257);

    auto g_0_zzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 258);

    auto g_0_zzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 259);

    auto g_0_zzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 260);

    auto g_0_zzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 261);

    auto g_0_zzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 262);

    auto g_0_zzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 263);

    auto g_0_zzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 264);

    auto g_0_zzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 265);

    auto g_0_zzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 266);

    auto g_0_zzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 267);

    auto g_0_zzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 268);

    auto g_0_zzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 269);

    auto g_0_zzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 270);

    auto g_0_zzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 271);

    auto g_0_zzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 272);

    auto g_0_zzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 273);

    auto g_0_zzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 274);

    auto g_0_zzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 275);

    auto g_0_zzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 276);

    auto g_0_zzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 277);

    auto g_0_zzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 278);

    auto g_0_zzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 279);

    /// Set up components of auxilary buffer : SFSI

    auto g_0_xxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sfsi);

    auto g_0_xxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sfsi + 1);

    auto g_0_xxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sfsi + 2);

    auto g_0_xxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sfsi + 3);

    auto g_0_xxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sfsi + 4);

    auto g_0_xxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sfsi + 5);

    auto g_0_xxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sfsi + 6);

    auto g_0_xxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sfsi + 7);

    auto g_0_xxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sfsi + 8);

    auto g_0_xxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sfsi + 9);

    auto g_0_xxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 10);

    auto g_0_xxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 11);

    auto g_0_xxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 12);

    auto g_0_xxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 13);

    auto g_0_xxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 14);

    auto g_0_xxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 15);

    auto g_0_xxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 16);

    auto g_0_xxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 17);

    auto g_0_xxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 18);

    auto g_0_xxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 19);

    auto g_0_xxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 20);

    auto g_0_xxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 21);

    auto g_0_xxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 22);

    auto g_0_xxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 23);

    auto g_0_xxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 24);

    auto g_0_xxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 25);

    auto g_0_xxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 26);

    auto g_0_xxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 27);

    auto g_0_xxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sfsi + 28);

    auto g_0_xxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sfsi + 30);

    auto g_0_xxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sfsi + 33);

    auto g_0_xxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sfsi + 37);

    auto g_0_xxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 42);

    auto g_0_xxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 48);

    auto g_0_xxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sfsi + 56);

    auto g_0_xxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sfsi + 57);

    auto g_0_xxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sfsi + 59);

    auto g_0_xxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sfsi + 62);

    auto g_0_xxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 66);

    auto g_0_xxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 71);

    auto g_0_xyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sfsi + 85);

    auto g_0_xyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sfsi + 87);

    auto g_0_xyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sfsi + 88);

    auto g_0_xyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sfsi + 90);

    auto g_0_xyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sfsi + 91);

    auto g_0_xyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sfsi + 92);

    auto g_0_xyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 94);

    auto g_0_xyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 95);

    auto g_0_xyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 96);

    auto g_0_xyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 97);

    auto g_0_xyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 99);

    auto g_0_xyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 100);

    auto g_0_xyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 101);

    auto g_0_xyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 102);

    auto g_0_xyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 103);

    auto g_0_xyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 105);

    auto g_0_xyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 106);

    auto g_0_xyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 107);

    auto g_0_xyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 108);

    auto g_0_xyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 109);

    auto g_0_xyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 110);

    auto g_0_xyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 111);

    auto g_0_xzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sfsi + 142);

    auto g_0_xzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sfsi + 144);

    auto g_0_xzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sfsi + 145);

    auto g_0_xzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sfsi + 147);

    auto g_0_xzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sfsi + 148);

    auto g_0_xzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sfsi + 149);

    auto g_0_xzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 151);

    auto g_0_xzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 152);

    auto g_0_xzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 153);

    auto g_0_xzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 154);

    auto g_0_xzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 156);

    auto g_0_xzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 157);

    auto g_0_xzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 158);

    auto g_0_xzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 159);

    auto g_0_xzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 160);

    auto g_0_xzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 161);

    auto g_0_xzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 162);

    auto g_0_xzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 163);

    auto g_0_xzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 164);

    auto g_0_xzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 165);

    auto g_0_xzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 166);

    auto g_0_xzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 167);

    auto g_0_yyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sfsi + 168);

    auto g_0_yyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sfsi + 169);

    auto g_0_yyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sfsi + 170);

    auto g_0_yyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sfsi + 171);

    auto g_0_yyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sfsi + 172);

    auto g_0_yyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sfsi + 173);

    auto g_0_yyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sfsi + 174);

    auto g_0_yyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sfsi + 175);

    auto g_0_yyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sfsi + 176);

    auto g_0_yyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sfsi + 177);

    auto g_0_yyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 178);

    auto g_0_yyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 179);

    auto g_0_yyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 180);

    auto g_0_yyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 181);

    auto g_0_yyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 182);

    auto g_0_yyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 183);

    auto g_0_yyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 184);

    auto g_0_yyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 185);

    auto g_0_yyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 186);

    auto g_0_yyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 187);

    auto g_0_yyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 188);

    auto g_0_yyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 189);

    auto g_0_yyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 190);

    auto g_0_yyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 191);

    auto g_0_yyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 192);

    auto g_0_yyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 193);

    auto g_0_yyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 194);

    auto g_0_yyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 195);

    auto g_0_yyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sfsi + 197);

    auto g_0_yyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sfsi + 199);

    auto g_0_yyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sfsi + 202);

    auto g_0_yyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 206);

    auto g_0_yyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 211);

    auto g_0_yyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 217);

    auto g_0_yzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sfsi + 224);

    auto g_0_yzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sfsi + 226);

    auto g_0_yzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sfsi + 228);

    auto g_0_yzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sfsi + 229);

    auto g_0_yzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sfsi + 231);

    auto g_0_yzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sfsi + 232);

    auto g_0_yzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sfsi + 233);

    auto g_0_yzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 235);

    auto g_0_yzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 236);

    auto g_0_yzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 237);

    auto g_0_yzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 238);

    auto g_0_yzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 240);

    auto g_0_yzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 241);

    auto g_0_yzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 242);

    auto g_0_yzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 243);

    auto g_0_yzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 244);

    auto g_0_yzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 246);

    auto g_0_yzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 247);

    auto g_0_yzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 248);

    auto g_0_yzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 249);

    auto g_0_yzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 250);

    auto g_0_yzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 251);

    auto g_0_zzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sfsi + 252);

    auto g_0_zzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sfsi + 253);

    auto g_0_zzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sfsi + 254);

    auto g_0_zzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sfsi + 255);

    auto g_0_zzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sfsi + 256);

    auto g_0_zzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sfsi + 257);

    auto g_0_zzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sfsi + 258);

    auto g_0_zzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sfsi + 259);

    auto g_0_zzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sfsi + 260);

    auto g_0_zzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sfsi + 261);

    auto g_0_zzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 262);

    auto g_0_zzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 263);

    auto g_0_zzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 264);

    auto g_0_zzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 265);

    auto g_0_zzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 266);

    auto g_0_zzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 267);

    auto g_0_zzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 268);

    auto g_0_zzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 269);

    auto g_0_zzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 270);

    auto g_0_zzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 271);

    auto g_0_zzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 272);

    auto g_0_zzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sfsi + 273);

    auto g_0_zzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sfsi + 274);

    auto g_0_zzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sfsi + 275);

    auto g_0_zzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sfsi + 276);

    auto g_0_zzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 277);

    auto g_0_zzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 278);

    auto g_0_zzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sfsi + 279);

    /// Set up components of auxilary buffer : SGSH

    auto g_0_xxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh);

    auto g_0_xxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 1);

    auto g_0_xxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 2);

    auto g_0_xxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 3);

    auto g_0_xxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 4);

    auto g_0_xxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 5);

    auto g_0_xxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 6);

    auto g_0_xxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 7);

    auto g_0_xxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 8);

    auto g_0_xxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 9);

    auto g_0_xxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 10);

    auto g_0_xxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 11);

    auto g_0_xxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 12);

    auto g_0_xxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 13);

    auto g_0_xxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 14);

    auto g_0_xxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 15);

    auto g_0_xxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 16);

    auto g_0_xxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 17);

    auto g_0_xxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 18);

    auto g_0_xxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 19);

    auto g_0_xxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 20);

    auto g_0_xxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 44);

    auto g_0_xxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 46);

    auto g_0_xxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 47);

    auto g_0_xxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 49);

    auto g_0_xxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 50);

    auto g_0_xxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 51);

    auto g_0_xxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 53);

    auto g_0_xxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 54);

    auto g_0_xxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 55);

    auto g_0_xxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 56);

    auto g_0_xxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 58);

    auto g_0_xxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 59);

    auto g_0_xxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 60);

    auto g_0_xxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 61);

    auto g_0_xxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 62);

    auto g_0_xxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 63);

    auto g_0_xxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 64);

    auto g_0_xxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 65);

    auto g_0_xxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 66);

    auto g_0_xxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 67);

    auto g_0_xxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 68);

    auto g_0_xxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 69);

    auto g_0_xxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 70);

    auto g_0_xxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 71);

    auto g_0_xxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 72);

    auto g_0_xxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 73);

    auto g_0_xxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 74);

    auto g_0_xxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 75);

    auto g_0_xxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 76);

    auto g_0_xxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 77);

    auto g_0_xxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 78);

    auto g_0_xxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 79);

    auto g_0_xxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 80);

    auto g_0_xxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 81);

    auto g_0_xxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 82);

    auto g_0_xxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 83);

    auto g_0_xxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 105);

    auto g_0_xxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 106);

    auto g_0_xxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 107);

    auto g_0_xxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 108);

    auto g_0_xxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 109);

    auto g_0_xxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 110);

    auto g_0_xxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 111);

    auto g_0_xxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 112);

    auto g_0_xxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 113);

    auto g_0_xxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 114);

    auto g_0_xxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 115);

    auto g_0_xxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 116);

    auto g_0_xxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 117);

    auto g_0_xxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 118);

    auto g_0_xxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 119);

    auto g_0_xxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 120);

    auto g_0_xxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 121);

    auto g_0_xxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 122);

    auto g_0_xxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 123);

    auto g_0_xxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 124);

    auto g_0_xxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 125);

    auto g_0_xyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 127);

    auto g_0_xyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 129);

    auto g_0_xyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 130);

    auto g_0_xyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 132);

    auto g_0_xyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 133);

    auto g_0_xyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 134);

    auto g_0_xyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 136);

    auto g_0_xyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 137);

    auto g_0_xyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 138);

    auto g_0_xyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 139);

    auto g_0_xyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 141);

    auto g_0_xyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 142);

    auto g_0_xyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 143);

    auto g_0_xyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 144);

    auto g_0_xyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 145);

    auto g_0_xzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 191);

    auto g_0_xzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 193);

    auto g_0_xzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 194);

    auto g_0_xzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 196);

    auto g_0_xzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 197);

    auto g_0_xzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 198);

    auto g_0_xzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 200);

    auto g_0_xzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 201);

    auto g_0_xzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 202);

    auto g_0_xzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 203);

    auto g_0_xzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 205);

    auto g_0_xzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 206);

    auto g_0_xzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 207);

    auto g_0_xzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 208);

    auto g_0_xzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 209);

    auto g_0_yyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 210);

    auto g_0_yyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 211);

    auto g_0_yyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 212);

    auto g_0_yyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 213);

    auto g_0_yyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 214);

    auto g_0_yyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 215);

    auto g_0_yyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 216);

    auto g_0_yyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 217);

    auto g_0_yyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 218);

    auto g_0_yyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 219);

    auto g_0_yyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 220);

    auto g_0_yyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 221);

    auto g_0_yyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 222);

    auto g_0_yyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 223);

    auto g_0_yyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 224);

    auto g_0_yyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 225);

    auto g_0_yyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 226);

    auto g_0_yyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 227);

    auto g_0_yyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 228);

    auto g_0_yyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 229);

    auto g_0_yyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 230);

    auto g_0_yyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 233);

    auto g_0_yyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 235);

    auto g_0_yyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 236);

    auto g_0_yyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 238);

    auto g_0_yyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 239);

    auto g_0_yyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 240);

    auto g_0_yyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 242);

    auto g_0_yyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 243);

    auto g_0_yyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 244);

    auto g_0_yyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 245);

    auto g_0_yyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 247);

    auto g_0_yyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 248);

    auto g_0_yyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 249);

    auto g_0_yyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 250);

    auto g_0_yyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 251);

    auto g_0_yyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 252);

    auto g_0_yyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 253);

    auto g_0_yyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 254);

    auto g_0_yyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 255);

    auto g_0_yyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 256);

    auto g_0_yyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 257);

    auto g_0_yyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 258);

    auto g_0_yyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 259);

    auto g_0_yyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 260);

    auto g_0_yyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 261);

    auto g_0_yyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 262);

    auto g_0_yyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 263);

    auto g_0_yyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 264);

    auto g_0_yyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 265);

    auto g_0_yyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 266);

    auto g_0_yyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 267);

    auto g_0_yyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 268);

    auto g_0_yyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 269);

    auto g_0_yyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 270);

    auto g_0_yyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 271);

    auto g_0_yyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 272);

    auto g_0_yzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 274);

    auto g_0_yzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 275);

    auto g_0_yzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 276);

    auto g_0_yzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 277);

    auto g_0_yzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 278);

    auto g_0_yzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 279);

    auto g_0_yzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 280);

    auto g_0_yzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 281);

    auto g_0_yzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 282);

    auto g_0_yzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 283);

    auto g_0_yzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 284);

    auto g_0_yzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 285);

    auto g_0_yzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 286);

    auto g_0_yzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 287);

    auto g_0_yzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 288);

    auto g_0_yzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 289);

    auto g_0_yzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 290);

    auto g_0_yzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 291);

    auto g_0_yzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 292);

    auto g_0_yzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 293);

    auto g_0_zzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 294);

    auto g_0_zzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 295);

    auto g_0_zzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 296);

    auto g_0_zzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 297);

    auto g_0_zzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 298);

    auto g_0_zzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 299);

    auto g_0_zzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 300);

    auto g_0_zzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 301);

    auto g_0_zzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 302);

    auto g_0_zzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 303);

    auto g_0_zzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 304);

    auto g_0_zzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 305);

    auto g_0_zzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 306);

    auto g_0_zzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 307);

    auto g_0_zzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 308);

    auto g_0_zzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 309);

    auto g_0_zzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 310);

    auto g_0_zzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 311);

    auto g_0_zzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 312);

    auto g_0_zzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 313);

    auto g_0_zzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 314);

    /// Set up components of auxilary buffer : SGSI

    auto g_0_xxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi);

    auto g_0_xxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 1);

    auto g_0_xxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 2);

    auto g_0_xxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 3);

    auto g_0_xxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 4);

    auto g_0_xxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 5);

    auto g_0_xxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 6);

    auto g_0_xxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 7);

    auto g_0_xxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 8);

    auto g_0_xxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 9);

    auto g_0_xxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 10);

    auto g_0_xxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 11);

    auto g_0_xxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 12);

    auto g_0_xxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 13);

    auto g_0_xxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 14);

    auto g_0_xxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 15);

    auto g_0_xxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 16);

    auto g_0_xxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 17);

    auto g_0_xxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 18);

    auto g_0_xxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 19);

    auto g_0_xxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 20);

    auto g_0_xxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 21);

    auto g_0_xxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 22);

    auto g_0_xxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 23);

    auto g_0_xxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 24);

    auto g_0_xxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 25);

    auto g_0_xxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 26);

    auto g_0_xxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 27);

    auto g_0_xxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 28);

    auto g_0_xxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 29);

    auto g_0_xxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 30);

    auto g_0_xxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 31);

    auto g_0_xxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 33);

    auto g_0_xxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 34);

    auto g_0_xxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 37);

    auto g_0_xxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 38);

    auto g_0_xxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 42);

    auto g_0_xxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 43);

    auto g_0_xxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 48);

    auto g_0_xxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 49);

    auto g_0_xxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 56);

    auto g_0_xxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 57);

    auto g_0_xxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 58);

    auto g_0_xxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 59);

    auto g_0_xxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 60);

    auto g_0_xxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 61);

    auto g_0_xxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 62);

    auto g_0_xxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 63);

    auto g_0_xxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 64);

    auto g_0_xxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 65);

    auto g_0_xxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 66);

    auto g_0_xxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 67);

    auto g_0_xxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 68);

    auto g_0_xxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 69);

    auto g_0_xxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 70);

    auto g_0_xxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 71);

    auto g_0_xxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 72);

    auto g_0_xxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 73);

    auto g_0_xxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 74);

    auto g_0_xxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 75);

    auto g_0_xxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 76);

    auto g_0_xxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 78);

    auto g_0_xxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 79);

    auto g_0_xxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 80);

    auto g_0_xxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 81);

    auto g_0_xxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 82);

    auto g_0_xxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 83);

    auto g_0_xxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 84);

    auto g_0_xxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 85);

    auto g_0_xxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 86);

    auto g_0_xxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 87);

    auto g_0_xxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 88);

    auto g_0_xxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 89);

    auto g_0_xxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 90);

    auto g_0_xxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 91);

    auto g_0_xxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 92);

    auto g_0_xxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 93);

    auto g_0_xxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 94);

    auto g_0_xxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 95);

    auto g_0_xxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 96);

    auto g_0_xxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 97);

    auto g_0_xxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 98);

    auto g_0_xxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 99);

    auto g_0_xxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 100);

    auto g_0_xxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 101);

    auto g_0_xxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 102);

    auto g_0_xxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 103);

    auto g_0_xxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 104);

    auto g_0_xxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 105);

    auto g_0_xxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 106);

    auto g_0_xxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 107);

    auto g_0_xxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 108);

    auto g_0_xxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 109);

    auto g_0_xxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 110);

    auto g_0_xxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 111);

    auto g_0_xxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 140);

    auto g_0_xxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 141);

    auto g_0_xxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 142);

    auto g_0_xxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 143);

    auto g_0_xxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 144);

    auto g_0_xxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 145);

    auto g_0_xxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 146);

    auto g_0_xxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 147);

    auto g_0_xxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 148);

    auto g_0_xxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 149);

    auto g_0_xxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 150);

    auto g_0_xxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 151);

    auto g_0_xxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 152);

    auto g_0_xxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 153);

    auto g_0_xxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 154);

    auto g_0_xxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 155);

    auto g_0_xxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 156);

    auto g_0_xxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 157);

    auto g_0_xxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 158);

    auto g_0_xxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 159);

    auto g_0_xxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 160);

    auto g_0_xxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 161);

    auto g_0_xxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 162);

    auto g_0_xxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 163);

    auto g_0_xxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 164);

    auto g_0_xxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 165);

    auto g_0_xxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 166);

    auto g_0_xxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 167);

    auto g_0_xyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 168);

    auto g_0_xyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 169);

    auto g_0_xyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 171);

    auto g_0_xyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 172);

    auto g_0_xyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 174);

    auto g_0_xyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 175);

    auto g_0_xyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 176);

    auto g_0_xyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 178);

    auto g_0_xyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 179);

    auto g_0_xyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 180);

    auto g_0_xyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 181);

    auto g_0_xyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 183);

    auto g_0_xyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 184);

    auto g_0_xyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 185);

    auto g_0_xyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 186);

    auto g_0_xyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 187);

    auto g_0_xyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 189);

    auto g_0_xyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 190);

    auto g_0_xyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 191);

    auto g_0_xyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 192);

    auto g_0_xyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 193);

    auto g_0_xyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 194);

    auto g_0_xyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 195);

    auto g_0_xzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 252);

    auto g_0_xzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 254);

    auto g_0_xzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 256);

    auto g_0_xzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 257);

    auto g_0_xzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 259);

    auto g_0_xzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 260);

    auto g_0_xzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 261);

    auto g_0_xzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 263);

    auto g_0_xzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 264);

    auto g_0_xzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 265);

    auto g_0_xzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 266);

    auto g_0_xzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 268);

    auto g_0_xzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 269);

    auto g_0_xzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 270);

    auto g_0_xzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 271);

    auto g_0_xzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 272);

    auto g_0_xzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 273);

    auto g_0_xzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 274);

    auto g_0_xzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 275);

    auto g_0_xzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 276);

    auto g_0_xzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 277);

    auto g_0_xzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 278);

    auto g_0_xzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 279);

    auto g_0_yyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 280);

    auto g_0_yyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 281);

    auto g_0_yyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 282);

    auto g_0_yyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 283);

    auto g_0_yyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 284);

    auto g_0_yyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 285);

    auto g_0_yyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 286);

    auto g_0_yyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 287);

    auto g_0_yyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 288);

    auto g_0_yyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 289);

    auto g_0_yyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 290);

    auto g_0_yyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 291);

    auto g_0_yyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 292);

    auto g_0_yyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 293);

    auto g_0_yyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 294);

    auto g_0_yyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 295);

    auto g_0_yyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 296);

    auto g_0_yyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 297);

    auto g_0_yyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 298);

    auto g_0_yyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 299);

    auto g_0_yyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 300);

    auto g_0_yyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 301);

    auto g_0_yyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 302);

    auto g_0_yyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 303);

    auto g_0_yyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 304);

    auto g_0_yyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 305);

    auto g_0_yyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 306);

    auto g_0_yyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 307);

    auto g_0_yyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 309);

    auto g_0_yyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 310);

    auto g_0_yyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 311);

    auto g_0_yyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 312);

    auto g_0_yyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 313);

    auto g_0_yyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 314);

    auto g_0_yyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 315);

    auto g_0_yyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 316);

    auto g_0_yyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 317);

    auto g_0_yyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 318);

    auto g_0_yyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 319);

    auto g_0_yyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 320);

    auto g_0_yyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 321);

    auto g_0_yyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 322);

    auto g_0_yyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 323);

    auto g_0_yyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 324);

    auto g_0_yyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 325);

    auto g_0_yyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 326);

    auto g_0_yyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 327);

    auto g_0_yyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 328);

    auto g_0_yyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 329);

    auto g_0_yyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 330);

    auto g_0_yyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 331);

    auto g_0_yyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 332);

    auto g_0_yyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 333);

    auto g_0_yyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 334);

    auto g_0_yyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 335);

    auto g_0_yyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 336);

    auto g_0_yyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 337);

    auto g_0_yyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 338);

    auto g_0_yyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 339);

    auto g_0_yyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 340);

    auto g_0_yyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 341);

    auto g_0_yyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 342);

    auto g_0_yyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 343);

    auto g_0_yyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 344);

    auto g_0_yyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 345);

    auto g_0_yyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 346);

    auto g_0_yyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 347);

    auto g_0_yyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 348);

    auto g_0_yyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 349);

    auto g_0_yyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 350);

    auto g_0_yyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 351);

    auto g_0_yyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 352);

    auto g_0_yyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 353);

    auto g_0_yyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 354);

    auto g_0_yyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 355);

    auto g_0_yyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 356);

    auto g_0_yyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 357);

    auto g_0_yyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 358);

    auto g_0_yyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 359);

    auto g_0_yyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 360);

    auto g_0_yyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 361);

    auto g_0_yyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 362);

    auto g_0_yyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 363);

    auto g_0_yzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 364);

    auto g_0_yzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 365);

    auto g_0_yzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 366);

    auto g_0_yzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 367);

    auto g_0_yzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 368);

    auto g_0_yzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 369);

    auto g_0_yzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 370);

    auto g_0_yzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 371);

    auto g_0_yzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 372);

    auto g_0_yzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 373);

    auto g_0_yzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 374);

    auto g_0_yzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 375);

    auto g_0_yzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 376);

    auto g_0_yzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 377);

    auto g_0_yzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 378);

    auto g_0_yzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 379);

    auto g_0_yzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 380);

    auto g_0_yzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 381);

    auto g_0_yzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 382);

    auto g_0_yzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 383);

    auto g_0_yzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 384);

    auto g_0_yzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 385);

    auto g_0_yzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 386);

    auto g_0_yzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 387);

    auto g_0_yzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 388);

    auto g_0_yzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 389);

    auto g_0_yzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 390);

    auto g_0_yzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 391);

    auto g_0_zzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 392);

    auto g_0_zzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 393);

    auto g_0_zzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 394);

    auto g_0_zzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 395);

    auto g_0_zzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 396);

    auto g_0_zzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 397);

    auto g_0_zzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 398);

    auto g_0_zzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 399);

    auto g_0_zzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 400);

    auto g_0_zzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 401);

    auto g_0_zzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 402);

    auto g_0_zzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 403);

    auto g_0_zzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 404);

    auto g_0_zzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 405);

    auto g_0_zzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 406);

    auto g_0_zzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 407);

    auto g_0_zzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 408);

    auto g_0_zzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 409);

    auto g_0_zzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 410);

    auto g_0_zzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 411);

    auto g_0_zzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 412);

    auto g_0_zzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 413);

    auto g_0_zzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 414);

    auto g_0_zzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 415);

    auto g_0_zzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 416);

    auto g_0_zzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 417);

    auto g_0_zzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 418);

    auto g_0_zzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 419);

    /// Set up components of auxilary buffer : SGSI

    auto g_0_xxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi);

    auto g_0_xxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 1);

    auto g_0_xxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 2);

    auto g_0_xxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 3);

    auto g_0_xxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 4);

    auto g_0_xxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 5);

    auto g_0_xxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 6);

    auto g_0_xxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 7);

    auto g_0_xxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 8);

    auto g_0_xxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 9);

    auto g_0_xxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 10);

    auto g_0_xxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 11);

    auto g_0_xxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 12);

    auto g_0_xxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 13);

    auto g_0_xxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 14);

    auto g_0_xxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 15);

    auto g_0_xxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 16);

    auto g_0_xxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 17);

    auto g_0_xxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 18);

    auto g_0_xxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 19);

    auto g_0_xxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 20);

    auto g_0_xxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 21);

    auto g_0_xxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 22);

    auto g_0_xxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 23);

    auto g_0_xxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 24);

    auto g_0_xxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 25);

    auto g_0_xxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 26);

    auto g_0_xxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 27);

    auto g_0_xxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 28);

    auto g_0_xxxy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 29);

    auto g_0_xxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 30);

    auto g_0_xxxy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 31);

    auto g_0_xxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 33);

    auto g_0_xxxy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 34);

    auto g_0_xxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 37);

    auto g_0_xxxy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 38);

    auto g_0_xxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 42);

    auto g_0_xxxy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 43);

    auto g_0_xxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 48);

    auto g_0_xxxy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 49);

    auto g_0_xxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 56);

    auto g_0_xxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 57);

    auto g_0_xxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 58);

    auto g_0_xxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 59);

    auto g_0_xxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 60);

    auto g_0_xxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 61);

    auto g_0_xxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 62);

    auto g_0_xxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 63);

    auto g_0_xxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 64);

    auto g_0_xxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 65);

    auto g_0_xxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 66);

    auto g_0_xxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 67);

    auto g_0_xxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 68);

    auto g_0_xxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 69);

    auto g_0_xxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 70);

    auto g_0_xxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 71);

    auto g_0_xxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 72);

    auto g_0_xxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 73);

    auto g_0_xxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 74);

    auto g_0_xxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 75);

    auto g_0_xxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 76);

    auto g_0_xxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 78);

    auto g_0_xxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 79);

    auto g_0_xxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 80);

    auto g_0_xxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 81);

    auto g_0_xxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 82);

    auto g_0_xxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 83);

    auto g_0_xxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 84);

    auto g_0_xxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 85);

    auto g_0_xxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 86);

    auto g_0_xxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 87);

    auto g_0_xxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 88);

    auto g_0_xxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 89);

    auto g_0_xxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 90);

    auto g_0_xxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 91);

    auto g_0_xxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 92);

    auto g_0_xxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 93);

    auto g_0_xxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 94);

    auto g_0_xxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 95);

    auto g_0_xxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 96);

    auto g_0_xxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 97);

    auto g_0_xxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 98);

    auto g_0_xxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 99);

    auto g_0_xxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 100);

    auto g_0_xxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 101);

    auto g_0_xxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 102);

    auto g_0_xxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 103);

    auto g_0_xxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 104);

    auto g_0_xxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 105);

    auto g_0_xxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 106);

    auto g_0_xxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 107);

    auto g_0_xxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 108);

    auto g_0_xxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 109);

    auto g_0_xxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 110);

    auto g_0_xxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 111);

    auto g_0_xxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 140);

    auto g_0_xxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 141);

    auto g_0_xxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 142);

    auto g_0_xxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 143);

    auto g_0_xxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 144);

    auto g_0_xxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 145);

    auto g_0_xxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 146);

    auto g_0_xxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 147);

    auto g_0_xxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 148);

    auto g_0_xxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 149);

    auto g_0_xxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 150);

    auto g_0_xxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 151);

    auto g_0_xxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 152);

    auto g_0_xxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 153);

    auto g_0_xxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 154);

    auto g_0_xxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 155);

    auto g_0_xxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 156);

    auto g_0_xxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 157);

    auto g_0_xxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 158);

    auto g_0_xxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 159);

    auto g_0_xxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 160);

    auto g_0_xxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 161);

    auto g_0_xxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 162);

    auto g_0_xxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 163);

    auto g_0_xxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 164);

    auto g_0_xxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 165);

    auto g_0_xxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 166);

    auto g_0_xxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 167);

    auto g_0_xyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 168);

    auto g_0_xyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 169);

    auto g_0_xyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 171);

    auto g_0_xyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 172);

    auto g_0_xyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 174);

    auto g_0_xyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 175);

    auto g_0_xyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 176);

    auto g_0_xyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 178);

    auto g_0_xyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 179);

    auto g_0_xyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 180);

    auto g_0_xyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 181);

    auto g_0_xyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 183);

    auto g_0_xyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 184);

    auto g_0_xyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 185);

    auto g_0_xyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 186);

    auto g_0_xyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 187);

    auto g_0_xyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 189);

    auto g_0_xyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 190);

    auto g_0_xyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 191);

    auto g_0_xyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 192);

    auto g_0_xyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 193);

    auto g_0_xyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 194);

    auto g_0_xyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 195);

    auto g_0_xzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 252);

    auto g_0_xzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 254);

    auto g_0_xzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 256);

    auto g_0_xzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 257);

    auto g_0_xzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 259);

    auto g_0_xzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 260);

    auto g_0_xzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 261);

    auto g_0_xzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 263);

    auto g_0_xzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 264);

    auto g_0_xzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 265);

    auto g_0_xzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 266);

    auto g_0_xzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 268);

    auto g_0_xzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 269);

    auto g_0_xzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 270);

    auto g_0_xzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 271);

    auto g_0_xzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 272);

    auto g_0_xzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 273);

    auto g_0_xzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 274);

    auto g_0_xzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 275);

    auto g_0_xzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 276);

    auto g_0_xzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 277);

    auto g_0_xzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 278);

    auto g_0_xzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 279);

    auto g_0_yyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 280);

    auto g_0_yyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 281);

    auto g_0_yyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 282);

    auto g_0_yyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 283);

    auto g_0_yyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 284);

    auto g_0_yyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 285);

    auto g_0_yyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 286);

    auto g_0_yyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 287);

    auto g_0_yyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 288);

    auto g_0_yyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 289);

    auto g_0_yyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 290);

    auto g_0_yyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 291);

    auto g_0_yyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 292);

    auto g_0_yyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 293);

    auto g_0_yyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 294);

    auto g_0_yyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 295);

    auto g_0_yyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 296);

    auto g_0_yyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 297);

    auto g_0_yyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 298);

    auto g_0_yyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 299);

    auto g_0_yyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 300);

    auto g_0_yyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 301);

    auto g_0_yyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 302);

    auto g_0_yyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 303);

    auto g_0_yyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 304);

    auto g_0_yyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 305);

    auto g_0_yyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 306);

    auto g_0_yyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 307);

    auto g_0_yyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 309);

    auto g_0_yyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 310);

    auto g_0_yyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 311);

    auto g_0_yyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 312);

    auto g_0_yyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 313);

    auto g_0_yyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 314);

    auto g_0_yyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 315);

    auto g_0_yyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 316);

    auto g_0_yyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 317);

    auto g_0_yyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 318);

    auto g_0_yyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 319);

    auto g_0_yyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 320);

    auto g_0_yyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 321);

    auto g_0_yyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 322);

    auto g_0_yyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 323);

    auto g_0_yyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 324);

    auto g_0_yyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 325);

    auto g_0_yyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 326);

    auto g_0_yyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 327);

    auto g_0_yyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 328);

    auto g_0_yyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 329);

    auto g_0_yyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 330);

    auto g_0_yyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 331);

    auto g_0_yyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 332);

    auto g_0_yyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 333);

    auto g_0_yyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 334);

    auto g_0_yyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 335);

    auto g_0_yyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 336);

    auto g_0_yyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 337);

    auto g_0_yyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 338);

    auto g_0_yyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 339);

    auto g_0_yyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 340);

    auto g_0_yyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 341);

    auto g_0_yyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 342);

    auto g_0_yyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 343);

    auto g_0_yyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 344);

    auto g_0_yyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 345);

    auto g_0_yyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 346);

    auto g_0_yyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 347);

    auto g_0_yyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 348);

    auto g_0_yyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 349);

    auto g_0_yyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 350);

    auto g_0_yyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 351);

    auto g_0_yyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 352);

    auto g_0_yyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 353);

    auto g_0_yyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 354);

    auto g_0_yyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 355);

    auto g_0_yyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 356);

    auto g_0_yyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 357);

    auto g_0_yyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 358);

    auto g_0_yyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 359);

    auto g_0_yyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 360);

    auto g_0_yyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 361);

    auto g_0_yyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 362);

    auto g_0_yyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 363);

    auto g_0_yzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 364);

    auto g_0_yzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 365);

    auto g_0_yzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 366);

    auto g_0_yzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 367);

    auto g_0_yzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 368);

    auto g_0_yzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 369);

    auto g_0_yzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 370);

    auto g_0_yzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 371);

    auto g_0_yzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 372);

    auto g_0_yzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 373);

    auto g_0_yzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 374);

    auto g_0_yzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 375);

    auto g_0_yzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 376);

    auto g_0_yzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 377);

    auto g_0_yzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 378);

    auto g_0_yzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 379);

    auto g_0_yzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 380);

    auto g_0_yzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 381);

    auto g_0_yzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 382);

    auto g_0_yzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 383);

    auto g_0_yzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 384);

    auto g_0_yzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 385);

    auto g_0_yzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 386);

    auto g_0_yzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 387);

    auto g_0_yzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 388);

    auto g_0_yzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 389);

    auto g_0_yzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 390);

    auto g_0_yzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 391);

    auto g_0_zzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 392);

    auto g_0_zzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 393);

    auto g_0_zzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 394);

    auto g_0_zzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 395);

    auto g_0_zzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 396);

    auto g_0_zzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 397);

    auto g_0_zzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 398);

    auto g_0_zzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 399);

    auto g_0_zzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 400);

    auto g_0_zzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 401);

    auto g_0_zzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 402);

    auto g_0_zzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 403);

    auto g_0_zzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 404);

    auto g_0_zzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 405);

    auto g_0_zzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 406);

    auto g_0_zzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 407);

    auto g_0_zzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 408);

    auto g_0_zzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 409);

    auto g_0_zzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 410);

    auto g_0_zzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 411);

    auto g_0_zzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 412);

    auto g_0_zzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 413);

    auto g_0_zzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 414);

    auto g_0_zzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 415);

    auto g_0_zzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 416);

    auto g_0_zzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 417);

    auto g_0_zzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 418);

    auto g_0_zzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 419);

    /// Set up 0-28 components of targeted buffer : SHSI

    auto g_0_xxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi);

    auto g_0_xxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 1);

    auto g_0_xxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 2);

    auto g_0_xxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 3);

    auto g_0_xxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 4);

    auto g_0_xxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 5);

    auto g_0_xxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 6);

    auto g_0_xxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 7);

    auto g_0_xxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 8);

    auto g_0_xxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 9);

    auto g_0_xxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 10);

    auto g_0_xxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 11);

    auto g_0_xxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 12);

    auto g_0_xxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 13);

    auto g_0_xxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 14);

    auto g_0_xxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 15);

    auto g_0_xxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 16);

    auto g_0_xxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 17);

    auto g_0_xxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 18);

    auto g_0_xxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 19);

    auto g_0_xxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 20);

    auto g_0_xxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 21);

    auto g_0_xxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 22);

    auto g_0_xxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 23);

    auto g_0_xxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 24);

    auto g_0_xxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 25);

    auto g_0_xxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 26);

    auto g_0_xxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 27);

#pragma omp simd aligned(g_0_xxx_0_xxxxxx_0,       \
                             g_0_xxx_0_xxxxxx_1,   \
                             g_0_xxx_0_xxxxxy_0,   \
                             g_0_xxx_0_xxxxxy_1,   \
                             g_0_xxx_0_xxxxxz_0,   \
                             g_0_xxx_0_xxxxxz_1,   \
                             g_0_xxx_0_xxxxyy_0,   \
                             g_0_xxx_0_xxxxyy_1,   \
                             g_0_xxx_0_xxxxyz_0,   \
                             g_0_xxx_0_xxxxyz_1,   \
                             g_0_xxx_0_xxxxzz_0,   \
                             g_0_xxx_0_xxxxzz_1,   \
                             g_0_xxx_0_xxxyyy_0,   \
                             g_0_xxx_0_xxxyyy_1,   \
                             g_0_xxx_0_xxxyyz_0,   \
                             g_0_xxx_0_xxxyyz_1,   \
                             g_0_xxx_0_xxxyzz_0,   \
                             g_0_xxx_0_xxxyzz_1,   \
                             g_0_xxx_0_xxxzzz_0,   \
                             g_0_xxx_0_xxxzzz_1,   \
                             g_0_xxx_0_xxyyyy_0,   \
                             g_0_xxx_0_xxyyyy_1,   \
                             g_0_xxx_0_xxyyyz_0,   \
                             g_0_xxx_0_xxyyyz_1,   \
                             g_0_xxx_0_xxyyzz_0,   \
                             g_0_xxx_0_xxyyzz_1,   \
                             g_0_xxx_0_xxyzzz_0,   \
                             g_0_xxx_0_xxyzzz_1,   \
                             g_0_xxx_0_xxzzzz_0,   \
                             g_0_xxx_0_xxzzzz_1,   \
                             g_0_xxx_0_xyyyyy_0,   \
                             g_0_xxx_0_xyyyyy_1,   \
                             g_0_xxx_0_xyyyyz_0,   \
                             g_0_xxx_0_xyyyyz_1,   \
                             g_0_xxx_0_xyyyzz_0,   \
                             g_0_xxx_0_xyyyzz_1,   \
                             g_0_xxx_0_xyyzzz_0,   \
                             g_0_xxx_0_xyyzzz_1,   \
                             g_0_xxx_0_xyzzzz_0,   \
                             g_0_xxx_0_xyzzzz_1,   \
                             g_0_xxx_0_xzzzzz_0,   \
                             g_0_xxx_0_xzzzzz_1,   \
                             g_0_xxx_0_yyyyyy_0,   \
                             g_0_xxx_0_yyyyyy_1,   \
                             g_0_xxx_0_yyyyyz_0,   \
                             g_0_xxx_0_yyyyyz_1,   \
                             g_0_xxx_0_yyyyzz_0,   \
                             g_0_xxx_0_yyyyzz_1,   \
                             g_0_xxx_0_yyyzzz_0,   \
                             g_0_xxx_0_yyyzzz_1,   \
                             g_0_xxx_0_yyzzzz_0,   \
                             g_0_xxx_0_yyzzzz_1,   \
                             g_0_xxx_0_yzzzzz_0,   \
                             g_0_xxx_0_yzzzzz_1,   \
                             g_0_xxx_0_zzzzzz_0,   \
                             g_0_xxx_0_zzzzzz_1,   \
                             g_0_xxxx_0_xxxxx_1,   \
                             g_0_xxxx_0_xxxxxx_0,  \
                             g_0_xxxx_0_xxxxxx_1,  \
                             g_0_xxxx_0_xxxxxy_0,  \
                             g_0_xxxx_0_xxxxxy_1,  \
                             g_0_xxxx_0_xxxxxz_0,  \
                             g_0_xxxx_0_xxxxxz_1,  \
                             g_0_xxxx_0_xxxxy_1,   \
                             g_0_xxxx_0_xxxxyy_0,  \
                             g_0_xxxx_0_xxxxyy_1,  \
                             g_0_xxxx_0_xxxxyz_0,  \
                             g_0_xxxx_0_xxxxyz_1,  \
                             g_0_xxxx_0_xxxxz_1,   \
                             g_0_xxxx_0_xxxxzz_0,  \
                             g_0_xxxx_0_xxxxzz_1,  \
                             g_0_xxxx_0_xxxyy_1,   \
                             g_0_xxxx_0_xxxyyy_0,  \
                             g_0_xxxx_0_xxxyyy_1,  \
                             g_0_xxxx_0_xxxyyz_0,  \
                             g_0_xxxx_0_xxxyyz_1,  \
                             g_0_xxxx_0_xxxyz_1,   \
                             g_0_xxxx_0_xxxyzz_0,  \
                             g_0_xxxx_0_xxxyzz_1,  \
                             g_0_xxxx_0_xxxzz_1,   \
                             g_0_xxxx_0_xxxzzz_0,  \
                             g_0_xxxx_0_xxxzzz_1,  \
                             g_0_xxxx_0_xxyyy_1,   \
                             g_0_xxxx_0_xxyyyy_0,  \
                             g_0_xxxx_0_xxyyyy_1,  \
                             g_0_xxxx_0_xxyyyz_0,  \
                             g_0_xxxx_0_xxyyyz_1,  \
                             g_0_xxxx_0_xxyyz_1,   \
                             g_0_xxxx_0_xxyyzz_0,  \
                             g_0_xxxx_0_xxyyzz_1,  \
                             g_0_xxxx_0_xxyzz_1,   \
                             g_0_xxxx_0_xxyzzz_0,  \
                             g_0_xxxx_0_xxyzzz_1,  \
                             g_0_xxxx_0_xxzzz_1,   \
                             g_0_xxxx_0_xxzzzz_0,  \
                             g_0_xxxx_0_xxzzzz_1,  \
                             g_0_xxxx_0_xyyyy_1,   \
                             g_0_xxxx_0_xyyyyy_0,  \
                             g_0_xxxx_0_xyyyyy_1,  \
                             g_0_xxxx_0_xyyyyz_0,  \
                             g_0_xxxx_0_xyyyyz_1,  \
                             g_0_xxxx_0_xyyyz_1,   \
                             g_0_xxxx_0_xyyyzz_0,  \
                             g_0_xxxx_0_xyyyzz_1,  \
                             g_0_xxxx_0_xyyzz_1,   \
                             g_0_xxxx_0_xyyzzz_0,  \
                             g_0_xxxx_0_xyyzzz_1,  \
                             g_0_xxxx_0_xyzzz_1,   \
                             g_0_xxxx_0_xyzzzz_0,  \
                             g_0_xxxx_0_xyzzzz_1,  \
                             g_0_xxxx_0_xzzzz_1,   \
                             g_0_xxxx_0_xzzzzz_0,  \
                             g_0_xxxx_0_xzzzzz_1,  \
                             g_0_xxxx_0_yyyyy_1,   \
                             g_0_xxxx_0_yyyyyy_0,  \
                             g_0_xxxx_0_yyyyyy_1,  \
                             g_0_xxxx_0_yyyyyz_0,  \
                             g_0_xxxx_0_yyyyyz_1,  \
                             g_0_xxxx_0_yyyyz_1,   \
                             g_0_xxxx_0_yyyyzz_0,  \
                             g_0_xxxx_0_yyyyzz_1,  \
                             g_0_xxxx_0_yyyzz_1,   \
                             g_0_xxxx_0_yyyzzz_0,  \
                             g_0_xxxx_0_yyyzzz_1,  \
                             g_0_xxxx_0_yyzzz_1,   \
                             g_0_xxxx_0_yyzzzz_0,  \
                             g_0_xxxx_0_yyzzzz_1,  \
                             g_0_xxxx_0_yzzzz_1,   \
                             g_0_xxxx_0_yzzzzz_0,  \
                             g_0_xxxx_0_yzzzzz_1,  \
                             g_0_xxxx_0_zzzzz_1,   \
                             g_0_xxxx_0_zzzzzz_0,  \
                             g_0_xxxx_0_zzzzzz_1,  \
                             g_0_xxxxx_0_xxxxxx_0, \
                             g_0_xxxxx_0_xxxxxy_0, \
                             g_0_xxxxx_0_xxxxxz_0, \
                             g_0_xxxxx_0_xxxxyy_0, \
                             g_0_xxxxx_0_xxxxyz_0, \
                             g_0_xxxxx_0_xxxxzz_0, \
                             g_0_xxxxx_0_xxxyyy_0, \
                             g_0_xxxxx_0_xxxyyz_0, \
                             g_0_xxxxx_0_xxxyzz_0, \
                             g_0_xxxxx_0_xxxzzz_0, \
                             g_0_xxxxx_0_xxyyyy_0, \
                             g_0_xxxxx_0_xxyyyz_0, \
                             g_0_xxxxx_0_xxyyzz_0, \
                             g_0_xxxxx_0_xxyzzz_0, \
                             g_0_xxxxx_0_xxzzzz_0, \
                             g_0_xxxxx_0_xyyyyy_0, \
                             g_0_xxxxx_0_xyyyyz_0, \
                             g_0_xxxxx_0_xyyyzz_0, \
                             g_0_xxxxx_0_xyyzzz_0, \
                             g_0_xxxxx_0_xyzzzz_0, \
                             g_0_xxxxx_0_xzzzzz_0, \
                             g_0_xxxxx_0_yyyyyy_0, \
                             g_0_xxxxx_0_yyyyyz_0, \
                             g_0_xxxxx_0_yyyyzz_0, \
                             g_0_xxxxx_0_yyyzzz_0, \
                             g_0_xxxxx_0_yyzzzz_0, \
                             g_0_xxxxx_0_yzzzzz_0, \
                             g_0_xxxxx_0_zzzzzz_0, \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_xxxxxx_0[i] = 4.0 * g_0_xxx_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxx_1[i] * fti_ab_0 +
                                  6.0 * g_0_xxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxx_0[i] * pb_x + g_0_xxxx_0_xxxxxx_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxy_0[i] = 4.0 * g_0_xxx_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxy_1[i] * fti_ab_0 +
                                  5.0 * g_0_xxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxy_0[i] * pb_x + g_0_xxxx_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxz_0[i] = 4.0 * g_0_xxx_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxz_1[i] * fti_ab_0 +
                                  5.0 * g_0_xxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxz_0[i] * pb_x + g_0_xxxx_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyy_0[i] = 4.0 * g_0_xxx_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyy_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyy_0[i] * pb_x + g_0_xxxx_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyz_0[i] = 4.0 * g_0_xxx_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyz_0[i] * pb_x + g_0_xxxx_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxzz_0[i] = 4.0 * g_0_xxx_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxzz_0[i] * pb_x + g_0_xxxx_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyy_0[i] = 4.0 * g_0_xxx_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyy_0[i] * pb_x + g_0_xxxx_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyz_0[i] = 4.0 * g_0_xxx_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyz_0[i] * pb_x + g_0_xxxx_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyzz_0[i] = 4.0 * g_0_xxx_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzz_0[i] * pb_x + g_0_xxxx_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxzzz_0[i] = 4.0 * g_0_xxx_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxzzz_0[i] * pb_x + g_0_xxxx_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyy_0[i] = 4.0 * g_0_xxx_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyy_0[i] * pb_x + g_0_xxxx_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyz_0[i] = 4.0 * g_0_xxx_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyz_0[i] * pb_x + g_0_xxxx_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyzz_0[i] = 4.0 * g_0_xxx_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzz_0[i] * pb_x + g_0_xxxx_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyzzz_0[i] = 4.0 * g_0_xxx_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzz_0[i] * pb_x + g_0_xxxx_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxzzzz_0[i] = 4.0 * g_0_xxx_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzzzz_0[i] * pb_x + g_0_xxxx_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyy_0[i] = 4.0 * g_0_xxx_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyy_1[i] * fi_abcd_0 +
                                  g_0_xxxx_0_xyyyyy_0[i] * pb_x + g_0_xxxx_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyz_0[i] = 4.0 * g_0_xxx_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyz_1[i] * fi_abcd_0 +
                                  g_0_xxxx_0_xyyyyz_0[i] * pb_x + g_0_xxxx_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyzz_0[i] = 4.0 * g_0_xxx_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyzz_1[i] * fi_abcd_0 +
                                  g_0_xxxx_0_xyyyzz_0[i] * pb_x + g_0_xxxx_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyzzz_0[i] = 4.0 * g_0_xxx_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzzz_1[i] * fi_abcd_0 +
                                  g_0_xxxx_0_xyyzzz_0[i] * pb_x + g_0_xxxx_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyzzzz_0[i] = 4.0 * g_0_xxx_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzzz_1[i] * fi_abcd_0 +
                                  g_0_xxxx_0_xyzzzz_0[i] * pb_x + g_0_xxxx_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xzzzzz_0[i] = 4.0 * g_0_xxx_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzzz_1[i] * fi_abcd_0 +
                                  g_0_xxxx_0_xzzzzz_0[i] * pb_x + g_0_xxxx_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyy_0[i] = 4.0 * g_0_xxx_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyy_0[i] * pb_x +
                                  g_0_xxxx_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyz_0[i] = 4.0 * g_0_xxx_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyz_0[i] * pb_x +
                                  g_0_xxxx_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyzz_0[i] = 4.0 * g_0_xxx_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyzz_0[i] * pb_x +
                                  g_0_xxxx_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyzzz_0[i] = 4.0 * g_0_xxx_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyzzz_0[i] * pb_x +
                                  g_0_xxxx_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyzzzz_0[i] = 4.0 * g_0_xxx_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzzzz_0[i] * pb_x +
                                  g_0_xxxx_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yzzzzz_0[i] = 4.0 * g_0_xxx_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzzzz_0[i] * pb_x +
                                  g_0_xxxx_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_zzzzzz_0[i] = 4.0 * g_0_xxx_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzzzz_0[i] * pb_x +
                                  g_0_xxxx_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 28-56 components of targeted buffer : SHSI

    auto g_0_xxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 28);

    auto g_0_xxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 29);

    auto g_0_xxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 30);

    auto g_0_xxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 31);

    auto g_0_xxxxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 32);

    auto g_0_xxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 33);

    auto g_0_xxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 34);

    auto g_0_xxxxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 35);

    auto g_0_xxxxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 36);

    auto g_0_xxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 37);

    auto g_0_xxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 38);

    auto g_0_xxxxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 39);

    auto g_0_xxxxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 40);

    auto g_0_xxxxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 41);

    auto g_0_xxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 42);

    auto g_0_xxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 43);

    auto g_0_xxxxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 44);

    auto g_0_xxxxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 45);

    auto g_0_xxxxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 46);

    auto g_0_xxxxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 47);

    auto g_0_xxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 48);

    auto g_0_xxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 49);

    auto g_0_xxxxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 50);

    auto g_0_xxxxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 51);

    auto g_0_xxxxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 52);

    auto g_0_xxxxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 53);

    auto g_0_xxxxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 54);

    auto g_0_xxxxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 55);

#pragma omp simd aligned(g_0_xxxx_0_xxxxx_1,       \
                             g_0_xxxx_0_xxxxxx_0,  \
                             g_0_xxxx_0_xxxxxx_1,  \
                             g_0_xxxx_0_xxxxxy_0,  \
                             g_0_xxxx_0_xxxxxy_1,  \
                             g_0_xxxx_0_xxxxxz_0,  \
                             g_0_xxxx_0_xxxxxz_1,  \
                             g_0_xxxx_0_xxxxy_1,   \
                             g_0_xxxx_0_xxxxyy_0,  \
                             g_0_xxxx_0_xxxxyy_1,  \
                             g_0_xxxx_0_xxxxyz_0,  \
                             g_0_xxxx_0_xxxxyz_1,  \
                             g_0_xxxx_0_xxxxz_1,   \
                             g_0_xxxx_0_xxxxzz_0,  \
                             g_0_xxxx_0_xxxxzz_1,  \
                             g_0_xxxx_0_xxxyy_1,   \
                             g_0_xxxx_0_xxxyyy_0,  \
                             g_0_xxxx_0_xxxyyy_1,  \
                             g_0_xxxx_0_xxxyyz_0,  \
                             g_0_xxxx_0_xxxyyz_1,  \
                             g_0_xxxx_0_xxxyz_1,   \
                             g_0_xxxx_0_xxxyzz_0,  \
                             g_0_xxxx_0_xxxyzz_1,  \
                             g_0_xxxx_0_xxxzz_1,   \
                             g_0_xxxx_0_xxxzzz_0,  \
                             g_0_xxxx_0_xxxzzz_1,  \
                             g_0_xxxx_0_xxyyy_1,   \
                             g_0_xxxx_0_xxyyyy_0,  \
                             g_0_xxxx_0_xxyyyy_1,  \
                             g_0_xxxx_0_xxyyyz_0,  \
                             g_0_xxxx_0_xxyyyz_1,  \
                             g_0_xxxx_0_xxyyz_1,   \
                             g_0_xxxx_0_xxyyzz_0,  \
                             g_0_xxxx_0_xxyyzz_1,  \
                             g_0_xxxx_0_xxyzz_1,   \
                             g_0_xxxx_0_xxyzzz_0,  \
                             g_0_xxxx_0_xxyzzz_1,  \
                             g_0_xxxx_0_xxzzz_1,   \
                             g_0_xxxx_0_xxzzzz_0,  \
                             g_0_xxxx_0_xxzzzz_1,  \
                             g_0_xxxx_0_xyyyy_1,   \
                             g_0_xxxx_0_xyyyyy_0,  \
                             g_0_xxxx_0_xyyyyy_1,  \
                             g_0_xxxx_0_xyyyyz_0,  \
                             g_0_xxxx_0_xyyyyz_1,  \
                             g_0_xxxx_0_xyyyz_1,   \
                             g_0_xxxx_0_xyyyzz_0,  \
                             g_0_xxxx_0_xyyyzz_1,  \
                             g_0_xxxx_0_xyyzz_1,   \
                             g_0_xxxx_0_xyyzzz_0,  \
                             g_0_xxxx_0_xyyzzz_1,  \
                             g_0_xxxx_0_xyzzz_1,   \
                             g_0_xxxx_0_xyzzzz_0,  \
                             g_0_xxxx_0_xyzzzz_1,  \
                             g_0_xxxx_0_xzzzz_1,   \
                             g_0_xxxx_0_xzzzzz_0,  \
                             g_0_xxxx_0_xzzzzz_1,  \
                             g_0_xxxx_0_yyyyy_1,   \
                             g_0_xxxx_0_yyyyyy_0,  \
                             g_0_xxxx_0_yyyyyy_1,  \
                             g_0_xxxx_0_yyyyyz_0,  \
                             g_0_xxxx_0_yyyyyz_1,  \
                             g_0_xxxx_0_yyyyz_1,   \
                             g_0_xxxx_0_yyyyzz_0,  \
                             g_0_xxxx_0_yyyyzz_1,  \
                             g_0_xxxx_0_yyyzz_1,   \
                             g_0_xxxx_0_yyyzzz_0,  \
                             g_0_xxxx_0_yyyzzz_1,  \
                             g_0_xxxx_0_yyzzz_1,   \
                             g_0_xxxx_0_yyzzzz_0,  \
                             g_0_xxxx_0_yyzzzz_1,  \
                             g_0_xxxx_0_yzzzz_1,   \
                             g_0_xxxx_0_yzzzzz_0,  \
                             g_0_xxxx_0_yzzzzz_1,  \
                             g_0_xxxx_0_zzzzz_1,   \
                             g_0_xxxx_0_zzzzzz_0,  \
                             g_0_xxxx_0_zzzzzz_1,  \
                             g_0_xxxxy_0_xxxxxx_0, \
                             g_0_xxxxy_0_xxxxxy_0, \
                             g_0_xxxxy_0_xxxxxz_0, \
                             g_0_xxxxy_0_xxxxyy_0, \
                             g_0_xxxxy_0_xxxxyz_0, \
                             g_0_xxxxy_0_xxxxzz_0, \
                             g_0_xxxxy_0_xxxyyy_0, \
                             g_0_xxxxy_0_xxxyyz_0, \
                             g_0_xxxxy_0_xxxyzz_0, \
                             g_0_xxxxy_0_xxxzzz_0, \
                             g_0_xxxxy_0_xxyyyy_0, \
                             g_0_xxxxy_0_xxyyyz_0, \
                             g_0_xxxxy_0_xxyyzz_0, \
                             g_0_xxxxy_0_xxyzzz_0, \
                             g_0_xxxxy_0_xxzzzz_0, \
                             g_0_xxxxy_0_xyyyyy_0, \
                             g_0_xxxxy_0_xyyyyz_0, \
                             g_0_xxxxy_0_xyyyzz_0, \
                             g_0_xxxxy_0_xyyzzz_0, \
                             g_0_xxxxy_0_xyzzzz_0, \
                             g_0_xxxxy_0_xzzzzz_0, \
                             g_0_xxxxy_0_yyyyyy_0, \
                             g_0_xxxxy_0_yyyyyz_0, \
                             g_0_xxxxy_0_yyyyzz_0, \
                             g_0_xxxxy_0_yyyzzz_0, \
                             g_0_xxxxy_0_yyzzzz_0, \
                             g_0_xxxxy_0_yzzzzz_0, \
                             g_0_xxxxy_0_zzzzzz_0, \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxy_0_xxxxxx_0[i] = g_0_xxxx_0_xxxxxx_0[i] * pb_y + g_0_xxxx_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxy_0[i] = g_0_xxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxy_0[i] * pb_y + g_0_xxxx_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxz_0[i] = g_0_xxxx_0_xxxxxz_0[i] * pb_y + g_0_xxxx_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyy_0[i] = 2.0 * g_0_xxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyy_0[i] * pb_y + g_0_xxxx_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyz_0[i] = g_0_xxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyz_0[i] * pb_y + g_0_xxxx_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxzz_0[i] = g_0_xxxx_0_xxxxzz_0[i] * pb_y + g_0_xxxx_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyy_0[i] = 3.0 * g_0_xxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyy_0[i] * pb_y + g_0_xxxx_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyz_0[i] = 2.0 * g_0_xxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyz_0[i] * pb_y + g_0_xxxx_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyzz_0[i] = g_0_xxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzz_0[i] * pb_y + g_0_xxxx_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxzzz_0[i] = g_0_xxxx_0_xxxzzz_0[i] * pb_y + g_0_xxxx_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyy_0[i] = 4.0 * g_0_xxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyy_0[i] * pb_y + g_0_xxxx_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyz_0[i] = 3.0 * g_0_xxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyz_0[i] * pb_y + g_0_xxxx_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyzz_0[i] = 2.0 * g_0_xxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzz_0[i] * pb_y + g_0_xxxx_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyzzz_0[i] = g_0_xxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzz_0[i] * pb_y + g_0_xxxx_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxzzzz_0[i] = g_0_xxxx_0_xxzzzz_0[i] * pb_y + g_0_xxxx_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyy_0[i] = 5.0 * g_0_xxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyy_0[i] * pb_y + g_0_xxxx_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyz_0[i] = 4.0 * g_0_xxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyz_0[i] * pb_y + g_0_xxxx_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyzz_0[i] = 3.0 * g_0_xxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyzz_0[i] * pb_y + g_0_xxxx_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyzzz_0[i] = 2.0 * g_0_xxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzzz_0[i] * pb_y + g_0_xxxx_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyzzzz_0[i] = g_0_xxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzzz_0[i] * pb_y + g_0_xxxx_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xzzzzz_0[i] = g_0_xxxx_0_xzzzzz_0[i] * pb_y + g_0_xxxx_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyy_0[i] = 6.0 * g_0_xxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyy_0[i] * pb_y + g_0_xxxx_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyz_0[i] = 5.0 * g_0_xxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyz_0[i] * pb_y + g_0_xxxx_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyzz_0[i] = 4.0 * g_0_xxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyzz_0[i] * pb_y + g_0_xxxx_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyzzz_0[i] = 3.0 * g_0_xxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyzzz_0[i] * pb_y + g_0_xxxx_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyzzzz_0[i] = 2.0 * g_0_xxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzzzz_0[i] * pb_y + g_0_xxxx_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yzzzzz_0[i] = g_0_xxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzzzz_0[i] * pb_y + g_0_xxxx_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_zzzzzz_0[i] = g_0_xxxx_0_zzzzzz_0[i] * pb_y + g_0_xxxx_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 56-84 components of targeted buffer : SHSI

    auto g_0_xxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 56);

    auto g_0_xxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 57);

    auto g_0_xxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 58);

    auto g_0_xxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 59);

    auto g_0_xxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 60);

    auto g_0_xxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 61);

    auto g_0_xxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 62);

    auto g_0_xxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 63);

    auto g_0_xxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 64);

    auto g_0_xxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 65);

    auto g_0_xxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 66);

    auto g_0_xxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 67);

    auto g_0_xxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 68);

    auto g_0_xxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 69);

    auto g_0_xxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 70);

    auto g_0_xxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 71);

    auto g_0_xxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 72);

    auto g_0_xxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 73);

    auto g_0_xxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 74);

    auto g_0_xxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 75);

    auto g_0_xxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 76);

    auto g_0_xxxxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 77);

    auto g_0_xxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 78);

    auto g_0_xxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 79);

    auto g_0_xxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 80);

    auto g_0_xxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 81);

    auto g_0_xxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 82);

    auto g_0_xxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 83);

#pragma omp simd aligned(g_0_xxxx_0_xxxxx_1,       \
                             g_0_xxxx_0_xxxxxx_0,  \
                             g_0_xxxx_0_xxxxxx_1,  \
                             g_0_xxxx_0_xxxxxy_0,  \
                             g_0_xxxx_0_xxxxxy_1,  \
                             g_0_xxxx_0_xxxxxz_0,  \
                             g_0_xxxx_0_xxxxxz_1,  \
                             g_0_xxxx_0_xxxxy_1,   \
                             g_0_xxxx_0_xxxxyy_0,  \
                             g_0_xxxx_0_xxxxyy_1,  \
                             g_0_xxxx_0_xxxxyz_0,  \
                             g_0_xxxx_0_xxxxyz_1,  \
                             g_0_xxxx_0_xxxxz_1,   \
                             g_0_xxxx_0_xxxxzz_0,  \
                             g_0_xxxx_0_xxxxzz_1,  \
                             g_0_xxxx_0_xxxyy_1,   \
                             g_0_xxxx_0_xxxyyy_0,  \
                             g_0_xxxx_0_xxxyyy_1,  \
                             g_0_xxxx_0_xxxyyz_0,  \
                             g_0_xxxx_0_xxxyyz_1,  \
                             g_0_xxxx_0_xxxyz_1,   \
                             g_0_xxxx_0_xxxyzz_0,  \
                             g_0_xxxx_0_xxxyzz_1,  \
                             g_0_xxxx_0_xxxzz_1,   \
                             g_0_xxxx_0_xxxzzz_0,  \
                             g_0_xxxx_0_xxxzzz_1,  \
                             g_0_xxxx_0_xxyyy_1,   \
                             g_0_xxxx_0_xxyyyy_0,  \
                             g_0_xxxx_0_xxyyyy_1,  \
                             g_0_xxxx_0_xxyyyz_0,  \
                             g_0_xxxx_0_xxyyyz_1,  \
                             g_0_xxxx_0_xxyyz_1,   \
                             g_0_xxxx_0_xxyyzz_0,  \
                             g_0_xxxx_0_xxyyzz_1,  \
                             g_0_xxxx_0_xxyzz_1,   \
                             g_0_xxxx_0_xxyzzz_0,  \
                             g_0_xxxx_0_xxyzzz_1,  \
                             g_0_xxxx_0_xxzzz_1,   \
                             g_0_xxxx_0_xxzzzz_0,  \
                             g_0_xxxx_0_xxzzzz_1,  \
                             g_0_xxxx_0_xyyyy_1,   \
                             g_0_xxxx_0_xyyyyy_0,  \
                             g_0_xxxx_0_xyyyyy_1,  \
                             g_0_xxxx_0_xyyyyz_0,  \
                             g_0_xxxx_0_xyyyyz_1,  \
                             g_0_xxxx_0_xyyyz_1,   \
                             g_0_xxxx_0_xyyyzz_0,  \
                             g_0_xxxx_0_xyyyzz_1,  \
                             g_0_xxxx_0_xyyzz_1,   \
                             g_0_xxxx_0_xyyzzz_0,  \
                             g_0_xxxx_0_xyyzzz_1,  \
                             g_0_xxxx_0_xyzzz_1,   \
                             g_0_xxxx_0_xyzzzz_0,  \
                             g_0_xxxx_0_xyzzzz_1,  \
                             g_0_xxxx_0_xzzzz_1,   \
                             g_0_xxxx_0_xzzzzz_0,  \
                             g_0_xxxx_0_xzzzzz_1,  \
                             g_0_xxxx_0_yyyyy_1,   \
                             g_0_xxxx_0_yyyyyy_0,  \
                             g_0_xxxx_0_yyyyyy_1,  \
                             g_0_xxxx_0_yyyyyz_0,  \
                             g_0_xxxx_0_yyyyyz_1,  \
                             g_0_xxxx_0_yyyyz_1,   \
                             g_0_xxxx_0_yyyyzz_0,  \
                             g_0_xxxx_0_yyyyzz_1,  \
                             g_0_xxxx_0_yyyzz_1,   \
                             g_0_xxxx_0_yyyzzz_0,  \
                             g_0_xxxx_0_yyyzzz_1,  \
                             g_0_xxxx_0_yyzzz_1,   \
                             g_0_xxxx_0_yyzzzz_0,  \
                             g_0_xxxx_0_yyzzzz_1,  \
                             g_0_xxxx_0_yzzzz_1,   \
                             g_0_xxxx_0_yzzzzz_0,  \
                             g_0_xxxx_0_yzzzzz_1,  \
                             g_0_xxxx_0_zzzzz_1,   \
                             g_0_xxxx_0_zzzzzz_0,  \
                             g_0_xxxx_0_zzzzzz_1,  \
                             g_0_xxxxz_0_xxxxxx_0, \
                             g_0_xxxxz_0_xxxxxy_0, \
                             g_0_xxxxz_0_xxxxxz_0, \
                             g_0_xxxxz_0_xxxxyy_0, \
                             g_0_xxxxz_0_xxxxyz_0, \
                             g_0_xxxxz_0_xxxxzz_0, \
                             g_0_xxxxz_0_xxxyyy_0, \
                             g_0_xxxxz_0_xxxyyz_0, \
                             g_0_xxxxz_0_xxxyzz_0, \
                             g_0_xxxxz_0_xxxzzz_0, \
                             g_0_xxxxz_0_xxyyyy_0, \
                             g_0_xxxxz_0_xxyyyz_0, \
                             g_0_xxxxz_0_xxyyzz_0, \
                             g_0_xxxxz_0_xxyzzz_0, \
                             g_0_xxxxz_0_xxzzzz_0, \
                             g_0_xxxxz_0_xyyyyy_0, \
                             g_0_xxxxz_0_xyyyyz_0, \
                             g_0_xxxxz_0_xyyyzz_0, \
                             g_0_xxxxz_0_xyyzzz_0, \
                             g_0_xxxxz_0_xyzzzz_0, \
                             g_0_xxxxz_0_xzzzzz_0, \
                             g_0_xxxxz_0_yyyyyy_0, \
                             g_0_xxxxz_0_yyyyyz_0, \
                             g_0_xxxxz_0_yyyyzz_0, \
                             g_0_xxxxz_0_yyyzzz_0, \
                             g_0_xxxxz_0_yyzzzz_0, \
                             g_0_xxxxz_0_yzzzzz_0, \
                             g_0_xxxxz_0_zzzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxz_0_xxxxxx_0[i] = g_0_xxxx_0_xxxxxx_0[i] * pb_z + g_0_xxxx_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxy_0[i] = g_0_xxxx_0_xxxxxy_0[i] * pb_z + g_0_xxxx_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxz_0[i] = g_0_xxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxz_0[i] * pb_z + g_0_xxxx_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyy_0[i] = g_0_xxxx_0_xxxxyy_0[i] * pb_z + g_0_xxxx_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyz_0[i] = g_0_xxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyz_0[i] * pb_z + g_0_xxxx_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxzz_0[i] = 2.0 * g_0_xxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxzz_0[i] * pb_z + g_0_xxxx_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyy_0[i] = g_0_xxxx_0_xxxyyy_0[i] * pb_z + g_0_xxxx_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyz_0[i] = g_0_xxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyz_0[i] * pb_z + g_0_xxxx_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyzz_0[i] = 2.0 * g_0_xxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzz_0[i] * pb_z + g_0_xxxx_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxzzz_0[i] = 3.0 * g_0_xxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxzzz_0[i] * pb_z + g_0_xxxx_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyy_0[i] = g_0_xxxx_0_xxyyyy_0[i] * pb_z + g_0_xxxx_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyz_0[i] = g_0_xxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyz_0[i] * pb_z + g_0_xxxx_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyzz_0[i] = 2.0 * g_0_xxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzz_0[i] * pb_z + g_0_xxxx_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyzzz_0[i] = 3.0 * g_0_xxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzz_0[i] * pb_z + g_0_xxxx_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxzzzz_0[i] = 4.0 * g_0_xxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzzzz_0[i] * pb_z + g_0_xxxx_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyy_0[i] = g_0_xxxx_0_xyyyyy_0[i] * pb_z + g_0_xxxx_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyz_0[i] = g_0_xxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyz_0[i] * pb_z + g_0_xxxx_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyzz_0[i] = 2.0 * g_0_xxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyzz_0[i] * pb_z + g_0_xxxx_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyzzz_0[i] = 3.0 * g_0_xxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzzz_0[i] * pb_z + g_0_xxxx_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyzzzz_0[i] = 4.0 * g_0_xxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzzz_0[i] * pb_z + g_0_xxxx_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xzzzzz_0[i] = 5.0 * g_0_xxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzzzzz_0[i] * pb_z + g_0_xxxx_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyy_0[i] = g_0_xxxx_0_yyyyyy_0[i] * pb_z + g_0_xxxx_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyz_0[i] = g_0_xxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyz_0[i] * pb_z + g_0_xxxx_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyzz_0[i] = 2.0 * g_0_xxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyzz_0[i] * pb_z + g_0_xxxx_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyzzz_0[i] = 3.0 * g_0_xxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyzzz_0[i] * pb_z + g_0_xxxx_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyzzzz_0[i] = 4.0 * g_0_xxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzzzz_0[i] * pb_z + g_0_xxxx_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yzzzzz_0[i] = 5.0 * g_0_xxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzzzz_0[i] * pb_z + g_0_xxxx_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_zzzzzz_0[i] = 6.0 * g_0_xxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_zzzzzz_0[i] * pb_z + g_0_xxxx_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 84-112 components of targeted buffer : SHSI

    auto g_0_xxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 84);

    auto g_0_xxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 85);

    auto g_0_xxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 86);

    auto g_0_xxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 87);

    auto g_0_xxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 88);

    auto g_0_xxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 89);

    auto g_0_xxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 90);

    auto g_0_xxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 91);

    auto g_0_xxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 92);

    auto g_0_xxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 93);

    auto g_0_xxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 94);

    auto g_0_xxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 95);

    auto g_0_xxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 96);

    auto g_0_xxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 97);

    auto g_0_xxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 98);

    auto g_0_xxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 99);

    auto g_0_xxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 100);

    auto g_0_xxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 101);

    auto g_0_xxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 102);

    auto g_0_xxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 103);

    auto g_0_xxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 104);

    auto g_0_xxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 105);

    auto g_0_xxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 106);

    auto g_0_xxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 107);

    auto g_0_xxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 108);

    auto g_0_xxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 109);

    auto g_0_xxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 110);

    auto g_0_xxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 111);

#pragma omp simd aligned(g_0_xxx_0_xxxxxx_0,       \
                             g_0_xxx_0_xxxxxx_1,   \
                             g_0_xxx_0_xxxxxz_0,   \
                             g_0_xxx_0_xxxxxz_1,   \
                             g_0_xxx_0_xxxxzz_0,   \
                             g_0_xxx_0_xxxxzz_1,   \
                             g_0_xxx_0_xxxzzz_0,   \
                             g_0_xxx_0_xxxzzz_1,   \
                             g_0_xxx_0_xxzzzz_0,   \
                             g_0_xxx_0_xxzzzz_1,   \
                             g_0_xxx_0_xzzzzz_0,   \
                             g_0_xxx_0_xzzzzz_1,   \
                             g_0_xxxy_0_xxxxxx_0,  \
                             g_0_xxxy_0_xxxxxx_1,  \
                             g_0_xxxy_0_xxxxxz_0,  \
                             g_0_xxxy_0_xxxxxz_1,  \
                             g_0_xxxy_0_xxxxzz_0,  \
                             g_0_xxxy_0_xxxxzz_1,  \
                             g_0_xxxy_0_xxxzzz_0,  \
                             g_0_xxxy_0_xxxzzz_1,  \
                             g_0_xxxy_0_xxzzzz_0,  \
                             g_0_xxxy_0_xxzzzz_1,  \
                             g_0_xxxy_0_xzzzzz_0,  \
                             g_0_xxxy_0_xzzzzz_1,  \
                             g_0_xxxyy_0_xxxxxx_0, \
                             g_0_xxxyy_0_xxxxxy_0, \
                             g_0_xxxyy_0_xxxxxz_0, \
                             g_0_xxxyy_0_xxxxyy_0, \
                             g_0_xxxyy_0_xxxxyz_0, \
                             g_0_xxxyy_0_xxxxzz_0, \
                             g_0_xxxyy_0_xxxyyy_0, \
                             g_0_xxxyy_0_xxxyyz_0, \
                             g_0_xxxyy_0_xxxyzz_0, \
                             g_0_xxxyy_0_xxxzzz_0, \
                             g_0_xxxyy_0_xxyyyy_0, \
                             g_0_xxxyy_0_xxyyyz_0, \
                             g_0_xxxyy_0_xxyyzz_0, \
                             g_0_xxxyy_0_xxyzzz_0, \
                             g_0_xxxyy_0_xxzzzz_0, \
                             g_0_xxxyy_0_xyyyyy_0, \
                             g_0_xxxyy_0_xyyyyz_0, \
                             g_0_xxxyy_0_xyyyzz_0, \
                             g_0_xxxyy_0_xyyzzz_0, \
                             g_0_xxxyy_0_xyzzzz_0, \
                             g_0_xxxyy_0_xzzzzz_0, \
                             g_0_xxxyy_0_yyyyyy_0, \
                             g_0_xxxyy_0_yyyyyz_0, \
                             g_0_xxxyy_0_yyyyzz_0, \
                             g_0_xxxyy_0_yyyzzz_0, \
                             g_0_xxxyy_0_yyzzzz_0, \
                             g_0_xxxyy_0_yzzzzz_0, \
                             g_0_xxxyy_0_zzzzzz_0, \
                             g_0_xxyy_0_xxxxxy_0,  \
                             g_0_xxyy_0_xxxxxy_1,  \
                             g_0_xxyy_0_xxxxy_1,   \
                             g_0_xxyy_0_xxxxyy_0,  \
                             g_0_xxyy_0_xxxxyy_1,  \
                             g_0_xxyy_0_xxxxyz_0,  \
                             g_0_xxyy_0_xxxxyz_1,  \
                             g_0_xxyy_0_xxxyy_1,   \
                             g_0_xxyy_0_xxxyyy_0,  \
                             g_0_xxyy_0_xxxyyy_1,  \
                             g_0_xxyy_0_xxxyyz_0,  \
                             g_0_xxyy_0_xxxyyz_1,  \
                             g_0_xxyy_0_xxxyz_1,   \
                             g_0_xxyy_0_xxxyzz_0,  \
                             g_0_xxyy_0_xxxyzz_1,  \
                             g_0_xxyy_0_xxyyy_1,   \
                             g_0_xxyy_0_xxyyyy_0,  \
                             g_0_xxyy_0_xxyyyy_1,  \
                             g_0_xxyy_0_xxyyyz_0,  \
                             g_0_xxyy_0_xxyyyz_1,  \
                             g_0_xxyy_0_xxyyz_1,   \
                             g_0_xxyy_0_xxyyzz_0,  \
                             g_0_xxyy_0_xxyyzz_1,  \
                             g_0_xxyy_0_xxyzz_1,   \
                             g_0_xxyy_0_xxyzzz_0,  \
                             g_0_xxyy_0_xxyzzz_1,  \
                             g_0_xxyy_0_xyyyy_1,   \
                             g_0_xxyy_0_xyyyyy_0,  \
                             g_0_xxyy_0_xyyyyy_1,  \
                             g_0_xxyy_0_xyyyyz_0,  \
                             g_0_xxyy_0_xyyyyz_1,  \
                             g_0_xxyy_0_xyyyz_1,   \
                             g_0_xxyy_0_xyyyzz_0,  \
                             g_0_xxyy_0_xyyyzz_1,  \
                             g_0_xxyy_0_xyyzz_1,   \
                             g_0_xxyy_0_xyyzzz_0,  \
                             g_0_xxyy_0_xyyzzz_1,  \
                             g_0_xxyy_0_xyzzz_1,   \
                             g_0_xxyy_0_xyzzzz_0,  \
                             g_0_xxyy_0_xyzzzz_1,  \
                             g_0_xxyy_0_yyyyy_1,   \
                             g_0_xxyy_0_yyyyyy_0,  \
                             g_0_xxyy_0_yyyyyy_1,  \
                             g_0_xxyy_0_yyyyyz_0,  \
                             g_0_xxyy_0_yyyyyz_1,  \
                             g_0_xxyy_0_yyyyz_1,   \
                             g_0_xxyy_0_yyyyzz_0,  \
                             g_0_xxyy_0_yyyyzz_1,  \
                             g_0_xxyy_0_yyyzz_1,   \
                             g_0_xxyy_0_yyyzzz_0,  \
                             g_0_xxyy_0_yyyzzz_1,  \
                             g_0_xxyy_0_yyzzz_1,   \
                             g_0_xxyy_0_yyzzzz_0,  \
                             g_0_xxyy_0_yyzzzz_1,  \
                             g_0_xxyy_0_yzzzz_1,   \
                             g_0_xxyy_0_yzzzzz_0,  \
                             g_0_xxyy_0_yzzzzz_1,  \
                             g_0_xxyy_0_zzzzzz_0,  \
                             g_0_xxyy_0_zzzzzz_1,  \
                             g_0_xyy_0_xxxxxy_0,   \
                             g_0_xyy_0_xxxxxy_1,   \
                             g_0_xyy_0_xxxxyy_0,   \
                             g_0_xyy_0_xxxxyy_1,   \
                             g_0_xyy_0_xxxxyz_0,   \
                             g_0_xyy_0_xxxxyz_1,   \
                             g_0_xyy_0_xxxyyy_0,   \
                             g_0_xyy_0_xxxyyy_1,   \
                             g_0_xyy_0_xxxyyz_0,   \
                             g_0_xyy_0_xxxyyz_1,   \
                             g_0_xyy_0_xxxyzz_0,   \
                             g_0_xyy_0_xxxyzz_1,   \
                             g_0_xyy_0_xxyyyy_0,   \
                             g_0_xyy_0_xxyyyy_1,   \
                             g_0_xyy_0_xxyyyz_0,   \
                             g_0_xyy_0_xxyyyz_1,   \
                             g_0_xyy_0_xxyyzz_0,   \
                             g_0_xyy_0_xxyyzz_1,   \
                             g_0_xyy_0_xxyzzz_0,   \
                             g_0_xyy_0_xxyzzz_1,   \
                             g_0_xyy_0_xyyyyy_0,   \
                             g_0_xyy_0_xyyyyy_1,   \
                             g_0_xyy_0_xyyyyz_0,   \
                             g_0_xyy_0_xyyyyz_1,   \
                             g_0_xyy_0_xyyyzz_0,   \
                             g_0_xyy_0_xyyyzz_1,   \
                             g_0_xyy_0_xyyzzz_0,   \
                             g_0_xyy_0_xyyzzz_1,   \
                             g_0_xyy_0_xyzzzz_0,   \
                             g_0_xyy_0_xyzzzz_1,   \
                             g_0_xyy_0_yyyyyy_0,   \
                             g_0_xyy_0_yyyyyy_1,   \
                             g_0_xyy_0_yyyyyz_0,   \
                             g_0_xyy_0_yyyyyz_1,   \
                             g_0_xyy_0_yyyyzz_0,   \
                             g_0_xyy_0_yyyyzz_1,   \
                             g_0_xyy_0_yyyzzz_0,   \
                             g_0_xyy_0_yyyzzz_1,   \
                             g_0_xyy_0_yyzzzz_0,   \
                             g_0_xyy_0_yyzzzz_1,   \
                             g_0_xyy_0_yzzzzz_0,   \
                             g_0_xyy_0_yzzzzz_1,   \
                             g_0_xyy_0_zzzzzz_0,   \
                             g_0_xyy_0_zzzzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyy_0_xxxxxx_0[i] =
            g_0_xxx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxx_0[i] * pb_y + g_0_xxxy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxxy_0[i] = 2.0 * g_0_xyy_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxy_1[i] * fti_ab_0 +
                                  5.0 * g_0_xxyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxy_0[i] * pb_x + g_0_xxyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxz_0[i] =
            g_0_xxx_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxz_0[i] * pb_y + g_0_xxxy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxyy_0[i] = 2.0 * g_0_xyy_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyy_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyy_0[i] * pb_x + g_0_xxyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxyz_0[i] = 2.0 * g_0_xyy_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyz_0[i] * pb_x + g_0_xxyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxzz_0[i] =
            g_0_xxx_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxzz_0[i] * pb_y + g_0_xxxy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxyyy_0[i] = 2.0 * g_0_xyy_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyy_0[i] * pb_x + g_0_xxyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyyz_0[i] = 2.0 * g_0_xyy_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyz_0[i] * pb_x + g_0_xxyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyzz_0[i] = 2.0 * g_0_xyy_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyzz_0[i] * pb_x + g_0_xxyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxzzz_0[i] =
            g_0_xxx_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxzzz_0[i] * pb_y + g_0_xxxy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxyyyy_0[i] = 2.0 * g_0_xyy_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyy_0[i] * pb_x + g_0_xxyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyyz_0[i] = 2.0 * g_0_xyy_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyz_0[i] * pb_x + g_0_xxyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyzz_0[i] = 2.0 * g_0_xyy_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyzz_0[i] * pb_x + g_0_xxyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyzzz_0[i] = 2.0 * g_0_xyy_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyzzz_0[i] * pb_x + g_0_xxyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxzzzz_0[i] =
            g_0_xxx_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxzzzz_0[i] * pb_y + g_0_xxxy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xyyyyy_0[i] = 2.0 * g_0_xyy_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyy_1[i] * fi_abcd_0 +
                                  g_0_xxyy_0_xyyyyy_0[i] * pb_x + g_0_xxyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyyz_0[i] = 2.0 * g_0_xyy_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyz_1[i] * fi_abcd_0 +
                                  g_0_xxyy_0_xyyyyz_0[i] * pb_x + g_0_xxyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyzz_0[i] = 2.0 * g_0_xyy_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyzz_1[i] * fi_abcd_0 +
                                  g_0_xxyy_0_xyyyzz_0[i] * pb_x + g_0_xxyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyzzz_0[i] = 2.0 * g_0_xyy_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzzz_1[i] * fi_abcd_0 +
                                  g_0_xxyy_0_xyyzzz_0[i] * pb_x + g_0_xxyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyzzzz_0[i] = 2.0 * g_0_xyy_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzzz_1[i] * fi_abcd_0 +
                                  g_0_xxyy_0_xyzzzz_0[i] * pb_x + g_0_xxyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xzzzzz_0[i] =
            g_0_xxx_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xzzzzz_0[i] * pb_y + g_0_xxxy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_yyyyyy_0[i] = 2.0 * g_0_xyy_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyy_0[i] * pb_x +
                                  g_0_xxyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyyz_0[i] = 2.0 * g_0_xyy_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyz_0[i] * pb_x +
                                  g_0_xxyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyzz_0[i] = 2.0 * g_0_xyy_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyzz_0[i] * pb_x +
                                  g_0_xxyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyzzz_0[i] = 2.0 * g_0_xyy_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyzzz_0[i] * pb_x +
                                  g_0_xxyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyzzzz_0[i] = 2.0 * g_0_xyy_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzzzz_0[i] * pb_x +
                                  g_0_xxyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yzzzzz_0[i] = 2.0 * g_0_xyy_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzzzz_0[i] * pb_x +
                                  g_0_xxyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_zzzzzz_0[i] = 2.0 * g_0_xyy_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_zzzzzz_0[i] * pb_x +
                                  g_0_xxyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 112-140 components of targeted buffer : SHSI

    auto g_0_xxxyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 112);

    auto g_0_xxxyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 113);

    auto g_0_xxxyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 114);

    auto g_0_xxxyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 115);

    auto g_0_xxxyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 116);

    auto g_0_xxxyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 117);

    auto g_0_xxxyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 118);

    auto g_0_xxxyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 119);

    auto g_0_xxxyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 120);

    auto g_0_xxxyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 121);

    auto g_0_xxxyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 122);

    auto g_0_xxxyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 123);

    auto g_0_xxxyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 124);

    auto g_0_xxxyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 125);

    auto g_0_xxxyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 126);

    auto g_0_xxxyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 127);

    auto g_0_xxxyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 128);

    auto g_0_xxxyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 129);

    auto g_0_xxxyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 130);

    auto g_0_xxxyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 131);

    auto g_0_xxxyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 132);

    auto g_0_xxxyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 133);

    auto g_0_xxxyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 134);

    auto g_0_xxxyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 135);

    auto g_0_xxxyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 136);

    auto g_0_xxxyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 137);

    auto g_0_xxxyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 138);

    auto g_0_xxxyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 139);

#pragma omp simd aligned(g_0_xxxy_0_xxxxxy_0,      \
                             g_0_xxxy_0_xxxxxy_1,  \
                             g_0_xxxy_0_xxxxyy_0,  \
                             g_0_xxxy_0_xxxxyy_1,  \
                             g_0_xxxy_0_xxxyyy_0,  \
                             g_0_xxxy_0_xxxyyy_1,  \
                             g_0_xxxy_0_xxyyyy_0,  \
                             g_0_xxxy_0_xxyyyy_1,  \
                             g_0_xxxy_0_xyyyyy_0,  \
                             g_0_xxxy_0_xyyyyy_1,  \
                             g_0_xxxy_0_yyyyyy_0,  \
                             g_0_xxxy_0_yyyyyy_1,  \
                             g_0_xxxyz_0_xxxxxx_0, \
                             g_0_xxxyz_0_xxxxxy_0, \
                             g_0_xxxyz_0_xxxxxz_0, \
                             g_0_xxxyz_0_xxxxyy_0, \
                             g_0_xxxyz_0_xxxxyz_0, \
                             g_0_xxxyz_0_xxxxzz_0, \
                             g_0_xxxyz_0_xxxyyy_0, \
                             g_0_xxxyz_0_xxxyyz_0, \
                             g_0_xxxyz_0_xxxyzz_0, \
                             g_0_xxxyz_0_xxxzzz_0, \
                             g_0_xxxyz_0_xxyyyy_0, \
                             g_0_xxxyz_0_xxyyyz_0, \
                             g_0_xxxyz_0_xxyyzz_0, \
                             g_0_xxxyz_0_xxyzzz_0, \
                             g_0_xxxyz_0_xxzzzz_0, \
                             g_0_xxxyz_0_xyyyyy_0, \
                             g_0_xxxyz_0_xyyyyz_0, \
                             g_0_xxxyz_0_xyyyzz_0, \
                             g_0_xxxyz_0_xyyzzz_0, \
                             g_0_xxxyz_0_xyzzzz_0, \
                             g_0_xxxyz_0_xzzzzz_0, \
                             g_0_xxxyz_0_yyyyyy_0, \
                             g_0_xxxyz_0_yyyyyz_0, \
                             g_0_xxxyz_0_yyyyzz_0, \
                             g_0_xxxyz_0_yyyzzz_0, \
                             g_0_xxxyz_0_yyzzzz_0, \
                             g_0_xxxyz_0_yzzzzz_0, \
                             g_0_xxxyz_0_zzzzzz_0, \
                             g_0_xxxz_0_xxxxxx_0,  \
                             g_0_xxxz_0_xxxxxx_1,  \
                             g_0_xxxz_0_xxxxxz_0,  \
                             g_0_xxxz_0_xxxxxz_1,  \
                             g_0_xxxz_0_xxxxyz_0,  \
                             g_0_xxxz_0_xxxxyz_1,  \
                             g_0_xxxz_0_xxxxz_1,   \
                             g_0_xxxz_0_xxxxzz_0,  \
                             g_0_xxxz_0_xxxxzz_1,  \
                             g_0_xxxz_0_xxxyyz_0,  \
                             g_0_xxxz_0_xxxyyz_1,  \
                             g_0_xxxz_0_xxxyz_1,   \
                             g_0_xxxz_0_xxxyzz_0,  \
                             g_0_xxxz_0_xxxyzz_1,  \
                             g_0_xxxz_0_xxxzz_1,   \
                             g_0_xxxz_0_xxxzzz_0,  \
                             g_0_xxxz_0_xxxzzz_1,  \
                             g_0_xxxz_0_xxyyyz_0,  \
                             g_0_xxxz_0_xxyyyz_1,  \
                             g_0_xxxz_0_xxyyz_1,   \
                             g_0_xxxz_0_xxyyzz_0,  \
                             g_0_xxxz_0_xxyyzz_1,  \
                             g_0_xxxz_0_xxyzz_1,   \
                             g_0_xxxz_0_xxyzzz_0,  \
                             g_0_xxxz_0_xxyzzz_1,  \
                             g_0_xxxz_0_xxzzz_1,   \
                             g_0_xxxz_0_xxzzzz_0,  \
                             g_0_xxxz_0_xxzzzz_1,  \
                             g_0_xxxz_0_xyyyyz_0,  \
                             g_0_xxxz_0_xyyyyz_1,  \
                             g_0_xxxz_0_xyyyz_1,   \
                             g_0_xxxz_0_xyyyzz_0,  \
                             g_0_xxxz_0_xyyyzz_1,  \
                             g_0_xxxz_0_xyyzz_1,   \
                             g_0_xxxz_0_xyyzzz_0,  \
                             g_0_xxxz_0_xyyzzz_1,  \
                             g_0_xxxz_0_xyzzz_1,   \
                             g_0_xxxz_0_xyzzzz_0,  \
                             g_0_xxxz_0_xyzzzz_1,  \
                             g_0_xxxz_0_xzzzz_1,   \
                             g_0_xxxz_0_xzzzzz_0,  \
                             g_0_xxxz_0_xzzzzz_1,  \
                             g_0_xxxz_0_yyyyyz_0,  \
                             g_0_xxxz_0_yyyyyz_1,  \
                             g_0_xxxz_0_yyyyz_1,   \
                             g_0_xxxz_0_yyyyzz_0,  \
                             g_0_xxxz_0_yyyyzz_1,  \
                             g_0_xxxz_0_yyyzz_1,   \
                             g_0_xxxz_0_yyyzzz_0,  \
                             g_0_xxxz_0_yyyzzz_1,  \
                             g_0_xxxz_0_yyzzz_1,   \
                             g_0_xxxz_0_yyzzzz_0,  \
                             g_0_xxxz_0_yyzzzz_1,  \
                             g_0_xxxz_0_yzzzz_1,   \
                             g_0_xxxz_0_yzzzzz_0,  \
                             g_0_xxxz_0_yzzzzz_1,  \
                             g_0_xxxz_0_zzzzz_1,   \
                             g_0_xxxz_0_zzzzzz_0,  \
                             g_0_xxxz_0_zzzzzz_1,  \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyz_0_xxxxxx_0[i] = g_0_xxxz_0_xxxxxx_0[i] * pb_y + g_0_xxxz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxy_0[i] = g_0_xxxy_0_xxxxxy_0[i] * pb_z + g_0_xxxy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxxz_0[i] = g_0_xxxz_0_xxxxxz_0[i] * pb_y + g_0_xxxz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxyy_0[i] = g_0_xxxy_0_xxxxyy_0[i] * pb_z + g_0_xxxy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxyz_0[i] = g_0_xxxz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxyz_0[i] * pb_y + g_0_xxxz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxzz_0[i] = g_0_xxxz_0_xxxxzz_0[i] * pb_y + g_0_xxxz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyyy_0[i] = g_0_xxxy_0_xxxyyy_0[i] * pb_z + g_0_xxxy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxyyz_0[i] = 2.0 * g_0_xxxz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyyz_0[i] * pb_y + g_0_xxxz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyzz_0[i] = g_0_xxxz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyzz_0[i] * pb_y + g_0_xxxz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxzzz_0[i] = g_0_xxxz_0_xxxzzz_0[i] * pb_y + g_0_xxxz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyyy_0[i] = g_0_xxxy_0_xxyyyy_0[i] * pb_z + g_0_xxxy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxyyyz_0[i] = 3.0 * g_0_xxxz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyyz_0[i] * pb_y + g_0_xxxz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyzz_0[i] = 2.0 * g_0_xxxz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyzz_0[i] * pb_y + g_0_xxxz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyzzz_0[i] = g_0_xxxz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyzzz_0[i] * pb_y + g_0_xxxz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxzzzz_0[i] = g_0_xxxz_0_xxzzzz_0[i] * pb_y + g_0_xxxz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyyy_0[i] = g_0_xxxy_0_xyyyyy_0[i] * pb_z + g_0_xxxy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xyyyyz_0[i] = 4.0 * g_0_xxxz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyyz_0[i] * pb_y + g_0_xxxz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyzz_0[i] = 3.0 * g_0_xxxz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyzz_0[i] * pb_y + g_0_xxxz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyzzz_0[i] = 2.0 * g_0_xxxz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyzzz_0[i] * pb_y + g_0_xxxz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyzzzz_0[i] = g_0_xxxz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyzzzz_0[i] * pb_y + g_0_xxxz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xzzzzz_0[i] = g_0_xxxz_0_xzzzzz_0[i] * pb_y + g_0_xxxz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyyy_0[i] = g_0_xxxy_0_yyyyyy_0[i] * pb_z + g_0_xxxy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_yyyyyz_0[i] = 5.0 * g_0_xxxz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyyz_0[i] * pb_y + g_0_xxxz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyzz_0[i] = 4.0 * g_0_xxxz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyzz_0[i] * pb_y + g_0_xxxz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyzzz_0[i] = 3.0 * g_0_xxxz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyzzz_0[i] * pb_y + g_0_xxxz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyzzzz_0[i] = 2.0 * g_0_xxxz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyzzzz_0[i] * pb_y + g_0_xxxz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yzzzzz_0[i] = g_0_xxxz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yzzzzz_0[i] * pb_y + g_0_xxxz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_zzzzzz_0[i] = g_0_xxxz_0_zzzzzz_0[i] * pb_y + g_0_xxxz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 140-168 components of targeted buffer : SHSI

    auto g_0_xxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 140);

    auto g_0_xxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 141);

    auto g_0_xxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 142);

    auto g_0_xxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 143);

    auto g_0_xxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 144);

    auto g_0_xxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 145);

    auto g_0_xxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 146);

    auto g_0_xxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 147);

    auto g_0_xxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 148);

    auto g_0_xxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 149);

    auto g_0_xxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 150);

    auto g_0_xxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 151);

    auto g_0_xxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 152);

    auto g_0_xxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 153);

    auto g_0_xxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 154);

    auto g_0_xxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 155);

    auto g_0_xxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 156);

    auto g_0_xxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 157);

    auto g_0_xxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 158);

    auto g_0_xxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 159);

    auto g_0_xxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 160);

    auto g_0_xxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 161);

    auto g_0_xxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 162);

    auto g_0_xxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 163);

    auto g_0_xxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 164);

    auto g_0_xxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 165);

    auto g_0_xxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 166);

    auto g_0_xxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 167);

#pragma omp simd aligned(g_0_xxx_0_xxxxxx_0,       \
                             g_0_xxx_0_xxxxxx_1,   \
                             g_0_xxx_0_xxxxxy_0,   \
                             g_0_xxx_0_xxxxxy_1,   \
                             g_0_xxx_0_xxxxyy_0,   \
                             g_0_xxx_0_xxxxyy_1,   \
                             g_0_xxx_0_xxxyyy_0,   \
                             g_0_xxx_0_xxxyyy_1,   \
                             g_0_xxx_0_xxyyyy_0,   \
                             g_0_xxx_0_xxyyyy_1,   \
                             g_0_xxx_0_xyyyyy_0,   \
                             g_0_xxx_0_xyyyyy_1,   \
                             g_0_xxxz_0_xxxxxx_0,  \
                             g_0_xxxz_0_xxxxxx_1,  \
                             g_0_xxxz_0_xxxxxy_0,  \
                             g_0_xxxz_0_xxxxxy_1,  \
                             g_0_xxxz_0_xxxxyy_0,  \
                             g_0_xxxz_0_xxxxyy_1,  \
                             g_0_xxxz_0_xxxyyy_0,  \
                             g_0_xxxz_0_xxxyyy_1,  \
                             g_0_xxxz_0_xxyyyy_0,  \
                             g_0_xxxz_0_xxyyyy_1,  \
                             g_0_xxxz_0_xyyyyy_0,  \
                             g_0_xxxz_0_xyyyyy_1,  \
                             g_0_xxxzz_0_xxxxxx_0, \
                             g_0_xxxzz_0_xxxxxy_0, \
                             g_0_xxxzz_0_xxxxxz_0, \
                             g_0_xxxzz_0_xxxxyy_0, \
                             g_0_xxxzz_0_xxxxyz_0, \
                             g_0_xxxzz_0_xxxxzz_0, \
                             g_0_xxxzz_0_xxxyyy_0, \
                             g_0_xxxzz_0_xxxyyz_0, \
                             g_0_xxxzz_0_xxxyzz_0, \
                             g_0_xxxzz_0_xxxzzz_0, \
                             g_0_xxxzz_0_xxyyyy_0, \
                             g_0_xxxzz_0_xxyyyz_0, \
                             g_0_xxxzz_0_xxyyzz_0, \
                             g_0_xxxzz_0_xxyzzz_0, \
                             g_0_xxxzz_0_xxzzzz_0, \
                             g_0_xxxzz_0_xyyyyy_0, \
                             g_0_xxxzz_0_xyyyyz_0, \
                             g_0_xxxzz_0_xyyyzz_0, \
                             g_0_xxxzz_0_xyyzzz_0, \
                             g_0_xxxzz_0_xyzzzz_0, \
                             g_0_xxxzz_0_xzzzzz_0, \
                             g_0_xxxzz_0_yyyyyy_0, \
                             g_0_xxxzz_0_yyyyyz_0, \
                             g_0_xxxzz_0_yyyyzz_0, \
                             g_0_xxxzz_0_yyyzzz_0, \
                             g_0_xxxzz_0_yyzzzz_0, \
                             g_0_xxxzz_0_yzzzzz_0, \
                             g_0_xxxzz_0_zzzzzz_0, \
                             g_0_xxzz_0_xxxxxz_0,  \
                             g_0_xxzz_0_xxxxxz_1,  \
                             g_0_xxzz_0_xxxxyz_0,  \
                             g_0_xxzz_0_xxxxyz_1,  \
                             g_0_xxzz_0_xxxxz_1,   \
                             g_0_xxzz_0_xxxxzz_0,  \
                             g_0_xxzz_0_xxxxzz_1,  \
                             g_0_xxzz_0_xxxyyz_0,  \
                             g_0_xxzz_0_xxxyyz_1,  \
                             g_0_xxzz_0_xxxyz_1,   \
                             g_0_xxzz_0_xxxyzz_0,  \
                             g_0_xxzz_0_xxxyzz_1,  \
                             g_0_xxzz_0_xxxzz_1,   \
                             g_0_xxzz_0_xxxzzz_0,  \
                             g_0_xxzz_0_xxxzzz_1,  \
                             g_0_xxzz_0_xxyyyz_0,  \
                             g_0_xxzz_0_xxyyyz_1,  \
                             g_0_xxzz_0_xxyyz_1,   \
                             g_0_xxzz_0_xxyyzz_0,  \
                             g_0_xxzz_0_xxyyzz_1,  \
                             g_0_xxzz_0_xxyzz_1,   \
                             g_0_xxzz_0_xxyzzz_0,  \
                             g_0_xxzz_0_xxyzzz_1,  \
                             g_0_xxzz_0_xxzzz_1,   \
                             g_0_xxzz_0_xxzzzz_0,  \
                             g_0_xxzz_0_xxzzzz_1,  \
                             g_0_xxzz_0_xyyyyz_0,  \
                             g_0_xxzz_0_xyyyyz_1,  \
                             g_0_xxzz_0_xyyyz_1,   \
                             g_0_xxzz_0_xyyyzz_0,  \
                             g_0_xxzz_0_xyyyzz_1,  \
                             g_0_xxzz_0_xyyzz_1,   \
                             g_0_xxzz_0_xyyzzz_0,  \
                             g_0_xxzz_0_xyyzzz_1,  \
                             g_0_xxzz_0_xyzzz_1,   \
                             g_0_xxzz_0_xyzzzz_0,  \
                             g_0_xxzz_0_xyzzzz_1,  \
                             g_0_xxzz_0_xzzzz_1,   \
                             g_0_xxzz_0_xzzzzz_0,  \
                             g_0_xxzz_0_xzzzzz_1,  \
                             g_0_xxzz_0_yyyyyy_0,  \
                             g_0_xxzz_0_yyyyyy_1,  \
                             g_0_xxzz_0_yyyyyz_0,  \
                             g_0_xxzz_0_yyyyyz_1,  \
                             g_0_xxzz_0_yyyyz_1,   \
                             g_0_xxzz_0_yyyyzz_0,  \
                             g_0_xxzz_0_yyyyzz_1,  \
                             g_0_xxzz_0_yyyzz_1,   \
                             g_0_xxzz_0_yyyzzz_0,  \
                             g_0_xxzz_0_yyyzzz_1,  \
                             g_0_xxzz_0_yyzzz_1,   \
                             g_0_xxzz_0_yyzzzz_0,  \
                             g_0_xxzz_0_yyzzzz_1,  \
                             g_0_xxzz_0_yzzzz_1,   \
                             g_0_xxzz_0_yzzzzz_0,  \
                             g_0_xxzz_0_yzzzzz_1,  \
                             g_0_xxzz_0_zzzzz_1,   \
                             g_0_xxzz_0_zzzzzz_0,  \
                             g_0_xxzz_0_zzzzzz_1,  \
                             g_0_xzz_0_xxxxxz_0,   \
                             g_0_xzz_0_xxxxxz_1,   \
                             g_0_xzz_0_xxxxyz_0,   \
                             g_0_xzz_0_xxxxyz_1,   \
                             g_0_xzz_0_xxxxzz_0,   \
                             g_0_xzz_0_xxxxzz_1,   \
                             g_0_xzz_0_xxxyyz_0,   \
                             g_0_xzz_0_xxxyyz_1,   \
                             g_0_xzz_0_xxxyzz_0,   \
                             g_0_xzz_0_xxxyzz_1,   \
                             g_0_xzz_0_xxxzzz_0,   \
                             g_0_xzz_0_xxxzzz_1,   \
                             g_0_xzz_0_xxyyyz_0,   \
                             g_0_xzz_0_xxyyyz_1,   \
                             g_0_xzz_0_xxyyzz_0,   \
                             g_0_xzz_0_xxyyzz_1,   \
                             g_0_xzz_0_xxyzzz_0,   \
                             g_0_xzz_0_xxyzzz_1,   \
                             g_0_xzz_0_xxzzzz_0,   \
                             g_0_xzz_0_xxzzzz_1,   \
                             g_0_xzz_0_xyyyyz_0,   \
                             g_0_xzz_0_xyyyyz_1,   \
                             g_0_xzz_0_xyyyzz_0,   \
                             g_0_xzz_0_xyyyzz_1,   \
                             g_0_xzz_0_xyyzzz_0,   \
                             g_0_xzz_0_xyyzzz_1,   \
                             g_0_xzz_0_xyzzzz_0,   \
                             g_0_xzz_0_xyzzzz_1,   \
                             g_0_xzz_0_xzzzzz_0,   \
                             g_0_xzz_0_xzzzzz_1,   \
                             g_0_xzz_0_yyyyyy_0,   \
                             g_0_xzz_0_yyyyyy_1,   \
                             g_0_xzz_0_yyyyyz_0,   \
                             g_0_xzz_0_yyyyyz_1,   \
                             g_0_xzz_0_yyyyzz_0,   \
                             g_0_xzz_0_yyyyzz_1,   \
                             g_0_xzz_0_yyyzzz_0,   \
                             g_0_xzz_0_yyyzzz_1,   \
                             g_0_xzz_0_yyzzzz_0,   \
                             g_0_xzz_0_yyzzzz_1,   \
                             g_0_xzz_0_yzzzzz_0,   \
                             g_0_xzz_0_yzzzzz_1,   \
                             g_0_xzz_0_zzzzzz_0,   \
                             g_0_xzz_0_zzzzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzz_0_xxxxxx_0[i] =
            g_0_xxx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxx_0[i] * pb_z + g_0_xxxz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxy_0[i] =
            g_0_xxx_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxy_0[i] * pb_z + g_0_xxxz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxz_0[i] = 2.0 * g_0_xzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxz_1[i] * fti_ab_0 +
                                  5.0 * g_0_xxzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxz_0[i] * pb_x + g_0_xxzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxyy_0[i] =
            g_0_xxx_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxyy_0[i] * pb_z + g_0_xxxz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxyz_0[i] = 2.0 * g_0_xzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyz_0[i] * pb_x + g_0_xxzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxzz_0[i] = 2.0 * g_0_xzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_xxzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxzz_0[i] * pb_x + g_0_xxzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyyy_0[i] =
            g_0_xxx_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxyyy_0[i] * pb_z + g_0_xxxz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxyyz_0[i] = 2.0 * g_0_xzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyz_0[i] * pb_x + g_0_xxzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyzz_0[i] = 2.0 * g_0_xzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyzz_0[i] * pb_x + g_0_xxzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxzzz_0[i] = 2.0 * g_0_xzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_xxzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxzzz_0[i] * pb_x + g_0_xxzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyyy_0[i] =
            g_0_xxx_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxyyyy_0[i] * pb_z + g_0_xxxz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxyyyz_0[i] = 2.0 * g_0_xzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyz_0[i] * pb_x + g_0_xxzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyzz_0[i] = 2.0 * g_0_xzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyzz_0[i] * pb_x + g_0_xxzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyzzz_0[i] = 2.0 * g_0_xzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyzzz_0[i] * pb_x + g_0_xxzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxzzzz_0[i] = 2.0 * g_0_xzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_xxzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxzzzz_0[i] * pb_x + g_0_xxzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyyy_0[i] =
            g_0_xxx_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xyyyyy_0[i] * pb_z + g_0_xxxz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xyyyyz_0[i] = 2.0 * g_0_xzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyz_1[i] * fi_abcd_0 +
                                  g_0_xxzz_0_xyyyyz_0[i] * pb_x + g_0_xxzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyzz_0[i] = 2.0 * g_0_xzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyzz_1[i] * fi_abcd_0 +
                                  g_0_xxzz_0_xyyyzz_0[i] * pb_x + g_0_xxzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyzzz_0[i] = 2.0 * g_0_xzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzzz_1[i] * fi_abcd_0 +
                                  g_0_xxzz_0_xyyzzz_0[i] * pb_x + g_0_xxzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyzzzz_0[i] = 2.0 * g_0_xzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzzz_1[i] * fi_abcd_0 +
                                  g_0_xxzz_0_xyzzzz_0[i] * pb_x + g_0_xxzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xzzzzz_0[i] = 2.0 * g_0_xzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzzz_1[i] * fi_abcd_0 +
                                  g_0_xxzz_0_xzzzzz_0[i] * pb_x + g_0_xxzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyy_0[i] = 2.0 * g_0_xzz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyy_0[i] * pb_x +
                                  g_0_xxzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyz_0[i] = 2.0 * g_0_xzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyz_0[i] * pb_x +
                                  g_0_xxzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyzz_0[i] = 2.0 * g_0_xzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyzz_0[i] * pb_x +
                                  g_0_xxzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyzzz_0[i] = 2.0 * g_0_xzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyzzz_0[i] * pb_x +
                                  g_0_xxzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyzzzz_0[i] = 2.0 * g_0_xzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzzzz_0[i] * pb_x +
                                  g_0_xxzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yzzzzz_0[i] = 2.0 * g_0_xzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzzzz_0[i] * pb_x +
                                  g_0_xxzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_zzzzzz_0[i] = 2.0 * g_0_xzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzzzz_0[i] * pb_x +
                                  g_0_xxzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 168-196 components of targeted buffer : SHSI

    auto g_0_xxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 168);

    auto g_0_xxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 169);

    auto g_0_xxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 170);

    auto g_0_xxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 171);

    auto g_0_xxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 172);

    auto g_0_xxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 173);

    auto g_0_xxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 174);

    auto g_0_xxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 175);

    auto g_0_xxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 176);

    auto g_0_xxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 177);

    auto g_0_xxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 178);

    auto g_0_xxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 179);

    auto g_0_xxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 180);

    auto g_0_xxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 181);

    auto g_0_xxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 182);

    auto g_0_xxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 183);

    auto g_0_xxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 184);

    auto g_0_xxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 185);

    auto g_0_xxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 186);

    auto g_0_xxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 187);

    auto g_0_xxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 188);

    auto g_0_xxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 189);

    auto g_0_xxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 190);

    auto g_0_xxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 191);

    auto g_0_xxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 192);

    auto g_0_xxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 193);

    auto g_0_xxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 194);

    auto g_0_xxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 195);

#pragma omp simd aligned(g_0_xxy_0_xxxxxx_0,       \
                             g_0_xxy_0_xxxxxx_1,   \
                             g_0_xxy_0_xxxxxz_0,   \
                             g_0_xxy_0_xxxxxz_1,   \
                             g_0_xxy_0_xxxxzz_0,   \
                             g_0_xxy_0_xxxxzz_1,   \
                             g_0_xxy_0_xxxzzz_0,   \
                             g_0_xxy_0_xxxzzz_1,   \
                             g_0_xxy_0_xxzzzz_0,   \
                             g_0_xxy_0_xxzzzz_1,   \
                             g_0_xxy_0_xzzzzz_0,   \
                             g_0_xxy_0_xzzzzz_1,   \
                             g_0_xxyy_0_xxxxxx_0,  \
                             g_0_xxyy_0_xxxxxx_1,  \
                             g_0_xxyy_0_xxxxxz_0,  \
                             g_0_xxyy_0_xxxxxz_1,  \
                             g_0_xxyy_0_xxxxzz_0,  \
                             g_0_xxyy_0_xxxxzz_1,  \
                             g_0_xxyy_0_xxxzzz_0,  \
                             g_0_xxyy_0_xxxzzz_1,  \
                             g_0_xxyy_0_xxzzzz_0,  \
                             g_0_xxyy_0_xxzzzz_1,  \
                             g_0_xxyy_0_xzzzzz_0,  \
                             g_0_xxyy_0_xzzzzz_1,  \
                             g_0_xxyyy_0_xxxxxx_0, \
                             g_0_xxyyy_0_xxxxxy_0, \
                             g_0_xxyyy_0_xxxxxz_0, \
                             g_0_xxyyy_0_xxxxyy_0, \
                             g_0_xxyyy_0_xxxxyz_0, \
                             g_0_xxyyy_0_xxxxzz_0, \
                             g_0_xxyyy_0_xxxyyy_0, \
                             g_0_xxyyy_0_xxxyyz_0, \
                             g_0_xxyyy_0_xxxyzz_0, \
                             g_0_xxyyy_0_xxxzzz_0, \
                             g_0_xxyyy_0_xxyyyy_0, \
                             g_0_xxyyy_0_xxyyyz_0, \
                             g_0_xxyyy_0_xxyyzz_0, \
                             g_0_xxyyy_0_xxyzzz_0, \
                             g_0_xxyyy_0_xxzzzz_0, \
                             g_0_xxyyy_0_xyyyyy_0, \
                             g_0_xxyyy_0_xyyyyz_0, \
                             g_0_xxyyy_0_xyyyzz_0, \
                             g_0_xxyyy_0_xyyzzz_0, \
                             g_0_xxyyy_0_xyzzzz_0, \
                             g_0_xxyyy_0_xzzzzz_0, \
                             g_0_xxyyy_0_yyyyyy_0, \
                             g_0_xxyyy_0_yyyyyz_0, \
                             g_0_xxyyy_0_yyyyzz_0, \
                             g_0_xxyyy_0_yyyzzz_0, \
                             g_0_xxyyy_0_yyzzzz_0, \
                             g_0_xxyyy_0_yzzzzz_0, \
                             g_0_xxyyy_0_zzzzzz_0, \
                             g_0_xyyy_0_xxxxxy_0,  \
                             g_0_xyyy_0_xxxxxy_1,  \
                             g_0_xyyy_0_xxxxy_1,   \
                             g_0_xyyy_0_xxxxyy_0,  \
                             g_0_xyyy_0_xxxxyy_1,  \
                             g_0_xyyy_0_xxxxyz_0,  \
                             g_0_xyyy_0_xxxxyz_1,  \
                             g_0_xyyy_0_xxxyy_1,   \
                             g_0_xyyy_0_xxxyyy_0,  \
                             g_0_xyyy_0_xxxyyy_1,  \
                             g_0_xyyy_0_xxxyyz_0,  \
                             g_0_xyyy_0_xxxyyz_1,  \
                             g_0_xyyy_0_xxxyz_1,   \
                             g_0_xyyy_0_xxxyzz_0,  \
                             g_0_xyyy_0_xxxyzz_1,  \
                             g_0_xyyy_0_xxyyy_1,   \
                             g_0_xyyy_0_xxyyyy_0,  \
                             g_0_xyyy_0_xxyyyy_1,  \
                             g_0_xyyy_0_xxyyyz_0,  \
                             g_0_xyyy_0_xxyyyz_1,  \
                             g_0_xyyy_0_xxyyz_1,   \
                             g_0_xyyy_0_xxyyzz_0,  \
                             g_0_xyyy_0_xxyyzz_1,  \
                             g_0_xyyy_0_xxyzz_1,   \
                             g_0_xyyy_0_xxyzzz_0,  \
                             g_0_xyyy_0_xxyzzz_1,  \
                             g_0_xyyy_0_xyyyy_1,   \
                             g_0_xyyy_0_xyyyyy_0,  \
                             g_0_xyyy_0_xyyyyy_1,  \
                             g_0_xyyy_0_xyyyyz_0,  \
                             g_0_xyyy_0_xyyyyz_1,  \
                             g_0_xyyy_0_xyyyz_1,   \
                             g_0_xyyy_0_xyyyzz_0,  \
                             g_0_xyyy_0_xyyyzz_1,  \
                             g_0_xyyy_0_xyyzz_1,   \
                             g_0_xyyy_0_xyyzzz_0,  \
                             g_0_xyyy_0_xyyzzz_1,  \
                             g_0_xyyy_0_xyzzz_1,   \
                             g_0_xyyy_0_xyzzzz_0,  \
                             g_0_xyyy_0_xyzzzz_1,  \
                             g_0_xyyy_0_yyyyy_1,   \
                             g_0_xyyy_0_yyyyyy_0,  \
                             g_0_xyyy_0_yyyyyy_1,  \
                             g_0_xyyy_0_yyyyyz_0,  \
                             g_0_xyyy_0_yyyyyz_1,  \
                             g_0_xyyy_0_yyyyz_1,   \
                             g_0_xyyy_0_yyyyzz_0,  \
                             g_0_xyyy_0_yyyyzz_1,  \
                             g_0_xyyy_0_yyyzz_1,   \
                             g_0_xyyy_0_yyyzzz_0,  \
                             g_0_xyyy_0_yyyzzz_1,  \
                             g_0_xyyy_0_yyzzz_1,   \
                             g_0_xyyy_0_yyzzzz_0,  \
                             g_0_xyyy_0_yyzzzz_1,  \
                             g_0_xyyy_0_yzzzz_1,   \
                             g_0_xyyy_0_yzzzzz_0,  \
                             g_0_xyyy_0_yzzzzz_1,  \
                             g_0_xyyy_0_zzzzzz_0,  \
                             g_0_xyyy_0_zzzzzz_1,  \
                             g_0_yyy_0_xxxxxy_0,   \
                             g_0_yyy_0_xxxxxy_1,   \
                             g_0_yyy_0_xxxxyy_0,   \
                             g_0_yyy_0_xxxxyy_1,   \
                             g_0_yyy_0_xxxxyz_0,   \
                             g_0_yyy_0_xxxxyz_1,   \
                             g_0_yyy_0_xxxyyy_0,   \
                             g_0_yyy_0_xxxyyy_1,   \
                             g_0_yyy_0_xxxyyz_0,   \
                             g_0_yyy_0_xxxyyz_1,   \
                             g_0_yyy_0_xxxyzz_0,   \
                             g_0_yyy_0_xxxyzz_1,   \
                             g_0_yyy_0_xxyyyy_0,   \
                             g_0_yyy_0_xxyyyy_1,   \
                             g_0_yyy_0_xxyyyz_0,   \
                             g_0_yyy_0_xxyyyz_1,   \
                             g_0_yyy_0_xxyyzz_0,   \
                             g_0_yyy_0_xxyyzz_1,   \
                             g_0_yyy_0_xxyzzz_0,   \
                             g_0_yyy_0_xxyzzz_1,   \
                             g_0_yyy_0_xyyyyy_0,   \
                             g_0_yyy_0_xyyyyy_1,   \
                             g_0_yyy_0_xyyyyz_0,   \
                             g_0_yyy_0_xyyyyz_1,   \
                             g_0_yyy_0_xyyyzz_0,   \
                             g_0_yyy_0_xyyyzz_1,   \
                             g_0_yyy_0_xyyzzz_0,   \
                             g_0_yyy_0_xyyzzz_1,   \
                             g_0_yyy_0_xyzzzz_0,   \
                             g_0_yyy_0_xyzzzz_1,   \
                             g_0_yyy_0_yyyyyy_0,   \
                             g_0_yyy_0_yyyyyy_1,   \
                             g_0_yyy_0_yyyyyz_0,   \
                             g_0_yyy_0_yyyyyz_1,   \
                             g_0_yyy_0_yyyyzz_0,   \
                             g_0_yyy_0_yyyyzz_1,   \
                             g_0_yyy_0_yyyzzz_0,   \
                             g_0_yyy_0_yyyzzz_1,   \
                             g_0_yyy_0_yyzzzz_0,   \
                             g_0_yyy_0_yyzzzz_1,   \
                             g_0_yyy_0_yzzzzz_0,   \
                             g_0_yyy_0_yzzzzz_1,   \
                             g_0_yyy_0_zzzzzz_0,   \
                             g_0_yyy_0_zzzzzz_1,   \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyy_0_xxxxxx_0[i] = 2.0 * g_0_xxy_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxx_0[i] * pb_y +
                                  g_0_xxyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxxy_0[i] = g_0_yyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxy_1[i] * fti_ab_0 + 5.0 * g_0_xyyy_0_xxxxy_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxxxxy_0[i] * pb_x + g_0_xyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxz_0[i] = 2.0 * g_0_xxy_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxz_0[i] * pb_y +
                                  g_0_xxyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxyy_0[i] = g_0_yyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyy_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyy_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxxxyy_0[i] * pb_x + g_0_xyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxyz_0[i] = g_0_yyy_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxxxyz_0[i] * pb_x + g_0_xyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxzz_0[i] = 2.0 * g_0_xxy_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxzz_0[i] * pb_y +
                                  g_0_xxyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxyyy_0[i] = g_0_yyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyy_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxxyyy_0[i] * pb_x + g_0_xyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyyz_0[i] = g_0_yyy_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxxyyz_0[i] * pb_x + g_0_xyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyzz_0[i] = g_0_yyy_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyzz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxxyzz_0[i] * pb_x + g_0_xyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxzzz_0[i] = 2.0 * g_0_xxy_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxzzz_0[i] * pb_y +
                                  g_0_xxyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxyyyy_0[i] = g_0_yyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyy_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxyyyy_0[i] * pb_x + g_0_xyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyyz_0[i] = g_0_yyy_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxyyyz_0[i] * pb_x + g_0_xyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyzz_0[i] = g_0_yyy_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyzz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxyyzz_0[i] * pb_x + g_0_xyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyzzz_0[i] = g_0_yyy_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyzzz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xxyzzz_0[i] * pb_x + g_0_xyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxzzzz_0[i] = 2.0 * g_0_xxy_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxzzzz_0[i] * pb_y +
                                  g_0_xxyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xyyyyy_0[i] = g_0_yyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyy_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xyyyyy_0[i] * pb_x + g_0_xyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyyz_0[i] = g_0_yyy_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xyyyyz_0[i] * pb_x + g_0_xyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyzz_0[i] = g_0_yyy_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyzz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xyyyzz_0[i] * pb_x + g_0_xyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyzzz_0[i] = g_0_yyy_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzzz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xyyzzz_0[i] * pb_x + g_0_xyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyzzzz_0[i] = g_0_yyy_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzzz_1[i] * fi_abcd_0 +
                                  g_0_xyyy_0_xyzzzz_0[i] * pb_x + g_0_xyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xzzzzz_0[i] = 2.0 * g_0_xxy_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xzzzzz_0[i] * pb_y +
                                  g_0_xxyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_yyyyyy_0[i] =
            g_0_yyy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyy_0[i] * pb_x + g_0_xyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyyz_0[i] =
            g_0_yyy_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyz_0[i] * pb_x + g_0_xyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyzz_0[i] =
            g_0_yyy_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyzz_0[i] * pb_x + g_0_xyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyzzz_0[i] =
            g_0_yyy_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyzzz_0[i] * pb_x + g_0_xyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyzzzz_0[i] =
            g_0_yyy_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzzzz_0[i] * pb_x + g_0_xyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yzzzzz_0[i] =
            g_0_yyy_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzzzz_0[i] * pb_x + g_0_xyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_zzzzzz_0[i] =
            g_0_yyy_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_zzzzzz_0[i] * pb_x + g_0_xyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 196-224 components of targeted buffer : SHSI

    auto g_0_xxyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 196);

    auto g_0_xxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 197);

    auto g_0_xxyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 198);

    auto g_0_xxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 199);

    auto g_0_xxyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 200);

    auto g_0_xxyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 201);

    auto g_0_xxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 202);

    auto g_0_xxyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 203);

    auto g_0_xxyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 204);

    auto g_0_xxyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 205);

    auto g_0_xxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 206);

    auto g_0_xxyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 207);

    auto g_0_xxyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 208);

    auto g_0_xxyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 209);

    auto g_0_xxyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 210);

    auto g_0_xxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 211);

    auto g_0_xxyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 212);

    auto g_0_xxyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 213);

    auto g_0_xxyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 214);

    auto g_0_xxyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 215);

    auto g_0_xxyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 216);

    auto g_0_xxyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 217);

    auto g_0_xxyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 218);

    auto g_0_xxyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 219);

    auto g_0_xxyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 220);

    auto g_0_xxyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 221);

    auto g_0_xxyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 222);

    auto g_0_xxyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 223);

#pragma omp simd aligned(g_0_xxyy_0_xxxxx_1,       \
                             g_0_xxyy_0_xxxxxx_0,  \
                             g_0_xxyy_0_xxxxxx_1,  \
                             g_0_xxyy_0_xxxxxy_0,  \
                             g_0_xxyy_0_xxxxxy_1,  \
                             g_0_xxyy_0_xxxxxz_0,  \
                             g_0_xxyy_0_xxxxxz_1,  \
                             g_0_xxyy_0_xxxxy_1,   \
                             g_0_xxyy_0_xxxxyy_0,  \
                             g_0_xxyy_0_xxxxyy_1,  \
                             g_0_xxyy_0_xxxxyz_0,  \
                             g_0_xxyy_0_xxxxyz_1,  \
                             g_0_xxyy_0_xxxxz_1,   \
                             g_0_xxyy_0_xxxxzz_0,  \
                             g_0_xxyy_0_xxxxzz_1,  \
                             g_0_xxyy_0_xxxyy_1,   \
                             g_0_xxyy_0_xxxyyy_0,  \
                             g_0_xxyy_0_xxxyyy_1,  \
                             g_0_xxyy_0_xxxyyz_0,  \
                             g_0_xxyy_0_xxxyyz_1,  \
                             g_0_xxyy_0_xxxyz_1,   \
                             g_0_xxyy_0_xxxyzz_0,  \
                             g_0_xxyy_0_xxxyzz_1,  \
                             g_0_xxyy_0_xxxzz_1,   \
                             g_0_xxyy_0_xxxzzz_0,  \
                             g_0_xxyy_0_xxxzzz_1,  \
                             g_0_xxyy_0_xxyyy_1,   \
                             g_0_xxyy_0_xxyyyy_0,  \
                             g_0_xxyy_0_xxyyyy_1,  \
                             g_0_xxyy_0_xxyyyz_0,  \
                             g_0_xxyy_0_xxyyyz_1,  \
                             g_0_xxyy_0_xxyyz_1,   \
                             g_0_xxyy_0_xxyyzz_0,  \
                             g_0_xxyy_0_xxyyzz_1,  \
                             g_0_xxyy_0_xxyzz_1,   \
                             g_0_xxyy_0_xxyzzz_0,  \
                             g_0_xxyy_0_xxyzzz_1,  \
                             g_0_xxyy_0_xxzzz_1,   \
                             g_0_xxyy_0_xxzzzz_0,  \
                             g_0_xxyy_0_xxzzzz_1,  \
                             g_0_xxyy_0_xyyyy_1,   \
                             g_0_xxyy_0_xyyyyy_0,  \
                             g_0_xxyy_0_xyyyyy_1,  \
                             g_0_xxyy_0_xyyyyz_0,  \
                             g_0_xxyy_0_xyyyyz_1,  \
                             g_0_xxyy_0_xyyyz_1,   \
                             g_0_xxyy_0_xyyyzz_0,  \
                             g_0_xxyy_0_xyyyzz_1,  \
                             g_0_xxyy_0_xyyzz_1,   \
                             g_0_xxyy_0_xyyzzz_0,  \
                             g_0_xxyy_0_xyyzzz_1,  \
                             g_0_xxyy_0_xyzzz_1,   \
                             g_0_xxyy_0_xyzzzz_0,  \
                             g_0_xxyy_0_xyzzzz_1,  \
                             g_0_xxyy_0_xzzzz_1,   \
                             g_0_xxyy_0_xzzzzz_0,  \
                             g_0_xxyy_0_xzzzzz_1,  \
                             g_0_xxyy_0_yyyyy_1,   \
                             g_0_xxyy_0_yyyyyy_0,  \
                             g_0_xxyy_0_yyyyyy_1,  \
                             g_0_xxyy_0_yyyyyz_0,  \
                             g_0_xxyy_0_yyyyyz_1,  \
                             g_0_xxyy_0_yyyyz_1,   \
                             g_0_xxyy_0_yyyyzz_0,  \
                             g_0_xxyy_0_yyyyzz_1,  \
                             g_0_xxyy_0_yyyzz_1,   \
                             g_0_xxyy_0_yyyzzz_0,  \
                             g_0_xxyy_0_yyyzzz_1,  \
                             g_0_xxyy_0_yyzzz_1,   \
                             g_0_xxyy_0_yyzzzz_0,  \
                             g_0_xxyy_0_yyzzzz_1,  \
                             g_0_xxyy_0_yzzzz_1,   \
                             g_0_xxyy_0_yzzzzz_0,  \
                             g_0_xxyy_0_yzzzzz_1,  \
                             g_0_xxyy_0_zzzzz_1,   \
                             g_0_xxyy_0_zzzzzz_0,  \
                             g_0_xxyy_0_zzzzzz_1,  \
                             g_0_xxyyz_0_xxxxxx_0, \
                             g_0_xxyyz_0_xxxxxy_0, \
                             g_0_xxyyz_0_xxxxxz_0, \
                             g_0_xxyyz_0_xxxxyy_0, \
                             g_0_xxyyz_0_xxxxyz_0, \
                             g_0_xxyyz_0_xxxxzz_0, \
                             g_0_xxyyz_0_xxxyyy_0, \
                             g_0_xxyyz_0_xxxyyz_0, \
                             g_0_xxyyz_0_xxxyzz_0, \
                             g_0_xxyyz_0_xxxzzz_0, \
                             g_0_xxyyz_0_xxyyyy_0, \
                             g_0_xxyyz_0_xxyyyz_0, \
                             g_0_xxyyz_0_xxyyzz_0, \
                             g_0_xxyyz_0_xxyzzz_0, \
                             g_0_xxyyz_0_xxzzzz_0, \
                             g_0_xxyyz_0_xyyyyy_0, \
                             g_0_xxyyz_0_xyyyyz_0, \
                             g_0_xxyyz_0_xyyyzz_0, \
                             g_0_xxyyz_0_xyyzzz_0, \
                             g_0_xxyyz_0_xyzzzz_0, \
                             g_0_xxyyz_0_xzzzzz_0, \
                             g_0_xxyyz_0_yyyyyy_0, \
                             g_0_xxyyz_0_yyyyyz_0, \
                             g_0_xxyyz_0_yyyyzz_0, \
                             g_0_xxyyz_0_yyyzzz_0, \
                             g_0_xxyyz_0_yyzzzz_0, \
                             g_0_xxyyz_0_yzzzzz_0, \
                             g_0_xxyyz_0_zzzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyz_0_xxxxxx_0[i] = g_0_xxyy_0_xxxxxx_0[i] * pb_z + g_0_xxyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxy_0[i] = g_0_xxyy_0_xxxxxy_0[i] * pb_z + g_0_xxyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxz_0[i] = g_0_xxyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxz_0[i] * pb_z + g_0_xxyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyy_0[i] = g_0_xxyy_0_xxxxyy_0[i] * pb_z + g_0_xxyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyz_0[i] = g_0_xxyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyz_0[i] * pb_z + g_0_xxyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxzz_0[i] = 2.0 * g_0_xxyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxzz_0[i] * pb_z + g_0_xxyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyy_0[i] = g_0_xxyy_0_xxxyyy_0[i] * pb_z + g_0_xxyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyz_0[i] = g_0_xxyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyz_0[i] * pb_z + g_0_xxyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyzz_0[i] = 2.0 * g_0_xxyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyzz_0[i] * pb_z + g_0_xxyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxzzz_0[i] = 3.0 * g_0_xxyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxzzz_0[i] * pb_z + g_0_xxyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyy_0[i] = g_0_xxyy_0_xxyyyy_0[i] * pb_z + g_0_xxyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyz_0[i] = g_0_xxyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyz_0[i] * pb_z + g_0_xxyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyzz_0[i] = 2.0 * g_0_xxyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyzz_0[i] * pb_z + g_0_xxyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyzzz_0[i] = 3.0 * g_0_xxyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyzzz_0[i] * pb_z + g_0_xxyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxzzzz_0[i] = 4.0 * g_0_xxyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxzzzz_0[i] * pb_z + g_0_xxyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyy_0[i] = g_0_xxyy_0_xyyyyy_0[i] * pb_z + g_0_xxyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyz_0[i] = g_0_xxyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyz_0[i] * pb_z + g_0_xxyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyzz_0[i] = 2.0 * g_0_xxyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyzz_0[i] * pb_z + g_0_xxyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyzzz_0[i] = 3.0 * g_0_xxyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyzzz_0[i] * pb_z + g_0_xxyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyzzzz_0[i] = 4.0 * g_0_xxyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyzzzz_0[i] * pb_z + g_0_xxyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xzzzzz_0[i] = 5.0 * g_0_xxyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xzzzzz_0[i] * pb_z + g_0_xxyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyy_0[i] = g_0_xxyy_0_yyyyyy_0[i] * pb_z + g_0_xxyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyz_0[i] = g_0_xxyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyyz_0[i] * pb_z + g_0_xxyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyzz_0[i] = 2.0 * g_0_xxyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyzz_0[i] * pb_z + g_0_xxyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyzzz_0[i] = 3.0 * g_0_xxyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyzzz_0[i] * pb_z + g_0_xxyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyzzzz_0[i] = 4.0 * g_0_xxyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyzzzz_0[i] * pb_z + g_0_xxyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yzzzzz_0[i] = 5.0 * g_0_xxyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yzzzzz_0[i] * pb_z + g_0_xxyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_zzzzzz_0[i] = 6.0 * g_0_xxyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_zzzzzz_0[i] * pb_z + g_0_xxyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 224-252 components of targeted buffer : SHSI

    auto g_0_xxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 224);

    auto g_0_xxyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 225);

    auto g_0_xxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 226);

    auto g_0_xxyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 227);

    auto g_0_xxyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 228);

    auto g_0_xxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 229);

    auto g_0_xxyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 230);

    auto g_0_xxyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 231);

    auto g_0_xxyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 232);

    auto g_0_xxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 233);

    auto g_0_xxyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 234);

    auto g_0_xxyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 235);

    auto g_0_xxyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 236);

    auto g_0_xxyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 237);

    auto g_0_xxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 238);

    auto g_0_xxyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 239);

    auto g_0_xxyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 240);

    auto g_0_xxyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 241);

    auto g_0_xxyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 242);

    auto g_0_xxyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 243);

    auto g_0_xxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 244);

    auto g_0_xxyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 245);

    auto g_0_xxyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 246);

    auto g_0_xxyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 247);

    auto g_0_xxyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 248);

    auto g_0_xxyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 249);

    auto g_0_xxyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 250);

    auto g_0_xxyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 251);

#pragma omp simd aligned(g_0_xxyzz_0_xxxxxx_0,     \
                             g_0_xxyzz_0_xxxxxy_0, \
                             g_0_xxyzz_0_xxxxxz_0, \
                             g_0_xxyzz_0_xxxxyy_0, \
                             g_0_xxyzz_0_xxxxyz_0, \
                             g_0_xxyzz_0_xxxxzz_0, \
                             g_0_xxyzz_0_xxxyyy_0, \
                             g_0_xxyzz_0_xxxyyz_0, \
                             g_0_xxyzz_0_xxxyzz_0, \
                             g_0_xxyzz_0_xxxzzz_0, \
                             g_0_xxyzz_0_xxyyyy_0, \
                             g_0_xxyzz_0_xxyyyz_0, \
                             g_0_xxyzz_0_xxyyzz_0, \
                             g_0_xxyzz_0_xxyzzz_0, \
                             g_0_xxyzz_0_xxzzzz_0, \
                             g_0_xxyzz_0_xyyyyy_0, \
                             g_0_xxyzz_0_xyyyyz_0, \
                             g_0_xxyzz_0_xyyyzz_0, \
                             g_0_xxyzz_0_xyyzzz_0, \
                             g_0_xxyzz_0_xyzzzz_0, \
                             g_0_xxyzz_0_xzzzzz_0, \
                             g_0_xxyzz_0_yyyyyy_0, \
                             g_0_xxyzz_0_yyyyyz_0, \
                             g_0_xxyzz_0_yyyyzz_0, \
                             g_0_xxyzz_0_yyyzzz_0, \
                             g_0_xxyzz_0_yyzzzz_0, \
                             g_0_xxyzz_0_yzzzzz_0, \
                             g_0_xxyzz_0_zzzzzz_0, \
                             g_0_xxzz_0_xxxxx_1,   \
                             g_0_xxzz_0_xxxxxx_0,  \
                             g_0_xxzz_0_xxxxxx_1,  \
                             g_0_xxzz_0_xxxxxy_0,  \
                             g_0_xxzz_0_xxxxxy_1,  \
                             g_0_xxzz_0_xxxxxz_0,  \
                             g_0_xxzz_0_xxxxxz_1,  \
                             g_0_xxzz_0_xxxxy_1,   \
                             g_0_xxzz_0_xxxxyy_0,  \
                             g_0_xxzz_0_xxxxyy_1,  \
                             g_0_xxzz_0_xxxxyz_0,  \
                             g_0_xxzz_0_xxxxyz_1,  \
                             g_0_xxzz_0_xxxxz_1,   \
                             g_0_xxzz_0_xxxxzz_0,  \
                             g_0_xxzz_0_xxxxzz_1,  \
                             g_0_xxzz_0_xxxyy_1,   \
                             g_0_xxzz_0_xxxyyy_0,  \
                             g_0_xxzz_0_xxxyyy_1,  \
                             g_0_xxzz_0_xxxyyz_0,  \
                             g_0_xxzz_0_xxxyyz_1,  \
                             g_0_xxzz_0_xxxyz_1,   \
                             g_0_xxzz_0_xxxyzz_0,  \
                             g_0_xxzz_0_xxxyzz_1,  \
                             g_0_xxzz_0_xxxzz_1,   \
                             g_0_xxzz_0_xxxzzz_0,  \
                             g_0_xxzz_0_xxxzzz_1,  \
                             g_0_xxzz_0_xxyyy_1,   \
                             g_0_xxzz_0_xxyyyy_0,  \
                             g_0_xxzz_0_xxyyyy_1,  \
                             g_0_xxzz_0_xxyyyz_0,  \
                             g_0_xxzz_0_xxyyyz_1,  \
                             g_0_xxzz_0_xxyyz_1,   \
                             g_0_xxzz_0_xxyyzz_0,  \
                             g_0_xxzz_0_xxyyzz_1,  \
                             g_0_xxzz_0_xxyzz_1,   \
                             g_0_xxzz_0_xxyzzz_0,  \
                             g_0_xxzz_0_xxyzzz_1,  \
                             g_0_xxzz_0_xxzzz_1,   \
                             g_0_xxzz_0_xxzzzz_0,  \
                             g_0_xxzz_0_xxzzzz_1,  \
                             g_0_xxzz_0_xyyyy_1,   \
                             g_0_xxzz_0_xyyyyy_0,  \
                             g_0_xxzz_0_xyyyyy_1,  \
                             g_0_xxzz_0_xyyyyz_0,  \
                             g_0_xxzz_0_xyyyyz_1,  \
                             g_0_xxzz_0_xyyyz_1,   \
                             g_0_xxzz_0_xyyyzz_0,  \
                             g_0_xxzz_0_xyyyzz_1,  \
                             g_0_xxzz_0_xyyzz_1,   \
                             g_0_xxzz_0_xyyzzz_0,  \
                             g_0_xxzz_0_xyyzzz_1,  \
                             g_0_xxzz_0_xyzzz_1,   \
                             g_0_xxzz_0_xyzzzz_0,  \
                             g_0_xxzz_0_xyzzzz_1,  \
                             g_0_xxzz_0_xzzzz_1,   \
                             g_0_xxzz_0_xzzzzz_0,  \
                             g_0_xxzz_0_xzzzzz_1,  \
                             g_0_xxzz_0_yyyyy_1,   \
                             g_0_xxzz_0_yyyyyy_0,  \
                             g_0_xxzz_0_yyyyyy_1,  \
                             g_0_xxzz_0_yyyyyz_0,  \
                             g_0_xxzz_0_yyyyyz_1,  \
                             g_0_xxzz_0_yyyyz_1,   \
                             g_0_xxzz_0_yyyyzz_0,  \
                             g_0_xxzz_0_yyyyzz_1,  \
                             g_0_xxzz_0_yyyzz_1,   \
                             g_0_xxzz_0_yyyzzz_0,  \
                             g_0_xxzz_0_yyyzzz_1,  \
                             g_0_xxzz_0_yyzzz_1,   \
                             g_0_xxzz_0_yyzzzz_0,  \
                             g_0_xxzz_0_yyzzzz_1,  \
                             g_0_xxzz_0_yzzzz_1,   \
                             g_0_xxzz_0_yzzzzz_0,  \
                             g_0_xxzz_0_yzzzzz_1,  \
                             g_0_xxzz_0_zzzzz_1,   \
                             g_0_xxzz_0_zzzzzz_0,  \
                             g_0_xxzz_0_zzzzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzz_0_xxxxxx_0[i] = g_0_xxzz_0_xxxxxx_0[i] * pb_y + g_0_xxzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxy_0[i] = g_0_xxzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxy_0[i] * pb_y + g_0_xxzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxz_0[i] = g_0_xxzz_0_xxxxxz_0[i] * pb_y + g_0_xxzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyy_0[i] = 2.0 * g_0_xxzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyy_0[i] * pb_y + g_0_xxzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyz_0[i] = g_0_xxzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyz_0[i] * pb_y + g_0_xxzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxzz_0[i] = g_0_xxzz_0_xxxxzz_0[i] * pb_y + g_0_xxzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyy_0[i] = 3.0 * g_0_xxzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyy_0[i] * pb_y + g_0_xxzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyz_0[i] = 2.0 * g_0_xxzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyz_0[i] * pb_y + g_0_xxzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyzz_0[i] = g_0_xxzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyzz_0[i] * pb_y + g_0_xxzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxzzz_0[i] = g_0_xxzz_0_xxxzzz_0[i] * pb_y + g_0_xxzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyy_0[i] = 4.0 * g_0_xxzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyy_0[i] * pb_y + g_0_xxzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyz_0[i] = 3.0 * g_0_xxzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyz_0[i] * pb_y + g_0_xxzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyzz_0[i] = 2.0 * g_0_xxzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyzz_0[i] * pb_y + g_0_xxzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyzzz_0[i] = g_0_xxzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyzzz_0[i] * pb_y + g_0_xxzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxzzzz_0[i] = g_0_xxzz_0_xxzzzz_0[i] * pb_y + g_0_xxzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyy_0[i] = 5.0 * g_0_xxzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyy_0[i] * pb_y + g_0_xxzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyz_0[i] = 4.0 * g_0_xxzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyz_0[i] * pb_y + g_0_xxzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyzz_0[i] = 3.0 * g_0_xxzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyzz_0[i] * pb_y + g_0_xxzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyzzz_0[i] = 2.0 * g_0_xxzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyzzz_0[i] * pb_y + g_0_xxzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyzzzz_0[i] = g_0_xxzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyzzzz_0[i] * pb_y + g_0_xxzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xzzzzz_0[i] = g_0_xxzz_0_xzzzzz_0[i] * pb_y + g_0_xxzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyy_0[i] = 6.0 * g_0_xxzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyy_0[i] * pb_y + g_0_xxzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyz_0[i] = 5.0 * g_0_xxzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyz_0[i] * pb_y + g_0_xxzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyzz_0[i] = 4.0 * g_0_xxzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyzz_0[i] * pb_y + g_0_xxzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyzzz_0[i] = 3.0 * g_0_xxzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyzzz_0[i] * pb_y + g_0_xxzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyzzzz_0[i] = 2.0 * g_0_xxzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyzzzz_0[i] * pb_y + g_0_xxzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yzzzzz_0[i] = g_0_xxzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yzzzzz_0[i] * pb_y + g_0_xxzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_zzzzzz_0[i] = g_0_xxzz_0_zzzzzz_0[i] * pb_y + g_0_xxzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 252-280 components of targeted buffer : SHSI

    auto g_0_xxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 252);

    auto g_0_xxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 253);

    auto g_0_xxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 254);

    auto g_0_xxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 255);

    auto g_0_xxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 256);

    auto g_0_xxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 257);

    auto g_0_xxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 258);

    auto g_0_xxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 259);

    auto g_0_xxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 260);

    auto g_0_xxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 261);

    auto g_0_xxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 262);

    auto g_0_xxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 263);

    auto g_0_xxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 264);

    auto g_0_xxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 265);

    auto g_0_xxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 266);

    auto g_0_xxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 267);

    auto g_0_xxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 268);

    auto g_0_xxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 269);

    auto g_0_xxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 270);

    auto g_0_xxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 271);

    auto g_0_xxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 272);

    auto g_0_xxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 273);

    auto g_0_xxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 274);

    auto g_0_xxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 275);

    auto g_0_xxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 276);

    auto g_0_xxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 277);

    auto g_0_xxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 278);

    auto g_0_xxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 279);

#pragma omp simd aligned(g_0_xxz_0_xxxxxx_0,       \
                             g_0_xxz_0_xxxxxx_1,   \
                             g_0_xxz_0_xxxxxy_0,   \
                             g_0_xxz_0_xxxxxy_1,   \
                             g_0_xxz_0_xxxxyy_0,   \
                             g_0_xxz_0_xxxxyy_1,   \
                             g_0_xxz_0_xxxyyy_0,   \
                             g_0_xxz_0_xxxyyy_1,   \
                             g_0_xxz_0_xxyyyy_0,   \
                             g_0_xxz_0_xxyyyy_1,   \
                             g_0_xxz_0_xyyyyy_0,   \
                             g_0_xxz_0_xyyyyy_1,   \
                             g_0_xxzz_0_xxxxxx_0,  \
                             g_0_xxzz_0_xxxxxx_1,  \
                             g_0_xxzz_0_xxxxxy_0,  \
                             g_0_xxzz_0_xxxxxy_1,  \
                             g_0_xxzz_0_xxxxyy_0,  \
                             g_0_xxzz_0_xxxxyy_1,  \
                             g_0_xxzz_0_xxxyyy_0,  \
                             g_0_xxzz_0_xxxyyy_1,  \
                             g_0_xxzz_0_xxyyyy_0,  \
                             g_0_xxzz_0_xxyyyy_1,  \
                             g_0_xxzz_0_xyyyyy_0,  \
                             g_0_xxzz_0_xyyyyy_1,  \
                             g_0_xxzzz_0_xxxxxx_0, \
                             g_0_xxzzz_0_xxxxxy_0, \
                             g_0_xxzzz_0_xxxxxz_0, \
                             g_0_xxzzz_0_xxxxyy_0, \
                             g_0_xxzzz_0_xxxxyz_0, \
                             g_0_xxzzz_0_xxxxzz_0, \
                             g_0_xxzzz_0_xxxyyy_0, \
                             g_0_xxzzz_0_xxxyyz_0, \
                             g_0_xxzzz_0_xxxyzz_0, \
                             g_0_xxzzz_0_xxxzzz_0, \
                             g_0_xxzzz_0_xxyyyy_0, \
                             g_0_xxzzz_0_xxyyyz_0, \
                             g_0_xxzzz_0_xxyyzz_0, \
                             g_0_xxzzz_0_xxyzzz_0, \
                             g_0_xxzzz_0_xxzzzz_0, \
                             g_0_xxzzz_0_xyyyyy_0, \
                             g_0_xxzzz_0_xyyyyz_0, \
                             g_0_xxzzz_0_xyyyzz_0, \
                             g_0_xxzzz_0_xyyzzz_0, \
                             g_0_xxzzz_0_xyzzzz_0, \
                             g_0_xxzzz_0_xzzzzz_0, \
                             g_0_xxzzz_0_yyyyyy_0, \
                             g_0_xxzzz_0_yyyyyz_0, \
                             g_0_xxzzz_0_yyyyzz_0, \
                             g_0_xxzzz_0_yyyzzz_0, \
                             g_0_xxzzz_0_yyzzzz_0, \
                             g_0_xxzzz_0_yzzzzz_0, \
                             g_0_xxzzz_0_zzzzzz_0, \
                             g_0_xzzz_0_xxxxxz_0,  \
                             g_0_xzzz_0_xxxxxz_1,  \
                             g_0_xzzz_0_xxxxyz_0,  \
                             g_0_xzzz_0_xxxxyz_1,  \
                             g_0_xzzz_0_xxxxz_1,   \
                             g_0_xzzz_0_xxxxzz_0,  \
                             g_0_xzzz_0_xxxxzz_1,  \
                             g_0_xzzz_0_xxxyyz_0,  \
                             g_0_xzzz_0_xxxyyz_1,  \
                             g_0_xzzz_0_xxxyz_1,   \
                             g_0_xzzz_0_xxxyzz_0,  \
                             g_0_xzzz_0_xxxyzz_1,  \
                             g_0_xzzz_0_xxxzz_1,   \
                             g_0_xzzz_0_xxxzzz_0,  \
                             g_0_xzzz_0_xxxzzz_1,  \
                             g_0_xzzz_0_xxyyyz_0,  \
                             g_0_xzzz_0_xxyyyz_1,  \
                             g_0_xzzz_0_xxyyz_1,   \
                             g_0_xzzz_0_xxyyzz_0,  \
                             g_0_xzzz_0_xxyyzz_1,  \
                             g_0_xzzz_0_xxyzz_1,   \
                             g_0_xzzz_0_xxyzzz_0,  \
                             g_0_xzzz_0_xxyzzz_1,  \
                             g_0_xzzz_0_xxzzz_1,   \
                             g_0_xzzz_0_xxzzzz_0,  \
                             g_0_xzzz_0_xxzzzz_1,  \
                             g_0_xzzz_0_xyyyyz_0,  \
                             g_0_xzzz_0_xyyyyz_1,  \
                             g_0_xzzz_0_xyyyz_1,   \
                             g_0_xzzz_0_xyyyzz_0,  \
                             g_0_xzzz_0_xyyyzz_1,  \
                             g_0_xzzz_0_xyyzz_1,   \
                             g_0_xzzz_0_xyyzzz_0,  \
                             g_0_xzzz_0_xyyzzz_1,  \
                             g_0_xzzz_0_xyzzz_1,   \
                             g_0_xzzz_0_xyzzzz_0,  \
                             g_0_xzzz_0_xyzzzz_1,  \
                             g_0_xzzz_0_xzzzz_1,   \
                             g_0_xzzz_0_xzzzzz_0,  \
                             g_0_xzzz_0_xzzzzz_1,  \
                             g_0_xzzz_0_yyyyyy_0,  \
                             g_0_xzzz_0_yyyyyy_1,  \
                             g_0_xzzz_0_yyyyyz_0,  \
                             g_0_xzzz_0_yyyyyz_1,  \
                             g_0_xzzz_0_yyyyz_1,   \
                             g_0_xzzz_0_yyyyzz_0,  \
                             g_0_xzzz_0_yyyyzz_1,  \
                             g_0_xzzz_0_yyyzz_1,   \
                             g_0_xzzz_0_yyyzzz_0,  \
                             g_0_xzzz_0_yyyzzz_1,  \
                             g_0_xzzz_0_yyzzz_1,   \
                             g_0_xzzz_0_yyzzzz_0,  \
                             g_0_xzzz_0_yyzzzz_1,  \
                             g_0_xzzz_0_yzzzz_1,   \
                             g_0_xzzz_0_yzzzzz_0,  \
                             g_0_xzzz_0_yzzzzz_1,  \
                             g_0_xzzz_0_zzzzz_1,   \
                             g_0_xzzz_0_zzzzzz_0,  \
                             g_0_xzzz_0_zzzzzz_1,  \
                             g_0_zzz_0_xxxxxz_0,   \
                             g_0_zzz_0_xxxxxz_1,   \
                             g_0_zzz_0_xxxxyz_0,   \
                             g_0_zzz_0_xxxxyz_1,   \
                             g_0_zzz_0_xxxxzz_0,   \
                             g_0_zzz_0_xxxxzz_1,   \
                             g_0_zzz_0_xxxyyz_0,   \
                             g_0_zzz_0_xxxyyz_1,   \
                             g_0_zzz_0_xxxyzz_0,   \
                             g_0_zzz_0_xxxyzz_1,   \
                             g_0_zzz_0_xxxzzz_0,   \
                             g_0_zzz_0_xxxzzz_1,   \
                             g_0_zzz_0_xxyyyz_0,   \
                             g_0_zzz_0_xxyyyz_1,   \
                             g_0_zzz_0_xxyyzz_0,   \
                             g_0_zzz_0_xxyyzz_1,   \
                             g_0_zzz_0_xxyzzz_0,   \
                             g_0_zzz_0_xxyzzz_1,   \
                             g_0_zzz_0_xxzzzz_0,   \
                             g_0_zzz_0_xxzzzz_1,   \
                             g_0_zzz_0_xyyyyz_0,   \
                             g_0_zzz_0_xyyyyz_1,   \
                             g_0_zzz_0_xyyyzz_0,   \
                             g_0_zzz_0_xyyyzz_1,   \
                             g_0_zzz_0_xyyzzz_0,   \
                             g_0_zzz_0_xyyzzz_1,   \
                             g_0_zzz_0_xyzzzz_0,   \
                             g_0_zzz_0_xyzzzz_1,   \
                             g_0_zzz_0_xzzzzz_0,   \
                             g_0_zzz_0_xzzzzz_1,   \
                             g_0_zzz_0_yyyyyy_0,   \
                             g_0_zzz_0_yyyyyy_1,   \
                             g_0_zzz_0_yyyyyz_0,   \
                             g_0_zzz_0_yyyyyz_1,   \
                             g_0_zzz_0_yyyyzz_0,   \
                             g_0_zzz_0_yyyyzz_1,   \
                             g_0_zzz_0_yyyzzz_0,   \
                             g_0_zzz_0_yyyzzz_1,   \
                             g_0_zzz_0_yyzzzz_0,   \
                             g_0_zzz_0_yyzzzz_1,   \
                             g_0_zzz_0_yzzzzz_0,   \
                             g_0_zzz_0_yzzzzz_1,   \
                             g_0_zzz_0_zzzzzz_0,   \
                             g_0_zzz_0_zzzzzz_1,   \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzz_0_xxxxxx_0[i] = 2.0 * g_0_xxz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxx_0[i] * pb_z +
                                  g_0_xxzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxy_0[i] = 2.0 * g_0_xxz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxy_0[i] * pb_z +
                                  g_0_xxzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxz_0[i] = g_0_zzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxz_1[i] * fti_ab_0 + 5.0 * g_0_xzzz_0_xxxxz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxxxxz_0[i] * pb_x + g_0_xzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxyy_0[i] * pb_z +
                                  g_0_xxzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxyz_0[i] = g_0_zzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxyz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxxxyz_0[i] * pb_x + g_0_xzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxzz_0[i] = g_0_zzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxxxzz_0[i] * pb_x + g_0_xzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyyy_0[i] = 2.0 * g_0_xxz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxyyy_0[i] * pb_z +
                                  g_0_xxzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxyyz_0[i] = g_0_zzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyyz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxxyyz_0[i] * pb_x + g_0_xzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyzz_0[i] = g_0_zzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxxyzz_0[i] * pb_x + g_0_xzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxzzz_0[i] = g_0_zzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxxzzz_0[i] * pb_x + g_0_xzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyyy_0[i] = 2.0 * g_0_xxz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxyyyy_0[i] * pb_z +
                                  g_0_xxzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxyyyz_0[i] = g_0_zzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyyz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxyyyz_0[i] * pb_x + g_0_xzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyzz_0[i] = g_0_zzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxyyzz_0[i] * pb_x + g_0_xzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyzzz_0[i] = g_0_zzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxyzzz_0[i] * pb_x + g_0_xzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxzzzz_0[i] = g_0_zzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xzzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xxzzzz_0[i] * pb_x + g_0_xzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyyy_0[i] = 2.0 * g_0_xxz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xyyyyy_0[i] * pb_z +
                                  g_0_xxzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xyyyyz_0[i] = g_0_zzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xyyyyz_0[i] * pb_x + g_0_xzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyzz_0[i] = g_0_zzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xyyyzz_0[i] * pb_x + g_0_xzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyzzz_0[i] = g_0_zzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xyyzzz_0[i] * pb_x + g_0_xzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyzzzz_0[i] = g_0_zzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xyzzzz_0[i] * pb_x + g_0_xzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xzzzzz_0[i] = g_0_zzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzzz_1[i] * fi_abcd_0 +
                                  g_0_xzzz_0_xzzzzz_0[i] * pb_x + g_0_xzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyy_0[i] =
            g_0_zzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyy_0[i] * pb_x + g_0_xzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyz_0[i] =
            g_0_zzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyz_0[i] * pb_x + g_0_xzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyzz_0[i] =
            g_0_zzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyzz_0[i] * pb_x + g_0_xzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyzzz_0[i] =
            g_0_zzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyzzz_0[i] * pb_x + g_0_xzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyzzzz_0[i] =
            g_0_zzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzzzz_0[i] * pb_x + g_0_xzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yzzzzz_0[i] =
            g_0_zzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzzzz_0[i] * pb_x + g_0_xzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_zzzzzz_0[i] =
            g_0_zzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzzzz_0[i] * pb_x + g_0_xzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 280-308 components of targeted buffer : SHSI

    auto g_0_xyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 280);

    auto g_0_xyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 281);

    auto g_0_xyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 282);

    auto g_0_xyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 283);

    auto g_0_xyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 284);

    auto g_0_xyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 285);

    auto g_0_xyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 286);

    auto g_0_xyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 287);

    auto g_0_xyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 288);

    auto g_0_xyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 289);

    auto g_0_xyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 290);

    auto g_0_xyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 291);

    auto g_0_xyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 292);

    auto g_0_xyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 293);

    auto g_0_xyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 294);

    auto g_0_xyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 295);

    auto g_0_xyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 296);

    auto g_0_xyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 297);

    auto g_0_xyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 298);

    auto g_0_xyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 299);

    auto g_0_xyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 300);

    auto g_0_xyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 301);

    auto g_0_xyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 302);

    auto g_0_xyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 303);

    auto g_0_xyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 304);

    auto g_0_xyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 305);

    auto g_0_xyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 306);

    auto g_0_xyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 307);

#pragma omp simd aligned(g_0_xyyyy_0_xxxxxx_0,     \
                             g_0_xyyyy_0_xxxxxy_0, \
                             g_0_xyyyy_0_xxxxxz_0, \
                             g_0_xyyyy_0_xxxxyy_0, \
                             g_0_xyyyy_0_xxxxyz_0, \
                             g_0_xyyyy_0_xxxxzz_0, \
                             g_0_xyyyy_0_xxxyyy_0, \
                             g_0_xyyyy_0_xxxyyz_0, \
                             g_0_xyyyy_0_xxxyzz_0, \
                             g_0_xyyyy_0_xxxzzz_0, \
                             g_0_xyyyy_0_xxyyyy_0, \
                             g_0_xyyyy_0_xxyyyz_0, \
                             g_0_xyyyy_0_xxyyzz_0, \
                             g_0_xyyyy_0_xxyzzz_0, \
                             g_0_xyyyy_0_xxzzzz_0, \
                             g_0_xyyyy_0_xyyyyy_0, \
                             g_0_xyyyy_0_xyyyyz_0, \
                             g_0_xyyyy_0_xyyyzz_0, \
                             g_0_xyyyy_0_xyyzzz_0, \
                             g_0_xyyyy_0_xyzzzz_0, \
                             g_0_xyyyy_0_xzzzzz_0, \
                             g_0_xyyyy_0_yyyyyy_0, \
                             g_0_xyyyy_0_yyyyyz_0, \
                             g_0_xyyyy_0_yyyyzz_0, \
                             g_0_xyyyy_0_yyyzzz_0, \
                             g_0_xyyyy_0_yyzzzz_0, \
                             g_0_xyyyy_0_yzzzzz_0, \
                             g_0_xyyyy_0_zzzzzz_0, \
                             g_0_yyyy_0_xxxxx_1,   \
                             g_0_yyyy_0_xxxxxx_0,  \
                             g_0_yyyy_0_xxxxxx_1,  \
                             g_0_yyyy_0_xxxxxy_0,  \
                             g_0_yyyy_0_xxxxxy_1,  \
                             g_0_yyyy_0_xxxxxz_0,  \
                             g_0_yyyy_0_xxxxxz_1,  \
                             g_0_yyyy_0_xxxxy_1,   \
                             g_0_yyyy_0_xxxxyy_0,  \
                             g_0_yyyy_0_xxxxyy_1,  \
                             g_0_yyyy_0_xxxxyz_0,  \
                             g_0_yyyy_0_xxxxyz_1,  \
                             g_0_yyyy_0_xxxxz_1,   \
                             g_0_yyyy_0_xxxxzz_0,  \
                             g_0_yyyy_0_xxxxzz_1,  \
                             g_0_yyyy_0_xxxyy_1,   \
                             g_0_yyyy_0_xxxyyy_0,  \
                             g_0_yyyy_0_xxxyyy_1,  \
                             g_0_yyyy_0_xxxyyz_0,  \
                             g_0_yyyy_0_xxxyyz_1,  \
                             g_0_yyyy_0_xxxyz_1,   \
                             g_0_yyyy_0_xxxyzz_0,  \
                             g_0_yyyy_0_xxxyzz_1,  \
                             g_0_yyyy_0_xxxzz_1,   \
                             g_0_yyyy_0_xxxzzz_0,  \
                             g_0_yyyy_0_xxxzzz_1,  \
                             g_0_yyyy_0_xxyyy_1,   \
                             g_0_yyyy_0_xxyyyy_0,  \
                             g_0_yyyy_0_xxyyyy_1,  \
                             g_0_yyyy_0_xxyyyz_0,  \
                             g_0_yyyy_0_xxyyyz_1,  \
                             g_0_yyyy_0_xxyyz_1,   \
                             g_0_yyyy_0_xxyyzz_0,  \
                             g_0_yyyy_0_xxyyzz_1,  \
                             g_0_yyyy_0_xxyzz_1,   \
                             g_0_yyyy_0_xxyzzz_0,  \
                             g_0_yyyy_0_xxyzzz_1,  \
                             g_0_yyyy_0_xxzzz_1,   \
                             g_0_yyyy_0_xxzzzz_0,  \
                             g_0_yyyy_0_xxzzzz_1,  \
                             g_0_yyyy_0_xyyyy_1,   \
                             g_0_yyyy_0_xyyyyy_0,  \
                             g_0_yyyy_0_xyyyyy_1,  \
                             g_0_yyyy_0_xyyyyz_0,  \
                             g_0_yyyy_0_xyyyyz_1,  \
                             g_0_yyyy_0_xyyyz_1,   \
                             g_0_yyyy_0_xyyyzz_0,  \
                             g_0_yyyy_0_xyyyzz_1,  \
                             g_0_yyyy_0_xyyzz_1,   \
                             g_0_yyyy_0_xyyzzz_0,  \
                             g_0_yyyy_0_xyyzzz_1,  \
                             g_0_yyyy_0_xyzzz_1,   \
                             g_0_yyyy_0_xyzzzz_0,  \
                             g_0_yyyy_0_xyzzzz_1,  \
                             g_0_yyyy_0_xzzzz_1,   \
                             g_0_yyyy_0_xzzzzz_0,  \
                             g_0_yyyy_0_xzzzzz_1,  \
                             g_0_yyyy_0_yyyyy_1,   \
                             g_0_yyyy_0_yyyyyy_0,  \
                             g_0_yyyy_0_yyyyyy_1,  \
                             g_0_yyyy_0_yyyyyz_0,  \
                             g_0_yyyy_0_yyyyyz_1,  \
                             g_0_yyyy_0_yyyyz_1,   \
                             g_0_yyyy_0_yyyyzz_0,  \
                             g_0_yyyy_0_yyyyzz_1,  \
                             g_0_yyyy_0_yyyzz_1,   \
                             g_0_yyyy_0_yyyzzz_0,  \
                             g_0_yyyy_0_yyyzzz_1,  \
                             g_0_yyyy_0_yyzzz_1,   \
                             g_0_yyyy_0_yyzzzz_0,  \
                             g_0_yyyy_0_yyzzzz_1,  \
                             g_0_yyyy_0_yzzzz_1,   \
                             g_0_yyyy_0_yzzzzz_0,  \
                             g_0_yyyy_0_yzzzzz_1,  \
                             g_0_yyyy_0_zzzzz_1,   \
                             g_0_yyyy_0_zzzzzz_0,  \
                             g_0_yyyy_0_zzzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyy_0_xxxxxx_0[i] = 6.0 * g_0_yyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxx_0[i] * pb_x + g_0_yyyy_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxy_0[i] = 5.0 * g_0_yyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxy_0[i] * pb_x + g_0_yyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxz_0[i] = 5.0 * g_0_yyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxz_0[i] * pb_x + g_0_yyyy_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyy_0[i] = 4.0 * g_0_yyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyy_0[i] * pb_x + g_0_yyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyz_0[i] = 4.0 * g_0_yyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyz_0[i] * pb_x + g_0_yyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxzz_0[i] = 4.0 * g_0_yyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxzz_0[i] * pb_x + g_0_yyyy_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyy_0[i] = 3.0 * g_0_yyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyy_0[i] * pb_x + g_0_yyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyz_0[i] = 3.0 * g_0_yyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyz_0[i] * pb_x + g_0_yyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyzz_0[i] = 3.0 * g_0_yyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyzz_0[i] * pb_x + g_0_yyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxzzz_0[i] = 3.0 * g_0_yyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxzzz_0[i] * pb_x + g_0_yyyy_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyy_0[i] = 2.0 * g_0_yyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyy_0[i] * pb_x + g_0_yyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyz_0[i] = 2.0 * g_0_yyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyz_0[i] * pb_x + g_0_yyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyzz_0[i] = 2.0 * g_0_yyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzz_0[i] * pb_x + g_0_yyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyzzz_0[i] = 2.0 * g_0_yyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzzz_0[i] * pb_x + g_0_yyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxzzzz_0[i] = 2.0 * g_0_yyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzzzz_0[i] * pb_x + g_0_yyyy_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyy_0[i] = g_0_yyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyy_0[i] * pb_x + g_0_yyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyz_0[i] = g_0_yyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyz_0[i] * pb_x + g_0_yyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyzz_0[i] = g_0_yyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzz_0[i] * pb_x + g_0_yyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyzzz_0[i] = g_0_yyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzz_0[i] * pb_x + g_0_yyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyzzzz_0[i] = g_0_yyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzzz_0[i] * pb_x + g_0_yyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xzzzzz_0[i] = g_0_yyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzzzz_0[i] * pb_x + g_0_yyyy_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyy_0[i] = g_0_yyyy_0_yyyyyy_0[i] * pb_x + g_0_yyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyz_0[i] = g_0_yyyy_0_yyyyyz_0[i] * pb_x + g_0_yyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyzz_0[i] = g_0_yyyy_0_yyyyzz_0[i] * pb_x + g_0_yyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyzzz_0[i] = g_0_yyyy_0_yyyzzz_0[i] * pb_x + g_0_yyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyzzzz_0[i] = g_0_yyyy_0_yyzzzz_0[i] * pb_x + g_0_yyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yzzzzz_0[i] = g_0_yyyy_0_yzzzzz_0[i] * pb_x + g_0_yyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_zzzzzz_0[i] = g_0_yyyy_0_zzzzzz_0[i] * pb_x + g_0_yyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 308-336 components of targeted buffer : SHSI

    auto g_0_xyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 308);

    auto g_0_xyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 309);

    auto g_0_xyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 310);

    auto g_0_xyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 311);

    auto g_0_xyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 312);

    auto g_0_xyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 313);

    auto g_0_xyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 314);

    auto g_0_xyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 315);

    auto g_0_xyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 316);

    auto g_0_xyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 317);

    auto g_0_xyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 318);

    auto g_0_xyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 319);

    auto g_0_xyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 320);

    auto g_0_xyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 321);

    auto g_0_xyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 322);

    auto g_0_xyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 323);

    auto g_0_xyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 324);

    auto g_0_xyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 325);

    auto g_0_xyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 326);

    auto g_0_xyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 327);

    auto g_0_xyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 328);

    auto g_0_xyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 329);

    auto g_0_xyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 330);

    auto g_0_xyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 331);

    auto g_0_xyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 332);

    auto g_0_xyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 333);

    auto g_0_xyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 334);

    auto g_0_xyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 335);

#pragma omp simd aligned(g_0_xyyy_0_xxxxxx_0,      \
                             g_0_xyyy_0_xxxxxx_1,  \
                             g_0_xyyy_0_xxxxxy_0,  \
                             g_0_xyyy_0_xxxxxy_1,  \
                             g_0_xyyy_0_xxxxyy_0,  \
                             g_0_xyyy_0_xxxxyy_1,  \
                             g_0_xyyy_0_xxxyyy_0,  \
                             g_0_xyyy_0_xxxyyy_1,  \
                             g_0_xyyy_0_xxyyyy_0,  \
                             g_0_xyyy_0_xxyyyy_1,  \
                             g_0_xyyy_0_xyyyyy_0,  \
                             g_0_xyyy_0_xyyyyy_1,  \
                             g_0_xyyyz_0_xxxxxx_0, \
                             g_0_xyyyz_0_xxxxxy_0, \
                             g_0_xyyyz_0_xxxxxz_0, \
                             g_0_xyyyz_0_xxxxyy_0, \
                             g_0_xyyyz_0_xxxxyz_0, \
                             g_0_xyyyz_0_xxxxzz_0, \
                             g_0_xyyyz_0_xxxyyy_0, \
                             g_0_xyyyz_0_xxxyyz_0, \
                             g_0_xyyyz_0_xxxyzz_0, \
                             g_0_xyyyz_0_xxxzzz_0, \
                             g_0_xyyyz_0_xxyyyy_0, \
                             g_0_xyyyz_0_xxyyyz_0, \
                             g_0_xyyyz_0_xxyyzz_0, \
                             g_0_xyyyz_0_xxyzzz_0, \
                             g_0_xyyyz_0_xxzzzz_0, \
                             g_0_xyyyz_0_xyyyyy_0, \
                             g_0_xyyyz_0_xyyyyz_0, \
                             g_0_xyyyz_0_xyyyzz_0, \
                             g_0_xyyyz_0_xyyzzz_0, \
                             g_0_xyyyz_0_xyzzzz_0, \
                             g_0_xyyyz_0_xzzzzz_0, \
                             g_0_xyyyz_0_yyyyyy_0, \
                             g_0_xyyyz_0_yyyyyz_0, \
                             g_0_xyyyz_0_yyyyzz_0, \
                             g_0_xyyyz_0_yyyzzz_0, \
                             g_0_xyyyz_0_yyzzzz_0, \
                             g_0_xyyyz_0_yzzzzz_0, \
                             g_0_xyyyz_0_zzzzzz_0, \
                             g_0_yyyz_0_xxxxxz_0,  \
                             g_0_yyyz_0_xxxxxz_1,  \
                             g_0_yyyz_0_xxxxyz_0,  \
                             g_0_yyyz_0_xxxxyz_1,  \
                             g_0_yyyz_0_xxxxz_1,   \
                             g_0_yyyz_0_xxxxzz_0,  \
                             g_0_yyyz_0_xxxxzz_1,  \
                             g_0_yyyz_0_xxxyyz_0,  \
                             g_0_yyyz_0_xxxyyz_1,  \
                             g_0_yyyz_0_xxxyz_1,   \
                             g_0_yyyz_0_xxxyzz_0,  \
                             g_0_yyyz_0_xxxyzz_1,  \
                             g_0_yyyz_0_xxxzz_1,   \
                             g_0_yyyz_0_xxxzzz_0,  \
                             g_0_yyyz_0_xxxzzz_1,  \
                             g_0_yyyz_0_xxyyyz_0,  \
                             g_0_yyyz_0_xxyyyz_1,  \
                             g_0_yyyz_0_xxyyz_1,   \
                             g_0_yyyz_0_xxyyzz_0,  \
                             g_0_yyyz_0_xxyyzz_1,  \
                             g_0_yyyz_0_xxyzz_1,   \
                             g_0_yyyz_0_xxyzzz_0,  \
                             g_0_yyyz_0_xxyzzz_1,  \
                             g_0_yyyz_0_xxzzz_1,   \
                             g_0_yyyz_0_xxzzzz_0,  \
                             g_0_yyyz_0_xxzzzz_1,  \
                             g_0_yyyz_0_xyyyyz_0,  \
                             g_0_yyyz_0_xyyyyz_1,  \
                             g_0_yyyz_0_xyyyz_1,   \
                             g_0_yyyz_0_xyyyzz_0,  \
                             g_0_yyyz_0_xyyyzz_1,  \
                             g_0_yyyz_0_xyyzz_1,   \
                             g_0_yyyz_0_xyyzzz_0,  \
                             g_0_yyyz_0_xyyzzz_1,  \
                             g_0_yyyz_0_xyzzz_1,   \
                             g_0_yyyz_0_xyzzzz_0,  \
                             g_0_yyyz_0_xyzzzz_1,  \
                             g_0_yyyz_0_xzzzz_1,   \
                             g_0_yyyz_0_xzzzzz_0,  \
                             g_0_yyyz_0_xzzzzz_1,  \
                             g_0_yyyz_0_yyyyyy_0,  \
                             g_0_yyyz_0_yyyyyy_1,  \
                             g_0_yyyz_0_yyyyyz_0,  \
                             g_0_yyyz_0_yyyyyz_1,  \
                             g_0_yyyz_0_yyyyz_1,   \
                             g_0_yyyz_0_yyyyzz_0,  \
                             g_0_yyyz_0_yyyyzz_1,  \
                             g_0_yyyz_0_yyyzz_1,   \
                             g_0_yyyz_0_yyyzzz_0,  \
                             g_0_yyyz_0_yyyzzz_1,  \
                             g_0_yyyz_0_yyzzz_1,   \
                             g_0_yyyz_0_yyzzzz_0,  \
                             g_0_yyyz_0_yyzzzz_1,  \
                             g_0_yyyz_0_yzzzz_1,   \
                             g_0_yyyz_0_yzzzzz_0,  \
                             g_0_yyyz_0_yzzzzz_1,  \
                             g_0_yyyz_0_zzzzz_1,   \
                             g_0_yyyz_0_zzzzzz_0,  \
                             g_0_yyyz_0_zzzzzz_1,  \
                             wp_x,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyz_0_xxxxxx_0[i] = g_0_xyyy_0_xxxxxx_0[i] * pb_z + g_0_xyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxy_0[i] = g_0_xyyy_0_xxxxxy_0[i] * pb_z + g_0_xyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxz_0[i] = 5.0 * g_0_yyyz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxz_0[i] * pb_x + g_0_yyyz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxyy_0[i] = g_0_xyyy_0_xxxxyy_0[i] * pb_z + g_0_xyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxyz_0[i] = 4.0 * g_0_yyyz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxyz_0[i] * pb_x + g_0_yyyz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxzz_0[i] = 4.0 * g_0_yyyz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxzz_0[i] * pb_x + g_0_yyyz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyyy_0[i] = g_0_xyyy_0_xxxyyy_0[i] * pb_z + g_0_xyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxyyz_0[i] = 3.0 * g_0_yyyz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyyz_0[i] * pb_x + g_0_yyyz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyzz_0[i] = 3.0 * g_0_yyyz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyzz_0[i] * pb_x + g_0_yyyz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxzzz_0[i] = 3.0 * g_0_yyyz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxzzz_0[i] * pb_x + g_0_yyyz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyyy_0[i] = g_0_xyyy_0_xxyyyy_0[i] * pb_z + g_0_xyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxyyyz_0[i] = 2.0 * g_0_yyyz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyyz_0[i] * pb_x + g_0_yyyz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyzz_0[i] = 2.0 * g_0_yyyz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyzz_0[i] * pb_x + g_0_yyyz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyzzz_0[i] = 2.0 * g_0_yyyz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyzzz_0[i] * pb_x + g_0_yyyz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxzzzz_0[i] = 2.0 * g_0_yyyz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxzzzz_0[i] * pb_x + g_0_yyyz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyyy_0[i] = g_0_xyyy_0_xyyyyy_0[i] * pb_z + g_0_xyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xyyyyz_0[i] = g_0_yyyz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyyz_0[i] * pb_x + g_0_yyyz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyzz_0[i] = g_0_yyyz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyzz_0[i] * pb_x + g_0_yyyz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyzzz_0[i] = g_0_yyyz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyzzz_0[i] * pb_x + g_0_yyyz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyzzzz_0[i] = g_0_yyyz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyzzzz_0[i] * pb_x + g_0_yyyz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xzzzzz_0[i] = g_0_yyyz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xzzzzz_0[i] * pb_x + g_0_yyyz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyy_0[i] = g_0_yyyz_0_yyyyyy_0[i] * pb_x + g_0_yyyz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyz_0[i] = g_0_yyyz_0_yyyyyz_0[i] * pb_x + g_0_yyyz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyzz_0[i] = g_0_yyyz_0_yyyyzz_0[i] * pb_x + g_0_yyyz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyzzz_0[i] = g_0_yyyz_0_yyyzzz_0[i] * pb_x + g_0_yyyz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyzzzz_0[i] = g_0_yyyz_0_yyzzzz_0[i] * pb_x + g_0_yyyz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yzzzzz_0[i] = g_0_yyyz_0_yzzzzz_0[i] * pb_x + g_0_yyyz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_zzzzzz_0[i] = g_0_yyyz_0_zzzzzz_0[i] * pb_x + g_0_yyyz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 336-364 components of targeted buffer : SHSI

    auto g_0_xyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 336);

    auto g_0_xyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 337);

    auto g_0_xyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 338);

    auto g_0_xyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 339);

    auto g_0_xyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 340);

    auto g_0_xyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 341);

    auto g_0_xyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 342);

    auto g_0_xyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 343);

    auto g_0_xyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 344);

    auto g_0_xyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 345);

    auto g_0_xyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 346);

    auto g_0_xyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 347);

    auto g_0_xyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 348);

    auto g_0_xyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 349);

    auto g_0_xyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 350);

    auto g_0_xyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 351);

    auto g_0_xyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 352);

    auto g_0_xyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 353);

    auto g_0_xyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 354);

    auto g_0_xyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 355);

    auto g_0_xyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 356);

    auto g_0_xyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 357);

    auto g_0_xyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 358);

    auto g_0_xyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 359);

    auto g_0_xyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 360);

    auto g_0_xyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 361);

    auto g_0_xyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 362);

    auto g_0_xyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 363);

#pragma omp simd aligned(g_0_xyyzz_0_xxxxxx_0,     \
                             g_0_xyyzz_0_xxxxxy_0, \
                             g_0_xyyzz_0_xxxxxz_0, \
                             g_0_xyyzz_0_xxxxyy_0, \
                             g_0_xyyzz_0_xxxxyz_0, \
                             g_0_xyyzz_0_xxxxzz_0, \
                             g_0_xyyzz_0_xxxyyy_0, \
                             g_0_xyyzz_0_xxxyyz_0, \
                             g_0_xyyzz_0_xxxyzz_0, \
                             g_0_xyyzz_0_xxxzzz_0, \
                             g_0_xyyzz_0_xxyyyy_0, \
                             g_0_xyyzz_0_xxyyyz_0, \
                             g_0_xyyzz_0_xxyyzz_0, \
                             g_0_xyyzz_0_xxyzzz_0, \
                             g_0_xyyzz_0_xxzzzz_0, \
                             g_0_xyyzz_0_xyyyyy_0, \
                             g_0_xyyzz_0_xyyyyz_0, \
                             g_0_xyyzz_0_xyyyzz_0, \
                             g_0_xyyzz_0_xyyzzz_0, \
                             g_0_xyyzz_0_xyzzzz_0, \
                             g_0_xyyzz_0_xzzzzz_0, \
                             g_0_xyyzz_0_yyyyyy_0, \
                             g_0_xyyzz_0_yyyyyz_0, \
                             g_0_xyyzz_0_yyyyzz_0, \
                             g_0_xyyzz_0_yyyzzz_0, \
                             g_0_xyyzz_0_yyzzzz_0, \
                             g_0_xyyzz_0_yzzzzz_0, \
                             g_0_xyyzz_0_zzzzzz_0, \
                             g_0_yyzz_0_xxxxx_1,   \
                             g_0_yyzz_0_xxxxxx_0,  \
                             g_0_yyzz_0_xxxxxx_1,  \
                             g_0_yyzz_0_xxxxxy_0,  \
                             g_0_yyzz_0_xxxxxy_1,  \
                             g_0_yyzz_0_xxxxxz_0,  \
                             g_0_yyzz_0_xxxxxz_1,  \
                             g_0_yyzz_0_xxxxy_1,   \
                             g_0_yyzz_0_xxxxyy_0,  \
                             g_0_yyzz_0_xxxxyy_1,  \
                             g_0_yyzz_0_xxxxyz_0,  \
                             g_0_yyzz_0_xxxxyz_1,  \
                             g_0_yyzz_0_xxxxz_1,   \
                             g_0_yyzz_0_xxxxzz_0,  \
                             g_0_yyzz_0_xxxxzz_1,  \
                             g_0_yyzz_0_xxxyy_1,   \
                             g_0_yyzz_0_xxxyyy_0,  \
                             g_0_yyzz_0_xxxyyy_1,  \
                             g_0_yyzz_0_xxxyyz_0,  \
                             g_0_yyzz_0_xxxyyz_1,  \
                             g_0_yyzz_0_xxxyz_1,   \
                             g_0_yyzz_0_xxxyzz_0,  \
                             g_0_yyzz_0_xxxyzz_1,  \
                             g_0_yyzz_0_xxxzz_1,   \
                             g_0_yyzz_0_xxxzzz_0,  \
                             g_0_yyzz_0_xxxzzz_1,  \
                             g_0_yyzz_0_xxyyy_1,   \
                             g_0_yyzz_0_xxyyyy_0,  \
                             g_0_yyzz_0_xxyyyy_1,  \
                             g_0_yyzz_0_xxyyyz_0,  \
                             g_0_yyzz_0_xxyyyz_1,  \
                             g_0_yyzz_0_xxyyz_1,   \
                             g_0_yyzz_0_xxyyzz_0,  \
                             g_0_yyzz_0_xxyyzz_1,  \
                             g_0_yyzz_0_xxyzz_1,   \
                             g_0_yyzz_0_xxyzzz_0,  \
                             g_0_yyzz_0_xxyzzz_1,  \
                             g_0_yyzz_0_xxzzz_1,   \
                             g_0_yyzz_0_xxzzzz_0,  \
                             g_0_yyzz_0_xxzzzz_1,  \
                             g_0_yyzz_0_xyyyy_1,   \
                             g_0_yyzz_0_xyyyyy_0,  \
                             g_0_yyzz_0_xyyyyy_1,  \
                             g_0_yyzz_0_xyyyyz_0,  \
                             g_0_yyzz_0_xyyyyz_1,  \
                             g_0_yyzz_0_xyyyz_1,   \
                             g_0_yyzz_0_xyyyzz_0,  \
                             g_0_yyzz_0_xyyyzz_1,  \
                             g_0_yyzz_0_xyyzz_1,   \
                             g_0_yyzz_0_xyyzzz_0,  \
                             g_0_yyzz_0_xyyzzz_1,  \
                             g_0_yyzz_0_xyzzz_1,   \
                             g_0_yyzz_0_xyzzzz_0,  \
                             g_0_yyzz_0_xyzzzz_1,  \
                             g_0_yyzz_0_xzzzz_1,   \
                             g_0_yyzz_0_xzzzzz_0,  \
                             g_0_yyzz_0_xzzzzz_1,  \
                             g_0_yyzz_0_yyyyy_1,   \
                             g_0_yyzz_0_yyyyyy_0,  \
                             g_0_yyzz_0_yyyyyy_1,  \
                             g_0_yyzz_0_yyyyyz_0,  \
                             g_0_yyzz_0_yyyyyz_1,  \
                             g_0_yyzz_0_yyyyz_1,   \
                             g_0_yyzz_0_yyyyzz_0,  \
                             g_0_yyzz_0_yyyyzz_1,  \
                             g_0_yyzz_0_yyyzz_1,   \
                             g_0_yyzz_0_yyyzzz_0,  \
                             g_0_yyzz_0_yyyzzz_1,  \
                             g_0_yyzz_0_yyzzz_1,   \
                             g_0_yyzz_0_yyzzzz_0,  \
                             g_0_yyzz_0_yyzzzz_1,  \
                             g_0_yyzz_0_yzzzz_1,   \
                             g_0_yyzz_0_yzzzzz_0,  \
                             g_0_yyzz_0_yzzzzz_1,  \
                             g_0_yyzz_0_zzzzz_1,   \
                             g_0_yyzz_0_zzzzzz_0,  \
                             g_0_yyzz_0_zzzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzz_0_xxxxxx_0[i] = 6.0 * g_0_yyzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxx_0[i] * pb_x + g_0_yyzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxy_0[i] = 5.0 * g_0_yyzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxy_0[i] * pb_x + g_0_yyzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxz_0[i] = 5.0 * g_0_yyzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxz_0[i] * pb_x + g_0_yyzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyy_0[i] = 4.0 * g_0_yyzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyy_0[i] * pb_x + g_0_yyzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyz_0[i] = 4.0 * g_0_yyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyz_0[i] * pb_x + g_0_yyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxzz_0[i] = 4.0 * g_0_yyzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxzz_0[i] * pb_x + g_0_yyzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyy_0[i] = 3.0 * g_0_yyzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyy_0[i] * pb_x + g_0_yyzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyz_0[i] = 3.0 * g_0_yyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyz_0[i] * pb_x + g_0_yyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyzz_0[i] = 3.0 * g_0_yyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyzz_0[i] * pb_x + g_0_yyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxzzz_0[i] = 3.0 * g_0_yyzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxzzz_0[i] * pb_x + g_0_yyzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyy_0[i] = 2.0 * g_0_yyzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyy_0[i] * pb_x + g_0_yyzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyz_0[i] = 2.0 * g_0_yyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyz_0[i] * pb_x + g_0_yyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyzz_0[i] = 2.0 * g_0_yyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyzz_0[i] * pb_x + g_0_yyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyzzz_0[i] = 2.0 * g_0_yyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyzzz_0[i] * pb_x + g_0_yyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxzzzz_0[i] = 2.0 * g_0_yyzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxzzzz_0[i] * pb_x + g_0_yyzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyy_0[i] = g_0_yyzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyy_0[i] * pb_x + g_0_yyzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyz_0[i] = g_0_yyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyz_0[i] * pb_x + g_0_yyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyzz_0[i] = g_0_yyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyzz_0[i] * pb_x + g_0_yyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyzzz_0[i] = g_0_yyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyzzz_0[i] * pb_x + g_0_yyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyzzzz_0[i] = g_0_yyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyzzzz_0[i] * pb_x + g_0_yyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xzzzzz_0[i] = g_0_yyzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xzzzzz_0[i] * pb_x + g_0_yyzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyy_0[i] = g_0_yyzz_0_yyyyyy_0[i] * pb_x + g_0_yyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyz_0[i] = g_0_yyzz_0_yyyyyz_0[i] * pb_x + g_0_yyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyzz_0[i] = g_0_yyzz_0_yyyyzz_0[i] * pb_x + g_0_yyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyzzz_0[i] = g_0_yyzz_0_yyyzzz_0[i] * pb_x + g_0_yyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyzzzz_0[i] = g_0_yyzz_0_yyzzzz_0[i] * pb_x + g_0_yyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yzzzzz_0[i] = g_0_yyzz_0_yzzzzz_0[i] * pb_x + g_0_yyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_zzzzzz_0[i] = g_0_yyzz_0_zzzzzz_0[i] * pb_x + g_0_yyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 364-392 components of targeted buffer : SHSI

    auto g_0_xyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 364);

    auto g_0_xyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 365);

    auto g_0_xyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 366);

    auto g_0_xyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 367);

    auto g_0_xyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 368);

    auto g_0_xyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 369);

    auto g_0_xyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 370);

    auto g_0_xyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 371);

    auto g_0_xyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 372);

    auto g_0_xyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 373);

    auto g_0_xyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 374);

    auto g_0_xyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 375);

    auto g_0_xyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 376);

    auto g_0_xyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 377);

    auto g_0_xyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 378);

    auto g_0_xyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 379);

    auto g_0_xyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 380);

    auto g_0_xyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 381);

    auto g_0_xyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 382);

    auto g_0_xyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 383);

    auto g_0_xyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 384);

    auto g_0_xyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 385);

    auto g_0_xyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 386);

    auto g_0_xyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 387);

    auto g_0_xyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 388);

    auto g_0_xyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 389);

    auto g_0_xyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 390);

    auto g_0_xyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 391);

#pragma omp simd aligned(g_0_xyzzz_0_xxxxxx_0,     \
                             g_0_xyzzz_0_xxxxxy_0, \
                             g_0_xyzzz_0_xxxxxz_0, \
                             g_0_xyzzz_0_xxxxyy_0, \
                             g_0_xyzzz_0_xxxxyz_0, \
                             g_0_xyzzz_0_xxxxzz_0, \
                             g_0_xyzzz_0_xxxyyy_0, \
                             g_0_xyzzz_0_xxxyyz_0, \
                             g_0_xyzzz_0_xxxyzz_0, \
                             g_0_xyzzz_0_xxxzzz_0, \
                             g_0_xyzzz_0_xxyyyy_0, \
                             g_0_xyzzz_0_xxyyyz_0, \
                             g_0_xyzzz_0_xxyyzz_0, \
                             g_0_xyzzz_0_xxyzzz_0, \
                             g_0_xyzzz_0_xxzzzz_0, \
                             g_0_xyzzz_0_xyyyyy_0, \
                             g_0_xyzzz_0_xyyyyz_0, \
                             g_0_xyzzz_0_xyyyzz_0, \
                             g_0_xyzzz_0_xyyzzz_0, \
                             g_0_xyzzz_0_xyzzzz_0, \
                             g_0_xyzzz_0_xzzzzz_0, \
                             g_0_xyzzz_0_yyyyyy_0, \
                             g_0_xyzzz_0_yyyyyz_0, \
                             g_0_xyzzz_0_yyyyzz_0, \
                             g_0_xyzzz_0_yyyzzz_0, \
                             g_0_xyzzz_0_yyzzzz_0, \
                             g_0_xyzzz_0_yzzzzz_0, \
                             g_0_xyzzz_0_zzzzzz_0, \
                             g_0_xzzz_0_xxxxxx_0,  \
                             g_0_xzzz_0_xxxxxx_1,  \
                             g_0_xzzz_0_xxxxxz_0,  \
                             g_0_xzzz_0_xxxxxz_1,  \
                             g_0_xzzz_0_xxxxzz_0,  \
                             g_0_xzzz_0_xxxxzz_1,  \
                             g_0_xzzz_0_xxxzzz_0,  \
                             g_0_xzzz_0_xxxzzz_1,  \
                             g_0_xzzz_0_xxzzzz_0,  \
                             g_0_xzzz_0_xxzzzz_1,  \
                             g_0_xzzz_0_xzzzzz_0,  \
                             g_0_xzzz_0_xzzzzz_1,  \
                             g_0_yzzz_0_xxxxxy_0,  \
                             g_0_yzzz_0_xxxxxy_1,  \
                             g_0_yzzz_0_xxxxy_1,   \
                             g_0_yzzz_0_xxxxyy_0,  \
                             g_0_yzzz_0_xxxxyy_1,  \
                             g_0_yzzz_0_xxxxyz_0,  \
                             g_0_yzzz_0_xxxxyz_1,  \
                             g_0_yzzz_0_xxxyy_1,   \
                             g_0_yzzz_0_xxxyyy_0,  \
                             g_0_yzzz_0_xxxyyy_1,  \
                             g_0_yzzz_0_xxxyyz_0,  \
                             g_0_yzzz_0_xxxyyz_1,  \
                             g_0_yzzz_0_xxxyz_1,   \
                             g_0_yzzz_0_xxxyzz_0,  \
                             g_0_yzzz_0_xxxyzz_1,  \
                             g_0_yzzz_0_xxyyy_1,   \
                             g_0_yzzz_0_xxyyyy_0,  \
                             g_0_yzzz_0_xxyyyy_1,  \
                             g_0_yzzz_0_xxyyyz_0,  \
                             g_0_yzzz_0_xxyyyz_1,  \
                             g_0_yzzz_0_xxyyz_1,   \
                             g_0_yzzz_0_xxyyzz_0,  \
                             g_0_yzzz_0_xxyyzz_1,  \
                             g_0_yzzz_0_xxyzz_1,   \
                             g_0_yzzz_0_xxyzzz_0,  \
                             g_0_yzzz_0_xxyzzz_1,  \
                             g_0_yzzz_0_xyyyy_1,   \
                             g_0_yzzz_0_xyyyyy_0,  \
                             g_0_yzzz_0_xyyyyy_1,  \
                             g_0_yzzz_0_xyyyyz_0,  \
                             g_0_yzzz_0_xyyyyz_1,  \
                             g_0_yzzz_0_xyyyz_1,   \
                             g_0_yzzz_0_xyyyzz_0,  \
                             g_0_yzzz_0_xyyyzz_1,  \
                             g_0_yzzz_0_xyyzz_1,   \
                             g_0_yzzz_0_xyyzzz_0,  \
                             g_0_yzzz_0_xyyzzz_1,  \
                             g_0_yzzz_0_xyzzz_1,   \
                             g_0_yzzz_0_xyzzzz_0,  \
                             g_0_yzzz_0_xyzzzz_1,  \
                             g_0_yzzz_0_yyyyy_1,   \
                             g_0_yzzz_0_yyyyyy_0,  \
                             g_0_yzzz_0_yyyyyy_1,  \
                             g_0_yzzz_0_yyyyyz_0,  \
                             g_0_yzzz_0_yyyyyz_1,  \
                             g_0_yzzz_0_yyyyz_1,   \
                             g_0_yzzz_0_yyyyzz_0,  \
                             g_0_yzzz_0_yyyyzz_1,  \
                             g_0_yzzz_0_yyyzz_1,   \
                             g_0_yzzz_0_yyyzzz_0,  \
                             g_0_yzzz_0_yyyzzz_1,  \
                             g_0_yzzz_0_yyzzz_1,   \
                             g_0_yzzz_0_yyzzzz_0,  \
                             g_0_yzzz_0_yyzzzz_1,  \
                             g_0_yzzz_0_yzzzz_1,   \
                             g_0_yzzz_0_yzzzzz_0,  \
                             g_0_yzzz_0_yzzzzz_1,  \
                             g_0_yzzz_0_zzzzzz_0,  \
                             g_0_yzzz_0_zzzzzz_1,  \
                             wp_x,                 \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzz_0_xxxxxx_0[i] = g_0_xzzz_0_xxxxxx_0[i] * pb_y + g_0_xzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxxy_0[i] = 5.0 * g_0_yzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxy_0[i] * pb_x + g_0_yzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxz_0[i] = g_0_xzzz_0_xxxxxz_0[i] * pb_y + g_0_xzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxyy_0[i] = 4.0 * g_0_yzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyy_0[i] * pb_x + g_0_yzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxyz_0[i] = 4.0 * g_0_yzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyz_0[i] * pb_x + g_0_yzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxzz_0[i] = g_0_xzzz_0_xxxxzz_0[i] * pb_y + g_0_xzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxyyy_0[i] = 3.0 * g_0_yzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyy_0[i] * pb_x + g_0_yzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyyz_0[i] = 3.0 * g_0_yzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyz_0[i] * pb_x + g_0_yzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyzz_0[i] = 3.0 * g_0_yzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyzz_0[i] * pb_x + g_0_yzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxzzz_0[i] = g_0_xzzz_0_xxxzzz_0[i] * pb_y + g_0_xzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxyyyy_0[i] = 2.0 * g_0_yzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyy_0[i] * pb_x + g_0_yzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyyz_0[i] = 2.0 * g_0_yzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyz_0[i] * pb_x + g_0_yzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyzz_0[i] = 2.0 * g_0_yzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyzz_0[i] * pb_x + g_0_yzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyzzz_0[i] = 2.0 * g_0_yzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyzzz_0[i] * pb_x + g_0_yzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxzzzz_0[i] = g_0_xzzz_0_xxzzzz_0[i] * pb_y + g_0_xzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xyyyyy_0[i] = g_0_yzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyy_0[i] * pb_x + g_0_yzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyyz_0[i] = g_0_yzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyz_0[i] * pb_x + g_0_yzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyzz_0[i] = g_0_yzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyzz_0[i] * pb_x + g_0_yzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyzzz_0[i] = g_0_yzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyzzz_0[i] * pb_x + g_0_yzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyzzzz_0[i] = g_0_yzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyzzzz_0[i] * pb_x + g_0_yzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xzzzzz_0[i] = g_0_xzzz_0_xzzzzz_0[i] * pb_y + g_0_xzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_yyyyyy_0[i] = g_0_yzzz_0_yyyyyy_0[i] * pb_x + g_0_yzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyyz_0[i] = g_0_yzzz_0_yyyyyz_0[i] * pb_x + g_0_yzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyzz_0[i] = g_0_yzzz_0_yyyyzz_0[i] * pb_x + g_0_yzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyzzz_0[i] = g_0_yzzz_0_yyyzzz_0[i] * pb_x + g_0_yzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyzzzz_0[i] = g_0_yzzz_0_yyzzzz_0[i] * pb_x + g_0_yzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yzzzzz_0[i] = g_0_yzzz_0_yzzzzz_0[i] * pb_x + g_0_yzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_zzzzzz_0[i] = g_0_yzzz_0_zzzzzz_0[i] * pb_x + g_0_yzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 392-420 components of targeted buffer : SHSI

    auto g_0_xzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 392);

    auto g_0_xzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 393);

    auto g_0_xzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 394);

    auto g_0_xzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 395);

    auto g_0_xzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 396);

    auto g_0_xzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 397);

    auto g_0_xzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 398);

    auto g_0_xzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 399);

    auto g_0_xzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 400);

    auto g_0_xzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 401);

    auto g_0_xzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 402);

    auto g_0_xzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 403);

    auto g_0_xzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 404);

    auto g_0_xzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 405);

    auto g_0_xzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 406);

    auto g_0_xzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 407);

    auto g_0_xzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 408);

    auto g_0_xzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 409);

    auto g_0_xzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 410);

    auto g_0_xzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 411);

    auto g_0_xzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 412);

    auto g_0_xzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 413);

    auto g_0_xzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 414);

    auto g_0_xzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 415);

    auto g_0_xzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 416);

    auto g_0_xzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 417);

    auto g_0_xzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 418);

    auto g_0_xzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 419);

#pragma omp simd aligned(g_0_xzzzz_0_xxxxxx_0,     \
                             g_0_xzzzz_0_xxxxxy_0, \
                             g_0_xzzzz_0_xxxxxz_0, \
                             g_0_xzzzz_0_xxxxyy_0, \
                             g_0_xzzzz_0_xxxxyz_0, \
                             g_0_xzzzz_0_xxxxzz_0, \
                             g_0_xzzzz_0_xxxyyy_0, \
                             g_0_xzzzz_0_xxxyyz_0, \
                             g_0_xzzzz_0_xxxyzz_0, \
                             g_0_xzzzz_0_xxxzzz_0, \
                             g_0_xzzzz_0_xxyyyy_0, \
                             g_0_xzzzz_0_xxyyyz_0, \
                             g_0_xzzzz_0_xxyyzz_0, \
                             g_0_xzzzz_0_xxyzzz_0, \
                             g_0_xzzzz_0_xxzzzz_0, \
                             g_0_xzzzz_0_xyyyyy_0, \
                             g_0_xzzzz_0_xyyyyz_0, \
                             g_0_xzzzz_0_xyyyzz_0, \
                             g_0_xzzzz_0_xyyzzz_0, \
                             g_0_xzzzz_0_xyzzzz_0, \
                             g_0_xzzzz_0_xzzzzz_0, \
                             g_0_xzzzz_0_yyyyyy_0, \
                             g_0_xzzzz_0_yyyyyz_0, \
                             g_0_xzzzz_0_yyyyzz_0, \
                             g_0_xzzzz_0_yyyzzz_0, \
                             g_0_xzzzz_0_yyzzzz_0, \
                             g_0_xzzzz_0_yzzzzz_0, \
                             g_0_xzzzz_0_zzzzzz_0, \
                             g_0_zzzz_0_xxxxx_1,   \
                             g_0_zzzz_0_xxxxxx_0,  \
                             g_0_zzzz_0_xxxxxx_1,  \
                             g_0_zzzz_0_xxxxxy_0,  \
                             g_0_zzzz_0_xxxxxy_1,  \
                             g_0_zzzz_0_xxxxxz_0,  \
                             g_0_zzzz_0_xxxxxz_1,  \
                             g_0_zzzz_0_xxxxy_1,   \
                             g_0_zzzz_0_xxxxyy_0,  \
                             g_0_zzzz_0_xxxxyy_1,  \
                             g_0_zzzz_0_xxxxyz_0,  \
                             g_0_zzzz_0_xxxxyz_1,  \
                             g_0_zzzz_0_xxxxz_1,   \
                             g_0_zzzz_0_xxxxzz_0,  \
                             g_0_zzzz_0_xxxxzz_1,  \
                             g_0_zzzz_0_xxxyy_1,   \
                             g_0_zzzz_0_xxxyyy_0,  \
                             g_0_zzzz_0_xxxyyy_1,  \
                             g_0_zzzz_0_xxxyyz_0,  \
                             g_0_zzzz_0_xxxyyz_1,  \
                             g_0_zzzz_0_xxxyz_1,   \
                             g_0_zzzz_0_xxxyzz_0,  \
                             g_0_zzzz_0_xxxyzz_1,  \
                             g_0_zzzz_0_xxxzz_1,   \
                             g_0_zzzz_0_xxxzzz_0,  \
                             g_0_zzzz_0_xxxzzz_1,  \
                             g_0_zzzz_0_xxyyy_1,   \
                             g_0_zzzz_0_xxyyyy_0,  \
                             g_0_zzzz_0_xxyyyy_1,  \
                             g_0_zzzz_0_xxyyyz_0,  \
                             g_0_zzzz_0_xxyyyz_1,  \
                             g_0_zzzz_0_xxyyz_1,   \
                             g_0_zzzz_0_xxyyzz_0,  \
                             g_0_zzzz_0_xxyyzz_1,  \
                             g_0_zzzz_0_xxyzz_1,   \
                             g_0_zzzz_0_xxyzzz_0,  \
                             g_0_zzzz_0_xxyzzz_1,  \
                             g_0_zzzz_0_xxzzz_1,   \
                             g_0_zzzz_0_xxzzzz_0,  \
                             g_0_zzzz_0_xxzzzz_1,  \
                             g_0_zzzz_0_xyyyy_1,   \
                             g_0_zzzz_0_xyyyyy_0,  \
                             g_0_zzzz_0_xyyyyy_1,  \
                             g_0_zzzz_0_xyyyyz_0,  \
                             g_0_zzzz_0_xyyyyz_1,  \
                             g_0_zzzz_0_xyyyz_1,   \
                             g_0_zzzz_0_xyyyzz_0,  \
                             g_0_zzzz_0_xyyyzz_1,  \
                             g_0_zzzz_0_xyyzz_1,   \
                             g_0_zzzz_0_xyyzzz_0,  \
                             g_0_zzzz_0_xyyzzz_1,  \
                             g_0_zzzz_0_xyzzz_1,   \
                             g_0_zzzz_0_xyzzzz_0,  \
                             g_0_zzzz_0_xyzzzz_1,  \
                             g_0_zzzz_0_xzzzz_1,   \
                             g_0_zzzz_0_xzzzzz_0,  \
                             g_0_zzzz_0_xzzzzz_1,  \
                             g_0_zzzz_0_yyyyy_1,   \
                             g_0_zzzz_0_yyyyyy_0,  \
                             g_0_zzzz_0_yyyyyy_1,  \
                             g_0_zzzz_0_yyyyyz_0,  \
                             g_0_zzzz_0_yyyyyz_1,  \
                             g_0_zzzz_0_yyyyz_1,   \
                             g_0_zzzz_0_yyyyzz_0,  \
                             g_0_zzzz_0_yyyyzz_1,  \
                             g_0_zzzz_0_yyyzz_1,   \
                             g_0_zzzz_0_yyyzzz_0,  \
                             g_0_zzzz_0_yyyzzz_1,  \
                             g_0_zzzz_0_yyzzz_1,   \
                             g_0_zzzz_0_yyzzzz_0,  \
                             g_0_zzzz_0_yyzzzz_1,  \
                             g_0_zzzz_0_yzzzz_1,   \
                             g_0_zzzz_0_yzzzzz_0,  \
                             g_0_zzzz_0_yzzzzz_1,  \
                             g_0_zzzz_0_zzzzz_1,   \
                             g_0_zzzz_0_zzzzzz_0,  \
                             g_0_zzzz_0_zzzzzz_1,  \
                             wp_x,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzz_0_xxxxxx_0[i] = 6.0 * g_0_zzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxx_0[i] * pb_x + g_0_zzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxy_0[i] = 5.0 * g_0_zzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxy_0[i] * pb_x + g_0_zzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxz_0[i] = 5.0 * g_0_zzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxz_0[i] * pb_x + g_0_zzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyy_0[i] = 4.0 * g_0_zzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyy_0[i] * pb_x + g_0_zzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyz_0[i] = 4.0 * g_0_zzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyz_0[i] * pb_x + g_0_zzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxzz_0[i] = 4.0 * g_0_zzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxzz_0[i] * pb_x + g_0_zzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyy_0[i] = 3.0 * g_0_zzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyy_0[i] * pb_x + g_0_zzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyz_0[i] = 3.0 * g_0_zzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyz_0[i] * pb_x + g_0_zzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyzz_0[i] = 3.0 * g_0_zzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzz_0[i] * pb_x + g_0_zzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxzzz_0[i] = 3.0 * g_0_zzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxzzz_0[i] * pb_x + g_0_zzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyy_0[i] = 2.0 * g_0_zzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyy_0[i] * pb_x + g_0_zzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyz_0[i] = 2.0 * g_0_zzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyz_0[i] * pb_x + g_0_zzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyzz_0[i] = 2.0 * g_0_zzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzz_0[i] * pb_x + g_0_zzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyzzz_0[i] = 2.0 * g_0_zzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzz_0[i] * pb_x + g_0_zzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxzzzz_0[i] = 2.0 * g_0_zzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzzzz_0[i] * pb_x + g_0_zzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyy_0[i] = g_0_zzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyy_0[i] * pb_x + g_0_zzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyz_0[i] = g_0_zzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyz_0[i] * pb_x + g_0_zzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyzz_0[i] = g_0_zzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzz_0[i] * pb_x + g_0_zzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyzzz_0[i] = g_0_zzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzz_0[i] * pb_x + g_0_zzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyzzzz_0[i] = g_0_zzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzz_0[i] * pb_x + g_0_zzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xzzzzz_0[i] = g_0_zzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzzzz_0[i] * pb_x + g_0_zzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyy_0[i] = g_0_zzzz_0_yyyyyy_0[i] * pb_x + g_0_zzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyz_0[i] = g_0_zzzz_0_yyyyyz_0[i] * pb_x + g_0_zzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyzz_0[i] = g_0_zzzz_0_yyyyzz_0[i] * pb_x + g_0_zzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyzzz_0[i] = g_0_zzzz_0_yyyzzz_0[i] * pb_x + g_0_zzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyzzzz_0[i] = g_0_zzzz_0_yyzzzz_0[i] * pb_x + g_0_zzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yzzzzz_0[i] = g_0_zzzz_0_yzzzzz_0[i] * pb_x + g_0_zzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_zzzzzz_0[i] = g_0_zzzz_0_zzzzzz_0[i] * pb_x + g_0_zzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 420-448 components of targeted buffer : SHSI

    auto g_0_yyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 420);

    auto g_0_yyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 421);

    auto g_0_yyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 422);

    auto g_0_yyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 423);

    auto g_0_yyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 424);

    auto g_0_yyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 425);

    auto g_0_yyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 426);

    auto g_0_yyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 427);

    auto g_0_yyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 428);

    auto g_0_yyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 429);

    auto g_0_yyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 430);

    auto g_0_yyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 431);

    auto g_0_yyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 432);

    auto g_0_yyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 433);

    auto g_0_yyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 434);

    auto g_0_yyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 435);

    auto g_0_yyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 436);

    auto g_0_yyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 437);

    auto g_0_yyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 438);

    auto g_0_yyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 439);

    auto g_0_yyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 440);

    auto g_0_yyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 441);

    auto g_0_yyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 442);

    auto g_0_yyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 443);

    auto g_0_yyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 444);

    auto g_0_yyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 445);

    auto g_0_yyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 446);

    auto g_0_yyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 447);

#pragma omp simd aligned(g_0_yyy_0_xxxxxx_0,       \
                             g_0_yyy_0_xxxxxx_1,   \
                             g_0_yyy_0_xxxxxy_0,   \
                             g_0_yyy_0_xxxxxy_1,   \
                             g_0_yyy_0_xxxxxz_0,   \
                             g_0_yyy_0_xxxxxz_1,   \
                             g_0_yyy_0_xxxxyy_0,   \
                             g_0_yyy_0_xxxxyy_1,   \
                             g_0_yyy_0_xxxxyz_0,   \
                             g_0_yyy_0_xxxxyz_1,   \
                             g_0_yyy_0_xxxxzz_0,   \
                             g_0_yyy_0_xxxxzz_1,   \
                             g_0_yyy_0_xxxyyy_0,   \
                             g_0_yyy_0_xxxyyy_1,   \
                             g_0_yyy_0_xxxyyz_0,   \
                             g_0_yyy_0_xxxyyz_1,   \
                             g_0_yyy_0_xxxyzz_0,   \
                             g_0_yyy_0_xxxyzz_1,   \
                             g_0_yyy_0_xxxzzz_0,   \
                             g_0_yyy_0_xxxzzz_1,   \
                             g_0_yyy_0_xxyyyy_0,   \
                             g_0_yyy_0_xxyyyy_1,   \
                             g_0_yyy_0_xxyyyz_0,   \
                             g_0_yyy_0_xxyyyz_1,   \
                             g_0_yyy_0_xxyyzz_0,   \
                             g_0_yyy_0_xxyyzz_1,   \
                             g_0_yyy_0_xxyzzz_0,   \
                             g_0_yyy_0_xxyzzz_1,   \
                             g_0_yyy_0_xxzzzz_0,   \
                             g_0_yyy_0_xxzzzz_1,   \
                             g_0_yyy_0_xyyyyy_0,   \
                             g_0_yyy_0_xyyyyy_1,   \
                             g_0_yyy_0_xyyyyz_0,   \
                             g_0_yyy_0_xyyyyz_1,   \
                             g_0_yyy_0_xyyyzz_0,   \
                             g_0_yyy_0_xyyyzz_1,   \
                             g_0_yyy_0_xyyzzz_0,   \
                             g_0_yyy_0_xyyzzz_1,   \
                             g_0_yyy_0_xyzzzz_0,   \
                             g_0_yyy_0_xyzzzz_1,   \
                             g_0_yyy_0_xzzzzz_0,   \
                             g_0_yyy_0_xzzzzz_1,   \
                             g_0_yyy_0_yyyyyy_0,   \
                             g_0_yyy_0_yyyyyy_1,   \
                             g_0_yyy_0_yyyyyz_0,   \
                             g_0_yyy_0_yyyyyz_1,   \
                             g_0_yyy_0_yyyyzz_0,   \
                             g_0_yyy_0_yyyyzz_1,   \
                             g_0_yyy_0_yyyzzz_0,   \
                             g_0_yyy_0_yyyzzz_1,   \
                             g_0_yyy_0_yyzzzz_0,   \
                             g_0_yyy_0_yyzzzz_1,   \
                             g_0_yyy_0_yzzzzz_0,   \
                             g_0_yyy_0_yzzzzz_1,   \
                             g_0_yyy_0_zzzzzz_0,   \
                             g_0_yyy_0_zzzzzz_1,   \
                             g_0_yyyy_0_xxxxx_1,   \
                             g_0_yyyy_0_xxxxxx_0,  \
                             g_0_yyyy_0_xxxxxx_1,  \
                             g_0_yyyy_0_xxxxxy_0,  \
                             g_0_yyyy_0_xxxxxy_1,  \
                             g_0_yyyy_0_xxxxxz_0,  \
                             g_0_yyyy_0_xxxxxz_1,  \
                             g_0_yyyy_0_xxxxy_1,   \
                             g_0_yyyy_0_xxxxyy_0,  \
                             g_0_yyyy_0_xxxxyy_1,  \
                             g_0_yyyy_0_xxxxyz_0,  \
                             g_0_yyyy_0_xxxxyz_1,  \
                             g_0_yyyy_0_xxxxz_1,   \
                             g_0_yyyy_0_xxxxzz_0,  \
                             g_0_yyyy_0_xxxxzz_1,  \
                             g_0_yyyy_0_xxxyy_1,   \
                             g_0_yyyy_0_xxxyyy_0,  \
                             g_0_yyyy_0_xxxyyy_1,  \
                             g_0_yyyy_0_xxxyyz_0,  \
                             g_0_yyyy_0_xxxyyz_1,  \
                             g_0_yyyy_0_xxxyz_1,   \
                             g_0_yyyy_0_xxxyzz_0,  \
                             g_0_yyyy_0_xxxyzz_1,  \
                             g_0_yyyy_0_xxxzz_1,   \
                             g_0_yyyy_0_xxxzzz_0,  \
                             g_0_yyyy_0_xxxzzz_1,  \
                             g_0_yyyy_0_xxyyy_1,   \
                             g_0_yyyy_0_xxyyyy_0,  \
                             g_0_yyyy_0_xxyyyy_1,  \
                             g_0_yyyy_0_xxyyyz_0,  \
                             g_0_yyyy_0_xxyyyz_1,  \
                             g_0_yyyy_0_xxyyz_1,   \
                             g_0_yyyy_0_xxyyzz_0,  \
                             g_0_yyyy_0_xxyyzz_1,  \
                             g_0_yyyy_0_xxyzz_1,   \
                             g_0_yyyy_0_xxyzzz_0,  \
                             g_0_yyyy_0_xxyzzz_1,  \
                             g_0_yyyy_0_xxzzz_1,   \
                             g_0_yyyy_0_xxzzzz_0,  \
                             g_0_yyyy_0_xxzzzz_1,  \
                             g_0_yyyy_0_xyyyy_1,   \
                             g_0_yyyy_0_xyyyyy_0,  \
                             g_0_yyyy_0_xyyyyy_1,  \
                             g_0_yyyy_0_xyyyyz_0,  \
                             g_0_yyyy_0_xyyyyz_1,  \
                             g_0_yyyy_0_xyyyz_1,   \
                             g_0_yyyy_0_xyyyzz_0,  \
                             g_0_yyyy_0_xyyyzz_1,  \
                             g_0_yyyy_0_xyyzz_1,   \
                             g_0_yyyy_0_xyyzzz_0,  \
                             g_0_yyyy_0_xyyzzz_1,  \
                             g_0_yyyy_0_xyzzz_1,   \
                             g_0_yyyy_0_xyzzzz_0,  \
                             g_0_yyyy_0_xyzzzz_1,  \
                             g_0_yyyy_0_xzzzz_1,   \
                             g_0_yyyy_0_xzzzzz_0,  \
                             g_0_yyyy_0_xzzzzz_1,  \
                             g_0_yyyy_0_yyyyy_1,   \
                             g_0_yyyy_0_yyyyyy_0,  \
                             g_0_yyyy_0_yyyyyy_1,  \
                             g_0_yyyy_0_yyyyyz_0,  \
                             g_0_yyyy_0_yyyyyz_1,  \
                             g_0_yyyy_0_yyyyz_1,   \
                             g_0_yyyy_0_yyyyzz_0,  \
                             g_0_yyyy_0_yyyyzz_1,  \
                             g_0_yyyy_0_yyyzz_1,   \
                             g_0_yyyy_0_yyyzzz_0,  \
                             g_0_yyyy_0_yyyzzz_1,  \
                             g_0_yyyy_0_yyzzz_1,   \
                             g_0_yyyy_0_yyzzzz_0,  \
                             g_0_yyyy_0_yyzzzz_1,  \
                             g_0_yyyy_0_yzzzz_1,   \
                             g_0_yyyy_0_yzzzzz_0,  \
                             g_0_yyyy_0_yzzzzz_1,  \
                             g_0_yyyy_0_zzzzz_1,   \
                             g_0_yyyy_0_zzzzzz_0,  \
                             g_0_yyyy_0_zzzzzz_1,  \
                             g_0_yyyyy_0_xxxxxx_0, \
                             g_0_yyyyy_0_xxxxxy_0, \
                             g_0_yyyyy_0_xxxxxz_0, \
                             g_0_yyyyy_0_xxxxyy_0, \
                             g_0_yyyyy_0_xxxxyz_0, \
                             g_0_yyyyy_0_xxxxzz_0, \
                             g_0_yyyyy_0_xxxyyy_0, \
                             g_0_yyyyy_0_xxxyyz_0, \
                             g_0_yyyyy_0_xxxyzz_0, \
                             g_0_yyyyy_0_xxxzzz_0, \
                             g_0_yyyyy_0_xxyyyy_0, \
                             g_0_yyyyy_0_xxyyyz_0, \
                             g_0_yyyyy_0_xxyyzz_0, \
                             g_0_yyyyy_0_xxyzzz_0, \
                             g_0_yyyyy_0_xxzzzz_0, \
                             g_0_yyyyy_0_xyyyyy_0, \
                             g_0_yyyyy_0_xyyyyz_0, \
                             g_0_yyyyy_0_xyyyzz_0, \
                             g_0_yyyyy_0_xyyzzz_0, \
                             g_0_yyyyy_0_xyzzzz_0, \
                             g_0_yyyyy_0_xzzzzz_0, \
                             g_0_yyyyy_0_yyyyyy_0, \
                             g_0_yyyyy_0_yyyyyz_0, \
                             g_0_yyyyy_0_yyyyzz_0, \
                             g_0_yyyyy_0_yyyzzz_0, \
                             g_0_yyyyy_0_yyzzzz_0, \
                             g_0_yyyyy_0_yzzzzz_0, \
                             g_0_yyyyy_0_zzzzzz_0, \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyy_0_xxxxxx_0[i] = 4.0 * g_0_yyy_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxx_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxx_0[i] * pb_y +
                                  g_0_yyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxy_0[i] = 4.0 * g_0_yyy_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxx_1[i] * fi_abcd_0 +
                                  g_0_yyyy_0_xxxxxy_0[i] * pb_y + g_0_yyyy_0_xxxxxy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxz_0[i] = 4.0 * g_0_yyy_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxz_0[i] * pb_y +
                                  g_0_yyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyy_0[i] = 4.0 * g_0_yyy_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyy_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyy_0[i] * pb_y + g_0_yyyy_0_xxxxyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyz_0[i] = 4.0 * g_0_yyy_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxz_1[i] * fi_abcd_0 +
                                  g_0_yyyy_0_xxxxyz_0[i] * pb_y + g_0_yyyy_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxzz_0[i] = 4.0 * g_0_yyy_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxzz_0[i] * pb_y +
                                  g_0_yyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyy_0[i] = 4.0 * g_0_yyy_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyy_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyy_0[i] * pb_y + g_0_yyyy_0_xxxyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyz_0[i] = 4.0 * g_0_yyy_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyz_0[i] * pb_y + g_0_yyyy_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyzz_0[i] = 4.0 * g_0_yyy_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxzz_1[i] * fi_abcd_0 +
                                  g_0_yyyy_0_xxxyzz_0[i] * pb_y + g_0_yyyy_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxzzz_0[i] = 4.0 * g_0_yyy_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxzzz_0[i] * pb_y +
                                  g_0_yyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyy_0[i] = 4.0 * g_0_yyy_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyy_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyy_0[i] * pb_y + g_0_yyyy_0_xxyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyz_0[i] = 4.0 * g_0_yyy_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyz_0[i] * pb_y + g_0_yyyy_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyzz_0[i] = 4.0 * g_0_yyy_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzz_0[i] * pb_y + g_0_yyyy_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyzzz_0[i] = 4.0 * g_0_yyy_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzzz_1[i] * fi_abcd_0 +
                                  g_0_yyyy_0_xxyzzz_0[i] * pb_y + g_0_yyyy_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxzzzz_0[i] = 4.0 * g_0_yyy_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzzzz_0[i] * pb_y +
                                  g_0_yyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyy_0[i] = 4.0 * g_0_yyy_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyy_1[i] * fti_ab_0 +
                                  5.0 * g_0_yyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyy_0[i] * pb_y + g_0_yyyy_0_xyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyz_0[i] = 4.0 * g_0_yyy_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyz_0[i] * pb_y + g_0_yyyy_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyzz_0[i] = 4.0 * g_0_yyy_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzz_0[i] * pb_y + g_0_yyyy_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyzzz_0[i] = 4.0 * g_0_yyy_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzz_0[i] * pb_y + g_0_yyyy_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyzzzz_0[i] = 4.0 * g_0_yyy_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzzz_1[i] * fi_abcd_0 +
                                  g_0_yyyy_0_xyzzzz_0[i] * pb_y + g_0_yyyy_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xzzzzz_0[i] = 4.0 * g_0_yyy_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzzzz_0[i] * pb_y +
                                  g_0_yyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyy_0[i] = 4.0 * g_0_yyy_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyy_1[i] * fti_ab_0 +
                                  6.0 * g_0_yyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyy_0[i] * pb_y + g_0_yyyy_0_yyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyz_0[i] = 4.0 * g_0_yyy_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyz_1[i] * fti_ab_0 +
                                  5.0 * g_0_yyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyz_0[i] * pb_y + g_0_yyyy_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyzz_0[i] = 4.0 * g_0_yyy_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyzz_0[i] * pb_y + g_0_yyyy_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyzzz_0[i] = 4.0 * g_0_yyy_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyzzz_0[i] * pb_y + g_0_yyyy_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyzzzz_0[i] = 4.0 * g_0_yyy_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzzzz_0[i] * pb_y + g_0_yyyy_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yzzzzz_0[i] = 4.0 * g_0_yyy_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzzz_1[i] * fi_abcd_0 +
                                  g_0_yyyy_0_yzzzzz_0[i] * pb_y + g_0_yyyy_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_zzzzzz_0[i] = 4.0 * g_0_yyy_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzzzz_0[i] * pb_y +
                                  g_0_yyyy_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 448-476 components of targeted buffer : SHSI

    auto g_0_yyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 448);

    auto g_0_yyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 449);

    auto g_0_yyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 450);

    auto g_0_yyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 451);

    auto g_0_yyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 452);

    auto g_0_yyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 453);

    auto g_0_yyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 454);

    auto g_0_yyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 455);

    auto g_0_yyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 456);

    auto g_0_yyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 457);

    auto g_0_yyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 458);

    auto g_0_yyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 459);

    auto g_0_yyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 460);

    auto g_0_yyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 461);

    auto g_0_yyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 462);

    auto g_0_yyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 463);

    auto g_0_yyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 464);

    auto g_0_yyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 465);

    auto g_0_yyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 466);

    auto g_0_yyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 467);

    auto g_0_yyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 468);

    auto g_0_yyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 469);

    auto g_0_yyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 470);

    auto g_0_yyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 471);

    auto g_0_yyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 472);

    auto g_0_yyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 473);

    auto g_0_yyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 474);

    auto g_0_yyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 475);

#pragma omp simd aligned(g_0_yyyy_0_xxxxx_1,       \
                             g_0_yyyy_0_xxxxxx_0,  \
                             g_0_yyyy_0_xxxxxx_1,  \
                             g_0_yyyy_0_xxxxxy_0,  \
                             g_0_yyyy_0_xxxxxy_1,  \
                             g_0_yyyy_0_xxxxxz_0,  \
                             g_0_yyyy_0_xxxxxz_1,  \
                             g_0_yyyy_0_xxxxy_1,   \
                             g_0_yyyy_0_xxxxyy_0,  \
                             g_0_yyyy_0_xxxxyy_1,  \
                             g_0_yyyy_0_xxxxyz_0,  \
                             g_0_yyyy_0_xxxxyz_1,  \
                             g_0_yyyy_0_xxxxz_1,   \
                             g_0_yyyy_0_xxxxzz_0,  \
                             g_0_yyyy_0_xxxxzz_1,  \
                             g_0_yyyy_0_xxxyy_1,   \
                             g_0_yyyy_0_xxxyyy_0,  \
                             g_0_yyyy_0_xxxyyy_1,  \
                             g_0_yyyy_0_xxxyyz_0,  \
                             g_0_yyyy_0_xxxyyz_1,  \
                             g_0_yyyy_0_xxxyz_1,   \
                             g_0_yyyy_0_xxxyzz_0,  \
                             g_0_yyyy_0_xxxyzz_1,  \
                             g_0_yyyy_0_xxxzz_1,   \
                             g_0_yyyy_0_xxxzzz_0,  \
                             g_0_yyyy_0_xxxzzz_1,  \
                             g_0_yyyy_0_xxyyy_1,   \
                             g_0_yyyy_0_xxyyyy_0,  \
                             g_0_yyyy_0_xxyyyy_1,  \
                             g_0_yyyy_0_xxyyyz_0,  \
                             g_0_yyyy_0_xxyyyz_1,  \
                             g_0_yyyy_0_xxyyz_1,   \
                             g_0_yyyy_0_xxyyzz_0,  \
                             g_0_yyyy_0_xxyyzz_1,  \
                             g_0_yyyy_0_xxyzz_1,   \
                             g_0_yyyy_0_xxyzzz_0,  \
                             g_0_yyyy_0_xxyzzz_1,  \
                             g_0_yyyy_0_xxzzz_1,   \
                             g_0_yyyy_0_xxzzzz_0,  \
                             g_0_yyyy_0_xxzzzz_1,  \
                             g_0_yyyy_0_xyyyy_1,   \
                             g_0_yyyy_0_xyyyyy_0,  \
                             g_0_yyyy_0_xyyyyy_1,  \
                             g_0_yyyy_0_xyyyyz_0,  \
                             g_0_yyyy_0_xyyyyz_1,  \
                             g_0_yyyy_0_xyyyz_1,   \
                             g_0_yyyy_0_xyyyzz_0,  \
                             g_0_yyyy_0_xyyyzz_1,  \
                             g_0_yyyy_0_xyyzz_1,   \
                             g_0_yyyy_0_xyyzzz_0,  \
                             g_0_yyyy_0_xyyzzz_1,  \
                             g_0_yyyy_0_xyzzz_1,   \
                             g_0_yyyy_0_xyzzzz_0,  \
                             g_0_yyyy_0_xyzzzz_1,  \
                             g_0_yyyy_0_xzzzz_1,   \
                             g_0_yyyy_0_xzzzzz_0,  \
                             g_0_yyyy_0_xzzzzz_1,  \
                             g_0_yyyy_0_yyyyy_1,   \
                             g_0_yyyy_0_yyyyyy_0,  \
                             g_0_yyyy_0_yyyyyy_1,  \
                             g_0_yyyy_0_yyyyyz_0,  \
                             g_0_yyyy_0_yyyyyz_1,  \
                             g_0_yyyy_0_yyyyz_1,   \
                             g_0_yyyy_0_yyyyzz_0,  \
                             g_0_yyyy_0_yyyyzz_1,  \
                             g_0_yyyy_0_yyyzz_1,   \
                             g_0_yyyy_0_yyyzzz_0,  \
                             g_0_yyyy_0_yyyzzz_1,  \
                             g_0_yyyy_0_yyzzz_1,   \
                             g_0_yyyy_0_yyzzzz_0,  \
                             g_0_yyyy_0_yyzzzz_1,  \
                             g_0_yyyy_0_yzzzz_1,   \
                             g_0_yyyy_0_yzzzzz_0,  \
                             g_0_yyyy_0_yzzzzz_1,  \
                             g_0_yyyy_0_zzzzz_1,   \
                             g_0_yyyy_0_zzzzzz_0,  \
                             g_0_yyyy_0_zzzzzz_1,  \
                             g_0_yyyyz_0_xxxxxx_0, \
                             g_0_yyyyz_0_xxxxxy_0, \
                             g_0_yyyyz_0_xxxxxz_0, \
                             g_0_yyyyz_0_xxxxyy_0, \
                             g_0_yyyyz_0_xxxxyz_0, \
                             g_0_yyyyz_0_xxxxzz_0, \
                             g_0_yyyyz_0_xxxyyy_0, \
                             g_0_yyyyz_0_xxxyyz_0, \
                             g_0_yyyyz_0_xxxyzz_0, \
                             g_0_yyyyz_0_xxxzzz_0, \
                             g_0_yyyyz_0_xxyyyy_0, \
                             g_0_yyyyz_0_xxyyyz_0, \
                             g_0_yyyyz_0_xxyyzz_0, \
                             g_0_yyyyz_0_xxyzzz_0, \
                             g_0_yyyyz_0_xxzzzz_0, \
                             g_0_yyyyz_0_xyyyyy_0, \
                             g_0_yyyyz_0_xyyyyz_0, \
                             g_0_yyyyz_0_xyyyzz_0, \
                             g_0_yyyyz_0_xyyzzz_0, \
                             g_0_yyyyz_0_xyzzzz_0, \
                             g_0_yyyyz_0_xzzzzz_0, \
                             g_0_yyyyz_0_yyyyyy_0, \
                             g_0_yyyyz_0_yyyyyz_0, \
                             g_0_yyyyz_0_yyyyzz_0, \
                             g_0_yyyyz_0_yyyzzz_0, \
                             g_0_yyyyz_0_yyzzzz_0, \
                             g_0_yyyyz_0_yzzzzz_0, \
                             g_0_yyyyz_0_zzzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyz_0_xxxxxx_0[i] = g_0_yyyy_0_xxxxxx_0[i] * pb_z + g_0_yyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxy_0[i] = g_0_yyyy_0_xxxxxy_0[i] * pb_z + g_0_yyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxz_0[i] = g_0_yyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxz_0[i] * pb_z + g_0_yyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyy_0[i] = g_0_yyyy_0_xxxxyy_0[i] * pb_z + g_0_yyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyz_0[i] = g_0_yyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyz_0[i] * pb_z + g_0_yyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxzz_0[i] = 2.0 * g_0_yyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxzz_0[i] * pb_z + g_0_yyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyy_0[i] = g_0_yyyy_0_xxxyyy_0[i] * pb_z + g_0_yyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyz_0[i] = g_0_yyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyz_0[i] * pb_z + g_0_yyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyzz_0[i] = 2.0 * g_0_yyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyzz_0[i] * pb_z + g_0_yyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxzzz_0[i] = 3.0 * g_0_yyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxzzz_0[i] * pb_z + g_0_yyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyy_0[i] = g_0_yyyy_0_xxyyyy_0[i] * pb_z + g_0_yyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyz_0[i] = g_0_yyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyz_0[i] * pb_z + g_0_yyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyzz_0[i] = 2.0 * g_0_yyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzz_0[i] * pb_z + g_0_yyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyzzz_0[i] = 3.0 * g_0_yyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzzz_0[i] * pb_z + g_0_yyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxzzzz_0[i] = 4.0 * g_0_yyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzzzz_0[i] * pb_z + g_0_yyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyy_0[i] = g_0_yyyy_0_xyyyyy_0[i] * pb_z + g_0_yyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyz_0[i] = g_0_yyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyz_0[i] * pb_z + g_0_yyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyzz_0[i] = 2.0 * g_0_yyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzz_0[i] * pb_z + g_0_yyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyzzz_0[i] = 3.0 * g_0_yyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzz_0[i] * pb_z + g_0_yyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyzzzz_0[i] = 4.0 * g_0_yyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzzz_0[i] * pb_z + g_0_yyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xzzzzz_0[i] = 5.0 * g_0_yyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzzzz_0[i] * pb_z + g_0_yyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyy_0[i] = g_0_yyyy_0_yyyyyy_0[i] * pb_z + g_0_yyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyz_0[i] = g_0_yyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyz_0[i] * pb_z + g_0_yyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyzz_0[i] = 2.0 * g_0_yyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyzz_0[i] * pb_z + g_0_yyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyzzz_0[i] = 3.0 * g_0_yyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyzzz_0[i] * pb_z + g_0_yyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyzzzz_0[i] = 4.0 * g_0_yyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzzzz_0[i] * pb_z + g_0_yyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yzzzzz_0[i] = 5.0 * g_0_yyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzzzzz_0[i] * pb_z + g_0_yyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_zzzzzz_0[i] = 6.0 * g_0_yyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_zzzzzz_0[i] * pb_z + g_0_yyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 476-504 components of targeted buffer : SHSI

    auto g_0_yyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 476);

    auto g_0_yyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 477);

    auto g_0_yyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 478);

    auto g_0_yyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 479);

    auto g_0_yyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 480);

    auto g_0_yyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 481);

    auto g_0_yyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 482);

    auto g_0_yyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 483);

    auto g_0_yyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 484);

    auto g_0_yyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 485);

    auto g_0_yyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 486);

    auto g_0_yyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 487);

    auto g_0_yyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 488);

    auto g_0_yyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 489);

    auto g_0_yyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 490);

    auto g_0_yyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 491);

    auto g_0_yyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 492);

    auto g_0_yyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 493);

    auto g_0_yyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 494);

    auto g_0_yyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 495);

    auto g_0_yyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 496);

    auto g_0_yyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 497);

    auto g_0_yyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 498);

    auto g_0_yyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 499);

    auto g_0_yyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 500);

    auto g_0_yyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 501);

    auto g_0_yyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 502);

    auto g_0_yyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 503);

#pragma omp simd aligned(g_0_yyy_0_xxxxxy_0,       \
                             g_0_yyy_0_xxxxxy_1,   \
                             g_0_yyy_0_xxxxyy_0,   \
                             g_0_yyy_0_xxxxyy_1,   \
                             g_0_yyy_0_xxxyyy_0,   \
                             g_0_yyy_0_xxxyyy_1,   \
                             g_0_yyy_0_xxyyyy_0,   \
                             g_0_yyy_0_xxyyyy_1,   \
                             g_0_yyy_0_xyyyyy_0,   \
                             g_0_yyy_0_xyyyyy_1,   \
                             g_0_yyy_0_yyyyyy_0,   \
                             g_0_yyy_0_yyyyyy_1,   \
                             g_0_yyyz_0_xxxxxy_0,  \
                             g_0_yyyz_0_xxxxxy_1,  \
                             g_0_yyyz_0_xxxxyy_0,  \
                             g_0_yyyz_0_xxxxyy_1,  \
                             g_0_yyyz_0_xxxyyy_0,  \
                             g_0_yyyz_0_xxxyyy_1,  \
                             g_0_yyyz_0_xxyyyy_0,  \
                             g_0_yyyz_0_xxyyyy_1,  \
                             g_0_yyyz_0_xyyyyy_0,  \
                             g_0_yyyz_0_xyyyyy_1,  \
                             g_0_yyyz_0_yyyyyy_0,  \
                             g_0_yyyz_0_yyyyyy_1,  \
                             g_0_yyyzz_0_xxxxxx_0, \
                             g_0_yyyzz_0_xxxxxy_0, \
                             g_0_yyyzz_0_xxxxxz_0, \
                             g_0_yyyzz_0_xxxxyy_0, \
                             g_0_yyyzz_0_xxxxyz_0, \
                             g_0_yyyzz_0_xxxxzz_0, \
                             g_0_yyyzz_0_xxxyyy_0, \
                             g_0_yyyzz_0_xxxyyz_0, \
                             g_0_yyyzz_0_xxxyzz_0, \
                             g_0_yyyzz_0_xxxzzz_0, \
                             g_0_yyyzz_0_xxyyyy_0, \
                             g_0_yyyzz_0_xxyyyz_0, \
                             g_0_yyyzz_0_xxyyzz_0, \
                             g_0_yyyzz_0_xxyzzz_0, \
                             g_0_yyyzz_0_xxzzzz_0, \
                             g_0_yyyzz_0_xyyyyy_0, \
                             g_0_yyyzz_0_xyyyyz_0, \
                             g_0_yyyzz_0_xyyyzz_0, \
                             g_0_yyyzz_0_xyyzzz_0, \
                             g_0_yyyzz_0_xyzzzz_0, \
                             g_0_yyyzz_0_xzzzzz_0, \
                             g_0_yyyzz_0_yyyyyy_0, \
                             g_0_yyyzz_0_yyyyyz_0, \
                             g_0_yyyzz_0_yyyyzz_0, \
                             g_0_yyyzz_0_yyyzzz_0, \
                             g_0_yyyzz_0_yyzzzz_0, \
                             g_0_yyyzz_0_yzzzzz_0, \
                             g_0_yyyzz_0_zzzzzz_0, \
                             g_0_yyzz_0_xxxxxx_0,  \
                             g_0_yyzz_0_xxxxxx_1,  \
                             g_0_yyzz_0_xxxxxz_0,  \
                             g_0_yyzz_0_xxxxxz_1,  \
                             g_0_yyzz_0_xxxxyz_0,  \
                             g_0_yyzz_0_xxxxyz_1,  \
                             g_0_yyzz_0_xxxxz_1,   \
                             g_0_yyzz_0_xxxxzz_0,  \
                             g_0_yyzz_0_xxxxzz_1,  \
                             g_0_yyzz_0_xxxyyz_0,  \
                             g_0_yyzz_0_xxxyyz_1,  \
                             g_0_yyzz_0_xxxyz_1,   \
                             g_0_yyzz_0_xxxyzz_0,  \
                             g_0_yyzz_0_xxxyzz_1,  \
                             g_0_yyzz_0_xxxzz_1,   \
                             g_0_yyzz_0_xxxzzz_0,  \
                             g_0_yyzz_0_xxxzzz_1,  \
                             g_0_yyzz_0_xxyyyz_0,  \
                             g_0_yyzz_0_xxyyyz_1,  \
                             g_0_yyzz_0_xxyyz_1,   \
                             g_0_yyzz_0_xxyyzz_0,  \
                             g_0_yyzz_0_xxyyzz_1,  \
                             g_0_yyzz_0_xxyzz_1,   \
                             g_0_yyzz_0_xxyzzz_0,  \
                             g_0_yyzz_0_xxyzzz_1,  \
                             g_0_yyzz_0_xxzzz_1,   \
                             g_0_yyzz_0_xxzzzz_0,  \
                             g_0_yyzz_0_xxzzzz_1,  \
                             g_0_yyzz_0_xyyyyz_0,  \
                             g_0_yyzz_0_xyyyyz_1,  \
                             g_0_yyzz_0_xyyyz_1,   \
                             g_0_yyzz_0_xyyyzz_0,  \
                             g_0_yyzz_0_xyyyzz_1,  \
                             g_0_yyzz_0_xyyzz_1,   \
                             g_0_yyzz_0_xyyzzz_0,  \
                             g_0_yyzz_0_xyyzzz_1,  \
                             g_0_yyzz_0_xyzzz_1,   \
                             g_0_yyzz_0_xyzzzz_0,  \
                             g_0_yyzz_0_xyzzzz_1,  \
                             g_0_yyzz_0_xzzzz_1,   \
                             g_0_yyzz_0_xzzzzz_0,  \
                             g_0_yyzz_0_xzzzzz_1,  \
                             g_0_yyzz_0_yyyyyz_0,  \
                             g_0_yyzz_0_yyyyyz_1,  \
                             g_0_yyzz_0_yyyyz_1,   \
                             g_0_yyzz_0_yyyyzz_0,  \
                             g_0_yyzz_0_yyyyzz_1,  \
                             g_0_yyzz_0_yyyzz_1,   \
                             g_0_yyzz_0_yyyzzz_0,  \
                             g_0_yyzz_0_yyyzzz_1,  \
                             g_0_yyzz_0_yyzzz_1,   \
                             g_0_yyzz_0_yyzzzz_0,  \
                             g_0_yyzz_0_yyzzzz_1,  \
                             g_0_yyzz_0_yzzzz_1,   \
                             g_0_yyzz_0_yzzzzz_0,  \
                             g_0_yyzz_0_yzzzzz_1,  \
                             g_0_yyzz_0_zzzzz_1,   \
                             g_0_yyzz_0_zzzzzz_0,  \
                             g_0_yyzz_0_zzzzzz_1,  \
                             g_0_yzz_0_xxxxxx_0,   \
                             g_0_yzz_0_xxxxxx_1,   \
                             g_0_yzz_0_xxxxxz_0,   \
                             g_0_yzz_0_xxxxxz_1,   \
                             g_0_yzz_0_xxxxyz_0,   \
                             g_0_yzz_0_xxxxyz_1,   \
                             g_0_yzz_0_xxxxzz_0,   \
                             g_0_yzz_0_xxxxzz_1,   \
                             g_0_yzz_0_xxxyyz_0,   \
                             g_0_yzz_0_xxxyyz_1,   \
                             g_0_yzz_0_xxxyzz_0,   \
                             g_0_yzz_0_xxxyzz_1,   \
                             g_0_yzz_0_xxxzzz_0,   \
                             g_0_yzz_0_xxxzzz_1,   \
                             g_0_yzz_0_xxyyyz_0,   \
                             g_0_yzz_0_xxyyyz_1,   \
                             g_0_yzz_0_xxyyzz_0,   \
                             g_0_yzz_0_xxyyzz_1,   \
                             g_0_yzz_0_xxyzzz_0,   \
                             g_0_yzz_0_xxyzzz_1,   \
                             g_0_yzz_0_xxzzzz_0,   \
                             g_0_yzz_0_xxzzzz_1,   \
                             g_0_yzz_0_xyyyyz_0,   \
                             g_0_yzz_0_xyyyyz_1,   \
                             g_0_yzz_0_xyyyzz_0,   \
                             g_0_yzz_0_xyyyzz_1,   \
                             g_0_yzz_0_xyyzzz_0,   \
                             g_0_yzz_0_xyyzzz_1,   \
                             g_0_yzz_0_xyzzzz_0,   \
                             g_0_yzz_0_xyzzzz_1,   \
                             g_0_yzz_0_xzzzzz_0,   \
                             g_0_yzz_0_xzzzzz_1,   \
                             g_0_yzz_0_yyyyyz_0,   \
                             g_0_yzz_0_yyyyyz_1,   \
                             g_0_yzz_0_yyyyzz_0,   \
                             g_0_yzz_0_yyyyzz_1,   \
                             g_0_yzz_0_yyyzzz_0,   \
                             g_0_yzz_0_yyyzzz_1,   \
                             g_0_yzz_0_yyzzzz_0,   \
                             g_0_yzz_0_yyzzzz_1,   \
                             g_0_yzz_0_yzzzzz_0,   \
                             g_0_yzz_0_yzzzzz_1,   \
                             g_0_yzz_0_zzzzzz_0,   \
                             g_0_yzz_0_zzzzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzz_0_xxxxxx_0[i] = 2.0 * g_0_yzz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxx_0[i] * pb_y +
                                  g_0_yyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxy_0[i] =
            g_0_yyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxxy_0[i] * pb_z + g_0_yyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxxz_0[i] = 2.0 * g_0_yzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxz_0[i] * pb_y +
                                  g_0_yyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxyy_0[i] =
            g_0_yyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxyy_0[i] * pb_z + g_0_yyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxyz_0[i] = 2.0 * g_0_yzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxyz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxz_1[i] * fi_abcd_0 +
                                  g_0_yyzz_0_xxxxyz_0[i] * pb_y + g_0_yyzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxzz_0[i] = 2.0 * g_0_yzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxzz_0[i] * pb_y +
                                  g_0_yyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyyy_0[i] =
            g_0_yyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxyyy_0[i] * pb_z + g_0_yyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxyyz_0[i] = 2.0 * g_0_yzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyyz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyz_0[i] * pb_y + g_0_yyzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyzz_0[i] = 2.0 * g_0_yzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxzz_1[i] * fi_abcd_0 +
                                  g_0_yyzz_0_xxxyzz_0[i] * pb_y + g_0_yyzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxzzz_0[i] = 2.0 * g_0_yzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxzzz_0[i] * pb_y +
                                  g_0_yyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyyy_0[i] =
            g_0_yyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxyyyy_0[i] * pb_z + g_0_yyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxyyyz_0[i] = 2.0 * g_0_yzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyyz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyz_0[i] * pb_y + g_0_yyzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyzz_0[i] = 2.0 * g_0_yzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyzz_0[i] * pb_y + g_0_yyzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyzzz_0[i] = 2.0 * g_0_yzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzzz_1[i] * fi_abcd_0 +
                                  g_0_yyzz_0_xxyzzz_0[i] * pb_y + g_0_yyzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxzzzz_0[i] = 2.0 * g_0_yzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzzzz_0[i] * pb_y +
                                  g_0_yyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyyy_0[i] =
            g_0_yyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xyyyyy_0[i] * pb_z + g_0_yyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xyyyyz_0[i] = 2.0 * g_0_yzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyyz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyz_0[i] * pb_y + g_0_yyzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyzz_0[i] = 2.0 * g_0_yzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyzz_0[i] * pb_y + g_0_yyzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyzzz_0[i] = 2.0 * g_0_yzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyzzz_0[i] * pb_y + g_0_yyzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyzzzz_0[i] = 2.0 * g_0_yzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzzz_1[i] * fi_abcd_0 +
                                  g_0_yyzz_0_xyzzzz_0[i] * pb_y + g_0_yyzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xzzzzz_0[i] = 2.0 * g_0_yzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzzzz_0[i] * pb_y +
                                  g_0_yyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyyy_0[i] =
            g_0_yyy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_yyyyyy_0[i] * pb_z + g_0_yyyz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_yyyyyz_0[i] = 2.0 * g_0_yzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyyz_1[i] * fti_ab_0 +
                                  5.0 * g_0_yyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyyz_0[i] * pb_y + g_0_yyzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyzz_0[i] = 2.0 * g_0_yzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_yyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyzz_0[i] * pb_y + g_0_yyzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyzzz_0[i] = 2.0 * g_0_yzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_yyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyzzz_0[i] * pb_y + g_0_yyzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyzzzz_0[i] = 2.0 * g_0_yzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyzzzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_yyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyzzzz_0[i] * pb_y + g_0_yyzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yzzzzz_0[i] = 2.0 * g_0_yzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzzz_1[i] * fi_abcd_0 +
                                  g_0_yyzz_0_yzzzzz_0[i] * pb_y + g_0_yyzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_zzzzzz_0[i] = 2.0 * g_0_yzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzzzz_0[i] * pb_y +
                                  g_0_yyzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 504-532 components of targeted buffer : SHSI

    auto g_0_yyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 504);

    auto g_0_yyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 505);

    auto g_0_yyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 506);

    auto g_0_yyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 507);

    auto g_0_yyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 508);

    auto g_0_yyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 509);

    auto g_0_yyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 510);

    auto g_0_yyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 511);

    auto g_0_yyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 512);

    auto g_0_yyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 513);

    auto g_0_yyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 514);

    auto g_0_yyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 515);

    auto g_0_yyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 516);

    auto g_0_yyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 517);

    auto g_0_yyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 518);

    auto g_0_yyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 519);

    auto g_0_yyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 520);

    auto g_0_yyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 521);

    auto g_0_yyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 522);

    auto g_0_yyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 523);

    auto g_0_yyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 524);

    auto g_0_yyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 525);

    auto g_0_yyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 526);

    auto g_0_yyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 527);

    auto g_0_yyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 528);

    auto g_0_yyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 529);

    auto g_0_yyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 530);

    auto g_0_yyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 531);

#pragma omp simd aligned(g_0_yyz_0_xxxxxy_0,       \
                             g_0_yyz_0_xxxxxy_1,   \
                             g_0_yyz_0_xxxxyy_0,   \
                             g_0_yyz_0_xxxxyy_1,   \
                             g_0_yyz_0_xxxyyy_0,   \
                             g_0_yyz_0_xxxyyy_1,   \
                             g_0_yyz_0_xxyyyy_0,   \
                             g_0_yyz_0_xxyyyy_1,   \
                             g_0_yyz_0_xyyyyy_0,   \
                             g_0_yyz_0_xyyyyy_1,   \
                             g_0_yyz_0_yyyyyy_0,   \
                             g_0_yyz_0_yyyyyy_1,   \
                             g_0_yyzz_0_xxxxxy_0,  \
                             g_0_yyzz_0_xxxxxy_1,  \
                             g_0_yyzz_0_xxxxyy_0,  \
                             g_0_yyzz_0_xxxxyy_1,  \
                             g_0_yyzz_0_xxxyyy_0,  \
                             g_0_yyzz_0_xxxyyy_1,  \
                             g_0_yyzz_0_xxyyyy_0,  \
                             g_0_yyzz_0_xxyyyy_1,  \
                             g_0_yyzz_0_xyyyyy_0,  \
                             g_0_yyzz_0_xyyyyy_1,  \
                             g_0_yyzz_0_yyyyyy_0,  \
                             g_0_yyzz_0_yyyyyy_1,  \
                             g_0_yyzzz_0_xxxxxx_0, \
                             g_0_yyzzz_0_xxxxxy_0, \
                             g_0_yyzzz_0_xxxxxz_0, \
                             g_0_yyzzz_0_xxxxyy_0, \
                             g_0_yyzzz_0_xxxxyz_0, \
                             g_0_yyzzz_0_xxxxzz_0, \
                             g_0_yyzzz_0_xxxyyy_0, \
                             g_0_yyzzz_0_xxxyyz_0, \
                             g_0_yyzzz_0_xxxyzz_0, \
                             g_0_yyzzz_0_xxxzzz_0, \
                             g_0_yyzzz_0_xxyyyy_0, \
                             g_0_yyzzz_0_xxyyyz_0, \
                             g_0_yyzzz_0_xxyyzz_0, \
                             g_0_yyzzz_0_xxyzzz_0, \
                             g_0_yyzzz_0_xxzzzz_0, \
                             g_0_yyzzz_0_xyyyyy_0, \
                             g_0_yyzzz_0_xyyyyz_0, \
                             g_0_yyzzz_0_xyyyzz_0, \
                             g_0_yyzzz_0_xyyzzz_0, \
                             g_0_yyzzz_0_xyzzzz_0, \
                             g_0_yyzzz_0_xzzzzz_0, \
                             g_0_yyzzz_0_yyyyyy_0, \
                             g_0_yyzzz_0_yyyyyz_0, \
                             g_0_yyzzz_0_yyyyzz_0, \
                             g_0_yyzzz_0_yyyzzz_0, \
                             g_0_yyzzz_0_yyzzzz_0, \
                             g_0_yyzzz_0_yzzzzz_0, \
                             g_0_yyzzz_0_zzzzzz_0, \
                             g_0_yzzz_0_xxxxxx_0,  \
                             g_0_yzzz_0_xxxxxx_1,  \
                             g_0_yzzz_0_xxxxxz_0,  \
                             g_0_yzzz_0_xxxxxz_1,  \
                             g_0_yzzz_0_xxxxyz_0,  \
                             g_0_yzzz_0_xxxxyz_1,  \
                             g_0_yzzz_0_xxxxz_1,   \
                             g_0_yzzz_0_xxxxzz_0,  \
                             g_0_yzzz_0_xxxxzz_1,  \
                             g_0_yzzz_0_xxxyyz_0,  \
                             g_0_yzzz_0_xxxyyz_1,  \
                             g_0_yzzz_0_xxxyz_1,   \
                             g_0_yzzz_0_xxxyzz_0,  \
                             g_0_yzzz_0_xxxyzz_1,  \
                             g_0_yzzz_0_xxxzz_1,   \
                             g_0_yzzz_0_xxxzzz_0,  \
                             g_0_yzzz_0_xxxzzz_1,  \
                             g_0_yzzz_0_xxyyyz_0,  \
                             g_0_yzzz_0_xxyyyz_1,  \
                             g_0_yzzz_0_xxyyz_1,   \
                             g_0_yzzz_0_xxyyzz_0,  \
                             g_0_yzzz_0_xxyyzz_1,  \
                             g_0_yzzz_0_xxyzz_1,   \
                             g_0_yzzz_0_xxyzzz_0,  \
                             g_0_yzzz_0_xxyzzz_1,  \
                             g_0_yzzz_0_xxzzz_1,   \
                             g_0_yzzz_0_xxzzzz_0,  \
                             g_0_yzzz_0_xxzzzz_1,  \
                             g_0_yzzz_0_xyyyyz_0,  \
                             g_0_yzzz_0_xyyyyz_1,  \
                             g_0_yzzz_0_xyyyz_1,   \
                             g_0_yzzz_0_xyyyzz_0,  \
                             g_0_yzzz_0_xyyyzz_1,  \
                             g_0_yzzz_0_xyyzz_1,   \
                             g_0_yzzz_0_xyyzzz_0,  \
                             g_0_yzzz_0_xyyzzz_1,  \
                             g_0_yzzz_0_xyzzz_1,   \
                             g_0_yzzz_0_xyzzzz_0,  \
                             g_0_yzzz_0_xyzzzz_1,  \
                             g_0_yzzz_0_xzzzz_1,   \
                             g_0_yzzz_0_xzzzzz_0,  \
                             g_0_yzzz_0_xzzzzz_1,  \
                             g_0_yzzz_0_yyyyyz_0,  \
                             g_0_yzzz_0_yyyyyz_1,  \
                             g_0_yzzz_0_yyyyz_1,   \
                             g_0_yzzz_0_yyyyzz_0,  \
                             g_0_yzzz_0_yyyyzz_1,  \
                             g_0_yzzz_0_yyyzz_1,   \
                             g_0_yzzz_0_yyyzzz_0,  \
                             g_0_yzzz_0_yyyzzz_1,  \
                             g_0_yzzz_0_yyzzz_1,   \
                             g_0_yzzz_0_yyzzzz_0,  \
                             g_0_yzzz_0_yyzzzz_1,  \
                             g_0_yzzz_0_yzzzz_1,   \
                             g_0_yzzz_0_yzzzzz_0,  \
                             g_0_yzzz_0_yzzzzz_1,  \
                             g_0_yzzz_0_zzzzz_1,   \
                             g_0_yzzz_0_zzzzzz_0,  \
                             g_0_yzzz_0_zzzzzz_1,  \
                             g_0_zzz_0_xxxxxx_0,   \
                             g_0_zzz_0_xxxxxx_1,   \
                             g_0_zzz_0_xxxxxz_0,   \
                             g_0_zzz_0_xxxxxz_1,   \
                             g_0_zzz_0_xxxxyz_0,   \
                             g_0_zzz_0_xxxxyz_1,   \
                             g_0_zzz_0_xxxxzz_0,   \
                             g_0_zzz_0_xxxxzz_1,   \
                             g_0_zzz_0_xxxyyz_0,   \
                             g_0_zzz_0_xxxyyz_1,   \
                             g_0_zzz_0_xxxyzz_0,   \
                             g_0_zzz_0_xxxyzz_1,   \
                             g_0_zzz_0_xxxzzz_0,   \
                             g_0_zzz_0_xxxzzz_1,   \
                             g_0_zzz_0_xxyyyz_0,   \
                             g_0_zzz_0_xxyyyz_1,   \
                             g_0_zzz_0_xxyyzz_0,   \
                             g_0_zzz_0_xxyyzz_1,   \
                             g_0_zzz_0_xxyzzz_0,   \
                             g_0_zzz_0_xxyzzz_1,   \
                             g_0_zzz_0_xxzzzz_0,   \
                             g_0_zzz_0_xxzzzz_1,   \
                             g_0_zzz_0_xyyyyz_0,   \
                             g_0_zzz_0_xyyyyz_1,   \
                             g_0_zzz_0_xyyyzz_0,   \
                             g_0_zzz_0_xyyyzz_1,   \
                             g_0_zzz_0_xyyzzz_0,   \
                             g_0_zzz_0_xyyzzz_1,   \
                             g_0_zzz_0_xyzzzz_0,   \
                             g_0_zzz_0_xyzzzz_1,   \
                             g_0_zzz_0_xzzzzz_0,   \
                             g_0_zzz_0_xzzzzz_1,   \
                             g_0_zzz_0_yyyyyz_0,   \
                             g_0_zzz_0_yyyyyz_1,   \
                             g_0_zzz_0_yyyyzz_0,   \
                             g_0_zzz_0_yyyyzz_1,   \
                             g_0_zzz_0_yyyzzz_0,   \
                             g_0_zzz_0_yyyzzz_1,   \
                             g_0_zzz_0_yyzzzz_0,   \
                             g_0_zzz_0_yyzzzz_1,   \
                             g_0_zzz_0_yzzzzz_0,   \
                             g_0_zzz_0_yzzzzz_1,   \
                             g_0_zzz_0_zzzzzz_0,   \
                             g_0_zzz_0_zzzzzz_1,   \
                             wp_y,                 \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzz_0_xxxxxx_0[i] =
            g_0_zzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxx_0[i] * pb_y + g_0_yzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxy_0[i] = 2.0 * g_0_yyz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxy_0[i] * pb_z +
                                  g_0_yyzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxxz_0[i] =
            g_0_zzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxz_0[i] * pb_y + g_0_yzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxyy_0[i] = 2.0 * g_0_yyz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxyy_0[i] * pb_z +
                                  g_0_yyzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxyz_0[i] = g_0_zzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xxxxyz_0[i] * pb_y + g_0_yzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxzz_0[i] =
            g_0_zzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxzz_0[i] * pb_y + g_0_yzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyyy_0[i] = 2.0 * g_0_yyz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxyyy_0[i] * pb_z +
                                  g_0_yyzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxyyz_0[i] = g_0_zzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxxyz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xxxyyz_0[i] * pb_y + g_0_yzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyzz_0[i] = g_0_zzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xxxyzz_0[i] * pb_y + g_0_yzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxzzz_0[i] =
            g_0_zzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxzzz_0[i] * pb_y + g_0_yzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyyy_0[i] = 2.0 * g_0_yyz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxyyyy_0[i] * pb_z +
                                  g_0_yyzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxyyyz_0[i] = g_0_zzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xxyyz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xxyyyz_0[i] * pb_y + g_0_yzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyzz_0[i] = g_0_zzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxyzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xxyyzz_0[i] * pb_y + g_0_yzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyzzz_0[i] = g_0_zzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xxyzzz_0[i] * pb_y + g_0_yzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxzzzz_0[i] =
            g_0_zzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzzzz_0[i] * pb_y + g_0_yzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyyy_0[i] = 2.0 * g_0_yyz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xyyyyy_0[i] * pb_z +
                                  g_0_yyzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xyyyyz_0[i] = g_0_zzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_xyyyz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xyyyyz_0[i] * pb_y + g_0_yzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyzz_0[i] = g_0_zzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xyyzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xyyyzz_0[i] * pb_y + g_0_yzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyzzz_0[i] = g_0_zzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xyzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xyyzzz_0[i] * pb_y + g_0_yzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyzzzz_0[i] = g_0_zzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_xyzzzz_0[i] * pb_y + g_0_yzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xzzzzz_0[i] =
            g_0_zzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzzzz_0[i] * pb_y + g_0_yzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyyy_0[i] = 2.0 * g_0_yyz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_yyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_yyyyyy_0[i] * pb_z +
                                  g_0_yyzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_yyyyyz_0[i] = g_0_zzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_yyyyyz_0[i] * pb_y + g_0_yzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyzz_0[i] = g_0_zzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_yyyyzz_0[i] * pb_y + g_0_yzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyzzz_0[i] = g_0_zzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_yyyzzz_0[i] * pb_y + g_0_yzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyzzzz_0[i] = g_0_zzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_yyzzzz_0[i] * pb_y + g_0_yzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yzzzzz_0[i] = g_0_zzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzzz_1[i] * fi_abcd_0 +
                                  g_0_yzzz_0_yzzzzz_0[i] * pb_y + g_0_yzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_zzzzzz_0[i] =
            g_0_zzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzzzz_0[i] * pb_y + g_0_yzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 532-560 components of targeted buffer : SHSI

    auto g_0_yzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 532);

    auto g_0_yzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 533);

    auto g_0_yzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 534);

    auto g_0_yzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 535);

    auto g_0_yzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 536);

    auto g_0_yzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 537);

    auto g_0_yzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 538);

    auto g_0_yzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 539);

    auto g_0_yzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 540);

    auto g_0_yzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 541);

    auto g_0_yzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 542);

    auto g_0_yzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 543);

    auto g_0_yzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 544);

    auto g_0_yzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 545);

    auto g_0_yzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 546);

    auto g_0_yzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 547);

    auto g_0_yzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 548);

    auto g_0_yzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 549);

    auto g_0_yzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 550);

    auto g_0_yzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 551);

    auto g_0_yzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 552);

    auto g_0_yzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 553);

    auto g_0_yzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 554);

    auto g_0_yzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 555);

    auto g_0_yzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 556);

    auto g_0_yzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 557);

    auto g_0_yzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 558);

    auto g_0_yzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 559);

#pragma omp simd aligned(g_0_yzzzz_0_xxxxxx_0,     \
                             g_0_yzzzz_0_xxxxxy_0, \
                             g_0_yzzzz_0_xxxxxz_0, \
                             g_0_yzzzz_0_xxxxyy_0, \
                             g_0_yzzzz_0_xxxxyz_0, \
                             g_0_yzzzz_0_xxxxzz_0, \
                             g_0_yzzzz_0_xxxyyy_0, \
                             g_0_yzzzz_0_xxxyyz_0, \
                             g_0_yzzzz_0_xxxyzz_0, \
                             g_0_yzzzz_0_xxxzzz_0, \
                             g_0_yzzzz_0_xxyyyy_0, \
                             g_0_yzzzz_0_xxyyyz_0, \
                             g_0_yzzzz_0_xxyyzz_0, \
                             g_0_yzzzz_0_xxyzzz_0, \
                             g_0_yzzzz_0_xxzzzz_0, \
                             g_0_yzzzz_0_xyyyyy_0, \
                             g_0_yzzzz_0_xyyyyz_0, \
                             g_0_yzzzz_0_xyyyzz_0, \
                             g_0_yzzzz_0_xyyzzz_0, \
                             g_0_yzzzz_0_xyzzzz_0, \
                             g_0_yzzzz_0_xzzzzz_0, \
                             g_0_yzzzz_0_yyyyyy_0, \
                             g_0_yzzzz_0_yyyyyz_0, \
                             g_0_yzzzz_0_yyyyzz_0, \
                             g_0_yzzzz_0_yyyzzz_0, \
                             g_0_yzzzz_0_yyzzzz_0, \
                             g_0_yzzzz_0_yzzzzz_0, \
                             g_0_yzzzz_0_zzzzzz_0, \
                             g_0_zzzz_0_xxxxx_1,   \
                             g_0_zzzz_0_xxxxxx_0,  \
                             g_0_zzzz_0_xxxxxx_1,  \
                             g_0_zzzz_0_xxxxxy_0,  \
                             g_0_zzzz_0_xxxxxy_1,  \
                             g_0_zzzz_0_xxxxxz_0,  \
                             g_0_zzzz_0_xxxxxz_1,  \
                             g_0_zzzz_0_xxxxy_1,   \
                             g_0_zzzz_0_xxxxyy_0,  \
                             g_0_zzzz_0_xxxxyy_1,  \
                             g_0_zzzz_0_xxxxyz_0,  \
                             g_0_zzzz_0_xxxxyz_1,  \
                             g_0_zzzz_0_xxxxz_1,   \
                             g_0_zzzz_0_xxxxzz_0,  \
                             g_0_zzzz_0_xxxxzz_1,  \
                             g_0_zzzz_0_xxxyy_1,   \
                             g_0_zzzz_0_xxxyyy_0,  \
                             g_0_zzzz_0_xxxyyy_1,  \
                             g_0_zzzz_0_xxxyyz_0,  \
                             g_0_zzzz_0_xxxyyz_1,  \
                             g_0_zzzz_0_xxxyz_1,   \
                             g_0_zzzz_0_xxxyzz_0,  \
                             g_0_zzzz_0_xxxyzz_1,  \
                             g_0_zzzz_0_xxxzz_1,   \
                             g_0_zzzz_0_xxxzzz_0,  \
                             g_0_zzzz_0_xxxzzz_1,  \
                             g_0_zzzz_0_xxyyy_1,   \
                             g_0_zzzz_0_xxyyyy_0,  \
                             g_0_zzzz_0_xxyyyy_1,  \
                             g_0_zzzz_0_xxyyyz_0,  \
                             g_0_zzzz_0_xxyyyz_1,  \
                             g_0_zzzz_0_xxyyz_1,   \
                             g_0_zzzz_0_xxyyzz_0,  \
                             g_0_zzzz_0_xxyyzz_1,  \
                             g_0_zzzz_0_xxyzz_1,   \
                             g_0_zzzz_0_xxyzzz_0,  \
                             g_0_zzzz_0_xxyzzz_1,  \
                             g_0_zzzz_0_xxzzz_1,   \
                             g_0_zzzz_0_xxzzzz_0,  \
                             g_0_zzzz_0_xxzzzz_1,  \
                             g_0_zzzz_0_xyyyy_1,   \
                             g_0_zzzz_0_xyyyyy_0,  \
                             g_0_zzzz_0_xyyyyy_1,  \
                             g_0_zzzz_0_xyyyyz_0,  \
                             g_0_zzzz_0_xyyyyz_1,  \
                             g_0_zzzz_0_xyyyz_1,   \
                             g_0_zzzz_0_xyyyzz_0,  \
                             g_0_zzzz_0_xyyyzz_1,  \
                             g_0_zzzz_0_xyyzz_1,   \
                             g_0_zzzz_0_xyyzzz_0,  \
                             g_0_zzzz_0_xyyzzz_1,  \
                             g_0_zzzz_0_xyzzz_1,   \
                             g_0_zzzz_0_xyzzzz_0,  \
                             g_0_zzzz_0_xyzzzz_1,  \
                             g_0_zzzz_0_xzzzz_1,   \
                             g_0_zzzz_0_xzzzzz_0,  \
                             g_0_zzzz_0_xzzzzz_1,  \
                             g_0_zzzz_0_yyyyy_1,   \
                             g_0_zzzz_0_yyyyyy_0,  \
                             g_0_zzzz_0_yyyyyy_1,  \
                             g_0_zzzz_0_yyyyyz_0,  \
                             g_0_zzzz_0_yyyyyz_1,  \
                             g_0_zzzz_0_yyyyz_1,   \
                             g_0_zzzz_0_yyyyzz_0,  \
                             g_0_zzzz_0_yyyyzz_1,  \
                             g_0_zzzz_0_yyyzz_1,   \
                             g_0_zzzz_0_yyyzzz_0,  \
                             g_0_zzzz_0_yyyzzz_1,  \
                             g_0_zzzz_0_yyzzz_1,   \
                             g_0_zzzz_0_yyzzzz_0,  \
                             g_0_zzzz_0_yyzzzz_1,  \
                             g_0_zzzz_0_yzzzz_1,   \
                             g_0_zzzz_0_yzzzzz_0,  \
                             g_0_zzzz_0_yzzzzz_1,  \
                             g_0_zzzz_0_zzzzz_1,   \
                             g_0_zzzz_0_zzzzzz_0,  \
                             g_0_zzzz_0_zzzzzz_1,  \
                             wp_y,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzz_0_xxxxxx_0[i] = g_0_zzzz_0_xxxxxx_0[i] * pb_y + g_0_zzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxy_0[i] = g_0_zzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxy_0[i] * pb_y + g_0_zzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxz_0[i] = g_0_zzzz_0_xxxxxz_0[i] * pb_y + g_0_zzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyy_0[i] = 2.0 * g_0_zzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyy_0[i] * pb_y + g_0_zzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyz_0[i] = g_0_zzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyz_0[i] * pb_y + g_0_zzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxzz_0[i] = g_0_zzzz_0_xxxxzz_0[i] * pb_y + g_0_zzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyy_0[i] = 3.0 * g_0_zzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyy_0[i] * pb_y + g_0_zzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyz_0[i] = 2.0 * g_0_zzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyz_0[i] * pb_y + g_0_zzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyzz_0[i] = g_0_zzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzz_0[i] * pb_y + g_0_zzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxzzz_0[i] = g_0_zzzz_0_xxxzzz_0[i] * pb_y + g_0_zzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyy_0[i] = 4.0 * g_0_zzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyy_0[i] * pb_y + g_0_zzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyz_0[i] = 3.0 * g_0_zzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyz_0[i] * pb_y + g_0_zzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyzz_0[i] = 2.0 * g_0_zzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzz_0[i] * pb_y + g_0_zzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyzzz_0[i] = g_0_zzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzz_0[i] * pb_y + g_0_zzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxzzzz_0[i] = g_0_zzzz_0_xxzzzz_0[i] * pb_y + g_0_zzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyy_0[i] = 5.0 * g_0_zzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyy_0[i] * pb_y + g_0_zzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyz_0[i] = 4.0 * g_0_zzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyz_0[i] * pb_y + g_0_zzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyzz_0[i] = 3.0 * g_0_zzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzz_0[i] * pb_y + g_0_zzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyzzz_0[i] = 2.0 * g_0_zzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzz_0[i] * pb_y + g_0_zzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyzzzz_0[i] = g_0_zzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzz_0[i] * pb_y + g_0_zzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xzzzzz_0[i] = g_0_zzzz_0_xzzzzz_0[i] * pb_y + g_0_zzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyy_0[i] = 6.0 * g_0_zzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyy_0[i] * pb_y + g_0_zzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyz_0[i] = 5.0 * g_0_zzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyz_0[i] * pb_y + g_0_zzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyzz_0[i] = 4.0 * g_0_zzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyzz_0[i] * pb_y + g_0_zzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyzzz_0[i] = 3.0 * g_0_zzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyzzz_0[i] * pb_y + g_0_zzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyzzzz_0[i] = 2.0 * g_0_zzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzzzz_0[i] * pb_y + g_0_zzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yzzzzz_0[i] = g_0_zzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzzzz_0[i] * pb_y + g_0_zzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_zzzzzz_0[i] = g_0_zzzz_0_zzzzzz_0[i] * pb_y + g_0_zzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 560-588 components of targeted buffer : SHSI

    auto g_0_zzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 560);

    auto g_0_zzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 561);

    auto g_0_zzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 562);

    auto g_0_zzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 563);

    auto g_0_zzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 564);

    auto g_0_zzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 565);

    auto g_0_zzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 566);

    auto g_0_zzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 567);

    auto g_0_zzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 568);

    auto g_0_zzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 569);

    auto g_0_zzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 570);

    auto g_0_zzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 571);

    auto g_0_zzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 572);

    auto g_0_zzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 573);

    auto g_0_zzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 574);

    auto g_0_zzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 575);

    auto g_0_zzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 576);

    auto g_0_zzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 577);

    auto g_0_zzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 578);

    auto g_0_zzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 579);

    auto g_0_zzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 580);

    auto g_0_zzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 581);

    auto g_0_zzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 582);

    auto g_0_zzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 583);

    auto g_0_zzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 584);

    auto g_0_zzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 585);

    auto g_0_zzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 586);

    auto g_0_zzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 587);

#pragma omp simd aligned(g_0_zzz_0_xxxxxx_0,       \
                             g_0_zzz_0_xxxxxx_1,   \
                             g_0_zzz_0_xxxxxy_0,   \
                             g_0_zzz_0_xxxxxy_1,   \
                             g_0_zzz_0_xxxxxz_0,   \
                             g_0_zzz_0_xxxxxz_1,   \
                             g_0_zzz_0_xxxxyy_0,   \
                             g_0_zzz_0_xxxxyy_1,   \
                             g_0_zzz_0_xxxxyz_0,   \
                             g_0_zzz_0_xxxxyz_1,   \
                             g_0_zzz_0_xxxxzz_0,   \
                             g_0_zzz_0_xxxxzz_1,   \
                             g_0_zzz_0_xxxyyy_0,   \
                             g_0_zzz_0_xxxyyy_1,   \
                             g_0_zzz_0_xxxyyz_0,   \
                             g_0_zzz_0_xxxyyz_1,   \
                             g_0_zzz_0_xxxyzz_0,   \
                             g_0_zzz_0_xxxyzz_1,   \
                             g_0_zzz_0_xxxzzz_0,   \
                             g_0_zzz_0_xxxzzz_1,   \
                             g_0_zzz_0_xxyyyy_0,   \
                             g_0_zzz_0_xxyyyy_1,   \
                             g_0_zzz_0_xxyyyz_0,   \
                             g_0_zzz_0_xxyyyz_1,   \
                             g_0_zzz_0_xxyyzz_0,   \
                             g_0_zzz_0_xxyyzz_1,   \
                             g_0_zzz_0_xxyzzz_0,   \
                             g_0_zzz_0_xxyzzz_1,   \
                             g_0_zzz_0_xxzzzz_0,   \
                             g_0_zzz_0_xxzzzz_1,   \
                             g_0_zzz_0_xyyyyy_0,   \
                             g_0_zzz_0_xyyyyy_1,   \
                             g_0_zzz_0_xyyyyz_0,   \
                             g_0_zzz_0_xyyyyz_1,   \
                             g_0_zzz_0_xyyyzz_0,   \
                             g_0_zzz_0_xyyyzz_1,   \
                             g_0_zzz_0_xyyzzz_0,   \
                             g_0_zzz_0_xyyzzz_1,   \
                             g_0_zzz_0_xyzzzz_0,   \
                             g_0_zzz_0_xyzzzz_1,   \
                             g_0_zzz_0_xzzzzz_0,   \
                             g_0_zzz_0_xzzzzz_1,   \
                             g_0_zzz_0_yyyyyy_0,   \
                             g_0_zzz_0_yyyyyy_1,   \
                             g_0_zzz_0_yyyyyz_0,   \
                             g_0_zzz_0_yyyyyz_1,   \
                             g_0_zzz_0_yyyyzz_0,   \
                             g_0_zzz_0_yyyyzz_1,   \
                             g_0_zzz_0_yyyzzz_0,   \
                             g_0_zzz_0_yyyzzz_1,   \
                             g_0_zzz_0_yyzzzz_0,   \
                             g_0_zzz_0_yyzzzz_1,   \
                             g_0_zzz_0_yzzzzz_0,   \
                             g_0_zzz_0_yzzzzz_1,   \
                             g_0_zzz_0_zzzzzz_0,   \
                             g_0_zzz_0_zzzzzz_1,   \
                             g_0_zzzz_0_xxxxx_1,   \
                             g_0_zzzz_0_xxxxxx_0,  \
                             g_0_zzzz_0_xxxxxx_1,  \
                             g_0_zzzz_0_xxxxxy_0,  \
                             g_0_zzzz_0_xxxxxy_1,  \
                             g_0_zzzz_0_xxxxxz_0,  \
                             g_0_zzzz_0_xxxxxz_1,  \
                             g_0_zzzz_0_xxxxy_1,   \
                             g_0_zzzz_0_xxxxyy_0,  \
                             g_0_zzzz_0_xxxxyy_1,  \
                             g_0_zzzz_0_xxxxyz_0,  \
                             g_0_zzzz_0_xxxxyz_1,  \
                             g_0_zzzz_0_xxxxz_1,   \
                             g_0_zzzz_0_xxxxzz_0,  \
                             g_0_zzzz_0_xxxxzz_1,  \
                             g_0_zzzz_0_xxxyy_1,   \
                             g_0_zzzz_0_xxxyyy_0,  \
                             g_0_zzzz_0_xxxyyy_1,  \
                             g_0_zzzz_0_xxxyyz_0,  \
                             g_0_zzzz_0_xxxyyz_1,  \
                             g_0_zzzz_0_xxxyz_1,   \
                             g_0_zzzz_0_xxxyzz_0,  \
                             g_0_zzzz_0_xxxyzz_1,  \
                             g_0_zzzz_0_xxxzz_1,   \
                             g_0_zzzz_0_xxxzzz_0,  \
                             g_0_zzzz_0_xxxzzz_1,  \
                             g_0_zzzz_0_xxyyy_1,   \
                             g_0_zzzz_0_xxyyyy_0,  \
                             g_0_zzzz_0_xxyyyy_1,  \
                             g_0_zzzz_0_xxyyyz_0,  \
                             g_0_zzzz_0_xxyyyz_1,  \
                             g_0_zzzz_0_xxyyz_1,   \
                             g_0_zzzz_0_xxyyzz_0,  \
                             g_0_zzzz_0_xxyyzz_1,  \
                             g_0_zzzz_0_xxyzz_1,   \
                             g_0_zzzz_0_xxyzzz_0,  \
                             g_0_zzzz_0_xxyzzz_1,  \
                             g_0_zzzz_0_xxzzz_1,   \
                             g_0_zzzz_0_xxzzzz_0,  \
                             g_0_zzzz_0_xxzzzz_1,  \
                             g_0_zzzz_0_xyyyy_1,   \
                             g_0_zzzz_0_xyyyyy_0,  \
                             g_0_zzzz_0_xyyyyy_1,  \
                             g_0_zzzz_0_xyyyyz_0,  \
                             g_0_zzzz_0_xyyyyz_1,  \
                             g_0_zzzz_0_xyyyz_1,   \
                             g_0_zzzz_0_xyyyzz_0,  \
                             g_0_zzzz_0_xyyyzz_1,  \
                             g_0_zzzz_0_xyyzz_1,   \
                             g_0_zzzz_0_xyyzzz_0,  \
                             g_0_zzzz_0_xyyzzz_1,  \
                             g_0_zzzz_0_xyzzz_1,   \
                             g_0_zzzz_0_xyzzzz_0,  \
                             g_0_zzzz_0_xyzzzz_1,  \
                             g_0_zzzz_0_xzzzz_1,   \
                             g_0_zzzz_0_xzzzzz_0,  \
                             g_0_zzzz_0_xzzzzz_1,  \
                             g_0_zzzz_0_yyyyy_1,   \
                             g_0_zzzz_0_yyyyyy_0,  \
                             g_0_zzzz_0_yyyyyy_1,  \
                             g_0_zzzz_0_yyyyyz_0,  \
                             g_0_zzzz_0_yyyyyz_1,  \
                             g_0_zzzz_0_yyyyz_1,   \
                             g_0_zzzz_0_yyyyzz_0,  \
                             g_0_zzzz_0_yyyyzz_1,  \
                             g_0_zzzz_0_yyyzz_1,   \
                             g_0_zzzz_0_yyyzzz_0,  \
                             g_0_zzzz_0_yyyzzz_1,  \
                             g_0_zzzz_0_yyzzz_1,   \
                             g_0_zzzz_0_yyzzzz_0,  \
                             g_0_zzzz_0_yyzzzz_1,  \
                             g_0_zzzz_0_yzzzz_1,   \
                             g_0_zzzz_0_yzzzzz_0,  \
                             g_0_zzzz_0_yzzzzz_1,  \
                             g_0_zzzz_0_zzzzz_1,   \
                             g_0_zzzz_0_zzzzzz_0,  \
                             g_0_zzzz_0_zzzzzz_1,  \
                             g_0_zzzzz_0_xxxxxx_0, \
                             g_0_zzzzz_0_xxxxxy_0, \
                             g_0_zzzzz_0_xxxxxz_0, \
                             g_0_zzzzz_0_xxxxyy_0, \
                             g_0_zzzzz_0_xxxxyz_0, \
                             g_0_zzzzz_0_xxxxzz_0, \
                             g_0_zzzzz_0_xxxyyy_0, \
                             g_0_zzzzz_0_xxxyyz_0, \
                             g_0_zzzzz_0_xxxyzz_0, \
                             g_0_zzzzz_0_xxxzzz_0, \
                             g_0_zzzzz_0_xxyyyy_0, \
                             g_0_zzzzz_0_xxyyyz_0, \
                             g_0_zzzzz_0_xxyyzz_0, \
                             g_0_zzzzz_0_xxyzzz_0, \
                             g_0_zzzzz_0_xxzzzz_0, \
                             g_0_zzzzz_0_xyyyyy_0, \
                             g_0_zzzzz_0_xyyyyz_0, \
                             g_0_zzzzz_0_xyyyzz_0, \
                             g_0_zzzzz_0_xyyzzz_0, \
                             g_0_zzzzz_0_xyzzzz_0, \
                             g_0_zzzzz_0_xzzzzz_0, \
                             g_0_zzzzz_0_yyyyyy_0, \
                             g_0_zzzzz_0_yyyyyz_0, \
                             g_0_zzzzz_0_yyyyzz_0, \
                             g_0_zzzzz_0_yyyzzz_0, \
                             g_0_zzzzz_0_yyzzzz_0, \
                             g_0_zzzzz_0_yzzzzz_0, \
                             g_0_zzzzz_0_zzzzzz_0, \
                             wp_z,                 \
                             c_exps,               \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzz_0_xxxxxx_0[i] = 4.0 * g_0_zzz_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxx_0[i] * pb_z +
                                  g_0_zzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxy_0[i] = 4.0 * g_0_zzz_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxy_0[i] * pb_z +
                                  g_0_zzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxz_0[i] = 4.0 * g_0_zzz_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxx_1[i] * fi_abcd_0 +
                                  g_0_zzzz_0_xxxxxz_0[i] * pb_z + g_0_zzzz_0_xxxxxz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyy_0[i] = 4.0 * g_0_zzz_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxyy_0[i] * pb_z +
                                  g_0_zzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyz_0[i] = 4.0 * g_0_zzz_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxy_1[i] * fi_abcd_0 +
                                  g_0_zzzz_0_xxxxyz_0[i] * pb_z + g_0_zzzz_0_xxxxyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxzz_0[i] = 4.0 * g_0_zzz_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxzz_0[i] * pb_z + g_0_zzzz_0_xxxxzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyy_0[i] = 4.0 * g_0_zzz_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxyyy_0[i] * pb_z +
                                  g_0_zzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyz_0[i] = 4.0 * g_0_zzz_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxyy_1[i] * fi_abcd_0 +
                                  g_0_zzzz_0_xxxyyz_0[i] * pb_z + g_0_zzzz_0_xxxyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyzz_0[i] = 4.0 * g_0_zzz_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzz_0[i] * pb_z + g_0_zzzz_0_xxxyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxzzz_0[i] = 4.0 * g_0_zzz_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxzzz_0[i] * pb_z + g_0_zzzz_0_xxxzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyy_0[i] = 4.0 * g_0_zzz_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxyyyy_0[i] * pb_z +
                                  g_0_zzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyz_0[i] = 4.0 * g_0_zzz_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxyyy_1[i] * fi_abcd_0 +
                                  g_0_zzzz_0_xxyyyz_0[i] * pb_z + g_0_zzzz_0_xxyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyzz_0[i] = 4.0 * g_0_zzz_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzz_0[i] * pb_z + g_0_zzzz_0_xxyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyzzz_0[i] = 4.0 * g_0_zzz_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzz_0[i] * pb_z + g_0_zzzz_0_xxyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxzzzz_0[i] = 4.0 * g_0_zzz_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzzzz_0[i] * pb_z + g_0_zzzz_0_xxzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyy_0[i] = 4.0 * g_0_zzz_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xyyyyy_0[i] * pb_z +
                                  g_0_zzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyz_0[i] = 4.0 * g_0_zzz_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xyyyy_1[i] * fi_abcd_0 +
                                  g_0_zzzz_0_xyyyyz_0[i] * pb_z + g_0_zzzz_0_xyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyzz_0[i] = 4.0 * g_0_zzz_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzz_0[i] * pb_z + g_0_zzzz_0_xyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyzzz_0[i] = 4.0 * g_0_zzz_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzz_0[i] * pb_z + g_0_zzzz_0_xyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyzzzz_0[i] = 4.0 * g_0_zzz_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzz_0[i] * pb_z + g_0_zzzz_0_xyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xzzzzz_0[i] = 4.0 * g_0_zzz_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xzzzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_zzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzzzz_0[i] * pb_z + g_0_zzzz_0_xzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyy_0[i] = 4.0 * g_0_zzz_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_yyyyyy_0[i] * pb_z +
                                  g_0_zzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyz_0[i] = 4.0 * g_0_zzz_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_yyyyy_1[i] * fi_abcd_0 +
                                  g_0_zzzz_0_yyyyyz_0[i] * pb_z + g_0_zzzz_0_yyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyzz_0[i] = 4.0 * g_0_zzz_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyzz_1[i] * fti_ab_0 +
                                  2.0 * g_0_zzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyzz_0[i] * pb_z + g_0_zzzz_0_yyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyzzz_0[i] = 4.0 * g_0_zzz_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyzzz_1[i] * fti_ab_0 +
                                  3.0 * g_0_zzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyzzz_0[i] * pb_z + g_0_zzzz_0_yyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyzzzz_0[i] = 4.0 * g_0_zzz_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyzzzz_1[i] * fti_ab_0 +
                                  4.0 * g_0_zzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzzzz_0[i] * pb_z + g_0_zzzz_0_yyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yzzzzz_0[i] = 4.0 * g_0_zzz_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yzzzzz_1[i] * fti_ab_0 +
                                  5.0 * g_0_zzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzzzz_0[i] * pb_z + g_0_zzzz_0_yzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_zzzzzz_0[i] = 4.0 * g_0_zzz_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_zzzzzz_1[i] * fti_ab_0 +
                                  6.0 * g_0_zzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_zzzzzz_0[i] * pb_z + g_0_zzzz_0_zzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
