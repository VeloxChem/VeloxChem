#include "NuclearPotentialPrimRecHI.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_hi(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_hi,
                               const size_t              idx_npot_0_fi,
                               const size_t              idx_npot_1_fi,
                               const size_t              idx_npot_0_gh,
                               const size_t              idx_npot_1_gh,
                               const size_t              idx_npot_0_gi,
                               const size_t              idx_npot_1_gi,
                               const CSimdArray<double>& factors,
                               const size_t              idx_rpa,
                               const size_t              idx_rpc,
                               const double              a_exp) -> void
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

    auto ta_xxy_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 30);

    auto ta_xxy_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 33);

    auto ta_xxy_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 37);

    auto ta_xxy_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 42);

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

    auto ta_xxz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 61);

    auto ta_xxz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 62);

    auto ta_xxz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 65);

    auto ta_xxz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 66);

    auto ta_xxz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 70);

    auto ta_xxz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 71);

    auto ta_xxz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 76);

    auto ta_xxz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 78);

    auto ta_xxz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 79);

    auto ta_xxz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 80);

    auto ta_xxz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 81);

    auto ta_xxz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 82);

    auto ta_xxz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 83);

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

    auto ta_yyz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 201);

    auto ta_yyz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 202);

    auto ta_yyz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 205);

    auto ta_yyz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 206);

    auto ta_yyz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 210);

    auto ta_yyz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 211);

    auto ta_yyz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 216);

    auto ta_yyz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 217);

    auto ta_yyz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 218);

    auto ta_yyz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 219);

    auto ta_yyz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 220);

    auto ta_yyz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 221);

    auto ta_yyz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 222);

    auto ta_yyz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 223);

    auto ta_yzz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 224);

    auto ta_yzz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 226);

    auto ta_yzz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 228);

    auto ta_yzz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 229);

    auto ta_yzz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 231);

    auto ta_yzz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 232);

    auto ta_yzz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 233);

    auto ta_yzz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 235);

    auto ta_yzz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 236);

    auto ta_yzz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 237);

    auto ta_yzz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 238);

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

    auto ta_xxy_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 30);

    auto ta_xxy_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 33);

    auto ta_xxy_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 37);

    auto ta_xxy_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 42);

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

    auto ta_xxz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 61);

    auto ta_xxz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 62);

    auto ta_xxz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 65);

    auto ta_xxz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 66);

    auto ta_xxz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 70);

    auto ta_xxz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 71);

    auto ta_xxz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 76);

    auto ta_xxz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 78);

    auto ta_xxz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 79);

    auto ta_xxz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 80);

    auto ta_xxz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 81);

    auto ta_xxz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 82);

    auto ta_xxz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 83);

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

    auto ta_yyz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 201);

    auto ta_yyz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 202);

    auto ta_yyz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 205);

    auto ta_yyz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 206);

    auto ta_yyz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 210);

    auto ta_yyz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 211);

    auto ta_yyz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 216);

    auto ta_yyz_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 217);

    auto ta_yyz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 218);

    auto ta_yyz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 219);

    auto ta_yyz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 220);

    auto ta_yyz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 221);

    auto ta_yyz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 222);

    auto ta_yyz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 223);

    auto ta_yzz_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 224);

    auto ta_yzz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 226);

    auto ta_yzz_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 228);

    auto ta_yzz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 229);

    auto ta_yzz_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 231);

    auto ta_yzz_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 232);

    auto ta_yzz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 233);

    auto ta_yzz_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 235);

    auto ta_yzz_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 236);

    auto ta_yzz_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 237);

    auto ta_yzz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 238);

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

    // Set up components of auxiliary buffer : GH

    auto ta_xxxx_xxxxx_0 = pbuffer.data(idx_npot_0_gh);

    auto ta_xxxx_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 1);

    auto ta_xxxx_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 2);

    auto ta_xxxx_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 3);

    auto ta_xxxx_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 4);

    auto ta_xxxx_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 5);

    auto ta_xxxx_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 6);

    auto ta_xxxx_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 7);

    auto ta_xxxx_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 8);

    auto ta_xxxx_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 9);

    auto ta_xxxx_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 10);

    auto ta_xxxx_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 11);

    auto ta_xxxx_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 12);

    auto ta_xxxx_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 13);

    auto ta_xxxx_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 14);

    auto ta_xxxx_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 15);

    auto ta_xxxx_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 16);

    auto ta_xxxx_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 17);

    auto ta_xxxx_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 18);

    auto ta_xxxx_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 19);

    auto ta_xxxx_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 20);

    auto ta_xxxz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 44);

    auto ta_xxxz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 46);

    auto ta_xxxz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 47);

    auto ta_xxxz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 49);

    auto ta_xxxz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 50);

    auto ta_xxxz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 51);

    auto ta_xxxz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 53);

    auto ta_xxxz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 54);

    auto ta_xxxz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 55);

    auto ta_xxxz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 56);

    auto ta_xxyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 64);

    auto ta_xxyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 66);

    auto ta_xxyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 67);

    auto ta_xxyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 69);

    auto ta_xxyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 70);

    auto ta_xxyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 71);

    auto ta_xxyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 73);

    auto ta_xxyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 74);

    auto ta_xxyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 75);

    auto ta_xxyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 76);

    auto ta_xxyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 78);

    auto ta_xxyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 79);

    auto ta_xxyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 80);

    auto ta_xxyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 81);

    auto ta_xxyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 82);

    auto ta_xxzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 105);

    auto ta_xxzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 106);

    auto ta_xxzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 107);

    auto ta_xxzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 108);

    auto ta_xxzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 109);

    auto ta_xxzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 110);

    auto ta_xxzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 111);

    auto ta_xxzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 112);

    auto ta_xxzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 113);

    auto ta_xxzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 114);

    auto ta_xxzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 115);

    auto ta_xxzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 116);

    auto ta_xxzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 117);

    auto ta_xxzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 118);

    auto ta_xxzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 119);

    auto ta_xxzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 121);

    auto ta_xxzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 122);

    auto ta_xxzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 123);

    auto ta_xxzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 124);

    auto ta_xxzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 125);

    auto ta_xyyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 127);

    auto ta_xyyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 129);

    auto ta_xyyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 130);

    auto ta_xyyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 132);

    auto ta_xyyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 133);

    auto ta_xyyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 134);

    auto ta_xyyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 136);

    auto ta_xyyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 137);

    auto ta_xyyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 138);

    auto ta_xyyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 139);

    auto ta_xyyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 141);

    auto ta_xyyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 142);

    auto ta_xyyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 143);

    auto ta_xyyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 144);

    auto ta_xyyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 145);

    auto ta_xzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 191);

    auto ta_xzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 193);

    auto ta_xzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 194);

    auto ta_xzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 196);

    auto ta_xzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 197);

    auto ta_xzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 198);

    auto ta_xzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 200);

    auto ta_xzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 201);

    auto ta_xzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 202);

    auto ta_xzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 203);

    auto ta_xzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 205);

    auto ta_xzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 206);

    auto ta_xzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 207);

    auto ta_xzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 208);

    auto ta_xzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 209);

    auto ta_yyyy_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 210);

    auto ta_yyyy_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 211);

    auto ta_yyyy_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 212);

    auto ta_yyyy_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 213);

    auto ta_yyyy_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 214);

    auto ta_yyyy_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 215);

    auto ta_yyyy_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 216);

    auto ta_yyyy_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 217);

    auto ta_yyyy_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 218);

    auto ta_yyyy_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 219);

    auto ta_yyyy_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 220);

    auto ta_yyyy_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 221);

    auto ta_yyyy_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 222);

    auto ta_yyyy_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 223);

    auto ta_yyyy_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 224);

    auto ta_yyyy_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 225);

    auto ta_yyyy_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 226);

    auto ta_yyyy_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 227);

    auto ta_yyyy_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 228);

    auto ta_yyyy_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 229);

    auto ta_yyyy_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 230);

    auto ta_yyyz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 233);

    auto ta_yyyz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 235);

    auto ta_yyyz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 236);

    auto ta_yyyz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 238);

    auto ta_yyyz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 239);

    auto ta_yyyz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 240);

    auto ta_yyyz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 242);

    auto ta_yyyz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 243);

    auto ta_yyyz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 244);

    auto ta_yyyz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 245);

    auto ta_yyyz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 247);

    auto ta_yyyz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 248);

    auto ta_yyyz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 249);

    auto ta_yyyz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 250);

    auto ta_yyyz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 251);

    auto ta_yyzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 252);

    auto ta_yyzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 253);

    auto ta_yyzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 254);

    auto ta_yyzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 255);

    auto ta_yyzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 256);

    auto ta_yyzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 257);

    auto ta_yyzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 258);

    auto ta_yyzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 259);

    auto ta_yyzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 260);

    auto ta_yyzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 261);

    auto ta_yyzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 262);

    auto ta_yyzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 263);

    auto ta_yyzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 264);

    auto ta_yyzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 265);

    auto ta_yyzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 266);

    auto ta_yyzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 267);

    auto ta_yyzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 268);

    auto ta_yyzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 269);

    auto ta_yyzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 270);

    auto ta_yyzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 271);

    auto ta_yyzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 272);

    auto ta_yzzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 274);

    auto ta_yzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 275);

    auto ta_yzzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 276);

    auto ta_yzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 277);

    auto ta_yzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 278);

    auto ta_yzzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 279);

    auto ta_yzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 280);

    auto ta_yzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 281);

    auto ta_yzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 282);

    auto ta_yzzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 283);

    auto ta_yzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 284);

    auto ta_yzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 285);

    auto ta_yzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 286);

    auto ta_yzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 287);

    auto ta_yzzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 288);

    auto ta_yzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 289);

    auto ta_yzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 290);

    auto ta_yzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 291);

    auto ta_yzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 292);

    auto ta_yzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 293);

    auto ta_zzzz_xxxxx_0 = pbuffer.data(idx_npot_0_gh + 294);

    auto ta_zzzz_xxxxy_0 = pbuffer.data(idx_npot_0_gh + 295);

    auto ta_zzzz_xxxxz_0 = pbuffer.data(idx_npot_0_gh + 296);

    auto ta_zzzz_xxxyy_0 = pbuffer.data(idx_npot_0_gh + 297);

    auto ta_zzzz_xxxyz_0 = pbuffer.data(idx_npot_0_gh + 298);

    auto ta_zzzz_xxxzz_0 = pbuffer.data(idx_npot_0_gh + 299);

    auto ta_zzzz_xxyyy_0 = pbuffer.data(idx_npot_0_gh + 300);

    auto ta_zzzz_xxyyz_0 = pbuffer.data(idx_npot_0_gh + 301);

    auto ta_zzzz_xxyzz_0 = pbuffer.data(idx_npot_0_gh + 302);

    auto ta_zzzz_xxzzz_0 = pbuffer.data(idx_npot_0_gh + 303);

    auto ta_zzzz_xyyyy_0 = pbuffer.data(idx_npot_0_gh + 304);

    auto ta_zzzz_xyyyz_0 = pbuffer.data(idx_npot_0_gh + 305);

    auto ta_zzzz_xyyzz_0 = pbuffer.data(idx_npot_0_gh + 306);

    auto ta_zzzz_xyzzz_0 = pbuffer.data(idx_npot_0_gh + 307);

    auto ta_zzzz_xzzzz_0 = pbuffer.data(idx_npot_0_gh + 308);

    auto ta_zzzz_yyyyy_0 = pbuffer.data(idx_npot_0_gh + 309);

    auto ta_zzzz_yyyyz_0 = pbuffer.data(idx_npot_0_gh + 310);

    auto ta_zzzz_yyyzz_0 = pbuffer.data(idx_npot_0_gh + 311);

    auto ta_zzzz_yyzzz_0 = pbuffer.data(idx_npot_0_gh + 312);

    auto ta_zzzz_yzzzz_0 = pbuffer.data(idx_npot_0_gh + 313);

    auto ta_zzzz_zzzzz_0 = pbuffer.data(idx_npot_0_gh + 314);

    // Set up components of auxiliary buffer : GH

    auto ta_xxxx_xxxxx_1 = pbuffer.data(idx_npot_1_gh);

    auto ta_xxxx_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 1);

    auto ta_xxxx_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 2);

    auto ta_xxxx_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 3);

    auto ta_xxxx_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 4);

    auto ta_xxxx_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 5);

    auto ta_xxxx_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 6);

    auto ta_xxxx_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 7);

    auto ta_xxxx_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 8);

    auto ta_xxxx_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 9);

    auto ta_xxxx_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 10);

    auto ta_xxxx_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 11);

    auto ta_xxxx_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 12);

    auto ta_xxxx_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 13);

    auto ta_xxxx_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 14);

    auto ta_xxxx_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 15);

    auto ta_xxxx_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 16);

    auto ta_xxxx_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 17);

    auto ta_xxxx_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 18);

    auto ta_xxxx_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 19);

    auto ta_xxxx_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 20);

    auto ta_xxxz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 44);

    auto ta_xxxz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 46);

    auto ta_xxxz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 47);

    auto ta_xxxz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 49);

    auto ta_xxxz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 50);

    auto ta_xxxz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 51);

    auto ta_xxxz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 53);

    auto ta_xxxz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 54);

    auto ta_xxxz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 55);

    auto ta_xxxz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 56);

    auto ta_xxyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 64);

    auto ta_xxyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 66);

    auto ta_xxyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 67);

    auto ta_xxyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 69);

    auto ta_xxyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 70);

    auto ta_xxyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 71);

    auto ta_xxyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 73);

    auto ta_xxyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 74);

    auto ta_xxyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 75);

    auto ta_xxyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 76);

    auto ta_xxyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 78);

    auto ta_xxyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 79);

    auto ta_xxyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 80);

    auto ta_xxyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 81);

    auto ta_xxyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 82);

    auto ta_xxzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 105);

    auto ta_xxzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 106);

    auto ta_xxzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 107);

    auto ta_xxzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 108);

    auto ta_xxzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 109);

    auto ta_xxzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 110);

    auto ta_xxzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 111);

    auto ta_xxzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 112);

    auto ta_xxzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 113);

    auto ta_xxzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 114);

    auto ta_xxzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 115);

    auto ta_xxzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 116);

    auto ta_xxzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 117);

    auto ta_xxzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 118);

    auto ta_xxzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 119);

    auto ta_xxzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 121);

    auto ta_xxzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 122);

    auto ta_xxzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 123);

    auto ta_xxzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 124);

    auto ta_xxzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 125);

    auto ta_xyyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 127);

    auto ta_xyyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 129);

    auto ta_xyyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 130);

    auto ta_xyyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 132);

    auto ta_xyyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 133);

    auto ta_xyyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 134);

    auto ta_xyyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 136);

    auto ta_xyyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 137);

    auto ta_xyyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 138);

    auto ta_xyyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 139);

    auto ta_xyyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 141);

    auto ta_xyyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 142);

    auto ta_xyyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 143);

    auto ta_xyyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 144);

    auto ta_xyyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 145);

    auto ta_xzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 191);

    auto ta_xzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 193);

    auto ta_xzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 194);

    auto ta_xzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 196);

    auto ta_xzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 197);

    auto ta_xzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 198);

    auto ta_xzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 200);

    auto ta_xzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 201);

    auto ta_xzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 202);

    auto ta_xzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 203);

    auto ta_xzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 205);

    auto ta_xzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 206);

    auto ta_xzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 207);

    auto ta_xzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 208);

    auto ta_xzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 209);

    auto ta_yyyy_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 210);

    auto ta_yyyy_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 211);

    auto ta_yyyy_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 212);

    auto ta_yyyy_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 213);

    auto ta_yyyy_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 214);

    auto ta_yyyy_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 215);

    auto ta_yyyy_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 216);

    auto ta_yyyy_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 217);

    auto ta_yyyy_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 218);

    auto ta_yyyy_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 219);

    auto ta_yyyy_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 220);

    auto ta_yyyy_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 221);

    auto ta_yyyy_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 222);

    auto ta_yyyy_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 223);

    auto ta_yyyy_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 224);

    auto ta_yyyy_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 225);

    auto ta_yyyy_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 226);

    auto ta_yyyy_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 227);

    auto ta_yyyy_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 228);

    auto ta_yyyy_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 229);

    auto ta_yyyy_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 230);

    auto ta_yyyz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 233);

    auto ta_yyyz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 235);

    auto ta_yyyz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 236);

    auto ta_yyyz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 238);

    auto ta_yyyz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 239);

    auto ta_yyyz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 240);

    auto ta_yyyz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 242);

    auto ta_yyyz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 243);

    auto ta_yyyz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 244);

    auto ta_yyyz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 245);

    auto ta_yyyz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 247);

    auto ta_yyyz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 248);

    auto ta_yyyz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 249);

    auto ta_yyyz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 250);

    auto ta_yyyz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 251);

    auto ta_yyzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 252);

    auto ta_yyzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 253);

    auto ta_yyzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 254);

    auto ta_yyzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 255);

    auto ta_yyzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 256);

    auto ta_yyzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 257);

    auto ta_yyzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 258);

    auto ta_yyzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 259);

    auto ta_yyzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 260);

    auto ta_yyzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 261);

    auto ta_yyzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 262);

    auto ta_yyzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 263);

    auto ta_yyzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 264);

    auto ta_yyzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 265);

    auto ta_yyzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 266);

    auto ta_yyzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 267);

    auto ta_yyzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 268);

    auto ta_yyzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 269);

    auto ta_yyzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 270);

    auto ta_yyzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 271);

    auto ta_yyzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 272);

    auto ta_yzzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 274);

    auto ta_yzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 275);

    auto ta_yzzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 276);

    auto ta_yzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 277);

    auto ta_yzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 278);

    auto ta_yzzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 279);

    auto ta_yzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 280);

    auto ta_yzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 281);

    auto ta_yzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 282);

    auto ta_yzzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 283);

    auto ta_yzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 284);

    auto ta_yzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 285);

    auto ta_yzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 286);

    auto ta_yzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 287);

    auto ta_yzzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 288);

    auto ta_yzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 289);

    auto ta_yzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 290);

    auto ta_yzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 291);

    auto ta_yzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 292);

    auto ta_yzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 293);

    auto ta_zzzz_xxxxx_1 = pbuffer.data(idx_npot_1_gh + 294);

    auto ta_zzzz_xxxxy_1 = pbuffer.data(idx_npot_1_gh + 295);

    auto ta_zzzz_xxxxz_1 = pbuffer.data(idx_npot_1_gh + 296);

    auto ta_zzzz_xxxyy_1 = pbuffer.data(idx_npot_1_gh + 297);

    auto ta_zzzz_xxxyz_1 = pbuffer.data(idx_npot_1_gh + 298);

    auto ta_zzzz_xxxzz_1 = pbuffer.data(idx_npot_1_gh + 299);

    auto ta_zzzz_xxyyy_1 = pbuffer.data(idx_npot_1_gh + 300);

    auto ta_zzzz_xxyyz_1 = pbuffer.data(idx_npot_1_gh + 301);

    auto ta_zzzz_xxyzz_1 = pbuffer.data(idx_npot_1_gh + 302);

    auto ta_zzzz_xxzzz_1 = pbuffer.data(idx_npot_1_gh + 303);

    auto ta_zzzz_xyyyy_1 = pbuffer.data(idx_npot_1_gh + 304);

    auto ta_zzzz_xyyyz_1 = pbuffer.data(idx_npot_1_gh + 305);

    auto ta_zzzz_xyyzz_1 = pbuffer.data(idx_npot_1_gh + 306);

    auto ta_zzzz_xyzzz_1 = pbuffer.data(idx_npot_1_gh + 307);

    auto ta_zzzz_xzzzz_1 = pbuffer.data(idx_npot_1_gh + 308);

    auto ta_zzzz_yyyyy_1 = pbuffer.data(idx_npot_1_gh + 309);

    auto ta_zzzz_yyyyz_1 = pbuffer.data(idx_npot_1_gh + 310);

    auto ta_zzzz_yyyzz_1 = pbuffer.data(idx_npot_1_gh + 311);

    auto ta_zzzz_yyzzz_1 = pbuffer.data(idx_npot_1_gh + 312);

    auto ta_zzzz_yzzzz_1 = pbuffer.data(idx_npot_1_gh + 313);

    auto ta_zzzz_zzzzz_1 = pbuffer.data(idx_npot_1_gh + 314);

    // Set up components of auxiliary buffer : GI

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

    auto ta_xxxy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 28);

    auto ta_xxxy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 29);

    auto ta_xxxy_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 30);

    auto ta_xxxy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 31);

    auto ta_xxxy_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 33);

    auto ta_xxxy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 34);

    auto ta_xxxy_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 37);

    auto ta_xxxy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 38);

    auto ta_xxxy_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 42);

    auto ta_xxxy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 43);

    auto ta_xxxy_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 48);

    auto ta_xxxy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 49);

    auto ta_xxxy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 50);

    auto ta_xxxy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 51);

    auto ta_xxxy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 52);

    auto ta_xxxy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 53);

    auto ta_xxxy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 54);

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

    auto ta_xxxz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 78);

    auto ta_xxxz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 79);

    auto ta_xxxz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 80);

    auto ta_xxxz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 81);

    auto ta_xxxz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 82);

    auto ta_xxxz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 83);

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

    auto ta_xxyz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 114);

    auto ta_xxyz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 117);

    auto ta_xxyz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 121);

    auto ta_xxyz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 126);

    auto ta_xxyz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 132);

    auto ta_xxyz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 134);

    auto ta_xxyz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 135);

    auto ta_xxyz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 136);

    auto ta_xxyz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 137);

    auto ta_xxyz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 138);

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

    auto ta_xyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 168);

    auto ta_xyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 169);

    auto ta_xyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 171);

    auto ta_xyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 172);

    auto ta_xyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 174);

    auto ta_xyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 175);

    auto ta_xyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 176);

    auto ta_xyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 178);

    auto ta_xyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 179);

    auto ta_xyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 180);

    auto ta_xyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 181);

    auto ta_xyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 183);

    auto ta_xyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 184);

    auto ta_xyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 185);

    auto ta_xyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 186);

    auto ta_xyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 187);

    auto ta_xyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 189);

    auto ta_xyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 190);

    auto ta_xyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 191);

    auto ta_xyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 192);

    auto ta_xyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 193);

    auto ta_xyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 194);

    auto ta_xyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 195);

    auto ta_xyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 218);

    auto ta_xyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 219);

    auto ta_xyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 220);

    auto ta_xyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 221);

    auto ta_xyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 222);

    auto ta_xyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 223);

    auto ta_xyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 245);

    auto ta_xyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 246);

    auto ta_xyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 247);

    auto ta_xyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 248);

    auto ta_xyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 249);

    auto ta_xyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 250);

    auto ta_xzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 252);

    auto ta_xzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 254);

    auto ta_xzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 256);

    auto ta_xzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 257);

    auto ta_xzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 259);

    auto ta_xzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 260);

    auto ta_xzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 261);

    auto ta_xzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 263);

    auto ta_xzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 264);

    auto ta_xzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 265);

    auto ta_xzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 266);

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

    // Set up components of auxiliary buffer : GI

    auto ta_xxxx_xxxxxx_1 = pbuffer.data(idx_npot_1_gi);

    auto ta_xxxx_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 1);

    auto ta_xxxx_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 2);

    auto ta_xxxx_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 3);

    auto ta_xxxx_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 4);

    auto ta_xxxx_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 5);

    auto ta_xxxx_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 6);

    auto ta_xxxx_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 7);

    auto ta_xxxx_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 8);

    auto ta_xxxx_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 9);

    auto ta_xxxx_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 10);

    auto ta_xxxx_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 11);

    auto ta_xxxx_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 12);

    auto ta_xxxx_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 13);

    auto ta_xxxx_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 14);

    auto ta_xxxx_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 15);

    auto ta_xxxx_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 16);

    auto ta_xxxx_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 17);

    auto ta_xxxx_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 18);

    auto ta_xxxx_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 19);

    auto ta_xxxx_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 20);

    auto ta_xxxx_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 21);

    auto ta_xxxx_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 22);

    auto ta_xxxx_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 23);

    auto ta_xxxx_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 24);

    auto ta_xxxx_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 25);

    auto ta_xxxx_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 26);

    auto ta_xxxx_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 27);

    auto ta_xxxy_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 28);

    auto ta_xxxy_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 29);

    auto ta_xxxy_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 30);

    auto ta_xxxy_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 31);

    auto ta_xxxy_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 33);

    auto ta_xxxy_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 34);

    auto ta_xxxy_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 37);

    auto ta_xxxy_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 38);

    auto ta_xxxy_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 42);

    auto ta_xxxy_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 43);

    auto ta_xxxy_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 48);

    auto ta_xxxy_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 49);

    auto ta_xxxy_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 50);

    auto ta_xxxy_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 51);

    auto ta_xxxy_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 52);

    auto ta_xxxy_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 53);

    auto ta_xxxy_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 54);

    auto ta_xxxz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 56);

    auto ta_xxxz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 57);

    auto ta_xxxz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 58);

    auto ta_xxxz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 59);

    auto ta_xxxz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 60);

    auto ta_xxxz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 61);

    auto ta_xxxz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 62);

    auto ta_xxxz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 63);

    auto ta_xxxz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 64);

    auto ta_xxxz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 65);

    auto ta_xxxz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 66);

    auto ta_xxxz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 67);

    auto ta_xxxz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 68);

    auto ta_xxxz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 69);

    auto ta_xxxz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 70);

    auto ta_xxxz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 71);

    auto ta_xxxz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 72);

    auto ta_xxxz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 73);

    auto ta_xxxz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 74);

    auto ta_xxxz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 75);

    auto ta_xxxz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 76);

    auto ta_xxxz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 78);

    auto ta_xxxz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 79);

    auto ta_xxxz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 80);

    auto ta_xxxz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 81);

    auto ta_xxxz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 82);

    auto ta_xxxz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 83);

    auto ta_xxyy_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 84);

    auto ta_xxyy_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 85);

    auto ta_xxyy_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 86);

    auto ta_xxyy_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 87);

    auto ta_xxyy_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 88);

    auto ta_xxyy_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 89);

    auto ta_xxyy_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 90);

    auto ta_xxyy_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 91);

    auto ta_xxyy_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 92);

    auto ta_xxyy_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 93);

    auto ta_xxyy_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 94);

    auto ta_xxyy_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 95);

    auto ta_xxyy_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 96);

    auto ta_xxyy_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 97);

    auto ta_xxyy_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 98);

    auto ta_xxyy_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 99);

    auto ta_xxyy_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 100);

    auto ta_xxyy_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 101);

    auto ta_xxyy_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 102);

    auto ta_xxyy_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 103);

    auto ta_xxyy_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 104);

    auto ta_xxyy_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 105);

    auto ta_xxyy_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 106);

    auto ta_xxyy_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 107);

    auto ta_xxyy_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 108);

    auto ta_xxyy_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 109);

    auto ta_xxyy_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 110);

    auto ta_xxyy_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 111);

    auto ta_xxyz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 114);

    auto ta_xxyz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 117);

    auto ta_xxyz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 121);

    auto ta_xxyz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 126);

    auto ta_xxyz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 132);

    auto ta_xxyz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 134);

    auto ta_xxyz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 135);

    auto ta_xxyz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 136);

    auto ta_xxyz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 137);

    auto ta_xxyz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 138);

    auto ta_xxzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 140);

    auto ta_xxzz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 141);

    auto ta_xxzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 142);

    auto ta_xxzz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 143);

    auto ta_xxzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 144);

    auto ta_xxzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 145);

    auto ta_xxzz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 146);

    auto ta_xxzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 147);

    auto ta_xxzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 148);

    auto ta_xxzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 149);

    auto ta_xxzz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 150);

    auto ta_xxzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 151);

    auto ta_xxzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 152);

    auto ta_xxzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 153);

    auto ta_xxzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 154);

    auto ta_xxzz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 155);

    auto ta_xxzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 156);

    auto ta_xxzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 157);

    auto ta_xxzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 158);

    auto ta_xxzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 159);

    auto ta_xxzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 160);

    auto ta_xxzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 161);

    auto ta_xxzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 162);

    auto ta_xxzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 163);

    auto ta_xxzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 164);

    auto ta_xxzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 165);

    auto ta_xxzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 166);

    auto ta_xxzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 167);

    auto ta_xyyy_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 168);

    auto ta_xyyy_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 169);

    auto ta_xyyy_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 171);

    auto ta_xyyy_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 172);

    auto ta_xyyy_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 174);

    auto ta_xyyy_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 175);

    auto ta_xyyy_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 176);

    auto ta_xyyy_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 178);

    auto ta_xyyy_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 179);

    auto ta_xyyy_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 180);

    auto ta_xyyy_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 181);

    auto ta_xyyy_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 183);

    auto ta_xyyy_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 184);

    auto ta_xyyy_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 185);

    auto ta_xyyy_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 186);

    auto ta_xyyy_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 187);

    auto ta_xyyy_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 189);

    auto ta_xyyy_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 190);

    auto ta_xyyy_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 191);

    auto ta_xyyy_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 192);

    auto ta_xyyy_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 193);

    auto ta_xyyy_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 194);

    auto ta_xyyy_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 195);

    auto ta_xyyz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 218);

    auto ta_xyyz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 219);

    auto ta_xyyz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 220);

    auto ta_xyyz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 221);

    auto ta_xyyz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 222);

    auto ta_xyyz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 223);

    auto ta_xyzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 245);

    auto ta_xyzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 246);

    auto ta_xyzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 247);

    auto ta_xyzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 248);

    auto ta_xyzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 249);

    auto ta_xyzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 250);

    auto ta_xzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 252);

    auto ta_xzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 254);

    auto ta_xzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 256);

    auto ta_xzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 257);

    auto ta_xzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 259);

    auto ta_xzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 260);

    auto ta_xzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 261);

    auto ta_xzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 263);

    auto ta_xzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 264);

    auto ta_xzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 265);

    auto ta_xzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 266);

    auto ta_xzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 268);

    auto ta_xzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 269);

    auto ta_xzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 270);

    auto ta_xzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 271);

    auto ta_xzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 272);

    auto ta_xzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 273);

    auto ta_xzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 274);

    auto ta_xzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 275);

    auto ta_xzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 276);

    auto ta_xzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 277);

    auto ta_xzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 278);

    auto ta_xzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 279);

    auto ta_yyyy_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 280);

    auto ta_yyyy_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 281);

    auto ta_yyyy_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 282);

    auto ta_yyyy_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 283);

    auto ta_yyyy_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 284);

    auto ta_yyyy_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 285);

    auto ta_yyyy_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 286);

    auto ta_yyyy_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 287);

    auto ta_yyyy_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 288);

    auto ta_yyyy_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 289);

    auto ta_yyyy_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 290);

    auto ta_yyyy_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 291);

    auto ta_yyyy_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 292);

    auto ta_yyyy_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 293);

    auto ta_yyyy_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 294);

    auto ta_yyyy_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 295);

    auto ta_yyyy_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 296);

    auto ta_yyyy_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 297);

    auto ta_yyyy_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 298);

    auto ta_yyyy_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 299);

    auto ta_yyyy_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 300);

    auto ta_yyyy_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 301);

    auto ta_yyyy_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 302);

    auto ta_yyyy_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 303);

    auto ta_yyyy_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 304);

    auto ta_yyyy_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 305);

    auto ta_yyyy_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 306);

    auto ta_yyyy_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 307);

    auto ta_yyyz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 309);

    auto ta_yyyz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 310);

    auto ta_yyyz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 311);

    auto ta_yyyz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 312);

    auto ta_yyyz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 313);

    auto ta_yyyz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 314);

    auto ta_yyyz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 315);

    auto ta_yyyz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 316);

    auto ta_yyyz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 317);

    auto ta_yyyz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 318);

    auto ta_yyyz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 319);

    auto ta_yyyz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 320);

    auto ta_yyyz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 321);

    auto ta_yyyz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 322);

    auto ta_yyyz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 323);

    auto ta_yyyz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 324);

    auto ta_yyyz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 325);

    auto ta_yyyz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 326);

    auto ta_yyyz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 327);

    auto ta_yyyz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 328);

    auto ta_yyyz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 329);

    auto ta_yyyz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 330);

    auto ta_yyyz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 331);

    auto ta_yyyz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 332);

    auto ta_yyyz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 333);

    auto ta_yyyz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 334);

    auto ta_yyyz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 335);

    auto ta_yyzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 336);

    auto ta_yyzz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 337);

    auto ta_yyzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 338);

    auto ta_yyzz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 339);

    auto ta_yyzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 340);

    auto ta_yyzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 341);

    auto ta_yyzz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 342);

    auto ta_yyzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 343);

    auto ta_yyzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 344);

    auto ta_yyzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 345);

    auto ta_yyzz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 346);

    auto ta_yyzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 347);

    auto ta_yyzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 348);

    auto ta_yyzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 349);

    auto ta_yyzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 350);

    auto ta_yyzz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 351);

    auto ta_yyzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 352);

    auto ta_yyzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 353);

    auto ta_yyzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 354);

    auto ta_yyzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 355);

    auto ta_yyzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 356);

    auto ta_yyzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 357);

    auto ta_yyzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 358);

    auto ta_yyzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 359);

    auto ta_yyzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 360);

    auto ta_yyzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 361);

    auto ta_yyzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 362);

    auto ta_yyzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 363);

    auto ta_yzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 364);

    auto ta_yzzz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 365);

    auto ta_yzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 366);

    auto ta_yzzz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 367);

    auto ta_yzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 368);

    auto ta_yzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 369);

    auto ta_yzzz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 370);

    auto ta_yzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 371);

    auto ta_yzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 372);

    auto ta_yzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 373);

    auto ta_yzzz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 374);

    auto ta_yzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 375);

    auto ta_yzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 376);

    auto ta_yzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 377);

    auto ta_yzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 378);

    auto ta_yzzz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 379);

    auto ta_yzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 380);

    auto ta_yzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 381);

    auto ta_yzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 382);

    auto ta_yzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 383);

    auto ta_yzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 384);

    auto ta_yzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 385);

    auto ta_yzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 386);

    auto ta_yzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 387);

    auto ta_yzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 388);

    auto ta_yzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 389);

    auto ta_yzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 390);

    auto ta_yzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 391);

    auto ta_zzzz_xxxxxx_1 = pbuffer.data(idx_npot_1_gi + 392);

    auto ta_zzzz_xxxxxy_1 = pbuffer.data(idx_npot_1_gi + 393);

    auto ta_zzzz_xxxxxz_1 = pbuffer.data(idx_npot_1_gi + 394);

    auto ta_zzzz_xxxxyy_1 = pbuffer.data(idx_npot_1_gi + 395);

    auto ta_zzzz_xxxxyz_1 = pbuffer.data(idx_npot_1_gi + 396);

    auto ta_zzzz_xxxxzz_1 = pbuffer.data(idx_npot_1_gi + 397);

    auto ta_zzzz_xxxyyy_1 = pbuffer.data(idx_npot_1_gi + 398);

    auto ta_zzzz_xxxyyz_1 = pbuffer.data(idx_npot_1_gi + 399);

    auto ta_zzzz_xxxyzz_1 = pbuffer.data(idx_npot_1_gi + 400);

    auto ta_zzzz_xxxzzz_1 = pbuffer.data(idx_npot_1_gi + 401);

    auto ta_zzzz_xxyyyy_1 = pbuffer.data(idx_npot_1_gi + 402);

    auto ta_zzzz_xxyyyz_1 = pbuffer.data(idx_npot_1_gi + 403);

    auto ta_zzzz_xxyyzz_1 = pbuffer.data(idx_npot_1_gi + 404);

    auto ta_zzzz_xxyzzz_1 = pbuffer.data(idx_npot_1_gi + 405);

    auto ta_zzzz_xxzzzz_1 = pbuffer.data(idx_npot_1_gi + 406);

    auto ta_zzzz_xyyyyy_1 = pbuffer.data(idx_npot_1_gi + 407);

    auto ta_zzzz_xyyyyz_1 = pbuffer.data(idx_npot_1_gi + 408);

    auto ta_zzzz_xyyyzz_1 = pbuffer.data(idx_npot_1_gi + 409);

    auto ta_zzzz_xyyzzz_1 = pbuffer.data(idx_npot_1_gi + 410);

    auto ta_zzzz_xyzzzz_1 = pbuffer.data(idx_npot_1_gi + 411);

    auto ta_zzzz_xzzzzz_1 = pbuffer.data(idx_npot_1_gi + 412);

    auto ta_zzzz_yyyyyy_1 = pbuffer.data(idx_npot_1_gi + 413);

    auto ta_zzzz_yyyyyz_1 = pbuffer.data(idx_npot_1_gi + 414);

    auto ta_zzzz_yyyyzz_1 = pbuffer.data(idx_npot_1_gi + 415);

    auto ta_zzzz_yyyzzz_1 = pbuffer.data(idx_npot_1_gi + 416);

    auto ta_zzzz_yyzzzz_1 = pbuffer.data(idx_npot_1_gi + 417);

    auto ta_zzzz_yzzzzz_1 = pbuffer.data(idx_npot_1_gi + 418);

    auto ta_zzzz_zzzzzz_1 = pbuffer.data(idx_npot_1_gi + 419);

    // Set up 0-28 components of targeted buffer : HI

    auto ta_xxxxx_xxxxxx_0 = pbuffer.data(idx_npot_0_hi);

    auto ta_xxxxx_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 1);

    auto ta_xxxxx_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 2);

    auto ta_xxxxx_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 3);

    auto ta_xxxxx_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 4);

    auto ta_xxxxx_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 5);

    auto ta_xxxxx_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 6);

    auto ta_xxxxx_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 7);

    auto ta_xxxxx_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 8);

    auto ta_xxxxx_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 9);

    auto ta_xxxxx_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 10);

    auto ta_xxxxx_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 11);

    auto ta_xxxxx_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 12);

    auto ta_xxxxx_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 13);

    auto ta_xxxxx_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 14);

    auto ta_xxxxx_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 15);

    auto ta_xxxxx_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 16);

    auto ta_xxxxx_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 17);

    auto ta_xxxxx_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 18);

    auto ta_xxxxx_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 19);

    auto ta_xxxxx_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 20);

    auto ta_xxxxx_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 21);

    auto ta_xxxxx_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 22);

    auto ta_xxxxx_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 23);

    auto ta_xxxxx_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 24);

    auto ta_xxxxx_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 25);

    auto ta_xxxxx_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 26);

    auto ta_xxxxx_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 27);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xxx_xxxxxx_0,   \
                             ta_xxx_xxxxxx_1,   \
                             ta_xxx_xxxxxy_0,   \
                             ta_xxx_xxxxxy_1,   \
                             ta_xxx_xxxxxz_0,   \
                             ta_xxx_xxxxxz_1,   \
                             ta_xxx_xxxxyy_0,   \
                             ta_xxx_xxxxyy_1,   \
                             ta_xxx_xxxxyz_0,   \
                             ta_xxx_xxxxyz_1,   \
                             ta_xxx_xxxxzz_0,   \
                             ta_xxx_xxxxzz_1,   \
                             ta_xxx_xxxyyy_0,   \
                             ta_xxx_xxxyyy_1,   \
                             ta_xxx_xxxyyz_0,   \
                             ta_xxx_xxxyyz_1,   \
                             ta_xxx_xxxyzz_0,   \
                             ta_xxx_xxxyzz_1,   \
                             ta_xxx_xxxzzz_0,   \
                             ta_xxx_xxxzzz_1,   \
                             ta_xxx_xxyyyy_0,   \
                             ta_xxx_xxyyyy_1,   \
                             ta_xxx_xxyyyz_0,   \
                             ta_xxx_xxyyyz_1,   \
                             ta_xxx_xxyyzz_0,   \
                             ta_xxx_xxyyzz_1,   \
                             ta_xxx_xxyzzz_0,   \
                             ta_xxx_xxyzzz_1,   \
                             ta_xxx_xxzzzz_0,   \
                             ta_xxx_xxzzzz_1,   \
                             ta_xxx_xyyyyy_0,   \
                             ta_xxx_xyyyyy_1,   \
                             ta_xxx_xyyyyz_0,   \
                             ta_xxx_xyyyyz_1,   \
                             ta_xxx_xyyyzz_0,   \
                             ta_xxx_xyyyzz_1,   \
                             ta_xxx_xyyzzz_0,   \
                             ta_xxx_xyyzzz_1,   \
                             ta_xxx_xyzzzz_0,   \
                             ta_xxx_xyzzzz_1,   \
                             ta_xxx_xzzzzz_0,   \
                             ta_xxx_xzzzzz_1,   \
                             ta_xxx_yyyyyy_0,   \
                             ta_xxx_yyyyyy_1,   \
                             ta_xxx_yyyyyz_0,   \
                             ta_xxx_yyyyyz_1,   \
                             ta_xxx_yyyyzz_0,   \
                             ta_xxx_yyyyzz_1,   \
                             ta_xxx_yyyzzz_0,   \
                             ta_xxx_yyyzzz_1,   \
                             ta_xxx_yyzzzz_0,   \
                             ta_xxx_yyzzzz_1,   \
                             ta_xxx_yzzzzz_0,   \
                             ta_xxx_yzzzzz_1,   \
                             ta_xxx_zzzzzz_0,   \
                             ta_xxx_zzzzzz_1,   \
                             ta_xxxx_xxxxx_0,   \
                             ta_xxxx_xxxxx_1,   \
                             ta_xxxx_xxxxxx_0,  \
                             ta_xxxx_xxxxxx_1,  \
                             ta_xxxx_xxxxxy_0,  \
                             ta_xxxx_xxxxxy_1,  \
                             ta_xxxx_xxxxxz_0,  \
                             ta_xxxx_xxxxxz_1,  \
                             ta_xxxx_xxxxy_0,   \
                             ta_xxxx_xxxxy_1,   \
                             ta_xxxx_xxxxyy_0,  \
                             ta_xxxx_xxxxyy_1,  \
                             ta_xxxx_xxxxyz_0,  \
                             ta_xxxx_xxxxyz_1,  \
                             ta_xxxx_xxxxz_0,   \
                             ta_xxxx_xxxxz_1,   \
                             ta_xxxx_xxxxzz_0,  \
                             ta_xxxx_xxxxzz_1,  \
                             ta_xxxx_xxxyy_0,   \
                             ta_xxxx_xxxyy_1,   \
                             ta_xxxx_xxxyyy_0,  \
                             ta_xxxx_xxxyyy_1,  \
                             ta_xxxx_xxxyyz_0,  \
                             ta_xxxx_xxxyyz_1,  \
                             ta_xxxx_xxxyz_0,   \
                             ta_xxxx_xxxyz_1,   \
                             ta_xxxx_xxxyzz_0,  \
                             ta_xxxx_xxxyzz_1,  \
                             ta_xxxx_xxxzz_0,   \
                             ta_xxxx_xxxzz_1,   \
                             ta_xxxx_xxxzzz_0,  \
                             ta_xxxx_xxxzzz_1,  \
                             ta_xxxx_xxyyy_0,   \
                             ta_xxxx_xxyyy_1,   \
                             ta_xxxx_xxyyyy_0,  \
                             ta_xxxx_xxyyyy_1,  \
                             ta_xxxx_xxyyyz_0,  \
                             ta_xxxx_xxyyyz_1,  \
                             ta_xxxx_xxyyz_0,   \
                             ta_xxxx_xxyyz_1,   \
                             ta_xxxx_xxyyzz_0,  \
                             ta_xxxx_xxyyzz_1,  \
                             ta_xxxx_xxyzz_0,   \
                             ta_xxxx_xxyzz_1,   \
                             ta_xxxx_xxyzzz_0,  \
                             ta_xxxx_xxyzzz_1,  \
                             ta_xxxx_xxzzz_0,   \
                             ta_xxxx_xxzzz_1,   \
                             ta_xxxx_xxzzzz_0,  \
                             ta_xxxx_xxzzzz_1,  \
                             ta_xxxx_xyyyy_0,   \
                             ta_xxxx_xyyyy_1,   \
                             ta_xxxx_xyyyyy_0,  \
                             ta_xxxx_xyyyyy_1,  \
                             ta_xxxx_xyyyyz_0,  \
                             ta_xxxx_xyyyyz_1,  \
                             ta_xxxx_xyyyz_0,   \
                             ta_xxxx_xyyyz_1,   \
                             ta_xxxx_xyyyzz_0,  \
                             ta_xxxx_xyyyzz_1,  \
                             ta_xxxx_xyyzz_0,   \
                             ta_xxxx_xyyzz_1,   \
                             ta_xxxx_xyyzzz_0,  \
                             ta_xxxx_xyyzzz_1,  \
                             ta_xxxx_xyzzz_0,   \
                             ta_xxxx_xyzzz_1,   \
                             ta_xxxx_xyzzzz_0,  \
                             ta_xxxx_xyzzzz_1,  \
                             ta_xxxx_xzzzz_0,   \
                             ta_xxxx_xzzzz_1,   \
                             ta_xxxx_xzzzzz_0,  \
                             ta_xxxx_xzzzzz_1,  \
                             ta_xxxx_yyyyy_0,   \
                             ta_xxxx_yyyyy_1,   \
                             ta_xxxx_yyyyyy_0,  \
                             ta_xxxx_yyyyyy_1,  \
                             ta_xxxx_yyyyyz_0,  \
                             ta_xxxx_yyyyyz_1,  \
                             ta_xxxx_yyyyz_0,   \
                             ta_xxxx_yyyyz_1,   \
                             ta_xxxx_yyyyzz_0,  \
                             ta_xxxx_yyyyzz_1,  \
                             ta_xxxx_yyyzz_0,   \
                             ta_xxxx_yyyzz_1,   \
                             ta_xxxx_yyyzzz_0,  \
                             ta_xxxx_yyyzzz_1,  \
                             ta_xxxx_yyzzz_0,   \
                             ta_xxxx_yyzzz_1,   \
                             ta_xxxx_yyzzzz_0,  \
                             ta_xxxx_yyzzzz_1,  \
                             ta_xxxx_yzzzz_0,   \
                             ta_xxxx_yzzzz_1,   \
                             ta_xxxx_yzzzzz_0,  \
                             ta_xxxx_yzzzzz_1,  \
                             ta_xxxx_zzzzz_0,   \
                             ta_xxxx_zzzzz_1,   \
                             ta_xxxx_zzzzzz_0,  \
                             ta_xxxx_zzzzzz_1,  \
                             ta_xxxxx_xxxxxx_0, \
                             ta_xxxxx_xxxxxy_0, \
                             ta_xxxxx_xxxxxz_0, \
                             ta_xxxxx_xxxxyy_0, \
                             ta_xxxxx_xxxxyz_0, \
                             ta_xxxxx_xxxxzz_0, \
                             ta_xxxxx_xxxyyy_0, \
                             ta_xxxxx_xxxyyz_0, \
                             ta_xxxxx_xxxyzz_0, \
                             ta_xxxxx_xxxzzz_0, \
                             ta_xxxxx_xxyyyy_0, \
                             ta_xxxxx_xxyyyz_0, \
                             ta_xxxxx_xxyyzz_0, \
                             ta_xxxxx_xxyzzz_0, \
                             ta_xxxxx_xxzzzz_0, \
                             ta_xxxxx_xyyyyy_0, \
                             ta_xxxxx_xyyyyz_0, \
                             ta_xxxxx_xyyyzz_0, \
                             ta_xxxxx_xyyzzz_0, \
                             ta_xxxxx_xyzzzz_0, \
                             ta_xxxxx_xzzzzz_0, \
                             ta_xxxxx_yyyyyy_0, \
                             ta_xxxxx_yyyyyz_0, \
                             ta_xxxxx_yyyyzz_0, \
                             ta_xxxxx_yyyzzz_0, \
                             ta_xxxxx_yyzzzz_0, \
                             ta_xxxxx_yzzzzz_0, \
                             ta_xxxxx_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxx_xxxxxx_0[i] = 4.0 * ta_xxx_xxxxxx_0[i] * fe_0 - 4.0 * ta_xxx_xxxxxx_1[i] * fe_0 + 6.0 * ta_xxxx_xxxxx_0[i] * fe_0 -
                               6.0 * ta_xxxx_xxxxx_1[i] * fe_0 + ta_xxxx_xxxxxx_0[i] * pa_x[i] - ta_xxxx_xxxxxx_1[i] * pc_x[i];

        ta_xxxxx_xxxxxy_0[i] = 4.0 * ta_xxx_xxxxxy_0[i] * fe_0 - 4.0 * ta_xxx_xxxxxy_1[i] * fe_0 + 5.0 * ta_xxxx_xxxxy_0[i] * fe_0 -
                               5.0 * ta_xxxx_xxxxy_1[i] * fe_0 + ta_xxxx_xxxxxy_0[i] * pa_x[i] - ta_xxxx_xxxxxy_1[i] * pc_x[i];

        ta_xxxxx_xxxxxz_0[i] = 4.0 * ta_xxx_xxxxxz_0[i] * fe_0 - 4.0 * ta_xxx_xxxxxz_1[i] * fe_0 + 5.0 * ta_xxxx_xxxxz_0[i] * fe_0 -
                               5.0 * ta_xxxx_xxxxz_1[i] * fe_0 + ta_xxxx_xxxxxz_0[i] * pa_x[i] - ta_xxxx_xxxxxz_1[i] * pc_x[i];

        ta_xxxxx_xxxxyy_0[i] = 4.0 * ta_xxx_xxxxyy_0[i] * fe_0 - 4.0 * ta_xxx_xxxxyy_1[i] * fe_0 + 4.0 * ta_xxxx_xxxyy_0[i] * fe_0 -
                               4.0 * ta_xxxx_xxxyy_1[i] * fe_0 + ta_xxxx_xxxxyy_0[i] * pa_x[i] - ta_xxxx_xxxxyy_1[i] * pc_x[i];

        ta_xxxxx_xxxxyz_0[i] = 4.0 * ta_xxx_xxxxyz_0[i] * fe_0 - 4.0 * ta_xxx_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxxx_xxxyz_0[i] * fe_0 -
                               4.0 * ta_xxxx_xxxyz_1[i] * fe_0 + ta_xxxx_xxxxyz_0[i] * pa_x[i] - ta_xxxx_xxxxyz_1[i] * pc_x[i];

        ta_xxxxx_xxxxzz_0[i] = 4.0 * ta_xxx_xxxxzz_0[i] * fe_0 - 4.0 * ta_xxx_xxxxzz_1[i] * fe_0 + 4.0 * ta_xxxx_xxxzz_0[i] * fe_0 -
                               4.0 * ta_xxxx_xxxzz_1[i] * fe_0 + ta_xxxx_xxxxzz_0[i] * pa_x[i] - ta_xxxx_xxxxzz_1[i] * pc_x[i];

        ta_xxxxx_xxxyyy_0[i] = 4.0 * ta_xxx_xxxyyy_0[i] * fe_0 - 4.0 * ta_xxx_xxxyyy_1[i] * fe_0 + 3.0 * ta_xxxx_xxyyy_0[i] * fe_0 -
                               3.0 * ta_xxxx_xxyyy_1[i] * fe_0 + ta_xxxx_xxxyyy_0[i] * pa_x[i] - ta_xxxx_xxxyyy_1[i] * pc_x[i];

        ta_xxxxx_xxxyyz_0[i] = 4.0 * ta_xxx_xxxyyz_0[i] * fe_0 - 4.0 * ta_xxx_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxxx_xxyyz_0[i] * fe_0 -
                               3.0 * ta_xxxx_xxyyz_1[i] * fe_0 + ta_xxxx_xxxyyz_0[i] * pa_x[i] - ta_xxxx_xxxyyz_1[i] * pc_x[i];

        ta_xxxxx_xxxyzz_0[i] = 4.0 * ta_xxx_xxxyzz_0[i] * fe_0 - 4.0 * ta_xxx_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxxx_xxyzz_0[i] * fe_0 -
                               3.0 * ta_xxxx_xxyzz_1[i] * fe_0 + ta_xxxx_xxxyzz_0[i] * pa_x[i] - ta_xxxx_xxxyzz_1[i] * pc_x[i];

        ta_xxxxx_xxxzzz_0[i] = 4.0 * ta_xxx_xxxzzz_0[i] * fe_0 - 4.0 * ta_xxx_xxxzzz_1[i] * fe_0 + 3.0 * ta_xxxx_xxzzz_0[i] * fe_0 -
                               3.0 * ta_xxxx_xxzzz_1[i] * fe_0 + ta_xxxx_xxxzzz_0[i] * pa_x[i] - ta_xxxx_xxxzzz_1[i] * pc_x[i];

        ta_xxxxx_xxyyyy_0[i] = 4.0 * ta_xxx_xxyyyy_0[i] * fe_0 - 4.0 * ta_xxx_xxyyyy_1[i] * fe_0 + 2.0 * ta_xxxx_xyyyy_0[i] * fe_0 -
                               2.0 * ta_xxxx_xyyyy_1[i] * fe_0 + ta_xxxx_xxyyyy_0[i] * pa_x[i] - ta_xxxx_xxyyyy_1[i] * pc_x[i];

        ta_xxxxx_xxyyyz_0[i] = 4.0 * ta_xxx_xxyyyz_0[i] * fe_0 - 4.0 * ta_xxx_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxxx_xyyyz_0[i] * fe_0 -
                               2.0 * ta_xxxx_xyyyz_1[i] * fe_0 + ta_xxxx_xxyyyz_0[i] * pa_x[i] - ta_xxxx_xxyyyz_1[i] * pc_x[i];

        ta_xxxxx_xxyyzz_0[i] = 4.0 * ta_xxx_xxyyzz_0[i] * fe_0 - 4.0 * ta_xxx_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxxx_xyyzz_0[i] * fe_0 -
                               2.0 * ta_xxxx_xyyzz_1[i] * fe_0 + ta_xxxx_xxyyzz_0[i] * pa_x[i] - ta_xxxx_xxyyzz_1[i] * pc_x[i];

        ta_xxxxx_xxyzzz_0[i] = 4.0 * ta_xxx_xxyzzz_0[i] * fe_0 - 4.0 * ta_xxx_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxxx_xyzzz_0[i] * fe_0 -
                               2.0 * ta_xxxx_xyzzz_1[i] * fe_0 + ta_xxxx_xxyzzz_0[i] * pa_x[i] - ta_xxxx_xxyzzz_1[i] * pc_x[i];

        ta_xxxxx_xxzzzz_0[i] = 4.0 * ta_xxx_xxzzzz_0[i] * fe_0 - 4.0 * ta_xxx_xxzzzz_1[i] * fe_0 + 2.0 * ta_xxxx_xzzzz_0[i] * fe_0 -
                               2.0 * ta_xxxx_xzzzz_1[i] * fe_0 + ta_xxxx_xxzzzz_0[i] * pa_x[i] - ta_xxxx_xxzzzz_1[i] * pc_x[i];

        ta_xxxxx_xyyyyy_0[i] = 4.0 * ta_xxx_xyyyyy_0[i] * fe_0 - 4.0 * ta_xxx_xyyyyy_1[i] * fe_0 + ta_xxxx_yyyyy_0[i] * fe_0 -
                               ta_xxxx_yyyyy_1[i] * fe_0 + ta_xxxx_xyyyyy_0[i] * pa_x[i] - ta_xxxx_xyyyyy_1[i] * pc_x[i];

        ta_xxxxx_xyyyyz_0[i] = 4.0 * ta_xxx_xyyyyz_0[i] * fe_0 - 4.0 * ta_xxx_xyyyyz_1[i] * fe_0 + ta_xxxx_yyyyz_0[i] * fe_0 -
                               ta_xxxx_yyyyz_1[i] * fe_0 + ta_xxxx_xyyyyz_0[i] * pa_x[i] - ta_xxxx_xyyyyz_1[i] * pc_x[i];

        ta_xxxxx_xyyyzz_0[i] = 4.0 * ta_xxx_xyyyzz_0[i] * fe_0 - 4.0 * ta_xxx_xyyyzz_1[i] * fe_0 + ta_xxxx_yyyzz_0[i] * fe_0 -
                               ta_xxxx_yyyzz_1[i] * fe_0 + ta_xxxx_xyyyzz_0[i] * pa_x[i] - ta_xxxx_xyyyzz_1[i] * pc_x[i];

        ta_xxxxx_xyyzzz_0[i] = 4.0 * ta_xxx_xyyzzz_0[i] * fe_0 - 4.0 * ta_xxx_xyyzzz_1[i] * fe_0 + ta_xxxx_yyzzz_0[i] * fe_0 -
                               ta_xxxx_yyzzz_1[i] * fe_0 + ta_xxxx_xyyzzz_0[i] * pa_x[i] - ta_xxxx_xyyzzz_1[i] * pc_x[i];

        ta_xxxxx_xyzzzz_0[i] = 4.0 * ta_xxx_xyzzzz_0[i] * fe_0 - 4.0 * ta_xxx_xyzzzz_1[i] * fe_0 + ta_xxxx_yzzzz_0[i] * fe_0 -
                               ta_xxxx_yzzzz_1[i] * fe_0 + ta_xxxx_xyzzzz_0[i] * pa_x[i] - ta_xxxx_xyzzzz_1[i] * pc_x[i];

        ta_xxxxx_xzzzzz_0[i] = 4.0 * ta_xxx_xzzzzz_0[i] * fe_0 - 4.0 * ta_xxx_xzzzzz_1[i] * fe_0 + ta_xxxx_zzzzz_0[i] * fe_0 -
                               ta_xxxx_zzzzz_1[i] * fe_0 + ta_xxxx_xzzzzz_0[i] * pa_x[i] - ta_xxxx_xzzzzz_1[i] * pc_x[i];

        ta_xxxxx_yyyyyy_0[i] =
            4.0 * ta_xxx_yyyyyy_0[i] * fe_0 - 4.0 * ta_xxx_yyyyyy_1[i] * fe_0 + ta_xxxx_yyyyyy_0[i] * pa_x[i] - ta_xxxx_yyyyyy_1[i] * pc_x[i];

        ta_xxxxx_yyyyyz_0[i] =
            4.0 * ta_xxx_yyyyyz_0[i] * fe_0 - 4.0 * ta_xxx_yyyyyz_1[i] * fe_0 + ta_xxxx_yyyyyz_0[i] * pa_x[i] - ta_xxxx_yyyyyz_1[i] * pc_x[i];

        ta_xxxxx_yyyyzz_0[i] =
            4.0 * ta_xxx_yyyyzz_0[i] * fe_0 - 4.0 * ta_xxx_yyyyzz_1[i] * fe_0 + ta_xxxx_yyyyzz_0[i] * pa_x[i] - ta_xxxx_yyyyzz_1[i] * pc_x[i];

        ta_xxxxx_yyyzzz_0[i] =
            4.0 * ta_xxx_yyyzzz_0[i] * fe_0 - 4.0 * ta_xxx_yyyzzz_1[i] * fe_0 + ta_xxxx_yyyzzz_0[i] * pa_x[i] - ta_xxxx_yyyzzz_1[i] * pc_x[i];

        ta_xxxxx_yyzzzz_0[i] =
            4.0 * ta_xxx_yyzzzz_0[i] * fe_0 - 4.0 * ta_xxx_yyzzzz_1[i] * fe_0 + ta_xxxx_yyzzzz_0[i] * pa_x[i] - ta_xxxx_yyzzzz_1[i] * pc_x[i];

        ta_xxxxx_yzzzzz_0[i] =
            4.0 * ta_xxx_yzzzzz_0[i] * fe_0 - 4.0 * ta_xxx_yzzzzz_1[i] * fe_0 + ta_xxxx_yzzzzz_0[i] * pa_x[i] - ta_xxxx_yzzzzz_1[i] * pc_x[i];

        ta_xxxxx_zzzzzz_0[i] =
            4.0 * ta_xxx_zzzzzz_0[i] * fe_0 - 4.0 * ta_xxx_zzzzzz_1[i] * fe_0 + ta_xxxx_zzzzzz_0[i] * pa_x[i] - ta_xxxx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : HI

    auto ta_xxxxy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 28);

    auto ta_xxxxy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 29);

    auto ta_xxxxy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 30);

    auto ta_xxxxy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 31);

    auto ta_xxxxy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 32);

    auto ta_xxxxy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 33);

    auto ta_xxxxy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 34);

    auto ta_xxxxy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 35);

    auto ta_xxxxy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 36);

    auto ta_xxxxy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 37);

    auto ta_xxxxy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 38);

    auto ta_xxxxy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 39);

    auto ta_xxxxy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 40);

    auto ta_xxxxy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 41);

    auto ta_xxxxy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 42);

    auto ta_xxxxy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 43);

    auto ta_xxxxy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 44);

    auto ta_xxxxy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 45);

    auto ta_xxxxy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 46);

    auto ta_xxxxy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 47);

    auto ta_xxxxy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 48);

    auto ta_xxxxy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 49);

    auto ta_xxxxy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 50);

    auto ta_xxxxy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 51);

    auto ta_xxxxy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 52);

    auto ta_xxxxy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 53);

    auto ta_xxxxy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 54);

    auto ta_xxxxy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 55);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxxx_xxxxx_0,   \
                             ta_xxxx_xxxxx_1,   \
                             ta_xxxx_xxxxxx_0,  \
                             ta_xxxx_xxxxxx_1,  \
                             ta_xxxx_xxxxxy_0,  \
                             ta_xxxx_xxxxxy_1,  \
                             ta_xxxx_xxxxxz_0,  \
                             ta_xxxx_xxxxxz_1,  \
                             ta_xxxx_xxxxy_0,   \
                             ta_xxxx_xxxxy_1,   \
                             ta_xxxx_xxxxyy_0,  \
                             ta_xxxx_xxxxyy_1,  \
                             ta_xxxx_xxxxyz_0,  \
                             ta_xxxx_xxxxyz_1,  \
                             ta_xxxx_xxxxz_0,   \
                             ta_xxxx_xxxxz_1,   \
                             ta_xxxx_xxxxzz_0,  \
                             ta_xxxx_xxxxzz_1,  \
                             ta_xxxx_xxxyy_0,   \
                             ta_xxxx_xxxyy_1,   \
                             ta_xxxx_xxxyyy_0,  \
                             ta_xxxx_xxxyyy_1,  \
                             ta_xxxx_xxxyyz_0,  \
                             ta_xxxx_xxxyyz_1,  \
                             ta_xxxx_xxxyz_0,   \
                             ta_xxxx_xxxyz_1,   \
                             ta_xxxx_xxxyzz_0,  \
                             ta_xxxx_xxxyzz_1,  \
                             ta_xxxx_xxxzz_0,   \
                             ta_xxxx_xxxzz_1,   \
                             ta_xxxx_xxxzzz_0,  \
                             ta_xxxx_xxxzzz_1,  \
                             ta_xxxx_xxyyy_0,   \
                             ta_xxxx_xxyyy_1,   \
                             ta_xxxx_xxyyyy_0,  \
                             ta_xxxx_xxyyyy_1,  \
                             ta_xxxx_xxyyyz_0,  \
                             ta_xxxx_xxyyyz_1,  \
                             ta_xxxx_xxyyz_0,   \
                             ta_xxxx_xxyyz_1,   \
                             ta_xxxx_xxyyzz_0,  \
                             ta_xxxx_xxyyzz_1,  \
                             ta_xxxx_xxyzz_0,   \
                             ta_xxxx_xxyzz_1,   \
                             ta_xxxx_xxyzzz_0,  \
                             ta_xxxx_xxyzzz_1,  \
                             ta_xxxx_xxzzz_0,   \
                             ta_xxxx_xxzzz_1,   \
                             ta_xxxx_xxzzzz_0,  \
                             ta_xxxx_xxzzzz_1,  \
                             ta_xxxx_xyyyy_0,   \
                             ta_xxxx_xyyyy_1,   \
                             ta_xxxx_xyyyyy_0,  \
                             ta_xxxx_xyyyyy_1,  \
                             ta_xxxx_xyyyyz_0,  \
                             ta_xxxx_xyyyyz_1,  \
                             ta_xxxx_xyyyz_0,   \
                             ta_xxxx_xyyyz_1,   \
                             ta_xxxx_xyyyzz_0,  \
                             ta_xxxx_xyyyzz_1,  \
                             ta_xxxx_xyyzz_0,   \
                             ta_xxxx_xyyzz_1,   \
                             ta_xxxx_xyyzzz_0,  \
                             ta_xxxx_xyyzzz_1,  \
                             ta_xxxx_xyzzz_0,   \
                             ta_xxxx_xyzzz_1,   \
                             ta_xxxx_xyzzzz_0,  \
                             ta_xxxx_xyzzzz_1,  \
                             ta_xxxx_xzzzz_0,   \
                             ta_xxxx_xzzzz_1,   \
                             ta_xxxx_xzzzzz_0,  \
                             ta_xxxx_xzzzzz_1,  \
                             ta_xxxx_zzzzzz_0,  \
                             ta_xxxx_zzzzzz_1,  \
                             ta_xxxxy_xxxxxx_0, \
                             ta_xxxxy_xxxxxy_0, \
                             ta_xxxxy_xxxxxz_0, \
                             ta_xxxxy_xxxxyy_0, \
                             ta_xxxxy_xxxxyz_0, \
                             ta_xxxxy_xxxxzz_0, \
                             ta_xxxxy_xxxyyy_0, \
                             ta_xxxxy_xxxyyz_0, \
                             ta_xxxxy_xxxyzz_0, \
                             ta_xxxxy_xxxzzz_0, \
                             ta_xxxxy_xxyyyy_0, \
                             ta_xxxxy_xxyyyz_0, \
                             ta_xxxxy_xxyyzz_0, \
                             ta_xxxxy_xxyzzz_0, \
                             ta_xxxxy_xxzzzz_0, \
                             ta_xxxxy_xyyyyy_0, \
                             ta_xxxxy_xyyyyz_0, \
                             ta_xxxxy_xyyyzz_0, \
                             ta_xxxxy_xyyzzz_0, \
                             ta_xxxxy_xyzzzz_0, \
                             ta_xxxxy_xzzzzz_0, \
                             ta_xxxxy_yyyyyy_0, \
                             ta_xxxxy_yyyyyz_0, \
                             ta_xxxxy_yyyyzz_0, \
                             ta_xxxxy_yyyzzz_0, \
                             ta_xxxxy_yyzzzz_0, \
                             ta_xxxxy_yzzzzz_0, \
                             ta_xxxxy_zzzzzz_0, \
                             ta_xxxy_yyyyyy_0,  \
                             ta_xxxy_yyyyyy_1,  \
                             ta_xxxy_yyyyyz_0,  \
                             ta_xxxy_yyyyyz_1,  \
                             ta_xxxy_yyyyzz_0,  \
                             ta_xxxy_yyyyzz_1,  \
                             ta_xxxy_yyyzzz_0,  \
                             ta_xxxy_yyyzzz_1,  \
                             ta_xxxy_yyzzzz_0,  \
                             ta_xxxy_yyzzzz_1,  \
                             ta_xxxy_yzzzzz_0,  \
                             ta_xxxy_yzzzzz_1,  \
                             ta_xxy_yyyyyy_0,   \
                             ta_xxy_yyyyyy_1,   \
                             ta_xxy_yyyyyz_0,   \
                             ta_xxy_yyyyyz_1,   \
                             ta_xxy_yyyyzz_0,   \
                             ta_xxy_yyyyzz_1,   \
                             ta_xxy_yyyzzz_0,   \
                             ta_xxy_yyyzzz_1,   \
                             ta_xxy_yyzzzz_0,   \
                             ta_xxy_yyzzzz_1,   \
                             ta_xxy_yzzzzz_0,   \
                             ta_xxy_yzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxy_xxxxxx_0[i] = ta_xxxx_xxxxxx_0[i] * pa_y[i] - ta_xxxx_xxxxxx_1[i] * pc_y[i];

        ta_xxxxy_xxxxxy_0[i] = ta_xxxx_xxxxx_0[i] * fe_0 - ta_xxxx_xxxxx_1[i] * fe_0 + ta_xxxx_xxxxxy_0[i] * pa_y[i] - ta_xxxx_xxxxxy_1[i] * pc_y[i];

        ta_xxxxy_xxxxxz_0[i] = ta_xxxx_xxxxxz_0[i] * pa_y[i] - ta_xxxx_xxxxxz_1[i] * pc_y[i];

        ta_xxxxy_xxxxyy_0[i] =
            2.0 * ta_xxxx_xxxxy_0[i] * fe_0 - 2.0 * ta_xxxx_xxxxy_1[i] * fe_0 + ta_xxxx_xxxxyy_0[i] * pa_y[i] - ta_xxxx_xxxxyy_1[i] * pc_y[i];

        ta_xxxxy_xxxxyz_0[i] = ta_xxxx_xxxxz_0[i] * fe_0 - ta_xxxx_xxxxz_1[i] * fe_0 + ta_xxxx_xxxxyz_0[i] * pa_y[i] - ta_xxxx_xxxxyz_1[i] * pc_y[i];

        ta_xxxxy_xxxxzz_0[i] = ta_xxxx_xxxxzz_0[i] * pa_y[i] - ta_xxxx_xxxxzz_1[i] * pc_y[i];

        ta_xxxxy_xxxyyy_0[i] =
            3.0 * ta_xxxx_xxxyy_0[i] * fe_0 - 3.0 * ta_xxxx_xxxyy_1[i] * fe_0 + ta_xxxx_xxxyyy_0[i] * pa_y[i] - ta_xxxx_xxxyyy_1[i] * pc_y[i];

        ta_xxxxy_xxxyyz_0[i] =
            2.0 * ta_xxxx_xxxyz_0[i] * fe_0 - 2.0 * ta_xxxx_xxxyz_1[i] * fe_0 + ta_xxxx_xxxyyz_0[i] * pa_y[i] - ta_xxxx_xxxyyz_1[i] * pc_y[i];

        ta_xxxxy_xxxyzz_0[i] = ta_xxxx_xxxzz_0[i] * fe_0 - ta_xxxx_xxxzz_1[i] * fe_0 + ta_xxxx_xxxyzz_0[i] * pa_y[i] - ta_xxxx_xxxyzz_1[i] * pc_y[i];

        ta_xxxxy_xxxzzz_0[i] = ta_xxxx_xxxzzz_0[i] * pa_y[i] - ta_xxxx_xxxzzz_1[i] * pc_y[i];

        ta_xxxxy_xxyyyy_0[i] =
            4.0 * ta_xxxx_xxyyy_0[i] * fe_0 - 4.0 * ta_xxxx_xxyyy_1[i] * fe_0 + ta_xxxx_xxyyyy_0[i] * pa_y[i] - ta_xxxx_xxyyyy_1[i] * pc_y[i];

        ta_xxxxy_xxyyyz_0[i] =
            3.0 * ta_xxxx_xxyyz_0[i] * fe_0 - 3.0 * ta_xxxx_xxyyz_1[i] * fe_0 + ta_xxxx_xxyyyz_0[i] * pa_y[i] - ta_xxxx_xxyyyz_1[i] * pc_y[i];

        ta_xxxxy_xxyyzz_0[i] =
            2.0 * ta_xxxx_xxyzz_0[i] * fe_0 - 2.0 * ta_xxxx_xxyzz_1[i] * fe_0 + ta_xxxx_xxyyzz_0[i] * pa_y[i] - ta_xxxx_xxyyzz_1[i] * pc_y[i];

        ta_xxxxy_xxyzzz_0[i] = ta_xxxx_xxzzz_0[i] * fe_0 - ta_xxxx_xxzzz_1[i] * fe_0 + ta_xxxx_xxyzzz_0[i] * pa_y[i] - ta_xxxx_xxyzzz_1[i] * pc_y[i];

        ta_xxxxy_xxzzzz_0[i] = ta_xxxx_xxzzzz_0[i] * pa_y[i] - ta_xxxx_xxzzzz_1[i] * pc_y[i];

        ta_xxxxy_xyyyyy_0[i] =
            5.0 * ta_xxxx_xyyyy_0[i] * fe_0 - 5.0 * ta_xxxx_xyyyy_1[i] * fe_0 + ta_xxxx_xyyyyy_0[i] * pa_y[i] - ta_xxxx_xyyyyy_1[i] * pc_y[i];

        ta_xxxxy_xyyyyz_0[i] =
            4.0 * ta_xxxx_xyyyz_0[i] * fe_0 - 4.0 * ta_xxxx_xyyyz_1[i] * fe_0 + ta_xxxx_xyyyyz_0[i] * pa_y[i] - ta_xxxx_xyyyyz_1[i] * pc_y[i];

        ta_xxxxy_xyyyzz_0[i] =
            3.0 * ta_xxxx_xyyzz_0[i] * fe_0 - 3.0 * ta_xxxx_xyyzz_1[i] * fe_0 + ta_xxxx_xyyyzz_0[i] * pa_y[i] - ta_xxxx_xyyyzz_1[i] * pc_y[i];

        ta_xxxxy_xyyzzz_0[i] =
            2.0 * ta_xxxx_xyzzz_0[i] * fe_0 - 2.0 * ta_xxxx_xyzzz_1[i] * fe_0 + ta_xxxx_xyyzzz_0[i] * pa_y[i] - ta_xxxx_xyyzzz_1[i] * pc_y[i];

        ta_xxxxy_xyzzzz_0[i] = ta_xxxx_xzzzz_0[i] * fe_0 - ta_xxxx_xzzzz_1[i] * fe_0 + ta_xxxx_xyzzzz_0[i] * pa_y[i] - ta_xxxx_xyzzzz_1[i] * pc_y[i];

        ta_xxxxy_xzzzzz_0[i] = ta_xxxx_xzzzzz_0[i] * pa_y[i] - ta_xxxx_xzzzzz_1[i] * pc_y[i];

        ta_xxxxy_yyyyyy_0[i] =
            3.0 * ta_xxy_yyyyyy_0[i] * fe_0 - 3.0 * ta_xxy_yyyyyy_1[i] * fe_0 + ta_xxxy_yyyyyy_0[i] * pa_x[i] - ta_xxxy_yyyyyy_1[i] * pc_x[i];

        ta_xxxxy_yyyyyz_0[i] =
            3.0 * ta_xxy_yyyyyz_0[i] * fe_0 - 3.0 * ta_xxy_yyyyyz_1[i] * fe_0 + ta_xxxy_yyyyyz_0[i] * pa_x[i] - ta_xxxy_yyyyyz_1[i] * pc_x[i];

        ta_xxxxy_yyyyzz_0[i] =
            3.0 * ta_xxy_yyyyzz_0[i] * fe_0 - 3.0 * ta_xxy_yyyyzz_1[i] * fe_0 + ta_xxxy_yyyyzz_0[i] * pa_x[i] - ta_xxxy_yyyyzz_1[i] * pc_x[i];

        ta_xxxxy_yyyzzz_0[i] =
            3.0 * ta_xxy_yyyzzz_0[i] * fe_0 - 3.0 * ta_xxy_yyyzzz_1[i] * fe_0 + ta_xxxy_yyyzzz_0[i] * pa_x[i] - ta_xxxy_yyyzzz_1[i] * pc_x[i];

        ta_xxxxy_yyzzzz_0[i] =
            3.0 * ta_xxy_yyzzzz_0[i] * fe_0 - 3.0 * ta_xxy_yyzzzz_1[i] * fe_0 + ta_xxxy_yyzzzz_0[i] * pa_x[i] - ta_xxxy_yyzzzz_1[i] * pc_x[i];

        ta_xxxxy_yzzzzz_0[i] =
            3.0 * ta_xxy_yzzzzz_0[i] * fe_0 - 3.0 * ta_xxy_yzzzzz_1[i] * fe_0 + ta_xxxy_yzzzzz_0[i] * pa_x[i] - ta_xxxy_yzzzzz_1[i] * pc_x[i];

        ta_xxxxy_zzzzzz_0[i] = ta_xxxx_zzzzzz_0[i] * pa_y[i] - ta_xxxx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : HI

    auto ta_xxxxz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 56);

    auto ta_xxxxz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 57);

    auto ta_xxxxz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 58);

    auto ta_xxxxz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 59);

    auto ta_xxxxz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 60);

    auto ta_xxxxz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 61);

    auto ta_xxxxz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 62);

    auto ta_xxxxz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 63);

    auto ta_xxxxz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 64);

    auto ta_xxxxz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 65);

    auto ta_xxxxz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 66);

    auto ta_xxxxz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 67);

    auto ta_xxxxz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 68);

    auto ta_xxxxz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 69);

    auto ta_xxxxz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 70);

    auto ta_xxxxz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 71);

    auto ta_xxxxz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 72);

    auto ta_xxxxz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 73);

    auto ta_xxxxz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 74);

    auto ta_xxxxz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 75);

    auto ta_xxxxz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 76);

    auto ta_xxxxz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 77);

    auto ta_xxxxz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 78);

    auto ta_xxxxz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 79);

    auto ta_xxxxz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 80);

    auto ta_xxxxz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 81);

    auto ta_xxxxz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 82);

    auto ta_xxxxz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 83);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xxxx_xxxxx_0,   \
                             ta_xxxx_xxxxx_1,   \
                             ta_xxxx_xxxxxx_0,  \
                             ta_xxxx_xxxxxx_1,  \
                             ta_xxxx_xxxxxy_0,  \
                             ta_xxxx_xxxxxy_1,  \
                             ta_xxxx_xxxxxz_0,  \
                             ta_xxxx_xxxxxz_1,  \
                             ta_xxxx_xxxxy_0,   \
                             ta_xxxx_xxxxy_1,   \
                             ta_xxxx_xxxxyy_0,  \
                             ta_xxxx_xxxxyy_1,  \
                             ta_xxxx_xxxxyz_0,  \
                             ta_xxxx_xxxxyz_1,  \
                             ta_xxxx_xxxxz_0,   \
                             ta_xxxx_xxxxz_1,   \
                             ta_xxxx_xxxxzz_0,  \
                             ta_xxxx_xxxxzz_1,  \
                             ta_xxxx_xxxyy_0,   \
                             ta_xxxx_xxxyy_1,   \
                             ta_xxxx_xxxyyy_0,  \
                             ta_xxxx_xxxyyy_1,  \
                             ta_xxxx_xxxyyz_0,  \
                             ta_xxxx_xxxyyz_1,  \
                             ta_xxxx_xxxyz_0,   \
                             ta_xxxx_xxxyz_1,   \
                             ta_xxxx_xxxyzz_0,  \
                             ta_xxxx_xxxyzz_1,  \
                             ta_xxxx_xxxzz_0,   \
                             ta_xxxx_xxxzz_1,   \
                             ta_xxxx_xxxzzz_0,  \
                             ta_xxxx_xxxzzz_1,  \
                             ta_xxxx_xxyyy_0,   \
                             ta_xxxx_xxyyy_1,   \
                             ta_xxxx_xxyyyy_0,  \
                             ta_xxxx_xxyyyy_1,  \
                             ta_xxxx_xxyyyz_0,  \
                             ta_xxxx_xxyyyz_1,  \
                             ta_xxxx_xxyyz_0,   \
                             ta_xxxx_xxyyz_1,   \
                             ta_xxxx_xxyyzz_0,  \
                             ta_xxxx_xxyyzz_1,  \
                             ta_xxxx_xxyzz_0,   \
                             ta_xxxx_xxyzz_1,   \
                             ta_xxxx_xxyzzz_0,  \
                             ta_xxxx_xxyzzz_1,  \
                             ta_xxxx_xxzzz_0,   \
                             ta_xxxx_xxzzz_1,   \
                             ta_xxxx_xxzzzz_0,  \
                             ta_xxxx_xxzzzz_1,  \
                             ta_xxxx_xyyyy_0,   \
                             ta_xxxx_xyyyy_1,   \
                             ta_xxxx_xyyyyy_0,  \
                             ta_xxxx_xyyyyy_1,  \
                             ta_xxxx_xyyyyz_0,  \
                             ta_xxxx_xyyyyz_1,  \
                             ta_xxxx_xyyyz_0,   \
                             ta_xxxx_xyyyz_1,   \
                             ta_xxxx_xyyyzz_0,  \
                             ta_xxxx_xyyyzz_1,  \
                             ta_xxxx_xyyzz_0,   \
                             ta_xxxx_xyyzz_1,   \
                             ta_xxxx_xyyzzz_0,  \
                             ta_xxxx_xyyzzz_1,  \
                             ta_xxxx_xyzzz_0,   \
                             ta_xxxx_xyzzz_1,   \
                             ta_xxxx_xyzzzz_0,  \
                             ta_xxxx_xyzzzz_1,  \
                             ta_xxxx_xzzzz_0,   \
                             ta_xxxx_xzzzz_1,   \
                             ta_xxxx_xzzzzz_0,  \
                             ta_xxxx_xzzzzz_1,  \
                             ta_xxxx_yyyyyy_0,  \
                             ta_xxxx_yyyyyy_1,  \
                             ta_xxxxz_xxxxxx_0, \
                             ta_xxxxz_xxxxxy_0, \
                             ta_xxxxz_xxxxxz_0, \
                             ta_xxxxz_xxxxyy_0, \
                             ta_xxxxz_xxxxyz_0, \
                             ta_xxxxz_xxxxzz_0, \
                             ta_xxxxz_xxxyyy_0, \
                             ta_xxxxz_xxxyyz_0, \
                             ta_xxxxz_xxxyzz_0, \
                             ta_xxxxz_xxxzzz_0, \
                             ta_xxxxz_xxyyyy_0, \
                             ta_xxxxz_xxyyyz_0, \
                             ta_xxxxz_xxyyzz_0, \
                             ta_xxxxz_xxyzzz_0, \
                             ta_xxxxz_xxzzzz_0, \
                             ta_xxxxz_xyyyyy_0, \
                             ta_xxxxz_xyyyyz_0, \
                             ta_xxxxz_xyyyzz_0, \
                             ta_xxxxz_xyyzzz_0, \
                             ta_xxxxz_xyzzzz_0, \
                             ta_xxxxz_xzzzzz_0, \
                             ta_xxxxz_yyyyyy_0, \
                             ta_xxxxz_yyyyyz_0, \
                             ta_xxxxz_yyyyzz_0, \
                             ta_xxxxz_yyyzzz_0, \
                             ta_xxxxz_yyzzzz_0, \
                             ta_xxxxz_yzzzzz_0, \
                             ta_xxxxz_zzzzzz_0, \
                             ta_xxxz_yyyyyz_0,  \
                             ta_xxxz_yyyyyz_1,  \
                             ta_xxxz_yyyyzz_0,  \
                             ta_xxxz_yyyyzz_1,  \
                             ta_xxxz_yyyzzz_0,  \
                             ta_xxxz_yyyzzz_1,  \
                             ta_xxxz_yyzzzz_0,  \
                             ta_xxxz_yyzzzz_1,  \
                             ta_xxxz_yzzzzz_0,  \
                             ta_xxxz_yzzzzz_1,  \
                             ta_xxxz_zzzzzz_0,  \
                             ta_xxxz_zzzzzz_1,  \
                             ta_xxz_yyyyyz_0,   \
                             ta_xxz_yyyyyz_1,   \
                             ta_xxz_yyyyzz_0,   \
                             ta_xxz_yyyyzz_1,   \
                             ta_xxz_yyyzzz_0,   \
                             ta_xxz_yyyzzz_1,   \
                             ta_xxz_yyzzzz_0,   \
                             ta_xxz_yyzzzz_1,   \
                             ta_xxz_yzzzzz_0,   \
                             ta_xxz_yzzzzz_1,   \
                             ta_xxz_zzzzzz_0,   \
                             ta_xxz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxz_xxxxxx_0[i] = ta_xxxx_xxxxxx_0[i] * pa_z[i] - ta_xxxx_xxxxxx_1[i] * pc_z[i];

        ta_xxxxz_xxxxxy_0[i] = ta_xxxx_xxxxxy_0[i] * pa_z[i] - ta_xxxx_xxxxxy_1[i] * pc_z[i];

        ta_xxxxz_xxxxxz_0[i] = ta_xxxx_xxxxx_0[i] * fe_0 - ta_xxxx_xxxxx_1[i] * fe_0 + ta_xxxx_xxxxxz_0[i] * pa_z[i] - ta_xxxx_xxxxxz_1[i] * pc_z[i];

        ta_xxxxz_xxxxyy_0[i] = ta_xxxx_xxxxyy_0[i] * pa_z[i] - ta_xxxx_xxxxyy_1[i] * pc_z[i];

        ta_xxxxz_xxxxyz_0[i] = ta_xxxx_xxxxy_0[i] * fe_0 - ta_xxxx_xxxxy_1[i] * fe_0 + ta_xxxx_xxxxyz_0[i] * pa_z[i] - ta_xxxx_xxxxyz_1[i] * pc_z[i];

        ta_xxxxz_xxxxzz_0[i] =
            2.0 * ta_xxxx_xxxxz_0[i] * fe_0 - 2.0 * ta_xxxx_xxxxz_1[i] * fe_0 + ta_xxxx_xxxxzz_0[i] * pa_z[i] - ta_xxxx_xxxxzz_1[i] * pc_z[i];

        ta_xxxxz_xxxyyy_0[i] = ta_xxxx_xxxyyy_0[i] * pa_z[i] - ta_xxxx_xxxyyy_1[i] * pc_z[i];

        ta_xxxxz_xxxyyz_0[i] = ta_xxxx_xxxyy_0[i] * fe_0 - ta_xxxx_xxxyy_1[i] * fe_0 + ta_xxxx_xxxyyz_0[i] * pa_z[i] - ta_xxxx_xxxyyz_1[i] * pc_z[i];

        ta_xxxxz_xxxyzz_0[i] =
            2.0 * ta_xxxx_xxxyz_0[i] * fe_0 - 2.0 * ta_xxxx_xxxyz_1[i] * fe_0 + ta_xxxx_xxxyzz_0[i] * pa_z[i] - ta_xxxx_xxxyzz_1[i] * pc_z[i];

        ta_xxxxz_xxxzzz_0[i] =
            3.0 * ta_xxxx_xxxzz_0[i] * fe_0 - 3.0 * ta_xxxx_xxxzz_1[i] * fe_0 + ta_xxxx_xxxzzz_0[i] * pa_z[i] - ta_xxxx_xxxzzz_1[i] * pc_z[i];

        ta_xxxxz_xxyyyy_0[i] = ta_xxxx_xxyyyy_0[i] * pa_z[i] - ta_xxxx_xxyyyy_1[i] * pc_z[i];

        ta_xxxxz_xxyyyz_0[i] = ta_xxxx_xxyyy_0[i] * fe_0 - ta_xxxx_xxyyy_1[i] * fe_0 + ta_xxxx_xxyyyz_0[i] * pa_z[i] - ta_xxxx_xxyyyz_1[i] * pc_z[i];

        ta_xxxxz_xxyyzz_0[i] =
            2.0 * ta_xxxx_xxyyz_0[i] * fe_0 - 2.0 * ta_xxxx_xxyyz_1[i] * fe_0 + ta_xxxx_xxyyzz_0[i] * pa_z[i] - ta_xxxx_xxyyzz_1[i] * pc_z[i];

        ta_xxxxz_xxyzzz_0[i] =
            3.0 * ta_xxxx_xxyzz_0[i] * fe_0 - 3.0 * ta_xxxx_xxyzz_1[i] * fe_0 + ta_xxxx_xxyzzz_0[i] * pa_z[i] - ta_xxxx_xxyzzz_1[i] * pc_z[i];

        ta_xxxxz_xxzzzz_0[i] =
            4.0 * ta_xxxx_xxzzz_0[i] * fe_0 - 4.0 * ta_xxxx_xxzzz_1[i] * fe_0 + ta_xxxx_xxzzzz_0[i] * pa_z[i] - ta_xxxx_xxzzzz_1[i] * pc_z[i];

        ta_xxxxz_xyyyyy_0[i] = ta_xxxx_xyyyyy_0[i] * pa_z[i] - ta_xxxx_xyyyyy_1[i] * pc_z[i];

        ta_xxxxz_xyyyyz_0[i] = ta_xxxx_xyyyy_0[i] * fe_0 - ta_xxxx_xyyyy_1[i] * fe_0 + ta_xxxx_xyyyyz_0[i] * pa_z[i] - ta_xxxx_xyyyyz_1[i] * pc_z[i];

        ta_xxxxz_xyyyzz_0[i] =
            2.0 * ta_xxxx_xyyyz_0[i] * fe_0 - 2.0 * ta_xxxx_xyyyz_1[i] * fe_0 + ta_xxxx_xyyyzz_0[i] * pa_z[i] - ta_xxxx_xyyyzz_1[i] * pc_z[i];

        ta_xxxxz_xyyzzz_0[i] =
            3.0 * ta_xxxx_xyyzz_0[i] * fe_0 - 3.0 * ta_xxxx_xyyzz_1[i] * fe_0 + ta_xxxx_xyyzzz_0[i] * pa_z[i] - ta_xxxx_xyyzzz_1[i] * pc_z[i];

        ta_xxxxz_xyzzzz_0[i] =
            4.0 * ta_xxxx_xyzzz_0[i] * fe_0 - 4.0 * ta_xxxx_xyzzz_1[i] * fe_0 + ta_xxxx_xyzzzz_0[i] * pa_z[i] - ta_xxxx_xyzzzz_1[i] * pc_z[i];

        ta_xxxxz_xzzzzz_0[i] =
            5.0 * ta_xxxx_xzzzz_0[i] * fe_0 - 5.0 * ta_xxxx_xzzzz_1[i] * fe_0 + ta_xxxx_xzzzzz_0[i] * pa_z[i] - ta_xxxx_xzzzzz_1[i] * pc_z[i];

        ta_xxxxz_yyyyyy_0[i] = ta_xxxx_yyyyyy_0[i] * pa_z[i] - ta_xxxx_yyyyyy_1[i] * pc_z[i];

        ta_xxxxz_yyyyyz_0[i] =
            3.0 * ta_xxz_yyyyyz_0[i] * fe_0 - 3.0 * ta_xxz_yyyyyz_1[i] * fe_0 + ta_xxxz_yyyyyz_0[i] * pa_x[i] - ta_xxxz_yyyyyz_1[i] * pc_x[i];

        ta_xxxxz_yyyyzz_0[i] =
            3.0 * ta_xxz_yyyyzz_0[i] * fe_0 - 3.0 * ta_xxz_yyyyzz_1[i] * fe_0 + ta_xxxz_yyyyzz_0[i] * pa_x[i] - ta_xxxz_yyyyzz_1[i] * pc_x[i];

        ta_xxxxz_yyyzzz_0[i] =
            3.0 * ta_xxz_yyyzzz_0[i] * fe_0 - 3.0 * ta_xxz_yyyzzz_1[i] * fe_0 + ta_xxxz_yyyzzz_0[i] * pa_x[i] - ta_xxxz_yyyzzz_1[i] * pc_x[i];

        ta_xxxxz_yyzzzz_0[i] =
            3.0 * ta_xxz_yyzzzz_0[i] * fe_0 - 3.0 * ta_xxz_yyzzzz_1[i] * fe_0 + ta_xxxz_yyzzzz_0[i] * pa_x[i] - ta_xxxz_yyzzzz_1[i] * pc_x[i];

        ta_xxxxz_yzzzzz_0[i] =
            3.0 * ta_xxz_yzzzzz_0[i] * fe_0 - 3.0 * ta_xxz_yzzzzz_1[i] * fe_0 + ta_xxxz_yzzzzz_0[i] * pa_x[i] - ta_xxxz_yzzzzz_1[i] * pc_x[i];

        ta_xxxxz_zzzzzz_0[i] =
            3.0 * ta_xxz_zzzzzz_0[i] * fe_0 - 3.0 * ta_xxz_zzzzzz_1[i] * fe_0 + ta_xxxz_zzzzzz_0[i] * pa_x[i] - ta_xxxz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 84-112 components of targeted buffer : HI

    auto ta_xxxyy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 84);

    auto ta_xxxyy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 85);

    auto ta_xxxyy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 86);

    auto ta_xxxyy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 87);

    auto ta_xxxyy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 88);

    auto ta_xxxyy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 89);

    auto ta_xxxyy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 90);

    auto ta_xxxyy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 91);

    auto ta_xxxyy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 92);

    auto ta_xxxyy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 93);

    auto ta_xxxyy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 94);

    auto ta_xxxyy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 95);

    auto ta_xxxyy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 96);

    auto ta_xxxyy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 97);

    auto ta_xxxyy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 98);

    auto ta_xxxyy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 99);

    auto ta_xxxyy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 100);

    auto ta_xxxyy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 101);

    auto ta_xxxyy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 102);

    auto ta_xxxyy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 103);

    auto ta_xxxyy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 104);

    auto ta_xxxyy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 105);

    auto ta_xxxyy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 106);

    auto ta_xxxyy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 107);

    auto ta_xxxyy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 108);

    auto ta_xxxyy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 109);

    auto ta_xxxyy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 110);

    auto ta_xxxyy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 111);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxx_xxxxxx_0,   \
                             ta_xxx_xxxxxx_1,   \
                             ta_xxx_xxxxxz_0,   \
                             ta_xxx_xxxxxz_1,   \
                             ta_xxx_xxxxzz_0,   \
                             ta_xxx_xxxxzz_1,   \
                             ta_xxx_xxxzzz_0,   \
                             ta_xxx_xxxzzz_1,   \
                             ta_xxx_xxzzzz_0,   \
                             ta_xxx_xxzzzz_1,   \
                             ta_xxx_xzzzzz_0,   \
                             ta_xxx_xzzzzz_1,   \
                             ta_xxxy_xxxxxx_0,  \
                             ta_xxxy_xxxxxx_1,  \
                             ta_xxxy_xxxxxz_0,  \
                             ta_xxxy_xxxxxz_1,  \
                             ta_xxxy_xxxxzz_0,  \
                             ta_xxxy_xxxxzz_1,  \
                             ta_xxxy_xxxzzz_0,  \
                             ta_xxxy_xxxzzz_1,  \
                             ta_xxxy_xxzzzz_0,  \
                             ta_xxxy_xxzzzz_1,  \
                             ta_xxxy_xzzzzz_0,  \
                             ta_xxxy_xzzzzz_1,  \
                             ta_xxxyy_xxxxxx_0, \
                             ta_xxxyy_xxxxxy_0, \
                             ta_xxxyy_xxxxxz_0, \
                             ta_xxxyy_xxxxyy_0, \
                             ta_xxxyy_xxxxyz_0, \
                             ta_xxxyy_xxxxzz_0, \
                             ta_xxxyy_xxxyyy_0, \
                             ta_xxxyy_xxxyyz_0, \
                             ta_xxxyy_xxxyzz_0, \
                             ta_xxxyy_xxxzzz_0, \
                             ta_xxxyy_xxyyyy_0, \
                             ta_xxxyy_xxyyyz_0, \
                             ta_xxxyy_xxyyzz_0, \
                             ta_xxxyy_xxyzzz_0, \
                             ta_xxxyy_xxzzzz_0, \
                             ta_xxxyy_xyyyyy_0, \
                             ta_xxxyy_xyyyyz_0, \
                             ta_xxxyy_xyyyzz_0, \
                             ta_xxxyy_xyyzzz_0, \
                             ta_xxxyy_xyzzzz_0, \
                             ta_xxxyy_xzzzzz_0, \
                             ta_xxxyy_yyyyyy_0, \
                             ta_xxxyy_yyyyyz_0, \
                             ta_xxxyy_yyyyzz_0, \
                             ta_xxxyy_yyyzzz_0, \
                             ta_xxxyy_yyzzzz_0, \
                             ta_xxxyy_yzzzzz_0, \
                             ta_xxxyy_zzzzzz_0, \
                             ta_xxyy_xxxxxy_0,  \
                             ta_xxyy_xxxxxy_1,  \
                             ta_xxyy_xxxxy_0,   \
                             ta_xxyy_xxxxy_1,   \
                             ta_xxyy_xxxxyy_0,  \
                             ta_xxyy_xxxxyy_1,  \
                             ta_xxyy_xxxxyz_0,  \
                             ta_xxyy_xxxxyz_1,  \
                             ta_xxyy_xxxyy_0,   \
                             ta_xxyy_xxxyy_1,   \
                             ta_xxyy_xxxyyy_0,  \
                             ta_xxyy_xxxyyy_1,  \
                             ta_xxyy_xxxyyz_0,  \
                             ta_xxyy_xxxyyz_1,  \
                             ta_xxyy_xxxyz_0,   \
                             ta_xxyy_xxxyz_1,   \
                             ta_xxyy_xxxyzz_0,  \
                             ta_xxyy_xxxyzz_1,  \
                             ta_xxyy_xxyyy_0,   \
                             ta_xxyy_xxyyy_1,   \
                             ta_xxyy_xxyyyy_0,  \
                             ta_xxyy_xxyyyy_1,  \
                             ta_xxyy_xxyyyz_0,  \
                             ta_xxyy_xxyyyz_1,  \
                             ta_xxyy_xxyyz_0,   \
                             ta_xxyy_xxyyz_1,   \
                             ta_xxyy_xxyyzz_0,  \
                             ta_xxyy_xxyyzz_1,  \
                             ta_xxyy_xxyzz_0,   \
                             ta_xxyy_xxyzz_1,   \
                             ta_xxyy_xxyzzz_0,  \
                             ta_xxyy_xxyzzz_1,  \
                             ta_xxyy_xyyyy_0,   \
                             ta_xxyy_xyyyy_1,   \
                             ta_xxyy_xyyyyy_0,  \
                             ta_xxyy_xyyyyy_1,  \
                             ta_xxyy_xyyyyz_0,  \
                             ta_xxyy_xyyyyz_1,  \
                             ta_xxyy_xyyyz_0,   \
                             ta_xxyy_xyyyz_1,   \
                             ta_xxyy_xyyyzz_0,  \
                             ta_xxyy_xyyyzz_1,  \
                             ta_xxyy_xyyzz_0,   \
                             ta_xxyy_xyyzz_1,   \
                             ta_xxyy_xyyzzz_0,  \
                             ta_xxyy_xyyzzz_1,  \
                             ta_xxyy_xyzzz_0,   \
                             ta_xxyy_xyzzz_1,   \
                             ta_xxyy_xyzzzz_0,  \
                             ta_xxyy_xyzzzz_1,  \
                             ta_xxyy_yyyyy_0,   \
                             ta_xxyy_yyyyy_1,   \
                             ta_xxyy_yyyyyy_0,  \
                             ta_xxyy_yyyyyy_1,  \
                             ta_xxyy_yyyyyz_0,  \
                             ta_xxyy_yyyyyz_1,  \
                             ta_xxyy_yyyyz_0,   \
                             ta_xxyy_yyyyz_1,   \
                             ta_xxyy_yyyyzz_0,  \
                             ta_xxyy_yyyyzz_1,  \
                             ta_xxyy_yyyzz_0,   \
                             ta_xxyy_yyyzz_1,   \
                             ta_xxyy_yyyzzz_0,  \
                             ta_xxyy_yyyzzz_1,  \
                             ta_xxyy_yyzzz_0,   \
                             ta_xxyy_yyzzz_1,   \
                             ta_xxyy_yyzzzz_0,  \
                             ta_xxyy_yyzzzz_1,  \
                             ta_xxyy_yzzzz_0,   \
                             ta_xxyy_yzzzz_1,   \
                             ta_xxyy_yzzzzz_0,  \
                             ta_xxyy_yzzzzz_1,  \
                             ta_xxyy_zzzzzz_0,  \
                             ta_xxyy_zzzzzz_1,  \
                             ta_xyy_xxxxxy_0,   \
                             ta_xyy_xxxxxy_1,   \
                             ta_xyy_xxxxyy_0,   \
                             ta_xyy_xxxxyy_1,   \
                             ta_xyy_xxxxyz_0,   \
                             ta_xyy_xxxxyz_1,   \
                             ta_xyy_xxxyyy_0,   \
                             ta_xyy_xxxyyy_1,   \
                             ta_xyy_xxxyyz_0,   \
                             ta_xyy_xxxyyz_1,   \
                             ta_xyy_xxxyzz_0,   \
                             ta_xyy_xxxyzz_1,   \
                             ta_xyy_xxyyyy_0,   \
                             ta_xyy_xxyyyy_1,   \
                             ta_xyy_xxyyyz_0,   \
                             ta_xyy_xxyyyz_1,   \
                             ta_xyy_xxyyzz_0,   \
                             ta_xyy_xxyyzz_1,   \
                             ta_xyy_xxyzzz_0,   \
                             ta_xyy_xxyzzz_1,   \
                             ta_xyy_xyyyyy_0,   \
                             ta_xyy_xyyyyy_1,   \
                             ta_xyy_xyyyyz_0,   \
                             ta_xyy_xyyyyz_1,   \
                             ta_xyy_xyyyzz_0,   \
                             ta_xyy_xyyyzz_1,   \
                             ta_xyy_xyyzzz_0,   \
                             ta_xyy_xyyzzz_1,   \
                             ta_xyy_xyzzzz_0,   \
                             ta_xyy_xyzzzz_1,   \
                             ta_xyy_yyyyyy_0,   \
                             ta_xyy_yyyyyy_1,   \
                             ta_xyy_yyyyyz_0,   \
                             ta_xyy_yyyyyz_1,   \
                             ta_xyy_yyyyzz_0,   \
                             ta_xyy_yyyyzz_1,   \
                             ta_xyy_yyyzzz_0,   \
                             ta_xyy_yyyzzz_1,   \
                             ta_xyy_yyzzzz_0,   \
                             ta_xyy_yyzzzz_1,   \
                             ta_xyy_yzzzzz_0,   \
                             ta_xyy_yzzzzz_1,   \
                             ta_xyy_zzzzzz_0,   \
                             ta_xyy_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyy_xxxxxx_0[i] = ta_xxx_xxxxxx_0[i] * fe_0 - ta_xxx_xxxxxx_1[i] * fe_0 + ta_xxxy_xxxxxx_0[i] * pa_y[i] - ta_xxxy_xxxxxx_1[i] * pc_y[i];

        ta_xxxyy_xxxxxy_0[i] = 2.0 * ta_xyy_xxxxxy_0[i] * fe_0 - 2.0 * ta_xyy_xxxxxy_1[i] * fe_0 + 5.0 * ta_xxyy_xxxxy_0[i] * fe_0 -
                               5.0 * ta_xxyy_xxxxy_1[i] * fe_0 + ta_xxyy_xxxxxy_0[i] * pa_x[i] - ta_xxyy_xxxxxy_1[i] * pc_x[i];

        ta_xxxyy_xxxxxz_0[i] = ta_xxx_xxxxxz_0[i] * fe_0 - ta_xxx_xxxxxz_1[i] * fe_0 + ta_xxxy_xxxxxz_0[i] * pa_y[i] - ta_xxxy_xxxxxz_1[i] * pc_y[i];

        ta_xxxyy_xxxxyy_0[i] = 2.0 * ta_xyy_xxxxyy_0[i] * fe_0 - 2.0 * ta_xyy_xxxxyy_1[i] * fe_0 + 4.0 * ta_xxyy_xxxyy_0[i] * fe_0 -
                               4.0 * ta_xxyy_xxxyy_1[i] * fe_0 + ta_xxyy_xxxxyy_0[i] * pa_x[i] - ta_xxyy_xxxxyy_1[i] * pc_x[i];

        ta_xxxyy_xxxxyz_0[i] = 2.0 * ta_xyy_xxxxyz_0[i] * fe_0 - 2.0 * ta_xyy_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxyy_xxxyz_0[i] * fe_0 -
                               4.0 * ta_xxyy_xxxyz_1[i] * fe_0 + ta_xxyy_xxxxyz_0[i] * pa_x[i] - ta_xxyy_xxxxyz_1[i] * pc_x[i];

        ta_xxxyy_xxxxzz_0[i] = ta_xxx_xxxxzz_0[i] * fe_0 - ta_xxx_xxxxzz_1[i] * fe_0 + ta_xxxy_xxxxzz_0[i] * pa_y[i] - ta_xxxy_xxxxzz_1[i] * pc_y[i];

        ta_xxxyy_xxxyyy_0[i] = 2.0 * ta_xyy_xxxyyy_0[i] * fe_0 - 2.0 * ta_xyy_xxxyyy_1[i] * fe_0 + 3.0 * ta_xxyy_xxyyy_0[i] * fe_0 -
                               3.0 * ta_xxyy_xxyyy_1[i] * fe_0 + ta_xxyy_xxxyyy_0[i] * pa_x[i] - ta_xxyy_xxxyyy_1[i] * pc_x[i];

        ta_xxxyy_xxxyyz_0[i] = 2.0 * ta_xyy_xxxyyz_0[i] * fe_0 - 2.0 * ta_xyy_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxyy_xxyyz_0[i] * fe_0 -
                               3.0 * ta_xxyy_xxyyz_1[i] * fe_0 + ta_xxyy_xxxyyz_0[i] * pa_x[i] - ta_xxyy_xxxyyz_1[i] * pc_x[i];

        ta_xxxyy_xxxyzz_0[i] = 2.0 * ta_xyy_xxxyzz_0[i] * fe_0 - 2.0 * ta_xyy_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxyy_xxyzz_0[i] * fe_0 -
                               3.0 * ta_xxyy_xxyzz_1[i] * fe_0 + ta_xxyy_xxxyzz_0[i] * pa_x[i] - ta_xxyy_xxxyzz_1[i] * pc_x[i];

        ta_xxxyy_xxxzzz_0[i] = ta_xxx_xxxzzz_0[i] * fe_0 - ta_xxx_xxxzzz_1[i] * fe_0 + ta_xxxy_xxxzzz_0[i] * pa_y[i] - ta_xxxy_xxxzzz_1[i] * pc_y[i];

        ta_xxxyy_xxyyyy_0[i] = 2.0 * ta_xyy_xxyyyy_0[i] * fe_0 - 2.0 * ta_xyy_xxyyyy_1[i] * fe_0 + 2.0 * ta_xxyy_xyyyy_0[i] * fe_0 -
                               2.0 * ta_xxyy_xyyyy_1[i] * fe_0 + ta_xxyy_xxyyyy_0[i] * pa_x[i] - ta_xxyy_xxyyyy_1[i] * pc_x[i];

        ta_xxxyy_xxyyyz_0[i] = 2.0 * ta_xyy_xxyyyz_0[i] * fe_0 - 2.0 * ta_xyy_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxyy_xyyyz_0[i] * fe_0 -
                               2.0 * ta_xxyy_xyyyz_1[i] * fe_0 + ta_xxyy_xxyyyz_0[i] * pa_x[i] - ta_xxyy_xxyyyz_1[i] * pc_x[i];

        ta_xxxyy_xxyyzz_0[i] = 2.0 * ta_xyy_xxyyzz_0[i] * fe_0 - 2.0 * ta_xyy_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxyy_xyyzz_0[i] * fe_0 -
                               2.0 * ta_xxyy_xyyzz_1[i] * fe_0 + ta_xxyy_xxyyzz_0[i] * pa_x[i] - ta_xxyy_xxyyzz_1[i] * pc_x[i];

        ta_xxxyy_xxyzzz_0[i] = 2.0 * ta_xyy_xxyzzz_0[i] * fe_0 - 2.0 * ta_xyy_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxyy_xyzzz_0[i] * fe_0 -
                               2.0 * ta_xxyy_xyzzz_1[i] * fe_0 + ta_xxyy_xxyzzz_0[i] * pa_x[i] - ta_xxyy_xxyzzz_1[i] * pc_x[i];

        ta_xxxyy_xxzzzz_0[i] = ta_xxx_xxzzzz_0[i] * fe_0 - ta_xxx_xxzzzz_1[i] * fe_0 + ta_xxxy_xxzzzz_0[i] * pa_y[i] - ta_xxxy_xxzzzz_1[i] * pc_y[i];

        ta_xxxyy_xyyyyy_0[i] = 2.0 * ta_xyy_xyyyyy_0[i] * fe_0 - 2.0 * ta_xyy_xyyyyy_1[i] * fe_0 + ta_xxyy_yyyyy_0[i] * fe_0 -
                               ta_xxyy_yyyyy_1[i] * fe_0 + ta_xxyy_xyyyyy_0[i] * pa_x[i] - ta_xxyy_xyyyyy_1[i] * pc_x[i];

        ta_xxxyy_xyyyyz_0[i] = 2.0 * ta_xyy_xyyyyz_0[i] * fe_0 - 2.0 * ta_xyy_xyyyyz_1[i] * fe_0 + ta_xxyy_yyyyz_0[i] * fe_0 -
                               ta_xxyy_yyyyz_1[i] * fe_0 + ta_xxyy_xyyyyz_0[i] * pa_x[i] - ta_xxyy_xyyyyz_1[i] * pc_x[i];

        ta_xxxyy_xyyyzz_0[i] = 2.0 * ta_xyy_xyyyzz_0[i] * fe_0 - 2.0 * ta_xyy_xyyyzz_1[i] * fe_0 + ta_xxyy_yyyzz_0[i] * fe_0 -
                               ta_xxyy_yyyzz_1[i] * fe_0 + ta_xxyy_xyyyzz_0[i] * pa_x[i] - ta_xxyy_xyyyzz_1[i] * pc_x[i];

        ta_xxxyy_xyyzzz_0[i] = 2.0 * ta_xyy_xyyzzz_0[i] * fe_0 - 2.0 * ta_xyy_xyyzzz_1[i] * fe_0 + ta_xxyy_yyzzz_0[i] * fe_0 -
                               ta_xxyy_yyzzz_1[i] * fe_0 + ta_xxyy_xyyzzz_0[i] * pa_x[i] - ta_xxyy_xyyzzz_1[i] * pc_x[i];

        ta_xxxyy_xyzzzz_0[i] = 2.0 * ta_xyy_xyzzzz_0[i] * fe_0 - 2.0 * ta_xyy_xyzzzz_1[i] * fe_0 + ta_xxyy_yzzzz_0[i] * fe_0 -
                               ta_xxyy_yzzzz_1[i] * fe_0 + ta_xxyy_xyzzzz_0[i] * pa_x[i] - ta_xxyy_xyzzzz_1[i] * pc_x[i];

        ta_xxxyy_xzzzzz_0[i] = ta_xxx_xzzzzz_0[i] * fe_0 - ta_xxx_xzzzzz_1[i] * fe_0 + ta_xxxy_xzzzzz_0[i] * pa_y[i] - ta_xxxy_xzzzzz_1[i] * pc_y[i];

        ta_xxxyy_yyyyyy_0[i] =
            2.0 * ta_xyy_yyyyyy_0[i] * fe_0 - 2.0 * ta_xyy_yyyyyy_1[i] * fe_0 + ta_xxyy_yyyyyy_0[i] * pa_x[i] - ta_xxyy_yyyyyy_1[i] * pc_x[i];

        ta_xxxyy_yyyyyz_0[i] =
            2.0 * ta_xyy_yyyyyz_0[i] * fe_0 - 2.0 * ta_xyy_yyyyyz_1[i] * fe_0 + ta_xxyy_yyyyyz_0[i] * pa_x[i] - ta_xxyy_yyyyyz_1[i] * pc_x[i];

        ta_xxxyy_yyyyzz_0[i] =
            2.0 * ta_xyy_yyyyzz_0[i] * fe_0 - 2.0 * ta_xyy_yyyyzz_1[i] * fe_0 + ta_xxyy_yyyyzz_0[i] * pa_x[i] - ta_xxyy_yyyyzz_1[i] * pc_x[i];

        ta_xxxyy_yyyzzz_0[i] =
            2.0 * ta_xyy_yyyzzz_0[i] * fe_0 - 2.0 * ta_xyy_yyyzzz_1[i] * fe_0 + ta_xxyy_yyyzzz_0[i] * pa_x[i] - ta_xxyy_yyyzzz_1[i] * pc_x[i];

        ta_xxxyy_yyzzzz_0[i] =
            2.0 * ta_xyy_yyzzzz_0[i] * fe_0 - 2.0 * ta_xyy_yyzzzz_1[i] * fe_0 + ta_xxyy_yyzzzz_0[i] * pa_x[i] - ta_xxyy_yyzzzz_1[i] * pc_x[i];

        ta_xxxyy_yzzzzz_0[i] =
            2.0 * ta_xyy_yzzzzz_0[i] * fe_0 - 2.0 * ta_xyy_yzzzzz_1[i] * fe_0 + ta_xxyy_yzzzzz_0[i] * pa_x[i] - ta_xxyy_yzzzzz_1[i] * pc_x[i];

        ta_xxxyy_zzzzzz_0[i] =
            2.0 * ta_xyy_zzzzzz_0[i] * fe_0 - 2.0 * ta_xyy_zzzzzz_1[i] * fe_0 + ta_xxyy_zzzzzz_0[i] * pa_x[i] - ta_xxyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 112-140 components of targeted buffer : HI

    auto ta_xxxyz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 112);

    auto ta_xxxyz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 113);

    auto ta_xxxyz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 114);

    auto ta_xxxyz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 115);

    auto ta_xxxyz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 116);

    auto ta_xxxyz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 117);

    auto ta_xxxyz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 118);

    auto ta_xxxyz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 119);

    auto ta_xxxyz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 120);

    auto ta_xxxyz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 121);

    auto ta_xxxyz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 122);

    auto ta_xxxyz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 123);

    auto ta_xxxyz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 124);

    auto ta_xxxyz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 125);

    auto ta_xxxyz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 126);

    auto ta_xxxyz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 127);

    auto ta_xxxyz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 128);

    auto ta_xxxyz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 129);

    auto ta_xxxyz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 130);

    auto ta_xxxyz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 131);

    auto ta_xxxyz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 132);

    auto ta_xxxyz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 133);

    auto ta_xxxyz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 134);

    auto ta_xxxyz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 135);

    auto ta_xxxyz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 136);

    auto ta_xxxyz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 137);

    auto ta_xxxyz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 138);

    auto ta_xxxyz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 139);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta_xxxy_xxxxxy_0,  \
                             ta_xxxy_xxxxxy_1,  \
                             ta_xxxy_xxxxyy_0,  \
                             ta_xxxy_xxxxyy_1,  \
                             ta_xxxy_xxxyyy_0,  \
                             ta_xxxy_xxxyyy_1,  \
                             ta_xxxy_xxyyyy_0,  \
                             ta_xxxy_xxyyyy_1,  \
                             ta_xxxy_xyyyyy_0,  \
                             ta_xxxy_xyyyyy_1,  \
                             ta_xxxy_yyyyyy_0,  \
                             ta_xxxy_yyyyyy_1,  \
                             ta_xxxyz_xxxxxx_0, \
                             ta_xxxyz_xxxxxy_0, \
                             ta_xxxyz_xxxxxz_0, \
                             ta_xxxyz_xxxxyy_0, \
                             ta_xxxyz_xxxxyz_0, \
                             ta_xxxyz_xxxxzz_0, \
                             ta_xxxyz_xxxyyy_0, \
                             ta_xxxyz_xxxyyz_0, \
                             ta_xxxyz_xxxyzz_0, \
                             ta_xxxyz_xxxzzz_0, \
                             ta_xxxyz_xxyyyy_0, \
                             ta_xxxyz_xxyyyz_0, \
                             ta_xxxyz_xxyyzz_0, \
                             ta_xxxyz_xxyzzz_0, \
                             ta_xxxyz_xxzzzz_0, \
                             ta_xxxyz_xyyyyy_0, \
                             ta_xxxyz_xyyyyz_0, \
                             ta_xxxyz_xyyyzz_0, \
                             ta_xxxyz_xyyzzz_0, \
                             ta_xxxyz_xyzzzz_0, \
                             ta_xxxyz_xzzzzz_0, \
                             ta_xxxyz_yyyyyy_0, \
                             ta_xxxyz_yyyyyz_0, \
                             ta_xxxyz_yyyyzz_0, \
                             ta_xxxyz_yyyzzz_0, \
                             ta_xxxyz_yyzzzz_0, \
                             ta_xxxyz_yzzzzz_0, \
                             ta_xxxyz_zzzzzz_0, \
                             ta_xxxz_xxxxxx_0,  \
                             ta_xxxz_xxxxxx_1,  \
                             ta_xxxz_xxxxxz_0,  \
                             ta_xxxz_xxxxxz_1,  \
                             ta_xxxz_xxxxyz_0,  \
                             ta_xxxz_xxxxyz_1,  \
                             ta_xxxz_xxxxz_0,   \
                             ta_xxxz_xxxxz_1,   \
                             ta_xxxz_xxxxzz_0,  \
                             ta_xxxz_xxxxzz_1,  \
                             ta_xxxz_xxxyyz_0,  \
                             ta_xxxz_xxxyyz_1,  \
                             ta_xxxz_xxxyz_0,   \
                             ta_xxxz_xxxyz_1,   \
                             ta_xxxz_xxxyzz_0,  \
                             ta_xxxz_xxxyzz_1,  \
                             ta_xxxz_xxxzz_0,   \
                             ta_xxxz_xxxzz_1,   \
                             ta_xxxz_xxxzzz_0,  \
                             ta_xxxz_xxxzzz_1,  \
                             ta_xxxz_xxyyyz_0,  \
                             ta_xxxz_xxyyyz_1,  \
                             ta_xxxz_xxyyz_0,   \
                             ta_xxxz_xxyyz_1,   \
                             ta_xxxz_xxyyzz_0,  \
                             ta_xxxz_xxyyzz_1,  \
                             ta_xxxz_xxyzz_0,   \
                             ta_xxxz_xxyzz_1,   \
                             ta_xxxz_xxyzzz_0,  \
                             ta_xxxz_xxyzzz_1,  \
                             ta_xxxz_xxzzz_0,   \
                             ta_xxxz_xxzzz_1,   \
                             ta_xxxz_xxzzzz_0,  \
                             ta_xxxz_xxzzzz_1,  \
                             ta_xxxz_xyyyyz_0,  \
                             ta_xxxz_xyyyyz_1,  \
                             ta_xxxz_xyyyz_0,   \
                             ta_xxxz_xyyyz_1,   \
                             ta_xxxz_xyyyzz_0,  \
                             ta_xxxz_xyyyzz_1,  \
                             ta_xxxz_xyyzz_0,   \
                             ta_xxxz_xyyzz_1,   \
                             ta_xxxz_xyyzzz_0,  \
                             ta_xxxz_xyyzzz_1,  \
                             ta_xxxz_xyzzz_0,   \
                             ta_xxxz_xyzzz_1,   \
                             ta_xxxz_xyzzzz_0,  \
                             ta_xxxz_xyzzzz_1,  \
                             ta_xxxz_xzzzz_0,   \
                             ta_xxxz_xzzzz_1,   \
                             ta_xxxz_xzzzzz_0,  \
                             ta_xxxz_xzzzzz_1,  \
                             ta_xxxz_zzzzzz_0,  \
                             ta_xxxz_zzzzzz_1,  \
                             ta_xxyz_yyyyyz_0,  \
                             ta_xxyz_yyyyyz_1,  \
                             ta_xxyz_yyyyzz_0,  \
                             ta_xxyz_yyyyzz_1,  \
                             ta_xxyz_yyyzzz_0,  \
                             ta_xxyz_yyyzzz_1,  \
                             ta_xxyz_yyzzzz_0,  \
                             ta_xxyz_yyzzzz_1,  \
                             ta_xxyz_yzzzzz_0,  \
                             ta_xxyz_yzzzzz_1,  \
                             ta_xyz_yyyyyz_0,   \
                             ta_xyz_yyyyyz_1,   \
                             ta_xyz_yyyyzz_0,   \
                             ta_xyz_yyyyzz_1,   \
                             ta_xyz_yyyzzz_0,   \
                             ta_xyz_yyyzzz_1,   \
                             ta_xyz_yyzzzz_0,   \
                             ta_xyz_yyzzzz_1,   \
                             ta_xyz_yzzzzz_0,   \
                             ta_xyz_yzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyz_xxxxxx_0[i] = ta_xxxz_xxxxxx_0[i] * pa_y[i] - ta_xxxz_xxxxxx_1[i] * pc_y[i];

        ta_xxxyz_xxxxxy_0[i] = ta_xxxy_xxxxxy_0[i] * pa_z[i] - ta_xxxy_xxxxxy_1[i] * pc_z[i];

        ta_xxxyz_xxxxxz_0[i] = ta_xxxz_xxxxxz_0[i] * pa_y[i] - ta_xxxz_xxxxxz_1[i] * pc_y[i];

        ta_xxxyz_xxxxyy_0[i] = ta_xxxy_xxxxyy_0[i] * pa_z[i] - ta_xxxy_xxxxyy_1[i] * pc_z[i];

        ta_xxxyz_xxxxyz_0[i] = ta_xxxz_xxxxz_0[i] * fe_0 - ta_xxxz_xxxxz_1[i] * fe_0 + ta_xxxz_xxxxyz_0[i] * pa_y[i] - ta_xxxz_xxxxyz_1[i] * pc_y[i];

        ta_xxxyz_xxxxzz_0[i] = ta_xxxz_xxxxzz_0[i] * pa_y[i] - ta_xxxz_xxxxzz_1[i] * pc_y[i];

        ta_xxxyz_xxxyyy_0[i] = ta_xxxy_xxxyyy_0[i] * pa_z[i] - ta_xxxy_xxxyyy_1[i] * pc_z[i];

        ta_xxxyz_xxxyyz_0[i] =
            2.0 * ta_xxxz_xxxyz_0[i] * fe_0 - 2.0 * ta_xxxz_xxxyz_1[i] * fe_0 + ta_xxxz_xxxyyz_0[i] * pa_y[i] - ta_xxxz_xxxyyz_1[i] * pc_y[i];

        ta_xxxyz_xxxyzz_0[i] = ta_xxxz_xxxzz_0[i] * fe_0 - ta_xxxz_xxxzz_1[i] * fe_0 + ta_xxxz_xxxyzz_0[i] * pa_y[i] - ta_xxxz_xxxyzz_1[i] * pc_y[i];

        ta_xxxyz_xxxzzz_0[i] = ta_xxxz_xxxzzz_0[i] * pa_y[i] - ta_xxxz_xxxzzz_1[i] * pc_y[i];

        ta_xxxyz_xxyyyy_0[i] = ta_xxxy_xxyyyy_0[i] * pa_z[i] - ta_xxxy_xxyyyy_1[i] * pc_z[i];

        ta_xxxyz_xxyyyz_0[i] =
            3.0 * ta_xxxz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxxz_xxyyz_1[i] * fe_0 + ta_xxxz_xxyyyz_0[i] * pa_y[i] - ta_xxxz_xxyyyz_1[i] * pc_y[i];

        ta_xxxyz_xxyyzz_0[i] =
            2.0 * ta_xxxz_xxyzz_0[i] * fe_0 - 2.0 * ta_xxxz_xxyzz_1[i] * fe_0 + ta_xxxz_xxyyzz_0[i] * pa_y[i] - ta_xxxz_xxyyzz_1[i] * pc_y[i];

        ta_xxxyz_xxyzzz_0[i] = ta_xxxz_xxzzz_0[i] * fe_0 - ta_xxxz_xxzzz_1[i] * fe_0 + ta_xxxz_xxyzzz_0[i] * pa_y[i] - ta_xxxz_xxyzzz_1[i] * pc_y[i];

        ta_xxxyz_xxzzzz_0[i] = ta_xxxz_xxzzzz_0[i] * pa_y[i] - ta_xxxz_xxzzzz_1[i] * pc_y[i];

        ta_xxxyz_xyyyyy_0[i] = ta_xxxy_xyyyyy_0[i] * pa_z[i] - ta_xxxy_xyyyyy_1[i] * pc_z[i];

        ta_xxxyz_xyyyyz_0[i] =
            4.0 * ta_xxxz_xyyyz_0[i] * fe_0 - 4.0 * ta_xxxz_xyyyz_1[i] * fe_0 + ta_xxxz_xyyyyz_0[i] * pa_y[i] - ta_xxxz_xyyyyz_1[i] * pc_y[i];

        ta_xxxyz_xyyyzz_0[i] =
            3.0 * ta_xxxz_xyyzz_0[i] * fe_0 - 3.0 * ta_xxxz_xyyzz_1[i] * fe_0 + ta_xxxz_xyyyzz_0[i] * pa_y[i] - ta_xxxz_xyyyzz_1[i] * pc_y[i];

        ta_xxxyz_xyyzzz_0[i] =
            2.0 * ta_xxxz_xyzzz_0[i] * fe_0 - 2.0 * ta_xxxz_xyzzz_1[i] * fe_0 + ta_xxxz_xyyzzz_0[i] * pa_y[i] - ta_xxxz_xyyzzz_1[i] * pc_y[i];

        ta_xxxyz_xyzzzz_0[i] = ta_xxxz_xzzzz_0[i] * fe_0 - ta_xxxz_xzzzz_1[i] * fe_0 + ta_xxxz_xyzzzz_0[i] * pa_y[i] - ta_xxxz_xyzzzz_1[i] * pc_y[i];

        ta_xxxyz_xzzzzz_0[i] = ta_xxxz_xzzzzz_0[i] * pa_y[i] - ta_xxxz_xzzzzz_1[i] * pc_y[i];

        ta_xxxyz_yyyyyy_0[i] = ta_xxxy_yyyyyy_0[i] * pa_z[i] - ta_xxxy_yyyyyy_1[i] * pc_z[i];

        ta_xxxyz_yyyyyz_0[i] =
            2.0 * ta_xyz_yyyyyz_0[i] * fe_0 - 2.0 * ta_xyz_yyyyyz_1[i] * fe_0 + ta_xxyz_yyyyyz_0[i] * pa_x[i] - ta_xxyz_yyyyyz_1[i] * pc_x[i];

        ta_xxxyz_yyyyzz_0[i] =
            2.0 * ta_xyz_yyyyzz_0[i] * fe_0 - 2.0 * ta_xyz_yyyyzz_1[i] * fe_0 + ta_xxyz_yyyyzz_0[i] * pa_x[i] - ta_xxyz_yyyyzz_1[i] * pc_x[i];

        ta_xxxyz_yyyzzz_0[i] =
            2.0 * ta_xyz_yyyzzz_0[i] * fe_0 - 2.0 * ta_xyz_yyyzzz_1[i] * fe_0 + ta_xxyz_yyyzzz_0[i] * pa_x[i] - ta_xxyz_yyyzzz_1[i] * pc_x[i];

        ta_xxxyz_yyzzzz_0[i] =
            2.0 * ta_xyz_yyzzzz_0[i] * fe_0 - 2.0 * ta_xyz_yyzzzz_1[i] * fe_0 + ta_xxyz_yyzzzz_0[i] * pa_x[i] - ta_xxyz_yyzzzz_1[i] * pc_x[i];

        ta_xxxyz_yzzzzz_0[i] =
            2.0 * ta_xyz_yzzzzz_0[i] * fe_0 - 2.0 * ta_xyz_yzzzzz_1[i] * fe_0 + ta_xxyz_yzzzzz_0[i] * pa_x[i] - ta_xxyz_yzzzzz_1[i] * pc_x[i];

        ta_xxxyz_zzzzzz_0[i] = ta_xxxz_zzzzzz_0[i] * pa_y[i] - ta_xxxz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 140-168 components of targeted buffer : HI

    auto ta_xxxzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 140);

    auto ta_xxxzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 141);

    auto ta_xxxzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 142);

    auto ta_xxxzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 143);

    auto ta_xxxzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 144);

    auto ta_xxxzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 145);

    auto ta_xxxzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 146);

    auto ta_xxxzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 147);

    auto ta_xxxzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 148);

    auto ta_xxxzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 149);

    auto ta_xxxzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 150);

    auto ta_xxxzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 151);

    auto ta_xxxzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 152);

    auto ta_xxxzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 153);

    auto ta_xxxzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 154);

    auto ta_xxxzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 155);

    auto ta_xxxzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 156);

    auto ta_xxxzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 157);

    auto ta_xxxzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 158);

    auto ta_xxxzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 159);

    auto ta_xxxzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 160);

    auto ta_xxxzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 161);

    auto ta_xxxzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 162);

    auto ta_xxxzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 163);

    auto ta_xxxzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 164);

    auto ta_xxxzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 165);

    auto ta_xxxzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 166);

    auto ta_xxxzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 167);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xxx_xxxxxx_0,   \
                             ta_xxx_xxxxxx_1,   \
                             ta_xxx_xxxxxy_0,   \
                             ta_xxx_xxxxxy_1,   \
                             ta_xxx_xxxxyy_0,   \
                             ta_xxx_xxxxyy_1,   \
                             ta_xxx_xxxyyy_0,   \
                             ta_xxx_xxxyyy_1,   \
                             ta_xxx_xxyyyy_0,   \
                             ta_xxx_xxyyyy_1,   \
                             ta_xxx_xyyyyy_0,   \
                             ta_xxx_xyyyyy_1,   \
                             ta_xxxz_xxxxxx_0,  \
                             ta_xxxz_xxxxxx_1,  \
                             ta_xxxz_xxxxxy_0,  \
                             ta_xxxz_xxxxxy_1,  \
                             ta_xxxz_xxxxyy_0,  \
                             ta_xxxz_xxxxyy_1,  \
                             ta_xxxz_xxxyyy_0,  \
                             ta_xxxz_xxxyyy_1,  \
                             ta_xxxz_xxyyyy_0,  \
                             ta_xxxz_xxyyyy_1,  \
                             ta_xxxz_xyyyyy_0,  \
                             ta_xxxz_xyyyyy_1,  \
                             ta_xxxzz_xxxxxx_0, \
                             ta_xxxzz_xxxxxy_0, \
                             ta_xxxzz_xxxxxz_0, \
                             ta_xxxzz_xxxxyy_0, \
                             ta_xxxzz_xxxxyz_0, \
                             ta_xxxzz_xxxxzz_0, \
                             ta_xxxzz_xxxyyy_0, \
                             ta_xxxzz_xxxyyz_0, \
                             ta_xxxzz_xxxyzz_0, \
                             ta_xxxzz_xxxzzz_0, \
                             ta_xxxzz_xxyyyy_0, \
                             ta_xxxzz_xxyyyz_0, \
                             ta_xxxzz_xxyyzz_0, \
                             ta_xxxzz_xxyzzz_0, \
                             ta_xxxzz_xxzzzz_0, \
                             ta_xxxzz_xyyyyy_0, \
                             ta_xxxzz_xyyyyz_0, \
                             ta_xxxzz_xyyyzz_0, \
                             ta_xxxzz_xyyzzz_0, \
                             ta_xxxzz_xyzzzz_0, \
                             ta_xxxzz_xzzzzz_0, \
                             ta_xxxzz_yyyyyy_0, \
                             ta_xxxzz_yyyyyz_0, \
                             ta_xxxzz_yyyyzz_0, \
                             ta_xxxzz_yyyzzz_0, \
                             ta_xxxzz_yyzzzz_0, \
                             ta_xxxzz_yzzzzz_0, \
                             ta_xxxzz_zzzzzz_0, \
                             ta_xxzz_xxxxxz_0,  \
                             ta_xxzz_xxxxxz_1,  \
                             ta_xxzz_xxxxyz_0,  \
                             ta_xxzz_xxxxyz_1,  \
                             ta_xxzz_xxxxz_0,   \
                             ta_xxzz_xxxxz_1,   \
                             ta_xxzz_xxxxzz_0,  \
                             ta_xxzz_xxxxzz_1,  \
                             ta_xxzz_xxxyyz_0,  \
                             ta_xxzz_xxxyyz_1,  \
                             ta_xxzz_xxxyz_0,   \
                             ta_xxzz_xxxyz_1,   \
                             ta_xxzz_xxxyzz_0,  \
                             ta_xxzz_xxxyzz_1,  \
                             ta_xxzz_xxxzz_0,   \
                             ta_xxzz_xxxzz_1,   \
                             ta_xxzz_xxxzzz_0,  \
                             ta_xxzz_xxxzzz_1,  \
                             ta_xxzz_xxyyyz_0,  \
                             ta_xxzz_xxyyyz_1,  \
                             ta_xxzz_xxyyz_0,   \
                             ta_xxzz_xxyyz_1,   \
                             ta_xxzz_xxyyzz_0,  \
                             ta_xxzz_xxyyzz_1,  \
                             ta_xxzz_xxyzz_0,   \
                             ta_xxzz_xxyzz_1,   \
                             ta_xxzz_xxyzzz_0,  \
                             ta_xxzz_xxyzzz_1,  \
                             ta_xxzz_xxzzz_0,   \
                             ta_xxzz_xxzzz_1,   \
                             ta_xxzz_xxzzzz_0,  \
                             ta_xxzz_xxzzzz_1,  \
                             ta_xxzz_xyyyyz_0,  \
                             ta_xxzz_xyyyyz_1,  \
                             ta_xxzz_xyyyz_0,   \
                             ta_xxzz_xyyyz_1,   \
                             ta_xxzz_xyyyzz_0,  \
                             ta_xxzz_xyyyzz_1,  \
                             ta_xxzz_xyyzz_0,   \
                             ta_xxzz_xyyzz_1,   \
                             ta_xxzz_xyyzzz_0,  \
                             ta_xxzz_xyyzzz_1,  \
                             ta_xxzz_xyzzz_0,   \
                             ta_xxzz_xyzzz_1,   \
                             ta_xxzz_xyzzzz_0,  \
                             ta_xxzz_xyzzzz_1,  \
                             ta_xxzz_xzzzz_0,   \
                             ta_xxzz_xzzzz_1,   \
                             ta_xxzz_xzzzzz_0,  \
                             ta_xxzz_xzzzzz_1,  \
                             ta_xxzz_yyyyyy_0,  \
                             ta_xxzz_yyyyyy_1,  \
                             ta_xxzz_yyyyyz_0,  \
                             ta_xxzz_yyyyyz_1,  \
                             ta_xxzz_yyyyz_0,   \
                             ta_xxzz_yyyyz_1,   \
                             ta_xxzz_yyyyzz_0,  \
                             ta_xxzz_yyyyzz_1,  \
                             ta_xxzz_yyyzz_0,   \
                             ta_xxzz_yyyzz_1,   \
                             ta_xxzz_yyyzzz_0,  \
                             ta_xxzz_yyyzzz_1,  \
                             ta_xxzz_yyzzz_0,   \
                             ta_xxzz_yyzzz_1,   \
                             ta_xxzz_yyzzzz_0,  \
                             ta_xxzz_yyzzzz_1,  \
                             ta_xxzz_yzzzz_0,   \
                             ta_xxzz_yzzzz_1,   \
                             ta_xxzz_yzzzzz_0,  \
                             ta_xxzz_yzzzzz_1,  \
                             ta_xxzz_zzzzz_0,   \
                             ta_xxzz_zzzzz_1,   \
                             ta_xxzz_zzzzzz_0,  \
                             ta_xxzz_zzzzzz_1,  \
                             ta_xzz_xxxxxz_0,   \
                             ta_xzz_xxxxxz_1,   \
                             ta_xzz_xxxxyz_0,   \
                             ta_xzz_xxxxyz_1,   \
                             ta_xzz_xxxxzz_0,   \
                             ta_xzz_xxxxzz_1,   \
                             ta_xzz_xxxyyz_0,   \
                             ta_xzz_xxxyyz_1,   \
                             ta_xzz_xxxyzz_0,   \
                             ta_xzz_xxxyzz_1,   \
                             ta_xzz_xxxzzz_0,   \
                             ta_xzz_xxxzzz_1,   \
                             ta_xzz_xxyyyz_0,   \
                             ta_xzz_xxyyyz_1,   \
                             ta_xzz_xxyyzz_0,   \
                             ta_xzz_xxyyzz_1,   \
                             ta_xzz_xxyzzz_0,   \
                             ta_xzz_xxyzzz_1,   \
                             ta_xzz_xxzzzz_0,   \
                             ta_xzz_xxzzzz_1,   \
                             ta_xzz_xyyyyz_0,   \
                             ta_xzz_xyyyyz_1,   \
                             ta_xzz_xyyyzz_0,   \
                             ta_xzz_xyyyzz_1,   \
                             ta_xzz_xyyzzz_0,   \
                             ta_xzz_xyyzzz_1,   \
                             ta_xzz_xyzzzz_0,   \
                             ta_xzz_xyzzzz_1,   \
                             ta_xzz_xzzzzz_0,   \
                             ta_xzz_xzzzzz_1,   \
                             ta_xzz_yyyyyy_0,   \
                             ta_xzz_yyyyyy_1,   \
                             ta_xzz_yyyyyz_0,   \
                             ta_xzz_yyyyyz_1,   \
                             ta_xzz_yyyyzz_0,   \
                             ta_xzz_yyyyzz_1,   \
                             ta_xzz_yyyzzz_0,   \
                             ta_xzz_yyyzzz_1,   \
                             ta_xzz_yyzzzz_0,   \
                             ta_xzz_yyzzzz_1,   \
                             ta_xzz_yzzzzz_0,   \
                             ta_xzz_yzzzzz_1,   \
                             ta_xzz_zzzzzz_0,   \
                             ta_xzz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzz_xxxxxx_0[i] = ta_xxx_xxxxxx_0[i] * fe_0 - ta_xxx_xxxxxx_1[i] * fe_0 + ta_xxxz_xxxxxx_0[i] * pa_z[i] - ta_xxxz_xxxxxx_1[i] * pc_z[i];

        ta_xxxzz_xxxxxy_0[i] = ta_xxx_xxxxxy_0[i] * fe_0 - ta_xxx_xxxxxy_1[i] * fe_0 + ta_xxxz_xxxxxy_0[i] * pa_z[i] - ta_xxxz_xxxxxy_1[i] * pc_z[i];

        ta_xxxzz_xxxxxz_0[i] = 2.0 * ta_xzz_xxxxxz_0[i] * fe_0 - 2.0 * ta_xzz_xxxxxz_1[i] * fe_0 + 5.0 * ta_xxzz_xxxxz_0[i] * fe_0 -
                               5.0 * ta_xxzz_xxxxz_1[i] * fe_0 + ta_xxzz_xxxxxz_0[i] * pa_x[i] - ta_xxzz_xxxxxz_1[i] * pc_x[i];

        ta_xxxzz_xxxxyy_0[i] = ta_xxx_xxxxyy_0[i] * fe_0 - ta_xxx_xxxxyy_1[i] * fe_0 + ta_xxxz_xxxxyy_0[i] * pa_z[i] - ta_xxxz_xxxxyy_1[i] * pc_z[i];

        ta_xxxzz_xxxxyz_0[i] = 2.0 * ta_xzz_xxxxyz_0[i] * fe_0 - 2.0 * ta_xzz_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxzz_xxxyz_0[i] * fe_0 -
                               4.0 * ta_xxzz_xxxyz_1[i] * fe_0 + ta_xxzz_xxxxyz_0[i] * pa_x[i] - ta_xxzz_xxxxyz_1[i] * pc_x[i];

        ta_xxxzz_xxxxzz_0[i] = 2.0 * ta_xzz_xxxxzz_0[i] * fe_0 - 2.0 * ta_xzz_xxxxzz_1[i] * fe_0 + 4.0 * ta_xxzz_xxxzz_0[i] * fe_0 -
                               4.0 * ta_xxzz_xxxzz_1[i] * fe_0 + ta_xxzz_xxxxzz_0[i] * pa_x[i] - ta_xxzz_xxxxzz_1[i] * pc_x[i];

        ta_xxxzz_xxxyyy_0[i] = ta_xxx_xxxyyy_0[i] * fe_0 - ta_xxx_xxxyyy_1[i] * fe_0 + ta_xxxz_xxxyyy_0[i] * pa_z[i] - ta_xxxz_xxxyyy_1[i] * pc_z[i];

        ta_xxxzz_xxxyyz_0[i] = 2.0 * ta_xzz_xxxyyz_0[i] * fe_0 - 2.0 * ta_xzz_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxzz_xxyyz_0[i] * fe_0 -
                               3.0 * ta_xxzz_xxyyz_1[i] * fe_0 + ta_xxzz_xxxyyz_0[i] * pa_x[i] - ta_xxzz_xxxyyz_1[i] * pc_x[i];

        ta_xxxzz_xxxyzz_0[i] = 2.0 * ta_xzz_xxxyzz_0[i] * fe_0 - 2.0 * ta_xzz_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxzz_xxyzz_0[i] * fe_0 -
                               3.0 * ta_xxzz_xxyzz_1[i] * fe_0 + ta_xxzz_xxxyzz_0[i] * pa_x[i] - ta_xxzz_xxxyzz_1[i] * pc_x[i];

        ta_xxxzz_xxxzzz_0[i] = 2.0 * ta_xzz_xxxzzz_0[i] * fe_0 - 2.0 * ta_xzz_xxxzzz_1[i] * fe_0 + 3.0 * ta_xxzz_xxzzz_0[i] * fe_0 -
                               3.0 * ta_xxzz_xxzzz_1[i] * fe_0 + ta_xxzz_xxxzzz_0[i] * pa_x[i] - ta_xxzz_xxxzzz_1[i] * pc_x[i];

        ta_xxxzz_xxyyyy_0[i] = ta_xxx_xxyyyy_0[i] * fe_0 - ta_xxx_xxyyyy_1[i] * fe_0 + ta_xxxz_xxyyyy_0[i] * pa_z[i] - ta_xxxz_xxyyyy_1[i] * pc_z[i];

        ta_xxxzz_xxyyyz_0[i] = 2.0 * ta_xzz_xxyyyz_0[i] * fe_0 - 2.0 * ta_xzz_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxzz_xyyyz_0[i] * fe_0 -
                               2.0 * ta_xxzz_xyyyz_1[i] * fe_0 + ta_xxzz_xxyyyz_0[i] * pa_x[i] - ta_xxzz_xxyyyz_1[i] * pc_x[i];

        ta_xxxzz_xxyyzz_0[i] = 2.0 * ta_xzz_xxyyzz_0[i] * fe_0 - 2.0 * ta_xzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxzz_xyyzz_0[i] * fe_0 -
                               2.0 * ta_xxzz_xyyzz_1[i] * fe_0 + ta_xxzz_xxyyzz_0[i] * pa_x[i] - ta_xxzz_xxyyzz_1[i] * pc_x[i];

        ta_xxxzz_xxyzzz_0[i] = 2.0 * ta_xzz_xxyzzz_0[i] * fe_0 - 2.0 * ta_xzz_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxzz_xyzzz_0[i] * fe_0 -
                               2.0 * ta_xxzz_xyzzz_1[i] * fe_0 + ta_xxzz_xxyzzz_0[i] * pa_x[i] - ta_xxzz_xxyzzz_1[i] * pc_x[i];

        ta_xxxzz_xxzzzz_0[i] = 2.0 * ta_xzz_xxzzzz_0[i] * fe_0 - 2.0 * ta_xzz_xxzzzz_1[i] * fe_0 + 2.0 * ta_xxzz_xzzzz_0[i] * fe_0 -
                               2.0 * ta_xxzz_xzzzz_1[i] * fe_0 + ta_xxzz_xxzzzz_0[i] * pa_x[i] - ta_xxzz_xxzzzz_1[i] * pc_x[i];

        ta_xxxzz_xyyyyy_0[i] = ta_xxx_xyyyyy_0[i] * fe_0 - ta_xxx_xyyyyy_1[i] * fe_0 + ta_xxxz_xyyyyy_0[i] * pa_z[i] - ta_xxxz_xyyyyy_1[i] * pc_z[i];

        ta_xxxzz_xyyyyz_0[i] = 2.0 * ta_xzz_xyyyyz_0[i] * fe_0 - 2.0 * ta_xzz_xyyyyz_1[i] * fe_0 + ta_xxzz_yyyyz_0[i] * fe_0 -
                               ta_xxzz_yyyyz_1[i] * fe_0 + ta_xxzz_xyyyyz_0[i] * pa_x[i] - ta_xxzz_xyyyyz_1[i] * pc_x[i];

        ta_xxxzz_xyyyzz_0[i] = 2.0 * ta_xzz_xyyyzz_0[i] * fe_0 - 2.0 * ta_xzz_xyyyzz_1[i] * fe_0 + ta_xxzz_yyyzz_0[i] * fe_0 -
                               ta_xxzz_yyyzz_1[i] * fe_0 + ta_xxzz_xyyyzz_0[i] * pa_x[i] - ta_xxzz_xyyyzz_1[i] * pc_x[i];

        ta_xxxzz_xyyzzz_0[i] = 2.0 * ta_xzz_xyyzzz_0[i] * fe_0 - 2.0 * ta_xzz_xyyzzz_1[i] * fe_0 + ta_xxzz_yyzzz_0[i] * fe_0 -
                               ta_xxzz_yyzzz_1[i] * fe_0 + ta_xxzz_xyyzzz_0[i] * pa_x[i] - ta_xxzz_xyyzzz_1[i] * pc_x[i];

        ta_xxxzz_xyzzzz_0[i] = 2.0 * ta_xzz_xyzzzz_0[i] * fe_0 - 2.0 * ta_xzz_xyzzzz_1[i] * fe_0 + ta_xxzz_yzzzz_0[i] * fe_0 -
                               ta_xxzz_yzzzz_1[i] * fe_0 + ta_xxzz_xyzzzz_0[i] * pa_x[i] - ta_xxzz_xyzzzz_1[i] * pc_x[i];

        ta_xxxzz_xzzzzz_0[i] = 2.0 * ta_xzz_xzzzzz_0[i] * fe_0 - 2.0 * ta_xzz_xzzzzz_1[i] * fe_0 + ta_xxzz_zzzzz_0[i] * fe_0 -
                               ta_xxzz_zzzzz_1[i] * fe_0 + ta_xxzz_xzzzzz_0[i] * pa_x[i] - ta_xxzz_xzzzzz_1[i] * pc_x[i];

        ta_xxxzz_yyyyyy_0[i] =
            2.0 * ta_xzz_yyyyyy_0[i] * fe_0 - 2.0 * ta_xzz_yyyyyy_1[i] * fe_0 + ta_xxzz_yyyyyy_0[i] * pa_x[i] - ta_xxzz_yyyyyy_1[i] * pc_x[i];

        ta_xxxzz_yyyyyz_0[i] =
            2.0 * ta_xzz_yyyyyz_0[i] * fe_0 - 2.0 * ta_xzz_yyyyyz_1[i] * fe_0 + ta_xxzz_yyyyyz_0[i] * pa_x[i] - ta_xxzz_yyyyyz_1[i] * pc_x[i];

        ta_xxxzz_yyyyzz_0[i] =
            2.0 * ta_xzz_yyyyzz_0[i] * fe_0 - 2.0 * ta_xzz_yyyyzz_1[i] * fe_0 + ta_xxzz_yyyyzz_0[i] * pa_x[i] - ta_xxzz_yyyyzz_1[i] * pc_x[i];

        ta_xxxzz_yyyzzz_0[i] =
            2.0 * ta_xzz_yyyzzz_0[i] * fe_0 - 2.0 * ta_xzz_yyyzzz_1[i] * fe_0 + ta_xxzz_yyyzzz_0[i] * pa_x[i] - ta_xxzz_yyyzzz_1[i] * pc_x[i];

        ta_xxxzz_yyzzzz_0[i] =
            2.0 * ta_xzz_yyzzzz_0[i] * fe_0 - 2.0 * ta_xzz_yyzzzz_1[i] * fe_0 + ta_xxzz_yyzzzz_0[i] * pa_x[i] - ta_xxzz_yyzzzz_1[i] * pc_x[i];

        ta_xxxzz_yzzzzz_0[i] =
            2.0 * ta_xzz_yzzzzz_0[i] * fe_0 - 2.0 * ta_xzz_yzzzzz_1[i] * fe_0 + ta_xxzz_yzzzzz_0[i] * pa_x[i] - ta_xxzz_yzzzzz_1[i] * pc_x[i];

        ta_xxxzz_zzzzzz_0[i] =
            2.0 * ta_xzz_zzzzzz_0[i] * fe_0 - 2.0 * ta_xzz_zzzzzz_1[i] * fe_0 + ta_xxzz_zzzzzz_0[i] * pa_x[i] - ta_xxzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 168-196 components of targeted buffer : HI

    auto ta_xxyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 168);

    auto ta_xxyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 169);

    auto ta_xxyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 170);

    auto ta_xxyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 171);

    auto ta_xxyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 172);

    auto ta_xxyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 173);

    auto ta_xxyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 174);

    auto ta_xxyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 175);

    auto ta_xxyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 176);

    auto ta_xxyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 177);

    auto ta_xxyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 178);

    auto ta_xxyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 179);

    auto ta_xxyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 180);

    auto ta_xxyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 181);

    auto ta_xxyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 182);

    auto ta_xxyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 183);

    auto ta_xxyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 184);

    auto ta_xxyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 185);

    auto ta_xxyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 186);

    auto ta_xxyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 187);

    auto ta_xxyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 188);

    auto ta_xxyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 189);

    auto ta_xxyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 190);

    auto ta_xxyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 191);

    auto ta_xxyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 192);

    auto ta_xxyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 193);

    auto ta_xxyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 194);

    auto ta_xxyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 195);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxy_xxxxxx_0,   \
                             ta_xxy_xxxxxx_1,   \
                             ta_xxy_xxxxxz_0,   \
                             ta_xxy_xxxxxz_1,   \
                             ta_xxy_xxxxzz_0,   \
                             ta_xxy_xxxxzz_1,   \
                             ta_xxy_xxxzzz_0,   \
                             ta_xxy_xxxzzz_1,   \
                             ta_xxy_xxzzzz_0,   \
                             ta_xxy_xxzzzz_1,   \
                             ta_xxy_xzzzzz_0,   \
                             ta_xxy_xzzzzz_1,   \
                             ta_xxyy_xxxxxx_0,  \
                             ta_xxyy_xxxxxx_1,  \
                             ta_xxyy_xxxxxz_0,  \
                             ta_xxyy_xxxxxz_1,  \
                             ta_xxyy_xxxxzz_0,  \
                             ta_xxyy_xxxxzz_1,  \
                             ta_xxyy_xxxzzz_0,  \
                             ta_xxyy_xxxzzz_1,  \
                             ta_xxyy_xxzzzz_0,  \
                             ta_xxyy_xxzzzz_1,  \
                             ta_xxyy_xzzzzz_0,  \
                             ta_xxyy_xzzzzz_1,  \
                             ta_xxyyy_xxxxxx_0, \
                             ta_xxyyy_xxxxxy_0, \
                             ta_xxyyy_xxxxxz_0, \
                             ta_xxyyy_xxxxyy_0, \
                             ta_xxyyy_xxxxyz_0, \
                             ta_xxyyy_xxxxzz_0, \
                             ta_xxyyy_xxxyyy_0, \
                             ta_xxyyy_xxxyyz_0, \
                             ta_xxyyy_xxxyzz_0, \
                             ta_xxyyy_xxxzzz_0, \
                             ta_xxyyy_xxyyyy_0, \
                             ta_xxyyy_xxyyyz_0, \
                             ta_xxyyy_xxyyzz_0, \
                             ta_xxyyy_xxyzzz_0, \
                             ta_xxyyy_xxzzzz_0, \
                             ta_xxyyy_xyyyyy_0, \
                             ta_xxyyy_xyyyyz_0, \
                             ta_xxyyy_xyyyzz_0, \
                             ta_xxyyy_xyyzzz_0, \
                             ta_xxyyy_xyzzzz_0, \
                             ta_xxyyy_xzzzzz_0, \
                             ta_xxyyy_yyyyyy_0, \
                             ta_xxyyy_yyyyyz_0, \
                             ta_xxyyy_yyyyzz_0, \
                             ta_xxyyy_yyyzzz_0, \
                             ta_xxyyy_yyzzzz_0, \
                             ta_xxyyy_yzzzzz_0, \
                             ta_xxyyy_zzzzzz_0, \
                             ta_xyyy_xxxxxy_0,  \
                             ta_xyyy_xxxxxy_1,  \
                             ta_xyyy_xxxxy_0,   \
                             ta_xyyy_xxxxy_1,   \
                             ta_xyyy_xxxxyy_0,  \
                             ta_xyyy_xxxxyy_1,  \
                             ta_xyyy_xxxxyz_0,  \
                             ta_xyyy_xxxxyz_1,  \
                             ta_xyyy_xxxyy_0,   \
                             ta_xyyy_xxxyy_1,   \
                             ta_xyyy_xxxyyy_0,  \
                             ta_xyyy_xxxyyy_1,  \
                             ta_xyyy_xxxyyz_0,  \
                             ta_xyyy_xxxyyz_1,  \
                             ta_xyyy_xxxyz_0,   \
                             ta_xyyy_xxxyz_1,   \
                             ta_xyyy_xxxyzz_0,  \
                             ta_xyyy_xxxyzz_1,  \
                             ta_xyyy_xxyyy_0,   \
                             ta_xyyy_xxyyy_1,   \
                             ta_xyyy_xxyyyy_0,  \
                             ta_xyyy_xxyyyy_1,  \
                             ta_xyyy_xxyyyz_0,  \
                             ta_xyyy_xxyyyz_1,  \
                             ta_xyyy_xxyyz_0,   \
                             ta_xyyy_xxyyz_1,   \
                             ta_xyyy_xxyyzz_0,  \
                             ta_xyyy_xxyyzz_1,  \
                             ta_xyyy_xxyzz_0,   \
                             ta_xyyy_xxyzz_1,   \
                             ta_xyyy_xxyzzz_0,  \
                             ta_xyyy_xxyzzz_1,  \
                             ta_xyyy_xyyyy_0,   \
                             ta_xyyy_xyyyy_1,   \
                             ta_xyyy_xyyyyy_0,  \
                             ta_xyyy_xyyyyy_1,  \
                             ta_xyyy_xyyyyz_0,  \
                             ta_xyyy_xyyyyz_1,  \
                             ta_xyyy_xyyyz_0,   \
                             ta_xyyy_xyyyz_1,   \
                             ta_xyyy_xyyyzz_0,  \
                             ta_xyyy_xyyyzz_1,  \
                             ta_xyyy_xyyzz_0,   \
                             ta_xyyy_xyyzz_1,   \
                             ta_xyyy_xyyzzz_0,  \
                             ta_xyyy_xyyzzz_1,  \
                             ta_xyyy_xyzzz_0,   \
                             ta_xyyy_xyzzz_1,   \
                             ta_xyyy_xyzzzz_0,  \
                             ta_xyyy_xyzzzz_1,  \
                             ta_xyyy_yyyyy_0,   \
                             ta_xyyy_yyyyy_1,   \
                             ta_xyyy_yyyyyy_0,  \
                             ta_xyyy_yyyyyy_1,  \
                             ta_xyyy_yyyyyz_0,  \
                             ta_xyyy_yyyyyz_1,  \
                             ta_xyyy_yyyyz_0,   \
                             ta_xyyy_yyyyz_1,   \
                             ta_xyyy_yyyyzz_0,  \
                             ta_xyyy_yyyyzz_1,  \
                             ta_xyyy_yyyzz_0,   \
                             ta_xyyy_yyyzz_1,   \
                             ta_xyyy_yyyzzz_0,  \
                             ta_xyyy_yyyzzz_1,  \
                             ta_xyyy_yyzzz_0,   \
                             ta_xyyy_yyzzz_1,   \
                             ta_xyyy_yyzzzz_0,  \
                             ta_xyyy_yyzzzz_1,  \
                             ta_xyyy_yzzzz_0,   \
                             ta_xyyy_yzzzz_1,   \
                             ta_xyyy_yzzzzz_0,  \
                             ta_xyyy_yzzzzz_1,  \
                             ta_xyyy_zzzzzz_0,  \
                             ta_xyyy_zzzzzz_1,  \
                             ta_yyy_xxxxxy_0,   \
                             ta_yyy_xxxxxy_1,   \
                             ta_yyy_xxxxyy_0,   \
                             ta_yyy_xxxxyy_1,   \
                             ta_yyy_xxxxyz_0,   \
                             ta_yyy_xxxxyz_1,   \
                             ta_yyy_xxxyyy_0,   \
                             ta_yyy_xxxyyy_1,   \
                             ta_yyy_xxxyyz_0,   \
                             ta_yyy_xxxyyz_1,   \
                             ta_yyy_xxxyzz_0,   \
                             ta_yyy_xxxyzz_1,   \
                             ta_yyy_xxyyyy_0,   \
                             ta_yyy_xxyyyy_1,   \
                             ta_yyy_xxyyyz_0,   \
                             ta_yyy_xxyyyz_1,   \
                             ta_yyy_xxyyzz_0,   \
                             ta_yyy_xxyyzz_1,   \
                             ta_yyy_xxyzzz_0,   \
                             ta_yyy_xxyzzz_1,   \
                             ta_yyy_xyyyyy_0,   \
                             ta_yyy_xyyyyy_1,   \
                             ta_yyy_xyyyyz_0,   \
                             ta_yyy_xyyyyz_1,   \
                             ta_yyy_xyyyzz_0,   \
                             ta_yyy_xyyyzz_1,   \
                             ta_yyy_xyyzzz_0,   \
                             ta_yyy_xyyzzz_1,   \
                             ta_yyy_xyzzzz_0,   \
                             ta_yyy_xyzzzz_1,   \
                             ta_yyy_yyyyyy_0,   \
                             ta_yyy_yyyyyy_1,   \
                             ta_yyy_yyyyyz_0,   \
                             ta_yyy_yyyyyz_1,   \
                             ta_yyy_yyyyzz_0,   \
                             ta_yyy_yyyyzz_1,   \
                             ta_yyy_yyyzzz_0,   \
                             ta_yyy_yyyzzz_1,   \
                             ta_yyy_yyzzzz_0,   \
                             ta_yyy_yyzzzz_1,   \
                             ta_yyy_yzzzzz_0,   \
                             ta_yyy_yzzzzz_1,   \
                             ta_yyy_zzzzzz_0,   \
                             ta_yyy_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyy_xxxxxx_0[i] =
            2.0 * ta_xxy_xxxxxx_0[i] * fe_0 - 2.0 * ta_xxy_xxxxxx_1[i] * fe_0 + ta_xxyy_xxxxxx_0[i] * pa_y[i] - ta_xxyy_xxxxxx_1[i] * pc_y[i];

        ta_xxyyy_xxxxxy_0[i] = ta_yyy_xxxxxy_0[i] * fe_0 - ta_yyy_xxxxxy_1[i] * fe_0 + 5.0 * ta_xyyy_xxxxy_0[i] * fe_0 -
                               5.0 * ta_xyyy_xxxxy_1[i] * fe_0 + ta_xyyy_xxxxxy_0[i] * pa_x[i] - ta_xyyy_xxxxxy_1[i] * pc_x[i];

        ta_xxyyy_xxxxxz_0[i] =
            2.0 * ta_xxy_xxxxxz_0[i] * fe_0 - 2.0 * ta_xxy_xxxxxz_1[i] * fe_0 + ta_xxyy_xxxxxz_0[i] * pa_y[i] - ta_xxyy_xxxxxz_1[i] * pc_y[i];

        ta_xxyyy_xxxxyy_0[i] = ta_yyy_xxxxyy_0[i] * fe_0 - ta_yyy_xxxxyy_1[i] * fe_0 + 4.0 * ta_xyyy_xxxyy_0[i] * fe_0 -
                               4.0 * ta_xyyy_xxxyy_1[i] * fe_0 + ta_xyyy_xxxxyy_0[i] * pa_x[i] - ta_xyyy_xxxxyy_1[i] * pc_x[i];

        ta_xxyyy_xxxxyz_0[i] = ta_yyy_xxxxyz_0[i] * fe_0 - ta_yyy_xxxxyz_1[i] * fe_0 + 4.0 * ta_xyyy_xxxyz_0[i] * fe_0 -
                               4.0 * ta_xyyy_xxxyz_1[i] * fe_0 + ta_xyyy_xxxxyz_0[i] * pa_x[i] - ta_xyyy_xxxxyz_1[i] * pc_x[i];

        ta_xxyyy_xxxxzz_0[i] =
            2.0 * ta_xxy_xxxxzz_0[i] * fe_0 - 2.0 * ta_xxy_xxxxzz_1[i] * fe_0 + ta_xxyy_xxxxzz_0[i] * pa_y[i] - ta_xxyy_xxxxzz_1[i] * pc_y[i];

        ta_xxyyy_xxxyyy_0[i] = ta_yyy_xxxyyy_0[i] * fe_0 - ta_yyy_xxxyyy_1[i] * fe_0 + 3.0 * ta_xyyy_xxyyy_0[i] * fe_0 -
                               3.0 * ta_xyyy_xxyyy_1[i] * fe_0 + ta_xyyy_xxxyyy_0[i] * pa_x[i] - ta_xyyy_xxxyyy_1[i] * pc_x[i];

        ta_xxyyy_xxxyyz_0[i] = ta_yyy_xxxyyz_0[i] * fe_0 - ta_yyy_xxxyyz_1[i] * fe_0 + 3.0 * ta_xyyy_xxyyz_0[i] * fe_0 -
                               3.0 * ta_xyyy_xxyyz_1[i] * fe_0 + ta_xyyy_xxxyyz_0[i] * pa_x[i] - ta_xyyy_xxxyyz_1[i] * pc_x[i];

        ta_xxyyy_xxxyzz_0[i] = ta_yyy_xxxyzz_0[i] * fe_0 - ta_yyy_xxxyzz_1[i] * fe_0 + 3.0 * ta_xyyy_xxyzz_0[i] * fe_0 -
                               3.0 * ta_xyyy_xxyzz_1[i] * fe_0 + ta_xyyy_xxxyzz_0[i] * pa_x[i] - ta_xyyy_xxxyzz_1[i] * pc_x[i];

        ta_xxyyy_xxxzzz_0[i] =
            2.0 * ta_xxy_xxxzzz_0[i] * fe_0 - 2.0 * ta_xxy_xxxzzz_1[i] * fe_0 + ta_xxyy_xxxzzz_0[i] * pa_y[i] - ta_xxyy_xxxzzz_1[i] * pc_y[i];

        ta_xxyyy_xxyyyy_0[i] = ta_yyy_xxyyyy_0[i] * fe_0 - ta_yyy_xxyyyy_1[i] * fe_0 + 2.0 * ta_xyyy_xyyyy_0[i] * fe_0 -
                               2.0 * ta_xyyy_xyyyy_1[i] * fe_0 + ta_xyyy_xxyyyy_0[i] * pa_x[i] - ta_xyyy_xxyyyy_1[i] * pc_x[i];

        ta_xxyyy_xxyyyz_0[i] = ta_yyy_xxyyyz_0[i] * fe_0 - ta_yyy_xxyyyz_1[i] * fe_0 + 2.0 * ta_xyyy_xyyyz_0[i] * fe_0 -
                               2.0 * ta_xyyy_xyyyz_1[i] * fe_0 + ta_xyyy_xxyyyz_0[i] * pa_x[i] - ta_xyyy_xxyyyz_1[i] * pc_x[i];

        ta_xxyyy_xxyyzz_0[i] = ta_yyy_xxyyzz_0[i] * fe_0 - ta_yyy_xxyyzz_1[i] * fe_0 + 2.0 * ta_xyyy_xyyzz_0[i] * fe_0 -
                               2.0 * ta_xyyy_xyyzz_1[i] * fe_0 + ta_xyyy_xxyyzz_0[i] * pa_x[i] - ta_xyyy_xxyyzz_1[i] * pc_x[i];

        ta_xxyyy_xxyzzz_0[i] = ta_yyy_xxyzzz_0[i] * fe_0 - ta_yyy_xxyzzz_1[i] * fe_0 + 2.0 * ta_xyyy_xyzzz_0[i] * fe_0 -
                               2.0 * ta_xyyy_xyzzz_1[i] * fe_0 + ta_xyyy_xxyzzz_0[i] * pa_x[i] - ta_xyyy_xxyzzz_1[i] * pc_x[i];

        ta_xxyyy_xxzzzz_0[i] =
            2.0 * ta_xxy_xxzzzz_0[i] * fe_0 - 2.0 * ta_xxy_xxzzzz_1[i] * fe_0 + ta_xxyy_xxzzzz_0[i] * pa_y[i] - ta_xxyy_xxzzzz_1[i] * pc_y[i];

        ta_xxyyy_xyyyyy_0[i] = ta_yyy_xyyyyy_0[i] * fe_0 - ta_yyy_xyyyyy_1[i] * fe_0 + ta_xyyy_yyyyy_0[i] * fe_0 - ta_xyyy_yyyyy_1[i] * fe_0 +
                               ta_xyyy_xyyyyy_0[i] * pa_x[i] - ta_xyyy_xyyyyy_1[i] * pc_x[i];

        ta_xxyyy_xyyyyz_0[i] = ta_yyy_xyyyyz_0[i] * fe_0 - ta_yyy_xyyyyz_1[i] * fe_0 + ta_xyyy_yyyyz_0[i] * fe_0 - ta_xyyy_yyyyz_1[i] * fe_0 +
                               ta_xyyy_xyyyyz_0[i] * pa_x[i] - ta_xyyy_xyyyyz_1[i] * pc_x[i];

        ta_xxyyy_xyyyzz_0[i] = ta_yyy_xyyyzz_0[i] * fe_0 - ta_yyy_xyyyzz_1[i] * fe_0 + ta_xyyy_yyyzz_0[i] * fe_0 - ta_xyyy_yyyzz_1[i] * fe_0 +
                               ta_xyyy_xyyyzz_0[i] * pa_x[i] - ta_xyyy_xyyyzz_1[i] * pc_x[i];

        ta_xxyyy_xyyzzz_0[i] = ta_yyy_xyyzzz_0[i] * fe_0 - ta_yyy_xyyzzz_1[i] * fe_0 + ta_xyyy_yyzzz_0[i] * fe_0 - ta_xyyy_yyzzz_1[i] * fe_0 +
                               ta_xyyy_xyyzzz_0[i] * pa_x[i] - ta_xyyy_xyyzzz_1[i] * pc_x[i];

        ta_xxyyy_xyzzzz_0[i] = ta_yyy_xyzzzz_0[i] * fe_0 - ta_yyy_xyzzzz_1[i] * fe_0 + ta_xyyy_yzzzz_0[i] * fe_0 - ta_xyyy_yzzzz_1[i] * fe_0 +
                               ta_xyyy_xyzzzz_0[i] * pa_x[i] - ta_xyyy_xyzzzz_1[i] * pc_x[i];

        ta_xxyyy_xzzzzz_0[i] =
            2.0 * ta_xxy_xzzzzz_0[i] * fe_0 - 2.0 * ta_xxy_xzzzzz_1[i] * fe_0 + ta_xxyy_xzzzzz_0[i] * pa_y[i] - ta_xxyy_xzzzzz_1[i] * pc_y[i];

        ta_xxyyy_yyyyyy_0[i] = ta_yyy_yyyyyy_0[i] * fe_0 - ta_yyy_yyyyyy_1[i] * fe_0 + ta_xyyy_yyyyyy_0[i] * pa_x[i] - ta_xyyy_yyyyyy_1[i] * pc_x[i];

        ta_xxyyy_yyyyyz_0[i] = ta_yyy_yyyyyz_0[i] * fe_0 - ta_yyy_yyyyyz_1[i] * fe_0 + ta_xyyy_yyyyyz_0[i] * pa_x[i] - ta_xyyy_yyyyyz_1[i] * pc_x[i];

        ta_xxyyy_yyyyzz_0[i] = ta_yyy_yyyyzz_0[i] * fe_0 - ta_yyy_yyyyzz_1[i] * fe_0 + ta_xyyy_yyyyzz_0[i] * pa_x[i] - ta_xyyy_yyyyzz_1[i] * pc_x[i];

        ta_xxyyy_yyyzzz_0[i] = ta_yyy_yyyzzz_0[i] * fe_0 - ta_yyy_yyyzzz_1[i] * fe_0 + ta_xyyy_yyyzzz_0[i] * pa_x[i] - ta_xyyy_yyyzzz_1[i] * pc_x[i];

        ta_xxyyy_yyzzzz_0[i] = ta_yyy_yyzzzz_0[i] * fe_0 - ta_yyy_yyzzzz_1[i] * fe_0 + ta_xyyy_yyzzzz_0[i] * pa_x[i] - ta_xyyy_yyzzzz_1[i] * pc_x[i];

        ta_xxyyy_yzzzzz_0[i] = ta_yyy_yzzzzz_0[i] * fe_0 - ta_yyy_yzzzzz_1[i] * fe_0 + ta_xyyy_yzzzzz_0[i] * pa_x[i] - ta_xyyy_yzzzzz_1[i] * pc_x[i];

        ta_xxyyy_zzzzzz_0[i] = ta_yyy_zzzzzz_0[i] * fe_0 - ta_yyy_zzzzzz_1[i] * fe_0 + ta_xyyy_zzzzzz_0[i] * pa_x[i] - ta_xyyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 196-224 components of targeted buffer : HI

    auto ta_xxyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 196);

    auto ta_xxyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 197);

    auto ta_xxyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 198);

    auto ta_xxyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 199);

    auto ta_xxyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 200);

    auto ta_xxyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 201);

    auto ta_xxyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 202);

    auto ta_xxyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 203);

    auto ta_xxyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 204);

    auto ta_xxyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 205);

    auto ta_xxyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 206);

    auto ta_xxyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 207);

    auto ta_xxyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 208);

    auto ta_xxyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 209);

    auto ta_xxyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 210);

    auto ta_xxyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 211);

    auto ta_xxyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 212);

    auto ta_xxyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 213);

    auto ta_xxyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 214);

    auto ta_xxyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 215);

    auto ta_xxyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 216);

    auto ta_xxyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 217);

    auto ta_xxyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 218);

    auto ta_xxyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 219);

    auto ta_xxyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 220);

    auto ta_xxyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 221);

    auto ta_xxyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 222);

    auto ta_xxyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 223);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta_xxyy_xxxxxx_0,  \
                             ta_xxyy_xxxxxx_1,  \
                             ta_xxyy_xxxxxy_0,  \
                             ta_xxyy_xxxxxy_1,  \
                             ta_xxyy_xxxxy_0,   \
                             ta_xxyy_xxxxy_1,   \
                             ta_xxyy_xxxxyy_0,  \
                             ta_xxyy_xxxxyy_1,  \
                             ta_xxyy_xxxxyz_0,  \
                             ta_xxyy_xxxxyz_1,  \
                             ta_xxyy_xxxyy_0,   \
                             ta_xxyy_xxxyy_1,   \
                             ta_xxyy_xxxyyy_0,  \
                             ta_xxyy_xxxyyy_1,  \
                             ta_xxyy_xxxyyz_0,  \
                             ta_xxyy_xxxyyz_1,  \
                             ta_xxyy_xxxyz_0,   \
                             ta_xxyy_xxxyz_1,   \
                             ta_xxyy_xxxyzz_0,  \
                             ta_xxyy_xxxyzz_1,  \
                             ta_xxyy_xxyyy_0,   \
                             ta_xxyy_xxyyy_1,   \
                             ta_xxyy_xxyyyy_0,  \
                             ta_xxyy_xxyyyy_1,  \
                             ta_xxyy_xxyyyz_0,  \
                             ta_xxyy_xxyyyz_1,  \
                             ta_xxyy_xxyyz_0,   \
                             ta_xxyy_xxyyz_1,   \
                             ta_xxyy_xxyyzz_0,  \
                             ta_xxyy_xxyyzz_1,  \
                             ta_xxyy_xxyzz_0,   \
                             ta_xxyy_xxyzz_1,   \
                             ta_xxyy_xxyzzz_0,  \
                             ta_xxyy_xxyzzz_1,  \
                             ta_xxyy_xyyyy_0,   \
                             ta_xxyy_xyyyy_1,   \
                             ta_xxyy_xyyyyy_0,  \
                             ta_xxyy_xyyyyy_1,  \
                             ta_xxyy_xyyyyz_0,  \
                             ta_xxyy_xyyyyz_1,  \
                             ta_xxyy_xyyyz_0,   \
                             ta_xxyy_xyyyz_1,   \
                             ta_xxyy_xyyyzz_0,  \
                             ta_xxyy_xyyyzz_1,  \
                             ta_xxyy_xyyzz_0,   \
                             ta_xxyy_xyyzz_1,   \
                             ta_xxyy_xyyzzz_0,  \
                             ta_xxyy_xyyzzz_1,  \
                             ta_xxyy_xyzzz_0,   \
                             ta_xxyy_xyzzz_1,   \
                             ta_xxyy_xyzzzz_0,  \
                             ta_xxyy_xyzzzz_1,  \
                             ta_xxyy_yyyyyy_0,  \
                             ta_xxyy_yyyyyy_1,  \
                             ta_xxyyz_xxxxxx_0, \
                             ta_xxyyz_xxxxxy_0, \
                             ta_xxyyz_xxxxxz_0, \
                             ta_xxyyz_xxxxyy_0, \
                             ta_xxyyz_xxxxyz_0, \
                             ta_xxyyz_xxxxzz_0, \
                             ta_xxyyz_xxxyyy_0, \
                             ta_xxyyz_xxxyyz_0, \
                             ta_xxyyz_xxxyzz_0, \
                             ta_xxyyz_xxxzzz_0, \
                             ta_xxyyz_xxyyyy_0, \
                             ta_xxyyz_xxyyyz_0, \
                             ta_xxyyz_xxyyzz_0, \
                             ta_xxyyz_xxyzzz_0, \
                             ta_xxyyz_xxzzzz_0, \
                             ta_xxyyz_xyyyyy_0, \
                             ta_xxyyz_xyyyyz_0, \
                             ta_xxyyz_xyyyzz_0, \
                             ta_xxyyz_xyyzzz_0, \
                             ta_xxyyz_xyzzzz_0, \
                             ta_xxyyz_xzzzzz_0, \
                             ta_xxyyz_yyyyyy_0, \
                             ta_xxyyz_yyyyyz_0, \
                             ta_xxyyz_yyyyzz_0, \
                             ta_xxyyz_yyyzzz_0, \
                             ta_xxyyz_yyzzzz_0, \
                             ta_xxyyz_yzzzzz_0, \
                             ta_xxyyz_zzzzzz_0, \
                             ta_xxyz_xxxxxz_0,  \
                             ta_xxyz_xxxxxz_1,  \
                             ta_xxyz_xxxxzz_0,  \
                             ta_xxyz_xxxxzz_1,  \
                             ta_xxyz_xxxzzz_0,  \
                             ta_xxyz_xxxzzz_1,  \
                             ta_xxyz_xxzzzz_0,  \
                             ta_xxyz_xxzzzz_1,  \
                             ta_xxyz_xzzzzz_0,  \
                             ta_xxyz_xzzzzz_1,  \
                             ta_xxz_xxxxxz_0,   \
                             ta_xxz_xxxxxz_1,   \
                             ta_xxz_xxxxzz_0,   \
                             ta_xxz_xxxxzz_1,   \
                             ta_xxz_xxxzzz_0,   \
                             ta_xxz_xxxzzz_1,   \
                             ta_xxz_xxzzzz_0,   \
                             ta_xxz_xxzzzz_1,   \
                             ta_xxz_xzzzzz_0,   \
                             ta_xxz_xzzzzz_1,   \
                             ta_xyyz_yyyyyz_0,  \
                             ta_xyyz_yyyyyz_1,  \
                             ta_xyyz_yyyyzz_0,  \
                             ta_xyyz_yyyyzz_1,  \
                             ta_xyyz_yyyzzz_0,  \
                             ta_xyyz_yyyzzz_1,  \
                             ta_xyyz_yyzzzz_0,  \
                             ta_xyyz_yyzzzz_1,  \
                             ta_xyyz_yzzzzz_0,  \
                             ta_xyyz_yzzzzz_1,  \
                             ta_xyyz_zzzzzz_0,  \
                             ta_xyyz_zzzzzz_1,  \
                             ta_yyz_yyyyyz_0,   \
                             ta_yyz_yyyyyz_1,   \
                             ta_yyz_yyyyzz_0,   \
                             ta_yyz_yyyyzz_1,   \
                             ta_yyz_yyyzzz_0,   \
                             ta_yyz_yyyzzz_1,   \
                             ta_yyz_yyzzzz_0,   \
                             ta_yyz_yyzzzz_1,   \
                             ta_yyz_yzzzzz_0,   \
                             ta_yyz_yzzzzz_1,   \
                             ta_yyz_zzzzzz_0,   \
                             ta_yyz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyz_xxxxxx_0[i] = ta_xxyy_xxxxxx_0[i] * pa_z[i] - ta_xxyy_xxxxxx_1[i] * pc_z[i];

        ta_xxyyz_xxxxxy_0[i] = ta_xxyy_xxxxxy_0[i] * pa_z[i] - ta_xxyy_xxxxxy_1[i] * pc_z[i];

        ta_xxyyz_xxxxxz_0[i] = ta_xxz_xxxxxz_0[i] * fe_0 - ta_xxz_xxxxxz_1[i] * fe_0 + ta_xxyz_xxxxxz_0[i] * pa_y[i] - ta_xxyz_xxxxxz_1[i] * pc_y[i];

        ta_xxyyz_xxxxyy_0[i] = ta_xxyy_xxxxyy_0[i] * pa_z[i] - ta_xxyy_xxxxyy_1[i] * pc_z[i];

        ta_xxyyz_xxxxyz_0[i] = ta_xxyy_xxxxy_0[i] * fe_0 - ta_xxyy_xxxxy_1[i] * fe_0 + ta_xxyy_xxxxyz_0[i] * pa_z[i] - ta_xxyy_xxxxyz_1[i] * pc_z[i];

        ta_xxyyz_xxxxzz_0[i] = ta_xxz_xxxxzz_0[i] * fe_0 - ta_xxz_xxxxzz_1[i] * fe_0 + ta_xxyz_xxxxzz_0[i] * pa_y[i] - ta_xxyz_xxxxzz_1[i] * pc_y[i];

        ta_xxyyz_xxxyyy_0[i] = ta_xxyy_xxxyyy_0[i] * pa_z[i] - ta_xxyy_xxxyyy_1[i] * pc_z[i];

        ta_xxyyz_xxxyyz_0[i] = ta_xxyy_xxxyy_0[i] * fe_0 - ta_xxyy_xxxyy_1[i] * fe_0 + ta_xxyy_xxxyyz_0[i] * pa_z[i] - ta_xxyy_xxxyyz_1[i] * pc_z[i];

        ta_xxyyz_xxxyzz_0[i] =
            2.0 * ta_xxyy_xxxyz_0[i] * fe_0 - 2.0 * ta_xxyy_xxxyz_1[i] * fe_0 + ta_xxyy_xxxyzz_0[i] * pa_z[i] - ta_xxyy_xxxyzz_1[i] * pc_z[i];

        ta_xxyyz_xxxzzz_0[i] = ta_xxz_xxxzzz_0[i] * fe_0 - ta_xxz_xxxzzz_1[i] * fe_0 + ta_xxyz_xxxzzz_0[i] * pa_y[i] - ta_xxyz_xxxzzz_1[i] * pc_y[i];

        ta_xxyyz_xxyyyy_0[i] = ta_xxyy_xxyyyy_0[i] * pa_z[i] - ta_xxyy_xxyyyy_1[i] * pc_z[i];

        ta_xxyyz_xxyyyz_0[i] = ta_xxyy_xxyyy_0[i] * fe_0 - ta_xxyy_xxyyy_1[i] * fe_0 + ta_xxyy_xxyyyz_0[i] * pa_z[i] - ta_xxyy_xxyyyz_1[i] * pc_z[i];

        ta_xxyyz_xxyyzz_0[i] =
            2.0 * ta_xxyy_xxyyz_0[i] * fe_0 - 2.0 * ta_xxyy_xxyyz_1[i] * fe_0 + ta_xxyy_xxyyzz_0[i] * pa_z[i] - ta_xxyy_xxyyzz_1[i] * pc_z[i];

        ta_xxyyz_xxyzzz_0[i] =
            3.0 * ta_xxyy_xxyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xxyzz_1[i] * fe_0 + ta_xxyy_xxyzzz_0[i] * pa_z[i] - ta_xxyy_xxyzzz_1[i] * pc_z[i];

        ta_xxyyz_xxzzzz_0[i] = ta_xxz_xxzzzz_0[i] * fe_0 - ta_xxz_xxzzzz_1[i] * fe_0 + ta_xxyz_xxzzzz_0[i] * pa_y[i] - ta_xxyz_xxzzzz_1[i] * pc_y[i];

        ta_xxyyz_xyyyyy_0[i] = ta_xxyy_xyyyyy_0[i] * pa_z[i] - ta_xxyy_xyyyyy_1[i] * pc_z[i];

        ta_xxyyz_xyyyyz_0[i] = ta_xxyy_xyyyy_0[i] * fe_0 - ta_xxyy_xyyyy_1[i] * fe_0 + ta_xxyy_xyyyyz_0[i] * pa_z[i] - ta_xxyy_xyyyyz_1[i] * pc_z[i];

        ta_xxyyz_xyyyzz_0[i] =
            2.0 * ta_xxyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xxyy_xyyyz_1[i] * fe_0 + ta_xxyy_xyyyzz_0[i] * pa_z[i] - ta_xxyy_xyyyzz_1[i] * pc_z[i];

        ta_xxyyz_xyyzzz_0[i] =
            3.0 * ta_xxyy_xyyzz_0[i] * fe_0 - 3.0 * ta_xxyy_xyyzz_1[i] * fe_0 + ta_xxyy_xyyzzz_0[i] * pa_z[i] - ta_xxyy_xyyzzz_1[i] * pc_z[i];

        ta_xxyyz_xyzzzz_0[i] =
            4.0 * ta_xxyy_xyzzz_0[i] * fe_0 - 4.0 * ta_xxyy_xyzzz_1[i] * fe_0 + ta_xxyy_xyzzzz_0[i] * pa_z[i] - ta_xxyy_xyzzzz_1[i] * pc_z[i];

        ta_xxyyz_xzzzzz_0[i] = ta_xxz_xzzzzz_0[i] * fe_0 - ta_xxz_xzzzzz_1[i] * fe_0 + ta_xxyz_xzzzzz_0[i] * pa_y[i] - ta_xxyz_xzzzzz_1[i] * pc_y[i];

        ta_xxyyz_yyyyyy_0[i] = ta_xxyy_yyyyyy_0[i] * pa_z[i] - ta_xxyy_yyyyyy_1[i] * pc_z[i];

        ta_xxyyz_yyyyyz_0[i] = ta_yyz_yyyyyz_0[i] * fe_0 - ta_yyz_yyyyyz_1[i] * fe_0 + ta_xyyz_yyyyyz_0[i] * pa_x[i] - ta_xyyz_yyyyyz_1[i] * pc_x[i];

        ta_xxyyz_yyyyzz_0[i] = ta_yyz_yyyyzz_0[i] * fe_0 - ta_yyz_yyyyzz_1[i] * fe_0 + ta_xyyz_yyyyzz_0[i] * pa_x[i] - ta_xyyz_yyyyzz_1[i] * pc_x[i];

        ta_xxyyz_yyyzzz_0[i] = ta_yyz_yyyzzz_0[i] * fe_0 - ta_yyz_yyyzzz_1[i] * fe_0 + ta_xyyz_yyyzzz_0[i] * pa_x[i] - ta_xyyz_yyyzzz_1[i] * pc_x[i];

        ta_xxyyz_yyzzzz_0[i] = ta_yyz_yyzzzz_0[i] * fe_0 - ta_yyz_yyzzzz_1[i] * fe_0 + ta_xyyz_yyzzzz_0[i] * pa_x[i] - ta_xyyz_yyzzzz_1[i] * pc_x[i];

        ta_xxyyz_yzzzzz_0[i] = ta_yyz_yzzzzz_0[i] * fe_0 - ta_yyz_yzzzzz_1[i] * fe_0 + ta_xyyz_yzzzzz_0[i] * pa_x[i] - ta_xyyz_yzzzzz_1[i] * pc_x[i];

        ta_xxyyz_zzzzzz_0[i] = ta_yyz_zzzzzz_0[i] * fe_0 - ta_yyz_zzzzzz_1[i] * fe_0 + ta_xyyz_zzzzzz_0[i] * pa_x[i] - ta_xyyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 224-252 components of targeted buffer : HI

    auto ta_xxyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 224);

    auto ta_xxyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 225);

    auto ta_xxyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 226);

    auto ta_xxyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 227);

    auto ta_xxyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 228);

    auto ta_xxyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 229);

    auto ta_xxyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 230);

    auto ta_xxyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 231);

    auto ta_xxyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 232);

    auto ta_xxyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 233);

    auto ta_xxyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 234);

    auto ta_xxyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 235);

    auto ta_xxyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 236);

    auto ta_xxyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 237);

    auto ta_xxyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 238);

    auto ta_xxyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 239);

    auto ta_xxyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 240);

    auto ta_xxyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 241);

    auto ta_xxyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 242);

    auto ta_xxyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 243);

    auto ta_xxyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 244);

    auto ta_xxyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 245);

    auto ta_xxyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 246);

    auto ta_xxyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 247);

    auto ta_xxyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 248);

    auto ta_xxyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 249);

    auto ta_xxyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 250);

    auto ta_xxyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 251);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xxyzz_xxxxxx_0, \
                             ta_xxyzz_xxxxxy_0, \
                             ta_xxyzz_xxxxxz_0, \
                             ta_xxyzz_xxxxyy_0, \
                             ta_xxyzz_xxxxyz_0, \
                             ta_xxyzz_xxxxzz_0, \
                             ta_xxyzz_xxxyyy_0, \
                             ta_xxyzz_xxxyyz_0, \
                             ta_xxyzz_xxxyzz_0, \
                             ta_xxyzz_xxxzzz_0, \
                             ta_xxyzz_xxyyyy_0, \
                             ta_xxyzz_xxyyyz_0, \
                             ta_xxyzz_xxyyzz_0, \
                             ta_xxyzz_xxyzzz_0, \
                             ta_xxyzz_xxzzzz_0, \
                             ta_xxyzz_xyyyyy_0, \
                             ta_xxyzz_xyyyyz_0, \
                             ta_xxyzz_xyyyzz_0, \
                             ta_xxyzz_xyyzzz_0, \
                             ta_xxyzz_xyzzzz_0, \
                             ta_xxyzz_xzzzzz_0, \
                             ta_xxyzz_yyyyyy_0, \
                             ta_xxyzz_yyyyyz_0, \
                             ta_xxyzz_yyyyzz_0, \
                             ta_xxyzz_yyyzzz_0, \
                             ta_xxyzz_yyzzzz_0, \
                             ta_xxyzz_yzzzzz_0, \
                             ta_xxyzz_zzzzzz_0, \
                             ta_xxzz_xxxxx_0,   \
                             ta_xxzz_xxxxx_1,   \
                             ta_xxzz_xxxxxx_0,  \
                             ta_xxzz_xxxxxx_1,  \
                             ta_xxzz_xxxxxy_0,  \
                             ta_xxzz_xxxxxy_1,  \
                             ta_xxzz_xxxxxz_0,  \
                             ta_xxzz_xxxxxz_1,  \
                             ta_xxzz_xxxxy_0,   \
                             ta_xxzz_xxxxy_1,   \
                             ta_xxzz_xxxxyy_0,  \
                             ta_xxzz_xxxxyy_1,  \
                             ta_xxzz_xxxxyz_0,  \
                             ta_xxzz_xxxxyz_1,  \
                             ta_xxzz_xxxxz_0,   \
                             ta_xxzz_xxxxz_1,   \
                             ta_xxzz_xxxxzz_0,  \
                             ta_xxzz_xxxxzz_1,  \
                             ta_xxzz_xxxyy_0,   \
                             ta_xxzz_xxxyy_1,   \
                             ta_xxzz_xxxyyy_0,  \
                             ta_xxzz_xxxyyy_1,  \
                             ta_xxzz_xxxyyz_0,  \
                             ta_xxzz_xxxyyz_1,  \
                             ta_xxzz_xxxyz_0,   \
                             ta_xxzz_xxxyz_1,   \
                             ta_xxzz_xxxyzz_0,  \
                             ta_xxzz_xxxyzz_1,  \
                             ta_xxzz_xxxzz_0,   \
                             ta_xxzz_xxxzz_1,   \
                             ta_xxzz_xxxzzz_0,  \
                             ta_xxzz_xxxzzz_1,  \
                             ta_xxzz_xxyyy_0,   \
                             ta_xxzz_xxyyy_1,   \
                             ta_xxzz_xxyyyy_0,  \
                             ta_xxzz_xxyyyy_1,  \
                             ta_xxzz_xxyyyz_0,  \
                             ta_xxzz_xxyyyz_1,  \
                             ta_xxzz_xxyyz_0,   \
                             ta_xxzz_xxyyz_1,   \
                             ta_xxzz_xxyyzz_0,  \
                             ta_xxzz_xxyyzz_1,  \
                             ta_xxzz_xxyzz_0,   \
                             ta_xxzz_xxyzz_1,   \
                             ta_xxzz_xxyzzz_0,  \
                             ta_xxzz_xxyzzz_1,  \
                             ta_xxzz_xxzzz_0,   \
                             ta_xxzz_xxzzz_1,   \
                             ta_xxzz_xxzzzz_0,  \
                             ta_xxzz_xxzzzz_1,  \
                             ta_xxzz_xyyyy_0,   \
                             ta_xxzz_xyyyy_1,   \
                             ta_xxzz_xyyyyy_0,  \
                             ta_xxzz_xyyyyy_1,  \
                             ta_xxzz_xyyyyz_0,  \
                             ta_xxzz_xyyyyz_1,  \
                             ta_xxzz_xyyyz_0,   \
                             ta_xxzz_xyyyz_1,   \
                             ta_xxzz_xyyyzz_0,  \
                             ta_xxzz_xyyyzz_1,  \
                             ta_xxzz_xyyzz_0,   \
                             ta_xxzz_xyyzz_1,   \
                             ta_xxzz_xyyzzz_0,  \
                             ta_xxzz_xyyzzz_1,  \
                             ta_xxzz_xyzzz_0,   \
                             ta_xxzz_xyzzz_1,   \
                             ta_xxzz_xyzzzz_0,  \
                             ta_xxzz_xyzzzz_1,  \
                             ta_xxzz_xzzzz_0,   \
                             ta_xxzz_xzzzz_1,   \
                             ta_xxzz_xzzzzz_0,  \
                             ta_xxzz_xzzzzz_1,  \
                             ta_xxzz_zzzzzz_0,  \
                             ta_xxzz_zzzzzz_1,  \
                             ta_xyzz_yyyyyy_0,  \
                             ta_xyzz_yyyyyy_1,  \
                             ta_xyzz_yyyyyz_0,  \
                             ta_xyzz_yyyyyz_1,  \
                             ta_xyzz_yyyyzz_0,  \
                             ta_xyzz_yyyyzz_1,  \
                             ta_xyzz_yyyzzz_0,  \
                             ta_xyzz_yyyzzz_1,  \
                             ta_xyzz_yyzzzz_0,  \
                             ta_xyzz_yyzzzz_1,  \
                             ta_xyzz_yzzzzz_0,  \
                             ta_xyzz_yzzzzz_1,  \
                             ta_yzz_yyyyyy_0,   \
                             ta_yzz_yyyyyy_1,   \
                             ta_yzz_yyyyyz_0,   \
                             ta_yzz_yyyyyz_1,   \
                             ta_yzz_yyyyzz_0,   \
                             ta_yzz_yyyyzz_1,   \
                             ta_yzz_yyyzzz_0,   \
                             ta_yzz_yyyzzz_1,   \
                             ta_yzz_yyzzzz_0,   \
                             ta_yzz_yyzzzz_1,   \
                             ta_yzz_yzzzzz_0,   \
                             ta_yzz_yzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzz_xxxxxx_0[i] = ta_xxzz_xxxxxx_0[i] * pa_y[i] - ta_xxzz_xxxxxx_1[i] * pc_y[i];

        ta_xxyzz_xxxxxy_0[i] = ta_xxzz_xxxxx_0[i] * fe_0 - ta_xxzz_xxxxx_1[i] * fe_0 + ta_xxzz_xxxxxy_0[i] * pa_y[i] - ta_xxzz_xxxxxy_1[i] * pc_y[i];

        ta_xxyzz_xxxxxz_0[i] = ta_xxzz_xxxxxz_0[i] * pa_y[i] - ta_xxzz_xxxxxz_1[i] * pc_y[i];

        ta_xxyzz_xxxxyy_0[i] =
            2.0 * ta_xxzz_xxxxy_0[i] * fe_0 - 2.0 * ta_xxzz_xxxxy_1[i] * fe_0 + ta_xxzz_xxxxyy_0[i] * pa_y[i] - ta_xxzz_xxxxyy_1[i] * pc_y[i];

        ta_xxyzz_xxxxyz_0[i] = ta_xxzz_xxxxz_0[i] * fe_0 - ta_xxzz_xxxxz_1[i] * fe_0 + ta_xxzz_xxxxyz_0[i] * pa_y[i] - ta_xxzz_xxxxyz_1[i] * pc_y[i];

        ta_xxyzz_xxxxzz_0[i] = ta_xxzz_xxxxzz_0[i] * pa_y[i] - ta_xxzz_xxxxzz_1[i] * pc_y[i];

        ta_xxyzz_xxxyyy_0[i] =
            3.0 * ta_xxzz_xxxyy_0[i] * fe_0 - 3.0 * ta_xxzz_xxxyy_1[i] * fe_0 + ta_xxzz_xxxyyy_0[i] * pa_y[i] - ta_xxzz_xxxyyy_1[i] * pc_y[i];

        ta_xxyzz_xxxyyz_0[i] =
            2.0 * ta_xxzz_xxxyz_0[i] * fe_0 - 2.0 * ta_xxzz_xxxyz_1[i] * fe_0 + ta_xxzz_xxxyyz_0[i] * pa_y[i] - ta_xxzz_xxxyyz_1[i] * pc_y[i];

        ta_xxyzz_xxxyzz_0[i] = ta_xxzz_xxxzz_0[i] * fe_0 - ta_xxzz_xxxzz_1[i] * fe_0 + ta_xxzz_xxxyzz_0[i] * pa_y[i] - ta_xxzz_xxxyzz_1[i] * pc_y[i];

        ta_xxyzz_xxxzzz_0[i] = ta_xxzz_xxxzzz_0[i] * pa_y[i] - ta_xxzz_xxxzzz_1[i] * pc_y[i];

        ta_xxyzz_xxyyyy_0[i] =
            4.0 * ta_xxzz_xxyyy_0[i] * fe_0 - 4.0 * ta_xxzz_xxyyy_1[i] * fe_0 + ta_xxzz_xxyyyy_0[i] * pa_y[i] - ta_xxzz_xxyyyy_1[i] * pc_y[i];

        ta_xxyzz_xxyyyz_0[i] =
            3.0 * ta_xxzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxzz_xxyyz_1[i] * fe_0 + ta_xxzz_xxyyyz_0[i] * pa_y[i] - ta_xxzz_xxyyyz_1[i] * pc_y[i];

        ta_xxyzz_xxyyzz_0[i] =
            2.0 * ta_xxzz_xxyzz_0[i] * fe_0 - 2.0 * ta_xxzz_xxyzz_1[i] * fe_0 + ta_xxzz_xxyyzz_0[i] * pa_y[i] - ta_xxzz_xxyyzz_1[i] * pc_y[i];

        ta_xxyzz_xxyzzz_0[i] = ta_xxzz_xxzzz_0[i] * fe_0 - ta_xxzz_xxzzz_1[i] * fe_0 + ta_xxzz_xxyzzz_0[i] * pa_y[i] - ta_xxzz_xxyzzz_1[i] * pc_y[i];

        ta_xxyzz_xxzzzz_0[i] = ta_xxzz_xxzzzz_0[i] * pa_y[i] - ta_xxzz_xxzzzz_1[i] * pc_y[i];

        ta_xxyzz_xyyyyy_0[i] =
            5.0 * ta_xxzz_xyyyy_0[i] * fe_0 - 5.0 * ta_xxzz_xyyyy_1[i] * fe_0 + ta_xxzz_xyyyyy_0[i] * pa_y[i] - ta_xxzz_xyyyyy_1[i] * pc_y[i];

        ta_xxyzz_xyyyyz_0[i] =
            4.0 * ta_xxzz_xyyyz_0[i] * fe_0 - 4.0 * ta_xxzz_xyyyz_1[i] * fe_0 + ta_xxzz_xyyyyz_0[i] * pa_y[i] - ta_xxzz_xyyyyz_1[i] * pc_y[i];

        ta_xxyzz_xyyyzz_0[i] =
            3.0 * ta_xxzz_xyyzz_0[i] * fe_0 - 3.0 * ta_xxzz_xyyzz_1[i] * fe_0 + ta_xxzz_xyyyzz_0[i] * pa_y[i] - ta_xxzz_xyyyzz_1[i] * pc_y[i];

        ta_xxyzz_xyyzzz_0[i] =
            2.0 * ta_xxzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xxzz_xyzzz_1[i] * fe_0 + ta_xxzz_xyyzzz_0[i] * pa_y[i] - ta_xxzz_xyyzzz_1[i] * pc_y[i];

        ta_xxyzz_xyzzzz_0[i] = ta_xxzz_xzzzz_0[i] * fe_0 - ta_xxzz_xzzzz_1[i] * fe_0 + ta_xxzz_xyzzzz_0[i] * pa_y[i] - ta_xxzz_xyzzzz_1[i] * pc_y[i];

        ta_xxyzz_xzzzzz_0[i] = ta_xxzz_xzzzzz_0[i] * pa_y[i] - ta_xxzz_xzzzzz_1[i] * pc_y[i];

        ta_xxyzz_yyyyyy_0[i] = ta_yzz_yyyyyy_0[i] * fe_0 - ta_yzz_yyyyyy_1[i] * fe_0 + ta_xyzz_yyyyyy_0[i] * pa_x[i] - ta_xyzz_yyyyyy_1[i] * pc_x[i];

        ta_xxyzz_yyyyyz_0[i] = ta_yzz_yyyyyz_0[i] * fe_0 - ta_yzz_yyyyyz_1[i] * fe_0 + ta_xyzz_yyyyyz_0[i] * pa_x[i] - ta_xyzz_yyyyyz_1[i] * pc_x[i];

        ta_xxyzz_yyyyzz_0[i] = ta_yzz_yyyyzz_0[i] * fe_0 - ta_yzz_yyyyzz_1[i] * fe_0 + ta_xyzz_yyyyzz_0[i] * pa_x[i] - ta_xyzz_yyyyzz_1[i] * pc_x[i];

        ta_xxyzz_yyyzzz_0[i] = ta_yzz_yyyzzz_0[i] * fe_0 - ta_yzz_yyyzzz_1[i] * fe_0 + ta_xyzz_yyyzzz_0[i] * pa_x[i] - ta_xyzz_yyyzzz_1[i] * pc_x[i];

        ta_xxyzz_yyzzzz_0[i] = ta_yzz_yyzzzz_0[i] * fe_0 - ta_yzz_yyzzzz_1[i] * fe_0 + ta_xyzz_yyzzzz_0[i] * pa_x[i] - ta_xyzz_yyzzzz_1[i] * pc_x[i];

        ta_xxyzz_yzzzzz_0[i] = ta_yzz_yzzzzz_0[i] * fe_0 - ta_yzz_yzzzzz_1[i] * fe_0 + ta_xyzz_yzzzzz_0[i] * pa_x[i] - ta_xyzz_yzzzzz_1[i] * pc_x[i];

        ta_xxyzz_zzzzzz_0[i] = ta_xxzz_zzzzzz_0[i] * pa_y[i] - ta_xxzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 252-280 components of targeted buffer : HI

    auto ta_xxzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 252);

    auto ta_xxzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 253);

    auto ta_xxzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 254);

    auto ta_xxzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 255);

    auto ta_xxzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 256);

    auto ta_xxzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 257);

    auto ta_xxzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 258);

    auto ta_xxzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 259);

    auto ta_xxzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 260);

    auto ta_xxzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 261);

    auto ta_xxzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 262);

    auto ta_xxzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 263);

    auto ta_xxzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 264);

    auto ta_xxzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 265);

    auto ta_xxzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 266);

    auto ta_xxzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 267);

    auto ta_xxzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 268);

    auto ta_xxzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 269);

    auto ta_xxzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 270);

    auto ta_xxzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 271);

    auto ta_xxzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 272);

    auto ta_xxzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 273);

    auto ta_xxzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 274);

    auto ta_xxzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 275);

    auto ta_xxzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 276);

    auto ta_xxzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 277);

    auto ta_xxzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 278);

    auto ta_xxzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 279);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xxz_xxxxxx_0,   \
                             ta_xxz_xxxxxx_1,   \
                             ta_xxz_xxxxxy_0,   \
                             ta_xxz_xxxxxy_1,   \
                             ta_xxz_xxxxyy_0,   \
                             ta_xxz_xxxxyy_1,   \
                             ta_xxz_xxxyyy_0,   \
                             ta_xxz_xxxyyy_1,   \
                             ta_xxz_xxyyyy_0,   \
                             ta_xxz_xxyyyy_1,   \
                             ta_xxz_xyyyyy_0,   \
                             ta_xxz_xyyyyy_1,   \
                             ta_xxzz_xxxxxx_0,  \
                             ta_xxzz_xxxxxx_1,  \
                             ta_xxzz_xxxxxy_0,  \
                             ta_xxzz_xxxxxy_1,  \
                             ta_xxzz_xxxxyy_0,  \
                             ta_xxzz_xxxxyy_1,  \
                             ta_xxzz_xxxyyy_0,  \
                             ta_xxzz_xxxyyy_1,  \
                             ta_xxzz_xxyyyy_0,  \
                             ta_xxzz_xxyyyy_1,  \
                             ta_xxzz_xyyyyy_0,  \
                             ta_xxzz_xyyyyy_1,  \
                             ta_xxzzz_xxxxxx_0, \
                             ta_xxzzz_xxxxxy_0, \
                             ta_xxzzz_xxxxxz_0, \
                             ta_xxzzz_xxxxyy_0, \
                             ta_xxzzz_xxxxyz_0, \
                             ta_xxzzz_xxxxzz_0, \
                             ta_xxzzz_xxxyyy_0, \
                             ta_xxzzz_xxxyyz_0, \
                             ta_xxzzz_xxxyzz_0, \
                             ta_xxzzz_xxxzzz_0, \
                             ta_xxzzz_xxyyyy_0, \
                             ta_xxzzz_xxyyyz_0, \
                             ta_xxzzz_xxyyzz_0, \
                             ta_xxzzz_xxyzzz_0, \
                             ta_xxzzz_xxzzzz_0, \
                             ta_xxzzz_xyyyyy_0, \
                             ta_xxzzz_xyyyyz_0, \
                             ta_xxzzz_xyyyzz_0, \
                             ta_xxzzz_xyyzzz_0, \
                             ta_xxzzz_xyzzzz_0, \
                             ta_xxzzz_xzzzzz_0, \
                             ta_xxzzz_yyyyyy_0, \
                             ta_xxzzz_yyyyyz_0, \
                             ta_xxzzz_yyyyzz_0, \
                             ta_xxzzz_yyyzzz_0, \
                             ta_xxzzz_yyzzzz_0, \
                             ta_xxzzz_yzzzzz_0, \
                             ta_xxzzz_zzzzzz_0, \
                             ta_xzzz_xxxxxz_0,  \
                             ta_xzzz_xxxxxz_1,  \
                             ta_xzzz_xxxxyz_0,  \
                             ta_xzzz_xxxxyz_1,  \
                             ta_xzzz_xxxxz_0,   \
                             ta_xzzz_xxxxz_1,   \
                             ta_xzzz_xxxxzz_0,  \
                             ta_xzzz_xxxxzz_1,  \
                             ta_xzzz_xxxyyz_0,  \
                             ta_xzzz_xxxyyz_1,  \
                             ta_xzzz_xxxyz_0,   \
                             ta_xzzz_xxxyz_1,   \
                             ta_xzzz_xxxyzz_0,  \
                             ta_xzzz_xxxyzz_1,  \
                             ta_xzzz_xxxzz_0,   \
                             ta_xzzz_xxxzz_1,   \
                             ta_xzzz_xxxzzz_0,  \
                             ta_xzzz_xxxzzz_1,  \
                             ta_xzzz_xxyyyz_0,  \
                             ta_xzzz_xxyyyz_1,  \
                             ta_xzzz_xxyyz_0,   \
                             ta_xzzz_xxyyz_1,   \
                             ta_xzzz_xxyyzz_0,  \
                             ta_xzzz_xxyyzz_1,  \
                             ta_xzzz_xxyzz_0,   \
                             ta_xzzz_xxyzz_1,   \
                             ta_xzzz_xxyzzz_0,  \
                             ta_xzzz_xxyzzz_1,  \
                             ta_xzzz_xxzzz_0,   \
                             ta_xzzz_xxzzz_1,   \
                             ta_xzzz_xxzzzz_0,  \
                             ta_xzzz_xxzzzz_1,  \
                             ta_xzzz_xyyyyz_0,  \
                             ta_xzzz_xyyyyz_1,  \
                             ta_xzzz_xyyyz_0,   \
                             ta_xzzz_xyyyz_1,   \
                             ta_xzzz_xyyyzz_0,  \
                             ta_xzzz_xyyyzz_1,  \
                             ta_xzzz_xyyzz_0,   \
                             ta_xzzz_xyyzz_1,   \
                             ta_xzzz_xyyzzz_0,  \
                             ta_xzzz_xyyzzz_1,  \
                             ta_xzzz_xyzzz_0,   \
                             ta_xzzz_xyzzz_1,   \
                             ta_xzzz_xyzzzz_0,  \
                             ta_xzzz_xyzzzz_1,  \
                             ta_xzzz_xzzzz_0,   \
                             ta_xzzz_xzzzz_1,   \
                             ta_xzzz_xzzzzz_0,  \
                             ta_xzzz_xzzzzz_1,  \
                             ta_xzzz_yyyyyy_0,  \
                             ta_xzzz_yyyyyy_1,  \
                             ta_xzzz_yyyyyz_0,  \
                             ta_xzzz_yyyyyz_1,  \
                             ta_xzzz_yyyyz_0,   \
                             ta_xzzz_yyyyz_1,   \
                             ta_xzzz_yyyyzz_0,  \
                             ta_xzzz_yyyyzz_1,  \
                             ta_xzzz_yyyzz_0,   \
                             ta_xzzz_yyyzz_1,   \
                             ta_xzzz_yyyzzz_0,  \
                             ta_xzzz_yyyzzz_1,  \
                             ta_xzzz_yyzzz_0,   \
                             ta_xzzz_yyzzz_1,   \
                             ta_xzzz_yyzzzz_0,  \
                             ta_xzzz_yyzzzz_1,  \
                             ta_xzzz_yzzzz_0,   \
                             ta_xzzz_yzzzz_1,   \
                             ta_xzzz_yzzzzz_0,  \
                             ta_xzzz_yzzzzz_1,  \
                             ta_xzzz_zzzzz_0,   \
                             ta_xzzz_zzzzz_1,   \
                             ta_xzzz_zzzzzz_0,  \
                             ta_xzzz_zzzzzz_1,  \
                             ta_zzz_xxxxxz_0,   \
                             ta_zzz_xxxxxz_1,   \
                             ta_zzz_xxxxyz_0,   \
                             ta_zzz_xxxxyz_1,   \
                             ta_zzz_xxxxzz_0,   \
                             ta_zzz_xxxxzz_1,   \
                             ta_zzz_xxxyyz_0,   \
                             ta_zzz_xxxyyz_1,   \
                             ta_zzz_xxxyzz_0,   \
                             ta_zzz_xxxyzz_1,   \
                             ta_zzz_xxxzzz_0,   \
                             ta_zzz_xxxzzz_1,   \
                             ta_zzz_xxyyyz_0,   \
                             ta_zzz_xxyyyz_1,   \
                             ta_zzz_xxyyzz_0,   \
                             ta_zzz_xxyyzz_1,   \
                             ta_zzz_xxyzzz_0,   \
                             ta_zzz_xxyzzz_1,   \
                             ta_zzz_xxzzzz_0,   \
                             ta_zzz_xxzzzz_1,   \
                             ta_zzz_xyyyyz_0,   \
                             ta_zzz_xyyyyz_1,   \
                             ta_zzz_xyyyzz_0,   \
                             ta_zzz_xyyyzz_1,   \
                             ta_zzz_xyyzzz_0,   \
                             ta_zzz_xyyzzz_1,   \
                             ta_zzz_xyzzzz_0,   \
                             ta_zzz_xyzzzz_1,   \
                             ta_zzz_xzzzzz_0,   \
                             ta_zzz_xzzzzz_1,   \
                             ta_zzz_yyyyyy_0,   \
                             ta_zzz_yyyyyy_1,   \
                             ta_zzz_yyyyyz_0,   \
                             ta_zzz_yyyyyz_1,   \
                             ta_zzz_yyyyzz_0,   \
                             ta_zzz_yyyyzz_1,   \
                             ta_zzz_yyyzzz_0,   \
                             ta_zzz_yyyzzz_1,   \
                             ta_zzz_yyzzzz_0,   \
                             ta_zzz_yyzzzz_1,   \
                             ta_zzz_yzzzzz_0,   \
                             ta_zzz_yzzzzz_1,   \
                             ta_zzz_zzzzzz_0,   \
                             ta_zzz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzz_xxxxxx_0[i] =
            2.0 * ta_xxz_xxxxxx_0[i] * fe_0 - 2.0 * ta_xxz_xxxxxx_1[i] * fe_0 + ta_xxzz_xxxxxx_0[i] * pa_z[i] - ta_xxzz_xxxxxx_1[i] * pc_z[i];

        ta_xxzzz_xxxxxy_0[i] =
            2.0 * ta_xxz_xxxxxy_0[i] * fe_0 - 2.0 * ta_xxz_xxxxxy_1[i] * fe_0 + ta_xxzz_xxxxxy_0[i] * pa_z[i] - ta_xxzz_xxxxxy_1[i] * pc_z[i];

        ta_xxzzz_xxxxxz_0[i] = ta_zzz_xxxxxz_0[i] * fe_0 - ta_zzz_xxxxxz_1[i] * fe_0 + 5.0 * ta_xzzz_xxxxz_0[i] * fe_0 -
                               5.0 * ta_xzzz_xxxxz_1[i] * fe_0 + ta_xzzz_xxxxxz_0[i] * pa_x[i] - ta_xzzz_xxxxxz_1[i] * pc_x[i];

        ta_xxzzz_xxxxyy_0[i] =
            2.0 * ta_xxz_xxxxyy_0[i] * fe_0 - 2.0 * ta_xxz_xxxxyy_1[i] * fe_0 + ta_xxzz_xxxxyy_0[i] * pa_z[i] - ta_xxzz_xxxxyy_1[i] * pc_z[i];

        ta_xxzzz_xxxxyz_0[i] = ta_zzz_xxxxyz_0[i] * fe_0 - ta_zzz_xxxxyz_1[i] * fe_0 + 4.0 * ta_xzzz_xxxyz_0[i] * fe_0 -
                               4.0 * ta_xzzz_xxxyz_1[i] * fe_0 + ta_xzzz_xxxxyz_0[i] * pa_x[i] - ta_xzzz_xxxxyz_1[i] * pc_x[i];

        ta_xxzzz_xxxxzz_0[i] = ta_zzz_xxxxzz_0[i] * fe_0 - ta_zzz_xxxxzz_1[i] * fe_0 + 4.0 * ta_xzzz_xxxzz_0[i] * fe_0 -
                               4.0 * ta_xzzz_xxxzz_1[i] * fe_0 + ta_xzzz_xxxxzz_0[i] * pa_x[i] - ta_xzzz_xxxxzz_1[i] * pc_x[i];

        ta_xxzzz_xxxyyy_0[i] =
            2.0 * ta_xxz_xxxyyy_0[i] * fe_0 - 2.0 * ta_xxz_xxxyyy_1[i] * fe_0 + ta_xxzz_xxxyyy_0[i] * pa_z[i] - ta_xxzz_xxxyyy_1[i] * pc_z[i];

        ta_xxzzz_xxxyyz_0[i] = ta_zzz_xxxyyz_0[i] * fe_0 - ta_zzz_xxxyyz_1[i] * fe_0 + 3.0 * ta_xzzz_xxyyz_0[i] * fe_0 -
                               3.0 * ta_xzzz_xxyyz_1[i] * fe_0 + ta_xzzz_xxxyyz_0[i] * pa_x[i] - ta_xzzz_xxxyyz_1[i] * pc_x[i];

        ta_xxzzz_xxxyzz_0[i] = ta_zzz_xxxyzz_0[i] * fe_0 - ta_zzz_xxxyzz_1[i] * fe_0 + 3.0 * ta_xzzz_xxyzz_0[i] * fe_0 -
                               3.0 * ta_xzzz_xxyzz_1[i] * fe_0 + ta_xzzz_xxxyzz_0[i] * pa_x[i] - ta_xzzz_xxxyzz_1[i] * pc_x[i];

        ta_xxzzz_xxxzzz_0[i] = ta_zzz_xxxzzz_0[i] * fe_0 - ta_zzz_xxxzzz_1[i] * fe_0 + 3.0 * ta_xzzz_xxzzz_0[i] * fe_0 -
                               3.0 * ta_xzzz_xxzzz_1[i] * fe_0 + ta_xzzz_xxxzzz_0[i] * pa_x[i] - ta_xzzz_xxxzzz_1[i] * pc_x[i];

        ta_xxzzz_xxyyyy_0[i] =
            2.0 * ta_xxz_xxyyyy_0[i] * fe_0 - 2.0 * ta_xxz_xxyyyy_1[i] * fe_0 + ta_xxzz_xxyyyy_0[i] * pa_z[i] - ta_xxzz_xxyyyy_1[i] * pc_z[i];

        ta_xxzzz_xxyyyz_0[i] = ta_zzz_xxyyyz_0[i] * fe_0 - ta_zzz_xxyyyz_1[i] * fe_0 + 2.0 * ta_xzzz_xyyyz_0[i] * fe_0 -
                               2.0 * ta_xzzz_xyyyz_1[i] * fe_0 + ta_xzzz_xxyyyz_0[i] * pa_x[i] - ta_xzzz_xxyyyz_1[i] * pc_x[i];

        ta_xxzzz_xxyyzz_0[i] = ta_zzz_xxyyzz_0[i] * fe_0 - ta_zzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_xzzz_xyyzz_0[i] * fe_0 -
                               2.0 * ta_xzzz_xyyzz_1[i] * fe_0 + ta_xzzz_xxyyzz_0[i] * pa_x[i] - ta_xzzz_xxyyzz_1[i] * pc_x[i];

        ta_xxzzz_xxyzzz_0[i] = ta_zzz_xxyzzz_0[i] * fe_0 - ta_zzz_xxyzzz_1[i] * fe_0 + 2.0 * ta_xzzz_xyzzz_0[i] * fe_0 -
                               2.0 * ta_xzzz_xyzzz_1[i] * fe_0 + ta_xzzz_xxyzzz_0[i] * pa_x[i] - ta_xzzz_xxyzzz_1[i] * pc_x[i];

        ta_xxzzz_xxzzzz_0[i] = ta_zzz_xxzzzz_0[i] * fe_0 - ta_zzz_xxzzzz_1[i] * fe_0 + 2.0 * ta_xzzz_xzzzz_0[i] * fe_0 -
                               2.0 * ta_xzzz_xzzzz_1[i] * fe_0 + ta_xzzz_xxzzzz_0[i] * pa_x[i] - ta_xzzz_xxzzzz_1[i] * pc_x[i];

        ta_xxzzz_xyyyyy_0[i] =
            2.0 * ta_xxz_xyyyyy_0[i] * fe_0 - 2.0 * ta_xxz_xyyyyy_1[i] * fe_0 + ta_xxzz_xyyyyy_0[i] * pa_z[i] - ta_xxzz_xyyyyy_1[i] * pc_z[i];

        ta_xxzzz_xyyyyz_0[i] = ta_zzz_xyyyyz_0[i] * fe_0 - ta_zzz_xyyyyz_1[i] * fe_0 + ta_xzzz_yyyyz_0[i] * fe_0 - ta_xzzz_yyyyz_1[i] * fe_0 +
                               ta_xzzz_xyyyyz_0[i] * pa_x[i] - ta_xzzz_xyyyyz_1[i] * pc_x[i];

        ta_xxzzz_xyyyzz_0[i] = ta_zzz_xyyyzz_0[i] * fe_0 - ta_zzz_xyyyzz_1[i] * fe_0 + ta_xzzz_yyyzz_0[i] * fe_0 - ta_xzzz_yyyzz_1[i] * fe_0 +
                               ta_xzzz_xyyyzz_0[i] * pa_x[i] - ta_xzzz_xyyyzz_1[i] * pc_x[i];

        ta_xxzzz_xyyzzz_0[i] = ta_zzz_xyyzzz_0[i] * fe_0 - ta_zzz_xyyzzz_1[i] * fe_0 + ta_xzzz_yyzzz_0[i] * fe_0 - ta_xzzz_yyzzz_1[i] * fe_0 +
                               ta_xzzz_xyyzzz_0[i] * pa_x[i] - ta_xzzz_xyyzzz_1[i] * pc_x[i];

        ta_xxzzz_xyzzzz_0[i] = ta_zzz_xyzzzz_0[i] * fe_0 - ta_zzz_xyzzzz_1[i] * fe_0 + ta_xzzz_yzzzz_0[i] * fe_0 - ta_xzzz_yzzzz_1[i] * fe_0 +
                               ta_xzzz_xyzzzz_0[i] * pa_x[i] - ta_xzzz_xyzzzz_1[i] * pc_x[i];

        ta_xxzzz_xzzzzz_0[i] = ta_zzz_xzzzzz_0[i] * fe_0 - ta_zzz_xzzzzz_1[i] * fe_0 + ta_xzzz_zzzzz_0[i] * fe_0 - ta_xzzz_zzzzz_1[i] * fe_0 +
                               ta_xzzz_xzzzzz_0[i] * pa_x[i] - ta_xzzz_xzzzzz_1[i] * pc_x[i];

        ta_xxzzz_yyyyyy_0[i] = ta_zzz_yyyyyy_0[i] * fe_0 - ta_zzz_yyyyyy_1[i] * fe_0 + ta_xzzz_yyyyyy_0[i] * pa_x[i] - ta_xzzz_yyyyyy_1[i] * pc_x[i];

        ta_xxzzz_yyyyyz_0[i] = ta_zzz_yyyyyz_0[i] * fe_0 - ta_zzz_yyyyyz_1[i] * fe_0 + ta_xzzz_yyyyyz_0[i] * pa_x[i] - ta_xzzz_yyyyyz_1[i] * pc_x[i];

        ta_xxzzz_yyyyzz_0[i] = ta_zzz_yyyyzz_0[i] * fe_0 - ta_zzz_yyyyzz_1[i] * fe_0 + ta_xzzz_yyyyzz_0[i] * pa_x[i] - ta_xzzz_yyyyzz_1[i] * pc_x[i];

        ta_xxzzz_yyyzzz_0[i] = ta_zzz_yyyzzz_0[i] * fe_0 - ta_zzz_yyyzzz_1[i] * fe_0 + ta_xzzz_yyyzzz_0[i] * pa_x[i] - ta_xzzz_yyyzzz_1[i] * pc_x[i];

        ta_xxzzz_yyzzzz_0[i] = ta_zzz_yyzzzz_0[i] * fe_0 - ta_zzz_yyzzzz_1[i] * fe_0 + ta_xzzz_yyzzzz_0[i] * pa_x[i] - ta_xzzz_yyzzzz_1[i] * pc_x[i];

        ta_xxzzz_yzzzzz_0[i] = ta_zzz_yzzzzz_0[i] * fe_0 - ta_zzz_yzzzzz_1[i] * fe_0 + ta_xzzz_yzzzzz_0[i] * pa_x[i] - ta_xzzz_yzzzzz_1[i] * pc_x[i];

        ta_xxzzz_zzzzzz_0[i] = ta_zzz_zzzzzz_0[i] * fe_0 - ta_zzz_zzzzzz_1[i] * fe_0 + ta_xzzz_zzzzzz_0[i] * pa_x[i] - ta_xzzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 280-308 components of targeted buffer : HI

    auto ta_xyyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 280);

    auto ta_xyyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 281);

    auto ta_xyyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 282);

    auto ta_xyyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 283);

    auto ta_xyyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 284);

    auto ta_xyyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 285);

    auto ta_xyyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 286);

    auto ta_xyyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 287);

    auto ta_xyyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 288);

    auto ta_xyyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 289);

    auto ta_xyyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 290);

    auto ta_xyyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 291);

    auto ta_xyyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 292);

    auto ta_xyyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 293);

    auto ta_xyyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 294);

    auto ta_xyyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 295);

    auto ta_xyyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 296);

    auto ta_xyyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 297);

    auto ta_xyyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 298);

    auto ta_xyyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 299);

    auto ta_xyyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 300);

    auto ta_xyyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 301);

    auto ta_xyyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 302);

    auto ta_xyyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 303);

    auto ta_xyyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 304);

    auto ta_xyyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 305);

    auto ta_xyyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 306);

    auto ta_xyyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 307);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xyyyy_xxxxxx_0, \
                             ta_xyyyy_xxxxxy_0, \
                             ta_xyyyy_xxxxxz_0, \
                             ta_xyyyy_xxxxyy_0, \
                             ta_xyyyy_xxxxyz_0, \
                             ta_xyyyy_xxxxzz_0, \
                             ta_xyyyy_xxxyyy_0, \
                             ta_xyyyy_xxxyyz_0, \
                             ta_xyyyy_xxxyzz_0, \
                             ta_xyyyy_xxxzzz_0, \
                             ta_xyyyy_xxyyyy_0, \
                             ta_xyyyy_xxyyyz_0, \
                             ta_xyyyy_xxyyzz_0, \
                             ta_xyyyy_xxyzzz_0, \
                             ta_xyyyy_xxzzzz_0, \
                             ta_xyyyy_xyyyyy_0, \
                             ta_xyyyy_xyyyyz_0, \
                             ta_xyyyy_xyyyzz_0, \
                             ta_xyyyy_xyyzzz_0, \
                             ta_xyyyy_xyzzzz_0, \
                             ta_xyyyy_xzzzzz_0, \
                             ta_xyyyy_yyyyyy_0, \
                             ta_xyyyy_yyyyyz_0, \
                             ta_xyyyy_yyyyzz_0, \
                             ta_xyyyy_yyyzzz_0, \
                             ta_xyyyy_yyzzzz_0, \
                             ta_xyyyy_yzzzzz_0, \
                             ta_xyyyy_zzzzzz_0, \
                             ta_yyyy_xxxxx_0,   \
                             ta_yyyy_xxxxx_1,   \
                             ta_yyyy_xxxxxx_0,  \
                             ta_yyyy_xxxxxx_1,  \
                             ta_yyyy_xxxxxy_0,  \
                             ta_yyyy_xxxxxy_1,  \
                             ta_yyyy_xxxxxz_0,  \
                             ta_yyyy_xxxxxz_1,  \
                             ta_yyyy_xxxxy_0,   \
                             ta_yyyy_xxxxy_1,   \
                             ta_yyyy_xxxxyy_0,  \
                             ta_yyyy_xxxxyy_1,  \
                             ta_yyyy_xxxxyz_0,  \
                             ta_yyyy_xxxxyz_1,  \
                             ta_yyyy_xxxxz_0,   \
                             ta_yyyy_xxxxz_1,   \
                             ta_yyyy_xxxxzz_0,  \
                             ta_yyyy_xxxxzz_1,  \
                             ta_yyyy_xxxyy_0,   \
                             ta_yyyy_xxxyy_1,   \
                             ta_yyyy_xxxyyy_0,  \
                             ta_yyyy_xxxyyy_1,  \
                             ta_yyyy_xxxyyz_0,  \
                             ta_yyyy_xxxyyz_1,  \
                             ta_yyyy_xxxyz_0,   \
                             ta_yyyy_xxxyz_1,   \
                             ta_yyyy_xxxyzz_0,  \
                             ta_yyyy_xxxyzz_1,  \
                             ta_yyyy_xxxzz_0,   \
                             ta_yyyy_xxxzz_1,   \
                             ta_yyyy_xxxzzz_0,  \
                             ta_yyyy_xxxzzz_1,  \
                             ta_yyyy_xxyyy_0,   \
                             ta_yyyy_xxyyy_1,   \
                             ta_yyyy_xxyyyy_0,  \
                             ta_yyyy_xxyyyy_1,  \
                             ta_yyyy_xxyyyz_0,  \
                             ta_yyyy_xxyyyz_1,  \
                             ta_yyyy_xxyyz_0,   \
                             ta_yyyy_xxyyz_1,   \
                             ta_yyyy_xxyyzz_0,  \
                             ta_yyyy_xxyyzz_1,  \
                             ta_yyyy_xxyzz_0,   \
                             ta_yyyy_xxyzz_1,   \
                             ta_yyyy_xxyzzz_0,  \
                             ta_yyyy_xxyzzz_1,  \
                             ta_yyyy_xxzzz_0,   \
                             ta_yyyy_xxzzz_1,   \
                             ta_yyyy_xxzzzz_0,  \
                             ta_yyyy_xxzzzz_1,  \
                             ta_yyyy_xyyyy_0,   \
                             ta_yyyy_xyyyy_1,   \
                             ta_yyyy_xyyyyy_0,  \
                             ta_yyyy_xyyyyy_1,  \
                             ta_yyyy_xyyyyz_0,  \
                             ta_yyyy_xyyyyz_1,  \
                             ta_yyyy_xyyyz_0,   \
                             ta_yyyy_xyyyz_1,   \
                             ta_yyyy_xyyyzz_0,  \
                             ta_yyyy_xyyyzz_1,  \
                             ta_yyyy_xyyzz_0,   \
                             ta_yyyy_xyyzz_1,   \
                             ta_yyyy_xyyzzz_0,  \
                             ta_yyyy_xyyzzz_1,  \
                             ta_yyyy_xyzzz_0,   \
                             ta_yyyy_xyzzz_1,   \
                             ta_yyyy_xyzzzz_0,  \
                             ta_yyyy_xyzzzz_1,  \
                             ta_yyyy_xzzzz_0,   \
                             ta_yyyy_xzzzz_1,   \
                             ta_yyyy_xzzzzz_0,  \
                             ta_yyyy_xzzzzz_1,  \
                             ta_yyyy_yyyyy_0,   \
                             ta_yyyy_yyyyy_1,   \
                             ta_yyyy_yyyyyy_0,  \
                             ta_yyyy_yyyyyy_1,  \
                             ta_yyyy_yyyyyz_0,  \
                             ta_yyyy_yyyyyz_1,  \
                             ta_yyyy_yyyyz_0,   \
                             ta_yyyy_yyyyz_1,   \
                             ta_yyyy_yyyyzz_0,  \
                             ta_yyyy_yyyyzz_1,  \
                             ta_yyyy_yyyzz_0,   \
                             ta_yyyy_yyyzz_1,   \
                             ta_yyyy_yyyzzz_0,  \
                             ta_yyyy_yyyzzz_1,  \
                             ta_yyyy_yyzzz_0,   \
                             ta_yyyy_yyzzz_1,   \
                             ta_yyyy_yyzzzz_0,  \
                             ta_yyyy_yyzzzz_1,  \
                             ta_yyyy_yzzzz_0,   \
                             ta_yyyy_yzzzz_1,   \
                             ta_yyyy_yzzzzz_0,  \
                             ta_yyyy_yzzzzz_1,  \
                             ta_yyyy_zzzzz_0,   \
                             ta_yyyy_zzzzz_1,   \
                             ta_yyyy_zzzzzz_0,  \
                             ta_yyyy_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyy_xxxxxx_0[i] =
            6.0 * ta_yyyy_xxxxx_0[i] * fe_0 - 6.0 * ta_yyyy_xxxxx_1[i] * fe_0 + ta_yyyy_xxxxxx_0[i] * pa_x[i] - ta_yyyy_xxxxxx_1[i] * pc_x[i];

        ta_xyyyy_xxxxxy_0[i] =
            5.0 * ta_yyyy_xxxxy_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxy_1[i] * fe_0 + ta_yyyy_xxxxxy_0[i] * pa_x[i] - ta_yyyy_xxxxxy_1[i] * pc_x[i];

        ta_xyyyy_xxxxxz_0[i] =
            5.0 * ta_yyyy_xxxxz_0[i] * fe_0 - 5.0 * ta_yyyy_xxxxz_1[i] * fe_0 + ta_yyyy_xxxxxz_0[i] * pa_x[i] - ta_yyyy_xxxxxz_1[i] * pc_x[i];

        ta_xyyyy_xxxxyy_0[i] =
            4.0 * ta_yyyy_xxxyy_0[i] * fe_0 - 4.0 * ta_yyyy_xxxyy_1[i] * fe_0 + ta_yyyy_xxxxyy_0[i] * pa_x[i] - ta_yyyy_xxxxyy_1[i] * pc_x[i];

        ta_xyyyy_xxxxyz_0[i] =
            4.0 * ta_yyyy_xxxyz_0[i] * fe_0 - 4.0 * ta_yyyy_xxxyz_1[i] * fe_0 + ta_yyyy_xxxxyz_0[i] * pa_x[i] - ta_yyyy_xxxxyz_1[i] * pc_x[i];

        ta_xyyyy_xxxxzz_0[i] =
            4.0 * ta_yyyy_xxxzz_0[i] * fe_0 - 4.0 * ta_yyyy_xxxzz_1[i] * fe_0 + ta_yyyy_xxxxzz_0[i] * pa_x[i] - ta_yyyy_xxxxzz_1[i] * pc_x[i];

        ta_xyyyy_xxxyyy_0[i] =
            3.0 * ta_yyyy_xxyyy_0[i] * fe_0 - 3.0 * ta_yyyy_xxyyy_1[i] * fe_0 + ta_yyyy_xxxyyy_0[i] * pa_x[i] - ta_yyyy_xxxyyy_1[i] * pc_x[i];

        ta_xyyyy_xxxyyz_0[i] =
            3.0 * ta_yyyy_xxyyz_0[i] * fe_0 - 3.0 * ta_yyyy_xxyyz_1[i] * fe_0 + ta_yyyy_xxxyyz_0[i] * pa_x[i] - ta_yyyy_xxxyyz_1[i] * pc_x[i];

        ta_xyyyy_xxxyzz_0[i] =
            3.0 * ta_yyyy_xxyzz_0[i] * fe_0 - 3.0 * ta_yyyy_xxyzz_1[i] * fe_0 + ta_yyyy_xxxyzz_0[i] * pa_x[i] - ta_yyyy_xxxyzz_1[i] * pc_x[i];

        ta_xyyyy_xxxzzz_0[i] =
            3.0 * ta_yyyy_xxzzz_0[i] * fe_0 - 3.0 * ta_yyyy_xxzzz_1[i] * fe_0 + ta_yyyy_xxxzzz_0[i] * pa_x[i] - ta_yyyy_xxxzzz_1[i] * pc_x[i];

        ta_xyyyy_xxyyyy_0[i] =
            2.0 * ta_yyyy_xyyyy_0[i] * fe_0 - 2.0 * ta_yyyy_xyyyy_1[i] * fe_0 + ta_yyyy_xxyyyy_0[i] * pa_x[i] - ta_yyyy_xxyyyy_1[i] * pc_x[i];

        ta_xyyyy_xxyyyz_0[i] =
            2.0 * ta_yyyy_xyyyz_0[i] * fe_0 - 2.0 * ta_yyyy_xyyyz_1[i] * fe_0 + ta_yyyy_xxyyyz_0[i] * pa_x[i] - ta_yyyy_xxyyyz_1[i] * pc_x[i];

        ta_xyyyy_xxyyzz_0[i] =
            2.0 * ta_yyyy_xyyzz_0[i] * fe_0 - 2.0 * ta_yyyy_xyyzz_1[i] * fe_0 + ta_yyyy_xxyyzz_0[i] * pa_x[i] - ta_yyyy_xxyyzz_1[i] * pc_x[i];

        ta_xyyyy_xxyzzz_0[i] =
            2.0 * ta_yyyy_xyzzz_0[i] * fe_0 - 2.0 * ta_yyyy_xyzzz_1[i] * fe_0 + ta_yyyy_xxyzzz_0[i] * pa_x[i] - ta_yyyy_xxyzzz_1[i] * pc_x[i];

        ta_xyyyy_xxzzzz_0[i] =
            2.0 * ta_yyyy_xzzzz_0[i] * fe_0 - 2.0 * ta_yyyy_xzzzz_1[i] * fe_0 + ta_yyyy_xxzzzz_0[i] * pa_x[i] - ta_yyyy_xxzzzz_1[i] * pc_x[i];

        ta_xyyyy_xyyyyy_0[i] = ta_yyyy_yyyyy_0[i] * fe_0 - ta_yyyy_yyyyy_1[i] * fe_0 + ta_yyyy_xyyyyy_0[i] * pa_x[i] - ta_yyyy_xyyyyy_1[i] * pc_x[i];

        ta_xyyyy_xyyyyz_0[i] = ta_yyyy_yyyyz_0[i] * fe_0 - ta_yyyy_yyyyz_1[i] * fe_0 + ta_yyyy_xyyyyz_0[i] * pa_x[i] - ta_yyyy_xyyyyz_1[i] * pc_x[i];

        ta_xyyyy_xyyyzz_0[i] = ta_yyyy_yyyzz_0[i] * fe_0 - ta_yyyy_yyyzz_1[i] * fe_0 + ta_yyyy_xyyyzz_0[i] * pa_x[i] - ta_yyyy_xyyyzz_1[i] * pc_x[i];

        ta_xyyyy_xyyzzz_0[i] = ta_yyyy_yyzzz_0[i] * fe_0 - ta_yyyy_yyzzz_1[i] * fe_0 + ta_yyyy_xyyzzz_0[i] * pa_x[i] - ta_yyyy_xyyzzz_1[i] * pc_x[i];

        ta_xyyyy_xyzzzz_0[i] = ta_yyyy_yzzzz_0[i] * fe_0 - ta_yyyy_yzzzz_1[i] * fe_0 + ta_yyyy_xyzzzz_0[i] * pa_x[i] - ta_yyyy_xyzzzz_1[i] * pc_x[i];

        ta_xyyyy_xzzzzz_0[i] = ta_yyyy_zzzzz_0[i] * fe_0 - ta_yyyy_zzzzz_1[i] * fe_0 + ta_yyyy_xzzzzz_0[i] * pa_x[i] - ta_yyyy_xzzzzz_1[i] * pc_x[i];

        ta_xyyyy_yyyyyy_0[i] = ta_yyyy_yyyyyy_0[i] * pa_x[i] - ta_yyyy_yyyyyy_1[i] * pc_x[i];

        ta_xyyyy_yyyyyz_0[i] = ta_yyyy_yyyyyz_0[i] * pa_x[i] - ta_yyyy_yyyyyz_1[i] * pc_x[i];

        ta_xyyyy_yyyyzz_0[i] = ta_yyyy_yyyyzz_0[i] * pa_x[i] - ta_yyyy_yyyyzz_1[i] * pc_x[i];

        ta_xyyyy_yyyzzz_0[i] = ta_yyyy_yyyzzz_0[i] * pa_x[i] - ta_yyyy_yyyzzz_1[i] * pc_x[i];

        ta_xyyyy_yyzzzz_0[i] = ta_yyyy_yyzzzz_0[i] * pa_x[i] - ta_yyyy_yyzzzz_1[i] * pc_x[i];

        ta_xyyyy_yzzzzz_0[i] = ta_yyyy_yzzzzz_0[i] * pa_x[i] - ta_yyyy_yzzzzz_1[i] * pc_x[i];

        ta_xyyyy_zzzzzz_0[i] = ta_yyyy_zzzzzz_0[i] * pa_x[i] - ta_yyyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 308-336 components of targeted buffer : HI

    auto ta_xyyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 308);

    auto ta_xyyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 309);

    auto ta_xyyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 310);

    auto ta_xyyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 311);

    auto ta_xyyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 312);

    auto ta_xyyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 313);

    auto ta_xyyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 314);

    auto ta_xyyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 315);

    auto ta_xyyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 316);

    auto ta_xyyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 317);

    auto ta_xyyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 318);

    auto ta_xyyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 319);

    auto ta_xyyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 320);

    auto ta_xyyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 321);

    auto ta_xyyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 322);

    auto ta_xyyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 323);

    auto ta_xyyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 324);

    auto ta_xyyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 325);

    auto ta_xyyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 326);

    auto ta_xyyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 327);

    auto ta_xyyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 328);

    auto ta_xyyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 329);

    auto ta_xyyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 330);

    auto ta_xyyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 331);

    auto ta_xyyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 332);

    auto ta_xyyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 333);

    auto ta_xyyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 334);

    auto ta_xyyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 335);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta_xyyy_xxxxxx_0,  \
                             ta_xyyy_xxxxxx_1,  \
                             ta_xyyy_xxxxxy_0,  \
                             ta_xyyy_xxxxxy_1,  \
                             ta_xyyy_xxxxyy_0,  \
                             ta_xyyy_xxxxyy_1,  \
                             ta_xyyy_xxxyyy_0,  \
                             ta_xyyy_xxxyyy_1,  \
                             ta_xyyy_xxyyyy_0,  \
                             ta_xyyy_xxyyyy_1,  \
                             ta_xyyy_xyyyyy_0,  \
                             ta_xyyy_xyyyyy_1,  \
                             ta_xyyyz_xxxxxx_0, \
                             ta_xyyyz_xxxxxy_0, \
                             ta_xyyyz_xxxxxz_0, \
                             ta_xyyyz_xxxxyy_0, \
                             ta_xyyyz_xxxxyz_0, \
                             ta_xyyyz_xxxxzz_0, \
                             ta_xyyyz_xxxyyy_0, \
                             ta_xyyyz_xxxyyz_0, \
                             ta_xyyyz_xxxyzz_0, \
                             ta_xyyyz_xxxzzz_0, \
                             ta_xyyyz_xxyyyy_0, \
                             ta_xyyyz_xxyyyz_0, \
                             ta_xyyyz_xxyyzz_0, \
                             ta_xyyyz_xxyzzz_0, \
                             ta_xyyyz_xxzzzz_0, \
                             ta_xyyyz_xyyyyy_0, \
                             ta_xyyyz_xyyyyz_0, \
                             ta_xyyyz_xyyyzz_0, \
                             ta_xyyyz_xyyzzz_0, \
                             ta_xyyyz_xyzzzz_0, \
                             ta_xyyyz_xzzzzz_0, \
                             ta_xyyyz_yyyyyy_0, \
                             ta_xyyyz_yyyyyz_0, \
                             ta_xyyyz_yyyyzz_0, \
                             ta_xyyyz_yyyzzz_0, \
                             ta_xyyyz_yyzzzz_0, \
                             ta_xyyyz_yzzzzz_0, \
                             ta_xyyyz_zzzzzz_0, \
                             ta_yyyz_xxxxxz_0,  \
                             ta_yyyz_xxxxxz_1,  \
                             ta_yyyz_xxxxyz_0,  \
                             ta_yyyz_xxxxyz_1,  \
                             ta_yyyz_xxxxz_0,   \
                             ta_yyyz_xxxxz_1,   \
                             ta_yyyz_xxxxzz_0,  \
                             ta_yyyz_xxxxzz_1,  \
                             ta_yyyz_xxxyyz_0,  \
                             ta_yyyz_xxxyyz_1,  \
                             ta_yyyz_xxxyz_0,   \
                             ta_yyyz_xxxyz_1,   \
                             ta_yyyz_xxxyzz_0,  \
                             ta_yyyz_xxxyzz_1,  \
                             ta_yyyz_xxxzz_0,   \
                             ta_yyyz_xxxzz_1,   \
                             ta_yyyz_xxxzzz_0,  \
                             ta_yyyz_xxxzzz_1,  \
                             ta_yyyz_xxyyyz_0,  \
                             ta_yyyz_xxyyyz_1,  \
                             ta_yyyz_xxyyz_0,   \
                             ta_yyyz_xxyyz_1,   \
                             ta_yyyz_xxyyzz_0,  \
                             ta_yyyz_xxyyzz_1,  \
                             ta_yyyz_xxyzz_0,   \
                             ta_yyyz_xxyzz_1,   \
                             ta_yyyz_xxyzzz_0,  \
                             ta_yyyz_xxyzzz_1,  \
                             ta_yyyz_xxzzz_0,   \
                             ta_yyyz_xxzzz_1,   \
                             ta_yyyz_xxzzzz_0,  \
                             ta_yyyz_xxzzzz_1,  \
                             ta_yyyz_xyyyyz_0,  \
                             ta_yyyz_xyyyyz_1,  \
                             ta_yyyz_xyyyz_0,   \
                             ta_yyyz_xyyyz_1,   \
                             ta_yyyz_xyyyzz_0,  \
                             ta_yyyz_xyyyzz_1,  \
                             ta_yyyz_xyyzz_0,   \
                             ta_yyyz_xyyzz_1,   \
                             ta_yyyz_xyyzzz_0,  \
                             ta_yyyz_xyyzzz_1,  \
                             ta_yyyz_xyzzz_0,   \
                             ta_yyyz_xyzzz_1,   \
                             ta_yyyz_xyzzzz_0,  \
                             ta_yyyz_xyzzzz_1,  \
                             ta_yyyz_xzzzz_0,   \
                             ta_yyyz_xzzzz_1,   \
                             ta_yyyz_xzzzzz_0,  \
                             ta_yyyz_xzzzzz_1,  \
                             ta_yyyz_yyyyyy_0,  \
                             ta_yyyz_yyyyyy_1,  \
                             ta_yyyz_yyyyyz_0,  \
                             ta_yyyz_yyyyyz_1,  \
                             ta_yyyz_yyyyz_0,   \
                             ta_yyyz_yyyyz_1,   \
                             ta_yyyz_yyyyzz_0,  \
                             ta_yyyz_yyyyzz_1,  \
                             ta_yyyz_yyyzz_0,   \
                             ta_yyyz_yyyzz_1,   \
                             ta_yyyz_yyyzzz_0,  \
                             ta_yyyz_yyyzzz_1,  \
                             ta_yyyz_yyzzz_0,   \
                             ta_yyyz_yyzzz_1,   \
                             ta_yyyz_yyzzzz_0,  \
                             ta_yyyz_yyzzzz_1,  \
                             ta_yyyz_yzzzz_0,   \
                             ta_yyyz_yzzzz_1,   \
                             ta_yyyz_yzzzzz_0,  \
                             ta_yyyz_yzzzzz_1,  \
                             ta_yyyz_zzzzz_0,   \
                             ta_yyyz_zzzzz_1,   \
                             ta_yyyz_zzzzzz_0,  \
                             ta_yyyz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyz_xxxxxx_0[i] = ta_xyyy_xxxxxx_0[i] * pa_z[i] - ta_xyyy_xxxxxx_1[i] * pc_z[i];

        ta_xyyyz_xxxxxy_0[i] = ta_xyyy_xxxxxy_0[i] * pa_z[i] - ta_xyyy_xxxxxy_1[i] * pc_z[i];

        ta_xyyyz_xxxxxz_0[i] =
            5.0 * ta_yyyz_xxxxz_0[i] * fe_0 - 5.0 * ta_yyyz_xxxxz_1[i] * fe_0 + ta_yyyz_xxxxxz_0[i] * pa_x[i] - ta_yyyz_xxxxxz_1[i] * pc_x[i];

        ta_xyyyz_xxxxyy_0[i] = ta_xyyy_xxxxyy_0[i] * pa_z[i] - ta_xyyy_xxxxyy_1[i] * pc_z[i];

        ta_xyyyz_xxxxyz_0[i] =
            4.0 * ta_yyyz_xxxyz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxyz_1[i] * fe_0 + ta_yyyz_xxxxyz_0[i] * pa_x[i] - ta_yyyz_xxxxyz_1[i] * pc_x[i];

        ta_xyyyz_xxxxzz_0[i] =
            4.0 * ta_yyyz_xxxzz_0[i] * fe_0 - 4.0 * ta_yyyz_xxxzz_1[i] * fe_0 + ta_yyyz_xxxxzz_0[i] * pa_x[i] - ta_yyyz_xxxxzz_1[i] * pc_x[i];

        ta_xyyyz_xxxyyy_0[i] = ta_xyyy_xxxyyy_0[i] * pa_z[i] - ta_xyyy_xxxyyy_1[i] * pc_z[i];

        ta_xyyyz_xxxyyz_0[i] =
            3.0 * ta_yyyz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyyz_xxyyz_1[i] * fe_0 + ta_yyyz_xxxyyz_0[i] * pa_x[i] - ta_yyyz_xxxyyz_1[i] * pc_x[i];

        ta_xyyyz_xxxyzz_0[i] =
            3.0 * ta_yyyz_xxyzz_0[i] * fe_0 - 3.0 * ta_yyyz_xxyzz_1[i] * fe_0 + ta_yyyz_xxxyzz_0[i] * pa_x[i] - ta_yyyz_xxxyzz_1[i] * pc_x[i];

        ta_xyyyz_xxxzzz_0[i] =
            3.0 * ta_yyyz_xxzzz_0[i] * fe_0 - 3.0 * ta_yyyz_xxzzz_1[i] * fe_0 + ta_yyyz_xxxzzz_0[i] * pa_x[i] - ta_yyyz_xxxzzz_1[i] * pc_x[i];

        ta_xyyyz_xxyyyy_0[i] = ta_xyyy_xxyyyy_0[i] * pa_z[i] - ta_xyyy_xxyyyy_1[i] * pc_z[i];

        ta_xyyyz_xxyyyz_0[i] =
            2.0 * ta_yyyz_xyyyz_0[i] * fe_0 - 2.0 * ta_yyyz_xyyyz_1[i] * fe_0 + ta_yyyz_xxyyyz_0[i] * pa_x[i] - ta_yyyz_xxyyyz_1[i] * pc_x[i];

        ta_xyyyz_xxyyzz_0[i] =
            2.0 * ta_yyyz_xyyzz_0[i] * fe_0 - 2.0 * ta_yyyz_xyyzz_1[i] * fe_0 + ta_yyyz_xxyyzz_0[i] * pa_x[i] - ta_yyyz_xxyyzz_1[i] * pc_x[i];

        ta_xyyyz_xxyzzz_0[i] =
            2.0 * ta_yyyz_xyzzz_0[i] * fe_0 - 2.0 * ta_yyyz_xyzzz_1[i] * fe_0 + ta_yyyz_xxyzzz_0[i] * pa_x[i] - ta_yyyz_xxyzzz_1[i] * pc_x[i];

        ta_xyyyz_xxzzzz_0[i] =
            2.0 * ta_yyyz_xzzzz_0[i] * fe_0 - 2.0 * ta_yyyz_xzzzz_1[i] * fe_0 + ta_yyyz_xxzzzz_0[i] * pa_x[i] - ta_yyyz_xxzzzz_1[i] * pc_x[i];

        ta_xyyyz_xyyyyy_0[i] = ta_xyyy_xyyyyy_0[i] * pa_z[i] - ta_xyyy_xyyyyy_1[i] * pc_z[i];

        ta_xyyyz_xyyyyz_0[i] = ta_yyyz_yyyyz_0[i] * fe_0 - ta_yyyz_yyyyz_1[i] * fe_0 + ta_yyyz_xyyyyz_0[i] * pa_x[i] - ta_yyyz_xyyyyz_1[i] * pc_x[i];

        ta_xyyyz_xyyyzz_0[i] = ta_yyyz_yyyzz_0[i] * fe_0 - ta_yyyz_yyyzz_1[i] * fe_0 + ta_yyyz_xyyyzz_0[i] * pa_x[i] - ta_yyyz_xyyyzz_1[i] * pc_x[i];

        ta_xyyyz_xyyzzz_0[i] = ta_yyyz_yyzzz_0[i] * fe_0 - ta_yyyz_yyzzz_1[i] * fe_0 + ta_yyyz_xyyzzz_0[i] * pa_x[i] - ta_yyyz_xyyzzz_1[i] * pc_x[i];

        ta_xyyyz_xyzzzz_0[i] = ta_yyyz_yzzzz_0[i] * fe_0 - ta_yyyz_yzzzz_1[i] * fe_0 + ta_yyyz_xyzzzz_0[i] * pa_x[i] - ta_yyyz_xyzzzz_1[i] * pc_x[i];

        ta_xyyyz_xzzzzz_0[i] = ta_yyyz_zzzzz_0[i] * fe_0 - ta_yyyz_zzzzz_1[i] * fe_0 + ta_yyyz_xzzzzz_0[i] * pa_x[i] - ta_yyyz_xzzzzz_1[i] * pc_x[i];

        ta_xyyyz_yyyyyy_0[i] = ta_yyyz_yyyyyy_0[i] * pa_x[i] - ta_yyyz_yyyyyy_1[i] * pc_x[i];

        ta_xyyyz_yyyyyz_0[i] = ta_yyyz_yyyyyz_0[i] * pa_x[i] - ta_yyyz_yyyyyz_1[i] * pc_x[i];

        ta_xyyyz_yyyyzz_0[i] = ta_yyyz_yyyyzz_0[i] * pa_x[i] - ta_yyyz_yyyyzz_1[i] * pc_x[i];

        ta_xyyyz_yyyzzz_0[i] = ta_yyyz_yyyzzz_0[i] * pa_x[i] - ta_yyyz_yyyzzz_1[i] * pc_x[i];

        ta_xyyyz_yyzzzz_0[i] = ta_yyyz_yyzzzz_0[i] * pa_x[i] - ta_yyyz_yyzzzz_1[i] * pc_x[i];

        ta_xyyyz_yzzzzz_0[i] = ta_yyyz_yzzzzz_0[i] * pa_x[i] - ta_yyyz_yzzzzz_1[i] * pc_x[i];

        ta_xyyyz_zzzzzz_0[i] = ta_yyyz_zzzzzz_0[i] * pa_x[i] - ta_yyyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 336-364 components of targeted buffer : HI

    auto ta_xyyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 336);

    auto ta_xyyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 337);

    auto ta_xyyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 338);

    auto ta_xyyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 339);

    auto ta_xyyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 340);

    auto ta_xyyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 341);

    auto ta_xyyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 342);

    auto ta_xyyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 343);

    auto ta_xyyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 344);

    auto ta_xyyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 345);

    auto ta_xyyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 346);

    auto ta_xyyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 347);

    auto ta_xyyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 348);

    auto ta_xyyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 349);

    auto ta_xyyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 350);

    auto ta_xyyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 351);

    auto ta_xyyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 352);

    auto ta_xyyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 353);

    auto ta_xyyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 354);

    auto ta_xyyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 355);

    auto ta_xyyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 356);

    auto ta_xyyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 357);

    auto ta_xyyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 358);

    auto ta_xyyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 359);

    auto ta_xyyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 360);

    auto ta_xyyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 361);

    auto ta_xyyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 362);

    auto ta_xyyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 363);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xyyzz_xxxxxx_0, \
                             ta_xyyzz_xxxxxy_0, \
                             ta_xyyzz_xxxxxz_0, \
                             ta_xyyzz_xxxxyy_0, \
                             ta_xyyzz_xxxxyz_0, \
                             ta_xyyzz_xxxxzz_0, \
                             ta_xyyzz_xxxyyy_0, \
                             ta_xyyzz_xxxyyz_0, \
                             ta_xyyzz_xxxyzz_0, \
                             ta_xyyzz_xxxzzz_0, \
                             ta_xyyzz_xxyyyy_0, \
                             ta_xyyzz_xxyyyz_0, \
                             ta_xyyzz_xxyyzz_0, \
                             ta_xyyzz_xxyzzz_0, \
                             ta_xyyzz_xxzzzz_0, \
                             ta_xyyzz_xyyyyy_0, \
                             ta_xyyzz_xyyyyz_0, \
                             ta_xyyzz_xyyyzz_0, \
                             ta_xyyzz_xyyzzz_0, \
                             ta_xyyzz_xyzzzz_0, \
                             ta_xyyzz_xzzzzz_0, \
                             ta_xyyzz_yyyyyy_0, \
                             ta_xyyzz_yyyyyz_0, \
                             ta_xyyzz_yyyyzz_0, \
                             ta_xyyzz_yyyzzz_0, \
                             ta_xyyzz_yyzzzz_0, \
                             ta_xyyzz_yzzzzz_0, \
                             ta_xyyzz_zzzzzz_0, \
                             ta_yyzz_xxxxx_0,   \
                             ta_yyzz_xxxxx_1,   \
                             ta_yyzz_xxxxxx_0,  \
                             ta_yyzz_xxxxxx_1,  \
                             ta_yyzz_xxxxxy_0,  \
                             ta_yyzz_xxxxxy_1,  \
                             ta_yyzz_xxxxxz_0,  \
                             ta_yyzz_xxxxxz_1,  \
                             ta_yyzz_xxxxy_0,   \
                             ta_yyzz_xxxxy_1,   \
                             ta_yyzz_xxxxyy_0,  \
                             ta_yyzz_xxxxyy_1,  \
                             ta_yyzz_xxxxyz_0,  \
                             ta_yyzz_xxxxyz_1,  \
                             ta_yyzz_xxxxz_0,   \
                             ta_yyzz_xxxxz_1,   \
                             ta_yyzz_xxxxzz_0,  \
                             ta_yyzz_xxxxzz_1,  \
                             ta_yyzz_xxxyy_0,   \
                             ta_yyzz_xxxyy_1,   \
                             ta_yyzz_xxxyyy_0,  \
                             ta_yyzz_xxxyyy_1,  \
                             ta_yyzz_xxxyyz_0,  \
                             ta_yyzz_xxxyyz_1,  \
                             ta_yyzz_xxxyz_0,   \
                             ta_yyzz_xxxyz_1,   \
                             ta_yyzz_xxxyzz_0,  \
                             ta_yyzz_xxxyzz_1,  \
                             ta_yyzz_xxxzz_0,   \
                             ta_yyzz_xxxzz_1,   \
                             ta_yyzz_xxxzzz_0,  \
                             ta_yyzz_xxxzzz_1,  \
                             ta_yyzz_xxyyy_0,   \
                             ta_yyzz_xxyyy_1,   \
                             ta_yyzz_xxyyyy_0,  \
                             ta_yyzz_xxyyyy_1,  \
                             ta_yyzz_xxyyyz_0,  \
                             ta_yyzz_xxyyyz_1,  \
                             ta_yyzz_xxyyz_0,   \
                             ta_yyzz_xxyyz_1,   \
                             ta_yyzz_xxyyzz_0,  \
                             ta_yyzz_xxyyzz_1,  \
                             ta_yyzz_xxyzz_0,   \
                             ta_yyzz_xxyzz_1,   \
                             ta_yyzz_xxyzzz_0,  \
                             ta_yyzz_xxyzzz_1,  \
                             ta_yyzz_xxzzz_0,   \
                             ta_yyzz_xxzzz_1,   \
                             ta_yyzz_xxzzzz_0,  \
                             ta_yyzz_xxzzzz_1,  \
                             ta_yyzz_xyyyy_0,   \
                             ta_yyzz_xyyyy_1,   \
                             ta_yyzz_xyyyyy_0,  \
                             ta_yyzz_xyyyyy_1,  \
                             ta_yyzz_xyyyyz_0,  \
                             ta_yyzz_xyyyyz_1,  \
                             ta_yyzz_xyyyz_0,   \
                             ta_yyzz_xyyyz_1,   \
                             ta_yyzz_xyyyzz_0,  \
                             ta_yyzz_xyyyzz_1,  \
                             ta_yyzz_xyyzz_0,   \
                             ta_yyzz_xyyzz_1,   \
                             ta_yyzz_xyyzzz_0,  \
                             ta_yyzz_xyyzzz_1,  \
                             ta_yyzz_xyzzz_0,   \
                             ta_yyzz_xyzzz_1,   \
                             ta_yyzz_xyzzzz_0,  \
                             ta_yyzz_xyzzzz_1,  \
                             ta_yyzz_xzzzz_0,   \
                             ta_yyzz_xzzzz_1,   \
                             ta_yyzz_xzzzzz_0,  \
                             ta_yyzz_xzzzzz_1,  \
                             ta_yyzz_yyyyy_0,   \
                             ta_yyzz_yyyyy_1,   \
                             ta_yyzz_yyyyyy_0,  \
                             ta_yyzz_yyyyyy_1,  \
                             ta_yyzz_yyyyyz_0,  \
                             ta_yyzz_yyyyyz_1,  \
                             ta_yyzz_yyyyz_0,   \
                             ta_yyzz_yyyyz_1,   \
                             ta_yyzz_yyyyzz_0,  \
                             ta_yyzz_yyyyzz_1,  \
                             ta_yyzz_yyyzz_0,   \
                             ta_yyzz_yyyzz_1,   \
                             ta_yyzz_yyyzzz_0,  \
                             ta_yyzz_yyyzzz_1,  \
                             ta_yyzz_yyzzz_0,   \
                             ta_yyzz_yyzzz_1,   \
                             ta_yyzz_yyzzzz_0,  \
                             ta_yyzz_yyzzzz_1,  \
                             ta_yyzz_yzzzz_0,   \
                             ta_yyzz_yzzzz_1,   \
                             ta_yyzz_yzzzzz_0,  \
                             ta_yyzz_yzzzzz_1,  \
                             ta_yyzz_zzzzz_0,   \
                             ta_yyzz_zzzzz_1,   \
                             ta_yyzz_zzzzzz_0,  \
                             ta_yyzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzz_xxxxxx_0[i] =
            6.0 * ta_yyzz_xxxxx_0[i] * fe_0 - 6.0 * ta_yyzz_xxxxx_1[i] * fe_0 + ta_yyzz_xxxxxx_0[i] * pa_x[i] - ta_yyzz_xxxxxx_1[i] * pc_x[i];

        ta_xyyzz_xxxxxy_0[i] =
            5.0 * ta_yyzz_xxxxy_0[i] * fe_0 - 5.0 * ta_yyzz_xxxxy_1[i] * fe_0 + ta_yyzz_xxxxxy_0[i] * pa_x[i] - ta_yyzz_xxxxxy_1[i] * pc_x[i];

        ta_xyyzz_xxxxxz_0[i] =
            5.0 * ta_yyzz_xxxxz_0[i] * fe_0 - 5.0 * ta_yyzz_xxxxz_1[i] * fe_0 + ta_yyzz_xxxxxz_0[i] * pa_x[i] - ta_yyzz_xxxxxz_1[i] * pc_x[i];

        ta_xyyzz_xxxxyy_0[i] =
            4.0 * ta_yyzz_xxxyy_0[i] * fe_0 - 4.0 * ta_yyzz_xxxyy_1[i] * fe_0 + ta_yyzz_xxxxyy_0[i] * pa_x[i] - ta_yyzz_xxxxyy_1[i] * pc_x[i];

        ta_xyyzz_xxxxyz_0[i] =
            4.0 * ta_yyzz_xxxyz_0[i] * fe_0 - 4.0 * ta_yyzz_xxxyz_1[i] * fe_0 + ta_yyzz_xxxxyz_0[i] * pa_x[i] - ta_yyzz_xxxxyz_1[i] * pc_x[i];

        ta_xyyzz_xxxxzz_0[i] =
            4.0 * ta_yyzz_xxxzz_0[i] * fe_0 - 4.0 * ta_yyzz_xxxzz_1[i] * fe_0 + ta_yyzz_xxxxzz_0[i] * pa_x[i] - ta_yyzz_xxxxzz_1[i] * pc_x[i];

        ta_xyyzz_xxxyyy_0[i] =
            3.0 * ta_yyzz_xxyyy_0[i] * fe_0 - 3.0 * ta_yyzz_xxyyy_1[i] * fe_0 + ta_yyzz_xxxyyy_0[i] * pa_x[i] - ta_yyzz_xxxyyy_1[i] * pc_x[i];

        ta_xyyzz_xxxyyz_0[i] =
            3.0 * ta_yyzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyyz_1[i] * fe_0 + ta_yyzz_xxxyyz_0[i] * pa_x[i] - ta_yyzz_xxxyyz_1[i] * pc_x[i];

        ta_xyyzz_xxxyzz_0[i] =
            3.0 * ta_yyzz_xxyzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxyzz_1[i] * fe_0 + ta_yyzz_xxxyzz_0[i] * pa_x[i] - ta_yyzz_xxxyzz_1[i] * pc_x[i];

        ta_xyyzz_xxxzzz_0[i] =
            3.0 * ta_yyzz_xxzzz_0[i] * fe_0 - 3.0 * ta_yyzz_xxzzz_1[i] * fe_0 + ta_yyzz_xxxzzz_0[i] * pa_x[i] - ta_yyzz_xxxzzz_1[i] * pc_x[i];

        ta_xyyzz_xxyyyy_0[i] =
            2.0 * ta_yyzz_xyyyy_0[i] * fe_0 - 2.0 * ta_yyzz_xyyyy_1[i] * fe_0 + ta_yyzz_xxyyyy_0[i] * pa_x[i] - ta_yyzz_xxyyyy_1[i] * pc_x[i];

        ta_xyyzz_xxyyyz_0[i] =
            2.0 * ta_yyzz_xyyyz_0[i] * fe_0 - 2.0 * ta_yyzz_xyyyz_1[i] * fe_0 + ta_yyzz_xxyyyz_0[i] * pa_x[i] - ta_yyzz_xxyyyz_1[i] * pc_x[i];

        ta_xyyzz_xxyyzz_0[i] =
            2.0 * ta_yyzz_xyyzz_0[i] * fe_0 - 2.0 * ta_yyzz_xyyzz_1[i] * fe_0 + ta_yyzz_xxyyzz_0[i] * pa_x[i] - ta_yyzz_xxyyzz_1[i] * pc_x[i];

        ta_xyyzz_xxyzzz_0[i] =
            2.0 * ta_yyzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yyzz_xyzzz_1[i] * fe_0 + ta_yyzz_xxyzzz_0[i] * pa_x[i] - ta_yyzz_xxyzzz_1[i] * pc_x[i];

        ta_xyyzz_xxzzzz_0[i] =
            2.0 * ta_yyzz_xzzzz_0[i] * fe_0 - 2.0 * ta_yyzz_xzzzz_1[i] * fe_0 + ta_yyzz_xxzzzz_0[i] * pa_x[i] - ta_yyzz_xxzzzz_1[i] * pc_x[i];

        ta_xyyzz_xyyyyy_0[i] = ta_yyzz_yyyyy_0[i] * fe_0 - ta_yyzz_yyyyy_1[i] * fe_0 + ta_yyzz_xyyyyy_0[i] * pa_x[i] - ta_yyzz_xyyyyy_1[i] * pc_x[i];

        ta_xyyzz_xyyyyz_0[i] = ta_yyzz_yyyyz_0[i] * fe_0 - ta_yyzz_yyyyz_1[i] * fe_0 + ta_yyzz_xyyyyz_0[i] * pa_x[i] - ta_yyzz_xyyyyz_1[i] * pc_x[i];

        ta_xyyzz_xyyyzz_0[i] = ta_yyzz_yyyzz_0[i] * fe_0 - ta_yyzz_yyyzz_1[i] * fe_0 + ta_yyzz_xyyyzz_0[i] * pa_x[i] - ta_yyzz_xyyyzz_1[i] * pc_x[i];

        ta_xyyzz_xyyzzz_0[i] = ta_yyzz_yyzzz_0[i] * fe_0 - ta_yyzz_yyzzz_1[i] * fe_0 + ta_yyzz_xyyzzz_0[i] * pa_x[i] - ta_yyzz_xyyzzz_1[i] * pc_x[i];

        ta_xyyzz_xyzzzz_0[i] = ta_yyzz_yzzzz_0[i] * fe_0 - ta_yyzz_yzzzz_1[i] * fe_0 + ta_yyzz_xyzzzz_0[i] * pa_x[i] - ta_yyzz_xyzzzz_1[i] * pc_x[i];

        ta_xyyzz_xzzzzz_0[i] = ta_yyzz_zzzzz_0[i] * fe_0 - ta_yyzz_zzzzz_1[i] * fe_0 + ta_yyzz_xzzzzz_0[i] * pa_x[i] - ta_yyzz_xzzzzz_1[i] * pc_x[i];

        ta_xyyzz_yyyyyy_0[i] = ta_yyzz_yyyyyy_0[i] * pa_x[i] - ta_yyzz_yyyyyy_1[i] * pc_x[i];

        ta_xyyzz_yyyyyz_0[i] = ta_yyzz_yyyyyz_0[i] * pa_x[i] - ta_yyzz_yyyyyz_1[i] * pc_x[i];

        ta_xyyzz_yyyyzz_0[i] = ta_yyzz_yyyyzz_0[i] * pa_x[i] - ta_yyzz_yyyyzz_1[i] * pc_x[i];

        ta_xyyzz_yyyzzz_0[i] = ta_yyzz_yyyzzz_0[i] * pa_x[i] - ta_yyzz_yyyzzz_1[i] * pc_x[i];

        ta_xyyzz_yyzzzz_0[i] = ta_yyzz_yyzzzz_0[i] * pa_x[i] - ta_yyzz_yyzzzz_1[i] * pc_x[i];

        ta_xyyzz_yzzzzz_0[i] = ta_yyzz_yzzzzz_0[i] * pa_x[i] - ta_yyzz_yzzzzz_1[i] * pc_x[i];

        ta_xyyzz_zzzzzz_0[i] = ta_yyzz_zzzzzz_0[i] * pa_x[i] - ta_yyzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 364-392 components of targeted buffer : HI

    auto ta_xyzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 364);

    auto ta_xyzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 365);

    auto ta_xyzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 366);

    auto ta_xyzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 367);

    auto ta_xyzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 368);

    auto ta_xyzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 369);

    auto ta_xyzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 370);

    auto ta_xyzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 371);

    auto ta_xyzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 372);

    auto ta_xyzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 373);

    auto ta_xyzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 374);

    auto ta_xyzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 375);

    auto ta_xyzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 376);

    auto ta_xyzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 377);

    auto ta_xyzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 378);

    auto ta_xyzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 379);

    auto ta_xyzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 380);

    auto ta_xyzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 381);

    auto ta_xyzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 382);

    auto ta_xyzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 383);

    auto ta_xyzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 384);

    auto ta_xyzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 385);

    auto ta_xyzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 386);

    auto ta_xyzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 387);

    auto ta_xyzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 388);

    auto ta_xyzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 389);

    auto ta_xyzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 390);

    auto ta_xyzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 391);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta_xyzzz_xxxxxx_0, \
                             ta_xyzzz_xxxxxy_0, \
                             ta_xyzzz_xxxxxz_0, \
                             ta_xyzzz_xxxxyy_0, \
                             ta_xyzzz_xxxxyz_0, \
                             ta_xyzzz_xxxxzz_0, \
                             ta_xyzzz_xxxyyy_0, \
                             ta_xyzzz_xxxyyz_0, \
                             ta_xyzzz_xxxyzz_0, \
                             ta_xyzzz_xxxzzz_0, \
                             ta_xyzzz_xxyyyy_0, \
                             ta_xyzzz_xxyyyz_0, \
                             ta_xyzzz_xxyyzz_0, \
                             ta_xyzzz_xxyzzz_0, \
                             ta_xyzzz_xxzzzz_0, \
                             ta_xyzzz_xyyyyy_0, \
                             ta_xyzzz_xyyyyz_0, \
                             ta_xyzzz_xyyyzz_0, \
                             ta_xyzzz_xyyzzz_0, \
                             ta_xyzzz_xyzzzz_0, \
                             ta_xyzzz_xzzzzz_0, \
                             ta_xyzzz_yyyyyy_0, \
                             ta_xyzzz_yyyyyz_0, \
                             ta_xyzzz_yyyyzz_0, \
                             ta_xyzzz_yyyzzz_0, \
                             ta_xyzzz_yyzzzz_0, \
                             ta_xyzzz_yzzzzz_0, \
                             ta_xyzzz_zzzzzz_0, \
                             ta_xzzz_xxxxxx_0,  \
                             ta_xzzz_xxxxxx_1,  \
                             ta_xzzz_xxxxxz_0,  \
                             ta_xzzz_xxxxxz_1,  \
                             ta_xzzz_xxxxzz_0,  \
                             ta_xzzz_xxxxzz_1,  \
                             ta_xzzz_xxxzzz_0,  \
                             ta_xzzz_xxxzzz_1,  \
                             ta_xzzz_xxzzzz_0,  \
                             ta_xzzz_xxzzzz_1,  \
                             ta_xzzz_xzzzzz_0,  \
                             ta_xzzz_xzzzzz_1,  \
                             ta_yzzz_xxxxxy_0,  \
                             ta_yzzz_xxxxxy_1,  \
                             ta_yzzz_xxxxy_0,   \
                             ta_yzzz_xxxxy_1,   \
                             ta_yzzz_xxxxyy_0,  \
                             ta_yzzz_xxxxyy_1,  \
                             ta_yzzz_xxxxyz_0,  \
                             ta_yzzz_xxxxyz_1,  \
                             ta_yzzz_xxxyy_0,   \
                             ta_yzzz_xxxyy_1,   \
                             ta_yzzz_xxxyyy_0,  \
                             ta_yzzz_xxxyyy_1,  \
                             ta_yzzz_xxxyyz_0,  \
                             ta_yzzz_xxxyyz_1,  \
                             ta_yzzz_xxxyz_0,   \
                             ta_yzzz_xxxyz_1,   \
                             ta_yzzz_xxxyzz_0,  \
                             ta_yzzz_xxxyzz_1,  \
                             ta_yzzz_xxyyy_0,   \
                             ta_yzzz_xxyyy_1,   \
                             ta_yzzz_xxyyyy_0,  \
                             ta_yzzz_xxyyyy_1,  \
                             ta_yzzz_xxyyyz_0,  \
                             ta_yzzz_xxyyyz_1,  \
                             ta_yzzz_xxyyz_0,   \
                             ta_yzzz_xxyyz_1,   \
                             ta_yzzz_xxyyzz_0,  \
                             ta_yzzz_xxyyzz_1,  \
                             ta_yzzz_xxyzz_0,   \
                             ta_yzzz_xxyzz_1,   \
                             ta_yzzz_xxyzzz_0,  \
                             ta_yzzz_xxyzzz_1,  \
                             ta_yzzz_xyyyy_0,   \
                             ta_yzzz_xyyyy_1,   \
                             ta_yzzz_xyyyyy_0,  \
                             ta_yzzz_xyyyyy_1,  \
                             ta_yzzz_xyyyyz_0,  \
                             ta_yzzz_xyyyyz_1,  \
                             ta_yzzz_xyyyz_0,   \
                             ta_yzzz_xyyyz_1,   \
                             ta_yzzz_xyyyzz_0,  \
                             ta_yzzz_xyyyzz_1,  \
                             ta_yzzz_xyyzz_0,   \
                             ta_yzzz_xyyzz_1,   \
                             ta_yzzz_xyyzzz_0,  \
                             ta_yzzz_xyyzzz_1,  \
                             ta_yzzz_xyzzz_0,   \
                             ta_yzzz_xyzzz_1,   \
                             ta_yzzz_xyzzzz_0,  \
                             ta_yzzz_xyzzzz_1,  \
                             ta_yzzz_yyyyy_0,   \
                             ta_yzzz_yyyyy_1,   \
                             ta_yzzz_yyyyyy_0,  \
                             ta_yzzz_yyyyyy_1,  \
                             ta_yzzz_yyyyyz_0,  \
                             ta_yzzz_yyyyyz_1,  \
                             ta_yzzz_yyyyz_0,   \
                             ta_yzzz_yyyyz_1,   \
                             ta_yzzz_yyyyzz_0,  \
                             ta_yzzz_yyyyzz_1,  \
                             ta_yzzz_yyyzz_0,   \
                             ta_yzzz_yyyzz_1,   \
                             ta_yzzz_yyyzzz_0,  \
                             ta_yzzz_yyyzzz_1,  \
                             ta_yzzz_yyzzz_0,   \
                             ta_yzzz_yyzzz_1,   \
                             ta_yzzz_yyzzzz_0,  \
                             ta_yzzz_yyzzzz_1,  \
                             ta_yzzz_yzzzz_0,   \
                             ta_yzzz_yzzzz_1,   \
                             ta_yzzz_yzzzzz_0,  \
                             ta_yzzz_yzzzzz_1,  \
                             ta_yzzz_zzzzzz_0,  \
                             ta_yzzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzz_xxxxxx_0[i] = ta_xzzz_xxxxxx_0[i] * pa_y[i] - ta_xzzz_xxxxxx_1[i] * pc_y[i];

        ta_xyzzz_xxxxxy_0[i] =
            5.0 * ta_yzzz_xxxxy_0[i] * fe_0 - 5.0 * ta_yzzz_xxxxy_1[i] * fe_0 + ta_yzzz_xxxxxy_0[i] * pa_x[i] - ta_yzzz_xxxxxy_1[i] * pc_x[i];

        ta_xyzzz_xxxxxz_0[i] = ta_xzzz_xxxxxz_0[i] * pa_y[i] - ta_xzzz_xxxxxz_1[i] * pc_y[i];

        ta_xyzzz_xxxxyy_0[i] =
            4.0 * ta_yzzz_xxxyy_0[i] * fe_0 - 4.0 * ta_yzzz_xxxyy_1[i] * fe_0 + ta_yzzz_xxxxyy_0[i] * pa_x[i] - ta_yzzz_xxxxyy_1[i] * pc_x[i];

        ta_xyzzz_xxxxyz_0[i] =
            4.0 * ta_yzzz_xxxyz_0[i] * fe_0 - 4.0 * ta_yzzz_xxxyz_1[i] * fe_0 + ta_yzzz_xxxxyz_0[i] * pa_x[i] - ta_yzzz_xxxxyz_1[i] * pc_x[i];

        ta_xyzzz_xxxxzz_0[i] = ta_xzzz_xxxxzz_0[i] * pa_y[i] - ta_xzzz_xxxxzz_1[i] * pc_y[i];

        ta_xyzzz_xxxyyy_0[i] =
            3.0 * ta_yzzz_xxyyy_0[i] * fe_0 - 3.0 * ta_yzzz_xxyyy_1[i] * fe_0 + ta_yzzz_xxxyyy_0[i] * pa_x[i] - ta_yzzz_xxxyyy_1[i] * pc_x[i];

        ta_xyzzz_xxxyyz_0[i] =
            3.0 * ta_yzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yzzz_xxyyz_1[i] * fe_0 + ta_yzzz_xxxyyz_0[i] * pa_x[i] - ta_yzzz_xxxyyz_1[i] * pc_x[i];

        ta_xyzzz_xxxyzz_0[i] =
            3.0 * ta_yzzz_xxyzz_0[i] * fe_0 - 3.0 * ta_yzzz_xxyzz_1[i] * fe_0 + ta_yzzz_xxxyzz_0[i] * pa_x[i] - ta_yzzz_xxxyzz_1[i] * pc_x[i];

        ta_xyzzz_xxxzzz_0[i] = ta_xzzz_xxxzzz_0[i] * pa_y[i] - ta_xzzz_xxxzzz_1[i] * pc_y[i];

        ta_xyzzz_xxyyyy_0[i] =
            2.0 * ta_yzzz_xyyyy_0[i] * fe_0 - 2.0 * ta_yzzz_xyyyy_1[i] * fe_0 + ta_yzzz_xxyyyy_0[i] * pa_x[i] - ta_yzzz_xxyyyy_1[i] * pc_x[i];

        ta_xyzzz_xxyyyz_0[i] =
            2.0 * ta_yzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyyz_1[i] * fe_0 + ta_yzzz_xxyyyz_0[i] * pa_x[i] - ta_yzzz_xxyyyz_1[i] * pc_x[i];

        ta_xyzzz_xxyyzz_0[i] =
            2.0 * ta_yzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyyzz_1[i] * fe_0 + ta_yzzz_xxyyzz_0[i] * pa_x[i] - ta_yzzz_xxyyzz_1[i] * pc_x[i];

        ta_xyzzz_xxyzzz_0[i] =
            2.0 * ta_yzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yzzz_xyzzz_1[i] * fe_0 + ta_yzzz_xxyzzz_0[i] * pa_x[i] - ta_yzzz_xxyzzz_1[i] * pc_x[i];

        ta_xyzzz_xxzzzz_0[i] = ta_xzzz_xxzzzz_0[i] * pa_y[i] - ta_xzzz_xxzzzz_1[i] * pc_y[i];

        ta_xyzzz_xyyyyy_0[i] = ta_yzzz_yyyyy_0[i] * fe_0 - ta_yzzz_yyyyy_1[i] * fe_0 + ta_yzzz_xyyyyy_0[i] * pa_x[i] - ta_yzzz_xyyyyy_1[i] * pc_x[i];

        ta_xyzzz_xyyyyz_0[i] = ta_yzzz_yyyyz_0[i] * fe_0 - ta_yzzz_yyyyz_1[i] * fe_0 + ta_yzzz_xyyyyz_0[i] * pa_x[i] - ta_yzzz_xyyyyz_1[i] * pc_x[i];

        ta_xyzzz_xyyyzz_0[i] = ta_yzzz_yyyzz_0[i] * fe_0 - ta_yzzz_yyyzz_1[i] * fe_0 + ta_yzzz_xyyyzz_0[i] * pa_x[i] - ta_yzzz_xyyyzz_1[i] * pc_x[i];

        ta_xyzzz_xyyzzz_0[i] = ta_yzzz_yyzzz_0[i] * fe_0 - ta_yzzz_yyzzz_1[i] * fe_0 + ta_yzzz_xyyzzz_0[i] * pa_x[i] - ta_yzzz_xyyzzz_1[i] * pc_x[i];

        ta_xyzzz_xyzzzz_0[i] = ta_yzzz_yzzzz_0[i] * fe_0 - ta_yzzz_yzzzz_1[i] * fe_0 + ta_yzzz_xyzzzz_0[i] * pa_x[i] - ta_yzzz_xyzzzz_1[i] * pc_x[i];

        ta_xyzzz_xzzzzz_0[i] = ta_xzzz_xzzzzz_0[i] * pa_y[i] - ta_xzzz_xzzzzz_1[i] * pc_y[i];

        ta_xyzzz_yyyyyy_0[i] = ta_yzzz_yyyyyy_0[i] * pa_x[i] - ta_yzzz_yyyyyy_1[i] * pc_x[i];

        ta_xyzzz_yyyyyz_0[i] = ta_yzzz_yyyyyz_0[i] * pa_x[i] - ta_yzzz_yyyyyz_1[i] * pc_x[i];

        ta_xyzzz_yyyyzz_0[i] = ta_yzzz_yyyyzz_0[i] * pa_x[i] - ta_yzzz_yyyyzz_1[i] * pc_x[i];

        ta_xyzzz_yyyzzz_0[i] = ta_yzzz_yyyzzz_0[i] * pa_x[i] - ta_yzzz_yyyzzz_1[i] * pc_x[i];

        ta_xyzzz_yyzzzz_0[i] = ta_yzzz_yyzzzz_0[i] * pa_x[i] - ta_yzzz_yyzzzz_1[i] * pc_x[i];

        ta_xyzzz_yzzzzz_0[i] = ta_yzzz_yzzzzz_0[i] * pa_x[i] - ta_yzzz_yzzzzz_1[i] * pc_x[i];

        ta_xyzzz_zzzzzz_0[i] = ta_yzzz_zzzzzz_0[i] * pa_x[i] - ta_yzzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 392-420 components of targeted buffer : HI

    auto ta_xzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 392);

    auto ta_xzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 393);

    auto ta_xzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 394);

    auto ta_xzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 395);

    auto ta_xzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 396);

    auto ta_xzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 397);

    auto ta_xzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 398);

    auto ta_xzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 399);

    auto ta_xzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 400);

    auto ta_xzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 401);

    auto ta_xzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 402);

    auto ta_xzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 403);

    auto ta_xzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 404);

    auto ta_xzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 405);

    auto ta_xzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 406);

    auto ta_xzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 407);

    auto ta_xzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 408);

    auto ta_xzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 409);

    auto ta_xzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 410);

    auto ta_xzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 411);

    auto ta_xzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 412);

    auto ta_xzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 413);

    auto ta_xzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 414);

    auto ta_xzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 415);

    auto ta_xzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 416);

    auto ta_xzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 417);

    auto ta_xzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 418);

    auto ta_xzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 419);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta_xzzzz_xxxxxx_0, \
                             ta_xzzzz_xxxxxy_0, \
                             ta_xzzzz_xxxxxz_0, \
                             ta_xzzzz_xxxxyy_0, \
                             ta_xzzzz_xxxxyz_0, \
                             ta_xzzzz_xxxxzz_0, \
                             ta_xzzzz_xxxyyy_0, \
                             ta_xzzzz_xxxyyz_0, \
                             ta_xzzzz_xxxyzz_0, \
                             ta_xzzzz_xxxzzz_0, \
                             ta_xzzzz_xxyyyy_0, \
                             ta_xzzzz_xxyyyz_0, \
                             ta_xzzzz_xxyyzz_0, \
                             ta_xzzzz_xxyzzz_0, \
                             ta_xzzzz_xxzzzz_0, \
                             ta_xzzzz_xyyyyy_0, \
                             ta_xzzzz_xyyyyz_0, \
                             ta_xzzzz_xyyyzz_0, \
                             ta_xzzzz_xyyzzz_0, \
                             ta_xzzzz_xyzzzz_0, \
                             ta_xzzzz_xzzzzz_0, \
                             ta_xzzzz_yyyyyy_0, \
                             ta_xzzzz_yyyyyz_0, \
                             ta_xzzzz_yyyyzz_0, \
                             ta_xzzzz_yyyzzz_0, \
                             ta_xzzzz_yyzzzz_0, \
                             ta_xzzzz_yzzzzz_0, \
                             ta_xzzzz_zzzzzz_0, \
                             ta_zzzz_xxxxx_0,   \
                             ta_zzzz_xxxxx_1,   \
                             ta_zzzz_xxxxxx_0,  \
                             ta_zzzz_xxxxxx_1,  \
                             ta_zzzz_xxxxxy_0,  \
                             ta_zzzz_xxxxxy_1,  \
                             ta_zzzz_xxxxxz_0,  \
                             ta_zzzz_xxxxxz_1,  \
                             ta_zzzz_xxxxy_0,   \
                             ta_zzzz_xxxxy_1,   \
                             ta_zzzz_xxxxyy_0,  \
                             ta_zzzz_xxxxyy_1,  \
                             ta_zzzz_xxxxyz_0,  \
                             ta_zzzz_xxxxyz_1,  \
                             ta_zzzz_xxxxz_0,   \
                             ta_zzzz_xxxxz_1,   \
                             ta_zzzz_xxxxzz_0,  \
                             ta_zzzz_xxxxzz_1,  \
                             ta_zzzz_xxxyy_0,   \
                             ta_zzzz_xxxyy_1,   \
                             ta_zzzz_xxxyyy_0,  \
                             ta_zzzz_xxxyyy_1,  \
                             ta_zzzz_xxxyyz_0,  \
                             ta_zzzz_xxxyyz_1,  \
                             ta_zzzz_xxxyz_0,   \
                             ta_zzzz_xxxyz_1,   \
                             ta_zzzz_xxxyzz_0,  \
                             ta_zzzz_xxxyzz_1,  \
                             ta_zzzz_xxxzz_0,   \
                             ta_zzzz_xxxzz_1,   \
                             ta_zzzz_xxxzzz_0,  \
                             ta_zzzz_xxxzzz_1,  \
                             ta_zzzz_xxyyy_0,   \
                             ta_zzzz_xxyyy_1,   \
                             ta_zzzz_xxyyyy_0,  \
                             ta_zzzz_xxyyyy_1,  \
                             ta_zzzz_xxyyyz_0,  \
                             ta_zzzz_xxyyyz_1,  \
                             ta_zzzz_xxyyz_0,   \
                             ta_zzzz_xxyyz_1,   \
                             ta_zzzz_xxyyzz_0,  \
                             ta_zzzz_xxyyzz_1,  \
                             ta_zzzz_xxyzz_0,   \
                             ta_zzzz_xxyzz_1,   \
                             ta_zzzz_xxyzzz_0,  \
                             ta_zzzz_xxyzzz_1,  \
                             ta_zzzz_xxzzz_0,   \
                             ta_zzzz_xxzzz_1,   \
                             ta_zzzz_xxzzzz_0,  \
                             ta_zzzz_xxzzzz_1,  \
                             ta_zzzz_xyyyy_0,   \
                             ta_zzzz_xyyyy_1,   \
                             ta_zzzz_xyyyyy_0,  \
                             ta_zzzz_xyyyyy_1,  \
                             ta_zzzz_xyyyyz_0,  \
                             ta_zzzz_xyyyyz_1,  \
                             ta_zzzz_xyyyz_0,   \
                             ta_zzzz_xyyyz_1,   \
                             ta_zzzz_xyyyzz_0,  \
                             ta_zzzz_xyyyzz_1,  \
                             ta_zzzz_xyyzz_0,   \
                             ta_zzzz_xyyzz_1,   \
                             ta_zzzz_xyyzzz_0,  \
                             ta_zzzz_xyyzzz_1,  \
                             ta_zzzz_xyzzz_0,   \
                             ta_zzzz_xyzzz_1,   \
                             ta_zzzz_xyzzzz_0,  \
                             ta_zzzz_xyzzzz_1,  \
                             ta_zzzz_xzzzz_0,   \
                             ta_zzzz_xzzzz_1,   \
                             ta_zzzz_xzzzzz_0,  \
                             ta_zzzz_xzzzzz_1,  \
                             ta_zzzz_yyyyy_0,   \
                             ta_zzzz_yyyyy_1,   \
                             ta_zzzz_yyyyyy_0,  \
                             ta_zzzz_yyyyyy_1,  \
                             ta_zzzz_yyyyyz_0,  \
                             ta_zzzz_yyyyyz_1,  \
                             ta_zzzz_yyyyz_0,   \
                             ta_zzzz_yyyyz_1,   \
                             ta_zzzz_yyyyzz_0,  \
                             ta_zzzz_yyyyzz_1,  \
                             ta_zzzz_yyyzz_0,   \
                             ta_zzzz_yyyzz_1,   \
                             ta_zzzz_yyyzzz_0,  \
                             ta_zzzz_yyyzzz_1,  \
                             ta_zzzz_yyzzz_0,   \
                             ta_zzzz_yyzzz_1,   \
                             ta_zzzz_yyzzzz_0,  \
                             ta_zzzz_yyzzzz_1,  \
                             ta_zzzz_yzzzz_0,   \
                             ta_zzzz_yzzzz_1,   \
                             ta_zzzz_yzzzzz_0,  \
                             ta_zzzz_yzzzzz_1,  \
                             ta_zzzz_zzzzz_0,   \
                             ta_zzzz_zzzzz_1,   \
                             ta_zzzz_zzzzzz_0,  \
                             ta_zzzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzz_xxxxxx_0[i] =
            6.0 * ta_zzzz_xxxxx_0[i] * fe_0 - 6.0 * ta_zzzz_xxxxx_1[i] * fe_0 + ta_zzzz_xxxxxx_0[i] * pa_x[i] - ta_zzzz_xxxxxx_1[i] * pc_x[i];

        ta_xzzzz_xxxxxy_0[i] =
            5.0 * ta_zzzz_xxxxy_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxy_1[i] * fe_0 + ta_zzzz_xxxxxy_0[i] * pa_x[i] - ta_zzzz_xxxxxy_1[i] * pc_x[i];

        ta_xzzzz_xxxxxz_0[i] =
            5.0 * ta_zzzz_xxxxz_0[i] * fe_0 - 5.0 * ta_zzzz_xxxxz_1[i] * fe_0 + ta_zzzz_xxxxxz_0[i] * pa_x[i] - ta_zzzz_xxxxxz_1[i] * pc_x[i];

        ta_xzzzz_xxxxyy_0[i] =
            4.0 * ta_zzzz_xxxyy_0[i] * fe_0 - 4.0 * ta_zzzz_xxxyy_1[i] * fe_0 + ta_zzzz_xxxxyy_0[i] * pa_x[i] - ta_zzzz_xxxxyy_1[i] * pc_x[i];

        ta_xzzzz_xxxxyz_0[i] =
            4.0 * ta_zzzz_xxxyz_0[i] * fe_0 - 4.0 * ta_zzzz_xxxyz_1[i] * fe_0 + ta_zzzz_xxxxyz_0[i] * pa_x[i] - ta_zzzz_xxxxyz_1[i] * pc_x[i];

        ta_xzzzz_xxxxzz_0[i] =
            4.0 * ta_zzzz_xxxzz_0[i] * fe_0 - 4.0 * ta_zzzz_xxxzz_1[i] * fe_0 + ta_zzzz_xxxxzz_0[i] * pa_x[i] - ta_zzzz_xxxxzz_1[i] * pc_x[i];

        ta_xzzzz_xxxyyy_0[i] =
            3.0 * ta_zzzz_xxyyy_0[i] * fe_0 - 3.0 * ta_zzzz_xxyyy_1[i] * fe_0 + ta_zzzz_xxxyyy_0[i] * pa_x[i] - ta_zzzz_xxxyyy_1[i] * pc_x[i];

        ta_xzzzz_xxxyyz_0[i] =
            3.0 * ta_zzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_zzzz_xxyyz_1[i] * fe_0 + ta_zzzz_xxxyyz_0[i] * pa_x[i] - ta_zzzz_xxxyyz_1[i] * pc_x[i];

        ta_xzzzz_xxxyzz_0[i] =
            3.0 * ta_zzzz_xxyzz_0[i] * fe_0 - 3.0 * ta_zzzz_xxyzz_1[i] * fe_0 + ta_zzzz_xxxyzz_0[i] * pa_x[i] - ta_zzzz_xxxyzz_1[i] * pc_x[i];

        ta_xzzzz_xxxzzz_0[i] =
            3.0 * ta_zzzz_xxzzz_0[i] * fe_0 - 3.0 * ta_zzzz_xxzzz_1[i] * fe_0 + ta_zzzz_xxxzzz_0[i] * pa_x[i] - ta_zzzz_xxxzzz_1[i] * pc_x[i];

        ta_xzzzz_xxyyyy_0[i] =
            2.0 * ta_zzzz_xyyyy_0[i] * fe_0 - 2.0 * ta_zzzz_xyyyy_1[i] * fe_0 + ta_zzzz_xxyyyy_0[i] * pa_x[i] - ta_zzzz_xxyyyy_1[i] * pc_x[i];

        ta_xzzzz_xxyyyz_0[i] =
            2.0 * ta_zzzz_xyyyz_0[i] * fe_0 - 2.0 * ta_zzzz_xyyyz_1[i] * fe_0 + ta_zzzz_xxyyyz_0[i] * pa_x[i] - ta_zzzz_xxyyyz_1[i] * pc_x[i];

        ta_xzzzz_xxyyzz_0[i] =
            2.0 * ta_zzzz_xyyzz_0[i] * fe_0 - 2.0 * ta_zzzz_xyyzz_1[i] * fe_0 + ta_zzzz_xxyyzz_0[i] * pa_x[i] - ta_zzzz_xxyyzz_1[i] * pc_x[i];

        ta_xzzzz_xxyzzz_0[i] =
            2.0 * ta_zzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_zzzz_xyzzz_1[i] * fe_0 + ta_zzzz_xxyzzz_0[i] * pa_x[i] - ta_zzzz_xxyzzz_1[i] * pc_x[i];

        ta_xzzzz_xxzzzz_0[i] =
            2.0 * ta_zzzz_xzzzz_0[i] * fe_0 - 2.0 * ta_zzzz_xzzzz_1[i] * fe_0 + ta_zzzz_xxzzzz_0[i] * pa_x[i] - ta_zzzz_xxzzzz_1[i] * pc_x[i];

        ta_xzzzz_xyyyyy_0[i] = ta_zzzz_yyyyy_0[i] * fe_0 - ta_zzzz_yyyyy_1[i] * fe_0 + ta_zzzz_xyyyyy_0[i] * pa_x[i] - ta_zzzz_xyyyyy_1[i] * pc_x[i];

        ta_xzzzz_xyyyyz_0[i] = ta_zzzz_yyyyz_0[i] * fe_0 - ta_zzzz_yyyyz_1[i] * fe_0 + ta_zzzz_xyyyyz_0[i] * pa_x[i] - ta_zzzz_xyyyyz_1[i] * pc_x[i];

        ta_xzzzz_xyyyzz_0[i] = ta_zzzz_yyyzz_0[i] * fe_0 - ta_zzzz_yyyzz_1[i] * fe_0 + ta_zzzz_xyyyzz_0[i] * pa_x[i] - ta_zzzz_xyyyzz_1[i] * pc_x[i];

        ta_xzzzz_xyyzzz_0[i] = ta_zzzz_yyzzz_0[i] * fe_0 - ta_zzzz_yyzzz_1[i] * fe_0 + ta_zzzz_xyyzzz_0[i] * pa_x[i] - ta_zzzz_xyyzzz_1[i] * pc_x[i];

        ta_xzzzz_xyzzzz_0[i] = ta_zzzz_yzzzz_0[i] * fe_0 - ta_zzzz_yzzzz_1[i] * fe_0 + ta_zzzz_xyzzzz_0[i] * pa_x[i] - ta_zzzz_xyzzzz_1[i] * pc_x[i];

        ta_xzzzz_xzzzzz_0[i] = ta_zzzz_zzzzz_0[i] * fe_0 - ta_zzzz_zzzzz_1[i] * fe_0 + ta_zzzz_xzzzzz_0[i] * pa_x[i] - ta_zzzz_xzzzzz_1[i] * pc_x[i];

        ta_xzzzz_yyyyyy_0[i] = ta_zzzz_yyyyyy_0[i] * pa_x[i] - ta_zzzz_yyyyyy_1[i] * pc_x[i];

        ta_xzzzz_yyyyyz_0[i] = ta_zzzz_yyyyyz_0[i] * pa_x[i] - ta_zzzz_yyyyyz_1[i] * pc_x[i];

        ta_xzzzz_yyyyzz_0[i] = ta_zzzz_yyyyzz_0[i] * pa_x[i] - ta_zzzz_yyyyzz_1[i] * pc_x[i];

        ta_xzzzz_yyyzzz_0[i] = ta_zzzz_yyyzzz_0[i] * pa_x[i] - ta_zzzz_yyyzzz_1[i] * pc_x[i];

        ta_xzzzz_yyzzzz_0[i] = ta_zzzz_yyzzzz_0[i] * pa_x[i] - ta_zzzz_yyzzzz_1[i] * pc_x[i];

        ta_xzzzz_yzzzzz_0[i] = ta_zzzz_yzzzzz_0[i] * pa_x[i] - ta_zzzz_yzzzzz_1[i] * pc_x[i];

        ta_xzzzz_zzzzzz_0[i] = ta_zzzz_zzzzzz_0[i] * pa_x[i] - ta_zzzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 420-448 components of targeted buffer : HI

    auto ta_yyyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 420);

    auto ta_yyyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 421);

    auto ta_yyyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 422);

    auto ta_yyyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 423);

    auto ta_yyyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 424);

    auto ta_yyyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 425);

    auto ta_yyyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 426);

    auto ta_yyyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 427);

    auto ta_yyyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 428);

    auto ta_yyyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 429);

    auto ta_yyyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 430);

    auto ta_yyyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 431);

    auto ta_yyyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 432);

    auto ta_yyyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 433);

    auto ta_yyyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 434);

    auto ta_yyyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 435);

    auto ta_yyyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 436);

    auto ta_yyyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 437);

    auto ta_yyyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 438);

    auto ta_yyyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 439);

    auto ta_yyyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 440);

    auto ta_yyyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 441);

    auto ta_yyyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 442);

    auto ta_yyyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 443);

    auto ta_yyyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 444);

    auto ta_yyyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 445);

    auto ta_yyyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 446);

    auto ta_yyyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 447);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta_yyy_xxxxxx_0,   \
                             ta_yyy_xxxxxx_1,   \
                             ta_yyy_xxxxxy_0,   \
                             ta_yyy_xxxxxy_1,   \
                             ta_yyy_xxxxxz_0,   \
                             ta_yyy_xxxxxz_1,   \
                             ta_yyy_xxxxyy_0,   \
                             ta_yyy_xxxxyy_1,   \
                             ta_yyy_xxxxyz_0,   \
                             ta_yyy_xxxxyz_1,   \
                             ta_yyy_xxxxzz_0,   \
                             ta_yyy_xxxxzz_1,   \
                             ta_yyy_xxxyyy_0,   \
                             ta_yyy_xxxyyy_1,   \
                             ta_yyy_xxxyyz_0,   \
                             ta_yyy_xxxyyz_1,   \
                             ta_yyy_xxxyzz_0,   \
                             ta_yyy_xxxyzz_1,   \
                             ta_yyy_xxxzzz_0,   \
                             ta_yyy_xxxzzz_1,   \
                             ta_yyy_xxyyyy_0,   \
                             ta_yyy_xxyyyy_1,   \
                             ta_yyy_xxyyyz_0,   \
                             ta_yyy_xxyyyz_1,   \
                             ta_yyy_xxyyzz_0,   \
                             ta_yyy_xxyyzz_1,   \
                             ta_yyy_xxyzzz_0,   \
                             ta_yyy_xxyzzz_1,   \
                             ta_yyy_xxzzzz_0,   \
                             ta_yyy_xxzzzz_1,   \
                             ta_yyy_xyyyyy_0,   \
                             ta_yyy_xyyyyy_1,   \
                             ta_yyy_xyyyyz_0,   \
                             ta_yyy_xyyyyz_1,   \
                             ta_yyy_xyyyzz_0,   \
                             ta_yyy_xyyyzz_1,   \
                             ta_yyy_xyyzzz_0,   \
                             ta_yyy_xyyzzz_1,   \
                             ta_yyy_xyzzzz_0,   \
                             ta_yyy_xyzzzz_1,   \
                             ta_yyy_xzzzzz_0,   \
                             ta_yyy_xzzzzz_1,   \
                             ta_yyy_yyyyyy_0,   \
                             ta_yyy_yyyyyy_1,   \
                             ta_yyy_yyyyyz_0,   \
                             ta_yyy_yyyyyz_1,   \
                             ta_yyy_yyyyzz_0,   \
                             ta_yyy_yyyyzz_1,   \
                             ta_yyy_yyyzzz_0,   \
                             ta_yyy_yyyzzz_1,   \
                             ta_yyy_yyzzzz_0,   \
                             ta_yyy_yyzzzz_1,   \
                             ta_yyy_yzzzzz_0,   \
                             ta_yyy_yzzzzz_1,   \
                             ta_yyy_zzzzzz_0,   \
                             ta_yyy_zzzzzz_1,   \
                             ta_yyyy_xxxxx_0,   \
                             ta_yyyy_xxxxx_1,   \
                             ta_yyyy_xxxxxx_0,  \
                             ta_yyyy_xxxxxx_1,  \
                             ta_yyyy_xxxxxy_0,  \
                             ta_yyyy_xxxxxy_1,  \
                             ta_yyyy_xxxxxz_0,  \
                             ta_yyyy_xxxxxz_1,  \
                             ta_yyyy_xxxxy_0,   \
                             ta_yyyy_xxxxy_1,   \
                             ta_yyyy_xxxxyy_0,  \
                             ta_yyyy_xxxxyy_1,  \
                             ta_yyyy_xxxxyz_0,  \
                             ta_yyyy_xxxxyz_1,  \
                             ta_yyyy_xxxxz_0,   \
                             ta_yyyy_xxxxz_1,   \
                             ta_yyyy_xxxxzz_0,  \
                             ta_yyyy_xxxxzz_1,  \
                             ta_yyyy_xxxyy_0,   \
                             ta_yyyy_xxxyy_1,   \
                             ta_yyyy_xxxyyy_0,  \
                             ta_yyyy_xxxyyy_1,  \
                             ta_yyyy_xxxyyz_0,  \
                             ta_yyyy_xxxyyz_1,  \
                             ta_yyyy_xxxyz_0,   \
                             ta_yyyy_xxxyz_1,   \
                             ta_yyyy_xxxyzz_0,  \
                             ta_yyyy_xxxyzz_1,  \
                             ta_yyyy_xxxzz_0,   \
                             ta_yyyy_xxxzz_1,   \
                             ta_yyyy_xxxzzz_0,  \
                             ta_yyyy_xxxzzz_1,  \
                             ta_yyyy_xxyyy_0,   \
                             ta_yyyy_xxyyy_1,   \
                             ta_yyyy_xxyyyy_0,  \
                             ta_yyyy_xxyyyy_1,  \
                             ta_yyyy_xxyyyz_0,  \
                             ta_yyyy_xxyyyz_1,  \
                             ta_yyyy_xxyyz_0,   \
                             ta_yyyy_xxyyz_1,   \
                             ta_yyyy_xxyyzz_0,  \
                             ta_yyyy_xxyyzz_1,  \
                             ta_yyyy_xxyzz_0,   \
                             ta_yyyy_xxyzz_1,   \
                             ta_yyyy_xxyzzz_0,  \
                             ta_yyyy_xxyzzz_1,  \
                             ta_yyyy_xxzzz_0,   \
                             ta_yyyy_xxzzz_1,   \
                             ta_yyyy_xxzzzz_0,  \
                             ta_yyyy_xxzzzz_1,  \
                             ta_yyyy_xyyyy_0,   \
                             ta_yyyy_xyyyy_1,   \
                             ta_yyyy_xyyyyy_0,  \
                             ta_yyyy_xyyyyy_1,  \
                             ta_yyyy_xyyyyz_0,  \
                             ta_yyyy_xyyyyz_1,  \
                             ta_yyyy_xyyyz_0,   \
                             ta_yyyy_xyyyz_1,   \
                             ta_yyyy_xyyyzz_0,  \
                             ta_yyyy_xyyyzz_1,  \
                             ta_yyyy_xyyzz_0,   \
                             ta_yyyy_xyyzz_1,   \
                             ta_yyyy_xyyzzz_0,  \
                             ta_yyyy_xyyzzz_1,  \
                             ta_yyyy_xyzzz_0,   \
                             ta_yyyy_xyzzz_1,   \
                             ta_yyyy_xyzzzz_0,  \
                             ta_yyyy_xyzzzz_1,  \
                             ta_yyyy_xzzzz_0,   \
                             ta_yyyy_xzzzz_1,   \
                             ta_yyyy_xzzzzz_0,  \
                             ta_yyyy_xzzzzz_1,  \
                             ta_yyyy_yyyyy_0,   \
                             ta_yyyy_yyyyy_1,   \
                             ta_yyyy_yyyyyy_0,  \
                             ta_yyyy_yyyyyy_1,  \
                             ta_yyyy_yyyyyz_0,  \
                             ta_yyyy_yyyyyz_1,  \
                             ta_yyyy_yyyyz_0,   \
                             ta_yyyy_yyyyz_1,   \
                             ta_yyyy_yyyyzz_0,  \
                             ta_yyyy_yyyyzz_1,  \
                             ta_yyyy_yyyzz_0,   \
                             ta_yyyy_yyyzz_1,   \
                             ta_yyyy_yyyzzz_0,  \
                             ta_yyyy_yyyzzz_1,  \
                             ta_yyyy_yyzzz_0,   \
                             ta_yyyy_yyzzz_1,   \
                             ta_yyyy_yyzzzz_0,  \
                             ta_yyyy_yyzzzz_1,  \
                             ta_yyyy_yzzzz_0,   \
                             ta_yyyy_yzzzz_1,   \
                             ta_yyyy_yzzzzz_0,  \
                             ta_yyyy_yzzzzz_1,  \
                             ta_yyyy_zzzzz_0,   \
                             ta_yyyy_zzzzz_1,   \
                             ta_yyyy_zzzzzz_0,  \
                             ta_yyyy_zzzzzz_1,  \
                             ta_yyyyy_xxxxxx_0, \
                             ta_yyyyy_xxxxxy_0, \
                             ta_yyyyy_xxxxxz_0, \
                             ta_yyyyy_xxxxyy_0, \
                             ta_yyyyy_xxxxyz_0, \
                             ta_yyyyy_xxxxzz_0, \
                             ta_yyyyy_xxxyyy_0, \
                             ta_yyyyy_xxxyyz_0, \
                             ta_yyyyy_xxxyzz_0, \
                             ta_yyyyy_xxxzzz_0, \
                             ta_yyyyy_xxyyyy_0, \
                             ta_yyyyy_xxyyyz_0, \
                             ta_yyyyy_xxyyzz_0, \
                             ta_yyyyy_xxyzzz_0, \
                             ta_yyyyy_xxzzzz_0, \
                             ta_yyyyy_xyyyyy_0, \
                             ta_yyyyy_xyyyyz_0, \
                             ta_yyyyy_xyyyzz_0, \
                             ta_yyyyy_xyyzzz_0, \
                             ta_yyyyy_xyzzzz_0, \
                             ta_yyyyy_xzzzzz_0, \
                             ta_yyyyy_yyyyyy_0, \
                             ta_yyyyy_yyyyyz_0, \
                             ta_yyyyy_yyyyzz_0, \
                             ta_yyyyy_yyyzzz_0, \
                             ta_yyyyy_yyzzzz_0, \
                             ta_yyyyy_yzzzzz_0, \
                             ta_yyyyy_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyy_xxxxxx_0[i] =
            4.0 * ta_yyy_xxxxxx_0[i] * fe_0 - 4.0 * ta_yyy_xxxxxx_1[i] * fe_0 + ta_yyyy_xxxxxx_0[i] * pa_y[i] - ta_yyyy_xxxxxx_1[i] * pc_y[i];

        ta_yyyyy_xxxxxy_0[i] = 4.0 * ta_yyy_xxxxxy_0[i] * fe_0 - 4.0 * ta_yyy_xxxxxy_1[i] * fe_0 + ta_yyyy_xxxxx_0[i] * fe_0 -
                               ta_yyyy_xxxxx_1[i] * fe_0 + ta_yyyy_xxxxxy_0[i] * pa_y[i] - ta_yyyy_xxxxxy_1[i] * pc_y[i];

        ta_yyyyy_xxxxxz_0[i] =
            4.0 * ta_yyy_xxxxxz_0[i] * fe_0 - 4.0 * ta_yyy_xxxxxz_1[i] * fe_0 + ta_yyyy_xxxxxz_0[i] * pa_y[i] - ta_yyyy_xxxxxz_1[i] * pc_y[i];

        ta_yyyyy_xxxxyy_0[i] = 4.0 * ta_yyy_xxxxyy_0[i] * fe_0 - 4.0 * ta_yyy_xxxxyy_1[i] * fe_0 + 2.0 * ta_yyyy_xxxxy_0[i] * fe_0 -
                               2.0 * ta_yyyy_xxxxy_1[i] * fe_0 + ta_yyyy_xxxxyy_0[i] * pa_y[i] - ta_yyyy_xxxxyy_1[i] * pc_y[i];

        ta_yyyyy_xxxxyz_0[i] = 4.0 * ta_yyy_xxxxyz_0[i] * fe_0 - 4.0 * ta_yyy_xxxxyz_1[i] * fe_0 + ta_yyyy_xxxxz_0[i] * fe_0 -
                               ta_yyyy_xxxxz_1[i] * fe_0 + ta_yyyy_xxxxyz_0[i] * pa_y[i] - ta_yyyy_xxxxyz_1[i] * pc_y[i];

        ta_yyyyy_xxxxzz_0[i] =
            4.0 * ta_yyy_xxxxzz_0[i] * fe_0 - 4.0 * ta_yyy_xxxxzz_1[i] * fe_0 + ta_yyyy_xxxxzz_0[i] * pa_y[i] - ta_yyyy_xxxxzz_1[i] * pc_y[i];

        ta_yyyyy_xxxyyy_0[i] = 4.0 * ta_yyy_xxxyyy_0[i] * fe_0 - 4.0 * ta_yyy_xxxyyy_1[i] * fe_0 + 3.0 * ta_yyyy_xxxyy_0[i] * fe_0 -
                               3.0 * ta_yyyy_xxxyy_1[i] * fe_0 + ta_yyyy_xxxyyy_0[i] * pa_y[i] - ta_yyyy_xxxyyy_1[i] * pc_y[i];

        ta_yyyyy_xxxyyz_0[i] = 4.0 * ta_yyy_xxxyyz_0[i] * fe_0 - 4.0 * ta_yyy_xxxyyz_1[i] * fe_0 + 2.0 * ta_yyyy_xxxyz_0[i] * fe_0 -
                               2.0 * ta_yyyy_xxxyz_1[i] * fe_0 + ta_yyyy_xxxyyz_0[i] * pa_y[i] - ta_yyyy_xxxyyz_1[i] * pc_y[i];

        ta_yyyyy_xxxyzz_0[i] = 4.0 * ta_yyy_xxxyzz_0[i] * fe_0 - 4.0 * ta_yyy_xxxyzz_1[i] * fe_0 + ta_yyyy_xxxzz_0[i] * fe_0 -
                               ta_yyyy_xxxzz_1[i] * fe_0 + ta_yyyy_xxxyzz_0[i] * pa_y[i] - ta_yyyy_xxxyzz_1[i] * pc_y[i];

        ta_yyyyy_xxxzzz_0[i] =
            4.0 * ta_yyy_xxxzzz_0[i] * fe_0 - 4.0 * ta_yyy_xxxzzz_1[i] * fe_0 + ta_yyyy_xxxzzz_0[i] * pa_y[i] - ta_yyyy_xxxzzz_1[i] * pc_y[i];

        ta_yyyyy_xxyyyy_0[i] = 4.0 * ta_yyy_xxyyyy_0[i] * fe_0 - 4.0 * ta_yyy_xxyyyy_1[i] * fe_0 + 4.0 * ta_yyyy_xxyyy_0[i] * fe_0 -
                               4.0 * ta_yyyy_xxyyy_1[i] * fe_0 + ta_yyyy_xxyyyy_0[i] * pa_y[i] - ta_yyyy_xxyyyy_1[i] * pc_y[i];

        ta_yyyyy_xxyyyz_0[i] = 4.0 * ta_yyy_xxyyyz_0[i] * fe_0 - 4.0 * ta_yyy_xxyyyz_1[i] * fe_0 + 3.0 * ta_yyyy_xxyyz_0[i] * fe_0 -
                               3.0 * ta_yyyy_xxyyz_1[i] * fe_0 + ta_yyyy_xxyyyz_0[i] * pa_y[i] - ta_yyyy_xxyyyz_1[i] * pc_y[i];

        ta_yyyyy_xxyyzz_0[i] = 4.0 * ta_yyy_xxyyzz_0[i] * fe_0 - 4.0 * ta_yyy_xxyyzz_1[i] * fe_0 + 2.0 * ta_yyyy_xxyzz_0[i] * fe_0 -
                               2.0 * ta_yyyy_xxyzz_1[i] * fe_0 + ta_yyyy_xxyyzz_0[i] * pa_y[i] - ta_yyyy_xxyyzz_1[i] * pc_y[i];

        ta_yyyyy_xxyzzz_0[i] = 4.0 * ta_yyy_xxyzzz_0[i] * fe_0 - 4.0 * ta_yyy_xxyzzz_1[i] * fe_0 + ta_yyyy_xxzzz_0[i] * fe_0 -
                               ta_yyyy_xxzzz_1[i] * fe_0 + ta_yyyy_xxyzzz_0[i] * pa_y[i] - ta_yyyy_xxyzzz_1[i] * pc_y[i];

        ta_yyyyy_xxzzzz_0[i] =
            4.0 * ta_yyy_xxzzzz_0[i] * fe_0 - 4.0 * ta_yyy_xxzzzz_1[i] * fe_0 + ta_yyyy_xxzzzz_0[i] * pa_y[i] - ta_yyyy_xxzzzz_1[i] * pc_y[i];

        ta_yyyyy_xyyyyy_0[i] = 4.0 * ta_yyy_xyyyyy_0[i] * fe_0 - 4.0 * ta_yyy_xyyyyy_1[i] * fe_0 + 5.0 * ta_yyyy_xyyyy_0[i] * fe_0 -
                               5.0 * ta_yyyy_xyyyy_1[i] * fe_0 + ta_yyyy_xyyyyy_0[i] * pa_y[i] - ta_yyyy_xyyyyy_1[i] * pc_y[i];

        ta_yyyyy_xyyyyz_0[i] = 4.0 * ta_yyy_xyyyyz_0[i] * fe_0 - 4.0 * ta_yyy_xyyyyz_1[i] * fe_0 + 4.0 * ta_yyyy_xyyyz_0[i] * fe_0 -
                               4.0 * ta_yyyy_xyyyz_1[i] * fe_0 + ta_yyyy_xyyyyz_0[i] * pa_y[i] - ta_yyyy_xyyyyz_1[i] * pc_y[i];

        ta_yyyyy_xyyyzz_0[i] = 4.0 * ta_yyy_xyyyzz_0[i] * fe_0 - 4.0 * ta_yyy_xyyyzz_1[i] * fe_0 + 3.0 * ta_yyyy_xyyzz_0[i] * fe_0 -
                               3.0 * ta_yyyy_xyyzz_1[i] * fe_0 + ta_yyyy_xyyyzz_0[i] * pa_y[i] - ta_yyyy_xyyyzz_1[i] * pc_y[i];

        ta_yyyyy_xyyzzz_0[i] = 4.0 * ta_yyy_xyyzzz_0[i] * fe_0 - 4.0 * ta_yyy_xyyzzz_1[i] * fe_0 + 2.0 * ta_yyyy_xyzzz_0[i] * fe_0 -
                               2.0 * ta_yyyy_xyzzz_1[i] * fe_0 + ta_yyyy_xyyzzz_0[i] * pa_y[i] - ta_yyyy_xyyzzz_1[i] * pc_y[i];

        ta_yyyyy_xyzzzz_0[i] = 4.0 * ta_yyy_xyzzzz_0[i] * fe_0 - 4.0 * ta_yyy_xyzzzz_1[i] * fe_0 + ta_yyyy_xzzzz_0[i] * fe_0 -
                               ta_yyyy_xzzzz_1[i] * fe_0 + ta_yyyy_xyzzzz_0[i] * pa_y[i] - ta_yyyy_xyzzzz_1[i] * pc_y[i];

        ta_yyyyy_xzzzzz_0[i] =
            4.0 * ta_yyy_xzzzzz_0[i] * fe_0 - 4.0 * ta_yyy_xzzzzz_1[i] * fe_0 + ta_yyyy_xzzzzz_0[i] * pa_y[i] - ta_yyyy_xzzzzz_1[i] * pc_y[i];

        ta_yyyyy_yyyyyy_0[i] = 4.0 * ta_yyy_yyyyyy_0[i] * fe_0 - 4.0 * ta_yyy_yyyyyy_1[i] * fe_0 + 6.0 * ta_yyyy_yyyyy_0[i] * fe_0 -
                               6.0 * ta_yyyy_yyyyy_1[i] * fe_0 + ta_yyyy_yyyyyy_0[i] * pa_y[i] - ta_yyyy_yyyyyy_1[i] * pc_y[i];

        ta_yyyyy_yyyyyz_0[i] = 4.0 * ta_yyy_yyyyyz_0[i] * fe_0 - 4.0 * ta_yyy_yyyyyz_1[i] * fe_0 + 5.0 * ta_yyyy_yyyyz_0[i] * fe_0 -
                               5.0 * ta_yyyy_yyyyz_1[i] * fe_0 + ta_yyyy_yyyyyz_0[i] * pa_y[i] - ta_yyyy_yyyyyz_1[i] * pc_y[i];

        ta_yyyyy_yyyyzz_0[i] = 4.0 * ta_yyy_yyyyzz_0[i] * fe_0 - 4.0 * ta_yyy_yyyyzz_1[i] * fe_0 + 4.0 * ta_yyyy_yyyzz_0[i] * fe_0 -
                               4.0 * ta_yyyy_yyyzz_1[i] * fe_0 + ta_yyyy_yyyyzz_0[i] * pa_y[i] - ta_yyyy_yyyyzz_1[i] * pc_y[i];

        ta_yyyyy_yyyzzz_0[i] = 4.0 * ta_yyy_yyyzzz_0[i] * fe_0 - 4.0 * ta_yyy_yyyzzz_1[i] * fe_0 + 3.0 * ta_yyyy_yyzzz_0[i] * fe_0 -
                               3.0 * ta_yyyy_yyzzz_1[i] * fe_0 + ta_yyyy_yyyzzz_0[i] * pa_y[i] - ta_yyyy_yyyzzz_1[i] * pc_y[i];

        ta_yyyyy_yyzzzz_0[i] = 4.0 * ta_yyy_yyzzzz_0[i] * fe_0 - 4.0 * ta_yyy_yyzzzz_1[i] * fe_0 + 2.0 * ta_yyyy_yzzzz_0[i] * fe_0 -
                               2.0 * ta_yyyy_yzzzz_1[i] * fe_0 + ta_yyyy_yyzzzz_0[i] * pa_y[i] - ta_yyyy_yyzzzz_1[i] * pc_y[i];

        ta_yyyyy_yzzzzz_0[i] = 4.0 * ta_yyy_yzzzzz_0[i] * fe_0 - 4.0 * ta_yyy_yzzzzz_1[i] * fe_0 + ta_yyyy_zzzzz_0[i] * fe_0 -
                               ta_yyyy_zzzzz_1[i] * fe_0 + ta_yyyy_yzzzzz_0[i] * pa_y[i] - ta_yyyy_yzzzzz_1[i] * pc_y[i];

        ta_yyyyy_zzzzzz_0[i] =
            4.0 * ta_yyy_zzzzzz_0[i] * fe_0 - 4.0 * ta_yyy_zzzzzz_1[i] * fe_0 + ta_yyyy_zzzzzz_0[i] * pa_y[i] - ta_yyyy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 448-476 components of targeted buffer : HI

    auto ta_yyyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 448);

    auto ta_yyyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 449);

    auto ta_yyyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 450);

    auto ta_yyyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 451);

    auto ta_yyyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 452);

    auto ta_yyyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 453);

    auto ta_yyyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 454);

    auto ta_yyyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 455);

    auto ta_yyyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 456);

    auto ta_yyyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 457);

    auto ta_yyyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 458);

    auto ta_yyyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 459);

    auto ta_yyyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 460);

    auto ta_yyyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 461);

    auto ta_yyyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 462);

    auto ta_yyyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 463);

    auto ta_yyyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 464);

    auto ta_yyyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 465);

    auto ta_yyyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 466);

    auto ta_yyyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 467);

    auto ta_yyyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 468);

    auto ta_yyyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 469);

    auto ta_yyyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 470);

    auto ta_yyyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 471);

    auto ta_yyyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 472);

    auto ta_yyyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 473);

    auto ta_yyyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 474);

    auto ta_yyyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 475);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta_yyyy_xxxxxx_0,  \
                             ta_yyyy_xxxxxx_1,  \
                             ta_yyyy_xxxxxy_0,  \
                             ta_yyyy_xxxxxy_1,  \
                             ta_yyyy_xxxxy_0,   \
                             ta_yyyy_xxxxy_1,   \
                             ta_yyyy_xxxxyy_0,  \
                             ta_yyyy_xxxxyy_1,  \
                             ta_yyyy_xxxxyz_0,  \
                             ta_yyyy_xxxxyz_1,  \
                             ta_yyyy_xxxyy_0,   \
                             ta_yyyy_xxxyy_1,   \
                             ta_yyyy_xxxyyy_0,  \
                             ta_yyyy_xxxyyy_1,  \
                             ta_yyyy_xxxyyz_0,  \
                             ta_yyyy_xxxyyz_1,  \
                             ta_yyyy_xxxyz_0,   \
                             ta_yyyy_xxxyz_1,   \
                             ta_yyyy_xxxyzz_0,  \
                             ta_yyyy_xxxyzz_1,  \
                             ta_yyyy_xxyyy_0,   \
                             ta_yyyy_xxyyy_1,   \
                             ta_yyyy_xxyyyy_0,  \
                             ta_yyyy_xxyyyy_1,  \
                             ta_yyyy_xxyyyz_0,  \
                             ta_yyyy_xxyyyz_1,  \
                             ta_yyyy_xxyyz_0,   \
                             ta_yyyy_xxyyz_1,   \
                             ta_yyyy_xxyyzz_0,  \
                             ta_yyyy_xxyyzz_1,  \
                             ta_yyyy_xxyzz_0,   \
                             ta_yyyy_xxyzz_1,   \
                             ta_yyyy_xxyzzz_0,  \
                             ta_yyyy_xxyzzz_1,  \
                             ta_yyyy_xyyyy_0,   \
                             ta_yyyy_xyyyy_1,   \
                             ta_yyyy_xyyyyy_0,  \
                             ta_yyyy_xyyyyy_1,  \
                             ta_yyyy_xyyyyz_0,  \
                             ta_yyyy_xyyyyz_1,  \
                             ta_yyyy_xyyyz_0,   \
                             ta_yyyy_xyyyz_1,   \
                             ta_yyyy_xyyyzz_0,  \
                             ta_yyyy_xyyyzz_1,  \
                             ta_yyyy_xyyzz_0,   \
                             ta_yyyy_xyyzz_1,   \
                             ta_yyyy_xyyzzz_0,  \
                             ta_yyyy_xyyzzz_1,  \
                             ta_yyyy_xyzzz_0,   \
                             ta_yyyy_xyzzz_1,   \
                             ta_yyyy_xyzzzz_0,  \
                             ta_yyyy_xyzzzz_1,  \
                             ta_yyyy_yyyyy_0,   \
                             ta_yyyy_yyyyy_1,   \
                             ta_yyyy_yyyyyy_0,  \
                             ta_yyyy_yyyyyy_1,  \
                             ta_yyyy_yyyyyz_0,  \
                             ta_yyyy_yyyyyz_1,  \
                             ta_yyyy_yyyyz_0,   \
                             ta_yyyy_yyyyz_1,   \
                             ta_yyyy_yyyyzz_0,  \
                             ta_yyyy_yyyyzz_1,  \
                             ta_yyyy_yyyzz_0,   \
                             ta_yyyy_yyyzz_1,   \
                             ta_yyyy_yyyzzz_0,  \
                             ta_yyyy_yyyzzz_1,  \
                             ta_yyyy_yyzzz_0,   \
                             ta_yyyy_yyzzz_1,   \
                             ta_yyyy_yyzzzz_0,  \
                             ta_yyyy_yyzzzz_1,  \
                             ta_yyyy_yzzzz_0,   \
                             ta_yyyy_yzzzz_1,   \
                             ta_yyyy_yzzzzz_0,  \
                             ta_yyyy_yzzzzz_1,  \
                             ta_yyyyz_xxxxxx_0, \
                             ta_yyyyz_xxxxxy_0, \
                             ta_yyyyz_xxxxxz_0, \
                             ta_yyyyz_xxxxyy_0, \
                             ta_yyyyz_xxxxyz_0, \
                             ta_yyyyz_xxxxzz_0, \
                             ta_yyyyz_xxxyyy_0, \
                             ta_yyyyz_xxxyyz_0, \
                             ta_yyyyz_xxxyzz_0, \
                             ta_yyyyz_xxxzzz_0, \
                             ta_yyyyz_xxyyyy_0, \
                             ta_yyyyz_xxyyyz_0, \
                             ta_yyyyz_xxyyzz_0, \
                             ta_yyyyz_xxyzzz_0, \
                             ta_yyyyz_xxzzzz_0, \
                             ta_yyyyz_xyyyyy_0, \
                             ta_yyyyz_xyyyyz_0, \
                             ta_yyyyz_xyyyzz_0, \
                             ta_yyyyz_xyyzzz_0, \
                             ta_yyyyz_xyzzzz_0, \
                             ta_yyyyz_xzzzzz_0, \
                             ta_yyyyz_yyyyyy_0, \
                             ta_yyyyz_yyyyyz_0, \
                             ta_yyyyz_yyyyzz_0, \
                             ta_yyyyz_yyyzzz_0, \
                             ta_yyyyz_yyzzzz_0, \
                             ta_yyyyz_yzzzzz_0, \
                             ta_yyyyz_zzzzzz_0, \
                             ta_yyyz_xxxxxz_0,  \
                             ta_yyyz_xxxxxz_1,  \
                             ta_yyyz_xxxxzz_0,  \
                             ta_yyyz_xxxxzz_1,  \
                             ta_yyyz_xxxzzz_0,  \
                             ta_yyyz_xxxzzz_1,  \
                             ta_yyyz_xxzzzz_0,  \
                             ta_yyyz_xxzzzz_1,  \
                             ta_yyyz_xzzzzz_0,  \
                             ta_yyyz_xzzzzz_1,  \
                             ta_yyyz_zzzzzz_0,  \
                             ta_yyyz_zzzzzz_1,  \
                             ta_yyz_xxxxxz_0,   \
                             ta_yyz_xxxxxz_1,   \
                             ta_yyz_xxxxzz_0,   \
                             ta_yyz_xxxxzz_1,   \
                             ta_yyz_xxxzzz_0,   \
                             ta_yyz_xxxzzz_1,   \
                             ta_yyz_xxzzzz_0,   \
                             ta_yyz_xxzzzz_1,   \
                             ta_yyz_xzzzzz_0,   \
                             ta_yyz_xzzzzz_1,   \
                             ta_yyz_zzzzzz_0,   \
                             ta_yyz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyz_xxxxxx_0[i] = ta_yyyy_xxxxxx_0[i] * pa_z[i] - ta_yyyy_xxxxxx_1[i] * pc_z[i];

        ta_yyyyz_xxxxxy_0[i] = ta_yyyy_xxxxxy_0[i] * pa_z[i] - ta_yyyy_xxxxxy_1[i] * pc_z[i];

        ta_yyyyz_xxxxxz_0[i] =
            3.0 * ta_yyz_xxxxxz_0[i] * fe_0 - 3.0 * ta_yyz_xxxxxz_1[i] * fe_0 + ta_yyyz_xxxxxz_0[i] * pa_y[i] - ta_yyyz_xxxxxz_1[i] * pc_y[i];

        ta_yyyyz_xxxxyy_0[i] = ta_yyyy_xxxxyy_0[i] * pa_z[i] - ta_yyyy_xxxxyy_1[i] * pc_z[i];

        ta_yyyyz_xxxxyz_0[i] = ta_yyyy_xxxxy_0[i] * fe_0 - ta_yyyy_xxxxy_1[i] * fe_0 + ta_yyyy_xxxxyz_0[i] * pa_z[i] - ta_yyyy_xxxxyz_1[i] * pc_z[i];

        ta_yyyyz_xxxxzz_0[i] =
            3.0 * ta_yyz_xxxxzz_0[i] * fe_0 - 3.0 * ta_yyz_xxxxzz_1[i] * fe_0 + ta_yyyz_xxxxzz_0[i] * pa_y[i] - ta_yyyz_xxxxzz_1[i] * pc_y[i];

        ta_yyyyz_xxxyyy_0[i] = ta_yyyy_xxxyyy_0[i] * pa_z[i] - ta_yyyy_xxxyyy_1[i] * pc_z[i];

        ta_yyyyz_xxxyyz_0[i] = ta_yyyy_xxxyy_0[i] * fe_0 - ta_yyyy_xxxyy_1[i] * fe_0 + ta_yyyy_xxxyyz_0[i] * pa_z[i] - ta_yyyy_xxxyyz_1[i] * pc_z[i];

        ta_yyyyz_xxxyzz_0[i] =
            2.0 * ta_yyyy_xxxyz_0[i] * fe_0 - 2.0 * ta_yyyy_xxxyz_1[i] * fe_0 + ta_yyyy_xxxyzz_0[i] * pa_z[i] - ta_yyyy_xxxyzz_1[i] * pc_z[i];

        ta_yyyyz_xxxzzz_0[i] =
            3.0 * ta_yyz_xxxzzz_0[i] * fe_0 - 3.0 * ta_yyz_xxxzzz_1[i] * fe_0 + ta_yyyz_xxxzzz_0[i] * pa_y[i] - ta_yyyz_xxxzzz_1[i] * pc_y[i];

        ta_yyyyz_xxyyyy_0[i] = ta_yyyy_xxyyyy_0[i] * pa_z[i] - ta_yyyy_xxyyyy_1[i] * pc_z[i];

        ta_yyyyz_xxyyyz_0[i] = ta_yyyy_xxyyy_0[i] * fe_0 - ta_yyyy_xxyyy_1[i] * fe_0 + ta_yyyy_xxyyyz_0[i] * pa_z[i] - ta_yyyy_xxyyyz_1[i] * pc_z[i];

        ta_yyyyz_xxyyzz_0[i] =
            2.0 * ta_yyyy_xxyyz_0[i] * fe_0 - 2.0 * ta_yyyy_xxyyz_1[i] * fe_0 + ta_yyyy_xxyyzz_0[i] * pa_z[i] - ta_yyyy_xxyyzz_1[i] * pc_z[i];

        ta_yyyyz_xxyzzz_0[i] =
            3.0 * ta_yyyy_xxyzz_0[i] * fe_0 - 3.0 * ta_yyyy_xxyzz_1[i] * fe_0 + ta_yyyy_xxyzzz_0[i] * pa_z[i] - ta_yyyy_xxyzzz_1[i] * pc_z[i];

        ta_yyyyz_xxzzzz_0[i] =
            3.0 * ta_yyz_xxzzzz_0[i] * fe_0 - 3.0 * ta_yyz_xxzzzz_1[i] * fe_0 + ta_yyyz_xxzzzz_0[i] * pa_y[i] - ta_yyyz_xxzzzz_1[i] * pc_y[i];

        ta_yyyyz_xyyyyy_0[i] = ta_yyyy_xyyyyy_0[i] * pa_z[i] - ta_yyyy_xyyyyy_1[i] * pc_z[i];

        ta_yyyyz_xyyyyz_0[i] = ta_yyyy_xyyyy_0[i] * fe_0 - ta_yyyy_xyyyy_1[i] * fe_0 + ta_yyyy_xyyyyz_0[i] * pa_z[i] - ta_yyyy_xyyyyz_1[i] * pc_z[i];

        ta_yyyyz_xyyyzz_0[i] =
            2.0 * ta_yyyy_xyyyz_0[i] * fe_0 - 2.0 * ta_yyyy_xyyyz_1[i] * fe_0 + ta_yyyy_xyyyzz_0[i] * pa_z[i] - ta_yyyy_xyyyzz_1[i] * pc_z[i];

        ta_yyyyz_xyyzzz_0[i] =
            3.0 * ta_yyyy_xyyzz_0[i] * fe_0 - 3.0 * ta_yyyy_xyyzz_1[i] * fe_0 + ta_yyyy_xyyzzz_0[i] * pa_z[i] - ta_yyyy_xyyzzz_1[i] * pc_z[i];

        ta_yyyyz_xyzzzz_0[i] =
            4.0 * ta_yyyy_xyzzz_0[i] * fe_0 - 4.0 * ta_yyyy_xyzzz_1[i] * fe_0 + ta_yyyy_xyzzzz_0[i] * pa_z[i] - ta_yyyy_xyzzzz_1[i] * pc_z[i];

        ta_yyyyz_xzzzzz_0[i] =
            3.0 * ta_yyz_xzzzzz_0[i] * fe_0 - 3.0 * ta_yyz_xzzzzz_1[i] * fe_0 + ta_yyyz_xzzzzz_0[i] * pa_y[i] - ta_yyyz_xzzzzz_1[i] * pc_y[i];

        ta_yyyyz_yyyyyy_0[i] = ta_yyyy_yyyyyy_0[i] * pa_z[i] - ta_yyyy_yyyyyy_1[i] * pc_z[i];

        ta_yyyyz_yyyyyz_0[i] = ta_yyyy_yyyyy_0[i] * fe_0 - ta_yyyy_yyyyy_1[i] * fe_0 + ta_yyyy_yyyyyz_0[i] * pa_z[i] - ta_yyyy_yyyyyz_1[i] * pc_z[i];

        ta_yyyyz_yyyyzz_0[i] =
            2.0 * ta_yyyy_yyyyz_0[i] * fe_0 - 2.0 * ta_yyyy_yyyyz_1[i] * fe_0 + ta_yyyy_yyyyzz_0[i] * pa_z[i] - ta_yyyy_yyyyzz_1[i] * pc_z[i];

        ta_yyyyz_yyyzzz_0[i] =
            3.0 * ta_yyyy_yyyzz_0[i] * fe_0 - 3.0 * ta_yyyy_yyyzz_1[i] * fe_0 + ta_yyyy_yyyzzz_0[i] * pa_z[i] - ta_yyyy_yyyzzz_1[i] * pc_z[i];

        ta_yyyyz_yyzzzz_0[i] =
            4.0 * ta_yyyy_yyzzz_0[i] * fe_0 - 4.0 * ta_yyyy_yyzzz_1[i] * fe_0 + ta_yyyy_yyzzzz_0[i] * pa_z[i] - ta_yyyy_yyzzzz_1[i] * pc_z[i];

        ta_yyyyz_yzzzzz_0[i] =
            5.0 * ta_yyyy_yzzzz_0[i] * fe_0 - 5.0 * ta_yyyy_yzzzz_1[i] * fe_0 + ta_yyyy_yzzzzz_0[i] * pa_z[i] - ta_yyyy_yzzzzz_1[i] * pc_z[i];

        ta_yyyyz_zzzzzz_0[i] =
            3.0 * ta_yyz_zzzzzz_0[i] * fe_0 - 3.0 * ta_yyz_zzzzzz_1[i] * fe_0 + ta_yyyz_zzzzzz_0[i] * pa_y[i] - ta_yyyz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 476-504 components of targeted buffer : HI

    auto ta_yyyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 476);

    auto ta_yyyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 477);

    auto ta_yyyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 478);

    auto ta_yyyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 479);

    auto ta_yyyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 480);

    auto ta_yyyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 481);

    auto ta_yyyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 482);

    auto ta_yyyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 483);

    auto ta_yyyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 484);

    auto ta_yyyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 485);

    auto ta_yyyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 486);

    auto ta_yyyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 487);

    auto ta_yyyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 488);

    auto ta_yyyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 489);

    auto ta_yyyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 490);

    auto ta_yyyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 491);

    auto ta_yyyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 492);

    auto ta_yyyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 493);

    auto ta_yyyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 494);

    auto ta_yyyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 495);

    auto ta_yyyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 496);

    auto ta_yyyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 497);

    auto ta_yyyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 498);

    auto ta_yyyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 499);

    auto ta_yyyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 500);

    auto ta_yyyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 501);

    auto ta_yyyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 502);

    auto ta_yyyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 503);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta_yyy_xxxxxy_0,   \
                             ta_yyy_xxxxxy_1,   \
                             ta_yyy_xxxxyy_0,   \
                             ta_yyy_xxxxyy_1,   \
                             ta_yyy_xxxyyy_0,   \
                             ta_yyy_xxxyyy_1,   \
                             ta_yyy_xxyyyy_0,   \
                             ta_yyy_xxyyyy_1,   \
                             ta_yyy_xyyyyy_0,   \
                             ta_yyy_xyyyyy_1,   \
                             ta_yyy_yyyyyy_0,   \
                             ta_yyy_yyyyyy_1,   \
                             ta_yyyz_xxxxxy_0,  \
                             ta_yyyz_xxxxxy_1,  \
                             ta_yyyz_xxxxyy_0,  \
                             ta_yyyz_xxxxyy_1,  \
                             ta_yyyz_xxxyyy_0,  \
                             ta_yyyz_xxxyyy_1,  \
                             ta_yyyz_xxyyyy_0,  \
                             ta_yyyz_xxyyyy_1,  \
                             ta_yyyz_xyyyyy_0,  \
                             ta_yyyz_xyyyyy_1,  \
                             ta_yyyz_yyyyyy_0,  \
                             ta_yyyz_yyyyyy_1,  \
                             ta_yyyzz_xxxxxx_0, \
                             ta_yyyzz_xxxxxy_0, \
                             ta_yyyzz_xxxxxz_0, \
                             ta_yyyzz_xxxxyy_0, \
                             ta_yyyzz_xxxxyz_0, \
                             ta_yyyzz_xxxxzz_0, \
                             ta_yyyzz_xxxyyy_0, \
                             ta_yyyzz_xxxyyz_0, \
                             ta_yyyzz_xxxyzz_0, \
                             ta_yyyzz_xxxzzz_0, \
                             ta_yyyzz_xxyyyy_0, \
                             ta_yyyzz_xxyyyz_0, \
                             ta_yyyzz_xxyyzz_0, \
                             ta_yyyzz_xxyzzz_0, \
                             ta_yyyzz_xxzzzz_0, \
                             ta_yyyzz_xyyyyy_0, \
                             ta_yyyzz_xyyyyz_0, \
                             ta_yyyzz_xyyyzz_0, \
                             ta_yyyzz_xyyzzz_0, \
                             ta_yyyzz_xyzzzz_0, \
                             ta_yyyzz_xzzzzz_0, \
                             ta_yyyzz_yyyyyy_0, \
                             ta_yyyzz_yyyyyz_0, \
                             ta_yyyzz_yyyyzz_0, \
                             ta_yyyzz_yyyzzz_0, \
                             ta_yyyzz_yyzzzz_0, \
                             ta_yyyzz_yzzzzz_0, \
                             ta_yyyzz_zzzzzz_0, \
                             ta_yyzz_xxxxxx_0,  \
                             ta_yyzz_xxxxxx_1,  \
                             ta_yyzz_xxxxxz_0,  \
                             ta_yyzz_xxxxxz_1,  \
                             ta_yyzz_xxxxyz_0,  \
                             ta_yyzz_xxxxyz_1,  \
                             ta_yyzz_xxxxz_0,   \
                             ta_yyzz_xxxxz_1,   \
                             ta_yyzz_xxxxzz_0,  \
                             ta_yyzz_xxxxzz_1,  \
                             ta_yyzz_xxxyyz_0,  \
                             ta_yyzz_xxxyyz_1,  \
                             ta_yyzz_xxxyz_0,   \
                             ta_yyzz_xxxyz_1,   \
                             ta_yyzz_xxxyzz_0,  \
                             ta_yyzz_xxxyzz_1,  \
                             ta_yyzz_xxxzz_0,   \
                             ta_yyzz_xxxzz_1,   \
                             ta_yyzz_xxxzzz_0,  \
                             ta_yyzz_xxxzzz_1,  \
                             ta_yyzz_xxyyyz_0,  \
                             ta_yyzz_xxyyyz_1,  \
                             ta_yyzz_xxyyz_0,   \
                             ta_yyzz_xxyyz_1,   \
                             ta_yyzz_xxyyzz_0,  \
                             ta_yyzz_xxyyzz_1,  \
                             ta_yyzz_xxyzz_0,   \
                             ta_yyzz_xxyzz_1,   \
                             ta_yyzz_xxyzzz_0,  \
                             ta_yyzz_xxyzzz_1,  \
                             ta_yyzz_xxzzz_0,   \
                             ta_yyzz_xxzzz_1,   \
                             ta_yyzz_xxzzzz_0,  \
                             ta_yyzz_xxzzzz_1,  \
                             ta_yyzz_xyyyyz_0,  \
                             ta_yyzz_xyyyyz_1,  \
                             ta_yyzz_xyyyz_0,   \
                             ta_yyzz_xyyyz_1,   \
                             ta_yyzz_xyyyzz_0,  \
                             ta_yyzz_xyyyzz_1,  \
                             ta_yyzz_xyyzz_0,   \
                             ta_yyzz_xyyzz_1,   \
                             ta_yyzz_xyyzzz_0,  \
                             ta_yyzz_xyyzzz_1,  \
                             ta_yyzz_xyzzz_0,   \
                             ta_yyzz_xyzzz_1,   \
                             ta_yyzz_xyzzzz_0,  \
                             ta_yyzz_xyzzzz_1,  \
                             ta_yyzz_xzzzz_0,   \
                             ta_yyzz_xzzzz_1,   \
                             ta_yyzz_xzzzzz_0,  \
                             ta_yyzz_xzzzzz_1,  \
                             ta_yyzz_yyyyyz_0,  \
                             ta_yyzz_yyyyyz_1,  \
                             ta_yyzz_yyyyz_0,   \
                             ta_yyzz_yyyyz_1,   \
                             ta_yyzz_yyyyzz_0,  \
                             ta_yyzz_yyyyzz_1,  \
                             ta_yyzz_yyyzz_0,   \
                             ta_yyzz_yyyzz_1,   \
                             ta_yyzz_yyyzzz_0,  \
                             ta_yyzz_yyyzzz_1,  \
                             ta_yyzz_yyzzz_0,   \
                             ta_yyzz_yyzzz_1,   \
                             ta_yyzz_yyzzzz_0,  \
                             ta_yyzz_yyzzzz_1,  \
                             ta_yyzz_yzzzz_0,   \
                             ta_yyzz_yzzzz_1,   \
                             ta_yyzz_yzzzzz_0,  \
                             ta_yyzz_yzzzzz_1,  \
                             ta_yyzz_zzzzz_0,   \
                             ta_yyzz_zzzzz_1,   \
                             ta_yyzz_zzzzzz_0,  \
                             ta_yyzz_zzzzzz_1,  \
                             ta_yzz_xxxxxx_0,   \
                             ta_yzz_xxxxxx_1,   \
                             ta_yzz_xxxxxz_0,   \
                             ta_yzz_xxxxxz_1,   \
                             ta_yzz_xxxxyz_0,   \
                             ta_yzz_xxxxyz_1,   \
                             ta_yzz_xxxxzz_0,   \
                             ta_yzz_xxxxzz_1,   \
                             ta_yzz_xxxyyz_0,   \
                             ta_yzz_xxxyyz_1,   \
                             ta_yzz_xxxyzz_0,   \
                             ta_yzz_xxxyzz_1,   \
                             ta_yzz_xxxzzz_0,   \
                             ta_yzz_xxxzzz_1,   \
                             ta_yzz_xxyyyz_0,   \
                             ta_yzz_xxyyyz_1,   \
                             ta_yzz_xxyyzz_0,   \
                             ta_yzz_xxyyzz_1,   \
                             ta_yzz_xxyzzz_0,   \
                             ta_yzz_xxyzzz_1,   \
                             ta_yzz_xxzzzz_0,   \
                             ta_yzz_xxzzzz_1,   \
                             ta_yzz_xyyyyz_0,   \
                             ta_yzz_xyyyyz_1,   \
                             ta_yzz_xyyyzz_0,   \
                             ta_yzz_xyyyzz_1,   \
                             ta_yzz_xyyzzz_0,   \
                             ta_yzz_xyyzzz_1,   \
                             ta_yzz_xyzzzz_0,   \
                             ta_yzz_xyzzzz_1,   \
                             ta_yzz_xzzzzz_0,   \
                             ta_yzz_xzzzzz_1,   \
                             ta_yzz_yyyyyz_0,   \
                             ta_yzz_yyyyyz_1,   \
                             ta_yzz_yyyyzz_0,   \
                             ta_yzz_yyyyzz_1,   \
                             ta_yzz_yyyzzz_0,   \
                             ta_yzz_yyyzzz_1,   \
                             ta_yzz_yyzzzz_0,   \
                             ta_yzz_yyzzzz_1,   \
                             ta_yzz_yzzzzz_0,   \
                             ta_yzz_yzzzzz_1,   \
                             ta_yzz_zzzzzz_0,   \
                             ta_yzz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzz_xxxxxx_0[i] =
            2.0 * ta_yzz_xxxxxx_0[i] * fe_0 - 2.0 * ta_yzz_xxxxxx_1[i] * fe_0 + ta_yyzz_xxxxxx_0[i] * pa_y[i] - ta_yyzz_xxxxxx_1[i] * pc_y[i];

        ta_yyyzz_xxxxxy_0[i] = ta_yyy_xxxxxy_0[i] * fe_0 - ta_yyy_xxxxxy_1[i] * fe_0 + ta_yyyz_xxxxxy_0[i] * pa_z[i] - ta_yyyz_xxxxxy_1[i] * pc_z[i];

        ta_yyyzz_xxxxxz_0[i] =
            2.0 * ta_yzz_xxxxxz_0[i] * fe_0 - 2.0 * ta_yzz_xxxxxz_1[i] * fe_0 + ta_yyzz_xxxxxz_0[i] * pa_y[i] - ta_yyzz_xxxxxz_1[i] * pc_y[i];

        ta_yyyzz_xxxxyy_0[i] = ta_yyy_xxxxyy_0[i] * fe_0 - ta_yyy_xxxxyy_1[i] * fe_0 + ta_yyyz_xxxxyy_0[i] * pa_z[i] - ta_yyyz_xxxxyy_1[i] * pc_z[i];

        ta_yyyzz_xxxxyz_0[i] = 2.0 * ta_yzz_xxxxyz_0[i] * fe_0 - 2.0 * ta_yzz_xxxxyz_1[i] * fe_0 + ta_yyzz_xxxxz_0[i] * fe_0 -
                               ta_yyzz_xxxxz_1[i] * fe_0 + ta_yyzz_xxxxyz_0[i] * pa_y[i] - ta_yyzz_xxxxyz_1[i] * pc_y[i];

        ta_yyyzz_xxxxzz_0[i] =
            2.0 * ta_yzz_xxxxzz_0[i] * fe_0 - 2.0 * ta_yzz_xxxxzz_1[i] * fe_0 + ta_yyzz_xxxxzz_0[i] * pa_y[i] - ta_yyzz_xxxxzz_1[i] * pc_y[i];

        ta_yyyzz_xxxyyy_0[i] = ta_yyy_xxxyyy_0[i] * fe_0 - ta_yyy_xxxyyy_1[i] * fe_0 + ta_yyyz_xxxyyy_0[i] * pa_z[i] - ta_yyyz_xxxyyy_1[i] * pc_z[i];

        ta_yyyzz_xxxyyz_0[i] = 2.0 * ta_yzz_xxxyyz_0[i] * fe_0 - 2.0 * ta_yzz_xxxyyz_1[i] * fe_0 + 2.0 * ta_yyzz_xxxyz_0[i] * fe_0 -
                               2.0 * ta_yyzz_xxxyz_1[i] * fe_0 + ta_yyzz_xxxyyz_0[i] * pa_y[i] - ta_yyzz_xxxyyz_1[i] * pc_y[i];

        ta_yyyzz_xxxyzz_0[i] = 2.0 * ta_yzz_xxxyzz_0[i] * fe_0 - 2.0 * ta_yzz_xxxyzz_1[i] * fe_0 + ta_yyzz_xxxzz_0[i] * fe_0 -
                               ta_yyzz_xxxzz_1[i] * fe_0 + ta_yyzz_xxxyzz_0[i] * pa_y[i] - ta_yyzz_xxxyzz_1[i] * pc_y[i];

        ta_yyyzz_xxxzzz_0[i] =
            2.0 * ta_yzz_xxxzzz_0[i] * fe_0 - 2.0 * ta_yzz_xxxzzz_1[i] * fe_0 + ta_yyzz_xxxzzz_0[i] * pa_y[i] - ta_yyzz_xxxzzz_1[i] * pc_y[i];

        ta_yyyzz_xxyyyy_0[i] = ta_yyy_xxyyyy_0[i] * fe_0 - ta_yyy_xxyyyy_1[i] * fe_0 + ta_yyyz_xxyyyy_0[i] * pa_z[i] - ta_yyyz_xxyyyy_1[i] * pc_z[i];

        ta_yyyzz_xxyyyz_0[i] = 2.0 * ta_yzz_xxyyyz_0[i] * fe_0 - 2.0 * ta_yzz_xxyyyz_1[i] * fe_0 + 3.0 * ta_yyzz_xxyyz_0[i] * fe_0 -
                               3.0 * ta_yyzz_xxyyz_1[i] * fe_0 + ta_yyzz_xxyyyz_0[i] * pa_y[i] - ta_yyzz_xxyyyz_1[i] * pc_y[i];

        ta_yyyzz_xxyyzz_0[i] = 2.0 * ta_yzz_xxyyzz_0[i] * fe_0 - 2.0 * ta_yzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_yyzz_xxyzz_0[i] * fe_0 -
                               2.0 * ta_yyzz_xxyzz_1[i] * fe_0 + ta_yyzz_xxyyzz_0[i] * pa_y[i] - ta_yyzz_xxyyzz_1[i] * pc_y[i];

        ta_yyyzz_xxyzzz_0[i] = 2.0 * ta_yzz_xxyzzz_0[i] * fe_0 - 2.0 * ta_yzz_xxyzzz_1[i] * fe_0 + ta_yyzz_xxzzz_0[i] * fe_0 -
                               ta_yyzz_xxzzz_1[i] * fe_0 + ta_yyzz_xxyzzz_0[i] * pa_y[i] - ta_yyzz_xxyzzz_1[i] * pc_y[i];

        ta_yyyzz_xxzzzz_0[i] =
            2.0 * ta_yzz_xxzzzz_0[i] * fe_0 - 2.0 * ta_yzz_xxzzzz_1[i] * fe_0 + ta_yyzz_xxzzzz_0[i] * pa_y[i] - ta_yyzz_xxzzzz_1[i] * pc_y[i];

        ta_yyyzz_xyyyyy_0[i] = ta_yyy_xyyyyy_0[i] * fe_0 - ta_yyy_xyyyyy_1[i] * fe_0 + ta_yyyz_xyyyyy_0[i] * pa_z[i] - ta_yyyz_xyyyyy_1[i] * pc_z[i];

        ta_yyyzz_xyyyyz_0[i] = 2.0 * ta_yzz_xyyyyz_0[i] * fe_0 - 2.0 * ta_yzz_xyyyyz_1[i] * fe_0 + 4.0 * ta_yyzz_xyyyz_0[i] * fe_0 -
                               4.0 * ta_yyzz_xyyyz_1[i] * fe_0 + ta_yyzz_xyyyyz_0[i] * pa_y[i] - ta_yyzz_xyyyyz_1[i] * pc_y[i];

        ta_yyyzz_xyyyzz_0[i] = 2.0 * ta_yzz_xyyyzz_0[i] * fe_0 - 2.0 * ta_yzz_xyyyzz_1[i] * fe_0 + 3.0 * ta_yyzz_xyyzz_0[i] * fe_0 -
                               3.0 * ta_yyzz_xyyzz_1[i] * fe_0 + ta_yyzz_xyyyzz_0[i] * pa_y[i] - ta_yyzz_xyyyzz_1[i] * pc_y[i];

        ta_yyyzz_xyyzzz_0[i] = 2.0 * ta_yzz_xyyzzz_0[i] * fe_0 - 2.0 * ta_yzz_xyyzzz_1[i] * fe_0 + 2.0 * ta_yyzz_xyzzz_0[i] * fe_0 -
                               2.0 * ta_yyzz_xyzzz_1[i] * fe_0 + ta_yyzz_xyyzzz_0[i] * pa_y[i] - ta_yyzz_xyyzzz_1[i] * pc_y[i];

        ta_yyyzz_xyzzzz_0[i] = 2.0 * ta_yzz_xyzzzz_0[i] * fe_0 - 2.0 * ta_yzz_xyzzzz_1[i] * fe_0 + ta_yyzz_xzzzz_0[i] * fe_0 -
                               ta_yyzz_xzzzz_1[i] * fe_0 + ta_yyzz_xyzzzz_0[i] * pa_y[i] - ta_yyzz_xyzzzz_1[i] * pc_y[i];

        ta_yyyzz_xzzzzz_0[i] =
            2.0 * ta_yzz_xzzzzz_0[i] * fe_0 - 2.0 * ta_yzz_xzzzzz_1[i] * fe_0 + ta_yyzz_xzzzzz_0[i] * pa_y[i] - ta_yyzz_xzzzzz_1[i] * pc_y[i];

        ta_yyyzz_yyyyyy_0[i] = ta_yyy_yyyyyy_0[i] * fe_0 - ta_yyy_yyyyyy_1[i] * fe_0 + ta_yyyz_yyyyyy_0[i] * pa_z[i] - ta_yyyz_yyyyyy_1[i] * pc_z[i];

        ta_yyyzz_yyyyyz_0[i] = 2.0 * ta_yzz_yyyyyz_0[i] * fe_0 - 2.0 * ta_yzz_yyyyyz_1[i] * fe_0 + 5.0 * ta_yyzz_yyyyz_0[i] * fe_0 -
                               5.0 * ta_yyzz_yyyyz_1[i] * fe_0 + ta_yyzz_yyyyyz_0[i] * pa_y[i] - ta_yyzz_yyyyyz_1[i] * pc_y[i];

        ta_yyyzz_yyyyzz_0[i] = 2.0 * ta_yzz_yyyyzz_0[i] * fe_0 - 2.0 * ta_yzz_yyyyzz_1[i] * fe_0 + 4.0 * ta_yyzz_yyyzz_0[i] * fe_0 -
                               4.0 * ta_yyzz_yyyzz_1[i] * fe_0 + ta_yyzz_yyyyzz_0[i] * pa_y[i] - ta_yyzz_yyyyzz_1[i] * pc_y[i];

        ta_yyyzz_yyyzzz_0[i] = 2.0 * ta_yzz_yyyzzz_0[i] * fe_0 - 2.0 * ta_yzz_yyyzzz_1[i] * fe_0 + 3.0 * ta_yyzz_yyzzz_0[i] * fe_0 -
                               3.0 * ta_yyzz_yyzzz_1[i] * fe_0 + ta_yyzz_yyyzzz_0[i] * pa_y[i] - ta_yyzz_yyyzzz_1[i] * pc_y[i];

        ta_yyyzz_yyzzzz_0[i] = 2.0 * ta_yzz_yyzzzz_0[i] * fe_0 - 2.0 * ta_yzz_yyzzzz_1[i] * fe_0 + 2.0 * ta_yyzz_yzzzz_0[i] * fe_0 -
                               2.0 * ta_yyzz_yzzzz_1[i] * fe_0 + ta_yyzz_yyzzzz_0[i] * pa_y[i] - ta_yyzz_yyzzzz_1[i] * pc_y[i];

        ta_yyyzz_yzzzzz_0[i] = 2.0 * ta_yzz_yzzzzz_0[i] * fe_0 - 2.0 * ta_yzz_yzzzzz_1[i] * fe_0 + ta_yyzz_zzzzz_0[i] * fe_0 -
                               ta_yyzz_zzzzz_1[i] * fe_0 + ta_yyzz_yzzzzz_0[i] * pa_y[i] - ta_yyzz_yzzzzz_1[i] * pc_y[i];

        ta_yyyzz_zzzzzz_0[i] =
            2.0 * ta_yzz_zzzzzz_0[i] * fe_0 - 2.0 * ta_yzz_zzzzzz_1[i] * fe_0 + ta_yyzz_zzzzzz_0[i] * pa_y[i] - ta_yyzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 504-532 components of targeted buffer : HI

    auto ta_yyzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 504);

    auto ta_yyzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 505);

    auto ta_yyzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 506);

    auto ta_yyzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 507);

    auto ta_yyzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 508);

    auto ta_yyzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 509);

    auto ta_yyzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 510);

    auto ta_yyzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 511);

    auto ta_yyzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 512);

    auto ta_yyzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 513);

    auto ta_yyzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 514);

    auto ta_yyzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 515);

    auto ta_yyzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 516);

    auto ta_yyzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 517);

    auto ta_yyzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 518);

    auto ta_yyzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 519);

    auto ta_yyzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 520);

    auto ta_yyzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 521);

    auto ta_yyzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 522);

    auto ta_yyzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 523);

    auto ta_yyzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 524);

    auto ta_yyzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 525);

    auto ta_yyzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 526);

    auto ta_yyzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 527);

    auto ta_yyzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 528);

    auto ta_yyzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 529);

    auto ta_yyzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 530);

    auto ta_yyzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 531);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta_yyz_xxxxxy_0,   \
                             ta_yyz_xxxxxy_1,   \
                             ta_yyz_xxxxyy_0,   \
                             ta_yyz_xxxxyy_1,   \
                             ta_yyz_xxxyyy_0,   \
                             ta_yyz_xxxyyy_1,   \
                             ta_yyz_xxyyyy_0,   \
                             ta_yyz_xxyyyy_1,   \
                             ta_yyz_xyyyyy_0,   \
                             ta_yyz_xyyyyy_1,   \
                             ta_yyz_yyyyyy_0,   \
                             ta_yyz_yyyyyy_1,   \
                             ta_yyzz_xxxxxy_0,  \
                             ta_yyzz_xxxxxy_1,  \
                             ta_yyzz_xxxxyy_0,  \
                             ta_yyzz_xxxxyy_1,  \
                             ta_yyzz_xxxyyy_0,  \
                             ta_yyzz_xxxyyy_1,  \
                             ta_yyzz_xxyyyy_0,  \
                             ta_yyzz_xxyyyy_1,  \
                             ta_yyzz_xyyyyy_0,  \
                             ta_yyzz_xyyyyy_1,  \
                             ta_yyzz_yyyyyy_0,  \
                             ta_yyzz_yyyyyy_1,  \
                             ta_yyzzz_xxxxxx_0, \
                             ta_yyzzz_xxxxxy_0, \
                             ta_yyzzz_xxxxxz_0, \
                             ta_yyzzz_xxxxyy_0, \
                             ta_yyzzz_xxxxyz_0, \
                             ta_yyzzz_xxxxzz_0, \
                             ta_yyzzz_xxxyyy_0, \
                             ta_yyzzz_xxxyyz_0, \
                             ta_yyzzz_xxxyzz_0, \
                             ta_yyzzz_xxxzzz_0, \
                             ta_yyzzz_xxyyyy_0, \
                             ta_yyzzz_xxyyyz_0, \
                             ta_yyzzz_xxyyzz_0, \
                             ta_yyzzz_xxyzzz_0, \
                             ta_yyzzz_xxzzzz_0, \
                             ta_yyzzz_xyyyyy_0, \
                             ta_yyzzz_xyyyyz_0, \
                             ta_yyzzz_xyyyzz_0, \
                             ta_yyzzz_xyyzzz_0, \
                             ta_yyzzz_xyzzzz_0, \
                             ta_yyzzz_xzzzzz_0, \
                             ta_yyzzz_yyyyyy_0, \
                             ta_yyzzz_yyyyyz_0, \
                             ta_yyzzz_yyyyzz_0, \
                             ta_yyzzz_yyyzzz_0, \
                             ta_yyzzz_yyzzzz_0, \
                             ta_yyzzz_yzzzzz_0, \
                             ta_yyzzz_zzzzzz_0, \
                             ta_yzzz_xxxxxx_0,  \
                             ta_yzzz_xxxxxx_1,  \
                             ta_yzzz_xxxxxz_0,  \
                             ta_yzzz_xxxxxz_1,  \
                             ta_yzzz_xxxxyz_0,  \
                             ta_yzzz_xxxxyz_1,  \
                             ta_yzzz_xxxxz_0,   \
                             ta_yzzz_xxxxz_1,   \
                             ta_yzzz_xxxxzz_0,  \
                             ta_yzzz_xxxxzz_1,  \
                             ta_yzzz_xxxyyz_0,  \
                             ta_yzzz_xxxyyz_1,  \
                             ta_yzzz_xxxyz_0,   \
                             ta_yzzz_xxxyz_1,   \
                             ta_yzzz_xxxyzz_0,  \
                             ta_yzzz_xxxyzz_1,  \
                             ta_yzzz_xxxzz_0,   \
                             ta_yzzz_xxxzz_1,   \
                             ta_yzzz_xxxzzz_0,  \
                             ta_yzzz_xxxzzz_1,  \
                             ta_yzzz_xxyyyz_0,  \
                             ta_yzzz_xxyyyz_1,  \
                             ta_yzzz_xxyyz_0,   \
                             ta_yzzz_xxyyz_1,   \
                             ta_yzzz_xxyyzz_0,  \
                             ta_yzzz_xxyyzz_1,  \
                             ta_yzzz_xxyzz_0,   \
                             ta_yzzz_xxyzz_1,   \
                             ta_yzzz_xxyzzz_0,  \
                             ta_yzzz_xxyzzz_1,  \
                             ta_yzzz_xxzzz_0,   \
                             ta_yzzz_xxzzz_1,   \
                             ta_yzzz_xxzzzz_0,  \
                             ta_yzzz_xxzzzz_1,  \
                             ta_yzzz_xyyyyz_0,  \
                             ta_yzzz_xyyyyz_1,  \
                             ta_yzzz_xyyyz_0,   \
                             ta_yzzz_xyyyz_1,   \
                             ta_yzzz_xyyyzz_0,  \
                             ta_yzzz_xyyyzz_1,  \
                             ta_yzzz_xyyzz_0,   \
                             ta_yzzz_xyyzz_1,   \
                             ta_yzzz_xyyzzz_0,  \
                             ta_yzzz_xyyzzz_1,  \
                             ta_yzzz_xyzzz_0,   \
                             ta_yzzz_xyzzz_1,   \
                             ta_yzzz_xyzzzz_0,  \
                             ta_yzzz_xyzzzz_1,  \
                             ta_yzzz_xzzzz_0,   \
                             ta_yzzz_xzzzz_1,   \
                             ta_yzzz_xzzzzz_0,  \
                             ta_yzzz_xzzzzz_1,  \
                             ta_yzzz_yyyyyz_0,  \
                             ta_yzzz_yyyyyz_1,  \
                             ta_yzzz_yyyyz_0,   \
                             ta_yzzz_yyyyz_1,   \
                             ta_yzzz_yyyyzz_0,  \
                             ta_yzzz_yyyyzz_1,  \
                             ta_yzzz_yyyzz_0,   \
                             ta_yzzz_yyyzz_1,   \
                             ta_yzzz_yyyzzz_0,  \
                             ta_yzzz_yyyzzz_1,  \
                             ta_yzzz_yyzzz_0,   \
                             ta_yzzz_yyzzz_1,   \
                             ta_yzzz_yyzzzz_0,  \
                             ta_yzzz_yyzzzz_1,  \
                             ta_yzzz_yzzzz_0,   \
                             ta_yzzz_yzzzz_1,   \
                             ta_yzzz_yzzzzz_0,  \
                             ta_yzzz_yzzzzz_1,  \
                             ta_yzzz_zzzzz_0,   \
                             ta_yzzz_zzzzz_1,   \
                             ta_yzzz_zzzzzz_0,  \
                             ta_yzzz_zzzzzz_1,  \
                             ta_zzz_xxxxxx_0,   \
                             ta_zzz_xxxxxx_1,   \
                             ta_zzz_xxxxxz_0,   \
                             ta_zzz_xxxxxz_1,   \
                             ta_zzz_xxxxyz_0,   \
                             ta_zzz_xxxxyz_1,   \
                             ta_zzz_xxxxzz_0,   \
                             ta_zzz_xxxxzz_1,   \
                             ta_zzz_xxxyyz_0,   \
                             ta_zzz_xxxyyz_1,   \
                             ta_zzz_xxxyzz_0,   \
                             ta_zzz_xxxyzz_1,   \
                             ta_zzz_xxxzzz_0,   \
                             ta_zzz_xxxzzz_1,   \
                             ta_zzz_xxyyyz_0,   \
                             ta_zzz_xxyyyz_1,   \
                             ta_zzz_xxyyzz_0,   \
                             ta_zzz_xxyyzz_1,   \
                             ta_zzz_xxyzzz_0,   \
                             ta_zzz_xxyzzz_1,   \
                             ta_zzz_xxzzzz_0,   \
                             ta_zzz_xxzzzz_1,   \
                             ta_zzz_xyyyyz_0,   \
                             ta_zzz_xyyyyz_1,   \
                             ta_zzz_xyyyzz_0,   \
                             ta_zzz_xyyyzz_1,   \
                             ta_zzz_xyyzzz_0,   \
                             ta_zzz_xyyzzz_1,   \
                             ta_zzz_xyzzzz_0,   \
                             ta_zzz_xyzzzz_1,   \
                             ta_zzz_xzzzzz_0,   \
                             ta_zzz_xzzzzz_1,   \
                             ta_zzz_yyyyyz_0,   \
                             ta_zzz_yyyyyz_1,   \
                             ta_zzz_yyyyzz_0,   \
                             ta_zzz_yyyyzz_1,   \
                             ta_zzz_yyyzzz_0,   \
                             ta_zzz_yyyzzz_1,   \
                             ta_zzz_yyzzzz_0,   \
                             ta_zzz_yyzzzz_1,   \
                             ta_zzz_yzzzzz_0,   \
                             ta_zzz_yzzzzz_1,   \
                             ta_zzz_zzzzzz_0,   \
                             ta_zzz_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzz_xxxxxx_0[i] = ta_zzz_xxxxxx_0[i] * fe_0 - ta_zzz_xxxxxx_1[i] * fe_0 + ta_yzzz_xxxxxx_0[i] * pa_y[i] - ta_yzzz_xxxxxx_1[i] * pc_y[i];

        ta_yyzzz_xxxxxy_0[i] =
            2.0 * ta_yyz_xxxxxy_0[i] * fe_0 - 2.0 * ta_yyz_xxxxxy_1[i] * fe_0 + ta_yyzz_xxxxxy_0[i] * pa_z[i] - ta_yyzz_xxxxxy_1[i] * pc_z[i];

        ta_yyzzz_xxxxxz_0[i] = ta_zzz_xxxxxz_0[i] * fe_0 - ta_zzz_xxxxxz_1[i] * fe_0 + ta_yzzz_xxxxxz_0[i] * pa_y[i] - ta_yzzz_xxxxxz_1[i] * pc_y[i];

        ta_yyzzz_xxxxyy_0[i] =
            2.0 * ta_yyz_xxxxyy_0[i] * fe_0 - 2.0 * ta_yyz_xxxxyy_1[i] * fe_0 + ta_yyzz_xxxxyy_0[i] * pa_z[i] - ta_yyzz_xxxxyy_1[i] * pc_z[i];

        ta_yyzzz_xxxxyz_0[i] = ta_zzz_xxxxyz_0[i] * fe_0 - ta_zzz_xxxxyz_1[i] * fe_0 + ta_yzzz_xxxxz_0[i] * fe_0 - ta_yzzz_xxxxz_1[i] * fe_0 +
                               ta_yzzz_xxxxyz_0[i] * pa_y[i] - ta_yzzz_xxxxyz_1[i] * pc_y[i];

        ta_yyzzz_xxxxzz_0[i] = ta_zzz_xxxxzz_0[i] * fe_0 - ta_zzz_xxxxzz_1[i] * fe_0 + ta_yzzz_xxxxzz_0[i] * pa_y[i] - ta_yzzz_xxxxzz_1[i] * pc_y[i];

        ta_yyzzz_xxxyyy_0[i] =
            2.0 * ta_yyz_xxxyyy_0[i] * fe_0 - 2.0 * ta_yyz_xxxyyy_1[i] * fe_0 + ta_yyzz_xxxyyy_0[i] * pa_z[i] - ta_yyzz_xxxyyy_1[i] * pc_z[i];

        ta_yyzzz_xxxyyz_0[i] = ta_zzz_xxxyyz_0[i] * fe_0 - ta_zzz_xxxyyz_1[i] * fe_0 + 2.0 * ta_yzzz_xxxyz_0[i] * fe_0 -
                               2.0 * ta_yzzz_xxxyz_1[i] * fe_0 + ta_yzzz_xxxyyz_0[i] * pa_y[i] - ta_yzzz_xxxyyz_1[i] * pc_y[i];

        ta_yyzzz_xxxyzz_0[i] = ta_zzz_xxxyzz_0[i] * fe_0 - ta_zzz_xxxyzz_1[i] * fe_0 + ta_yzzz_xxxzz_0[i] * fe_0 - ta_yzzz_xxxzz_1[i] * fe_0 +
                               ta_yzzz_xxxyzz_0[i] * pa_y[i] - ta_yzzz_xxxyzz_1[i] * pc_y[i];

        ta_yyzzz_xxxzzz_0[i] = ta_zzz_xxxzzz_0[i] * fe_0 - ta_zzz_xxxzzz_1[i] * fe_0 + ta_yzzz_xxxzzz_0[i] * pa_y[i] - ta_yzzz_xxxzzz_1[i] * pc_y[i];

        ta_yyzzz_xxyyyy_0[i] =
            2.0 * ta_yyz_xxyyyy_0[i] * fe_0 - 2.0 * ta_yyz_xxyyyy_1[i] * fe_0 + ta_yyzz_xxyyyy_0[i] * pa_z[i] - ta_yyzz_xxyyyy_1[i] * pc_z[i];

        ta_yyzzz_xxyyyz_0[i] = ta_zzz_xxyyyz_0[i] * fe_0 - ta_zzz_xxyyyz_1[i] * fe_0 + 3.0 * ta_yzzz_xxyyz_0[i] * fe_0 -
                               3.0 * ta_yzzz_xxyyz_1[i] * fe_0 + ta_yzzz_xxyyyz_0[i] * pa_y[i] - ta_yzzz_xxyyyz_1[i] * pc_y[i];

        ta_yyzzz_xxyyzz_0[i] = ta_zzz_xxyyzz_0[i] * fe_0 - ta_zzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_yzzz_xxyzz_0[i] * fe_0 -
                               2.0 * ta_yzzz_xxyzz_1[i] * fe_0 + ta_yzzz_xxyyzz_0[i] * pa_y[i] - ta_yzzz_xxyyzz_1[i] * pc_y[i];

        ta_yyzzz_xxyzzz_0[i] = ta_zzz_xxyzzz_0[i] * fe_0 - ta_zzz_xxyzzz_1[i] * fe_0 + ta_yzzz_xxzzz_0[i] * fe_0 - ta_yzzz_xxzzz_1[i] * fe_0 +
                               ta_yzzz_xxyzzz_0[i] * pa_y[i] - ta_yzzz_xxyzzz_1[i] * pc_y[i];

        ta_yyzzz_xxzzzz_0[i] = ta_zzz_xxzzzz_0[i] * fe_0 - ta_zzz_xxzzzz_1[i] * fe_0 + ta_yzzz_xxzzzz_0[i] * pa_y[i] - ta_yzzz_xxzzzz_1[i] * pc_y[i];

        ta_yyzzz_xyyyyy_0[i] =
            2.0 * ta_yyz_xyyyyy_0[i] * fe_0 - 2.0 * ta_yyz_xyyyyy_1[i] * fe_0 + ta_yyzz_xyyyyy_0[i] * pa_z[i] - ta_yyzz_xyyyyy_1[i] * pc_z[i];

        ta_yyzzz_xyyyyz_0[i] = ta_zzz_xyyyyz_0[i] * fe_0 - ta_zzz_xyyyyz_1[i] * fe_0 + 4.0 * ta_yzzz_xyyyz_0[i] * fe_0 -
                               4.0 * ta_yzzz_xyyyz_1[i] * fe_0 + ta_yzzz_xyyyyz_0[i] * pa_y[i] - ta_yzzz_xyyyyz_1[i] * pc_y[i];

        ta_yyzzz_xyyyzz_0[i] = ta_zzz_xyyyzz_0[i] * fe_0 - ta_zzz_xyyyzz_1[i] * fe_0 + 3.0 * ta_yzzz_xyyzz_0[i] * fe_0 -
                               3.0 * ta_yzzz_xyyzz_1[i] * fe_0 + ta_yzzz_xyyyzz_0[i] * pa_y[i] - ta_yzzz_xyyyzz_1[i] * pc_y[i];

        ta_yyzzz_xyyzzz_0[i] = ta_zzz_xyyzzz_0[i] * fe_0 - ta_zzz_xyyzzz_1[i] * fe_0 + 2.0 * ta_yzzz_xyzzz_0[i] * fe_0 -
                               2.0 * ta_yzzz_xyzzz_1[i] * fe_0 + ta_yzzz_xyyzzz_0[i] * pa_y[i] - ta_yzzz_xyyzzz_1[i] * pc_y[i];

        ta_yyzzz_xyzzzz_0[i] = ta_zzz_xyzzzz_0[i] * fe_0 - ta_zzz_xyzzzz_1[i] * fe_0 + ta_yzzz_xzzzz_0[i] * fe_0 - ta_yzzz_xzzzz_1[i] * fe_0 +
                               ta_yzzz_xyzzzz_0[i] * pa_y[i] - ta_yzzz_xyzzzz_1[i] * pc_y[i];

        ta_yyzzz_xzzzzz_0[i] = ta_zzz_xzzzzz_0[i] * fe_0 - ta_zzz_xzzzzz_1[i] * fe_0 + ta_yzzz_xzzzzz_0[i] * pa_y[i] - ta_yzzz_xzzzzz_1[i] * pc_y[i];

        ta_yyzzz_yyyyyy_0[i] =
            2.0 * ta_yyz_yyyyyy_0[i] * fe_0 - 2.0 * ta_yyz_yyyyyy_1[i] * fe_0 + ta_yyzz_yyyyyy_0[i] * pa_z[i] - ta_yyzz_yyyyyy_1[i] * pc_z[i];

        ta_yyzzz_yyyyyz_0[i] = ta_zzz_yyyyyz_0[i] * fe_0 - ta_zzz_yyyyyz_1[i] * fe_0 + 5.0 * ta_yzzz_yyyyz_0[i] * fe_0 -
                               5.0 * ta_yzzz_yyyyz_1[i] * fe_0 + ta_yzzz_yyyyyz_0[i] * pa_y[i] - ta_yzzz_yyyyyz_1[i] * pc_y[i];

        ta_yyzzz_yyyyzz_0[i] = ta_zzz_yyyyzz_0[i] * fe_0 - ta_zzz_yyyyzz_1[i] * fe_0 + 4.0 * ta_yzzz_yyyzz_0[i] * fe_0 -
                               4.0 * ta_yzzz_yyyzz_1[i] * fe_0 + ta_yzzz_yyyyzz_0[i] * pa_y[i] - ta_yzzz_yyyyzz_1[i] * pc_y[i];

        ta_yyzzz_yyyzzz_0[i] = ta_zzz_yyyzzz_0[i] * fe_0 - ta_zzz_yyyzzz_1[i] * fe_0 + 3.0 * ta_yzzz_yyzzz_0[i] * fe_0 -
                               3.0 * ta_yzzz_yyzzz_1[i] * fe_0 + ta_yzzz_yyyzzz_0[i] * pa_y[i] - ta_yzzz_yyyzzz_1[i] * pc_y[i];

        ta_yyzzz_yyzzzz_0[i] = ta_zzz_yyzzzz_0[i] * fe_0 - ta_zzz_yyzzzz_1[i] * fe_0 + 2.0 * ta_yzzz_yzzzz_0[i] * fe_0 -
                               2.0 * ta_yzzz_yzzzz_1[i] * fe_0 + ta_yzzz_yyzzzz_0[i] * pa_y[i] - ta_yzzz_yyzzzz_1[i] * pc_y[i];

        ta_yyzzz_yzzzzz_0[i] = ta_zzz_yzzzzz_0[i] * fe_0 - ta_zzz_yzzzzz_1[i] * fe_0 + ta_yzzz_zzzzz_0[i] * fe_0 - ta_yzzz_zzzzz_1[i] * fe_0 +
                               ta_yzzz_yzzzzz_0[i] * pa_y[i] - ta_yzzz_yzzzzz_1[i] * pc_y[i];

        ta_yyzzz_zzzzzz_0[i] = ta_zzz_zzzzzz_0[i] * fe_0 - ta_zzz_zzzzzz_1[i] * fe_0 + ta_yzzz_zzzzzz_0[i] * pa_y[i] - ta_yzzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 532-560 components of targeted buffer : HI

    auto ta_yzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 532);

    auto ta_yzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 533);

    auto ta_yzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 534);

    auto ta_yzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 535);

    auto ta_yzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 536);

    auto ta_yzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 537);

    auto ta_yzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 538);

    auto ta_yzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 539);

    auto ta_yzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 540);

    auto ta_yzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 541);

    auto ta_yzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 542);

    auto ta_yzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 543);

    auto ta_yzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 544);

    auto ta_yzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 545);

    auto ta_yzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 546);

    auto ta_yzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 547);

    auto ta_yzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 548);

    auto ta_yzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 549);

    auto ta_yzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 550);

    auto ta_yzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 551);

    auto ta_yzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 552);

    auto ta_yzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 553);

    auto ta_yzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 554);

    auto ta_yzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 555);

    auto ta_yzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 556);

    auto ta_yzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 557);

    auto ta_yzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 558);

    auto ta_yzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 559);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta_yzzzz_xxxxxx_0, \
                             ta_yzzzz_xxxxxy_0, \
                             ta_yzzzz_xxxxxz_0, \
                             ta_yzzzz_xxxxyy_0, \
                             ta_yzzzz_xxxxyz_0, \
                             ta_yzzzz_xxxxzz_0, \
                             ta_yzzzz_xxxyyy_0, \
                             ta_yzzzz_xxxyyz_0, \
                             ta_yzzzz_xxxyzz_0, \
                             ta_yzzzz_xxxzzz_0, \
                             ta_yzzzz_xxyyyy_0, \
                             ta_yzzzz_xxyyyz_0, \
                             ta_yzzzz_xxyyzz_0, \
                             ta_yzzzz_xxyzzz_0, \
                             ta_yzzzz_xxzzzz_0, \
                             ta_yzzzz_xyyyyy_0, \
                             ta_yzzzz_xyyyyz_0, \
                             ta_yzzzz_xyyyzz_0, \
                             ta_yzzzz_xyyzzz_0, \
                             ta_yzzzz_xyzzzz_0, \
                             ta_yzzzz_xzzzzz_0, \
                             ta_yzzzz_yyyyyy_0, \
                             ta_yzzzz_yyyyyz_0, \
                             ta_yzzzz_yyyyzz_0, \
                             ta_yzzzz_yyyzzz_0, \
                             ta_yzzzz_yyzzzz_0, \
                             ta_yzzzz_yzzzzz_0, \
                             ta_yzzzz_zzzzzz_0, \
                             ta_zzzz_xxxxx_0,   \
                             ta_zzzz_xxxxx_1,   \
                             ta_zzzz_xxxxxx_0,  \
                             ta_zzzz_xxxxxx_1,  \
                             ta_zzzz_xxxxxy_0,  \
                             ta_zzzz_xxxxxy_1,  \
                             ta_zzzz_xxxxxz_0,  \
                             ta_zzzz_xxxxxz_1,  \
                             ta_zzzz_xxxxy_0,   \
                             ta_zzzz_xxxxy_1,   \
                             ta_zzzz_xxxxyy_0,  \
                             ta_zzzz_xxxxyy_1,  \
                             ta_zzzz_xxxxyz_0,  \
                             ta_zzzz_xxxxyz_1,  \
                             ta_zzzz_xxxxz_0,   \
                             ta_zzzz_xxxxz_1,   \
                             ta_zzzz_xxxxzz_0,  \
                             ta_zzzz_xxxxzz_1,  \
                             ta_zzzz_xxxyy_0,   \
                             ta_zzzz_xxxyy_1,   \
                             ta_zzzz_xxxyyy_0,  \
                             ta_zzzz_xxxyyy_1,  \
                             ta_zzzz_xxxyyz_0,  \
                             ta_zzzz_xxxyyz_1,  \
                             ta_zzzz_xxxyz_0,   \
                             ta_zzzz_xxxyz_1,   \
                             ta_zzzz_xxxyzz_0,  \
                             ta_zzzz_xxxyzz_1,  \
                             ta_zzzz_xxxzz_0,   \
                             ta_zzzz_xxxzz_1,   \
                             ta_zzzz_xxxzzz_0,  \
                             ta_zzzz_xxxzzz_1,  \
                             ta_zzzz_xxyyy_0,   \
                             ta_zzzz_xxyyy_1,   \
                             ta_zzzz_xxyyyy_0,  \
                             ta_zzzz_xxyyyy_1,  \
                             ta_zzzz_xxyyyz_0,  \
                             ta_zzzz_xxyyyz_1,  \
                             ta_zzzz_xxyyz_0,   \
                             ta_zzzz_xxyyz_1,   \
                             ta_zzzz_xxyyzz_0,  \
                             ta_zzzz_xxyyzz_1,  \
                             ta_zzzz_xxyzz_0,   \
                             ta_zzzz_xxyzz_1,   \
                             ta_zzzz_xxyzzz_0,  \
                             ta_zzzz_xxyzzz_1,  \
                             ta_zzzz_xxzzz_0,   \
                             ta_zzzz_xxzzz_1,   \
                             ta_zzzz_xxzzzz_0,  \
                             ta_zzzz_xxzzzz_1,  \
                             ta_zzzz_xyyyy_0,   \
                             ta_zzzz_xyyyy_1,   \
                             ta_zzzz_xyyyyy_0,  \
                             ta_zzzz_xyyyyy_1,  \
                             ta_zzzz_xyyyyz_0,  \
                             ta_zzzz_xyyyyz_1,  \
                             ta_zzzz_xyyyz_0,   \
                             ta_zzzz_xyyyz_1,   \
                             ta_zzzz_xyyyzz_0,  \
                             ta_zzzz_xyyyzz_1,  \
                             ta_zzzz_xyyzz_0,   \
                             ta_zzzz_xyyzz_1,   \
                             ta_zzzz_xyyzzz_0,  \
                             ta_zzzz_xyyzzz_1,  \
                             ta_zzzz_xyzzz_0,   \
                             ta_zzzz_xyzzz_1,   \
                             ta_zzzz_xyzzzz_0,  \
                             ta_zzzz_xyzzzz_1,  \
                             ta_zzzz_xzzzz_0,   \
                             ta_zzzz_xzzzz_1,   \
                             ta_zzzz_xzzzzz_0,  \
                             ta_zzzz_xzzzzz_1,  \
                             ta_zzzz_yyyyy_0,   \
                             ta_zzzz_yyyyy_1,   \
                             ta_zzzz_yyyyyy_0,  \
                             ta_zzzz_yyyyyy_1,  \
                             ta_zzzz_yyyyyz_0,  \
                             ta_zzzz_yyyyyz_1,  \
                             ta_zzzz_yyyyz_0,   \
                             ta_zzzz_yyyyz_1,   \
                             ta_zzzz_yyyyzz_0,  \
                             ta_zzzz_yyyyzz_1,  \
                             ta_zzzz_yyyzz_0,   \
                             ta_zzzz_yyyzz_1,   \
                             ta_zzzz_yyyzzz_0,  \
                             ta_zzzz_yyyzzz_1,  \
                             ta_zzzz_yyzzz_0,   \
                             ta_zzzz_yyzzz_1,   \
                             ta_zzzz_yyzzzz_0,  \
                             ta_zzzz_yyzzzz_1,  \
                             ta_zzzz_yzzzz_0,   \
                             ta_zzzz_yzzzz_1,   \
                             ta_zzzz_yzzzzz_0,  \
                             ta_zzzz_yzzzzz_1,  \
                             ta_zzzz_zzzzz_0,   \
                             ta_zzzz_zzzzz_1,   \
                             ta_zzzz_zzzzzz_0,  \
                             ta_zzzz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzz_xxxxxx_0[i] = ta_zzzz_xxxxxx_0[i] * pa_y[i] - ta_zzzz_xxxxxx_1[i] * pc_y[i];

        ta_yzzzz_xxxxxy_0[i] = ta_zzzz_xxxxx_0[i] * fe_0 - ta_zzzz_xxxxx_1[i] * fe_0 + ta_zzzz_xxxxxy_0[i] * pa_y[i] - ta_zzzz_xxxxxy_1[i] * pc_y[i];

        ta_yzzzz_xxxxxz_0[i] = ta_zzzz_xxxxxz_0[i] * pa_y[i] - ta_zzzz_xxxxxz_1[i] * pc_y[i];

        ta_yzzzz_xxxxyy_0[i] =
            2.0 * ta_zzzz_xxxxy_0[i] * fe_0 - 2.0 * ta_zzzz_xxxxy_1[i] * fe_0 + ta_zzzz_xxxxyy_0[i] * pa_y[i] - ta_zzzz_xxxxyy_1[i] * pc_y[i];

        ta_yzzzz_xxxxyz_0[i] = ta_zzzz_xxxxz_0[i] * fe_0 - ta_zzzz_xxxxz_1[i] * fe_0 + ta_zzzz_xxxxyz_0[i] * pa_y[i] - ta_zzzz_xxxxyz_1[i] * pc_y[i];

        ta_yzzzz_xxxxzz_0[i] = ta_zzzz_xxxxzz_0[i] * pa_y[i] - ta_zzzz_xxxxzz_1[i] * pc_y[i];

        ta_yzzzz_xxxyyy_0[i] =
            3.0 * ta_zzzz_xxxyy_0[i] * fe_0 - 3.0 * ta_zzzz_xxxyy_1[i] * fe_0 + ta_zzzz_xxxyyy_0[i] * pa_y[i] - ta_zzzz_xxxyyy_1[i] * pc_y[i];

        ta_yzzzz_xxxyyz_0[i] =
            2.0 * ta_zzzz_xxxyz_0[i] * fe_0 - 2.0 * ta_zzzz_xxxyz_1[i] * fe_0 + ta_zzzz_xxxyyz_0[i] * pa_y[i] - ta_zzzz_xxxyyz_1[i] * pc_y[i];

        ta_yzzzz_xxxyzz_0[i] = ta_zzzz_xxxzz_0[i] * fe_0 - ta_zzzz_xxxzz_1[i] * fe_0 + ta_zzzz_xxxyzz_0[i] * pa_y[i] - ta_zzzz_xxxyzz_1[i] * pc_y[i];

        ta_yzzzz_xxxzzz_0[i] = ta_zzzz_xxxzzz_0[i] * pa_y[i] - ta_zzzz_xxxzzz_1[i] * pc_y[i];

        ta_yzzzz_xxyyyy_0[i] =
            4.0 * ta_zzzz_xxyyy_0[i] * fe_0 - 4.0 * ta_zzzz_xxyyy_1[i] * fe_0 + ta_zzzz_xxyyyy_0[i] * pa_y[i] - ta_zzzz_xxyyyy_1[i] * pc_y[i];

        ta_yzzzz_xxyyyz_0[i] =
            3.0 * ta_zzzz_xxyyz_0[i] * fe_0 - 3.0 * ta_zzzz_xxyyz_1[i] * fe_0 + ta_zzzz_xxyyyz_0[i] * pa_y[i] - ta_zzzz_xxyyyz_1[i] * pc_y[i];

        ta_yzzzz_xxyyzz_0[i] =
            2.0 * ta_zzzz_xxyzz_0[i] * fe_0 - 2.0 * ta_zzzz_xxyzz_1[i] * fe_0 + ta_zzzz_xxyyzz_0[i] * pa_y[i] - ta_zzzz_xxyyzz_1[i] * pc_y[i];

        ta_yzzzz_xxyzzz_0[i] = ta_zzzz_xxzzz_0[i] * fe_0 - ta_zzzz_xxzzz_1[i] * fe_0 + ta_zzzz_xxyzzz_0[i] * pa_y[i] - ta_zzzz_xxyzzz_1[i] * pc_y[i];

        ta_yzzzz_xxzzzz_0[i] = ta_zzzz_xxzzzz_0[i] * pa_y[i] - ta_zzzz_xxzzzz_1[i] * pc_y[i];

        ta_yzzzz_xyyyyy_0[i] =
            5.0 * ta_zzzz_xyyyy_0[i] * fe_0 - 5.0 * ta_zzzz_xyyyy_1[i] * fe_0 + ta_zzzz_xyyyyy_0[i] * pa_y[i] - ta_zzzz_xyyyyy_1[i] * pc_y[i];

        ta_yzzzz_xyyyyz_0[i] =
            4.0 * ta_zzzz_xyyyz_0[i] * fe_0 - 4.0 * ta_zzzz_xyyyz_1[i] * fe_0 + ta_zzzz_xyyyyz_0[i] * pa_y[i] - ta_zzzz_xyyyyz_1[i] * pc_y[i];

        ta_yzzzz_xyyyzz_0[i] =
            3.0 * ta_zzzz_xyyzz_0[i] * fe_0 - 3.0 * ta_zzzz_xyyzz_1[i] * fe_0 + ta_zzzz_xyyyzz_0[i] * pa_y[i] - ta_zzzz_xyyyzz_1[i] * pc_y[i];

        ta_yzzzz_xyyzzz_0[i] =
            2.0 * ta_zzzz_xyzzz_0[i] * fe_0 - 2.0 * ta_zzzz_xyzzz_1[i] * fe_0 + ta_zzzz_xyyzzz_0[i] * pa_y[i] - ta_zzzz_xyyzzz_1[i] * pc_y[i];

        ta_yzzzz_xyzzzz_0[i] = ta_zzzz_xzzzz_0[i] * fe_0 - ta_zzzz_xzzzz_1[i] * fe_0 + ta_zzzz_xyzzzz_0[i] * pa_y[i] - ta_zzzz_xyzzzz_1[i] * pc_y[i];

        ta_yzzzz_xzzzzz_0[i] = ta_zzzz_xzzzzz_0[i] * pa_y[i] - ta_zzzz_xzzzzz_1[i] * pc_y[i];

        ta_yzzzz_yyyyyy_0[i] =
            6.0 * ta_zzzz_yyyyy_0[i] * fe_0 - 6.0 * ta_zzzz_yyyyy_1[i] * fe_0 + ta_zzzz_yyyyyy_0[i] * pa_y[i] - ta_zzzz_yyyyyy_1[i] * pc_y[i];

        ta_yzzzz_yyyyyz_0[i] =
            5.0 * ta_zzzz_yyyyz_0[i] * fe_0 - 5.0 * ta_zzzz_yyyyz_1[i] * fe_0 + ta_zzzz_yyyyyz_0[i] * pa_y[i] - ta_zzzz_yyyyyz_1[i] * pc_y[i];

        ta_yzzzz_yyyyzz_0[i] =
            4.0 * ta_zzzz_yyyzz_0[i] * fe_0 - 4.0 * ta_zzzz_yyyzz_1[i] * fe_0 + ta_zzzz_yyyyzz_0[i] * pa_y[i] - ta_zzzz_yyyyzz_1[i] * pc_y[i];

        ta_yzzzz_yyyzzz_0[i] =
            3.0 * ta_zzzz_yyzzz_0[i] * fe_0 - 3.0 * ta_zzzz_yyzzz_1[i] * fe_0 + ta_zzzz_yyyzzz_0[i] * pa_y[i] - ta_zzzz_yyyzzz_1[i] * pc_y[i];

        ta_yzzzz_yyzzzz_0[i] =
            2.0 * ta_zzzz_yzzzz_0[i] * fe_0 - 2.0 * ta_zzzz_yzzzz_1[i] * fe_0 + ta_zzzz_yyzzzz_0[i] * pa_y[i] - ta_zzzz_yyzzzz_1[i] * pc_y[i];

        ta_yzzzz_yzzzzz_0[i] = ta_zzzz_zzzzz_0[i] * fe_0 - ta_zzzz_zzzzz_1[i] * fe_0 + ta_zzzz_yzzzzz_0[i] * pa_y[i] - ta_zzzz_yzzzzz_1[i] * pc_y[i];

        ta_yzzzz_zzzzzz_0[i] = ta_zzzz_zzzzzz_0[i] * pa_y[i] - ta_zzzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 560-588 components of targeted buffer : HI

    auto ta_zzzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_hi + 560);

    auto ta_zzzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_hi + 561);

    auto ta_zzzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_hi + 562);

    auto ta_zzzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_hi + 563);

    auto ta_zzzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_hi + 564);

    auto ta_zzzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_hi + 565);

    auto ta_zzzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_hi + 566);

    auto ta_zzzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_hi + 567);

    auto ta_zzzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_hi + 568);

    auto ta_zzzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_hi + 569);

    auto ta_zzzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_hi + 570);

    auto ta_zzzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_hi + 571);

    auto ta_zzzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_hi + 572);

    auto ta_zzzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_hi + 573);

    auto ta_zzzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_hi + 574);

    auto ta_zzzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_hi + 575);

    auto ta_zzzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_hi + 576);

    auto ta_zzzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_hi + 577);

    auto ta_zzzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_hi + 578);

    auto ta_zzzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_hi + 579);

    auto ta_zzzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_hi + 580);

    auto ta_zzzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_hi + 581);

    auto ta_zzzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_hi + 582);

    auto ta_zzzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_hi + 583);

    auto ta_zzzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_hi + 584);

    auto ta_zzzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_hi + 585);

    auto ta_zzzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_hi + 586);

    auto ta_zzzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_hi + 587);

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta_zzz_xxxxxx_0,   \
                             ta_zzz_xxxxxx_1,   \
                             ta_zzz_xxxxxy_0,   \
                             ta_zzz_xxxxxy_1,   \
                             ta_zzz_xxxxxz_0,   \
                             ta_zzz_xxxxxz_1,   \
                             ta_zzz_xxxxyy_0,   \
                             ta_zzz_xxxxyy_1,   \
                             ta_zzz_xxxxyz_0,   \
                             ta_zzz_xxxxyz_1,   \
                             ta_zzz_xxxxzz_0,   \
                             ta_zzz_xxxxzz_1,   \
                             ta_zzz_xxxyyy_0,   \
                             ta_zzz_xxxyyy_1,   \
                             ta_zzz_xxxyyz_0,   \
                             ta_zzz_xxxyyz_1,   \
                             ta_zzz_xxxyzz_0,   \
                             ta_zzz_xxxyzz_1,   \
                             ta_zzz_xxxzzz_0,   \
                             ta_zzz_xxxzzz_1,   \
                             ta_zzz_xxyyyy_0,   \
                             ta_zzz_xxyyyy_1,   \
                             ta_zzz_xxyyyz_0,   \
                             ta_zzz_xxyyyz_1,   \
                             ta_zzz_xxyyzz_0,   \
                             ta_zzz_xxyyzz_1,   \
                             ta_zzz_xxyzzz_0,   \
                             ta_zzz_xxyzzz_1,   \
                             ta_zzz_xxzzzz_0,   \
                             ta_zzz_xxzzzz_1,   \
                             ta_zzz_xyyyyy_0,   \
                             ta_zzz_xyyyyy_1,   \
                             ta_zzz_xyyyyz_0,   \
                             ta_zzz_xyyyyz_1,   \
                             ta_zzz_xyyyzz_0,   \
                             ta_zzz_xyyyzz_1,   \
                             ta_zzz_xyyzzz_0,   \
                             ta_zzz_xyyzzz_1,   \
                             ta_zzz_xyzzzz_0,   \
                             ta_zzz_xyzzzz_1,   \
                             ta_zzz_xzzzzz_0,   \
                             ta_zzz_xzzzzz_1,   \
                             ta_zzz_yyyyyy_0,   \
                             ta_zzz_yyyyyy_1,   \
                             ta_zzz_yyyyyz_0,   \
                             ta_zzz_yyyyyz_1,   \
                             ta_zzz_yyyyzz_0,   \
                             ta_zzz_yyyyzz_1,   \
                             ta_zzz_yyyzzz_0,   \
                             ta_zzz_yyyzzz_1,   \
                             ta_zzz_yyzzzz_0,   \
                             ta_zzz_yyzzzz_1,   \
                             ta_zzz_yzzzzz_0,   \
                             ta_zzz_yzzzzz_1,   \
                             ta_zzz_zzzzzz_0,   \
                             ta_zzz_zzzzzz_1,   \
                             ta_zzzz_xxxxx_0,   \
                             ta_zzzz_xxxxx_1,   \
                             ta_zzzz_xxxxxx_0,  \
                             ta_zzzz_xxxxxx_1,  \
                             ta_zzzz_xxxxxy_0,  \
                             ta_zzzz_xxxxxy_1,  \
                             ta_zzzz_xxxxxz_0,  \
                             ta_zzzz_xxxxxz_1,  \
                             ta_zzzz_xxxxy_0,   \
                             ta_zzzz_xxxxy_1,   \
                             ta_zzzz_xxxxyy_0,  \
                             ta_zzzz_xxxxyy_1,  \
                             ta_zzzz_xxxxyz_0,  \
                             ta_zzzz_xxxxyz_1,  \
                             ta_zzzz_xxxxz_0,   \
                             ta_zzzz_xxxxz_1,   \
                             ta_zzzz_xxxxzz_0,  \
                             ta_zzzz_xxxxzz_1,  \
                             ta_zzzz_xxxyy_0,   \
                             ta_zzzz_xxxyy_1,   \
                             ta_zzzz_xxxyyy_0,  \
                             ta_zzzz_xxxyyy_1,  \
                             ta_zzzz_xxxyyz_0,  \
                             ta_zzzz_xxxyyz_1,  \
                             ta_zzzz_xxxyz_0,   \
                             ta_zzzz_xxxyz_1,   \
                             ta_zzzz_xxxyzz_0,  \
                             ta_zzzz_xxxyzz_1,  \
                             ta_zzzz_xxxzz_0,   \
                             ta_zzzz_xxxzz_1,   \
                             ta_zzzz_xxxzzz_0,  \
                             ta_zzzz_xxxzzz_1,  \
                             ta_zzzz_xxyyy_0,   \
                             ta_zzzz_xxyyy_1,   \
                             ta_zzzz_xxyyyy_0,  \
                             ta_zzzz_xxyyyy_1,  \
                             ta_zzzz_xxyyyz_0,  \
                             ta_zzzz_xxyyyz_1,  \
                             ta_zzzz_xxyyz_0,   \
                             ta_zzzz_xxyyz_1,   \
                             ta_zzzz_xxyyzz_0,  \
                             ta_zzzz_xxyyzz_1,  \
                             ta_zzzz_xxyzz_0,   \
                             ta_zzzz_xxyzz_1,   \
                             ta_zzzz_xxyzzz_0,  \
                             ta_zzzz_xxyzzz_1,  \
                             ta_zzzz_xxzzz_0,   \
                             ta_zzzz_xxzzz_1,   \
                             ta_zzzz_xxzzzz_0,  \
                             ta_zzzz_xxzzzz_1,  \
                             ta_zzzz_xyyyy_0,   \
                             ta_zzzz_xyyyy_1,   \
                             ta_zzzz_xyyyyy_0,  \
                             ta_zzzz_xyyyyy_1,  \
                             ta_zzzz_xyyyyz_0,  \
                             ta_zzzz_xyyyyz_1,  \
                             ta_zzzz_xyyyz_0,   \
                             ta_zzzz_xyyyz_1,   \
                             ta_zzzz_xyyyzz_0,  \
                             ta_zzzz_xyyyzz_1,  \
                             ta_zzzz_xyyzz_0,   \
                             ta_zzzz_xyyzz_1,   \
                             ta_zzzz_xyyzzz_0,  \
                             ta_zzzz_xyyzzz_1,  \
                             ta_zzzz_xyzzz_0,   \
                             ta_zzzz_xyzzz_1,   \
                             ta_zzzz_xyzzzz_0,  \
                             ta_zzzz_xyzzzz_1,  \
                             ta_zzzz_xzzzz_0,   \
                             ta_zzzz_xzzzz_1,   \
                             ta_zzzz_xzzzzz_0,  \
                             ta_zzzz_xzzzzz_1,  \
                             ta_zzzz_yyyyy_0,   \
                             ta_zzzz_yyyyy_1,   \
                             ta_zzzz_yyyyyy_0,  \
                             ta_zzzz_yyyyyy_1,  \
                             ta_zzzz_yyyyyz_0,  \
                             ta_zzzz_yyyyyz_1,  \
                             ta_zzzz_yyyyz_0,   \
                             ta_zzzz_yyyyz_1,   \
                             ta_zzzz_yyyyzz_0,  \
                             ta_zzzz_yyyyzz_1,  \
                             ta_zzzz_yyyzz_0,   \
                             ta_zzzz_yyyzz_1,   \
                             ta_zzzz_yyyzzz_0,  \
                             ta_zzzz_yyyzzz_1,  \
                             ta_zzzz_yyzzz_0,   \
                             ta_zzzz_yyzzz_1,   \
                             ta_zzzz_yyzzzz_0,  \
                             ta_zzzz_yyzzzz_1,  \
                             ta_zzzz_yzzzz_0,   \
                             ta_zzzz_yzzzz_1,   \
                             ta_zzzz_yzzzzz_0,  \
                             ta_zzzz_yzzzzz_1,  \
                             ta_zzzz_zzzzz_0,   \
                             ta_zzzz_zzzzz_1,   \
                             ta_zzzz_zzzzzz_0,  \
                             ta_zzzz_zzzzzz_1,  \
                             ta_zzzzz_xxxxxx_0, \
                             ta_zzzzz_xxxxxy_0, \
                             ta_zzzzz_xxxxxz_0, \
                             ta_zzzzz_xxxxyy_0, \
                             ta_zzzzz_xxxxyz_0, \
                             ta_zzzzz_xxxxzz_0, \
                             ta_zzzzz_xxxyyy_0, \
                             ta_zzzzz_xxxyyz_0, \
                             ta_zzzzz_xxxyzz_0, \
                             ta_zzzzz_xxxzzz_0, \
                             ta_zzzzz_xxyyyy_0, \
                             ta_zzzzz_xxyyyz_0, \
                             ta_zzzzz_xxyyzz_0, \
                             ta_zzzzz_xxyzzz_0, \
                             ta_zzzzz_xxzzzz_0, \
                             ta_zzzzz_xyyyyy_0, \
                             ta_zzzzz_xyyyyz_0, \
                             ta_zzzzz_xyyyzz_0, \
                             ta_zzzzz_xyyzzz_0, \
                             ta_zzzzz_xyzzzz_0, \
                             ta_zzzzz_xzzzzz_0, \
                             ta_zzzzz_yyyyyy_0, \
                             ta_zzzzz_yyyyyz_0, \
                             ta_zzzzz_yyyyzz_0, \
                             ta_zzzzz_yyyzzz_0, \
                             ta_zzzzz_yyzzzz_0, \
                             ta_zzzzz_yzzzzz_0, \
                             ta_zzzzz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzz_xxxxxx_0[i] =
            4.0 * ta_zzz_xxxxxx_0[i] * fe_0 - 4.0 * ta_zzz_xxxxxx_1[i] * fe_0 + ta_zzzz_xxxxxx_0[i] * pa_z[i] - ta_zzzz_xxxxxx_1[i] * pc_z[i];

        ta_zzzzz_xxxxxy_0[i] =
            4.0 * ta_zzz_xxxxxy_0[i] * fe_0 - 4.0 * ta_zzz_xxxxxy_1[i] * fe_0 + ta_zzzz_xxxxxy_0[i] * pa_z[i] - ta_zzzz_xxxxxy_1[i] * pc_z[i];

        ta_zzzzz_xxxxxz_0[i] = 4.0 * ta_zzz_xxxxxz_0[i] * fe_0 - 4.0 * ta_zzz_xxxxxz_1[i] * fe_0 + ta_zzzz_xxxxx_0[i] * fe_0 -
                               ta_zzzz_xxxxx_1[i] * fe_0 + ta_zzzz_xxxxxz_0[i] * pa_z[i] - ta_zzzz_xxxxxz_1[i] * pc_z[i];

        ta_zzzzz_xxxxyy_0[i] =
            4.0 * ta_zzz_xxxxyy_0[i] * fe_0 - 4.0 * ta_zzz_xxxxyy_1[i] * fe_0 + ta_zzzz_xxxxyy_0[i] * pa_z[i] - ta_zzzz_xxxxyy_1[i] * pc_z[i];

        ta_zzzzz_xxxxyz_0[i] = 4.0 * ta_zzz_xxxxyz_0[i] * fe_0 - 4.0 * ta_zzz_xxxxyz_1[i] * fe_0 + ta_zzzz_xxxxy_0[i] * fe_0 -
                               ta_zzzz_xxxxy_1[i] * fe_0 + ta_zzzz_xxxxyz_0[i] * pa_z[i] - ta_zzzz_xxxxyz_1[i] * pc_z[i];

        ta_zzzzz_xxxxzz_0[i] = 4.0 * ta_zzz_xxxxzz_0[i] * fe_0 - 4.0 * ta_zzz_xxxxzz_1[i] * fe_0 + 2.0 * ta_zzzz_xxxxz_0[i] * fe_0 -
                               2.0 * ta_zzzz_xxxxz_1[i] * fe_0 + ta_zzzz_xxxxzz_0[i] * pa_z[i] - ta_zzzz_xxxxzz_1[i] * pc_z[i];

        ta_zzzzz_xxxyyy_0[i] =
            4.0 * ta_zzz_xxxyyy_0[i] * fe_0 - 4.0 * ta_zzz_xxxyyy_1[i] * fe_0 + ta_zzzz_xxxyyy_0[i] * pa_z[i] - ta_zzzz_xxxyyy_1[i] * pc_z[i];

        ta_zzzzz_xxxyyz_0[i] = 4.0 * ta_zzz_xxxyyz_0[i] * fe_0 - 4.0 * ta_zzz_xxxyyz_1[i] * fe_0 + ta_zzzz_xxxyy_0[i] * fe_0 -
                               ta_zzzz_xxxyy_1[i] * fe_0 + ta_zzzz_xxxyyz_0[i] * pa_z[i] - ta_zzzz_xxxyyz_1[i] * pc_z[i];

        ta_zzzzz_xxxyzz_0[i] = 4.0 * ta_zzz_xxxyzz_0[i] * fe_0 - 4.0 * ta_zzz_xxxyzz_1[i] * fe_0 + 2.0 * ta_zzzz_xxxyz_0[i] * fe_0 -
                               2.0 * ta_zzzz_xxxyz_1[i] * fe_0 + ta_zzzz_xxxyzz_0[i] * pa_z[i] - ta_zzzz_xxxyzz_1[i] * pc_z[i];

        ta_zzzzz_xxxzzz_0[i] = 4.0 * ta_zzz_xxxzzz_0[i] * fe_0 - 4.0 * ta_zzz_xxxzzz_1[i] * fe_0 + 3.0 * ta_zzzz_xxxzz_0[i] * fe_0 -
                               3.0 * ta_zzzz_xxxzz_1[i] * fe_0 + ta_zzzz_xxxzzz_0[i] * pa_z[i] - ta_zzzz_xxxzzz_1[i] * pc_z[i];

        ta_zzzzz_xxyyyy_0[i] =
            4.0 * ta_zzz_xxyyyy_0[i] * fe_0 - 4.0 * ta_zzz_xxyyyy_1[i] * fe_0 + ta_zzzz_xxyyyy_0[i] * pa_z[i] - ta_zzzz_xxyyyy_1[i] * pc_z[i];

        ta_zzzzz_xxyyyz_0[i] = 4.0 * ta_zzz_xxyyyz_0[i] * fe_0 - 4.0 * ta_zzz_xxyyyz_1[i] * fe_0 + ta_zzzz_xxyyy_0[i] * fe_0 -
                               ta_zzzz_xxyyy_1[i] * fe_0 + ta_zzzz_xxyyyz_0[i] * pa_z[i] - ta_zzzz_xxyyyz_1[i] * pc_z[i];

        ta_zzzzz_xxyyzz_0[i] = 4.0 * ta_zzz_xxyyzz_0[i] * fe_0 - 4.0 * ta_zzz_xxyyzz_1[i] * fe_0 + 2.0 * ta_zzzz_xxyyz_0[i] * fe_0 -
                               2.0 * ta_zzzz_xxyyz_1[i] * fe_0 + ta_zzzz_xxyyzz_0[i] * pa_z[i] - ta_zzzz_xxyyzz_1[i] * pc_z[i];

        ta_zzzzz_xxyzzz_0[i] = 4.0 * ta_zzz_xxyzzz_0[i] * fe_0 - 4.0 * ta_zzz_xxyzzz_1[i] * fe_0 + 3.0 * ta_zzzz_xxyzz_0[i] * fe_0 -
                               3.0 * ta_zzzz_xxyzz_1[i] * fe_0 + ta_zzzz_xxyzzz_0[i] * pa_z[i] - ta_zzzz_xxyzzz_1[i] * pc_z[i];

        ta_zzzzz_xxzzzz_0[i] = 4.0 * ta_zzz_xxzzzz_0[i] * fe_0 - 4.0 * ta_zzz_xxzzzz_1[i] * fe_0 + 4.0 * ta_zzzz_xxzzz_0[i] * fe_0 -
                               4.0 * ta_zzzz_xxzzz_1[i] * fe_0 + ta_zzzz_xxzzzz_0[i] * pa_z[i] - ta_zzzz_xxzzzz_1[i] * pc_z[i];

        ta_zzzzz_xyyyyy_0[i] =
            4.0 * ta_zzz_xyyyyy_0[i] * fe_0 - 4.0 * ta_zzz_xyyyyy_1[i] * fe_0 + ta_zzzz_xyyyyy_0[i] * pa_z[i] - ta_zzzz_xyyyyy_1[i] * pc_z[i];

        ta_zzzzz_xyyyyz_0[i] = 4.0 * ta_zzz_xyyyyz_0[i] * fe_0 - 4.0 * ta_zzz_xyyyyz_1[i] * fe_0 + ta_zzzz_xyyyy_0[i] * fe_0 -
                               ta_zzzz_xyyyy_1[i] * fe_0 + ta_zzzz_xyyyyz_0[i] * pa_z[i] - ta_zzzz_xyyyyz_1[i] * pc_z[i];

        ta_zzzzz_xyyyzz_0[i] = 4.0 * ta_zzz_xyyyzz_0[i] * fe_0 - 4.0 * ta_zzz_xyyyzz_1[i] * fe_0 + 2.0 * ta_zzzz_xyyyz_0[i] * fe_0 -
                               2.0 * ta_zzzz_xyyyz_1[i] * fe_0 + ta_zzzz_xyyyzz_0[i] * pa_z[i] - ta_zzzz_xyyyzz_1[i] * pc_z[i];

        ta_zzzzz_xyyzzz_0[i] = 4.0 * ta_zzz_xyyzzz_0[i] * fe_0 - 4.0 * ta_zzz_xyyzzz_1[i] * fe_0 + 3.0 * ta_zzzz_xyyzz_0[i] * fe_0 -
                               3.0 * ta_zzzz_xyyzz_1[i] * fe_0 + ta_zzzz_xyyzzz_0[i] * pa_z[i] - ta_zzzz_xyyzzz_1[i] * pc_z[i];

        ta_zzzzz_xyzzzz_0[i] = 4.0 * ta_zzz_xyzzzz_0[i] * fe_0 - 4.0 * ta_zzz_xyzzzz_1[i] * fe_0 + 4.0 * ta_zzzz_xyzzz_0[i] * fe_0 -
                               4.0 * ta_zzzz_xyzzz_1[i] * fe_0 + ta_zzzz_xyzzzz_0[i] * pa_z[i] - ta_zzzz_xyzzzz_1[i] * pc_z[i];

        ta_zzzzz_xzzzzz_0[i] = 4.0 * ta_zzz_xzzzzz_0[i] * fe_0 - 4.0 * ta_zzz_xzzzzz_1[i] * fe_0 + 5.0 * ta_zzzz_xzzzz_0[i] * fe_0 -
                               5.0 * ta_zzzz_xzzzz_1[i] * fe_0 + ta_zzzz_xzzzzz_0[i] * pa_z[i] - ta_zzzz_xzzzzz_1[i] * pc_z[i];

        ta_zzzzz_yyyyyy_0[i] =
            4.0 * ta_zzz_yyyyyy_0[i] * fe_0 - 4.0 * ta_zzz_yyyyyy_1[i] * fe_0 + ta_zzzz_yyyyyy_0[i] * pa_z[i] - ta_zzzz_yyyyyy_1[i] * pc_z[i];

        ta_zzzzz_yyyyyz_0[i] = 4.0 * ta_zzz_yyyyyz_0[i] * fe_0 - 4.0 * ta_zzz_yyyyyz_1[i] * fe_0 + ta_zzzz_yyyyy_0[i] * fe_0 -
                               ta_zzzz_yyyyy_1[i] * fe_0 + ta_zzzz_yyyyyz_0[i] * pa_z[i] - ta_zzzz_yyyyyz_1[i] * pc_z[i];

        ta_zzzzz_yyyyzz_0[i] = 4.0 * ta_zzz_yyyyzz_0[i] * fe_0 - 4.0 * ta_zzz_yyyyzz_1[i] * fe_0 + 2.0 * ta_zzzz_yyyyz_0[i] * fe_0 -
                               2.0 * ta_zzzz_yyyyz_1[i] * fe_0 + ta_zzzz_yyyyzz_0[i] * pa_z[i] - ta_zzzz_yyyyzz_1[i] * pc_z[i];

        ta_zzzzz_yyyzzz_0[i] = 4.0 * ta_zzz_yyyzzz_0[i] * fe_0 - 4.0 * ta_zzz_yyyzzz_1[i] * fe_0 + 3.0 * ta_zzzz_yyyzz_0[i] * fe_0 -
                               3.0 * ta_zzzz_yyyzz_1[i] * fe_0 + ta_zzzz_yyyzzz_0[i] * pa_z[i] - ta_zzzz_yyyzzz_1[i] * pc_z[i];

        ta_zzzzz_yyzzzz_0[i] = 4.0 * ta_zzz_yyzzzz_0[i] * fe_0 - 4.0 * ta_zzz_yyzzzz_1[i] * fe_0 + 4.0 * ta_zzzz_yyzzz_0[i] * fe_0 -
                               4.0 * ta_zzzz_yyzzz_1[i] * fe_0 + ta_zzzz_yyzzzz_0[i] * pa_z[i] - ta_zzzz_yyzzzz_1[i] * pc_z[i];

        ta_zzzzz_yzzzzz_0[i] = 4.0 * ta_zzz_yzzzzz_0[i] * fe_0 - 4.0 * ta_zzz_yzzzzz_1[i] * fe_0 + 5.0 * ta_zzzz_yzzzz_0[i] * fe_0 -
                               5.0 * ta_zzzz_yzzzz_1[i] * fe_0 + ta_zzzz_yzzzzz_0[i] * pa_z[i] - ta_zzzz_yzzzzz_1[i] * pc_z[i];

        ta_zzzzz_zzzzzz_0[i] = 4.0 * ta_zzz_zzzzzz_0[i] * fe_0 - 4.0 * ta_zzz_zzzzzz_1[i] * fe_0 + 6.0 * ta_zzzz_zzzzz_0[i] * fe_0 -
                               6.0 * ta_zzzz_zzzzz_1[i] * fe_0 + ta_zzzz_zzzzzz_0[i] * pa_z[i] - ta_zzzz_zzzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
