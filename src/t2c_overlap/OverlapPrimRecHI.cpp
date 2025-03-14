#include "OverlapPrimRecHI.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_hi(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_hi,
                     const size_t idx_ovl_fi,
                     const size_t idx_ovl_gh,
                     const size_t idx_ovl_gi,
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

    auto ts_xxy_xxxxxz = pbuffer.data(idx_ovl_fi + 30);

    auto ts_xxy_xxxxzz = pbuffer.data(idx_ovl_fi + 33);

    auto ts_xxy_xxxzzz = pbuffer.data(idx_ovl_fi + 37);

    auto ts_xxy_xxzzzz = pbuffer.data(idx_ovl_fi + 42);

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

    auto ts_xxz_xxxxzz = pbuffer.data(idx_ovl_fi + 61);

    auto ts_xxz_xxxyyy = pbuffer.data(idx_ovl_fi + 62);

    auto ts_xxz_xxxzzz = pbuffer.data(idx_ovl_fi + 65);

    auto ts_xxz_xxyyyy = pbuffer.data(idx_ovl_fi + 66);

    auto ts_xxz_xxzzzz = pbuffer.data(idx_ovl_fi + 70);

    auto ts_xxz_xyyyyy = pbuffer.data(idx_ovl_fi + 71);

    auto ts_xxz_xzzzzz = pbuffer.data(idx_ovl_fi + 76);

    auto ts_xxz_yyyyyz = pbuffer.data(idx_ovl_fi + 78);

    auto ts_xxz_yyyyzz = pbuffer.data(idx_ovl_fi + 79);

    auto ts_xxz_yyyzzz = pbuffer.data(idx_ovl_fi + 80);

    auto ts_xxz_yyzzzz = pbuffer.data(idx_ovl_fi + 81);

    auto ts_xxz_yzzzzz = pbuffer.data(idx_ovl_fi + 82);

    auto ts_xxz_zzzzzz = pbuffer.data(idx_ovl_fi + 83);

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

    auto ts_yyz_xxxxzz = pbuffer.data(idx_ovl_fi + 201);

    auto ts_yyz_xxxyyy = pbuffer.data(idx_ovl_fi + 202);

    auto ts_yyz_xxxzzz = pbuffer.data(idx_ovl_fi + 205);

    auto ts_yyz_xxyyyy = pbuffer.data(idx_ovl_fi + 206);

    auto ts_yyz_xxzzzz = pbuffer.data(idx_ovl_fi + 210);

    auto ts_yyz_xyyyyy = pbuffer.data(idx_ovl_fi + 211);

    auto ts_yyz_xzzzzz = pbuffer.data(idx_ovl_fi + 216);

    auto ts_yyz_yyyyyy = pbuffer.data(idx_ovl_fi + 217);

    auto ts_yyz_yyyyyz = pbuffer.data(idx_ovl_fi + 218);

    auto ts_yyz_yyyyzz = pbuffer.data(idx_ovl_fi + 219);

    auto ts_yyz_yyyzzz = pbuffer.data(idx_ovl_fi + 220);

    auto ts_yyz_yyzzzz = pbuffer.data(idx_ovl_fi + 221);

    auto ts_yyz_yzzzzz = pbuffer.data(idx_ovl_fi + 222);

    auto ts_yyz_zzzzzz = pbuffer.data(idx_ovl_fi + 223);

    auto ts_yzz_xxxxxx = pbuffer.data(idx_ovl_fi + 224);

    auto ts_yzz_xxxxxz = pbuffer.data(idx_ovl_fi + 226);

    auto ts_yzz_xxxxyz = pbuffer.data(idx_ovl_fi + 228);

    auto ts_yzz_xxxxzz = pbuffer.data(idx_ovl_fi + 229);

    auto ts_yzz_xxxyyz = pbuffer.data(idx_ovl_fi + 231);

    auto ts_yzz_xxxyzz = pbuffer.data(idx_ovl_fi + 232);

    auto ts_yzz_xxxzzz = pbuffer.data(idx_ovl_fi + 233);

    auto ts_yzz_xxyyyz = pbuffer.data(idx_ovl_fi + 235);

    auto ts_yzz_xxyyzz = pbuffer.data(idx_ovl_fi + 236);

    auto ts_yzz_xxyzzz = pbuffer.data(idx_ovl_fi + 237);

    auto ts_yzz_xxzzzz = pbuffer.data(idx_ovl_fi + 238);

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

    // Set up components of auxiliary buffer : GH

    auto ts_xxxx_xxxxx = pbuffer.data(idx_ovl_gh);

    auto ts_xxxx_xxxxy = pbuffer.data(idx_ovl_gh + 1);

    auto ts_xxxx_xxxxz = pbuffer.data(idx_ovl_gh + 2);

    auto ts_xxxx_xxxyy = pbuffer.data(idx_ovl_gh + 3);

    auto ts_xxxx_xxxyz = pbuffer.data(idx_ovl_gh + 4);

    auto ts_xxxx_xxxzz = pbuffer.data(idx_ovl_gh + 5);

    auto ts_xxxx_xxyyy = pbuffer.data(idx_ovl_gh + 6);

    auto ts_xxxx_xxyyz = pbuffer.data(idx_ovl_gh + 7);

    auto ts_xxxx_xxyzz = pbuffer.data(idx_ovl_gh + 8);

    auto ts_xxxx_xxzzz = pbuffer.data(idx_ovl_gh + 9);

    auto ts_xxxx_xyyyy = pbuffer.data(idx_ovl_gh + 10);

    auto ts_xxxx_xyyyz = pbuffer.data(idx_ovl_gh + 11);

    auto ts_xxxx_xyyzz = pbuffer.data(idx_ovl_gh + 12);

    auto ts_xxxx_xyzzz = pbuffer.data(idx_ovl_gh + 13);

    auto ts_xxxx_xzzzz = pbuffer.data(idx_ovl_gh + 14);

    auto ts_xxxx_yyyyy = pbuffer.data(idx_ovl_gh + 15);

    auto ts_xxxx_yyyyz = pbuffer.data(idx_ovl_gh + 16);

    auto ts_xxxx_yyyzz = pbuffer.data(idx_ovl_gh + 17);

    auto ts_xxxx_yyzzz = pbuffer.data(idx_ovl_gh + 18);

    auto ts_xxxx_yzzzz = pbuffer.data(idx_ovl_gh + 19);

    auto ts_xxxx_zzzzz = pbuffer.data(idx_ovl_gh + 20);

    auto ts_xxxz_xxxxz = pbuffer.data(idx_ovl_gh + 44);

    auto ts_xxxz_xxxyz = pbuffer.data(idx_ovl_gh + 46);

    auto ts_xxxz_xxxzz = pbuffer.data(idx_ovl_gh + 47);

    auto ts_xxxz_xxyyz = pbuffer.data(idx_ovl_gh + 49);

    auto ts_xxxz_xxyzz = pbuffer.data(idx_ovl_gh + 50);

    auto ts_xxxz_xxzzz = pbuffer.data(idx_ovl_gh + 51);

    auto ts_xxxz_xyyyz = pbuffer.data(idx_ovl_gh + 53);

    auto ts_xxxz_xyyzz = pbuffer.data(idx_ovl_gh + 54);

    auto ts_xxxz_xyzzz = pbuffer.data(idx_ovl_gh + 55);

    auto ts_xxxz_xzzzz = pbuffer.data(idx_ovl_gh + 56);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_ovl_gh + 64);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_ovl_gh + 66);

    auto ts_xxyy_xxxyz = pbuffer.data(idx_ovl_gh + 67);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_ovl_gh + 69);

    auto ts_xxyy_xxyyz = pbuffer.data(idx_ovl_gh + 70);

    auto ts_xxyy_xxyzz = pbuffer.data(idx_ovl_gh + 71);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_ovl_gh + 73);

    auto ts_xxyy_xyyyz = pbuffer.data(idx_ovl_gh + 74);

    auto ts_xxyy_xyyzz = pbuffer.data(idx_ovl_gh + 75);

    auto ts_xxyy_xyzzz = pbuffer.data(idx_ovl_gh + 76);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_ovl_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_ovl_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_ovl_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_ovl_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_ovl_gh + 82);

    auto ts_xxzz_xxxxx = pbuffer.data(idx_ovl_gh + 105);

    auto ts_xxzz_xxxxy = pbuffer.data(idx_ovl_gh + 106);

    auto ts_xxzz_xxxxz = pbuffer.data(idx_ovl_gh + 107);

    auto ts_xxzz_xxxyy = pbuffer.data(idx_ovl_gh + 108);

    auto ts_xxzz_xxxyz = pbuffer.data(idx_ovl_gh + 109);

    auto ts_xxzz_xxxzz = pbuffer.data(idx_ovl_gh + 110);

    auto ts_xxzz_xxyyy = pbuffer.data(idx_ovl_gh + 111);

    auto ts_xxzz_xxyyz = pbuffer.data(idx_ovl_gh + 112);

    auto ts_xxzz_xxyzz = pbuffer.data(idx_ovl_gh + 113);

    auto ts_xxzz_xxzzz = pbuffer.data(idx_ovl_gh + 114);

    auto ts_xxzz_xyyyy = pbuffer.data(idx_ovl_gh + 115);

    auto ts_xxzz_xyyyz = pbuffer.data(idx_ovl_gh + 116);

    auto ts_xxzz_xyyzz = pbuffer.data(idx_ovl_gh + 117);

    auto ts_xxzz_xyzzz = pbuffer.data(idx_ovl_gh + 118);

    auto ts_xxzz_xzzzz = pbuffer.data(idx_ovl_gh + 119);

    auto ts_xxzz_yyyyz = pbuffer.data(idx_ovl_gh + 121);

    auto ts_xxzz_yyyzz = pbuffer.data(idx_ovl_gh + 122);

    auto ts_xxzz_yyzzz = pbuffer.data(idx_ovl_gh + 123);

    auto ts_xxzz_yzzzz = pbuffer.data(idx_ovl_gh + 124);

    auto ts_xxzz_zzzzz = pbuffer.data(idx_ovl_gh + 125);

    auto ts_xyyy_xxxxy = pbuffer.data(idx_ovl_gh + 127);

    auto ts_xyyy_xxxyy = pbuffer.data(idx_ovl_gh + 129);

    auto ts_xyyy_xxxyz = pbuffer.data(idx_ovl_gh + 130);

    auto ts_xyyy_xxyyy = pbuffer.data(idx_ovl_gh + 132);

    auto ts_xyyy_xxyyz = pbuffer.data(idx_ovl_gh + 133);

    auto ts_xyyy_xxyzz = pbuffer.data(idx_ovl_gh + 134);

    auto ts_xyyy_xyyyy = pbuffer.data(idx_ovl_gh + 136);

    auto ts_xyyy_xyyyz = pbuffer.data(idx_ovl_gh + 137);

    auto ts_xyyy_xyyzz = pbuffer.data(idx_ovl_gh + 138);

    auto ts_xyyy_xyzzz = pbuffer.data(idx_ovl_gh + 139);

    auto ts_xyyy_yyyyy = pbuffer.data(idx_ovl_gh + 141);

    auto ts_xyyy_yyyyz = pbuffer.data(idx_ovl_gh + 142);

    auto ts_xyyy_yyyzz = pbuffer.data(idx_ovl_gh + 143);

    auto ts_xyyy_yyzzz = pbuffer.data(idx_ovl_gh + 144);

    auto ts_xyyy_yzzzz = pbuffer.data(idx_ovl_gh + 145);

    auto ts_xzzz_xxxxz = pbuffer.data(idx_ovl_gh + 191);

    auto ts_xzzz_xxxyz = pbuffer.data(idx_ovl_gh + 193);

    auto ts_xzzz_xxxzz = pbuffer.data(idx_ovl_gh + 194);

    auto ts_xzzz_xxyyz = pbuffer.data(idx_ovl_gh + 196);

    auto ts_xzzz_xxyzz = pbuffer.data(idx_ovl_gh + 197);

    auto ts_xzzz_xxzzz = pbuffer.data(idx_ovl_gh + 198);

    auto ts_xzzz_xyyyz = pbuffer.data(idx_ovl_gh + 200);

    auto ts_xzzz_xyyzz = pbuffer.data(idx_ovl_gh + 201);

    auto ts_xzzz_xyzzz = pbuffer.data(idx_ovl_gh + 202);

    auto ts_xzzz_xzzzz = pbuffer.data(idx_ovl_gh + 203);

    auto ts_xzzz_yyyyz = pbuffer.data(idx_ovl_gh + 205);

    auto ts_xzzz_yyyzz = pbuffer.data(idx_ovl_gh + 206);

    auto ts_xzzz_yyzzz = pbuffer.data(idx_ovl_gh + 207);

    auto ts_xzzz_yzzzz = pbuffer.data(idx_ovl_gh + 208);

    auto ts_xzzz_zzzzz = pbuffer.data(idx_ovl_gh + 209);

    auto ts_yyyy_xxxxx = pbuffer.data(idx_ovl_gh + 210);

    auto ts_yyyy_xxxxy = pbuffer.data(idx_ovl_gh + 211);

    auto ts_yyyy_xxxxz = pbuffer.data(idx_ovl_gh + 212);

    auto ts_yyyy_xxxyy = pbuffer.data(idx_ovl_gh + 213);

    auto ts_yyyy_xxxyz = pbuffer.data(idx_ovl_gh + 214);

    auto ts_yyyy_xxxzz = pbuffer.data(idx_ovl_gh + 215);

    auto ts_yyyy_xxyyy = pbuffer.data(idx_ovl_gh + 216);

    auto ts_yyyy_xxyyz = pbuffer.data(idx_ovl_gh + 217);

    auto ts_yyyy_xxyzz = pbuffer.data(idx_ovl_gh + 218);

    auto ts_yyyy_xxzzz = pbuffer.data(idx_ovl_gh + 219);

    auto ts_yyyy_xyyyy = pbuffer.data(idx_ovl_gh + 220);

    auto ts_yyyy_xyyyz = pbuffer.data(idx_ovl_gh + 221);

    auto ts_yyyy_xyyzz = pbuffer.data(idx_ovl_gh + 222);

    auto ts_yyyy_xyzzz = pbuffer.data(idx_ovl_gh + 223);

    auto ts_yyyy_xzzzz = pbuffer.data(idx_ovl_gh + 224);

    auto ts_yyyy_yyyyy = pbuffer.data(idx_ovl_gh + 225);

    auto ts_yyyy_yyyyz = pbuffer.data(idx_ovl_gh + 226);

    auto ts_yyyy_yyyzz = pbuffer.data(idx_ovl_gh + 227);

    auto ts_yyyy_yyzzz = pbuffer.data(idx_ovl_gh + 228);

    auto ts_yyyy_yzzzz = pbuffer.data(idx_ovl_gh + 229);

    auto ts_yyyy_zzzzz = pbuffer.data(idx_ovl_gh + 230);

    auto ts_yyyz_xxxxz = pbuffer.data(idx_ovl_gh + 233);

    auto ts_yyyz_xxxyz = pbuffer.data(idx_ovl_gh + 235);

    auto ts_yyyz_xxxzz = pbuffer.data(idx_ovl_gh + 236);

    auto ts_yyyz_xxyyz = pbuffer.data(idx_ovl_gh + 238);

    auto ts_yyyz_xxyzz = pbuffer.data(idx_ovl_gh + 239);

    auto ts_yyyz_xxzzz = pbuffer.data(idx_ovl_gh + 240);

    auto ts_yyyz_xyyyz = pbuffer.data(idx_ovl_gh + 242);

    auto ts_yyyz_xyyzz = pbuffer.data(idx_ovl_gh + 243);

    auto ts_yyyz_xyzzz = pbuffer.data(idx_ovl_gh + 244);

    auto ts_yyyz_xzzzz = pbuffer.data(idx_ovl_gh + 245);

    auto ts_yyyz_yyyyz = pbuffer.data(idx_ovl_gh + 247);

    auto ts_yyyz_yyyzz = pbuffer.data(idx_ovl_gh + 248);

    auto ts_yyyz_yyzzz = pbuffer.data(idx_ovl_gh + 249);

    auto ts_yyyz_yzzzz = pbuffer.data(idx_ovl_gh + 250);

    auto ts_yyyz_zzzzz = pbuffer.data(idx_ovl_gh + 251);

    auto ts_yyzz_xxxxx = pbuffer.data(idx_ovl_gh + 252);

    auto ts_yyzz_xxxxy = pbuffer.data(idx_ovl_gh + 253);

    auto ts_yyzz_xxxxz = pbuffer.data(idx_ovl_gh + 254);

    auto ts_yyzz_xxxyy = pbuffer.data(idx_ovl_gh + 255);

    auto ts_yyzz_xxxyz = pbuffer.data(idx_ovl_gh + 256);

    auto ts_yyzz_xxxzz = pbuffer.data(idx_ovl_gh + 257);

    auto ts_yyzz_xxyyy = pbuffer.data(idx_ovl_gh + 258);

    auto ts_yyzz_xxyyz = pbuffer.data(idx_ovl_gh + 259);

    auto ts_yyzz_xxyzz = pbuffer.data(idx_ovl_gh + 260);

    auto ts_yyzz_xxzzz = pbuffer.data(idx_ovl_gh + 261);

    auto ts_yyzz_xyyyy = pbuffer.data(idx_ovl_gh + 262);

    auto ts_yyzz_xyyyz = pbuffer.data(idx_ovl_gh + 263);

    auto ts_yyzz_xyyzz = pbuffer.data(idx_ovl_gh + 264);

    auto ts_yyzz_xyzzz = pbuffer.data(idx_ovl_gh + 265);

    auto ts_yyzz_xzzzz = pbuffer.data(idx_ovl_gh + 266);

    auto ts_yyzz_yyyyy = pbuffer.data(idx_ovl_gh + 267);

    auto ts_yyzz_yyyyz = pbuffer.data(idx_ovl_gh + 268);

    auto ts_yyzz_yyyzz = pbuffer.data(idx_ovl_gh + 269);

    auto ts_yyzz_yyzzz = pbuffer.data(idx_ovl_gh + 270);

    auto ts_yyzz_yzzzz = pbuffer.data(idx_ovl_gh + 271);

    auto ts_yyzz_zzzzz = pbuffer.data(idx_ovl_gh + 272);

    auto ts_yzzz_xxxxy = pbuffer.data(idx_ovl_gh + 274);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_ovl_gh + 275);

    auto ts_yzzz_xxxyy = pbuffer.data(idx_ovl_gh + 276);

    auto ts_yzzz_xxxyz = pbuffer.data(idx_ovl_gh + 277);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_ovl_gh + 278);

    auto ts_yzzz_xxyyy = pbuffer.data(idx_ovl_gh + 279);

    auto ts_yzzz_xxyyz = pbuffer.data(idx_ovl_gh + 280);

    auto ts_yzzz_xxyzz = pbuffer.data(idx_ovl_gh + 281);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_ovl_gh + 282);

    auto ts_yzzz_xyyyy = pbuffer.data(idx_ovl_gh + 283);

    auto ts_yzzz_xyyyz = pbuffer.data(idx_ovl_gh + 284);

    auto ts_yzzz_xyyzz = pbuffer.data(idx_ovl_gh + 285);

    auto ts_yzzz_xyzzz = pbuffer.data(idx_ovl_gh + 286);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_ovl_gh + 287);

    auto ts_yzzz_yyyyy = pbuffer.data(idx_ovl_gh + 288);

    auto ts_yzzz_yyyyz = pbuffer.data(idx_ovl_gh + 289);

    auto ts_yzzz_yyyzz = pbuffer.data(idx_ovl_gh + 290);

    auto ts_yzzz_yyzzz = pbuffer.data(idx_ovl_gh + 291);

    auto ts_yzzz_yzzzz = pbuffer.data(idx_ovl_gh + 292);

    auto ts_yzzz_zzzzz = pbuffer.data(idx_ovl_gh + 293);

    auto ts_zzzz_xxxxx = pbuffer.data(idx_ovl_gh + 294);

    auto ts_zzzz_xxxxy = pbuffer.data(idx_ovl_gh + 295);

    auto ts_zzzz_xxxxz = pbuffer.data(idx_ovl_gh + 296);

    auto ts_zzzz_xxxyy = pbuffer.data(idx_ovl_gh + 297);

    auto ts_zzzz_xxxyz = pbuffer.data(idx_ovl_gh + 298);

    auto ts_zzzz_xxxzz = pbuffer.data(idx_ovl_gh + 299);

    auto ts_zzzz_xxyyy = pbuffer.data(idx_ovl_gh + 300);

    auto ts_zzzz_xxyyz = pbuffer.data(idx_ovl_gh + 301);

    auto ts_zzzz_xxyzz = pbuffer.data(idx_ovl_gh + 302);

    auto ts_zzzz_xxzzz = pbuffer.data(idx_ovl_gh + 303);

    auto ts_zzzz_xyyyy = pbuffer.data(idx_ovl_gh + 304);

    auto ts_zzzz_xyyyz = pbuffer.data(idx_ovl_gh + 305);

    auto ts_zzzz_xyyzz = pbuffer.data(idx_ovl_gh + 306);

    auto ts_zzzz_xyzzz = pbuffer.data(idx_ovl_gh + 307);

    auto ts_zzzz_xzzzz = pbuffer.data(idx_ovl_gh + 308);

    auto ts_zzzz_yyyyy = pbuffer.data(idx_ovl_gh + 309);

    auto ts_zzzz_yyyyz = pbuffer.data(idx_ovl_gh + 310);

    auto ts_zzzz_yyyzz = pbuffer.data(idx_ovl_gh + 311);

    auto ts_zzzz_yyzzz = pbuffer.data(idx_ovl_gh + 312);

    auto ts_zzzz_yzzzz = pbuffer.data(idx_ovl_gh + 313);

    auto ts_zzzz_zzzzz = pbuffer.data(idx_ovl_gh + 314);

    // Set up components of auxiliary buffer : GI

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

    auto ts_xxxy_xxxxxx = pbuffer.data(idx_ovl_gi + 28);

    auto ts_xxxy_xxxxxy = pbuffer.data(idx_ovl_gi + 29);

    auto ts_xxxy_xxxxxz = pbuffer.data(idx_ovl_gi + 30);

    auto ts_xxxy_xxxxyy = pbuffer.data(idx_ovl_gi + 31);

    auto ts_xxxy_xxxxzz = pbuffer.data(idx_ovl_gi + 33);

    auto ts_xxxy_xxxyyy = pbuffer.data(idx_ovl_gi + 34);

    auto ts_xxxy_xxxzzz = pbuffer.data(idx_ovl_gi + 37);

    auto ts_xxxy_xxyyyy = pbuffer.data(idx_ovl_gi + 38);

    auto ts_xxxy_xxzzzz = pbuffer.data(idx_ovl_gi + 42);

    auto ts_xxxy_xyyyyy = pbuffer.data(idx_ovl_gi + 43);

    auto ts_xxxy_xzzzzz = pbuffer.data(idx_ovl_gi + 48);

    auto ts_xxxy_yyyyyy = pbuffer.data(idx_ovl_gi + 49);

    auto ts_xxxy_yyyyyz = pbuffer.data(idx_ovl_gi + 50);

    auto ts_xxxy_yyyyzz = pbuffer.data(idx_ovl_gi + 51);

    auto ts_xxxy_yyyzzz = pbuffer.data(idx_ovl_gi + 52);

    auto ts_xxxy_yyzzzz = pbuffer.data(idx_ovl_gi + 53);

    auto ts_xxxy_yzzzzz = pbuffer.data(idx_ovl_gi + 54);

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

    auto ts_xxxz_yyyyyz = pbuffer.data(idx_ovl_gi + 78);

    auto ts_xxxz_yyyyzz = pbuffer.data(idx_ovl_gi + 79);

    auto ts_xxxz_yyyzzz = pbuffer.data(idx_ovl_gi + 80);

    auto ts_xxxz_yyzzzz = pbuffer.data(idx_ovl_gi + 81);

    auto ts_xxxz_yzzzzz = pbuffer.data(idx_ovl_gi + 82);

    auto ts_xxxz_zzzzzz = pbuffer.data(idx_ovl_gi + 83);

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

    auto ts_xxyz_xxxxxz = pbuffer.data(idx_ovl_gi + 114);

    auto ts_xxyz_xxxxzz = pbuffer.data(idx_ovl_gi + 117);

    auto ts_xxyz_xxxzzz = pbuffer.data(idx_ovl_gi + 121);

    auto ts_xxyz_xxzzzz = pbuffer.data(idx_ovl_gi + 126);

    auto ts_xxyz_xzzzzz = pbuffer.data(idx_ovl_gi + 132);

    auto ts_xxyz_yyyyyz = pbuffer.data(idx_ovl_gi + 134);

    auto ts_xxyz_yyyyzz = pbuffer.data(idx_ovl_gi + 135);

    auto ts_xxyz_yyyzzz = pbuffer.data(idx_ovl_gi + 136);

    auto ts_xxyz_yyzzzz = pbuffer.data(idx_ovl_gi + 137);

    auto ts_xxyz_yzzzzz = pbuffer.data(idx_ovl_gi + 138);

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

    auto ts_xyyy_xxxxxx = pbuffer.data(idx_ovl_gi + 168);

    auto ts_xyyy_xxxxxy = pbuffer.data(idx_ovl_gi + 169);

    auto ts_xyyy_xxxxyy = pbuffer.data(idx_ovl_gi + 171);

    auto ts_xyyy_xxxxyz = pbuffer.data(idx_ovl_gi + 172);

    auto ts_xyyy_xxxyyy = pbuffer.data(idx_ovl_gi + 174);

    auto ts_xyyy_xxxyyz = pbuffer.data(idx_ovl_gi + 175);

    auto ts_xyyy_xxxyzz = pbuffer.data(idx_ovl_gi + 176);

    auto ts_xyyy_xxyyyy = pbuffer.data(idx_ovl_gi + 178);

    auto ts_xyyy_xxyyyz = pbuffer.data(idx_ovl_gi + 179);

    auto ts_xyyy_xxyyzz = pbuffer.data(idx_ovl_gi + 180);

    auto ts_xyyy_xxyzzz = pbuffer.data(idx_ovl_gi + 181);

    auto ts_xyyy_xyyyyy = pbuffer.data(idx_ovl_gi + 183);

    auto ts_xyyy_xyyyyz = pbuffer.data(idx_ovl_gi + 184);

    auto ts_xyyy_xyyyzz = pbuffer.data(idx_ovl_gi + 185);

    auto ts_xyyy_xyyzzz = pbuffer.data(idx_ovl_gi + 186);

    auto ts_xyyy_xyzzzz = pbuffer.data(idx_ovl_gi + 187);

    auto ts_xyyy_yyyyyy = pbuffer.data(idx_ovl_gi + 189);

    auto ts_xyyy_yyyyyz = pbuffer.data(idx_ovl_gi + 190);

    auto ts_xyyy_yyyyzz = pbuffer.data(idx_ovl_gi + 191);

    auto ts_xyyy_yyyzzz = pbuffer.data(idx_ovl_gi + 192);

    auto ts_xyyy_yyzzzz = pbuffer.data(idx_ovl_gi + 193);

    auto ts_xyyy_yzzzzz = pbuffer.data(idx_ovl_gi + 194);

    auto ts_xyyy_zzzzzz = pbuffer.data(idx_ovl_gi + 195);

    auto ts_xyyz_yyyyyz = pbuffer.data(idx_ovl_gi + 218);

    auto ts_xyyz_yyyyzz = pbuffer.data(idx_ovl_gi + 219);

    auto ts_xyyz_yyyzzz = pbuffer.data(idx_ovl_gi + 220);

    auto ts_xyyz_yyzzzz = pbuffer.data(idx_ovl_gi + 221);

    auto ts_xyyz_yzzzzz = pbuffer.data(idx_ovl_gi + 222);

    auto ts_xyyz_zzzzzz = pbuffer.data(idx_ovl_gi + 223);

    auto ts_xyzz_yyyyyy = pbuffer.data(idx_ovl_gi + 245);

    auto ts_xyzz_yyyyyz = pbuffer.data(idx_ovl_gi + 246);

    auto ts_xyzz_yyyyzz = pbuffer.data(idx_ovl_gi + 247);

    auto ts_xyzz_yyyzzz = pbuffer.data(idx_ovl_gi + 248);

    auto ts_xyzz_yyzzzz = pbuffer.data(idx_ovl_gi + 249);

    auto ts_xyzz_yzzzzz = pbuffer.data(idx_ovl_gi + 250);

    auto ts_xzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 252);

    auto ts_xzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 254);

    auto ts_xzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 256);

    auto ts_xzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 257);

    auto ts_xzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 259);

    auto ts_xzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 260);

    auto ts_xzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 261);

    auto ts_xzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 263);

    auto ts_xzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 264);

    auto ts_xzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 265);

    auto ts_xzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 266);

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

    // Set up 0-28 components of targeted buffer : HI

    auto ts_xxxxx_xxxxxx = pbuffer.data(idx_ovl_hi);

    auto ts_xxxxx_xxxxxy = pbuffer.data(idx_ovl_hi + 1);

    auto ts_xxxxx_xxxxxz = pbuffer.data(idx_ovl_hi + 2);

    auto ts_xxxxx_xxxxyy = pbuffer.data(idx_ovl_hi + 3);

    auto ts_xxxxx_xxxxyz = pbuffer.data(idx_ovl_hi + 4);

    auto ts_xxxxx_xxxxzz = pbuffer.data(idx_ovl_hi + 5);

    auto ts_xxxxx_xxxyyy = pbuffer.data(idx_ovl_hi + 6);

    auto ts_xxxxx_xxxyyz = pbuffer.data(idx_ovl_hi + 7);

    auto ts_xxxxx_xxxyzz = pbuffer.data(idx_ovl_hi + 8);

    auto ts_xxxxx_xxxzzz = pbuffer.data(idx_ovl_hi + 9);

    auto ts_xxxxx_xxyyyy = pbuffer.data(idx_ovl_hi + 10);

    auto ts_xxxxx_xxyyyz = pbuffer.data(idx_ovl_hi + 11);

    auto ts_xxxxx_xxyyzz = pbuffer.data(idx_ovl_hi + 12);

    auto ts_xxxxx_xxyzzz = pbuffer.data(idx_ovl_hi + 13);

    auto ts_xxxxx_xxzzzz = pbuffer.data(idx_ovl_hi + 14);

    auto ts_xxxxx_xyyyyy = pbuffer.data(idx_ovl_hi + 15);

    auto ts_xxxxx_xyyyyz = pbuffer.data(idx_ovl_hi + 16);

    auto ts_xxxxx_xyyyzz = pbuffer.data(idx_ovl_hi + 17);

    auto ts_xxxxx_xyyzzz = pbuffer.data(idx_ovl_hi + 18);

    auto ts_xxxxx_xyzzzz = pbuffer.data(idx_ovl_hi + 19);

    auto ts_xxxxx_xzzzzz = pbuffer.data(idx_ovl_hi + 20);

    auto ts_xxxxx_yyyyyy = pbuffer.data(idx_ovl_hi + 21);

    auto ts_xxxxx_yyyyyz = pbuffer.data(idx_ovl_hi + 22);

    auto ts_xxxxx_yyyyzz = pbuffer.data(idx_ovl_hi + 23);

    auto ts_xxxxx_yyyzzz = pbuffer.data(idx_ovl_hi + 24);

    auto ts_xxxxx_yyzzzz = pbuffer.data(idx_ovl_hi + 25);

    auto ts_xxxxx_yzzzzz = pbuffer.data(idx_ovl_hi + 26);

    auto ts_xxxxx_zzzzzz = pbuffer.data(idx_ovl_hi + 27);

    #pragma omp simd aligned(pa_x, ts_xxx_xxxxxx, ts_xxx_xxxxxy, ts_xxx_xxxxxz, ts_xxx_xxxxyy, ts_xxx_xxxxyz, ts_xxx_xxxxzz, ts_xxx_xxxyyy, ts_xxx_xxxyyz, ts_xxx_xxxyzz, ts_xxx_xxxzzz, ts_xxx_xxyyyy, ts_xxx_xxyyyz, ts_xxx_xxyyzz, ts_xxx_xxyzzz, ts_xxx_xxzzzz, ts_xxx_xyyyyy, ts_xxx_xyyyyz, ts_xxx_xyyyzz, ts_xxx_xyyzzz, ts_xxx_xyzzzz, ts_xxx_xzzzzz, ts_xxx_yyyyyy, ts_xxx_yyyyyz, ts_xxx_yyyyzz, ts_xxx_yyyzzz, ts_xxx_yyzzzz, ts_xxx_yzzzzz, ts_xxx_zzzzzz, ts_xxxx_xxxxx, ts_xxxx_xxxxxx, ts_xxxx_xxxxxy, ts_xxxx_xxxxxz, ts_xxxx_xxxxy, ts_xxxx_xxxxyy, ts_xxxx_xxxxyz, ts_xxxx_xxxxz, ts_xxxx_xxxxzz, ts_xxxx_xxxyy, ts_xxxx_xxxyyy, ts_xxxx_xxxyyz, ts_xxxx_xxxyz, ts_xxxx_xxxyzz, ts_xxxx_xxxzz, ts_xxxx_xxxzzz, ts_xxxx_xxyyy, ts_xxxx_xxyyyy, ts_xxxx_xxyyyz, ts_xxxx_xxyyz, ts_xxxx_xxyyzz, ts_xxxx_xxyzz, ts_xxxx_xxyzzz, ts_xxxx_xxzzz, ts_xxxx_xxzzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyyy, ts_xxxx_xyyyyz, ts_xxxx_xyyyz, ts_xxxx_xyyyzz, ts_xxxx_xyyzz, ts_xxxx_xyyzzz, ts_xxxx_xyzzz, ts_xxxx_xyzzzz, ts_xxxx_xzzzz, ts_xxxx_xzzzzz, ts_xxxx_yyyyy, ts_xxxx_yyyyyy, ts_xxxx_yyyyyz, ts_xxxx_yyyyz, ts_xxxx_yyyyzz, ts_xxxx_yyyzz, ts_xxxx_yyyzzz, ts_xxxx_yyzzz, ts_xxxx_yyzzzz, ts_xxxx_yzzzz, ts_xxxx_yzzzzz, ts_xxxx_zzzzz, ts_xxxx_zzzzzz, ts_xxxxx_xxxxxx, ts_xxxxx_xxxxxy, ts_xxxxx_xxxxxz, ts_xxxxx_xxxxyy, ts_xxxxx_xxxxyz, ts_xxxxx_xxxxzz, ts_xxxxx_xxxyyy, ts_xxxxx_xxxyyz, ts_xxxxx_xxxyzz, ts_xxxxx_xxxzzz, ts_xxxxx_xxyyyy, ts_xxxxx_xxyyyz, ts_xxxxx_xxyyzz, ts_xxxxx_xxyzzz, ts_xxxxx_xxzzzz, ts_xxxxx_xyyyyy, ts_xxxxx_xyyyyz, ts_xxxxx_xyyyzz, ts_xxxxx_xyyzzz, ts_xxxxx_xyzzzz, ts_xxxxx_xzzzzz, ts_xxxxx_yyyyyy, ts_xxxxx_yyyyyz, ts_xxxxx_yyyyzz, ts_xxxxx_yyyzzz, ts_xxxxx_yyzzzz, ts_xxxxx_yzzzzz, ts_xxxxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxx_xxxxxx[i] = 4.0 * ts_xxx_xxxxxx[i] * fe_0 + 6.0 * ts_xxxx_xxxxx[i] * fe_0 + ts_xxxx_xxxxxx[i] * pa_x[i];

        ts_xxxxx_xxxxxy[i] = 4.0 * ts_xxx_xxxxxy[i] * fe_0 + 5.0 * ts_xxxx_xxxxy[i] * fe_0 + ts_xxxx_xxxxxy[i] * pa_x[i];

        ts_xxxxx_xxxxxz[i] = 4.0 * ts_xxx_xxxxxz[i] * fe_0 + 5.0 * ts_xxxx_xxxxz[i] * fe_0 + ts_xxxx_xxxxxz[i] * pa_x[i];

        ts_xxxxx_xxxxyy[i] = 4.0 * ts_xxx_xxxxyy[i] * fe_0 + 4.0 * ts_xxxx_xxxyy[i] * fe_0 + ts_xxxx_xxxxyy[i] * pa_x[i];

        ts_xxxxx_xxxxyz[i] = 4.0 * ts_xxx_xxxxyz[i] * fe_0 + 4.0 * ts_xxxx_xxxyz[i] * fe_0 + ts_xxxx_xxxxyz[i] * pa_x[i];

        ts_xxxxx_xxxxzz[i] = 4.0 * ts_xxx_xxxxzz[i] * fe_0 + 4.0 * ts_xxxx_xxxzz[i] * fe_0 + ts_xxxx_xxxxzz[i] * pa_x[i];

        ts_xxxxx_xxxyyy[i] = 4.0 * ts_xxx_xxxyyy[i] * fe_0 + 3.0 * ts_xxxx_xxyyy[i] * fe_0 + ts_xxxx_xxxyyy[i] * pa_x[i];

        ts_xxxxx_xxxyyz[i] = 4.0 * ts_xxx_xxxyyz[i] * fe_0 + 3.0 * ts_xxxx_xxyyz[i] * fe_0 + ts_xxxx_xxxyyz[i] * pa_x[i];

        ts_xxxxx_xxxyzz[i] = 4.0 * ts_xxx_xxxyzz[i] * fe_0 + 3.0 * ts_xxxx_xxyzz[i] * fe_0 + ts_xxxx_xxxyzz[i] * pa_x[i];

        ts_xxxxx_xxxzzz[i] = 4.0 * ts_xxx_xxxzzz[i] * fe_0 + 3.0 * ts_xxxx_xxzzz[i] * fe_0 + ts_xxxx_xxxzzz[i] * pa_x[i];

        ts_xxxxx_xxyyyy[i] = 4.0 * ts_xxx_xxyyyy[i] * fe_0 + 2.0 * ts_xxxx_xyyyy[i] * fe_0 + ts_xxxx_xxyyyy[i] * pa_x[i];

        ts_xxxxx_xxyyyz[i] = 4.0 * ts_xxx_xxyyyz[i] * fe_0 + 2.0 * ts_xxxx_xyyyz[i] * fe_0 + ts_xxxx_xxyyyz[i] * pa_x[i];

        ts_xxxxx_xxyyzz[i] = 4.0 * ts_xxx_xxyyzz[i] * fe_0 + 2.0 * ts_xxxx_xyyzz[i] * fe_0 + ts_xxxx_xxyyzz[i] * pa_x[i];

        ts_xxxxx_xxyzzz[i] = 4.0 * ts_xxx_xxyzzz[i] * fe_0 + 2.0 * ts_xxxx_xyzzz[i] * fe_0 + ts_xxxx_xxyzzz[i] * pa_x[i];

        ts_xxxxx_xxzzzz[i] = 4.0 * ts_xxx_xxzzzz[i] * fe_0 + 2.0 * ts_xxxx_xzzzz[i] * fe_0 + ts_xxxx_xxzzzz[i] * pa_x[i];

        ts_xxxxx_xyyyyy[i] = 4.0 * ts_xxx_xyyyyy[i] * fe_0 + ts_xxxx_yyyyy[i] * fe_0 + ts_xxxx_xyyyyy[i] * pa_x[i];

        ts_xxxxx_xyyyyz[i] = 4.0 * ts_xxx_xyyyyz[i] * fe_0 + ts_xxxx_yyyyz[i] * fe_0 + ts_xxxx_xyyyyz[i] * pa_x[i];

        ts_xxxxx_xyyyzz[i] = 4.0 * ts_xxx_xyyyzz[i] * fe_0 + ts_xxxx_yyyzz[i] * fe_0 + ts_xxxx_xyyyzz[i] * pa_x[i];

        ts_xxxxx_xyyzzz[i] = 4.0 * ts_xxx_xyyzzz[i] * fe_0 + ts_xxxx_yyzzz[i] * fe_0 + ts_xxxx_xyyzzz[i] * pa_x[i];

        ts_xxxxx_xyzzzz[i] = 4.0 * ts_xxx_xyzzzz[i] * fe_0 + ts_xxxx_yzzzz[i] * fe_0 + ts_xxxx_xyzzzz[i] * pa_x[i];

        ts_xxxxx_xzzzzz[i] = 4.0 * ts_xxx_xzzzzz[i] * fe_0 + ts_xxxx_zzzzz[i] * fe_0 + ts_xxxx_xzzzzz[i] * pa_x[i];

        ts_xxxxx_yyyyyy[i] = 4.0 * ts_xxx_yyyyyy[i] * fe_0 + ts_xxxx_yyyyyy[i] * pa_x[i];

        ts_xxxxx_yyyyyz[i] = 4.0 * ts_xxx_yyyyyz[i] * fe_0 + ts_xxxx_yyyyyz[i] * pa_x[i];

        ts_xxxxx_yyyyzz[i] = 4.0 * ts_xxx_yyyyzz[i] * fe_0 + ts_xxxx_yyyyzz[i] * pa_x[i];

        ts_xxxxx_yyyzzz[i] = 4.0 * ts_xxx_yyyzzz[i] * fe_0 + ts_xxxx_yyyzzz[i] * pa_x[i];

        ts_xxxxx_yyzzzz[i] = 4.0 * ts_xxx_yyzzzz[i] * fe_0 + ts_xxxx_yyzzzz[i] * pa_x[i];

        ts_xxxxx_yzzzzz[i] = 4.0 * ts_xxx_yzzzzz[i] * fe_0 + ts_xxxx_yzzzzz[i] * pa_x[i];

        ts_xxxxx_zzzzzz[i] = 4.0 * ts_xxx_zzzzzz[i] * fe_0 + ts_xxxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : HI

    auto ts_xxxxy_xxxxxx = pbuffer.data(idx_ovl_hi + 28);

    auto ts_xxxxy_xxxxxy = pbuffer.data(idx_ovl_hi + 29);

    auto ts_xxxxy_xxxxxz = pbuffer.data(idx_ovl_hi + 30);

    auto ts_xxxxy_xxxxyy = pbuffer.data(idx_ovl_hi + 31);

    auto ts_xxxxy_xxxxyz = pbuffer.data(idx_ovl_hi + 32);

    auto ts_xxxxy_xxxxzz = pbuffer.data(idx_ovl_hi + 33);

    auto ts_xxxxy_xxxyyy = pbuffer.data(idx_ovl_hi + 34);

    auto ts_xxxxy_xxxyyz = pbuffer.data(idx_ovl_hi + 35);

    auto ts_xxxxy_xxxyzz = pbuffer.data(idx_ovl_hi + 36);

    auto ts_xxxxy_xxxzzz = pbuffer.data(idx_ovl_hi + 37);

    auto ts_xxxxy_xxyyyy = pbuffer.data(idx_ovl_hi + 38);

    auto ts_xxxxy_xxyyyz = pbuffer.data(idx_ovl_hi + 39);

    auto ts_xxxxy_xxyyzz = pbuffer.data(idx_ovl_hi + 40);

    auto ts_xxxxy_xxyzzz = pbuffer.data(idx_ovl_hi + 41);

    auto ts_xxxxy_xxzzzz = pbuffer.data(idx_ovl_hi + 42);

    auto ts_xxxxy_xyyyyy = pbuffer.data(idx_ovl_hi + 43);

    auto ts_xxxxy_xyyyyz = pbuffer.data(idx_ovl_hi + 44);

    auto ts_xxxxy_xyyyzz = pbuffer.data(idx_ovl_hi + 45);

    auto ts_xxxxy_xyyzzz = pbuffer.data(idx_ovl_hi + 46);

    auto ts_xxxxy_xyzzzz = pbuffer.data(idx_ovl_hi + 47);

    auto ts_xxxxy_xzzzzz = pbuffer.data(idx_ovl_hi + 48);

    auto ts_xxxxy_yyyyyy = pbuffer.data(idx_ovl_hi + 49);

    auto ts_xxxxy_yyyyyz = pbuffer.data(idx_ovl_hi + 50);

    auto ts_xxxxy_yyyyzz = pbuffer.data(idx_ovl_hi + 51);

    auto ts_xxxxy_yyyzzz = pbuffer.data(idx_ovl_hi + 52);

    auto ts_xxxxy_yyzzzz = pbuffer.data(idx_ovl_hi + 53);

    auto ts_xxxxy_yzzzzz = pbuffer.data(idx_ovl_hi + 54);

    auto ts_xxxxy_zzzzzz = pbuffer.data(idx_ovl_hi + 55);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_xxxxx, ts_xxxx_xxxxxx, ts_xxxx_xxxxxy, ts_xxxx_xxxxxz, ts_xxxx_xxxxy, ts_xxxx_xxxxyy, ts_xxxx_xxxxyz, ts_xxxx_xxxxz, ts_xxxx_xxxxzz, ts_xxxx_xxxyy, ts_xxxx_xxxyyy, ts_xxxx_xxxyyz, ts_xxxx_xxxyz, ts_xxxx_xxxyzz, ts_xxxx_xxxzz, ts_xxxx_xxxzzz, ts_xxxx_xxyyy, ts_xxxx_xxyyyy, ts_xxxx_xxyyyz, ts_xxxx_xxyyz, ts_xxxx_xxyyzz, ts_xxxx_xxyzz, ts_xxxx_xxyzzz, ts_xxxx_xxzzz, ts_xxxx_xxzzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyyy, ts_xxxx_xyyyyz, ts_xxxx_xyyyz, ts_xxxx_xyyyzz, ts_xxxx_xyyzz, ts_xxxx_xyyzzz, ts_xxxx_xyzzz, ts_xxxx_xyzzzz, ts_xxxx_xzzzz, ts_xxxx_xzzzzz, ts_xxxx_zzzzzz, ts_xxxxy_xxxxxx, ts_xxxxy_xxxxxy, ts_xxxxy_xxxxxz, ts_xxxxy_xxxxyy, ts_xxxxy_xxxxyz, ts_xxxxy_xxxxzz, ts_xxxxy_xxxyyy, ts_xxxxy_xxxyyz, ts_xxxxy_xxxyzz, ts_xxxxy_xxxzzz, ts_xxxxy_xxyyyy, ts_xxxxy_xxyyyz, ts_xxxxy_xxyyzz, ts_xxxxy_xxyzzz, ts_xxxxy_xxzzzz, ts_xxxxy_xyyyyy, ts_xxxxy_xyyyyz, ts_xxxxy_xyyyzz, ts_xxxxy_xyyzzz, ts_xxxxy_xyzzzz, ts_xxxxy_xzzzzz, ts_xxxxy_yyyyyy, ts_xxxxy_yyyyyz, ts_xxxxy_yyyyzz, ts_xxxxy_yyyzzz, ts_xxxxy_yyzzzz, ts_xxxxy_yzzzzz, ts_xxxxy_zzzzzz, ts_xxxy_yyyyyy, ts_xxxy_yyyyyz, ts_xxxy_yyyyzz, ts_xxxy_yyyzzz, ts_xxxy_yyzzzz, ts_xxxy_yzzzzz, ts_xxy_yyyyyy, ts_xxy_yyyyyz, ts_xxy_yyyyzz, ts_xxy_yyyzzz, ts_xxy_yyzzzz, ts_xxy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxy_xxxxxx[i] = ts_xxxx_xxxxxx[i] * pa_y[i];

        ts_xxxxy_xxxxxy[i] = ts_xxxx_xxxxx[i] * fe_0 + ts_xxxx_xxxxxy[i] * pa_y[i];

        ts_xxxxy_xxxxxz[i] = ts_xxxx_xxxxxz[i] * pa_y[i];

        ts_xxxxy_xxxxyy[i] = 2.0 * ts_xxxx_xxxxy[i] * fe_0 + ts_xxxx_xxxxyy[i] * pa_y[i];

        ts_xxxxy_xxxxyz[i] = ts_xxxx_xxxxz[i] * fe_0 + ts_xxxx_xxxxyz[i] * pa_y[i];

        ts_xxxxy_xxxxzz[i] = ts_xxxx_xxxxzz[i] * pa_y[i];

        ts_xxxxy_xxxyyy[i] = 3.0 * ts_xxxx_xxxyy[i] * fe_0 + ts_xxxx_xxxyyy[i] * pa_y[i];

        ts_xxxxy_xxxyyz[i] = 2.0 * ts_xxxx_xxxyz[i] * fe_0 + ts_xxxx_xxxyyz[i] * pa_y[i];

        ts_xxxxy_xxxyzz[i] = ts_xxxx_xxxzz[i] * fe_0 + ts_xxxx_xxxyzz[i] * pa_y[i];

        ts_xxxxy_xxxzzz[i] = ts_xxxx_xxxzzz[i] * pa_y[i];

        ts_xxxxy_xxyyyy[i] = 4.0 * ts_xxxx_xxyyy[i] * fe_0 + ts_xxxx_xxyyyy[i] * pa_y[i];

        ts_xxxxy_xxyyyz[i] = 3.0 * ts_xxxx_xxyyz[i] * fe_0 + ts_xxxx_xxyyyz[i] * pa_y[i];

        ts_xxxxy_xxyyzz[i] = 2.0 * ts_xxxx_xxyzz[i] * fe_0 + ts_xxxx_xxyyzz[i] * pa_y[i];

        ts_xxxxy_xxyzzz[i] = ts_xxxx_xxzzz[i] * fe_0 + ts_xxxx_xxyzzz[i] * pa_y[i];

        ts_xxxxy_xxzzzz[i] = ts_xxxx_xxzzzz[i] * pa_y[i];

        ts_xxxxy_xyyyyy[i] = 5.0 * ts_xxxx_xyyyy[i] * fe_0 + ts_xxxx_xyyyyy[i] * pa_y[i];

        ts_xxxxy_xyyyyz[i] = 4.0 * ts_xxxx_xyyyz[i] * fe_0 + ts_xxxx_xyyyyz[i] * pa_y[i];

        ts_xxxxy_xyyyzz[i] = 3.0 * ts_xxxx_xyyzz[i] * fe_0 + ts_xxxx_xyyyzz[i] * pa_y[i];

        ts_xxxxy_xyyzzz[i] = 2.0 * ts_xxxx_xyzzz[i] * fe_0 + ts_xxxx_xyyzzz[i] * pa_y[i];

        ts_xxxxy_xyzzzz[i] = ts_xxxx_xzzzz[i] * fe_0 + ts_xxxx_xyzzzz[i] * pa_y[i];

        ts_xxxxy_xzzzzz[i] = ts_xxxx_xzzzzz[i] * pa_y[i];

        ts_xxxxy_yyyyyy[i] = 3.0 * ts_xxy_yyyyyy[i] * fe_0 + ts_xxxy_yyyyyy[i] * pa_x[i];

        ts_xxxxy_yyyyyz[i] = 3.0 * ts_xxy_yyyyyz[i] * fe_0 + ts_xxxy_yyyyyz[i] * pa_x[i];

        ts_xxxxy_yyyyzz[i] = 3.0 * ts_xxy_yyyyzz[i] * fe_0 + ts_xxxy_yyyyzz[i] * pa_x[i];

        ts_xxxxy_yyyzzz[i] = 3.0 * ts_xxy_yyyzzz[i] * fe_0 + ts_xxxy_yyyzzz[i] * pa_x[i];

        ts_xxxxy_yyzzzz[i] = 3.0 * ts_xxy_yyzzzz[i] * fe_0 + ts_xxxy_yyzzzz[i] * pa_x[i];

        ts_xxxxy_yzzzzz[i] = 3.0 * ts_xxy_yzzzzz[i] * fe_0 + ts_xxxy_yzzzzz[i] * pa_x[i];

        ts_xxxxy_zzzzzz[i] = ts_xxxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : HI

    auto ts_xxxxz_xxxxxx = pbuffer.data(idx_ovl_hi + 56);

    auto ts_xxxxz_xxxxxy = pbuffer.data(idx_ovl_hi + 57);

    auto ts_xxxxz_xxxxxz = pbuffer.data(idx_ovl_hi + 58);

    auto ts_xxxxz_xxxxyy = pbuffer.data(idx_ovl_hi + 59);

    auto ts_xxxxz_xxxxyz = pbuffer.data(idx_ovl_hi + 60);

    auto ts_xxxxz_xxxxzz = pbuffer.data(idx_ovl_hi + 61);

    auto ts_xxxxz_xxxyyy = pbuffer.data(idx_ovl_hi + 62);

    auto ts_xxxxz_xxxyyz = pbuffer.data(idx_ovl_hi + 63);

    auto ts_xxxxz_xxxyzz = pbuffer.data(idx_ovl_hi + 64);

    auto ts_xxxxz_xxxzzz = pbuffer.data(idx_ovl_hi + 65);

    auto ts_xxxxz_xxyyyy = pbuffer.data(idx_ovl_hi + 66);

    auto ts_xxxxz_xxyyyz = pbuffer.data(idx_ovl_hi + 67);

    auto ts_xxxxz_xxyyzz = pbuffer.data(idx_ovl_hi + 68);

    auto ts_xxxxz_xxyzzz = pbuffer.data(idx_ovl_hi + 69);

    auto ts_xxxxz_xxzzzz = pbuffer.data(idx_ovl_hi + 70);

    auto ts_xxxxz_xyyyyy = pbuffer.data(idx_ovl_hi + 71);

    auto ts_xxxxz_xyyyyz = pbuffer.data(idx_ovl_hi + 72);

    auto ts_xxxxz_xyyyzz = pbuffer.data(idx_ovl_hi + 73);

    auto ts_xxxxz_xyyzzz = pbuffer.data(idx_ovl_hi + 74);

    auto ts_xxxxz_xyzzzz = pbuffer.data(idx_ovl_hi + 75);

    auto ts_xxxxz_xzzzzz = pbuffer.data(idx_ovl_hi + 76);

    auto ts_xxxxz_yyyyyy = pbuffer.data(idx_ovl_hi + 77);

    auto ts_xxxxz_yyyyyz = pbuffer.data(idx_ovl_hi + 78);

    auto ts_xxxxz_yyyyzz = pbuffer.data(idx_ovl_hi + 79);

    auto ts_xxxxz_yyyzzz = pbuffer.data(idx_ovl_hi + 80);

    auto ts_xxxxz_yyzzzz = pbuffer.data(idx_ovl_hi + 81);

    auto ts_xxxxz_yzzzzz = pbuffer.data(idx_ovl_hi + 82);

    auto ts_xxxxz_zzzzzz = pbuffer.data(idx_ovl_hi + 83);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_xxxxx, ts_xxxx_xxxxxx, ts_xxxx_xxxxxy, ts_xxxx_xxxxxz, ts_xxxx_xxxxy, ts_xxxx_xxxxyy, ts_xxxx_xxxxyz, ts_xxxx_xxxxz, ts_xxxx_xxxxzz, ts_xxxx_xxxyy, ts_xxxx_xxxyyy, ts_xxxx_xxxyyz, ts_xxxx_xxxyz, ts_xxxx_xxxyzz, ts_xxxx_xxxzz, ts_xxxx_xxxzzz, ts_xxxx_xxyyy, ts_xxxx_xxyyyy, ts_xxxx_xxyyyz, ts_xxxx_xxyyz, ts_xxxx_xxyyzz, ts_xxxx_xxyzz, ts_xxxx_xxyzzz, ts_xxxx_xxzzz, ts_xxxx_xxzzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyyy, ts_xxxx_xyyyyz, ts_xxxx_xyyyz, ts_xxxx_xyyyzz, ts_xxxx_xyyzz, ts_xxxx_xyyzzz, ts_xxxx_xyzzz, ts_xxxx_xyzzzz, ts_xxxx_xzzzz, ts_xxxx_xzzzzz, ts_xxxx_yyyyyy, ts_xxxxz_xxxxxx, ts_xxxxz_xxxxxy, ts_xxxxz_xxxxxz, ts_xxxxz_xxxxyy, ts_xxxxz_xxxxyz, ts_xxxxz_xxxxzz, ts_xxxxz_xxxyyy, ts_xxxxz_xxxyyz, ts_xxxxz_xxxyzz, ts_xxxxz_xxxzzz, ts_xxxxz_xxyyyy, ts_xxxxz_xxyyyz, ts_xxxxz_xxyyzz, ts_xxxxz_xxyzzz, ts_xxxxz_xxzzzz, ts_xxxxz_xyyyyy, ts_xxxxz_xyyyyz, ts_xxxxz_xyyyzz, ts_xxxxz_xyyzzz, ts_xxxxz_xyzzzz, ts_xxxxz_xzzzzz, ts_xxxxz_yyyyyy, ts_xxxxz_yyyyyz, ts_xxxxz_yyyyzz, ts_xxxxz_yyyzzz, ts_xxxxz_yyzzzz, ts_xxxxz_yzzzzz, ts_xxxxz_zzzzzz, ts_xxxz_yyyyyz, ts_xxxz_yyyyzz, ts_xxxz_yyyzzz, ts_xxxz_yyzzzz, ts_xxxz_yzzzzz, ts_xxxz_zzzzzz, ts_xxz_yyyyyz, ts_xxz_yyyyzz, ts_xxz_yyyzzz, ts_xxz_yyzzzz, ts_xxz_yzzzzz, ts_xxz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxz_xxxxxx[i] = ts_xxxx_xxxxxx[i] * pa_z[i];

        ts_xxxxz_xxxxxy[i] = ts_xxxx_xxxxxy[i] * pa_z[i];

        ts_xxxxz_xxxxxz[i] = ts_xxxx_xxxxx[i] * fe_0 + ts_xxxx_xxxxxz[i] * pa_z[i];

        ts_xxxxz_xxxxyy[i] = ts_xxxx_xxxxyy[i] * pa_z[i];

        ts_xxxxz_xxxxyz[i] = ts_xxxx_xxxxy[i] * fe_0 + ts_xxxx_xxxxyz[i] * pa_z[i];

        ts_xxxxz_xxxxzz[i] = 2.0 * ts_xxxx_xxxxz[i] * fe_0 + ts_xxxx_xxxxzz[i] * pa_z[i];

        ts_xxxxz_xxxyyy[i] = ts_xxxx_xxxyyy[i] * pa_z[i];

        ts_xxxxz_xxxyyz[i] = ts_xxxx_xxxyy[i] * fe_0 + ts_xxxx_xxxyyz[i] * pa_z[i];

        ts_xxxxz_xxxyzz[i] = 2.0 * ts_xxxx_xxxyz[i] * fe_0 + ts_xxxx_xxxyzz[i] * pa_z[i];

        ts_xxxxz_xxxzzz[i] = 3.0 * ts_xxxx_xxxzz[i] * fe_0 + ts_xxxx_xxxzzz[i] * pa_z[i];

        ts_xxxxz_xxyyyy[i] = ts_xxxx_xxyyyy[i] * pa_z[i];

        ts_xxxxz_xxyyyz[i] = ts_xxxx_xxyyy[i] * fe_0 + ts_xxxx_xxyyyz[i] * pa_z[i];

        ts_xxxxz_xxyyzz[i] = 2.0 * ts_xxxx_xxyyz[i] * fe_0 + ts_xxxx_xxyyzz[i] * pa_z[i];

        ts_xxxxz_xxyzzz[i] = 3.0 * ts_xxxx_xxyzz[i] * fe_0 + ts_xxxx_xxyzzz[i] * pa_z[i];

        ts_xxxxz_xxzzzz[i] = 4.0 * ts_xxxx_xxzzz[i] * fe_0 + ts_xxxx_xxzzzz[i] * pa_z[i];

        ts_xxxxz_xyyyyy[i] = ts_xxxx_xyyyyy[i] * pa_z[i];

        ts_xxxxz_xyyyyz[i] = ts_xxxx_xyyyy[i] * fe_0 + ts_xxxx_xyyyyz[i] * pa_z[i];

        ts_xxxxz_xyyyzz[i] = 2.0 * ts_xxxx_xyyyz[i] * fe_0 + ts_xxxx_xyyyzz[i] * pa_z[i];

        ts_xxxxz_xyyzzz[i] = 3.0 * ts_xxxx_xyyzz[i] * fe_0 + ts_xxxx_xyyzzz[i] * pa_z[i];

        ts_xxxxz_xyzzzz[i] = 4.0 * ts_xxxx_xyzzz[i] * fe_0 + ts_xxxx_xyzzzz[i] * pa_z[i];

        ts_xxxxz_xzzzzz[i] = 5.0 * ts_xxxx_xzzzz[i] * fe_0 + ts_xxxx_xzzzzz[i] * pa_z[i];

        ts_xxxxz_yyyyyy[i] = ts_xxxx_yyyyyy[i] * pa_z[i];

        ts_xxxxz_yyyyyz[i] = 3.0 * ts_xxz_yyyyyz[i] * fe_0 + ts_xxxz_yyyyyz[i] * pa_x[i];

        ts_xxxxz_yyyyzz[i] = 3.0 * ts_xxz_yyyyzz[i] * fe_0 + ts_xxxz_yyyyzz[i] * pa_x[i];

        ts_xxxxz_yyyzzz[i] = 3.0 * ts_xxz_yyyzzz[i] * fe_0 + ts_xxxz_yyyzzz[i] * pa_x[i];

        ts_xxxxz_yyzzzz[i] = 3.0 * ts_xxz_yyzzzz[i] * fe_0 + ts_xxxz_yyzzzz[i] * pa_x[i];

        ts_xxxxz_yzzzzz[i] = 3.0 * ts_xxz_yzzzzz[i] * fe_0 + ts_xxxz_yzzzzz[i] * pa_x[i];

        ts_xxxxz_zzzzzz[i] = 3.0 * ts_xxz_zzzzzz[i] * fe_0 + ts_xxxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 84-112 components of targeted buffer : HI

    auto ts_xxxyy_xxxxxx = pbuffer.data(idx_ovl_hi + 84);

    auto ts_xxxyy_xxxxxy = pbuffer.data(idx_ovl_hi + 85);

    auto ts_xxxyy_xxxxxz = pbuffer.data(idx_ovl_hi + 86);

    auto ts_xxxyy_xxxxyy = pbuffer.data(idx_ovl_hi + 87);

    auto ts_xxxyy_xxxxyz = pbuffer.data(idx_ovl_hi + 88);

    auto ts_xxxyy_xxxxzz = pbuffer.data(idx_ovl_hi + 89);

    auto ts_xxxyy_xxxyyy = pbuffer.data(idx_ovl_hi + 90);

    auto ts_xxxyy_xxxyyz = pbuffer.data(idx_ovl_hi + 91);

    auto ts_xxxyy_xxxyzz = pbuffer.data(idx_ovl_hi + 92);

    auto ts_xxxyy_xxxzzz = pbuffer.data(idx_ovl_hi + 93);

    auto ts_xxxyy_xxyyyy = pbuffer.data(idx_ovl_hi + 94);

    auto ts_xxxyy_xxyyyz = pbuffer.data(idx_ovl_hi + 95);

    auto ts_xxxyy_xxyyzz = pbuffer.data(idx_ovl_hi + 96);

    auto ts_xxxyy_xxyzzz = pbuffer.data(idx_ovl_hi + 97);

    auto ts_xxxyy_xxzzzz = pbuffer.data(idx_ovl_hi + 98);

    auto ts_xxxyy_xyyyyy = pbuffer.data(idx_ovl_hi + 99);

    auto ts_xxxyy_xyyyyz = pbuffer.data(idx_ovl_hi + 100);

    auto ts_xxxyy_xyyyzz = pbuffer.data(idx_ovl_hi + 101);

    auto ts_xxxyy_xyyzzz = pbuffer.data(idx_ovl_hi + 102);

    auto ts_xxxyy_xyzzzz = pbuffer.data(idx_ovl_hi + 103);

    auto ts_xxxyy_xzzzzz = pbuffer.data(idx_ovl_hi + 104);

    auto ts_xxxyy_yyyyyy = pbuffer.data(idx_ovl_hi + 105);

    auto ts_xxxyy_yyyyyz = pbuffer.data(idx_ovl_hi + 106);

    auto ts_xxxyy_yyyyzz = pbuffer.data(idx_ovl_hi + 107);

    auto ts_xxxyy_yyyzzz = pbuffer.data(idx_ovl_hi + 108);

    auto ts_xxxyy_yyzzzz = pbuffer.data(idx_ovl_hi + 109);

    auto ts_xxxyy_yzzzzz = pbuffer.data(idx_ovl_hi + 110);

    auto ts_xxxyy_zzzzzz = pbuffer.data(idx_ovl_hi + 111);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxx_xxxxxx, ts_xxx_xxxxxz, ts_xxx_xxxxzz, ts_xxx_xxxzzz, ts_xxx_xxzzzz, ts_xxx_xzzzzz, ts_xxxy_xxxxxx, ts_xxxy_xxxxxz, ts_xxxy_xxxxzz, ts_xxxy_xxxzzz, ts_xxxy_xxzzzz, ts_xxxy_xzzzzz, ts_xxxyy_xxxxxx, ts_xxxyy_xxxxxy, ts_xxxyy_xxxxxz, ts_xxxyy_xxxxyy, ts_xxxyy_xxxxyz, ts_xxxyy_xxxxzz, ts_xxxyy_xxxyyy, ts_xxxyy_xxxyyz, ts_xxxyy_xxxyzz, ts_xxxyy_xxxzzz, ts_xxxyy_xxyyyy, ts_xxxyy_xxyyyz, ts_xxxyy_xxyyzz, ts_xxxyy_xxyzzz, ts_xxxyy_xxzzzz, ts_xxxyy_xyyyyy, ts_xxxyy_xyyyyz, ts_xxxyy_xyyyzz, ts_xxxyy_xyyzzz, ts_xxxyy_xyzzzz, ts_xxxyy_xzzzzz, ts_xxxyy_yyyyyy, ts_xxxyy_yyyyyz, ts_xxxyy_yyyyzz, ts_xxxyy_yyyzzz, ts_xxxyy_yyzzzz, ts_xxxyy_yzzzzz, ts_xxxyy_zzzzzz, ts_xxyy_xxxxxy, ts_xxyy_xxxxy, ts_xxyy_xxxxyy, ts_xxyy_xxxxyz, ts_xxyy_xxxyy, ts_xxyy_xxxyyy, ts_xxyy_xxxyyz, ts_xxyy_xxxyz, ts_xxyy_xxxyzz, ts_xxyy_xxyyy, ts_xxyy_xxyyyy, ts_xxyy_xxyyyz, ts_xxyy_xxyyz, ts_xxyy_xxyyzz, ts_xxyy_xxyzz, ts_xxyy_xxyzzz, ts_xxyy_xyyyy, ts_xxyy_xyyyyy, ts_xxyy_xyyyyz, ts_xxyy_xyyyz, ts_xxyy_xyyyzz, ts_xxyy_xyyzz, ts_xxyy_xyyzzz, ts_xxyy_xyzzz, ts_xxyy_xyzzzz, ts_xxyy_yyyyy, ts_xxyy_yyyyyy, ts_xxyy_yyyyyz, ts_xxyy_yyyyz, ts_xxyy_yyyyzz, ts_xxyy_yyyzz, ts_xxyy_yyyzzz, ts_xxyy_yyzzz, ts_xxyy_yyzzzz, ts_xxyy_yzzzz, ts_xxyy_yzzzzz, ts_xxyy_zzzzzz, ts_xyy_xxxxxy, ts_xyy_xxxxyy, ts_xyy_xxxxyz, ts_xyy_xxxyyy, ts_xyy_xxxyyz, ts_xyy_xxxyzz, ts_xyy_xxyyyy, ts_xyy_xxyyyz, ts_xyy_xxyyzz, ts_xyy_xxyzzz, ts_xyy_xyyyyy, ts_xyy_xyyyyz, ts_xyy_xyyyzz, ts_xyy_xyyzzz, ts_xyy_xyzzzz, ts_xyy_yyyyyy, ts_xyy_yyyyyz, ts_xyy_yyyyzz, ts_xyy_yyyzzz, ts_xyy_yyzzzz, ts_xyy_yzzzzz, ts_xyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyy_xxxxxx[i] = ts_xxx_xxxxxx[i] * fe_0 + ts_xxxy_xxxxxx[i] * pa_y[i];

        ts_xxxyy_xxxxxy[i] = 2.0 * ts_xyy_xxxxxy[i] * fe_0 + 5.0 * ts_xxyy_xxxxy[i] * fe_0 + ts_xxyy_xxxxxy[i] * pa_x[i];

        ts_xxxyy_xxxxxz[i] = ts_xxx_xxxxxz[i] * fe_0 + ts_xxxy_xxxxxz[i] * pa_y[i];

        ts_xxxyy_xxxxyy[i] = 2.0 * ts_xyy_xxxxyy[i] * fe_0 + 4.0 * ts_xxyy_xxxyy[i] * fe_0 + ts_xxyy_xxxxyy[i] * pa_x[i];

        ts_xxxyy_xxxxyz[i] = 2.0 * ts_xyy_xxxxyz[i] * fe_0 + 4.0 * ts_xxyy_xxxyz[i] * fe_0 + ts_xxyy_xxxxyz[i] * pa_x[i];

        ts_xxxyy_xxxxzz[i] = ts_xxx_xxxxzz[i] * fe_0 + ts_xxxy_xxxxzz[i] * pa_y[i];

        ts_xxxyy_xxxyyy[i] = 2.0 * ts_xyy_xxxyyy[i] * fe_0 + 3.0 * ts_xxyy_xxyyy[i] * fe_0 + ts_xxyy_xxxyyy[i] * pa_x[i];

        ts_xxxyy_xxxyyz[i] = 2.0 * ts_xyy_xxxyyz[i] * fe_0 + 3.0 * ts_xxyy_xxyyz[i] * fe_0 + ts_xxyy_xxxyyz[i] * pa_x[i];

        ts_xxxyy_xxxyzz[i] = 2.0 * ts_xyy_xxxyzz[i] * fe_0 + 3.0 * ts_xxyy_xxyzz[i] * fe_0 + ts_xxyy_xxxyzz[i] * pa_x[i];

        ts_xxxyy_xxxzzz[i] = ts_xxx_xxxzzz[i] * fe_0 + ts_xxxy_xxxzzz[i] * pa_y[i];

        ts_xxxyy_xxyyyy[i] = 2.0 * ts_xyy_xxyyyy[i] * fe_0 + 2.0 * ts_xxyy_xyyyy[i] * fe_0 + ts_xxyy_xxyyyy[i] * pa_x[i];

        ts_xxxyy_xxyyyz[i] = 2.0 * ts_xyy_xxyyyz[i] * fe_0 + 2.0 * ts_xxyy_xyyyz[i] * fe_0 + ts_xxyy_xxyyyz[i] * pa_x[i];

        ts_xxxyy_xxyyzz[i] = 2.0 * ts_xyy_xxyyzz[i] * fe_0 + 2.0 * ts_xxyy_xyyzz[i] * fe_0 + ts_xxyy_xxyyzz[i] * pa_x[i];

        ts_xxxyy_xxyzzz[i] = 2.0 * ts_xyy_xxyzzz[i] * fe_0 + 2.0 * ts_xxyy_xyzzz[i] * fe_0 + ts_xxyy_xxyzzz[i] * pa_x[i];

        ts_xxxyy_xxzzzz[i] = ts_xxx_xxzzzz[i] * fe_0 + ts_xxxy_xxzzzz[i] * pa_y[i];

        ts_xxxyy_xyyyyy[i] = 2.0 * ts_xyy_xyyyyy[i] * fe_0 + ts_xxyy_yyyyy[i] * fe_0 + ts_xxyy_xyyyyy[i] * pa_x[i];

        ts_xxxyy_xyyyyz[i] = 2.0 * ts_xyy_xyyyyz[i] * fe_0 + ts_xxyy_yyyyz[i] * fe_0 + ts_xxyy_xyyyyz[i] * pa_x[i];

        ts_xxxyy_xyyyzz[i] = 2.0 * ts_xyy_xyyyzz[i] * fe_0 + ts_xxyy_yyyzz[i] * fe_0 + ts_xxyy_xyyyzz[i] * pa_x[i];

        ts_xxxyy_xyyzzz[i] = 2.0 * ts_xyy_xyyzzz[i] * fe_0 + ts_xxyy_yyzzz[i] * fe_0 + ts_xxyy_xyyzzz[i] * pa_x[i];

        ts_xxxyy_xyzzzz[i] = 2.0 * ts_xyy_xyzzzz[i] * fe_0 + ts_xxyy_yzzzz[i] * fe_0 + ts_xxyy_xyzzzz[i] * pa_x[i];

        ts_xxxyy_xzzzzz[i] = ts_xxx_xzzzzz[i] * fe_0 + ts_xxxy_xzzzzz[i] * pa_y[i];

        ts_xxxyy_yyyyyy[i] = 2.0 * ts_xyy_yyyyyy[i] * fe_0 + ts_xxyy_yyyyyy[i] * pa_x[i];

        ts_xxxyy_yyyyyz[i] = 2.0 * ts_xyy_yyyyyz[i] * fe_0 + ts_xxyy_yyyyyz[i] * pa_x[i];

        ts_xxxyy_yyyyzz[i] = 2.0 * ts_xyy_yyyyzz[i] * fe_0 + ts_xxyy_yyyyzz[i] * pa_x[i];

        ts_xxxyy_yyyzzz[i] = 2.0 * ts_xyy_yyyzzz[i] * fe_0 + ts_xxyy_yyyzzz[i] * pa_x[i];

        ts_xxxyy_yyzzzz[i] = 2.0 * ts_xyy_yyzzzz[i] * fe_0 + ts_xxyy_yyzzzz[i] * pa_x[i];

        ts_xxxyy_yzzzzz[i] = 2.0 * ts_xyy_yzzzzz[i] * fe_0 + ts_xxyy_yzzzzz[i] * pa_x[i];

        ts_xxxyy_zzzzzz[i] = 2.0 * ts_xyy_zzzzzz[i] * fe_0 + ts_xxyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : HI

    auto ts_xxxyz_xxxxxx = pbuffer.data(idx_ovl_hi + 112);

    auto ts_xxxyz_xxxxxy = pbuffer.data(idx_ovl_hi + 113);

    auto ts_xxxyz_xxxxxz = pbuffer.data(idx_ovl_hi + 114);

    auto ts_xxxyz_xxxxyy = pbuffer.data(idx_ovl_hi + 115);

    auto ts_xxxyz_xxxxyz = pbuffer.data(idx_ovl_hi + 116);

    auto ts_xxxyz_xxxxzz = pbuffer.data(idx_ovl_hi + 117);

    auto ts_xxxyz_xxxyyy = pbuffer.data(idx_ovl_hi + 118);

    auto ts_xxxyz_xxxyyz = pbuffer.data(idx_ovl_hi + 119);

    auto ts_xxxyz_xxxyzz = pbuffer.data(idx_ovl_hi + 120);

    auto ts_xxxyz_xxxzzz = pbuffer.data(idx_ovl_hi + 121);

    auto ts_xxxyz_xxyyyy = pbuffer.data(idx_ovl_hi + 122);

    auto ts_xxxyz_xxyyyz = pbuffer.data(idx_ovl_hi + 123);

    auto ts_xxxyz_xxyyzz = pbuffer.data(idx_ovl_hi + 124);

    auto ts_xxxyz_xxyzzz = pbuffer.data(idx_ovl_hi + 125);

    auto ts_xxxyz_xxzzzz = pbuffer.data(idx_ovl_hi + 126);

    auto ts_xxxyz_xyyyyy = pbuffer.data(idx_ovl_hi + 127);

    auto ts_xxxyz_xyyyyz = pbuffer.data(idx_ovl_hi + 128);

    auto ts_xxxyz_xyyyzz = pbuffer.data(idx_ovl_hi + 129);

    auto ts_xxxyz_xyyzzz = pbuffer.data(idx_ovl_hi + 130);

    auto ts_xxxyz_xyzzzz = pbuffer.data(idx_ovl_hi + 131);

    auto ts_xxxyz_xzzzzz = pbuffer.data(idx_ovl_hi + 132);

    auto ts_xxxyz_yyyyyy = pbuffer.data(idx_ovl_hi + 133);

    auto ts_xxxyz_yyyyyz = pbuffer.data(idx_ovl_hi + 134);

    auto ts_xxxyz_yyyyzz = pbuffer.data(idx_ovl_hi + 135);

    auto ts_xxxyz_yyyzzz = pbuffer.data(idx_ovl_hi + 136);

    auto ts_xxxyz_yyzzzz = pbuffer.data(idx_ovl_hi + 137);

    auto ts_xxxyz_yzzzzz = pbuffer.data(idx_ovl_hi + 138);

    auto ts_xxxyz_zzzzzz = pbuffer.data(idx_ovl_hi + 139);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxy_xxxxxy, ts_xxxy_xxxxyy, ts_xxxy_xxxyyy, ts_xxxy_xxyyyy, ts_xxxy_xyyyyy, ts_xxxy_yyyyyy, ts_xxxyz_xxxxxx, ts_xxxyz_xxxxxy, ts_xxxyz_xxxxxz, ts_xxxyz_xxxxyy, ts_xxxyz_xxxxyz, ts_xxxyz_xxxxzz, ts_xxxyz_xxxyyy, ts_xxxyz_xxxyyz, ts_xxxyz_xxxyzz, ts_xxxyz_xxxzzz, ts_xxxyz_xxyyyy, ts_xxxyz_xxyyyz, ts_xxxyz_xxyyzz, ts_xxxyz_xxyzzz, ts_xxxyz_xxzzzz, ts_xxxyz_xyyyyy, ts_xxxyz_xyyyyz, ts_xxxyz_xyyyzz, ts_xxxyz_xyyzzz, ts_xxxyz_xyzzzz, ts_xxxyz_xzzzzz, ts_xxxyz_yyyyyy, ts_xxxyz_yyyyyz, ts_xxxyz_yyyyzz, ts_xxxyz_yyyzzz, ts_xxxyz_yyzzzz, ts_xxxyz_yzzzzz, ts_xxxyz_zzzzzz, ts_xxxz_xxxxxx, ts_xxxz_xxxxxz, ts_xxxz_xxxxyz, ts_xxxz_xxxxz, ts_xxxz_xxxxzz, ts_xxxz_xxxyyz, ts_xxxz_xxxyz, ts_xxxz_xxxyzz, ts_xxxz_xxxzz, ts_xxxz_xxxzzz, ts_xxxz_xxyyyz, ts_xxxz_xxyyz, ts_xxxz_xxyyzz, ts_xxxz_xxyzz, ts_xxxz_xxyzzz, ts_xxxz_xxzzz, ts_xxxz_xxzzzz, ts_xxxz_xyyyyz, ts_xxxz_xyyyz, ts_xxxz_xyyyzz, ts_xxxz_xyyzz, ts_xxxz_xyyzzz, ts_xxxz_xyzzz, ts_xxxz_xyzzzz, ts_xxxz_xzzzz, ts_xxxz_xzzzzz, ts_xxxz_zzzzzz, ts_xxyz_yyyyyz, ts_xxyz_yyyyzz, ts_xxyz_yyyzzz, ts_xxyz_yyzzzz, ts_xxyz_yzzzzz, ts_xyz_yyyyyz, ts_xyz_yyyyzz, ts_xyz_yyyzzz, ts_xyz_yyzzzz, ts_xyz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyz_xxxxxx[i] = ts_xxxz_xxxxxx[i] * pa_y[i];

        ts_xxxyz_xxxxxy[i] = ts_xxxy_xxxxxy[i] * pa_z[i];

        ts_xxxyz_xxxxxz[i] = ts_xxxz_xxxxxz[i] * pa_y[i];

        ts_xxxyz_xxxxyy[i] = ts_xxxy_xxxxyy[i] * pa_z[i];

        ts_xxxyz_xxxxyz[i] = ts_xxxz_xxxxz[i] * fe_0 + ts_xxxz_xxxxyz[i] * pa_y[i];

        ts_xxxyz_xxxxzz[i] = ts_xxxz_xxxxzz[i] * pa_y[i];

        ts_xxxyz_xxxyyy[i] = ts_xxxy_xxxyyy[i] * pa_z[i];

        ts_xxxyz_xxxyyz[i] = 2.0 * ts_xxxz_xxxyz[i] * fe_0 + ts_xxxz_xxxyyz[i] * pa_y[i];

        ts_xxxyz_xxxyzz[i] = ts_xxxz_xxxzz[i] * fe_0 + ts_xxxz_xxxyzz[i] * pa_y[i];

        ts_xxxyz_xxxzzz[i] = ts_xxxz_xxxzzz[i] * pa_y[i];

        ts_xxxyz_xxyyyy[i] = ts_xxxy_xxyyyy[i] * pa_z[i];

        ts_xxxyz_xxyyyz[i] = 3.0 * ts_xxxz_xxyyz[i] * fe_0 + ts_xxxz_xxyyyz[i] * pa_y[i];

        ts_xxxyz_xxyyzz[i] = 2.0 * ts_xxxz_xxyzz[i] * fe_0 + ts_xxxz_xxyyzz[i] * pa_y[i];

        ts_xxxyz_xxyzzz[i] = ts_xxxz_xxzzz[i] * fe_0 + ts_xxxz_xxyzzz[i] * pa_y[i];

        ts_xxxyz_xxzzzz[i] = ts_xxxz_xxzzzz[i] * pa_y[i];

        ts_xxxyz_xyyyyy[i] = ts_xxxy_xyyyyy[i] * pa_z[i];

        ts_xxxyz_xyyyyz[i] = 4.0 * ts_xxxz_xyyyz[i] * fe_0 + ts_xxxz_xyyyyz[i] * pa_y[i];

        ts_xxxyz_xyyyzz[i] = 3.0 * ts_xxxz_xyyzz[i] * fe_0 + ts_xxxz_xyyyzz[i] * pa_y[i];

        ts_xxxyz_xyyzzz[i] = 2.0 * ts_xxxz_xyzzz[i] * fe_0 + ts_xxxz_xyyzzz[i] * pa_y[i];

        ts_xxxyz_xyzzzz[i] = ts_xxxz_xzzzz[i] * fe_0 + ts_xxxz_xyzzzz[i] * pa_y[i];

        ts_xxxyz_xzzzzz[i] = ts_xxxz_xzzzzz[i] * pa_y[i];

        ts_xxxyz_yyyyyy[i] = ts_xxxy_yyyyyy[i] * pa_z[i];

        ts_xxxyz_yyyyyz[i] = 2.0 * ts_xyz_yyyyyz[i] * fe_0 + ts_xxyz_yyyyyz[i] * pa_x[i];

        ts_xxxyz_yyyyzz[i] = 2.0 * ts_xyz_yyyyzz[i] * fe_0 + ts_xxyz_yyyyzz[i] * pa_x[i];

        ts_xxxyz_yyyzzz[i] = 2.0 * ts_xyz_yyyzzz[i] * fe_0 + ts_xxyz_yyyzzz[i] * pa_x[i];

        ts_xxxyz_yyzzzz[i] = 2.0 * ts_xyz_yyzzzz[i] * fe_0 + ts_xxyz_yyzzzz[i] * pa_x[i];

        ts_xxxyz_yzzzzz[i] = 2.0 * ts_xyz_yzzzzz[i] * fe_0 + ts_xxyz_yzzzzz[i] * pa_x[i];

        ts_xxxyz_zzzzzz[i] = ts_xxxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : HI

    auto ts_xxxzz_xxxxxx = pbuffer.data(idx_ovl_hi + 140);

    auto ts_xxxzz_xxxxxy = pbuffer.data(idx_ovl_hi + 141);

    auto ts_xxxzz_xxxxxz = pbuffer.data(idx_ovl_hi + 142);

    auto ts_xxxzz_xxxxyy = pbuffer.data(idx_ovl_hi + 143);

    auto ts_xxxzz_xxxxyz = pbuffer.data(idx_ovl_hi + 144);

    auto ts_xxxzz_xxxxzz = pbuffer.data(idx_ovl_hi + 145);

    auto ts_xxxzz_xxxyyy = pbuffer.data(idx_ovl_hi + 146);

    auto ts_xxxzz_xxxyyz = pbuffer.data(idx_ovl_hi + 147);

    auto ts_xxxzz_xxxyzz = pbuffer.data(idx_ovl_hi + 148);

    auto ts_xxxzz_xxxzzz = pbuffer.data(idx_ovl_hi + 149);

    auto ts_xxxzz_xxyyyy = pbuffer.data(idx_ovl_hi + 150);

    auto ts_xxxzz_xxyyyz = pbuffer.data(idx_ovl_hi + 151);

    auto ts_xxxzz_xxyyzz = pbuffer.data(idx_ovl_hi + 152);

    auto ts_xxxzz_xxyzzz = pbuffer.data(idx_ovl_hi + 153);

    auto ts_xxxzz_xxzzzz = pbuffer.data(idx_ovl_hi + 154);

    auto ts_xxxzz_xyyyyy = pbuffer.data(idx_ovl_hi + 155);

    auto ts_xxxzz_xyyyyz = pbuffer.data(idx_ovl_hi + 156);

    auto ts_xxxzz_xyyyzz = pbuffer.data(idx_ovl_hi + 157);

    auto ts_xxxzz_xyyzzz = pbuffer.data(idx_ovl_hi + 158);

    auto ts_xxxzz_xyzzzz = pbuffer.data(idx_ovl_hi + 159);

    auto ts_xxxzz_xzzzzz = pbuffer.data(idx_ovl_hi + 160);

    auto ts_xxxzz_yyyyyy = pbuffer.data(idx_ovl_hi + 161);

    auto ts_xxxzz_yyyyyz = pbuffer.data(idx_ovl_hi + 162);

    auto ts_xxxzz_yyyyzz = pbuffer.data(idx_ovl_hi + 163);

    auto ts_xxxzz_yyyzzz = pbuffer.data(idx_ovl_hi + 164);

    auto ts_xxxzz_yyzzzz = pbuffer.data(idx_ovl_hi + 165);

    auto ts_xxxzz_yzzzzz = pbuffer.data(idx_ovl_hi + 166);

    auto ts_xxxzz_zzzzzz = pbuffer.data(idx_ovl_hi + 167);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxx_xxxxxx, ts_xxx_xxxxxy, ts_xxx_xxxxyy, ts_xxx_xxxyyy, ts_xxx_xxyyyy, ts_xxx_xyyyyy, ts_xxxz_xxxxxx, ts_xxxz_xxxxxy, ts_xxxz_xxxxyy, ts_xxxz_xxxyyy, ts_xxxz_xxyyyy, ts_xxxz_xyyyyy, ts_xxxzz_xxxxxx, ts_xxxzz_xxxxxy, ts_xxxzz_xxxxxz, ts_xxxzz_xxxxyy, ts_xxxzz_xxxxyz, ts_xxxzz_xxxxzz, ts_xxxzz_xxxyyy, ts_xxxzz_xxxyyz, ts_xxxzz_xxxyzz, ts_xxxzz_xxxzzz, ts_xxxzz_xxyyyy, ts_xxxzz_xxyyyz, ts_xxxzz_xxyyzz, ts_xxxzz_xxyzzz, ts_xxxzz_xxzzzz, ts_xxxzz_xyyyyy, ts_xxxzz_xyyyyz, ts_xxxzz_xyyyzz, ts_xxxzz_xyyzzz, ts_xxxzz_xyzzzz, ts_xxxzz_xzzzzz, ts_xxxzz_yyyyyy, ts_xxxzz_yyyyyz, ts_xxxzz_yyyyzz, ts_xxxzz_yyyzzz, ts_xxxzz_yyzzzz, ts_xxxzz_yzzzzz, ts_xxxzz_zzzzzz, ts_xxzz_xxxxxz, ts_xxzz_xxxxyz, ts_xxzz_xxxxz, ts_xxzz_xxxxzz, ts_xxzz_xxxyyz, ts_xxzz_xxxyz, ts_xxzz_xxxyzz, ts_xxzz_xxxzz, ts_xxzz_xxxzzz, ts_xxzz_xxyyyz, ts_xxzz_xxyyz, ts_xxzz_xxyyzz, ts_xxzz_xxyzz, ts_xxzz_xxyzzz, ts_xxzz_xxzzz, ts_xxzz_xxzzzz, ts_xxzz_xyyyyz, ts_xxzz_xyyyz, ts_xxzz_xyyyzz, ts_xxzz_xyyzz, ts_xxzz_xyyzzz, ts_xxzz_xyzzz, ts_xxzz_xyzzzz, ts_xxzz_xzzzz, ts_xxzz_xzzzzz, ts_xxzz_yyyyyy, ts_xxzz_yyyyyz, ts_xxzz_yyyyz, ts_xxzz_yyyyzz, ts_xxzz_yyyzz, ts_xxzz_yyyzzz, ts_xxzz_yyzzz, ts_xxzz_yyzzzz, ts_xxzz_yzzzz, ts_xxzz_yzzzzz, ts_xxzz_zzzzz, ts_xxzz_zzzzzz, ts_xzz_xxxxxz, ts_xzz_xxxxyz, ts_xzz_xxxxzz, ts_xzz_xxxyyz, ts_xzz_xxxyzz, ts_xzz_xxxzzz, ts_xzz_xxyyyz, ts_xzz_xxyyzz, ts_xzz_xxyzzz, ts_xzz_xxzzzz, ts_xzz_xyyyyz, ts_xzz_xyyyzz, ts_xzz_xyyzzz, ts_xzz_xyzzzz, ts_xzz_xzzzzz, ts_xzz_yyyyyy, ts_xzz_yyyyyz, ts_xzz_yyyyzz, ts_xzz_yyyzzz, ts_xzz_yyzzzz, ts_xzz_yzzzzz, ts_xzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzz_xxxxxx[i] = ts_xxx_xxxxxx[i] * fe_0 + ts_xxxz_xxxxxx[i] * pa_z[i];

        ts_xxxzz_xxxxxy[i] = ts_xxx_xxxxxy[i] * fe_0 + ts_xxxz_xxxxxy[i] * pa_z[i];

        ts_xxxzz_xxxxxz[i] = 2.0 * ts_xzz_xxxxxz[i] * fe_0 + 5.0 * ts_xxzz_xxxxz[i] * fe_0 + ts_xxzz_xxxxxz[i] * pa_x[i];

        ts_xxxzz_xxxxyy[i] = ts_xxx_xxxxyy[i] * fe_0 + ts_xxxz_xxxxyy[i] * pa_z[i];

        ts_xxxzz_xxxxyz[i] = 2.0 * ts_xzz_xxxxyz[i] * fe_0 + 4.0 * ts_xxzz_xxxyz[i] * fe_0 + ts_xxzz_xxxxyz[i] * pa_x[i];

        ts_xxxzz_xxxxzz[i] = 2.0 * ts_xzz_xxxxzz[i] * fe_0 + 4.0 * ts_xxzz_xxxzz[i] * fe_0 + ts_xxzz_xxxxzz[i] * pa_x[i];

        ts_xxxzz_xxxyyy[i] = ts_xxx_xxxyyy[i] * fe_0 + ts_xxxz_xxxyyy[i] * pa_z[i];

        ts_xxxzz_xxxyyz[i] = 2.0 * ts_xzz_xxxyyz[i] * fe_0 + 3.0 * ts_xxzz_xxyyz[i] * fe_0 + ts_xxzz_xxxyyz[i] * pa_x[i];

        ts_xxxzz_xxxyzz[i] = 2.0 * ts_xzz_xxxyzz[i] * fe_0 + 3.0 * ts_xxzz_xxyzz[i] * fe_0 + ts_xxzz_xxxyzz[i] * pa_x[i];

        ts_xxxzz_xxxzzz[i] = 2.0 * ts_xzz_xxxzzz[i] * fe_0 + 3.0 * ts_xxzz_xxzzz[i] * fe_0 + ts_xxzz_xxxzzz[i] * pa_x[i];

        ts_xxxzz_xxyyyy[i] = ts_xxx_xxyyyy[i] * fe_0 + ts_xxxz_xxyyyy[i] * pa_z[i];

        ts_xxxzz_xxyyyz[i] = 2.0 * ts_xzz_xxyyyz[i] * fe_0 + 2.0 * ts_xxzz_xyyyz[i] * fe_0 + ts_xxzz_xxyyyz[i] * pa_x[i];

        ts_xxxzz_xxyyzz[i] = 2.0 * ts_xzz_xxyyzz[i] * fe_0 + 2.0 * ts_xxzz_xyyzz[i] * fe_0 + ts_xxzz_xxyyzz[i] * pa_x[i];

        ts_xxxzz_xxyzzz[i] = 2.0 * ts_xzz_xxyzzz[i] * fe_0 + 2.0 * ts_xxzz_xyzzz[i] * fe_0 + ts_xxzz_xxyzzz[i] * pa_x[i];

        ts_xxxzz_xxzzzz[i] = 2.0 * ts_xzz_xxzzzz[i] * fe_0 + 2.0 * ts_xxzz_xzzzz[i] * fe_0 + ts_xxzz_xxzzzz[i] * pa_x[i];

        ts_xxxzz_xyyyyy[i] = ts_xxx_xyyyyy[i] * fe_0 + ts_xxxz_xyyyyy[i] * pa_z[i];

        ts_xxxzz_xyyyyz[i] = 2.0 * ts_xzz_xyyyyz[i] * fe_0 + ts_xxzz_yyyyz[i] * fe_0 + ts_xxzz_xyyyyz[i] * pa_x[i];

        ts_xxxzz_xyyyzz[i] = 2.0 * ts_xzz_xyyyzz[i] * fe_0 + ts_xxzz_yyyzz[i] * fe_0 + ts_xxzz_xyyyzz[i] * pa_x[i];

        ts_xxxzz_xyyzzz[i] = 2.0 * ts_xzz_xyyzzz[i] * fe_0 + ts_xxzz_yyzzz[i] * fe_0 + ts_xxzz_xyyzzz[i] * pa_x[i];

        ts_xxxzz_xyzzzz[i] = 2.0 * ts_xzz_xyzzzz[i] * fe_0 + ts_xxzz_yzzzz[i] * fe_0 + ts_xxzz_xyzzzz[i] * pa_x[i];

        ts_xxxzz_xzzzzz[i] = 2.0 * ts_xzz_xzzzzz[i] * fe_0 + ts_xxzz_zzzzz[i] * fe_0 + ts_xxzz_xzzzzz[i] * pa_x[i];

        ts_xxxzz_yyyyyy[i] = 2.0 * ts_xzz_yyyyyy[i] * fe_0 + ts_xxzz_yyyyyy[i] * pa_x[i];

        ts_xxxzz_yyyyyz[i] = 2.0 * ts_xzz_yyyyyz[i] * fe_0 + ts_xxzz_yyyyyz[i] * pa_x[i];

        ts_xxxzz_yyyyzz[i] = 2.0 * ts_xzz_yyyyzz[i] * fe_0 + ts_xxzz_yyyyzz[i] * pa_x[i];

        ts_xxxzz_yyyzzz[i] = 2.0 * ts_xzz_yyyzzz[i] * fe_0 + ts_xxzz_yyyzzz[i] * pa_x[i];

        ts_xxxzz_yyzzzz[i] = 2.0 * ts_xzz_yyzzzz[i] * fe_0 + ts_xxzz_yyzzzz[i] * pa_x[i];

        ts_xxxzz_yzzzzz[i] = 2.0 * ts_xzz_yzzzzz[i] * fe_0 + ts_xxzz_yzzzzz[i] * pa_x[i];

        ts_xxxzz_zzzzzz[i] = 2.0 * ts_xzz_zzzzzz[i] * fe_0 + ts_xxzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : HI

    auto ts_xxyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 168);

    auto ts_xxyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 169);

    auto ts_xxyyy_xxxxxz = pbuffer.data(idx_ovl_hi + 170);

    auto ts_xxyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 171);

    auto ts_xxyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 172);

    auto ts_xxyyy_xxxxzz = pbuffer.data(idx_ovl_hi + 173);

    auto ts_xxyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 174);

    auto ts_xxyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 175);

    auto ts_xxyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 176);

    auto ts_xxyyy_xxxzzz = pbuffer.data(idx_ovl_hi + 177);

    auto ts_xxyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 178);

    auto ts_xxyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 179);

    auto ts_xxyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 180);

    auto ts_xxyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 181);

    auto ts_xxyyy_xxzzzz = pbuffer.data(idx_ovl_hi + 182);

    auto ts_xxyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 183);

    auto ts_xxyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 184);

    auto ts_xxyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 185);

    auto ts_xxyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 186);

    auto ts_xxyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 187);

    auto ts_xxyyy_xzzzzz = pbuffer.data(idx_ovl_hi + 188);

    auto ts_xxyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 189);

    auto ts_xxyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 190);

    auto ts_xxyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 191);

    auto ts_xxyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 192);

    auto ts_xxyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 193);

    auto ts_xxyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 194);

    auto ts_xxyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 195);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxy_xxxxxx, ts_xxy_xxxxxz, ts_xxy_xxxxzz, ts_xxy_xxxzzz, ts_xxy_xxzzzz, ts_xxy_xzzzzz, ts_xxyy_xxxxxx, ts_xxyy_xxxxxz, ts_xxyy_xxxxzz, ts_xxyy_xxxzzz, ts_xxyy_xxzzzz, ts_xxyy_xzzzzz, ts_xxyyy_xxxxxx, ts_xxyyy_xxxxxy, ts_xxyyy_xxxxxz, ts_xxyyy_xxxxyy, ts_xxyyy_xxxxyz, ts_xxyyy_xxxxzz, ts_xxyyy_xxxyyy, ts_xxyyy_xxxyyz, ts_xxyyy_xxxyzz, ts_xxyyy_xxxzzz, ts_xxyyy_xxyyyy, ts_xxyyy_xxyyyz, ts_xxyyy_xxyyzz, ts_xxyyy_xxyzzz, ts_xxyyy_xxzzzz, ts_xxyyy_xyyyyy, ts_xxyyy_xyyyyz, ts_xxyyy_xyyyzz, ts_xxyyy_xyyzzz, ts_xxyyy_xyzzzz, ts_xxyyy_xzzzzz, ts_xxyyy_yyyyyy, ts_xxyyy_yyyyyz, ts_xxyyy_yyyyzz, ts_xxyyy_yyyzzz, ts_xxyyy_yyzzzz, ts_xxyyy_yzzzzz, ts_xxyyy_zzzzzz, ts_xyyy_xxxxxy, ts_xyyy_xxxxy, ts_xyyy_xxxxyy, ts_xyyy_xxxxyz, ts_xyyy_xxxyy, ts_xyyy_xxxyyy, ts_xyyy_xxxyyz, ts_xyyy_xxxyz, ts_xyyy_xxxyzz, ts_xyyy_xxyyy, ts_xyyy_xxyyyy, ts_xyyy_xxyyyz, ts_xyyy_xxyyz, ts_xyyy_xxyyzz, ts_xyyy_xxyzz, ts_xyyy_xxyzzz, ts_xyyy_xyyyy, ts_xyyy_xyyyyy, ts_xyyy_xyyyyz, ts_xyyy_xyyyz, ts_xyyy_xyyyzz, ts_xyyy_xyyzz, ts_xyyy_xyyzzz, ts_xyyy_xyzzz, ts_xyyy_xyzzzz, ts_xyyy_yyyyy, ts_xyyy_yyyyyy, ts_xyyy_yyyyyz, ts_xyyy_yyyyz, ts_xyyy_yyyyzz, ts_xyyy_yyyzz, ts_xyyy_yyyzzz, ts_xyyy_yyzzz, ts_xyyy_yyzzzz, ts_xyyy_yzzzz, ts_xyyy_yzzzzz, ts_xyyy_zzzzzz, ts_yyy_xxxxxy, ts_yyy_xxxxyy, ts_yyy_xxxxyz, ts_yyy_xxxyyy, ts_yyy_xxxyyz, ts_yyy_xxxyzz, ts_yyy_xxyyyy, ts_yyy_xxyyyz, ts_yyy_xxyyzz, ts_yyy_xxyzzz, ts_yyy_xyyyyy, ts_yyy_xyyyyz, ts_yyy_xyyyzz, ts_yyy_xyyzzz, ts_yyy_xyzzzz, ts_yyy_yyyyyy, ts_yyy_yyyyyz, ts_yyy_yyyyzz, ts_yyy_yyyzzz, ts_yyy_yyzzzz, ts_yyy_yzzzzz, ts_yyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyy_xxxxxx[i] = 2.0 * ts_xxy_xxxxxx[i] * fe_0 + ts_xxyy_xxxxxx[i] * pa_y[i];

        ts_xxyyy_xxxxxy[i] = ts_yyy_xxxxxy[i] * fe_0 + 5.0 * ts_xyyy_xxxxy[i] * fe_0 + ts_xyyy_xxxxxy[i] * pa_x[i];

        ts_xxyyy_xxxxxz[i] = 2.0 * ts_xxy_xxxxxz[i] * fe_0 + ts_xxyy_xxxxxz[i] * pa_y[i];

        ts_xxyyy_xxxxyy[i] = ts_yyy_xxxxyy[i] * fe_0 + 4.0 * ts_xyyy_xxxyy[i] * fe_0 + ts_xyyy_xxxxyy[i] * pa_x[i];

        ts_xxyyy_xxxxyz[i] = ts_yyy_xxxxyz[i] * fe_0 + 4.0 * ts_xyyy_xxxyz[i] * fe_0 + ts_xyyy_xxxxyz[i] * pa_x[i];

        ts_xxyyy_xxxxzz[i] = 2.0 * ts_xxy_xxxxzz[i] * fe_0 + ts_xxyy_xxxxzz[i] * pa_y[i];

        ts_xxyyy_xxxyyy[i] = ts_yyy_xxxyyy[i] * fe_0 + 3.0 * ts_xyyy_xxyyy[i] * fe_0 + ts_xyyy_xxxyyy[i] * pa_x[i];

        ts_xxyyy_xxxyyz[i] = ts_yyy_xxxyyz[i] * fe_0 + 3.0 * ts_xyyy_xxyyz[i] * fe_0 + ts_xyyy_xxxyyz[i] * pa_x[i];

        ts_xxyyy_xxxyzz[i] = ts_yyy_xxxyzz[i] * fe_0 + 3.0 * ts_xyyy_xxyzz[i] * fe_0 + ts_xyyy_xxxyzz[i] * pa_x[i];

        ts_xxyyy_xxxzzz[i] = 2.0 * ts_xxy_xxxzzz[i] * fe_0 + ts_xxyy_xxxzzz[i] * pa_y[i];

        ts_xxyyy_xxyyyy[i] = ts_yyy_xxyyyy[i] * fe_0 + 2.0 * ts_xyyy_xyyyy[i] * fe_0 + ts_xyyy_xxyyyy[i] * pa_x[i];

        ts_xxyyy_xxyyyz[i] = ts_yyy_xxyyyz[i] * fe_0 + 2.0 * ts_xyyy_xyyyz[i] * fe_0 + ts_xyyy_xxyyyz[i] * pa_x[i];

        ts_xxyyy_xxyyzz[i] = ts_yyy_xxyyzz[i] * fe_0 + 2.0 * ts_xyyy_xyyzz[i] * fe_0 + ts_xyyy_xxyyzz[i] * pa_x[i];

        ts_xxyyy_xxyzzz[i] = ts_yyy_xxyzzz[i] * fe_0 + 2.0 * ts_xyyy_xyzzz[i] * fe_0 + ts_xyyy_xxyzzz[i] * pa_x[i];

        ts_xxyyy_xxzzzz[i] = 2.0 * ts_xxy_xxzzzz[i] * fe_0 + ts_xxyy_xxzzzz[i] * pa_y[i];

        ts_xxyyy_xyyyyy[i] = ts_yyy_xyyyyy[i] * fe_0 + ts_xyyy_yyyyy[i] * fe_0 + ts_xyyy_xyyyyy[i] * pa_x[i];

        ts_xxyyy_xyyyyz[i] = ts_yyy_xyyyyz[i] * fe_0 + ts_xyyy_yyyyz[i] * fe_0 + ts_xyyy_xyyyyz[i] * pa_x[i];

        ts_xxyyy_xyyyzz[i] = ts_yyy_xyyyzz[i] * fe_0 + ts_xyyy_yyyzz[i] * fe_0 + ts_xyyy_xyyyzz[i] * pa_x[i];

        ts_xxyyy_xyyzzz[i] = ts_yyy_xyyzzz[i] * fe_0 + ts_xyyy_yyzzz[i] * fe_0 + ts_xyyy_xyyzzz[i] * pa_x[i];

        ts_xxyyy_xyzzzz[i] = ts_yyy_xyzzzz[i] * fe_0 + ts_xyyy_yzzzz[i] * fe_0 + ts_xyyy_xyzzzz[i] * pa_x[i];

        ts_xxyyy_xzzzzz[i] = 2.0 * ts_xxy_xzzzzz[i] * fe_0 + ts_xxyy_xzzzzz[i] * pa_y[i];

        ts_xxyyy_yyyyyy[i] = ts_yyy_yyyyyy[i] * fe_0 + ts_xyyy_yyyyyy[i] * pa_x[i];

        ts_xxyyy_yyyyyz[i] = ts_yyy_yyyyyz[i] * fe_0 + ts_xyyy_yyyyyz[i] * pa_x[i];

        ts_xxyyy_yyyyzz[i] = ts_yyy_yyyyzz[i] * fe_0 + ts_xyyy_yyyyzz[i] * pa_x[i];

        ts_xxyyy_yyyzzz[i] = ts_yyy_yyyzzz[i] * fe_0 + ts_xyyy_yyyzzz[i] * pa_x[i];

        ts_xxyyy_yyzzzz[i] = ts_yyy_yyzzzz[i] * fe_0 + ts_xyyy_yyzzzz[i] * pa_x[i];

        ts_xxyyy_yzzzzz[i] = ts_yyy_yzzzzz[i] * fe_0 + ts_xyyy_yzzzzz[i] * pa_x[i];

        ts_xxyyy_zzzzzz[i] = ts_yyy_zzzzzz[i] * fe_0 + ts_xyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : HI

    auto ts_xxyyz_xxxxxx = pbuffer.data(idx_ovl_hi + 196);

    auto ts_xxyyz_xxxxxy = pbuffer.data(idx_ovl_hi + 197);

    auto ts_xxyyz_xxxxxz = pbuffer.data(idx_ovl_hi + 198);

    auto ts_xxyyz_xxxxyy = pbuffer.data(idx_ovl_hi + 199);

    auto ts_xxyyz_xxxxyz = pbuffer.data(idx_ovl_hi + 200);

    auto ts_xxyyz_xxxxzz = pbuffer.data(idx_ovl_hi + 201);

    auto ts_xxyyz_xxxyyy = pbuffer.data(idx_ovl_hi + 202);

    auto ts_xxyyz_xxxyyz = pbuffer.data(idx_ovl_hi + 203);

    auto ts_xxyyz_xxxyzz = pbuffer.data(idx_ovl_hi + 204);

    auto ts_xxyyz_xxxzzz = pbuffer.data(idx_ovl_hi + 205);

    auto ts_xxyyz_xxyyyy = pbuffer.data(idx_ovl_hi + 206);

    auto ts_xxyyz_xxyyyz = pbuffer.data(idx_ovl_hi + 207);

    auto ts_xxyyz_xxyyzz = pbuffer.data(idx_ovl_hi + 208);

    auto ts_xxyyz_xxyzzz = pbuffer.data(idx_ovl_hi + 209);

    auto ts_xxyyz_xxzzzz = pbuffer.data(idx_ovl_hi + 210);

    auto ts_xxyyz_xyyyyy = pbuffer.data(idx_ovl_hi + 211);

    auto ts_xxyyz_xyyyyz = pbuffer.data(idx_ovl_hi + 212);

    auto ts_xxyyz_xyyyzz = pbuffer.data(idx_ovl_hi + 213);

    auto ts_xxyyz_xyyzzz = pbuffer.data(idx_ovl_hi + 214);

    auto ts_xxyyz_xyzzzz = pbuffer.data(idx_ovl_hi + 215);

    auto ts_xxyyz_xzzzzz = pbuffer.data(idx_ovl_hi + 216);

    auto ts_xxyyz_yyyyyy = pbuffer.data(idx_ovl_hi + 217);

    auto ts_xxyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 218);

    auto ts_xxyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 219);

    auto ts_xxyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 220);

    auto ts_xxyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 221);

    auto ts_xxyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 222);

    auto ts_xxyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 223);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyy_xxxxxx, ts_xxyy_xxxxxy, ts_xxyy_xxxxy, ts_xxyy_xxxxyy, ts_xxyy_xxxxyz, ts_xxyy_xxxyy, ts_xxyy_xxxyyy, ts_xxyy_xxxyyz, ts_xxyy_xxxyz, ts_xxyy_xxxyzz, ts_xxyy_xxyyy, ts_xxyy_xxyyyy, ts_xxyy_xxyyyz, ts_xxyy_xxyyz, ts_xxyy_xxyyzz, ts_xxyy_xxyzz, ts_xxyy_xxyzzz, ts_xxyy_xyyyy, ts_xxyy_xyyyyy, ts_xxyy_xyyyyz, ts_xxyy_xyyyz, ts_xxyy_xyyyzz, ts_xxyy_xyyzz, ts_xxyy_xyyzzz, ts_xxyy_xyzzz, ts_xxyy_xyzzzz, ts_xxyy_yyyyyy, ts_xxyyz_xxxxxx, ts_xxyyz_xxxxxy, ts_xxyyz_xxxxxz, ts_xxyyz_xxxxyy, ts_xxyyz_xxxxyz, ts_xxyyz_xxxxzz, ts_xxyyz_xxxyyy, ts_xxyyz_xxxyyz, ts_xxyyz_xxxyzz, ts_xxyyz_xxxzzz, ts_xxyyz_xxyyyy, ts_xxyyz_xxyyyz, ts_xxyyz_xxyyzz, ts_xxyyz_xxyzzz, ts_xxyyz_xxzzzz, ts_xxyyz_xyyyyy, ts_xxyyz_xyyyyz, ts_xxyyz_xyyyzz, ts_xxyyz_xyyzzz, ts_xxyyz_xyzzzz, ts_xxyyz_xzzzzz, ts_xxyyz_yyyyyy, ts_xxyyz_yyyyyz, ts_xxyyz_yyyyzz, ts_xxyyz_yyyzzz, ts_xxyyz_yyzzzz, ts_xxyyz_yzzzzz, ts_xxyyz_zzzzzz, ts_xxyz_xxxxxz, ts_xxyz_xxxxzz, ts_xxyz_xxxzzz, ts_xxyz_xxzzzz, ts_xxyz_xzzzzz, ts_xxz_xxxxxz, ts_xxz_xxxxzz, ts_xxz_xxxzzz, ts_xxz_xxzzzz, ts_xxz_xzzzzz, ts_xyyz_yyyyyz, ts_xyyz_yyyyzz, ts_xyyz_yyyzzz, ts_xyyz_yyzzzz, ts_xyyz_yzzzzz, ts_xyyz_zzzzzz, ts_yyz_yyyyyz, ts_yyz_yyyyzz, ts_yyz_yyyzzz, ts_yyz_yyzzzz, ts_yyz_yzzzzz, ts_yyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyz_xxxxxx[i] = ts_xxyy_xxxxxx[i] * pa_z[i];

        ts_xxyyz_xxxxxy[i] = ts_xxyy_xxxxxy[i] * pa_z[i];

        ts_xxyyz_xxxxxz[i] = ts_xxz_xxxxxz[i] * fe_0 + ts_xxyz_xxxxxz[i] * pa_y[i];

        ts_xxyyz_xxxxyy[i] = ts_xxyy_xxxxyy[i] * pa_z[i];

        ts_xxyyz_xxxxyz[i] = ts_xxyy_xxxxy[i] * fe_0 + ts_xxyy_xxxxyz[i] * pa_z[i];

        ts_xxyyz_xxxxzz[i] = ts_xxz_xxxxzz[i] * fe_0 + ts_xxyz_xxxxzz[i] * pa_y[i];

        ts_xxyyz_xxxyyy[i] = ts_xxyy_xxxyyy[i] * pa_z[i];

        ts_xxyyz_xxxyyz[i] = ts_xxyy_xxxyy[i] * fe_0 + ts_xxyy_xxxyyz[i] * pa_z[i];

        ts_xxyyz_xxxyzz[i] = 2.0 * ts_xxyy_xxxyz[i] * fe_0 + ts_xxyy_xxxyzz[i] * pa_z[i];

        ts_xxyyz_xxxzzz[i] = ts_xxz_xxxzzz[i] * fe_0 + ts_xxyz_xxxzzz[i] * pa_y[i];

        ts_xxyyz_xxyyyy[i] = ts_xxyy_xxyyyy[i] * pa_z[i];

        ts_xxyyz_xxyyyz[i] = ts_xxyy_xxyyy[i] * fe_0 + ts_xxyy_xxyyyz[i] * pa_z[i];

        ts_xxyyz_xxyyzz[i] = 2.0 * ts_xxyy_xxyyz[i] * fe_0 + ts_xxyy_xxyyzz[i] * pa_z[i];

        ts_xxyyz_xxyzzz[i] = 3.0 * ts_xxyy_xxyzz[i] * fe_0 + ts_xxyy_xxyzzz[i] * pa_z[i];

        ts_xxyyz_xxzzzz[i] = ts_xxz_xxzzzz[i] * fe_0 + ts_xxyz_xxzzzz[i] * pa_y[i];

        ts_xxyyz_xyyyyy[i] = ts_xxyy_xyyyyy[i] * pa_z[i];

        ts_xxyyz_xyyyyz[i] = ts_xxyy_xyyyy[i] * fe_0 + ts_xxyy_xyyyyz[i] * pa_z[i];

        ts_xxyyz_xyyyzz[i] = 2.0 * ts_xxyy_xyyyz[i] * fe_0 + ts_xxyy_xyyyzz[i] * pa_z[i];

        ts_xxyyz_xyyzzz[i] = 3.0 * ts_xxyy_xyyzz[i] * fe_0 + ts_xxyy_xyyzzz[i] * pa_z[i];

        ts_xxyyz_xyzzzz[i] = 4.0 * ts_xxyy_xyzzz[i] * fe_0 + ts_xxyy_xyzzzz[i] * pa_z[i];

        ts_xxyyz_xzzzzz[i] = ts_xxz_xzzzzz[i] * fe_0 + ts_xxyz_xzzzzz[i] * pa_y[i];

        ts_xxyyz_yyyyyy[i] = ts_xxyy_yyyyyy[i] * pa_z[i];

        ts_xxyyz_yyyyyz[i] = ts_yyz_yyyyyz[i] * fe_0 + ts_xyyz_yyyyyz[i] * pa_x[i];

        ts_xxyyz_yyyyzz[i] = ts_yyz_yyyyzz[i] * fe_0 + ts_xyyz_yyyyzz[i] * pa_x[i];

        ts_xxyyz_yyyzzz[i] = ts_yyz_yyyzzz[i] * fe_0 + ts_xyyz_yyyzzz[i] * pa_x[i];

        ts_xxyyz_yyzzzz[i] = ts_yyz_yyzzzz[i] * fe_0 + ts_xyyz_yyzzzz[i] * pa_x[i];

        ts_xxyyz_yzzzzz[i] = ts_yyz_yzzzzz[i] * fe_0 + ts_xyyz_yzzzzz[i] * pa_x[i];

        ts_xxyyz_zzzzzz[i] = ts_yyz_zzzzzz[i] * fe_0 + ts_xyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 224-252 components of targeted buffer : HI

    auto ts_xxyzz_xxxxxx = pbuffer.data(idx_ovl_hi + 224);

    auto ts_xxyzz_xxxxxy = pbuffer.data(idx_ovl_hi + 225);

    auto ts_xxyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 226);

    auto ts_xxyzz_xxxxyy = pbuffer.data(idx_ovl_hi + 227);

    auto ts_xxyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 228);

    auto ts_xxyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 229);

    auto ts_xxyzz_xxxyyy = pbuffer.data(idx_ovl_hi + 230);

    auto ts_xxyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 231);

    auto ts_xxyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 232);

    auto ts_xxyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 233);

    auto ts_xxyzz_xxyyyy = pbuffer.data(idx_ovl_hi + 234);

    auto ts_xxyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 235);

    auto ts_xxyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 236);

    auto ts_xxyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 237);

    auto ts_xxyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 238);

    auto ts_xxyzz_xyyyyy = pbuffer.data(idx_ovl_hi + 239);

    auto ts_xxyzz_xyyyyz = pbuffer.data(idx_ovl_hi + 240);

    auto ts_xxyzz_xyyyzz = pbuffer.data(idx_ovl_hi + 241);

    auto ts_xxyzz_xyyzzz = pbuffer.data(idx_ovl_hi + 242);

    auto ts_xxyzz_xyzzzz = pbuffer.data(idx_ovl_hi + 243);

    auto ts_xxyzz_xzzzzz = pbuffer.data(idx_ovl_hi + 244);

    auto ts_xxyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 245);

    auto ts_xxyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 246);

    auto ts_xxyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 247);

    auto ts_xxyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 248);

    auto ts_xxyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 249);

    auto ts_xxyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 250);

    auto ts_xxyzz_zzzzzz = pbuffer.data(idx_ovl_hi + 251);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzz_xxxxxx, ts_xxyzz_xxxxxy, ts_xxyzz_xxxxxz, ts_xxyzz_xxxxyy, ts_xxyzz_xxxxyz, ts_xxyzz_xxxxzz, ts_xxyzz_xxxyyy, ts_xxyzz_xxxyyz, ts_xxyzz_xxxyzz, ts_xxyzz_xxxzzz, ts_xxyzz_xxyyyy, ts_xxyzz_xxyyyz, ts_xxyzz_xxyyzz, ts_xxyzz_xxyzzz, ts_xxyzz_xxzzzz, ts_xxyzz_xyyyyy, ts_xxyzz_xyyyyz, ts_xxyzz_xyyyzz, ts_xxyzz_xyyzzz, ts_xxyzz_xyzzzz, ts_xxyzz_xzzzzz, ts_xxyzz_yyyyyy, ts_xxyzz_yyyyyz, ts_xxyzz_yyyyzz, ts_xxyzz_yyyzzz, ts_xxyzz_yyzzzz, ts_xxyzz_yzzzzz, ts_xxyzz_zzzzzz, ts_xxzz_xxxxx, ts_xxzz_xxxxxx, ts_xxzz_xxxxxy, ts_xxzz_xxxxxz, ts_xxzz_xxxxy, ts_xxzz_xxxxyy, ts_xxzz_xxxxyz, ts_xxzz_xxxxz, ts_xxzz_xxxxzz, ts_xxzz_xxxyy, ts_xxzz_xxxyyy, ts_xxzz_xxxyyz, ts_xxzz_xxxyz, ts_xxzz_xxxyzz, ts_xxzz_xxxzz, ts_xxzz_xxxzzz, ts_xxzz_xxyyy, ts_xxzz_xxyyyy, ts_xxzz_xxyyyz, ts_xxzz_xxyyz, ts_xxzz_xxyyzz, ts_xxzz_xxyzz, ts_xxzz_xxyzzz, ts_xxzz_xxzzz, ts_xxzz_xxzzzz, ts_xxzz_xyyyy, ts_xxzz_xyyyyy, ts_xxzz_xyyyyz, ts_xxzz_xyyyz, ts_xxzz_xyyyzz, ts_xxzz_xyyzz, ts_xxzz_xyyzzz, ts_xxzz_xyzzz, ts_xxzz_xyzzzz, ts_xxzz_xzzzz, ts_xxzz_xzzzzz, ts_xxzz_zzzzzz, ts_xyzz_yyyyyy, ts_xyzz_yyyyyz, ts_xyzz_yyyyzz, ts_xyzz_yyyzzz, ts_xyzz_yyzzzz, ts_xyzz_yzzzzz, ts_yzz_yyyyyy, ts_yzz_yyyyyz, ts_yzz_yyyyzz, ts_yzz_yyyzzz, ts_yzz_yyzzzz, ts_yzz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzz_xxxxxx[i] = ts_xxzz_xxxxxx[i] * pa_y[i];

        ts_xxyzz_xxxxxy[i] = ts_xxzz_xxxxx[i] * fe_0 + ts_xxzz_xxxxxy[i] * pa_y[i];

        ts_xxyzz_xxxxxz[i] = ts_xxzz_xxxxxz[i] * pa_y[i];

        ts_xxyzz_xxxxyy[i] = 2.0 * ts_xxzz_xxxxy[i] * fe_0 + ts_xxzz_xxxxyy[i] * pa_y[i];

        ts_xxyzz_xxxxyz[i] = ts_xxzz_xxxxz[i] * fe_0 + ts_xxzz_xxxxyz[i] * pa_y[i];

        ts_xxyzz_xxxxzz[i] = ts_xxzz_xxxxzz[i] * pa_y[i];

        ts_xxyzz_xxxyyy[i] = 3.0 * ts_xxzz_xxxyy[i] * fe_0 + ts_xxzz_xxxyyy[i] * pa_y[i];

        ts_xxyzz_xxxyyz[i] = 2.0 * ts_xxzz_xxxyz[i] * fe_0 + ts_xxzz_xxxyyz[i] * pa_y[i];

        ts_xxyzz_xxxyzz[i] = ts_xxzz_xxxzz[i] * fe_0 + ts_xxzz_xxxyzz[i] * pa_y[i];

        ts_xxyzz_xxxzzz[i] = ts_xxzz_xxxzzz[i] * pa_y[i];

        ts_xxyzz_xxyyyy[i] = 4.0 * ts_xxzz_xxyyy[i] * fe_0 + ts_xxzz_xxyyyy[i] * pa_y[i];

        ts_xxyzz_xxyyyz[i] = 3.0 * ts_xxzz_xxyyz[i] * fe_0 + ts_xxzz_xxyyyz[i] * pa_y[i];

        ts_xxyzz_xxyyzz[i] = 2.0 * ts_xxzz_xxyzz[i] * fe_0 + ts_xxzz_xxyyzz[i] * pa_y[i];

        ts_xxyzz_xxyzzz[i] = ts_xxzz_xxzzz[i] * fe_0 + ts_xxzz_xxyzzz[i] * pa_y[i];

        ts_xxyzz_xxzzzz[i] = ts_xxzz_xxzzzz[i] * pa_y[i];

        ts_xxyzz_xyyyyy[i] = 5.0 * ts_xxzz_xyyyy[i] * fe_0 + ts_xxzz_xyyyyy[i] * pa_y[i];

        ts_xxyzz_xyyyyz[i] = 4.0 * ts_xxzz_xyyyz[i] * fe_0 + ts_xxzz_xyyyyz[i] * pa_y[i];

        ts_xxyzz_xyyyzz[i] = 3.0 * ts_xxzz_xyyzz[i] * fe_0 + ts_xxzz_xyyyzz[i] * pa_y[i];

        ts_xxyzz_xyyzzz[i] = 2.0 * ts_xxzz_xyzzz[i] * fe_0 + ts_xxzz_xyyzzz[i] * pa_y[i];

        ts_xxyzz_xyzzzz[i] = ts_xxzz_xzzzz[i] * fe_0 + ts_xxzz_xyzzzz[i] * pa_y[i];

        ts_xxyzz_xzzzzz[i] = ts_xxzz_xzzzzz[i] * pa_y[i];

        ts_xxyzz_yyyyyy[i] = ts_yzz_yyyyyy[i] * fe_0 + ts_xyzz_yyyyyy[i] * pa_x[i];

        ts_xxyzz_yyyyyz[i] = ts_yzz_yyyyyz[i] * fe_0 + ts_xyzz_yyyyyz[i] * pa_x[i];

        ts_xxyzz_yyyyzz[i] = ts_yzz_yyyyzz[i] * fe_0 + ts_xyzz_yyyyzz[i] * pa_x[i];

        ts_xxyzz_yyyzzz[i] = ts_yzz_yyyzzz[i] * fe_0 + ts_xyzz_yyyzzz[i] * pa_x[i];

        ts_xxyzz_yyzzzz[i] = ts_yzz_yyzzzz[i] * fe_0 + ts_xyzz_yyzzzz[i] * pa_x[i];

        ts_xxyzz_yzzzzz[i] = ts_yzz_yzzzzz[i] * fe_0 + ts_xyzz_yzzzzz[i] * pa_x[i];

        ts_xxyzz_zzzzzz[i] = ts_xxzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : HI

    auto ts_xxzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 252);

    auto ts_xxzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 253);

    auto ts_xxzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 254);

    auto ts_xxzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 255);

    auto ts_xxzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 256);

    auto ts_xxzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 257);

    auto ts_xxzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 258);

    auto ts_xxzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 259);

    auto ts_xxzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 260);

    auto ts_xxzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 261);

    auto ts_xxzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 262);

    auto ts_xxzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 263);

    auto ts_xxzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 264);

    auto ts_xxzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 265);

    auto ts_xxzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 266);

    auto ts_xxzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 267);

    auto ts_xxzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 268);

    auto ts_xxzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 269);

    auto ts_xxzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 270);

    auto ts_xxzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 271);

    auto ts_xxzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 272);

    auto ts_xxzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 273);

    auto ts_xxzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 274);

    auto ts_xxzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 275);

    auto ts_xxzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 276);

    auto ts_xxzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 277);

    auto ts_xxzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 278);

    auto ts_xxzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 279);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxz_xxxxxx, ts_xxz_xxxxxy, ts_xxz_xxxxyy, ts_xxz_xxxyyy, ts_xxz_xxyyyy, ts_xxz_xyyyyy, ts_xxzz_xxxxxx, ts_xxzz_xxxxxy, ts_xxzz_xxxxyy, ts_xxzz_xxxyyy, ts_xxzz_xxyyyy, ts_xxzz_xyyyyy, ts_xxzzz_xxxxxx, ts_xxzzz_xxxxxy, ts_xxzzz_xxxxxz, ts_xxzzz_xxxxyy, ts_xxzzz_xxxxyz, ts_xxzzz_xxxxzz, ts_xxzzz_xxxyyy, ts_xxzzz_xxxyyz, ts_xxzzz_xxxyzz, ts_xxzzz_xxxzzz, ts_xxzzz_xxyyyy, ts_xxzzz_xxyyyz, ts_xxzzz_xxyyzz, ts_xxzzz_xxyzzz, ts_xxzzz_xxzzzz, ts_xxzzz_xyyyyy, ts_xxzzz_xyyyyz, ts_xxzzz_xyyyzz, ts_xxzzz_xyyzzz, ts_xxzzz_xyzzzz, ts_xxzzz_xzzzzz, ts_xxzzz_yyyyyy, ts_xxzzz_yyyyyz, ts_xxzzz_yyyyzz, ts_xxzzz_yyyzzz, ts_xxzzz_yyzzzz, ts_xxzzz_yzzzzz, ts_xxzzz_zzzzzz, ts_xzzz_xxxxxz, ts_xzzz_xxxxyz, ts_xzzz_xxxxz, ts_xzzz_xxxxzz, ts_xzzz_xxxyyz, ts_xzzz_xxxyz, ts_xzzz_xxxyzz, ts_xzzz_xxxzz, ts_xzzz_xxxzzz, ts_xzzz_xxyyyz, ts_xzzz_xxyyz, ts_xzzz_xxyyzz, ts_xzzz_xxyzz, ts_xzzz_xxyzzz, ts_xzzz_xxzzz, ts_xzzz_xxzzzz, ts_xzzz_xyyyyz, ts_xzzz_xyyyz, ts_xzzz_xyyyzz, ts_xzzz_xyyzz, ts_xzzz_xyyzzz, ts_xzzz_xyzzz, ts_xzzz_xyzzzz, ts_xzzz_xzzzz, ts_xzzz_xzzzzz, ts_xzzz_yyyyyy, ts_xzzz_yyyyyz, ts_xzzz_yyyyz, ts_xzzz_yyyyzz, ts_xzzz_yyyzz, ts_xzzz_yyyzzz, ts_xzzz_yyzzz, ts_xzzz_yyzzzz, ts_xzzz_yzzzz, ts_xzzz_yzzzzz, ts_xzzz_zzzzz, ts_xzzz_zzzzzz, ts_zzz_xxxxxz, ts_zzz_xxxxyz, ts_zzz_xxxxzz, ts_zzz_xxxyyz, ts_zzz_xxxyzz, ts_zzz_xxxzzz, ts_zzz_xxyyyz, ts_zzz_xxyyzz, ts_zzz_xxyzzz, ts_zzz_xxzzzz, ts_zzz_xyyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyyy, ts_zzz_yyyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzz_xxxxxx[i] = 2.0 * ts_xxz_xxxxxx[i] * fe_0 + ts_xxzz_xxxxxx[i] * pa_z[i];

        ts_xxzzz_xxxxxy[i] = 2.0 * ts_xxz_xxxxxy[i] * fe_0 + ts_xxzz_xxxxxy[i] * pa_z[i];

        ts_xxzzz_xxxxxz[i] = ts_zzz_xxxxxz[i] * fe_0 + 5.0 * ts_xzzz_xxxxz[i] * fe_0 + ts_xzzz_xxxxxz[i] * pa_x[i];

        ts_xxzzz_xxxxyy[i] = 2.0 * ts_xxz_xxxxyy[i] * fe_0 + ts_xxzz_xxxxyy[i] * pa_z[i];

        ts_xxzzz_xxxxyz[i] = ts_zzz_xxxxyz[i] * fe_0 + 4.0 * ts_xzzz_xxxyz[i] * fe_0 + ts_xzzz_xxxxyz[i] * pa_x[i];

        ts_xxzzz_xxxxzz[i] = ts_zzz_xxxxzz[i] * fe_0 + 4.0 * ts_xzzz_xxxzz[i] * fe_0 + ts_xzzz_xxxxzz[i] * pa_x[i];

        ts_xxzzz_xxxyyy[i] = 2.0 * ts_xxz_xxxyyy[i] * fe_0 + ts_xxzz_xxxyyy[i] * pa_z[i];

        ts_xxzzz_xxxyyz[i] = ts_zzz_xxxyyz[i] * fe_0 + 3.0 * ts_xzzz_xxyyz[i] * fe_0 + ts_xzzz_xxxyyz[i] * pa_x[i];

        ts_xxzzz_xxxyzz[i] = ts_zzz_xxxyzz[i] * fe_0 + 3.0 * ts_xzzz_xxyzz[i] * fe_0 + ts_xzzz_xxxyzz[i] * pa_x[i];

        ts_xxzzz_xxxzzz[i] = ts_zzz_xxxzzz[i] * fe_0 + 3.0 * ts_xzzz_xxzzz[i] * fe_0 + ts_xzzz_xxxzzz[i] * pa_x[i];

        ts_xxzzz_xxyyyy[i] = 2.0 * ts_xxz_xxyyyy[i] * fe_0 + ts_xxzz_xxyyyy[i] * pa_z[i];

        ts_xxzzz_xxyyyz[i] = ts_zzz_xxyyyz[i] * fe_0 + 2.0 * ts_xzzz_xyyyz[i] * fe_0 + ts_xzzz_xxyyyz[i] * pa_x[i];

        ts_xxzzz_xxyyzz[i] = ts_zzz_xxyyzz[i] * fe_0 + 2.0 * ts_xzzz_xyyzz[i] * fe_0 + ts_xzzz_xxyyzz[i] * pa_x[i];

        ts_xxzzz_xxyzzz[i] = ts_zzz_xxyzzz[i] * fe_0 + 2.0 * ts_xzzz_xyzzz[i] * fe_0 + ts_xzzz_xxyzzz[i] * pa_x[i];

        ts_xxzzz_xxzzzz[i] = ts_zzz_xxzzzz[i] * fe_0 + 2.0 * ts_xzzz_xzzzz[i] * fe_0 + ts_xzzz_xxzzzz[i] * pa_x[i];

        ts_xxzzz_xyyyyy[i] = 2.0 * ts_xxz_xyyyyy[i] * fe_0 + ts_xxzz_xyyyyy[i] * pa_z[i];

        ts_xxzzz_xyyyyz[i] = ts_zzz_xyyyyz[i] * fe_0 + ts_xzzz_yyyyz[i] * fe_0 + ts_xzzz_xyyyyz[i] * pa_x[i];

        ts_xxzzz_xyyyzz[i] = ts_zzz_xyyyzz[i] * fe_0 + ts_xzzz_yyyzz[i] * fe_0 + ts_xzzz_xyyyzz[i] * pa_x[i];

        ts_xxzzz_xyyzzz[i] = ts_zzz_xyyzzz[i] * fe_0 + ts_xzzz_yyzzz[i] * fe_0 + ts_xzzz_xyyzzz[i] * pa_x[i];

        ts_xxzzz_xyzzzz[i] = ts_zzz_xyzzzz[i] * fe_0 + ts_xzzz_yzzzz[i] * fe_0 + ts_xzzz_xyzzzz[i] * pa_x[i];

        ts_xxzzz_xzzzzz[i] = ts_zzz_xzzzzz[i] * fe_0 + ts_xzzz_zzzzz[i] * fe_0 + ts_xzzz_xzzzzz[i] * pa_x[i];

        ts_xxzzz_yyyyyy[i] = ts_zzz_yyyyyy[i] * fe_0 + ts_xzzz_yyyyyy[i] * pa_x[i];

        ts_xxzzz_yyyyyz[i] = ts_zzz_yyyyyz[i] * fe_0 + ts_xzzz_yyyyyz[i] * pa_x[i];

        ts_xxzzz_yyyyzz[i] = ts_zzz_yyyyzz[i] * fe_0 + ts_xzzz_yyyyzz[i] * pa_x[i];

        ts_xxzzz_yyyzzz[i] = ts_zzz_yyyzzz[i] * fe_0 + ts_xzzz_yyyzzz[i] * pa_x[i];

        ts_xxzzz_yyzzzz[i] = ts_zzz_yyzzzz[i] * fe_0 + ts_xzzz_yyzzzz[i] * pa_x[i];

        ts_xxzzz_yzzzzz[i] = ts_zzz_yzzzzz[i] * fe_0 + ts_xzzz_yzzzzz[i] * pa_x[i];

        ts_xxzzz_zzzzzz[i] = ts_zzz_zzzzzz[i] * fe_0 + ts_xzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : HI

    auto ts_xyyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 280);

    auto ts_xyyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 281);

    auto ts_xyyyy_xxxxxz = pbuffer.data(idx_ovl_hi + 282);

    auto ts_xyyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 283);

    auto ts_xyyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 284);

    auto ts_xyyyy_xxxxzz = pbuffer.data(idx_ovl_hi + 285);

    auto ts_xyyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 286);

    auto ts_xyyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 287);

    auto ts_xyyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 288);

    auto ts_xyyyy_xxxzzz = pbuffer.data(idx_ovl_hi + 289);

    auto ts_xyyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 290);

    auto ts_xyyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 291);

    auto ts_xyyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 292);

    auto ts_xyyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 293);

    auto ts_xyyyy_xxzzzz = pbuffer.data(idx_ovl_hi + 294);

    auto ts_xyyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 295);

    auto ts_xyyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 296);

    auto ts_xyyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 297);

    auto ts_xyyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 298);

    auto ts_xyyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 299);

    auto ts_xyyyy_xzzzzz = pbuffer.data(idx_ovl_hi + 300);

    auto ts_xyyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 301);

    auto ts_xyyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 302);

    auto ts_xyyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 303);

    auto ts_xyyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 304);

    auto ts_xyyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 305);

    auto ts_xyyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 306);

    auto ts_xyyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 307);

    #pragma omp simd aligned(pa_x, ts_xyyyy_xxxxxx, ts_xyyyy_xxxxxy, ts_xyyyy_xxxxxz, ts_xyyyy_xxxxyy, ts_xyyyy_xxxxyz, ts_xyyyy_xxxxzz, ts_xyyyy_xxxyyy, ts_xyyyy_xxxyyz, ts_xyyyy_xxxyzz, ts_xyyyy_xxxzzz, ts_xyyyy_xxyyyy, ts_xyyyy_xxyyyz, ts_xyyyy_xxyyzz, ts_xyyyy_xxyzzz, ts_xyyyy_xxzzzz, ts_xyyyy_xyyyyy, ts_xyyyy_xyyyyz, ts_xyyyy_xyyyzz, ts_xyyyy_xyyzzz, ts_xyyyy_xyzzzz, ts_xyyyy_xzzzzz, ts_xyyyy_yyyyyy, ts_xyyyy_yyyyyz, ts_xyyyy_yyyyzz, ts_xyyyy_yyyzzz, ts_xyyyy_yyzzzz, ts_xyyyy_yzzzzz, ts_xyyyy_zzzzzz, ts_yyyy_xxxxx, ts_yyyy_xxxxxx, ts_yyyy_xxxxxy, ts_yyyy_xxxxxz, ts_yyyy_xxxxy, ts_yyyy_xxxxyy, ts_yyyy_xxxxyz, ts_yyyy_xxxxz, ts_yyyy_xxxxzz, ts_yyyy_xxxyy, ts_yyyy_xxxyyy, ts_yyyy_xxxyyz, ts_yyyy_xxxyz, ts_yyyy_xxxyzz, ts_yyyy_xxxzz, ts_yyyy_xxxzzz, ts_yyyy_xxyyy, ts_yyyy_xxyyyy, ts_yyyy_xxyyyz, ts_yyyy_xxyyz, ts_yyyy_xxyyzz, ts_yyyy_xxyzz, ts_yyyy_xxyzzz, ts_yyyy_xxzzz, ts_yyyy_xxzzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyyy, ts_yyyy_xyyyyz, ts_yyyy_xyyyz, ts_yyyy_xyyyzz, ts_yyyy_xyyzz, ts_yyyy_xyyzzz, ts_yyyy_xyzzz, ts_yyyy_xyzzzz, ts_yyyy_xzzzz, ts_yyyy_xzzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyyy, ts_yyyy_yyyyyz, ts_yyyy_yyyyz, ts_yyyy_yyyyzz, ts_yyyy_yyyzz, ts_yyyy_yyyzzz, ts_yyyy_yyzzz, ts_yyyy_yyzzzz, ts_yyyy_yzzzz, ts_yyyy_yzzzzz, ts_yyyy_zzzzz, ts_yyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyy_xxxxxx[i] = 6.0 * ts_yyyy_xxxxx[i] * fe_0 + ts_yyyy_xxxxxx[i] * pa_x[i];

        ts_xyyyy_xxxxxy[i] = 5.0 * ts_yyyy_xxxxy[i] * fe_0 + ts_yyyy_xxxxxy[i] * pa_x[i];

        ts_xyyyy_xxxxxz[i] = 5.0 * ts_yyyy_xxxxz[i] * fe_0 + ts_yyyy_xxxxxz[i] * pa_x[i];

        ts_xyyyy_xxxxyy[i] = 4.0 * ts_yyyy_xxxyy[i] * fe_0 + ts_yyyy_xxxxyy[i] * pa_x[i];

        ts_xyyyy_xxxxyz[i] = 4.0 * ts_yyyy_xxxyz[i] * fe_0 + ts_yyyy_xxxxyz[i] * pa_x[i];

        ts_xyyyy_xxxxzz[i] = 4.0 * ts_yyyy_xxxzz[i] * fe_0 + ts_yyyy_xxxxzz[i] * pa_x[i];

        ts_xyyyy_xxxyyy[i] = 3.0 * ts_yyyy_xxyyy[i] * fe_0 + ts_yyyy_xxxyyy[i] * pa_x[i];

        ts_xyyyy_xxxyyz[i] = 3.0 * ts_yyyy_xxyyz[i] * fe_0 + ts_yyyy_xxxyyz[i] * pa_x[i];

        ts_xyyyy_xxxyzz[i] = 3.0 * ts_yyyy_xxyzz[i] * fe_0 + ts_yyyy_xxxyzz[i] * pa_x[i];

        ts_xyyyy_xxxzzz[i] = 3.0 * ts_yyyy_xxzzz[i] * fe_0 + ts_yyyy_xxxzzz[i] * pa_x[i];

        ts_xyyyy_xxyyyy[i] = 2.0 * ts_yyyy_xyyyy[i] * fe_0 + ts_yyyy_xxyyyy[i] * pa_x[i];

        ts_xyyyy_xxyyyz[i] = 2.0 * ts_yyyy_xyyyz[i] * fe_0 + ts_yyyy_xxyyyz[i] * pa_x[i];

        ts_xyyyy_xxyyzz[i] = 2.0 * ts_yyyy_xyyzz[i] * fe_0 + ts_yyyy_xxyyzz[i] * pa_x[i];

        ts_xyyyy_xxyzzz[i] = 2.0 * ts_yyyy_xyzzz[i] * fe_0 + ts_yyyy_xxyzzz[i] * pa_x[i];

        ts_xyyyy_xxzzzz[i] = 2.0 * ts_yyyy_xzzzz[i] * fe_0 + ts_yyyy_xxzzzz[i] * pa_x[i];

        ts_xyyyy_xyyyyy[i] = ts_yyyy_yyyyy[i] * fe_0 + ts_yyyy_xyyyyy[i] * pa_x[i];

        ts_xyyyy_xyyyyz[i] = ts_yyyy_yyyyz[i] * fe_0 + ts_yyyy_xyyyyz[i] * pa_x[i];

        ts_xyyyy_xyyyzz[i] = ts_yyyy_yyyzz[i] * fe_0 + ts_yyyy_xyyyzz[i] * pa_x[i];

        ts_xyyyy_xyyzzz[i] = ts_yyyy_yyzzz[i] * fe_0 + ts_yyyy_xyyzzz[i] * pa_x[i];

        ts_xyyyy_xyzzzz[i] = ts_yyyy_yzzzz[i] * fe_0 + ts_yyyy_xyzzzz[i] * pa_x[i];

        ts_xyyyy_xzzzzz[i] = ts_yyyy_zzzzz[i] * fe_0 + ts_yyyy_xzzzzz[i] * pa_x[i];

        ts_xyyyy_yyyyyy[i] = ts_yyyy_yyyyyy[i] * pa_x[i];

        ts_xyyyy_yyyyyz[i] = ts_yyyy_yyyyyz[i] * pa_x[i];

        ts_xyyyy_yyyyzz[i] = ts_yyyy_yyyyzz[i] * pa_x[i];

        ts_xyyyy_yyyzzz[i] = ts_yyyy_yyyzzz[i] * pa_x[i];

        ts_xyyyy_yyzzzz[i] = ts_yyyy_yyzzzz[i] * pa_x[i];

        ts_xyyyy_yzzzzz[i] = ts_yyyy_yzzzzz[i] * pa_x[i];

        ts_xyyyy_zzzzzz[i] = ts_yyyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 308-336 components of targeted buffer : HI

    auto ts_xyyyz_xxxxxx = pbuffer.data(idx_ovl_hi + 308);

    auto ts_xyyyz_xxxxxy = pbuffer.data(idx_ovl_hi + 309);

    auto ts_xyyyz_xxxxxz = pbuffer.data(idx_ovl_hi + 310);

    auto ts_xyyyz_xxxxyy = pbuffer.data(idx_ovl_hi + 311);

    auto ts_xyyyz_xxxxyz = pbuffer.data(idx_ovl_hi + 312);

    auto ts_xyyyz_xxxxzz = pbuffer.data(idx_ovl_hi + 313);

    auto ts_xyyyz_xxxyyy = pbuffer.data(idx_ovl_hi + 314);

    auto ts_xyyyz_xxxyyz = pbuffer.data(idx_ovl_hi + 315);

    auto ts_xyyyz_xxxyzz = pbuffer.data(idx_ovl_hi + 316);

    auto ts_xyyyz_xxxzzz = pbuffer.data(idx_ovl_hi + 317);

    auto ts_xyyyz_xxyyyy = pbuffer.data(idx_ovl_hi + 318);

    auto ts_xyyyz_xxyyyz = pbuffer.data(idx_ovl_hi + 319);

    auto ts_xyyyz_xxyyzz = pbuffer.data(idx_ovl_hi + 320);

    auto ts_xyyyz_xxyzzz = pbuffer.data(idx_ovl_hi + 321);

    auto ts_xyyyz_xxzzzz = pbuffer.data(idx_ovl_hi + 322);

    auto ts_xyyyz_xyyyyy = pbuffer.data(idx_ovl_hi + 323);

    auto ts_xyyyz_xyyyyz = pbuffer.data(idx_ovl_hi + 324);

    auto ts_xyyyz_xyyyzz = pbuffer.data(idx_ovl_hi + 325);

    auto ts_xyyyz_xyyzzz = pbuffer.data(idx_ovl_hi + 326);

    auto ts_xyyyz_xyzzzz = pbuffer.data(idx_ovl_hi + 327);

    auto ts_xyyyz_xzzzzz = pbuffer.data(idx_ovl_hi + 328);

    auto ts_xyyyz_yyyyyy = pbuffer.data(idx_ovl_hi + 329);

    auto ts_xyyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 330);

    auto ts_xyyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 331);

    auto ts_xyyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 332);

    auto ts_xyyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 333);

    auto ts_xyyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 334);

    auto ts_xyyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 335);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyy_xxxxxx, ts_xyyy_xxxxxy, ts_xyyy_xxxxyy, ts_xyyy_xxxyyy, ts_xyyy_xxyyyy, ts_xyyy_xyyyyy, ts_xyyyz_xxxxxx, ts_xyyyz_xxxxxy, ts_xyyyz_xxxxxz, ts_xyyyz_xxxxyy, ts_xyyyz_xxxxyz, ts_xyyyz_xxxxzz, ts_xyyyz_xxxyyy, ts_xyyyz_xxxyyz, ts_xyyyz_xxxyzz, ts_xyyyz_xxxzzz, ts_xyyyz_xxyyyy, ts_xyyyz_xxyyyz, ts_xyyyz_xxyyzz, ts_xyyyz_xxyzzz, ts_xyyyz_xxzzzz, ts_xyyyz_xyyyyy, ts_xyyyz_xyyyyz, ts_xyyyz_xyyyzz, ts_xyyyz_xyyzzz, ts_xyyyz_xyzzzz, ts_xyyyz_xzzzzz, ts_xyyyz_yyyyyy, ts_xyyyz_yyyyyz, ts_xyyyz_yyyyzz, ts_xyyyz_yyyzzz, ts_xyyyz_yyzzzz, ts_xyyyz_yzzzzz, ts_xyyyz_zzzzzz, ts_yyyz_xxxxxz, ts_yyyz_xxxxyz, ts_yyyz_xxxxz, ts_yyyz_xxxxzz, ts_yyyz_xxxyyz, ts_yyyz_xxxyz, ts_yyyz_xxxyzz, ts_yyyz_xxxzz, ts_yyyz_xxxzzz, ts_yyyz_xxyyyz, ts_yyyz_xxyyz, ts_yyyz_xxyyzz, ts_yyyz_xxyzz, ts_yyyz_xxyzzz, ts_yyyz_xxzzz, ts_yyyz_xxzzzz, ts_yyyz_xyyyyz, ts_yyyz_xyyyz, ts_yyyz_xyyyzz, ts_yyyz_xyyzz, ts_yyyz_xyyzzz, ts_yyyz_xyzzz, ts_yyyz_xyzzzz, ts_yyyz_xzzzz, ts_yyyz_xzzzzz, ts_yyyz_yyyyyy, ts_yyyz_yyyyyz, ts_yyyz_yyyyz, ts_yyyz_yyyyzz, ts_yyyz_yyyzz, ts_yyyz_yyyzzz, ts_yyyz_yyzzz, ts_yyyz_yyzzzz, ts_yyyz_yzzzz, ts_yyyz_yzzzzz, ts_yyyz_zzzzz, ts_yyyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyz_xxxxxx[i] = ts_xyyy_xxxxxx[i] * pa_z[i];

        ts_xyyyz_xxxxxy[i] = ts_xyyy_xxxxxy[i] * pa_z[i];

        ts_xyyyz_xxxxxz[i] = 5.0 * ts_yyyz_xxxxz[i] * fe_0 + ts_yyyz_xxxxxz[i] * pa_x[i];

        ts_xyyyz_xxxxyy[i] = ts_xyyy_xxxxyy[i] * pa_z[i];

        ts_xyyyz_xxxxyz[i] = 4.0 * ts_yyyz_xxxyz[i] * fe_0 + ts_yyyz_xxxxyz[i] * pa_x[i];

        ts_xyyyz_xxxxzz[i] = 4.0 * ts_yyyz_xxxzz[i] * fe_0 + ts_yyyz_xxxxzz[i] * pa_x[i];

        ts_xyyyz_xxxyyy[i] = ts_xyyy_xxxyyy[i] * pa_z[i];

        ts_xyyyz_xxxyyz[i] = 3.0 * ts_yyyz_xxyyz[i] * fe_0 + ts_yyyz_xxxyyz[i] * pa_x[i];

        ts_xyyyz_xxxyzz[i] = 3.0 * ts_yyyz_xxyzz[i] * fe_0 + ts_yyyz_xxxyzz[i] * pa_x[i];

        ts_xyyyz_xxxzzz[i] = 3.0 * ts_yyyz_xxzzz[i] * fe_0 + ts_yyyz_xxxzzz[i] * pa_x[i];

        ts_xyyyz_xxyyyy[i] = ts_xyyy_xxyyyy[i] * pa_z[i];

        ts_xyyyz_xxyyyz[i] = 2.0 * ts_yyyz_xyyyz[i] * fe_0 + ts_yyyz_xxyyyz[i] * pa_x[i];

        ts_xyyyz_xxyyzz[i] = 2.0 * ts_yyyz_xyyzz[i] * fe_0 + ts_yyyz_xxyyzz[i] * pa_x[i];

        ts_xyyyz_xxyzzz[i] = 2.0 * ts_yyyz_xyzzz[i] * fe_0 + ts_yyyz_xxyzzz[i] * pa_x[i];

        ts_xyyyz_xxzzzz[i] = 2.0 * ts_yyyz_xzzzz[i] * fe_0 + ts_yyyz_xxzzzz[i] * pa_x[i];

        ts_xyyyz_xyyyyy[i] = ts_xyyy_xyyyyy[i] * pa_z[i];

        ts_xyyyz_xyyyyz[i] = ts_yyyz_yyyyz[i] * fe_0 + ts_yyyz_xyyyyz[i] * pa_x[i];

        ts_xyyyz_xyyyzz[i] = ts_yyyz_yyyzz[i] * fe_0 + ts_yyyz_xyyyzz[i] * pa_x[i];

        ts_xyyyz_xyyzzz[i] = ts_yyyz_yyzzz[i] * fe_0 + ts_yyyz_xyyzzz[i] * pa_x[i];

        ts_xyyyz_xyzzzz[i] = ts_yyyz_yzzzz[i] * fe_0 + ts_yyyz_xyzzzz[i] * pa_x[i];

        ts_xyyyz_xzzzzz[i] = ts_yyyz_zzzzz[i] * fe_0 + ts_yyyz_xzzzzz[i] * pa_x[i];

        ts_xyyyz_yyyyyy[i] = ts_yyyz_yyyyyy[i] * pa_x[i];

        ts_xyyyz_yyyyyz[i] = ts_yyyz_yyyyyz[i] * pa_x[i];

        ts_xyyyz_yyyyzz[i] = ts_yyyz_yyyyzz[i] * pa_x[i];

        ts_xyyyz_yyyzzz[i] = ts_yyyz_yyyzzz[i] * pa_x[i];

        ts_xyyyz_yyzzzz[i] = ts_yyyz_yyzzzz[i] * pa_x[i];

        ts_xyyyz_yzzzzz[i] = ts_yyyz_yzzzzz[i] * pa_x[i];

        ts_xyyyz_zzzzzz[i] = ts_yyyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 336-364 components of targeted buffer : HI

    auto ts_xyyzz_xxxxxx = pbuffer.data(idx_ovl_hi + 336);

    auto ts_xyyzz_xxxxxy = pbuffer.data(idx_ovl_hi + 337);

    auto ts_xyyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 338);

    auto ts_xyyzz_xxxxyy = pbuffer.data(idx_ovl_hi + 339);

    auto ts_xyyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 340);

    auto ts_xyyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 341);

    auto ts_xyyzz_xxxyyy = pbuffer.data(idx_ovl_hi + 342);

    auto ts_xyyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 343);

    auto ts_xyyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 344);

    auto ts_xyyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 345);

    auto ts_xyyzz_xxyyyy = pbuffer.data(idx_ovl_hi + 346);

    auto ts_xyyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 347);

    auto ts_xyyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 348);

    auto ts_xyyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 349);

    auto ts_xyyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 350);

    auto ts_xyyzz_xyyyyy = pbuffer.data(idx_ovl_hi + 351);

    auto ts_xyyzz_xyyyyz = pbuffer.data(idx_ovl_hi + 352);

    auto ts_xyyzz_xyyyzz = pbuffer.data(idx_ovl_hi + 353);

    auto ts_xyyzz_xyyzzz = pbuffer.data(idx_ovl_hi + 354);

    auto ts_xyyzz_xyzzzz = pbuffer.data(idx_ovl_hi + 355);

    auto ts_xyyzz_xzzzzz = pbuffer.data(idx_ovl_hi + 356);

    auto ts_xyyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 357);

    auto ts_xyyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 358);

    auto ts_xyyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 359);

    auto ts_xyyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 360);

    auto ts_xyyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 361);

    auto ts_xyyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 362);

    auto ts_xyyzz_zzzzzz = pbuffer.data(idx_ovl_hi + 363);

    #pragma omp simd aligned(pa_x, ts_xyyzz_xxxxxx, ts_xyyzz_xxxxxy, ts_xyyzz_xxxxxz, ts_xyyzz_xxxxyy, ts_xyyzz_xxxxyz, ts_xyyzz_xxxxzz, ts_xyyzz_xxxyyy, ts_xyyzz_xxxyyz, ts_xyyzz_xxxyzz, ts_xyyzz_xxxzzz, ts_xyyzz_xxyyyy, ts_xyyzz_xxyyyz, ts_xyyzz_xxyyzz, ts_xyyzz_xxyzzz, ts_xyyzz_xxzzzz, ts_xyyzz_xyyyyy, ts_xyyzz_xyyyyz, ts_xyyzz_xyyyzz, ts_xyyzz_xyyzzz, ts_xyyzz_xyzzzz, ts_xyyzz_xzzzzz, ts_xyyzz_yyyyyy, ts_xyyzz_yyyyyz, ts_xyyzz_yyyyzz, ts_xyyzz_yyyzzz, ts_xyyzz_yyzzzz, ts_xyyzz_yzzzzz, ts_xyyzz_zzzzzz, ts_yyzz_xxxxx, ts_yyzz_xxxxxx, ts_yyzz_xxxxxy, ts_yyzz_xxxxxz, ts_yyzz_xxxxy, ts_yyzz_xxxxyy, ts_yyzz_xxxxyz, ts_yyzz_xxxxz, ts_yyzz_xxxxzz, ts_yyzz_xxxyy, ts_yyzz_xxxyyy, ts_yyzz_xxxyyz, ts_yyzz_xxxyz, ts_yyzz_xxxyzz, ts_yyzz_xxxzz, ts_yyzz_xxxzzz, ts_yyzz_xxyyy, ts_yyzz_xxyyyy, ts_yyzz_xxyyyz, ts_yyzz_xxyyz, ts_yyzz_xxyyzz, ts_yyzz_xxyzz, ts_yyzz_xxyzzz, ts_yyzz_xxzzz, ts_yyzz_xxzzzz, ts_yyzz_xyyyy, ts_yyzz_xyyyyy, ts_yyzz_xyyyyz, ts_yyzz_xyyyz, ts_yyzz_xyyyzz, ts_yyzz_xyyzz, ts_yyzz_xyyzzz, ts_yyzz_xyzzz, ts_yyzz_xyzzzz, ts_yyzz_xzzzz, ts_yyzz_xzzzzz, ts_yyzz_yyyyy, ts_yyzz_yyyyyy, ts_yyzz_yyyyyz, ts_yyzz_yyyyz, ts_yyzz_yyyyzz, ts_yyzz_yyyzz, ts_yyzz_yyyzzz, ts_yyzz_yyzzz, ts_yyzz_yyzzzz, ts_yyzz_yzzzz, ts_yyzz_yzzzzz, ts_yyzz_zzzzz, ts_yyzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzz_xxxxxx[i] = 6.0 * ts_yyzz_xxxxx[i] * fe_0 + ts_yyzz_xxxxxx[i] * pa_x[i];

        ts_xyyzz_xxxxxy[i] = 5.0 * ts_yyzz_xxxxy[i] * fe_0 + ts_yyzz_xxxxxy[i] * pa_x[i];

        ts_xyyzz_xxxxxz[i] = 5.0 * ts_yyzz_xxxxz[i] * fe_0 + ts_yyzz_xxxxxz[i] * pa_x[i];

        ts_xyyzz_xxxxyy[i] = 4.0 * ts_yyzz_xxxyy[i] * fe_0 + ts_yyzz_xxxxyy[i] * pa_x[i];

        ts_xyyzz_xxxxyz[i] = 4.0 * ts_yyzz_xxxyz[i] * fe_0 + ts_yyzz_xxxxyz[i] * pa_x[i];

        ts_xyyzz_xxxxzz[i] = 4.0 * ts_yyzz_xxxzz[i] * fe_0 + ts_yyzz_xxxxzz[i] * pa_x[i];

        ts_xyyzz_xxxyyy[i] = 3.0 * ts_yyzz_xxyyy[i] * fe_0 + ts_yyzz_xxxyyy[i] * pa_x[i];

        ts_xyyzz_xxxyyz[i] = 3.0 * ts_yyzz_xxyyz[i] * fe_0 + ts_yyzz_xxxyyz[i] * pa_x[i];

        ts_xyyzz_xxxyzz[i] = 3.0 * ts_yyzz_xxyzz[i] * fe_0 + ts_yyzz_xxxyzz[i] * pa_x[i];

        ts_xyyzz_xxxzzz[i] = 3.0 * ts_yyzz_xxzzz[i] * fe_0 + ts_yyzz_xxxzzz[i] * pa_x[i];

        ts_xyyzz_xxyyyy[i] = 2.0 * ts_yyzz_xyyyy[i] * fe_0 + ts_yyzz_xxyyyy[i] * pa_x[i];

        ts_xyyzz_xxyyyz[i] = 2.0 * ts_yyzz_xyyyz[i] * fe_0 + ts_yyzz_xxyyyz[i] * pa_x[i];

        ts_xyyzz_xxyyzz[i] = 2.0 * ts_yyzz_xyyzz[i] * fe_0 + ts_yyzz_xxyyzz[i] * pa_x[i];

        ts_xyyzz_xxyzzz[i] = 2.0 * ts_yyzz_xyzzz[i] * fe_0 + ts_yyzz_xxyzzz[i] * pa_x[i];

        ts_xyyzz_xxzzzz[i] = 2.0 * ts_yyzz_xzzzz[i] * fe_0 + ts_yyzz_xxzzzz[i] * pa_x[i];

        ts_xyyzz_xyyyyy[i] = ts_yyzz_yyyyy[i] * fe_0 + ts_yyzz_xyyyyy[i] * pa_x[i];

        ts_xyyzz_xyyyyz[i] = ts_yyzz_yyyyz[i] * fe_0 + ts_yyzz_xyyyyz[i] * pa_x[i];

        ts_xyyzz_xyyyzz[i] = ts_yyzz_yyyzz[i] * fe_0 + ts_yyzz_xyyyzz[i] * pa_x[i];

        ts_xyyzz_xyyzzz[i] = ts_yyzz_yyzzz[i] * fe_0 + ts_yyzz_xyyzzz[i] * pa_x[i];

        ts_xyyzz_xyzzzz[i] = ts_yyzz_yzzzz[i] * fe_0 + ts_yyzz_xyzzzz[i] * pa_x[i];

        ts_xyyzz_xzzzzz[i] = ts_yyzz_zzzzz[i] * fe_0 + ts_yyzz_xzzzzz[i] * pa_x[i];

        ts_xyyzz_yyyyyy[i] = ts_yyzz_yyyyyy[i] * pa_x[i];

        ts_xyyzz_yyyyyz[i] = ts_yyzz_yyyyyz[i] * pa_x[i];

        ts_xyyzz_yyyyzz[i] = ts_yyzz_yyyyzz[i] * pa_x[i];

        ts_xyyzz_yyyzzz[i] = ts_yyzz_yyyzzz[i] * pa_x[i];

        ts_xyyzz_yyzzzz[i] = ts_yyzz_yyzzzz[i] * pa_x[i];

        ts_xyyzz_yzzzzz[i] = ts_yyzz_yzzzzz[i] * pa_x[i];

        ts_xyyzz_zzzzzz[i] = ts_yyzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 364-392 components of targeted buffer : HI

    auto ts_xyzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 364);

    auto ts_xyzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 365);

    auto ts_xyzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 366);

    auto ts_xyzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 367);

    auto ts_xyzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 368);

    auto ts_xyzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 369);

    auto ts_xyzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 370);

    auto ts_xyzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 371);

    auto ts_xyzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 372);

    auto ts_xyzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 373);

    auto ts_xyzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 374);

    auto ts_xyzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 375);

    auto ts_xyzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 376);

    auto ts_xyzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 377);

    auto ts_xyzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 378);

    auto ts_xyzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 379);

    auto ts_xyzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 380);

    auto ts_xyzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 381);

    auto ts_xyzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 382);

    auto ts_xyzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 383);

    auto ts_xyzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 384);

    auto ts_xyzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 385);

    auto ts_xyzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 386);

    auto ts_xyzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 387);

    auto ts_xyzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 388);

    auto ts_xyzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 389);

    auto ts_xyzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 390);

    auto ts_xyzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 391);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzz_xxxxxx, ts_xyzzz_xxxxxy, ts_xyzzz_xxxxxz, ts_xyzzz_xxxxyy, ts_xyzzz_xxxxyz, ts_xyzzz_xxxxzz, ts_xyzzz_xxxyyy, ts_xyzzz_xxxyyz, ts_xyzzz_xxxyzz, ts_xyzzz_xxxzzz, ts_xyzzz_xxyyyy, ts_xyzzz_xxyyyz, ts_xyzzz_xxyyzz, ts_xyzzz_xxyzzz, ts_xyzzz_xxzzzz, ts_xyzzz_xyyyyy, ts_xyzzz_xyyyyz, ts_xyzzz_xyyyzz, ts_xyzzz_xyyzzz, ts_xyzzz_xyzzzz, ts_xyzzz_xzzzzz, ts_xyzzz_yyyyyy, ts_xyzzz_yyyyyz, ts_xyzzz_yyyyzz, ts_xyzzz_yyyzzz, ts_xyzzz_yyzzzz, ts_xyzzz_yzzzzz, ts_xyzzz_zzzzzz, ts_xzzz_xxxxxx, ts_xzzz_xxxxxz, ts_xzzz_xxxxzz, ts_xzzz_xxxzzz, ts_xzzz_xxzzzz, ts_xzzz_xzzzzz, ts_yzzz_xxxxxy, ts_yzzz_xxxxy, ts_yzzz_xxxxyy, ts_yzzz_xxxxyz, ts_yzzz_xxxyy, ts_yzzz_xxxyyy, ts_yzzz_xxxyyz, ts_yzzz_xxxyz, ts_yzzz_xxxyzz, ts_yzzz_xxyyy, ts_yzzz_xxyyyy, ts_yzzz_xxyyyz, ts_yzzz_xxyyz, ts_yzzz_xxyyzz, ts_yzzz_xxyzz, ts_yzzz_xxyzzz, ts_yzzz_xyyyy, ts_yzzz_xyyyyy, ts_yzzz_xyyyyz, ts_yzzz_xyyyz, ts_yzzz_xyyyzz, ts_yzzz_xyyzz, ts_yzzz_xyyzzz, ts_yzzz_xyzzz, ts_yzzz_xyzzzz, ts_yzzz_yyyyy, ts_yzzz_yyyyyy, ts_yzzz_yyyyyz, ts_yzzz_yyyyz, ts_yzzz_yyyyzz, ts_yzzz_yyyzz, ts_yzzz_yyyzzz, ts_yzzz_yyzzz, ts_yzzz_yyzzzz, ts_yzzz_yzzzz, ts_yzzz_yzzzzz, ts_yzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzz_xxxxxx[i] = ts_xzzz_xxxxxx[i] * pa_y[i];

        ts_xyzzz_xxxxxy[i] = 5.0 * ts_yzzz_xxxxy[i] * fe_0 + ts_yzzz_xxxxxy[i] * pa_x[i];

        ts_xyzzz_xxxxxz[i] = ts_xzzz_xxxxxz[i] * pa_y[i];

        ts_xyzzz_xxxxyy[i] = 4.0 * ts_yzzz_xxxyy[i] * fe_0 + ts_yzzz_xxxxyy[i] * pa_x[i];

        ts_xyzzz_xxxxyz[i] = 4.0 * ts_yzzz_xxxyz[i] * fe_0 + ts_yzzz_xxxxyz[i] * pa_x[i];

        ts_xyzzz_xxxxzz[i] = ts_xzzz_xxxxzz[i] * pa_y[i];

        ts_xyzzz_xxxyyy[i] = 3.0 * ts_yzzz_xxyyy[i] * fe_0 + ts_yzzz_xxxyyy[i] * pa_x[i];

        ts_xyzzz_xxxyyz[i] = 3.0 * ts_yzzz_xxyyz[i] * fe_0 + ts_yzzz_xxxyyz[i] * pa_x[i];

        ts_xyzzz_xxxyzz[i] = 3.0 * ts_yzzz_xxyzz[i] * fe_0 + ts_yzzz_xxxyzz[i] * pa_x[i];

        ts_xyzzz_xxxzzz[i] = ts_xzzz_xxxzzz[i] * pa_y[i];

        ts_xyzzz_xxyyyy[i] = 2.0 * ts_yzzz_xyyyy[i] * fe_0 + ts_yzzz_xxyyyy[i] * pa_x[i];

        ts_xyzzz_xxyyyz[i] = 2.0 * ts_yzzz_xyyyz[i] * fe_0 + ts_yzzz_xxyyyz[i] * pa_x[i];

        ts_xyzzz_xxyyzz[i] = 2.0 * ts_yzzz_xyyzz[i] * fe_0 + ts_yzzz_xxyyzz[i] * pa_x[i];

        ts_xyzzz_xxyzzz[i] = 2.0 * ts_yzzz_xyzzz[i] * fe_0 + ts_yzzz_xxyzzz[i] * pa_x[i];

        ts_xyzzz_xxzzzz[i] = ts_xzzz_xxzzzz[i] * pa_y[i];

        ts_xyzzz_xyyyyy[i] = ts_yzzz_yyyyy[i] * fe_0 + ts_yzzz_xyyyyy[i] * pa_x[i];

        ts_xyzzz_xyyyyz[i] = ts_yzzz_yyyyz[i] * fe_0 + ts_yzzz_xyyyyz[i] * pa_x[i];

        ts_xyzzz_xyyyzz[i] = ts_yzzz_yyyzz[i] * fe_0 + ts_yzzz_xyyyzz[i] * pa_x[i];

        ts_xyzzz_xyyzzz[i] = ts_yzzz_yyzzz[i] * fe_0 + ts_yzzz_xyyzzz[i] * pa_x[i];

        ts_xyzzz_xyzzzz[i] = ts_yzzz_yzzzz[i] * fe_0 + ts_yzzz_xyzzzz[i] * pa_x[i];

        ts_xyzzz_xzzzzz[i] = ts_xzzz_xzzzzz[i] * pa_y[i];

        ts_xyzzz_yyyyyy[i] = ts_yzzz_yyyyyy[i] * pa_x[i];

        ts_xyzzz_yyyyyz[i] = ts_yzzz_yyyyyz[i] * pa_x[i];

        ts_xyzzz_yyyyzz[i] = ts_yzzz_yyyyzz[i] * pa_x[i];

        ts_xyzzz_yyyzzz[i] = ts_yzzz_yyyzzz[i] * pa_x[i];

        ts_xyzzz_yyzzzz[i] = ts_yzzz_yyzzzz[i] * pa_x[i];

        ts_xyzzz_yzzzzz[i] = ts_yzzz_yzzzzz[i] * pa_x[i];

        ts_xyzzz_zzzzzz[i] = ts_yzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 392-420 components of targeted buffer : HI

    auto ts_xzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 392);

    auto ts_xzzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 393);

    auto ts_xzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 394);

    auto ts_xzzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 395);

    auto ts_xzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 396);

    auto ts_xzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 397);

    auto ts_xzzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 398);

    auto ts_xzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 399);

    auto ts_xzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 400);

    auto ts_xzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 401);

    auto ts_xzzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 402);

    auto ts_xzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 403);

    auto ts_xzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 404);

    auto ts_xzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 405);

    auto ts_xzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 406);

    auto ts_xzzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 407);

    auto ts_xzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 408);

    auto ts_xzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 409);

    auto ts_xzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 410);

    auto ts_xzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 411);

    auto ts_xzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 412);

    auto ts_xzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 413);

    auto ts_xzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 414);

    auto ts_xzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 415);

    auto ts_xzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 416);

    auto ts_xzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 417);

    auto ts_xzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 418);

    auto ts_xzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 419);

    #pragma omp simd aligned(pa_x, ts_xzzzz_xxxxxx, ts_xzzzz_xxxxxy, ts_xzzzz_xxxxxz, ts_xzzzz_xxxxyy, ts_xzzzz_xxxxyz, ts_xzzzz_xxxxzz, ts_xzzzz_xxxyyy, ts_xzzzz_xxxyyz, ts_xzzzz_xxxyzz, ts_xzzzz_xxxzzz, ts_xzzzz_xxyyyy, ts_xzzzz_xxyyyz, ts_xzzzz_xxyyzz, ts_xzzzz_xxyzzz, ts_xzzzz_xxzzzz, ts_xzzzz_xyyyyy, ts_xzzzz_xyyyyz, ts_xzzzz_xyyyzz, ts_xzzzz_xyyzzz, ts_xzzzz_xyzzzz, ts_xzzzz_xzzzzz, ts_xzzzz_yyyyyy, ts_xzzzz_yyyyyz, ts_xzzzz_yyyyzz, ts_xzzzz_yyyzzz, ts_xzzzz_yyzzzz, ts_xzzzz_yzzzzz, ts_xzzzz_zzzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxxx, ts_zzzz_xxxxxy, ts_zzzz_xxxxxz, ts_zzzz_xxxxy, ts_zzzz_xxxxyy, ts_zzzz_xxxxyz, ts_zzzz_xxxxz, ts_zzzz_xxxxzz, ts_zzzz_xxxyy, ts_zzzz_xxxyyy, ts_zzzz_xxxyyz, ts_zzzz_xxxyz, ts_zzzz_xxxyzz, ts_zzzz_xxxzz, ts_zzzz_xxxzzz, ts_zzzz_xxyyy, ts_zzzz_xxyyyy, ts_zzzz_xxyyyz, ts_zzzz_xxyyz, ts_zzzz_xxyyzz, ts_zzzz_xxyzz, ts_zzzz_xxyzzz, ts_zzzz_xxzzz, ts_zzzz_xxzzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyyy, ts_zzzz_xyyyyz, ts_zzzz_xyyyz, ts_zzzz_xyyyzz, ts_zzzz_xyyzz, ts_zzzz_xyyzzz, ts_zzzz_xyzzz, ts_zzzz_xyzzzz, ts_zzzz_xzzzz, ts_zzzz_xzzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyyy, ts_zzzz_yyyyyz, ts_zzzz_yyyyz, ts_zzzz_yyyyzz, ts_zzzz_yyyzz, ts_zzzz_yyyzzz, ts_zzzz_yyzzz, ts_zzzz_yyzzzz, ts_zzzz_yzzzz, ts_zzzz_yzzzzz, ts_zzzz_zzzzz, ts_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzz_xxxxxx[i] = 6.0 * ts_zzzz_xxxxx[i] * fe_0 + ts_zzzz_xxxxxx[i] * pa_x[i];

        ts_xzzzz_xxxxxy[i] = 5.0 * ts_zzzz_xxxxy[i] * fe_0 + ts_zzzz_xxxxxy[i] * pa_x[i];

        ts_xzzzz_xxxxxz[i] = 5.0 * ts_zzzz_xxxxz[i] * fe_0 + ts_zzzz_xxxxxz[i] * pa_x[i];

        ts_xzzzz_xxxxyy[i] = 4.0 * ts_zzzz_xxxyy[i] * fe_0 + ts_zzzz_xxxxyy[i] * pa_x[i];

        ts_xzzzz_xxxxyz[i] = 4.0 * ts_zzzz_xxxyz[i] * fe_0 + ts_zzzz_xxxxyz[i] * pa_x[i];

        ts_xzzzz_xxxxzz[i] = 4.0 * ts_zzzz_xxxzz[i] * fe_0 + ts_zzzz_xxxxzz[i] * pa_x[i];

        ts_xzzzz_xxxyyy[i] = 3.0 * ts_zzzz_xxyyy[i] * fe_0 + ts_zzzz_xxxyyy[i] * pa_x[i];

        ts_xzzzz_xxxyyz[i] = 3.0 * ts_zzzz_xxyyz[i] * fe_0 + ts_zzzz_xxxyyz[i] * pa_x[i];

        ts_xzzzz_xxxyzz[i] = 3.0 * ts_zzzz_xxyzz[i] * fe_0 + ts_zzzz_xxxyzz[i] * pa_x[i];

        ts_xzzzz_xxxzzz[i] = 3.0 * ts_zzzz_xxzzz[i] * fe_0 + ts_zzzz_xxxzzz[i] * pa_x[i];

        ts_xzzzz_xxyyyy[i] = 2.0 * ts_zzzz_xyyyy[i] * fe_0 + ts_zzzz_xxyyyy[i] * pa_x[i];

        ts_xzzzz_xxyyyz[i] = 2.0 * ts_zzzz_xyyyz[i] * fe_0 + ts_zzzz_xxyyyz[i] * pa_x[i];

        ts_xzzzz_xxyyzz[i] = 2.0 * ts_zzzz_xyyzz[i] * fe_0 + ts_zzzz_xxyyzz[i] * pa_x[i];

        ts_xzzzz_xxyzzz[i] = 2.0 * ts_zzzz_xyzzz[i] * fe_0 + ts_zzzz_xxyzzz[i] * pa_x[i];

        ts_xzzzz_xxzzzz[i] = 2.0 * ts_zzzz_xzzzz[i] * fe_0 + ts_zzzz_xxzzzz[i] * pa_x[i];

        ts_xzzzz_xyyyyy[i] = ts_zzzz_yyyyy[i] * fe_0 + ts_zzzz_xyyyyy[i] * pa_x[i];

        ts_xzzzz_xyyyyz[i] = ts_zzzz_yyyyz[i] * fe_0 + ts_zzzz_xyyyyz[i] * pa_x[i];

        ts_xzzzz_xyyyzz[i] = ts_zzzz_yyyzz[i] * fe_0 + ts_zzzz_xyyyzz[i] * pa_x[i];

        ts_xzzzz_xyyzzz[i] = ts_zzzz_yyzzz[i] * fe_0 + ts_zzzz_xyyzzz[i] * pa_x[i];

        ts_xzzzz_xyzzzz[i] = ts_zzzz_yzzzz[i] * fe_0 + ts_zzzz_xyzzzz[i] * pa_x[i];

        ts_xzzzz_xzzzzz[i] = ts_zzzz_zzzzz[i] * fe_0 + ts_zzzz_xzzzzz[i] * pa_x[i];

        ts_xzzzz_yyyyyy[i] = ts_zzzz_yyyyyy[i] * pa_x[i];

        ts_xzzzz_yyyyyz[i] = ts_zzzz_yyyyyz[i] * pa_x[i];

        ts_xzzzz_yyyyzz[i] = ts_zzzz_yyyyzz[i] * pa_x[i];

        ts_xzzzz_yyyzzz[i] = ts_zzzz_yyyzzz[i] * pa_x[i];

        ts_xzzzz_yyzzzz[i] = ts_zzzz_yyzzzz[i] * pa_x[i];

        ts_xzzzz_yzzzzz[i] = ts_zzzz_yzzzzz[i] * pa_x[i];

        ts_xzzzz_zzzzzz[i] = ts_zzzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : HI

    auto ts_yyyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 420);

    auto ts_yyyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 421);

    auto ts_yyyyy_xxxxxz = pbuffer.data(idx_ovl_hi + 422);

    auto ts_yyyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 423);

    auto ts_yyyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 424);

    auto ts_yyyyy_xxxxzz = pbuffer.data(idx_ovl_hi + 425);

    auto ts_yyyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 426);

    auto ts_yyyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 427);

    auto ts_yyyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 428);

    auto ts_yyyyy_xxxzzz = pbuffer.data(idx_ovl_hi + 429);

    auto ts_yyyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 430);

    auto ts_yyyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 431);

    auto ts_yyyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 432);

    auto ts_yyyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 433);

    auto ts_yyyyy_xxzzzz = pbuffer.data(idx_ovl_hi + 434);

    auto ts_yyyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 435);

    auto ts_yyyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 436);

    auto ts_yyyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 437);

    auto ts_yyyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 438);

    auto ts_yyyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 439);

    auto ts_yyyyy_xzzzzz = pbuffer.data(idx_ovl_hi + 440);

    auto ts_yyyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 441);

    auto ts_yyyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 442);

    auto ts_yyyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 443);

    auto ts_yyyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 444);

    auto ts_yyyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 445);

    auto ts_yyyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 446);

    auto ts_yyyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 447);

    #pragma omp simd aligned(pa_y, ts_yyy_xxxxxx, ts_yyy_xxxxxy, ts_yyy_xxxxxz, ts_yyy_xxxxyy, ts_yyy_xxxxyz, ts_yyy_xxxxzz, ts_yyy_xxxyyy, ts_yyy_xxxyyz, ts_yyy_xxxyzz, ts_yyy_xxxzzz, ts_yyy_xxyyyy, ts_yyy_xxyyyz, ts_yyy_xxyyzz, ts_yyy_xxyzzz, ts_yyy_xxzzzz, ts_yyy_xyyyyy, ts_yyy_xyyyyz, ts_yyy_xyyyzz, ts_yyy_xyyzzz, ts_yyy_xyzzzz, ts_yyy_xzzzzz, ts_yyy_yyyyyy, ts_yyy_yyyyyz, ts_yyy_yyyyzz, ts_yyy_yyyzzz, ts_yyy_yyzzzz, ts_yyy_yzzzzz, ts_yyy_zzzzzz, ts_yyyy_xxxxx, ts_yyyy_xxxxxx, ts_yyyy_xxxxxy, ts_yyyy_xxxxxz, ts_yyyy_xxxxy, ts_yyyy_xxxxyy, ts_yyyy_xxxxyz, ts_yyyy_xxxxz, ts_yyyy_xxxxzz, ts_yyyy_xxxyy, ts_yyyy_xxxyyy, ts_yyyy_xxxyyz, ts_yyyy_xxxyz, ts_yyyy_xxxyzz, ts_yyyy_xxxzz, ts_yyyy_xxxzzz, ts_yyyy_xxyyy, ts_yyyy_xxyyyy, ts_yyyy_xxyyyz, ts_yyyy_xxyyz, ts_yyyy_xxyyzz, ts_yyyy_xxyzz, ts_yyyy_xxyzzz, ts_yyyy_xxzzz, ts_yyyy_xxzzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyyy, ts_yyyy_xyyyyz, ts_yyyy_xyyyz, ts_yyyy_xyyyzz, ts_yyyy_xyyzz, ts_yyyy_xyyzzz, ts_yyyy_xyzzz, ts_yyyy_xyzzzz, ts_yyyy_xzzzz, ts_yyyy_xzzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyyy, ts_yyyy_yyyyyz, ts_yyyy_yyyyz, ts_yyyy_yyyyzz, ts_yyyy_yyyzz, ts_yyyy_yyyzzz, ts_yyyy_yyzzz, ts_yyyy_yyzzzz, ts_yyyy_yzzzz, ts_yyyy_yzzzzz, ts_yyyy_zzzzz, ts_yyyy_zzzzzz, ts_yyyyy_xxxxxx, ts_yyyyy_xxxxxy, ts_yyyyy_xxxxxz, ts_yyyyy_xxxxyy, ts_yyyyy_xxxxyz, ts_yyyyy_xxxxzz, ts_yyyyy_xxxyyy, ts_yyyyy_xxxyyz, ts_yyyyy_xxxyzz, ts_yyyyy_xxxzzz, ts_yyyyy_xxyyyy, ts_yyyyy_xxyyyz, ts_yyyyy_xxyyzz, ts_yyyyy_xxyzzz, ts_yyyyy_xxzzzz, ts_yyyyy_xyyyyy, ts_yyyyy_xyyyyz, ts_yyyyy_xyyyzz, ts_yyyyy_xyyzzz, ts_yyyyy_xyzzzz, ts_yyyyy_xzzzzz, ts_yyyyy_yyyyyy, ts_yyyyy_yyyyyz, ts_yyyyy_yyyyzz, ts_yyyyy_yyyzzz, ts_yyyyy_yyzzzz, ts_yyyyy_yzzzzz, ts_yyyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyy_xxxxxx[i] = 4.0 * ts_yyy_xxxxxx[i] * fe_0 + ts_yyyy_xxxxxx[i] * pa_y[i];

        ts_yyyyy_xxxxxy[i] = 4.0 * ts_yyy_xxxxxy[i] * fe_0 + ts_yyyy_xxxxx[i] * fe_0 + ts_yyyy_xxxxxy[i] * pa_y[i];

        ts_yyyyy_xxxxxz[i] = 4.0 * ts_yyy_xxxxxz[i] * fe_0 + ts_yyyy_xxxxxz[i] * pa_y[i];

        ts_yyyyy_xxxxyy[i] = 4.0 * ts_yyy_xxxxyy[i] * fe_0 + 2.0 * ts_yyyy_xxxxy[i] * fe_0 + ts_yyyy_xxxxyy[i] * pa_y[i];

        ts_yyyyy_xxxxyz[i] = 4.0 * ts_yyy_xxxxyz[i] * fe_0 + ts_yyyy_xxxxz[i] * fe_0 + ts_yyyy_xxxxyz[i] * pa_y[i];

        ts_yyyyy_xxxxzz[i] = 4.0 * ts_yyy_xxxxzz[i] * fe_0 + ts_yyyy_xxxxzz[i] * pa_y[i];

        ts_yyyyy_xxxyyy[i] = 4.0 * ts_yyy_xxxyyy[i] * fe_0 + 3.0 * ts_yyyy_xxxyy[i] * fe_0 + ts_yyyy_xxxyyy[i] * pa_y[i];

        ts_yyyyy_xxxyyz[i] = 4.0 * ts_yyy_xxxyyz[i] * fe_0 + 2.0 * ts_yyyy_xxxyz[i] * fe_0 + ts_yyyy_xxxyyz[i] * pa_y[i];

        ts_yyyyy_xxxyzz[i] = 4.0 * ts_yyy_xxxyzz[i] * fe_0 + ts_yyyy_xxxzz[i] * fe_0 + ts_yyyy_xxxyzz[i] * pa_y[i];

        ts_yyyyy_xxxzzz[i] = 4.0 * ts_yyy_xxxzzz[i] * fe_0 + ts_yyyy_xxxzzz[i] * pa_y[i];

        ts_yyyyy_xxyyyy[i] = 4.0 * ts_yyy_xxyyyy[i] * fe_0 + 4.0 * ts_yyyy_xxyyy[i] * fe_0 + ts_yyyy_xxyyyy[i] * pa_y[i];

        ts_yyyyy_xxyyyz[i] = 4.0 * ts_yyy_xxyyyz[i] * fe_0 + 3.0 * ts_yyyy_xxyyz[i] * fe_0 + ts_yyyy_xxyyyz[i] * pa_y[i];

        ts_yyyyy_xxyyzz[i] = 4.0 * ts_yyy_xxyyzz[i] * fe_0 + 2.0 * ts_yyyy_xxyzz[i] * fe_0 + ts_yyyy_xxyyzz[i] * pa_y[i];

        ts_yyyyy_xxyzzz[i] = 4.0 * ts_yyy_xxyzzz[i] * fe_0 + ts_yyyy_xxzzz[i] * fe_0 + ts_yyyy_xxyzzz[i] * pa_y[i];

        ts_yyyyy_xxzzzz[i] = 4.0 * ts_yyy_xxzzzz[i] * fe_0 + ts_yyyy_xxzzzz[i] * pa_y[i];

        ts_yyyyy_xyyyyy[i] = 4.0 * ts_yyy_xyyyyy[i] * fe_0 + 5.0 * ts_yyyy_xyyyy[i] * fe_0 + ts_yyyy_xyyyyy[i] * pa_y[i];

        ts_yyyyy_xyyyyz[i] = 4.0 * ts_yyy_xyyyyz[i] * fe_0 + 4.0 * ts_yyyy_xyyyz[i] * fe_0 + ts_yyyy_xyyyyz[i] * pa_y[i];

        ts_yyyyy_xyyyzz[i] = 4.0 * ts_yyy_xyyyzz[i] * fe_0 + 3.0 * ts_yyyy_xyyzz[i] * fe_0 + ts_yyyy_xyyyzz[i] * pa_y[i];

        ts_yyyyy_xyyzzz[i] = 4.0 * ts_yyy_xyyzzz[i] * fe_0 + 2.0 * ts_yyyy_xyzzz[i] * fe_0 + ts_yyyy_xyyzzz[i] * pa_y[i];

        ts_yyyyy_xyzzzz[i] = 4.0 * ts_yyy_xyzzzz[i] * fe_0 + ts_yyyy_xzzzz[i] * fe_0 + ts_yyyy_xyzzzz[i] * pa_y[i];

        ts_yyyyy_xzzzzz[i] = 4.0 * ts_yyy_xzzzzz[i] * fe_0 + ts_yyyy_xzzzzz[i] * pa_y[i];

        ts_yyyyy_yyyyyy[i] = 4.0 * ts_yyy_yyyyyy[i] * fe_0 + 6.0 * ts_yyyy_yyyyy[i] * fe_0 + ts_yyyy_yyyyyy[i] * pa_y[i];

        ts_yyyyy_yyyyyz[i] = 4.0 * ts_yyy_yyyyyz[i] * fe_0 + 5.0 * ts_yyyy_yyyyz[i] * fe_0 + ts_yyyy_yyyyyz[i] * pa_y[i];

        ts_yyyyy_yyyyzz[i] = 4.0 * ts_yyy_yyyyzz[i] * fe_0 + 4.0 * ts_yyyy_yyyzz[i] * fe_0 + ts_yyyy_yyyyzz[i] * pa_y[i];

        ts_yyyyy_yyyzzz[i] = 4.0 * ts_yyy_yyyzzz[i] * fe_0 + 3.0 * ts_yyyy_yyzzz[i] * fe_0 + ts_yyyy_yyyzzz[i] * pa_y[i];

        ts_yyyyy_yyzzzz[i] = 4.0 * ts_yyy_yyzzzz[i] * fe_0 + 2.0 * ts_yyyy_yzzzz[i] * fe_0 + ts_yyyy_yyzzzz[i] * pa_y[i];

        ts_yyyyy_yzzzzz[i] = 4.0 * ts_yyy_yzzzzz[i] * fe_0 + ts_yyyy_zzzzz[i] * fe_0 + ts_yyyy_yzzzzz[i] * pa_y[i];

        ts_yyyyy_zzzzzz[i] = 4.0 * ts_yyy_zzzzzz[i] * fe_0 + ts_yyyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 448-476 components of targeted buffer : HI

    auto ts_yyyyz_xxxxxx = pbuffer.data(idx_ovl_hi + 448);

    auto ts_yyyyz_xxxxxy = pbuffer.data(idx_ovl_hi + 449);

    auto ts_yyyyz_xxxxxz = pbuffer.data(idx_ovl_hi + 450);

    auto ts_yyyyz_xxxxyy = pbuffer.data(idx_ovl_hi + 451);

    auto ts_yyyyz_xxxxyz = pbuffer.data(idx_ovl_hi + 452);

    auto ts_yyyyz_xxxxzz = pbuffer.data(idx_ovl_hi + 453);

    auto ts_yyyyz_xxxyyy = pbuffer.data(idx_ovl_hi + 454);

    auto ts_yyyyz_xxxyyz = pbuffer.data(idx_ovl_hi + 455);

    auto ts_yyyyz_xxxyzz = pbuffer.data(idx_ovl_hi + 456);

    auto ts_yyyyz_xxxzzz = pbuffer.data(idx_ovl_hi + 457);

    auto ts_yyyyz_xxyyyy = pbuffer.data(idx_ovl_hi + 458);

    auto ts_yyyyz_xxyyyz = pbuffer.data(idx_ovl_hi + 459);

    auto ts_yyyyz_xxyyzz = pbuffer.data(idx_ovl_hi + 460);

    auto ts_yyyyz_xxyzzz = pbuffer.data(idx_ovl_hi + 461);

    auto ts_yyyyz_xxzzzz = pbuffer.data(idx_ovl_hi + 462);

    auto ts_yyyyz_xyyyyy = pbuffer.data(idx_ovl_hi + 463);

    auto ts_yyyyz_xyyyyz = pbuffer.data(idx_ovl_hi + 464);

    auto ts_yyyyz_xyyyzz = pbuffer.data(idx_ovl_hi + 465);

    auto ts_yyyyz_xyyzzz = pbuffer.data(idx_ovl_hi + 466);

    auto ts_yyyyz_xyzzzz = pbuffer.data(idx_ovl_hi + 467);

    auto ts_yyyyz_xzzzzz = pbuffer.data(idx_ovl_hi + 468);

    auto ts_yyyyz_yyyyyy = pbuffer.data(idx_ovl_hi + 469);

    auto ts_yyyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 470);

    auto ts_yyyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 471);

    auto ts_yyyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 472);

    auto ts_yyyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 473);

    auto ts_yyyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 474);

    auto ts_yyyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 475);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_xxxxxx, ts_yyyy_xxxxxy, ts_yyyy_xxxxy, ts_yyyy_xxxxyy, ts_yyyy_xxxxyz, ts_yyyy_xxxyy, ts_yyyy_xxxyyy, ts_yyyy_xxxyyz, ts_yyyy_xxxyz, ts_yyyy_xxxyzz, ts_yyyy_xxyyy, ts_yyyy_xxyyyy, ts_yyyy_xxyyyz, ts_yyyy_xxyyz, ts_yyyy_xxyyzz, ts_yyyy_xxyzz, ts_yyyy_xxyzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyyy, ts_yyyy_xyyyyz, ts_yyyy_xyyyz, ts_yyyy_xyyyzz, ts_yyyy_xyyzz, ts_yyyy_xyyzzz, ts_yyyy_xyzzz, ts_yyyy_xyzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyyy, ts_yyyy_yyyyyz, ts_yyyy_yyyyz, ts_yyyy_yyyyzz, ts_yyyy_yyyzz, ts_yyyy_yyyzzz, ts_yyyy_yyzzz, ts_yyyy_yyzzzz, ts_yyyy_yzzzz, ts_yyyy_yzzzzz, ts_yyyyz_xxxxxx, ts_yyyyz_xxxxxy, ts_yyyyz_xxxxxz, ts_yyyyz_xxxxyy, ts_yyyyz_xxxxyz, ts_yyyyz_xxxxzz, ts_yyyyz_xxxyyy, ts_yyyyz_xxxyyz, ts_yyyyz_xxxyzz, ts_yyyyz_xxxzzz, ts_yyyyz_xxyyyy, ts_yyyyz_xxyyyz, ts_yyyyz_xxyyzz, ts_yyyyz_xxyzzz, ts_yyyyz_xxzzzz, ts_yyyyz_xyyyyy, ts_yyyyz_xyyyyz, ts_yyyyz_xyyyzz, ts_yyyyz_xyyzzz, ts_yyyyz_xyzzzz, ts_yyyyz_xzzzzz, ts_yyyyz_yyyyyy, ts_yyyyz_yyyyyz, ts_yyyyz_yyyyzz, ts_yyyyz_yyyzzz, ts_yyyyz_yyzzzz, ts_yyyyz_yzzzzz, ts_yyyyz_zzzzzz, ts_yyyz_xxxxxz, ts_yyyz_xxxxzz, ts_yyyz_xxxzzz, ts_yyyz_xxzzzz, ts_yyyz_xzzzzz, ts_yyyz_zzzzzz, ts_yyz_xxxxxz, ts_yyz_xxxxzz, ts_yyz_xxxzzz, ts_yyz_xxzzzz, ts_yyz_xzzzzz, ts_yyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyz_xxxxxx[i] = ts_yyyy_xxxxxx[i] * pa_z[i];

        ts_yyyyz_xxxxxy[i] = ts_yyyy_xxxxxy[i] * pa_z[i];

        ts_yyyyz_xxxxxz[i] = 3.0 * ts_yyz_xxxxxz[i] * fe_0 + ts_yyyz_xxxxxz[i] * pa_y[i];

        ts_yyyyz_xxxxyy[i] = ts_yyyy_xxxxyy[i] * pa_z[i];

        ts_yyyyz_xxxxyz[i] = ts_yyyy_xxxxy[i] * fe_0 + ts_yyyy_xxxxyz[i] * pa_z[i];

        ts_yyyyz_xxxxzz[i] = 3.0 * ts_yyz_xxxxzz[i] * fe_0 + ts_yyyz_xxxxzz[i] * pa_y[i];

        ts_yyyyz_xxxyyy[i] = ts_yyyy_xxxyyy[i] * pa_z[i];

        ts_yyyyz_xxxyyz[i] = ts_yyyy_xxxyy[i] * fe_0 + ts_yyyy_xxxyyz[i] * pa_z[i];

        ts_yyyyz_xxxyzz[i] = 2.0 * ts_yyyy_xxxyz[i] * fe_0 + ts_yyyy_xxxyzz[i] * pa_z[i];

        ts_yyyyz_xxxzzz[i] = 3.0 * ts_yyz_xxxzzz[i] * fe_0 + ts_yyyz_xxxzzz[i] * pa_y[i];

        ts_yyyyz_xxyyyy[i] = ts_yyyy_xxyyyy[i] * pa_z[i];

        ts_yyyyz_xxyyyz[i] = ts_yyyy_xxyyy[i] * fe_0 + ts_yyyy_xxyyyz[i] * pa_z[i];

        ts_yyyyz_xxyyzz[i] = 2.0 * ts_yyyy_xxyyz[i] * fe_0 + ts_yyyy_xxyyzz[i] * pa_z[i];

        ts_yyyyz_xxyzzz[i] = 3.0 * ts_yyyy_xxyzz[i] * fe_0 + ts_yyyy_xxyzzz[i] * pa_z[i];

        ts_yyyyz_xxzzzz[i] = 3.0 * ts_yyz_xxzzzz[i] * fe_0 + ts_yyyz_xxzzzz[i] * pa_y[i];

        ts_yyyyz_xyyyyy[i] = ts_yyyy_xyyyyy[i] * pa_z[i];

        ts_yyyyz_xyyyyz[i] = ts_yyyy_xyyyy[i] * fe_0 + ts_yyyy_xyyyyz[i] * pa_z[i];

        ts_yyyyz_xyyyzz[i] = 2.0 * ts_yyyy_xyyyz[i] * fe_0 + ts_yyyy_xyyyzz[i] * pa_z[i];

        ts_yyyyz_xyyzzz[i] = 3.0 * ts_yyyy_xyyzz[i] * fe_0 + ts_yyyy_xyyzzz[i] * pa_z[i];

        ts_yyyyz_xyzzzz[i] = 4.0 * ts_yyyy_xyzzz[i] * fe_0 + ts_yyyy_xyzzzz[i] * pa_z[i];

        ts_yyyyz_xzzzzz[i] = 3.0 * ts_yyz_xzzzzz[i] * fe_0 + ts_yyyz_xzzzzz[i] * pa_y[i];

        ts_yyyyz_yyyyyy[i] = ts_yyyy_yyyyyy[i] * pa_z[i];

        ts_yyyyz_yyyyyz[i] = ts_yyyy_yyyyy[i] * fe_0 + ts_yyyy_yyyyyz[i] * pa_z[i];

        ts_yyyyz_yyyyzz[i] = 2.0 * ts_yyyy_yyyyz[i] * fe_0 + ts_yyyy_yyyyzz[i] * pa_z[i];

        ts_yyyyz_yyyzzz[i] = 3.0 * ts_yyyy_yyyzz[i] * fe_0 + ts_yyyy_yyyzzz[i] * pa_z[i];

        ts_yyyyz_yyzzzz[i] = 4.0 * ts_yyyy_yyzzz[i] * fe_0 + ts_yyyy_yyzzzz[i] * pa_z[i];

        ts_yyyyz_yzzzzz[i] = 5.0 * ts_yyyy_yzzzz[i] * fe_0 + ts_yyyy_yzzzzz[i] * pa_z[i];

        ts_yyyyz_zzzzzz[i] = 3.0 * ts_yyz_zzzzzz[i] * fe_0 + ts_yyyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 476-504 components of targeted buffer : HI

    auto ts_yyyzz_xxxxxx = pbuffer.data(idx_ovl_hi + 476);

    auto ts_yyyzz_xxxxxy = pbuffer.data(idx_ovl_hi + 477);

    auto ts_yyyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 478);

    auto ts_yyyzz_xxxxyy = pbuffer.data(idx_ovl_hi + 479);

    auto ts_yyyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 480);

    auto ts_yyyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 481);

    auto ts_yyyzz_xxxyyy = pbuffer.data(idx_ovl_hi + 482);

    auto ts_yyyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 483);

    auto ts_yyyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 484);

    auto ts_yyyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 485);

    auto ts_yyyzz_xxyyyy = pbuffer.data(idx_ovl_hi + 486);

    auto ts_yyyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 487);

    auto ts_yyyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 488);

    auto ts_yyyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 489);

    auto ts_yyyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 490);

    auto ts_yyyzz_xyyyyy = pbuffer.data(idx_ovl_hi + 491);

    auto ts_yyyzz_xyyyyz = pbuffer.data(idx_ovl_hi + 492);

    auto ts_yyyzz_xyyyzz = pbuffer.data(idx_ovl_hi + 493);

    auto ts_yyyzz_xyyzzz = pbuffer.data(idx_ovl_hi + 494);

    auto ts_yyyzz_xyzzzz = pbuffer.data(idx_ovl_hi + 495);

    auto ts_yyyzz_xzzzzz = pbuffer.data(idx_ovl_hi + 496);

    auto ts_yyyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 497);

    auto ts_yyyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 498);

    auto ts_yyyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 499);

    auto ts_yyyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 500);

    auto ts_yyyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 501);

    auto ts_yyyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 502);

    auto ts_yyyzz_zzzzzz = pbuffer.data(idx_ovl_hi + 503);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyy_xxxxxy, ts_yyy_xxxxyy, ts_yyy_xxxyyy, ts_yyy_xxyyyy, ts_yyy_xyyyyy, ts_yyy_yyyyyy, ts_yyyz_xxxxxy, ts_yyyz_xxxxyy, ts_yyyz_xxxyyy, ts_yyyz_xxyyyy, ts_yyyz_xyyyyy, ts_yyyz_yyyyyy, ts_yyyzz_xxxxxx, ts_yyyzz_xxxxxy, ts_yyyzz_xxxxxz, ts_yyyzz_xxxxyy, ts_yyyzz_xxxxyz, ts_yyyzz_xxxxzz, ts_yyyzz_xxxyyy, ts_yyyzz_xxxyyz, ts_yyyzz_xxxyzz, ts_yyyzz_xxxzzz, ts_yyyzz_xxyyyy, ts_yyyzz_xxyyyz, ts_yyyzz_xxyyzz, ts_yyyzz_xxyzzz, ts_yyyzz_xxzzzz, ts_yyyzz_xyyyyy, ts_yyyzz_xyyyyz, ts_yyyzz_xyyyzz, ts_yyyzz_xyyzzz, ts_yyyzz_xyzzzz, ts_yyyzz_xzzzzz, ts_yyyzz_yyyyyy, ts_yyyzz_yyyyyz, ts_yyyzz_yyyyzz, ts_yyyzz_yyyzzz, ts_yyyzz_yyzzzz, ts_yyyzz_yzzzzz, ts_yyyzz_zzzzzz, ts_yyzz_xxxxxx, ts_yyzz_xxxxxz, ts_yyzz_xxxxyz, ts_yyzz_xxxxz, ts_yyzz_xxxxzz, ts_yyzz_xxxyyz, ts_yyzz_xxxyz, ts_yyzz_xxxyzz, ts_yyzz_xxxzz, ts_yyzz_xxxzzz, ts_yyzz_xxyyyz, ts_yyzz_xxyyz, ts_yyzz_xxyyzz, ts_yyzz_xxyzz, ts_yyzz_xxyzzz, ts_yyzz_xxzzz, ts_yyzz_xxzzzz, ts_yyzz_xyyyyz, ts_yyzz_xyyyz, ts_yyzz_xyyyzz, ts_yyzz_xyyzz, ts_yyzz_xyyzzz, ts_yyzz_xyzzz, ts_yyzz_xyzzzz, ts_yyzz_xzzzz, ts_yyzz_xzzzzz, ts_yyzz_yyyyyz, ts_yyzz_yyyyz, ts_yyzz_yyyyzz, ts_yyzz_yyyzz, ts_yyzz_yyyzzz, ts_yyzz_yyzzz, ts_yyzz_yyzzzz, ts_yyzz_yzzzz, ts_yyzz_yzzzzz, ts_yyzz_zzzzz, ts_yyzz_zzzzzz, ts_yzz_xxxxxx, ts_yzz_xxxxxz, ts_yzz_xxxxyz, ts_yzz_xxxxzz, ts_yzz_xxxyyz, ts_yzz_xxxyzz, ts_yzz_xxxzzz, ts_yzz_xxyyyz, ts_yzz_xxyyzz, ts_yzz_xxyzzz, ts_yzz_xxzzzz, ts_yzz_xyyyyz, ts_yzz_xyyyzz, ts_yzz_xyyzzz, ts_yzz_xyzzzz, ts_yzz_xzzzzz, ts_yzz_yyyyyz, ts_yzz_yyyyzz, ts_yzz_yyyzzz, ts_yzz_yyzzzz, ts_yzz_yzzzzz, ts_yzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzz_xxxxxx[i] = 2.0 * ts_yzz_xxxxxx[i] * fe_0 + ts_yyzz_xxxxxx[i] * pa_y[i];

        ts_yyyzz_xxxxxy[i] = ts_yyy_xxxxxy[i] * fe_0 + ts_yyyz_xxxxxy[i] * pa_z[i];

        ts_yyyzz_xxxxxz[i] = 2.0 * ts_yzz_xxxxxz[i] * fe_0 + ts_yyzz_xxxxxz[i] * pa_y[i];

        ts_yyyzz_xxxxyy[i] = ts_yyy_xxxxyy[i] * fe_0 + ts_yyyz_xxxxyy[i] * pa_z[i];

        ts_yyyzz_xxxxyz[i] = 2.0 * ts_yzz_xxxxyz[i] * fe_0 + ts_yyzz_xxxxz[i] * fe_0 + ts_yyzz_xxxxyz[i] * pa_y[i];

        ts_yyyzz_xxxxzz[i] = 2.0 * ts_yzz_xxxxzz[i] * fe_0 + ts_yyzz_xxxxzz[i] * pa_y[i];

        ts_yyyzz_xxxyyy[i] = ts_yyy_xxxyyy[i] * fe_0 + ts_yyyz_xxxyyy[i] * pa_z[i];

        ts_yyyzz_xxxyyz[i] = 2.0 * ts_yzz_xxxyyz[i] * fe_0 + 2.0 * ts_yyzz_xxxyz[i] * fe_0 + ts_yyzz_xxxyyz[i] * pa_y[i];

        ts_yyyzz_xxxyzz[i] = 2.0 * ts_yzz_xxxyzz[i] * fe_0 + ts_yyzz_xxxzz[i] * fe_0 + ts_yyzz_xxxyzz[i] * pa_y[i];

        ts_yyyzz_xxxzzz[i] = 2.0 * ts_yzz_xxxzzz[i] * fe_0 + ts_yyzz_xxxzzz[i] * pa_y[i];

        ts_yyyzz_xxyyyy[i] = ts_yyy_xxyyyy[i] * fe_0 + ts_yyyz_xxyyyy[i] * pa_z[i];

        ts_yyyzz_xxyyyz[i] = 2.0 * ts_yzz_xxyyyz[i] * fe_0 + 3.0 * ts_yyzz_xxyyz[i] * fe_0 + ts_yyzz_xxyyyz[i] * pa_y[i];

        ts_yyyzz_xxyyzz[i] = 2.0 * ts_yzz_xxyyzz[i] * fe_0 + 2.0 * ts_yyzz_xxyzz[i] * fe_0 + ts_yyzz_xxyyzz[i] * pa_y[i];

        ts_yyyzz_xxyzzz[i] = 2.0 * ts_yzz_xxyzzz[i] * fe_0 + ts_yyzz_xxzzz[i] * fe_0 + ts_yyzz_xxyzzz[i] * pa_y[i];

        ts_yyyzz_xxzzzz[i] = 2.0 * ts_yzz_xxzzzz[i] * fe_0 + ts_yyzz_xxzzzz[i] * pa_y[i];

        ts_yyyzz_xyyyyy[i] = ts_yyy_xyyyyy[i] * fe_0 + ts_yyyz_xyyyyy[i] * pa_z[i];

        ts_yyyzz_xyyyyz[i] = 2.0 * ts_yzz_xyyyyz[i] * fe_0 + 4.0 * ts_yyzz_xyyyz[i] * fe_0 + ts_yyzz_xyyyyz[i] * pa_y[i];

        ts_yyyzz_xyyyzz[i] = 2.0 * ts_yzz_xyyyzz[i] * fe_0 + 3.0 * ts_yyzz_xyyzz[i] * fe_0 + ts_yyzz_xyyyzz[i] * pa_y[i];

        ts_yyyzz_xyyzzz[i] = 2.0 * ts_yzz_xyyzzz[i] * fe_0 + 2.0 * ts_yyzz_xyzzz[i] * fe_0 + ts_yyzz_xyyzzz[i] * pa_y[i];

        ts_yyyzz_xyzzzz[i] = 2.0 * ts_yzz_xyzzzz[i] * fe_0 + ts_yyzz_xzzzz[i] * fe_0 + ts_yyzz_xyzzzz[i] * pa_y[i];

        ts_yyyzz_xzzzzz[i] = 2.0 * ts_yzz_xzzzzz[i] * fe_0 + ts_yyzz_xzzzzz[i] * pa_y[i];

        ts_yyyzz_yyyyyy[i] = ts_yyy_yyyyyy[i] * fe_0 + ts_yyyz_yyyyyy[i] * pa_z[i];

        ts_yyyzz_yyyyyz[i] = 2.0 * ts_yzz_yyyyyz[i] * fe_0 + 5.0 * ts_yyzz_yyyyz[i] * fe_0 + ts_yyzz_yyyyyz[i] * pa_y[i];

        ts_yyyzz_yyyyzz[i] = 2.0 * ts_yzz_yyyyzz[i] * fe_0 + 4.0 * ts_yyzz_yyyzz[i] * fe_0 + ts_yyzz_yyyyzz[i] * pa_y[i];

        ts_yyyzz_yyyzzz[i] = 2.0 * ts_yzz_yyyzzz[i] * fe_0 + 3.0 * ts_yyzz_yyzzz[i] * fe_0 + ts_yyzz_yyyzzz[i] * pa_y[i];

        ts_yyyzz_yyzzzz[i] = 2.0 * ts_yzz_yyzzzz[i] * fe_0 + 2.0 * ts_yyzz_yzzzz[i] * fe_0 + ts_yyzz_yyzzzz[i] * pa_y[i];

        ts_yyyzz_yzzzzz[i] = 2.0 * ts_yzz_yzzzzz[i] * fe_0 + ts_yyzz_zzzzz[i] * fe_0 + ts_yyzz_yzzzzz[i] * pa_y[i];

        ts_yyyzz_zzzzzz[i] = 2.0 * ts_yzz_zzzzzz[i] * fe_0 + ts_yyzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 504-532 components of targeted buffer : HI

    auto ts_yyzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 504);

    auto ts_yyzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 505);

    auto ts_yyzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 506);

    auto ts_yyzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 507);

    auto ts_yyzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 508);

    auto ts_yyzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 509);

    auto ts_yyzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 510);

    auto ts_yyzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 511);

    auto ts_yyzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 512);

    auto ts_yyzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 513);

    auto ts_yyzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 514);

    auto ts_yyzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 515);

    auto ts_yyzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 516);

    auto ts_yyzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 517);

    auto ts_yyzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 518);

    auto ts_yyzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 519);

    auto ts_yyzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 520);

    auto ts_yyzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 521);

    auto ts_yyzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 522);

    auto ts_yyzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 523);

    auto ts_yyzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 524);

    auto ts_yyzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 525);

    auto ts_yyzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 526);

    auto ts_yyzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 527);

    auto ts_yyzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 528);

    auto ts_yyzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 529);

    auto ts_yyzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 530);

    auto ts_yyzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 531);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyz_xxxxxy, ts_yyz_xxxxyy, ts_yyz_xxxyyy, ts_yyz_xxyyyy, ts_yyz_xyyyyy, ts_yyz_yyyyyy, ts_yyzz_xxxxxy, ts_yyzz_xxxxyy, ts_yyzz_xxxyyy, ts_yyzz_xxyyyy, ts_yyzz_xyyyyy, ts_yyzz_yyyyyy, ts_yyzzz_xxxxxx, ts_yyzzz_xxxxxy, ts_yyzzz_xxxxxz, ts_yyzzz_xxxxyy, ts_yyzzz_xxxxyz, ts_yyzzz_xxxxzz, ts_yyzzz_xxxyyy, ts_yyzzz_xxxyyz, ts_yyzzz_xxxyzz, ts_yyzzz_xxxzzz, ts_yyzzz_xxyyyy, ts_yyzzz_xxyyyz, ts_yyzzz_xxyyzz, ts_yyzzz_xxyzzz, ts_yyzzz_xxzzzz, ts_yyzzz_xyyyyy, ts_yyzzz_xyyyyz, ts_yyzzz_xyyyzz, ts_yyzzz_xyyzzz, ts_yyzzz_xyzzzz, ts_yyzzz_xzzzzz, ts_yyzzz_yyyyyy, ts_yyzzz_yyyyyz, ts_yyzzz_yyyyzz, ts_yyzzz_yyyzzz, ts_yyzzz_yyzzzz, ts_yyzzz_yzzzzz, ts_yyzzz_zzzzzz, ts_yzzz_xxxxxx, ts_yzzz_xxxxxz, ts_yzzz_xxxxyz, ts_yzzz_xxxxz, ts_yzzz_xxxxzz, ts_yzzz_xxxyyz, ts_yzzz_xxxyz, ts_yzzz_xxxyzz, ts_yzzz_xxxzz, ts_yzzz_xxxzzz, ts_yzzz_xxyyyz, ts_yzzz_xxyyz, ts_yzzz_xxyyzz, ts_yzzz_xxyzz, ts_yzzz_xxyzzz, ts_yzzz_xxzzz, ts_yzzz_xxzzzz, ts_yzzz_xyyyyz, ts_yzzz_xyyyz, ts_yzzz_xyyyzz, ts_yzzz_xyyzz, ts_yzzz_xyyzzz, ts_yzzz_xyzzz, ts_yzzz_xyzzzz, ts_yzzz_xzzzz, ts_yzzz_xzzzzz, ts_yzzz_yyyyyz, ts_yzzz_yyyyz, ts_yzzz_yyyyzz, ts_yzzz_yyyzz, ts_yzzz_yyyzzz, ts_yzzz_yyzzz, ts_yzzz_yyzzzz, ts_yzzz_yzzzz, ts_yzzz_yzzzzz, ts_yzzz_zzzzz, ts_yzzz_zzzzzz, ts_zzz_xxxxxx, ts_zzz_xxxxxz, ts_zzz_xxxxyz, ts_zzz_xxxxzz, ts_zzz_xxxyyz, ts_zzz_xxxyzz, ts_zzz_xxxzzz, ts_zzz_xxyyyz, ts_zzz_xxyyzz, ts_zzz_xxyzzz, ts_zzz_xxzzzz, ts_zzz_xyyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzz_xxxxxx[i] = ts_zzz_xxxxxx[i] * fe_0 + ts_yzzz_xxxxxx[i] * pa_y[i];

        ts_yyzzz_xxxxxy[i] = 2.0 * ts_yyz_xxxxxy[i] * fe_0 + ts_yyzz_xxxxxy[i] * pa_z[i];

        ts_yyzzz_xxxxxz[i] = ts_zzz_xxxxxz[i] * fe_0 + ts_yzzz_xxxxxz[i] * pa_y[i];

        ts_yyzzz_xxxxyy[i] = 2.0 * ts_yyz_xxxxyy[i] * fe_0 + ts_yyzz_xxxxyy[i] * pa_z[i];

        ts_yyzzz_xxxxyz[i] = ts_zzz_xxxxyz[i] * fe_0 + ts_yzzz_xxxxz[i] * fe_0 + ts_yzzz_xxxxyz[i] * pa_y[i];

        ts_yyzzz_xxxxzz[i] = ts_zzz_xxxxzz[i] * fe_0 + ts_yzzz_xxxxzz[i] * pa_y[i];

        ts_yyzzz_xxxyyy[i] = 2.0 * ts_yyz_xxxyyy[i] * fe_0 + ts_yyzz_xxxyyy[i] * pa_z[i];

        ts_yyzzz_xxxyyz[i] = ts_zzz_xxxyyz[i] * fe_0 + 2.0 * ts_yzzz_xxxyz[i] * fe_0 + ts_yzzz_xxxyyz[i] * pa_y[i];

        ts_yyzzz_xxxyzz[i] = ts_zzz_xxxyzz[i] * fe_0 + ts_yzzz_xxxzz[i] * fe_0 + ts_yzzz_xxxyzz[i] * pa_y[i];

        ts_yyzzz_xxxzzz[i] = ts_zzz_xxxzzz[i] * fe_0 + ts_yzzz_xxxzzz[i] * pa_y[i];

        ts_yyzzz_xxyyyy[i] = 2.0 * ts_yyz_xxyyyy[i] * fe_0 + ts_yyzz_xxyyyy[i] * pa_z[i];

        ts_yyzzz_xxyyyz[i] = ts_zzz_xxyyyz[i] * fe_0 + 3.0 * ts_yzzz_xxyyz[i] * fe_0 + ts_yzzz_xxyyyz[i] * pa_y[i];

        ts_yyzzz_xxyyzz[i] = ts_zzz_xxyyzz[i] * fe_0 + 2.0 * ts_yzzz_xxyzz[i] * fe_0 + ts_yzzz_xxyyzz[i] * pa_y[i];

        ts_yyzzz_xxyzzz[i] = ts_zzz_xxyzzz[i] * fe_0 + ts_yzzz_xxzzz[i] * fe_0 + ts_yzzz_xxyzzz[i] * pa_y[i];

        ts_yyzzz_xxzzzz[i] = ts_zzz_xxzzzz[i] * fe_0 + ts_yzzz_xxzzzz[i] * pa_y[i];

        ts_yyzzz_xyyyyy[i] = 2.0 * ts_yyz_xyyyyy[i] * fe_0 + ts_yyzz_xyyyyy[i] * pa_z[i];

        ts_yyzzz_xyyyyz[i] = ts_zzz_xyyyyz[i] * fe_0 + 4.0 * ts_yzzz_xyyyz[i] * fe_0 + ts_yzzz_xyyyyz[i] * pa_y[i];

        ts_yyzzz_xyyyzz[i] = ts_zzz_xyyyzz[i] * fe_0 + 3.0 * ts_yzzz_xyyzz[i] * fe_0 + ts_yzzz_xyyyzz[i] * pa_y[i];

        ts_yyzzz_xyyzzz[i] = ts_zzz_xyyzzz[i] * fe_0 + 2.0 * ts_yzzz_xyzzz[i] * fe_0 + ts_yzzz_xyyzzz[i] * pa_y[i];

        ts_yyzzz_xyzzzz[i] = ts_zzz_xyzzzz[i] * fe_0 + ts_yzzz_xzzzz[i] * fe_0 + ts_yzzz_xyzzzz[i] * pa_y[i];

        ts_yyzzz_xzzzzz[i] = ts_zzz_xzzzzz[i] * fe_0 + ts_yzzz_xzzzzz[i] * pa_y[i];

        ts_yyzzz_yyyyyy[i] = 2.0 * ts_yyz_yyyyyy[i] * fe_0 + ts_yyzz_yyyyyy[i] * pa_z[i];

        ts_yyzzz_yyyyyz[i] = ts_zzz_yyyyyz[i] * fe_0 + 5.0 * ts_yzzz_yyyyz[i] * fe_0 + ts_yzzz_yyyyyz[i] * pa_y[i];

        ts_yyzzz_yyyyzz[i] = ts_zzz_yyyyzz[i] * fe_0 + 4.0 * ts_yzzz_yyyzz[i] * fe_0 + ts_yzzz_yyyyzz[i] * pa_y[i];

        ts_yyzzz_yyyzzz[i] = ts_zzz_yyyzzz[i] * fe_0 + 3.0 * ts_yzzz_yyzzz[i] * fe_0 + ts_yzzz_yyyzzz[i] * pa_y[i];

        ts_yyzzz_yyzzzz[i] = ts_zzz_yyzzzz[i] * fe_0 + 2.0 * ts_yzzz_yzzzz[i] * fe_0 + ts_yzzz_yyzzzz[i] * pa_y[i];

        ts_yyzzz_yzzzzz[i] = ts_zzz_yzzzzz[i] * fe_0 + ts_yzzz_zzzzz[i] * fe_0 + ts_yzzz_yzzzzz[i] * pa_y[i];

        ts_yyzzz_zzzzzz[i] = ts_zzz_zzzzzz[i] * fe_0 + ts_yzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 532-560 components of targeted buffer : HI

    auto ts_yzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 532);

    auto ts_yzzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 533);

    auto ts_yzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 534);

    auto ts_yzzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 535);

    auto ts_yzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 536);

    auto ts_yzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 537);

    auto ts_yzzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 538);

    auto ts_yzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 539);

    auto ts_yzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 540);

    auto ts_yzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 541);

    auto ts_yzzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 542);

    auto ts_yzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 543);

    auto ts_yzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 544);

    auto ts_yzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 545);

    auto ts_yzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 546);

    auto ts_yzzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 547);

    auto ts_yzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 548);

    auto ts_yzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 549);

    auto ts_yzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 550);

    auto ts_yzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 551);

    auto ts_yzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 552);

    auto ts_yzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 553);

    auto ts_yzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 554);

    auto ts_yzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 555);

    auto ts_yzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 556);

    auto ts_yzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 557);

    auto ts_yzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 558);

    auto ts_yzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 559);

    #pragma omp simd aligned(pa_y, ts_yzzzz_xxxxxx, ts_yzzzz_xxxxxy, ts_yzzzz_xxxxxz, ts_yzzzz_xxxxyy, ts_yzzzz_xxxxyz, ts_yzzzz_xxxxzz, ts_yzzzz_xxxyyy, ts_yzzzz_xxxyyz, ts_yzzzz_xxxyzz, ts_yzzzz_xxxzzz, ts_yzzzz_xxyyyy, ts_yzzzz_xxyyyz, ts_yzzzz_xxyyzz, ts_yzzzz_xxyzzz, ts_yzzzz_xxzzzz, ts_yzzzz_xyyyyy, ts_yzzzz_xyyyyz, ts_yzzzz_xyyyzz, ts_yzzzz_xyyzzz, ts_yzzzz_xyzzzz, ts_yzzzz_xzzzzz, ts_yzzzz_yyyyyy, ts_yzzzz_yyyyyz, ts_yzzzz_yyyyzz, ts_yzzzz_yyyzzz, ts_yzzzz_yyzzzz, ts_yzzzz_yzzzzz, ts_yzzzz_zzzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxxx, ts_zzzz_xxxxxy, ts_zzzz_xxxxxz, ts_zzzz_xxxxy, ts_zzzz_xxxxyy, ts_zzzz_xxxxyz, ts_zzzz_xxxxz, ts_zzzz_xxxxzz, ts_zzzz_xxxyy, ts_zzzz_xxxyyy, ts_zzzz_xxxyyz, ts_zzzz_xxxyz, ts_zzzz_xxxyzz, ts_zzzz_xxxzz, ts_zzzz_xxxzzz, ts_zzzz_xxyyy, ts_zzzz_xxyyyy, ts_zzzz_xxyyyz, ts_zzzz_xxyyz, ts_zzzz_xxyyzz, ts_zzzz_xxyzz, ts_zzzz_xxyzzz, ts_zzzz_xxzzz, ts_zzzz_xxzzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyyy, ts_zzzz_xyyyyz, ts_zzzz_xyyyz, ts_zzzz_xyyyzz, ts_zzzz_xyyzz, ts_zzzz_xyyzzz, ts_zzzz_xyzzz, ts_zzzz_xyzzzz, ts_zzzz_xzzzz, ts_zzzz_xzzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyyy, ts_zzzz_yyyyyz, ts_zzzz_yyyyz, ts_zzzz_yyyyzz, ts_zzzz_yyyzz, ts_zzzz_yyyzzz, ts_zzzz_yyzzz, ts_zzzz_yyzzzz, ts_zzzz_yzzzz, ts_zzzz_yzzzzz, ts_zzzz_zzzzz, ts_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzz_xxxxxx[i] = ts_zzzz_xxxxxx[i] * pa_y[i];

        ts_yzzzz_xxxxxy[i] = ts_zzzz_xxxxx[i] * fe_0 + ts_zzzz_xxxxxy[i] * pa_y[i];

        ts_yzzzz_xxxxxz[i] = ts_zzzz_xxxxxz[i] * pa_y[i];

        ts_yzzzz_xxxxyy[i] = 2.0 * ts_zzzz_xxxxy[i] * fe_0 + ts_zzzz_xxxxyy[i] * pa_y[i];

        ts_yzzzz_xxxxyz[i] = ts_zzzz_xxxxz[i] * fe_0 + ts_zzzz_xxxxyz[i] * pa_y[i];

        ts_yzzzz_xxxxzz[i] = ts_zzzz_xxxxzz[i] * pa_y[i];

        ts_yzzzz_xxxyyy[i] = 3.0 * ts_zzzz_xxxyy[i] * fe_0 + ts_zzzz_xxxyyy[i] * pa_y[i];

        ts_yzzzz_xxxyyz[i] = 2.0 * ts_zzzz_xxxyz[i] * fe_0 + ts_zzzz_xxxyyz[i] * pa_y[i];

        ts_yzzzz_xxxyzz[i] = ts_zzzz_xxxzz[i] * fe_0 + ts_zzzz_xxxyzz[i] * pa_y[i];

        ts_yzzzz_xxxzzz[i] = ts_zzzz_xxxzzz[i] * pa_y[i];

        ts_yzzzz_xxyyyy[i] = 4.0 * ts_zzzz_xxyyy[i] * fe_0 + ts_zzzz_xxyyyy[i] * pa_y[i];

        ts_yzzzz_xxyyyz[i] = 3.0 * ts_zzzz_xxyyz[i] * fe_0 + ts_zzzz_xxyyyz[i] * pa_y[i];

        ts_yzzzz_xxyyzz[i] = 2.0 * ts_zzzz_xxyzz[i] * fe_0 + ts_zzzz_xxyyzz[i] * pa_y[i];

        ts_yzzzz_xxyzzz[i] = ts_zzzz_xxzzz[i] * fe_0 + ts_zzzz_xxyzzz[i] * pa_y[i];

        ts_yzzzz_xxzzzz[i] = ts_zzzz_xxzzzz[i] * pa_y[i];

        ts_yzzzz_xyyyyy[i] = 5.0 * ts_zzzz_xyyyy[i] * fe_0 + ts_zzzz_xyyyyy[i] * pa_y[i];

        ts_yzzzz_xyyyyz[i] = 4.0 * ts_zzzz_xyyyz[i] * fe_0 + ts_zzzz_xyyyyz[i] * pa_y[i];

        ts_yzzzz_xyyyzz[i] = 3.0 * ts_zzzz_xyyzz[i] * fe_0 + ts_zzzz_xyyyzz[i] * pa_y[i];

        ts_yzzzz_xyyzzz[i] = 2.0 * ts_zzzz_xyzzz[i] * fe_0 + ts_zzzz_xyyzzz[i] * pa_y[i];

        ts_yzzzz_xyzzzz[i] = ts_zzzz_xzzzz[i] * fe_0 + ts_zzzz_xyzzzz[i] * pa_y[i];

        ts_yzzzz_xzzzzz[i] = ts_zzzz_xzzzzz[i] * pa_y[i];

        ts_yzzzz_yyyyyy[i] = 6.0 * ts_zzzz_yyyyy[i] * fe_0 + ts_zzzz_yyyyyy[i] * pa_y[i];

        ts_yzzzz_yyyyyz[i] = 5.0 * ts_zzzz_yyyyz[i] * fe_0 + ts_zzzz_yyyyyz[i] * pa_y[i];

        ts_yzzzz_yyyyzz[i] = 4.0 * ts_zzzz_yyyzz[i] * fe_0 + ts_zzzz_yyyyzz[i] * pa_y[i];

        ts_yzzzz_yyyzzz[i] = 3.0 * ts_zzzz_yyzzz[i] * fe_0 + ts_zzzz_yyyzzz[i] * pa_y[i];

        ts_yzzzz_yyzzzz[i] = 2.0 * ts_zzzz_yzzzz[i] * fe_0 + ts_zzzz_yyzzzz[i] * pa_y[i];

        ts_yzzzz_yzzzzz[i] = ts_zzzz_zzzzz[i] * fe_0 + ts_zzzz_yzzzzz[i] * pa_y[i];

        ts_yzzzz_zzzzzz[i] = ts_zzzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 560-588 components of targeted buffer : HI

    auto ts_zzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 560);

    auto ts_zzzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 561);

    auto ts_zzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 562);

    auto ts_zzzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 563);

    auto ts_zzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 564);

    auto ts_zzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 565);

    auto ts_zzzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 566);

    auto ts_zzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 567);

    auto ts_zzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 568);

    auto ts_zzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 569);

    auto ts_zzzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 570);

    auto ts_zzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 571);

    auto ts_zzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 572);

    auto ts_zzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 573);

    auto ts_zzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 574);

    auto ts_zzzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 575);

    auto ts_zzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 576);

    auto ts_zzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 577);

    auto ts_zzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 578);

    auto ts_zzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 579);

    auto ts_zzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 580);

    auto ts_zzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 581);

    auto ts_zzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 582);

    auto ts_zzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 583);

    auto ts_zzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 584);

    auto ts_zzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 585);

    auto ts_zzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 586);

    auto ts_zzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 587);

    #pragma omp simd aligned(pa_z, ts_zzz_xxxxxx, ts_zzz_xxxxxy, ts_zzz_xxxxxz, ts_zzz_xxxxyy, ts_zzz_xxxxyz, ts_zzz_xxxxzz, ts_zzz_xxxyyy, ts_zzz_xxxyyz, ts_zzz_xxxyzz, ts_zzz_xxxzzz, ts_zzz_xxyyyy, ts_zzz_xxyyyz, ts_zzz_xxyyzz, ts_zzz_xxyzzz, ts_zzz_xxzzzz, ts_zzz_xyyyyy, ts_zzz_xyyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyyy, ts_zzz_yyyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxxx, ts_zzzz_xxxxxy, ts_zzzz_xxxxxz, ts_zzzz_xxxxy, ts_zzzz_xxxxyy, ts_zzzz_xxxxyz, ts_zzzz_xxxxz, ts_zzzz_xxxxzz, ts_zzzz_xxxyy, ts_zzzz_xxxyyy, ts_zzzz_xxxyyz, ts_zzzz_xxxyz, ts_zzzz_xxxyzz, ts_zzzz_xxxzz, ts_zzzz_xxxzzz, ts_zzzz_xxyyy, ts_zzzz_xxyyyy, ts_zzzz_xxyyyz, ts_zzzz_xxyyz, ts_zzzz_xxyyzz, ts_zzzz_xxyzz, ts_zzzz_xxyzzz, ts_zzzz_xxzzz, ts_zzzz_xxzzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyyy, ts_zzzz_xyyyyz, ts_zzzz_xyyyz, ts_zzzz_xyyyzz, ts_zzzz_xyyzz, ts_zzzz_xyyzzz, ts_zzzz_xyzzz, ts_zzzz_xyzzzz, ts_zzzz_xzzzz, ts_zzzz_xzzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyyy, ts_zzzz_yyyyyz, ts_zzzz_yyyyz, ts_zzzz_yyyyzz, ts_zzzz_yyyzz, ts_zzzz_yyyzzz, ts_zzzz_yyzzz, ts_zzzz_yyzzzz, ts_zzzz_yzzzz, ts_zzzz_yzzzzz, ts_zzzz_zzzzz, ts_zzzz_zzzzzz, ts_zzzzz_xxxxxx, ts_zzzzz_xxxxxy, ts_zzzzz_xxxxxz, ts_zzzzz_xxxxyy, ts_zzzzz_xxxxyz, ts_zzzzz_xxxxzz, ts_zzzzz_xxxyyy, ts_zzzzz_xxxyyz, ts_zzzzz_xxxyzz, ts_zzzzz_xxxzzz, ts_zzzzz_xxyyyy, ts_zzzzz_xxyyyz, ts_zzzzz_xxyyzz, ts_zzzzz_xxyzzz, ts_zzzzz_xxzzzz, ts_zzzzz_xyyyyy, ts_zzzzz_xyyyyz, ts_zzzzz_xyyyzz, ts_zzzzz_xyyzzz, ts_zzzzz_xyzzzz, ts_zzzzz_xzzzzz, ts_zzzzz_yyyyyy, ts_zzzzz_yyyyyz, ts_zzzzz_yyyyzz, ts_zzzzz_yyyzzz, ts_zzzzz_yyzzzz, ts_zzzzz_yzzzzz, ts_zzzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzz_xxxxxx[i] = 4.0 * ts_zzz_xxxxxx[i] * fe_0 + ts_zzzz_xxxxxx[i] * pa_z[i];

        ts_zzzzz_xxxxxy[i] = 4.0 * ts_zzz_xxxxxy[i] * fe_0 + ts_zzzz_xxxxxy[i] * pa_z[i];

        ts_zzzzz_xxxxxz[i] = 4.0 * ts_zzz_xxxxxz[i] * fe_0 + ts_zzzz_xxxxx[i] * fe_0 + ts_zzzz_xxxxxz[i] * pa_z[i];

        ts_zzzzz_xxxxyy[i] = 4.0 * ts_zzz_xxxxyy[i] * fe_0 + ts_zzzz_xxxxyy[i] * pa_z[i];

        ts_zzzzz_xxxxyz[i] = 4.0 * ts_zzz_xxxxyz[i] * fe_0 + ts_zzzz_xxxxy[i] * fe_0 + ts_zzzz_xxxxyz[i] * pa_z[i];

        ts_zzzzz_xxxxzz[i] = 4.0 * ts_zzz_xxxxzz[i] * fe_0 + 2.0 * ts_zzzz_xxxxz[i] * fe_0 + ts_zzzz_xxxxzz[i] * pa_z[i];

        ts_zzzzz_xxxyyy[i] = 4.0 * ts_zzz_xxxyyy[i] * fe_0 + ts_zzzz_xxxyyy[i] * pa_z[i];

        ts_zzzzz_xxxyyz[i] = 4.0 * ts_zzz_xxxyyz[i] * fe_0 + ts_zzzz_xxxyy[i] * fe_0 + ts_zzzz_xxxyyz[i] * pa_z[i];

        ts_zzzzz_xxxyzz[i] = 4.0 * ts_zzz_xxxyzz[i] * fe_0 + 2.0 * ts_zzzz_xxxyz[i] * fe_0 + ts_zzzz_xxxyzz[i] * pa_z[i];

        ts_zzzzz_xxxzzz[i] = 4.0 * ts_zzz_xxxzzz[i] * fe_0 + 3.0 * ts_zzzz_xxxzz[i] * fe_0 + ts_zzzz_xxxzzz[i] * pa_z[i];

        ts_zzzzz_xxyyyy[i] = 4.0 * ts_zzz_xxyyyy[i] * fe_0 + ts_zzzz_xxyyyy[i] * pa_z[i];

        ts_zzzzz_xxyyyz[i] = 4.0 * ts_zzz_xxyyyz[i] * fe_0 + ts_zzzz_xxyyy[i] * fe_0 + ts_zzzz_xxyyyz[i] * pa_z[i];

        ts_zzzzz_xxyyzz[i] = 4.0 * ts_zzz_xxyyzz[i] * fe_0 + 2.0 * ts_zzzz_xxyyz[i] * fe_0 + ts_zzzz_xxyyzz[i] * pa_z[i];

        ts_zzzzz_xxyzzz[i] = 4.0 * ts_zzz_xxyzzz[i] * fe_0 + 3.0 * ts_zzzz_xxyzz[i] * fe_0 + ts_zzzz_xxyzzz[i] * pa_z[i];

        ts_zzzzz_xxzzzz[i] = 4.0 * ts_zzz_xxzzzz[i] * fe_0 + 4.0 * ts_zzzz_xxzzz[i] * fe_0 + ts_zzzz_xxzzzz[i] * pa_z[i];

        ts_zzzzz_xyyyyy[i] = 4.0 * ts_zzz_xyyyyy[i] * fe_0 + ts_zzzz_xyyyyy[i] * pa_z[i];

        ts_zzzzz_xyyyyz[i] = 4.0 * ts_zzz_xyyyyz[i] * fe_0 + ts_zzzz_xyyyy[i] * fe_0 + ts_zzzz_xyyyyz[i] * pa_z[i];

        ts_zzzzz_xyyyzz[i] = 4.0 * ts_zzz_xyyyzz[i] * fe_0 + 2.0 * ts_zzzz_xyyyz[i] * fe_0 + ts_zzzz_xyyyzz[i] * pa_z[i];

        ts_zzzzz_xyyzzz[i] = 4.0 * ts_zzz_xyyzzz[i] * fe_0 + 3.0 * ts_zzzz_xyyzz[i] * fe_0 + ts_zzzz_xyyzzz[i] * pa_z[i];

        ts_zzzzz_xyzzzz[i] = 4.0 * ts_zzz_xyzzzz[i] * fe_0 + 4.0 * ts_zzzz_xyzzz[i] * fe_0 + ts_zzzz_xyzzzz[i] * pa_z[i];

        ts_zzzzz_xzzzzz[i] = 4.0 * ts_zzz_xzzzzz[i] * fe_0 + 5.0 * ts_zzzz_xzzzz[i] * fe_0 + ts_zzzz_xzzzzz[i] * pa_z[i];

        ts_zzzzz_yyyyyy[i] = 4.0 * ts_zzz_yyyyyy[i] * fe_0 + ts_zzzz_yyyyyy[i] * pa_z[i];

        ts_zzzzz_yyyyyz[i] = 4.0 * ts_zzz_yyyyyz[i] * fe_0 + ts_zzzz_yyyyy[i] * fe_0 + ts_zzzz_yyyyyz[i] * pa_z[i];

        ts_zzzzz_yyyyzz[i] = 4.0 * ts_zzz_yyyyzz[i] * fe_0 + 2.0 * ts_zzzz_yyyyz[i] * fe_0 + ts_zzzz_yyyyzz[i] * pa_z[i];

        ts_zzzzz_yyyzzz[i] = 4.0 * ts_zzz_yyyzzz[i] * fe_0 + 3.0 * ts_zzzz_yyyzz[i] * fe_0 + ts_zzzz_yyyzzz[i] * pa_z[i];

        ts_zzzzz_yyzzzz[i] = 4.0 * ts_zzz_yyzzzz[i] * fe_0 + 4.0 * ts_zzzz_yyzzz[i] * fe_0 + ts_zzzz_yyzzzz[i] * pa_z[i];

        ts_zzzzz_yzzzzz[i] = 4.0 * ts_zzz_yzzzzz[i] * fe_0 + 5.0 * ts_zzzz_yzzzz[i] * fe_0 + ts_zzzz_yzzzzz[i] * pa_z[i];

        ts_zzzzz_zzzzzz[i] = 4.0 * ts_zzz_zzzzzz[i] * fe_0 + 6.0 * ts_zzzz_zzzzz[i] * fe_0 + ts_zzzz_zzzzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

