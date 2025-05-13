#include "ThreeCenterOverlapGradientPrimRecHH.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_hh(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hh,
                              const size_t idx_gh,
                              const size_t idx_hg,
                              const size_t idx_hh,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : GH

    auto ts_xxxx_xxxxx = pbuffer.data(idx_gh);

    auto ts_xxxx_xxxxy = pbuffer.data(idx_gh + 1);

    auto ts_xxxx_xxxxz = pbuffer.data(idx_gh + 2);

    auto ts_xxxx_xxxyy = pbuffer.data(idx_gh + 3);

    auto ts_xxxx_xxxyz = pbuffer.data(idx_gh + 4);

    auto ts_xxxx_xxxzz = pbuffer.data(idx_gh + 5);

    auto ts_xxxx_xxyyy = pbuffer.data(idx_gh + 6);

    auto ts_xxxx_xxyyz = pbuffer.data(idx_gh + 7);

    auto ts_xxxx_xxyzz = pbuffer.data(idx_gh + 8);

    auto ts_xxxx_xxzzz = pbuffer.data(idx_gh + 9);

    auto ts_xxxx_xyyyy = pbuffer.data(idx_gh + 10);

    auto ts_xxxx_xyyyz = pbuffer.data(idx_gh + 11);

    auto ts_xxxx_xyyzz = pbuffer.data(idx_gh + 12);

    auto ts_xxxx_xyzzz = pbuffer.data(idx_gh + 13);

    auto ts_xxxx_xzzzz = pbuffer.data(idx_gh + 14);

    auto ts_xxxx_yyyyy = pbuffer.data(idx_gh + 15);

    auto ts_xxxx_yyyyz = pbuffer.data(idx_gh + 16);

    auto ts_xxxx_yyyzz = pbuffer.data(idx_gh + 17);

    auto ts_xxxx_yyzzz = pbuffer.data(idx_gh + 18);

    auto ts_xxxx_yzzzz = pbuffer.data(idx_gh + 19);

    auto ts_xxxx_zzzzz = pbuffer.data(idx_gh + 20);

    auto ts_xxxy_xxxxx = pbuffer.data(idx_gh + 21);

    auto ts_xxxy_xxxxy = pbuffer.data(idx_gh + 22);

    auto ts_xxxy_xxxxz = pbuffer.data(idx_gh + 23);

    auto ts_xxxy_xxxyy = pbuffer.data(idx_gh + 24);

    auto ts_xxxy_xxxyz = pbuffer.data(idx_gh + 25);

    auto ts_xxxy_xxxzz = pbuffer.data(idx_gh + 26);

    auto ts_xxxy_xxyyy = pbuffer.data(idx_gh + 27);

    auto ts_xxxy_xxyyz = pbuffer.data(idx_gh + 28);

    auto ts_xxxy_xxyzz = pbuffer.data(idx_gh + 29);

    auto ts_xxxy_xxzzz = pbuffer.data(idx_gh + 30);

    auto ts_xxxy_xyyyy = pbuffer.data(idx_gh + 31);

    auto ts_xxxy_xyyyz = pbuffer.data(idx_gh + 32);

    auto ts_xxxy_xyyzz = pbuffer.data(idx_gh + 33);

    auto ts_xxxy_xyzzz = pbuffer.data(idx_gh + 34);

    auto ts_xxxy_xzzzz = pbuffer.data(idx_gh + 35);

    auto ts_xxxy_yyyyy = pbuffer.data(idx_gh + 36);

    auto ts_xxxy_yyyyz = pbuffer.data(idx_gh + 37);

    auto ts_xxxy_yyyzz = pbuffer.data(idx_gh + 38);

    auto ts_xxxy_yyzzz = pbuffer.data(idx_gh + 39);

    auto ts_xxxy_yzzzz = pbuffer.data(idx_gh + 40);

    auto ts_xxxy_zzzzz = pbuffer.data(idx_gh + 41);

    auto ts_xxxz_xxxxx = pbuffer.data(idx_gh + 42);

    auto ts_xxxz_xxxxy = pbuffer.data(idx_gh + 43);

    auto ts_xxxz_xxxxz = pbuffer.data(idx_gh + 44);

    auto ts_xxxz_xxxyy = pbuffer.data(idx_gh + 45);

    auto ts_xxxz_xxxyz = pbuffer.data(idx_gh + 46);

    auto ts_xxxz_xxxzz = pbuffer.data(idx_gh + 47);

    auto ts_xxxz_xxyyy = pbuffer.data(idx_gh + 48);

    auto ts_xxxz_xxyyz = pbuffer.data(idx_gh + 49);

    auto ts_xxxz_xxyzz = pbuffer.data(idx_gh + 50);

    auto ts_xxxz_xxzzz = pbuffer.data(idx_gh + 51);

    auto ts_xxxz_xyyyy = pbuffer.data(idx_gh + 52);

    auto ts_xxxz_xyyyz = pbuffer.data(idx_gh + 53);

    auto ts_xxxz_xyyzz = pbuffer.data(idx_gh + 54);

    auto ts_xxxz_xyzzz = pbuffer.data(idx_gh + 55);

    auto ts_xxxz_xzzzz = pbuffer.data(idx_gh + 56);

    auto ts_xxxz_yyyyy = pbuffer.data(idx_gh + 57);

    auto ts_xxxz_yyyyz = pbuffer.data(idx_gh + 58);

    auto ts_xxxz_yyyzz = pbuffer.data(idx_gh + 59);

    auto ts_xxxz_yyzzz = pbuffer.data(idx_gh + 60);

    auto ts_xxxz_yzzzz = pbuffer.data(idx_gh + 61);

    auto ts_xxxz_zzzzz = pbuffer.data(idx_gh + 62);

    auto ts_xxyy_xxxxx = pbuffer.data(idx_gh + 63);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_gh + 64);

    auto ts_xxyy_xxxxz = pbuffer.data(idx_gh + 65);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_gh + 66);

    auto ts_xxyy_xxxyz = pbuffer.data(idx_gh + 67);

    auto ts_xxyy_xxxzz = pbuffer.data(idx_gh + 68);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_gh + 69);

    auto ts_xxyy_xxyyz = pbuffer.data(idx_gh + 70);

    auto ts_xxyy_xxyzz = pbuffer.data(idx_gh + 71);

    auto ts_xxyy_xxzzz = pbuffer.data(idx_gh + 72);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_gh + 73);

    auto ts_xxyy_xyyyz = pbuffer.data(idx_gh + 74);

    auto ts_xxyy_xyyzz = pbuffer.data(idx_gh + 75);

    auto ts_xxyy_xyzzz = pbuffer.data(idx_gh + 76);

    auto ts_xxyy_xzzzz = pbuffer.data(idx_gh + 77);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_gh + 82);

    auto ts_xxyy_zzzzz = pbuffer.data(idx_gh + 83);

    auto ts_xxyz_xxxxx = pbuffer.data(idx_gh + 84);

    auto ts_xxyz_xxxxy = pbuffer.data(idx_gh + 85);

    auto ts_xxyz_xxxxz = pbuffer.data(idx_gh + 86);

    auto ts_xxyz_xxxyy = pbuffer.data(idx_gh + 87);

    auto ts_xxyz_xxxyz = pbuffer.data(idx_gh + 88);

    auto ts_xxyz_xxxzz = pbuffer.data(idx_gh + 89);

    auto ts_xxyz_xxyyy = pbuffer.data(idx_gh + 90);

    auto ts_xxyz_xxyyz = pbuffer.data(idx_gh + 91);

    auto ts_xxyz_xxyzz = pbuffer.data(idx_gh + 92);

    auto ts_xxyz_xxzzz = pbuffer.data(idx_gh + 93);

    auto ts_xxyz_xyyyy = pbuffer.data(idx_gh + 94);

    auto ts_xxyz_xyyyz = pbuffer.data(idx_gh + 95);

    auto ts_xxyz_xyyzz = pbuffer.data(idx_gh + 96);

    auto ts_xxyz_xyzzz = pbuffer.data(idx_gh + 97);

    auto ts_xxyz_xzzzz = pbuffer.data(idx_gh + 98);

    auto ts_xxyz_yyyyy = pbuffer.data(idx_gh + 99);

    auto ts_xxyz_yyyyz = pbuffer.data(idx_gh + 100);

    auto ts_xxyz_yyyzz = pbuffer.data(idx_gh + 101);

    auto ts_xxyz_yyzzz = pbuffer.data(idx_gh + 102);

    auto ts_xxyz_yzzzz = pbuffer.data(idx_gh + 103);

    auto ts_xxyz_zzzzz = pbuffer.data(idx_gh + 104);

    auto ts_xxzz_xxxxx = pbuffer.data(idx_gh + 105);

    auto ts_xxzz_xxxxy = pbuffer.data(idx_gh + 106);

    auto ts_xxzz_xxxxz = pbuffer.data(idx_gh + 107);

    auto ts_xxzz_xxxyy = pbuffer.data(idx_gh + 108);

    auto ts_xxzz_xxxyz = pbuffer.data(idx_gh + 109);

    auto ts_xxzz_xxxzz = pbuffer.data(idx_gh + 110);

    auto ts_xxzz_xxyyy = pbuffer.data(idx_gh + 111);

    auto ts_xxzz_xxyyz = pbuffer.data(idx_gh + 112);

    auto ts_xxzz_xxyzz = pbuffer.data(idx_gh + 113);

    auto ts_xxzz_xxzzz = pbuffer.data(idx_gh + 114);

    auto ts_xxzz_xyyyy = pbuffer.data(idx_gh + 115);

    auto ts_xxzz_xyyyz = pbuffer.data(idx_gh + 116);

    auto ts_xxzz_xyyzz = pbuffer.data(idx_gh + 117);

    auto ts_xxzz_xyzzz = pbuffer.data(idx_gh + 118);

    auto ts_xxzz_xzzzz = pbuffer.data(idx_gh + 119);

    auto ts_xxzz_yyyyy = pbuffer.data(idx_gh + 120);

    auto ts_xxzz_yyyyz = pbuffer.data(idx_gh + 121);

    auto ts_xxzz_yyyzz = pbuffer.data(idx_gh + 122);

    auto ts_xxzz_yyzzz = pbuffer.data(idx_gh + 123);

    auto ts_xxzz_yzzzz = pbuffer.data(idx_gh + 124);

    auto ts_xxzz_zzzzz = pbuffer.data(idx_gh + 125);

    auto ts_xyyy_xxxxx = pbuffer.data(idx_gh + 126);

    auto ts_xyyy_xxxxy = pbuffer.data(idx_gh + 127);

    auto ts_xyyy_xxxxz = pbuffer.data(idx_gh + 128);

    auto ts_xyyy_xxxyy = pbuffer.data(idx_gh + 129);

    auto ts_xyyy_xxxyz = pbuffer.data(idx_gh + 130);

    auto ts_xyyy_xxxzz = pbuffer.data(idx_gh + 131);

    auto ts_xyyy_xxyyy = pbuffer.data(idx_gh + 132);

    auto ts_xyyy_xxyyz = pbuffer.data(idx_gh + 133);

    auto ts_xyyy_xxyzz = pbuffer.data(idx_gh + 134);

    auto ts_xyyy_xxzzz = pbuffer.data(idx_gh + 135);

    auto ts_xyyy_xyyyy = pbuffer.data(idx_gh + 136);

    auto ts_xyyy_xyyyz = pbuffer.data(idx_gh + 137);

    auto ts_xyyy_xyyzz = pbuffer.data(idx_gh + 138);

    auto ts_xyyy_xyzzz = pbuffer.data(idx_gh + 139);

    auto ts_xyyy_xzzzz = pbuffer.data(idx_gh + 140);

    auto ts_xyyy_yyyyy = pbuffer.data(idx_gh + 141);

    auto ts_xyyy_yyyyz = pbuffer.data(idx_gh + 142);

    auto ts_xyyy_yyyzz = pbuffer.data(idx_gh + 143);

    auto ts_xyyy_yyzzz = pbuffer.data(idx_gh + 144);

    auto ts_xyyy_yzzzz = pbuffer.data(idx_gh + 145);

    auto ts_xyyy_zzzzz = pbuffer.data(idx_gh + 146);

    auto ts_xyyz_xxxxx = pbuffer.data(idx_gh + 147);

    auto ts_xyyz_xxxxy = pbuffer.data(idx_gh + 148);

    auto ts_xyyz_xxxxz = pbuffer.data(idx_gh + 149);

    auto ts_xyyz_xxxyy = pbuffer.data(idx_gh + 150);

    auto ts_xyyz_xxxyz = pbuffer.data(idx_gh + 151);

    auto ts_xyyz_xxxzz = pbuffer.data(idx_gh + 152);

    auto ts_xyyz_xxyyy = pbuffer.data(idx_gh + 153);

    auto ts_xyyz_xxyyz = pbuffer.data(idx_gh + 154);

    auto ts_xyyz_xxyzz = pbuffer.data(idx_gh + 155);

    auto ts_xyyz_xxzzz = pbuffer.data(idx_gh + 156);

    auto ts_xyyz_xyyyy = pbuffer.data(idx_gh + 157);

    auto ts_xyyz_xyyyz = pbuffer.data(idx_gh + 158);

    auto ts_xyyz_xyyzz = pbuffer.data(idx_gh + 159);

    auto ts_xyyz_xyzzz = pbuffer.data(idx_gh + 160);

    auto ts_xyyz_xzzzz = pbuffer.data(idx_gh + 161);

    auto ts_xyyz_yyyyy = pbuffer.data(idx_gh + 162);

    auto ts_xyyz_yyyyz = pbuffer.data(idx_gh + 163);

    auto ts_xyyz_yyyzz = pbuffer.data(idx_gh + 164);

    auto ts_xyyz_yyzzz = pbuffer.data(idx_gh + 165);

    auto ts_xyyz_yzzzz = pbuffer.data(idx_gh + 166);

    auto ts_xyyz_zzzzz = pbuffer.data(idx_gh + 167);

    auto ts_xyzz_xxxxx = pbuffer.data(idx_gh + 168);

    auto ts_xyzz_xxxxy = pbuffer.data(idx_gh + 169);

    auto ts_xyzz_xxxxz = pbuffer.data(idx_gh + 170);

    auto ts_xyzz_xxxyy = pbuffer.data(idx_gh + 171);

    auto ts_xyzz_xxxyz = pbuffer.data(idx_gh + 172);

    auto ts_xyzz_xxxzz = pbuffer.data(idx_gh + 173);

    auto ts_xyzz_xxyyy = pbuffer.data(idx_gh + 174);

    auto ts_xyzz_xxyyz = pbuffer.data(idx_gh + 175);

    auto ts_xyzz_xxyzz = pbuffer.data(idx_gh + 176);

    auto ts_xyzz_xxzzz = pbuffer.data(idx_gh + 177);

    auto ts_xyzz_xyyyy = pbuffer.data(idx_gh + 178);

    auto ts_xyzz_xyyyz = pbuffer.data(idx_gh + 179);

    auto ts_xyzz_xyyzz = pbuffer.data(idx_gh + 180);

    auto ts_xyzz_xyzzz = pbuffer.data(idx_gh + 181);

    auto ts_xyzz_xzzzz = pbuffer.data(idx_gh + 182);

    auto ts_xyzz_yyyyy = pbuffer.data(idx_gh + 183);

    auto ts_xyzz_yyyyz = pbuffer.data(idx_gh + 184);

    auto ts_xyzz_yyyzz = pbuffer.data(idx_gh + 185);

    auto ts_xyzz_yyzzz = pbuffer.data(idx_gh + 186);

    auto ts_xyzz_yzzzz = pbuffer.data(idx_gh + 187);

    auto ts_xyzz_zzzzz = pbuffer.data(idx_gh + 188);

    auto ts_xzzz_xxxxx = pbuffer.data(idx_gh + 189);

    auto ts_xzzz_xxxxy = pbuffer.data(idx_gh + 190);

    auto ts_xzzz_xxxxz = pbuffer.data(idx_gh + 191);

    auto ts_xzzz_xxxyy = pbuffer.data(idx_gh + 192);

    auto ts_xzzz_xxxyz = pbuffer.data(idx_gh + 193);

    auto ts_xzzz_xxxzz = pbuffer.data(idx_gh + 194);

    auto ts_xzzz_xxyyy = pbuffer.data(idx_gh + 195);

    auto ts_xzzz_xxyyz = pbuffer.data(idx_gh + 196);

    auto ts_xzzz_xxyzz = pbuffer.data(idx_gh + 197);

    auto ts_xzzz_xxzzz = pbuffer.data(idx_gh + 198);

    auto ts_xzzz_xyyyy = pbuffer.data(idx_gh + 199);

    auto ts_xzzz_xyyyz = pbuffer.data(idx_gh + 200);

    auto ts_xzzz_xyyzz = pbuffer.data(idx_gh + 201);

    auto ts_xzzz_xyzzz = pbuffer.data(idx_gh + 202);

    auto ts_xzzz_xzzzz = pbuffer.data(idx_gh + 203);

    auto ts_xzzz_yyyyy = pbuffer.data(idx_gh + 204);

    auto ts_xzzz_yyyyz = pbuffer.data(idx_gh + 205);

    auto ts_xzzz_yyyzz = pbuffer.data(idx_gh + 206);

    auto ts_xzzz_yyzzz = pbuffer.data(idx_gh + 207);

    auto ts_xzzz_yzzzz = pbuffer.data(idx_gh + 208);

    auto ts_xzzz_zzzzz = pbuffer.data(idx_gh + 209);

    auto ts_yyyy_xxxxx = pbuffer.data(idx_gh + 210);

    auto ts_yyyy_xxxxy = pbuffer.data(idx_gh + 211);

    auto ts_yyyy_xxxxz = pbuffer.data(idx_gh + 212);

    auto ts_yyyy_xxxyy = pbuffer.data(idx_gh + 213);

    auto ts_yyyy_xxxyz = pbuffer.data(idx_gh + 214);

    auto ts_yyyy_xxxzz = pbuffer.data(idx_gh + 215);

    auto ts_yyyy_xxyyy = pbuffer.data(idx_gh + 216);

    auto ts_yyyy_xxyyz = pbuffer.data(idx_gh + 217);

    auto ts_yyyy_xxyzz = pbuffer.data(idx_gh + 218);

    auto ts_yyyy_xxzzz = pbuffer.data(idx_gh + 219);

    auto ts_yyyy_xyyyy = pbuffer.data(idx_gh + 220);

    auto ts_yyyy_xyyyz = pbuffer.data(idx_gh + 221);

    auto ts_yyyy_xyyzz = pbuffer.data(idx_gh + 222);

    auto ts_yyyy_xyzzz = pbuffer.data(idx_gh + 223);

    auto ts_yyyy_xzzzz = pbuffer.data(idx_gh + 224);

    auto ts_yyyy_yyyyy = pbuffer.data(idx_gh + 225);

    auto ts_yyyy_yyyyz = pbuffer.data(idx_gh + 226);

    auto ts_yyyy_yyyzz = pbuffer.data(idx_gh + 227);

    auto ts_yyyy_yyzzz = pbuffer.data(idx_gh + 228);

    auto ts_yyyy_yzzzz = pbuffer.data(idx_gh + 229);

    auto ts_yyyy_zzzzz = pbuffer.data(idx_gh + 230);

    auto ts_yyyz_xxxxx = pbuffer.data(idx_gh + 231);

    auto ts_yyyz_xxxxy = pbuffer.data(idx_gh + 232);

    auto ts_yyyz_xxxxz = pbuffer.data(idx_gh + 233);

    auto ts_yyyz_xxxyy = pbuffer.data(idx_gh + 234);

    auto ts_yyyz_xxxyz = pbuffer.data(idx_gh + 235);

    auto ts_yyyz_xxxzz = pbuffer.data(idx_gh + 236);

    auto ts_yyyz_xxyyy = pbuffer.data(idx_gh + 237);

    auto ts_yyyz_xxyyz = pbuffer.data(idx_gh + 238);

    auto ts_yyyz_xxyzz = pbuffer.data(idx_gh + 239);

    auto ts_yyyz_xxzzz = pbuffer.data(idx_gh + 240);

    auto ts_yyyz_xyyyy = pbuffer.data(idx_gh + 241);

    auto ts_yyyz_xyyyz = pbuffer.data(idx_gh + 242);

    auto ts_yyyz_xyyzz = pbuffer.data(idx_gh + 243);

    auto ts_yyyz_xyzzz = pbuffer.data(idx_gh + 244);

    auto ts_yyyz_xzzzz = pbuffer.data(idx_gh + 245);

    auto ts_yyyz_yyyyy = pbuffer.data(idx_gh + 246);

    auto ts_yyyz_yyyyz = pbuffer.data(idx_gh + 247);

    auto ts_yyyz_yyyzz = pbuffer.data(idx_gh + 248);

    auto ts_yyyz_yyzzz = pbuffer.data(idx_gh + 249);

    auto ts_yyyz_yzzzz = pbuffer.data(idx_gh + 250);

    auto ts_yyyz_zzzzz = pbuffer.data(idx_gh + 251);

    auto ts_yyzz_xxxxx = pbuffer.data(idx_gh + 252);

    auto ts_yyzz_xxxxy = pbuffer.data(idx_gh + 253);

    auto ts_yyzz_xxxxz = pbuffer.data(idx_gh + 254);

    auto ts_yyzz_xxxyy = pbuffer.data(idx_gh + 255);

    auto ts_yyzz_xxxyz = pbuffer.data(idx_gh + 256);

    auto ts_yyzz_xxxzz = pbuffer.data(idx_gh + 257);

    auto ts_yyzz_xxyyy = pbuffer.data(idx_gh + 258);

    auto ts_yyzz_xxyyz = pbuffer.data(idx_gh + 259);

    auto ts_yyzz_xxyzz = pbuffer.data(idx_gh + 260);

    auto ts_yyzz_xxzzz = pbuffer.data(idx_gh + 261);

    auto ts_yyzz_xyyyy = pbuffer.data(idx_gh + 262);

    auto ts_yyzz_xyyyz = pbuffer.data(idx_gh + 263);

    auto ts_yyzz_xyyzz = pbuffer.data(idx_gh + 264);

    auto ts_yyzz_xyzzz = pbuffer.data(idx_gh + 265);

    auto ts_yyzz_xzzzz = pbuffer.data(idx_gh + 266);

    auto ts_yyzz_yyyyy = pbuffer.data(idx_gh + 267);

    auto ts_yyzz_yyyyz = pbuffer.data(idx_gh + 268);

    auto ts_yyzz_yyyzz = pbuffer.data(idx_gh + 269);

    auto ts_yyzz_yyzzz = pbuffer.data(idx_gh + 270);

    auto ts_yyzz_yzzzz = pbuffer.data(idx_gh + 271);

    auto ts_yyzz_zzzzz = pbuffer.data(idx_gh + 272);

    auto ts_yzzz_xxxxx = pbuffer.data(idx_gh + 273);

    auto ts_yzzz_xxxxy = pbuffer.data(idx_gh + 274);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_gh + 275);

    auto ts_yzzz_xxxyy = pbuffer.data(idx_gh + 276);

    auto ts_yzzz_xxxyz = pbuffer.data(idx_gh + 277);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_gh + 278);

    auto ts_yzzz_xxyyy = pbuffer.data(idx_gh + 279);

    auto ts_yzzz_xxyyz = pbuffer.data(idx_gh + 280);

    auto ts_yzzz_xxyzz = pbuffer.data(idx_gh + 281);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_gh + 282);

    auto ts_yzzz_xyyyy = pbuffer.data(idx_gh + 283);

    auto ts_yzzz_xyyyz = pbuffer.data(idx_gh + 284);

    auto ts_yzzz_xyyzz = pbuffer.data(idx_gh + 285);

    auto ts_yzzz_xyzzz = pbuffer.data(idx_gh + 286);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_gh + 287);

    auto ts_yzzz_yyyyy = pbuffer.data(idx_gh + 288);

    auto ts_yzzz_yyyyz = pbuffer.data(idx_gh + 289);

    auto ts_yzzz_yyyzz = pbuffer.data(idx_gh + 290);

    auto ts_yzzz_yyzzz = pbuffer.data(idx_gh + 291);

    auto ts_yzzz_yzzzz = pbuffer.data(idx_gh + 292);

    auto ts_yzzz_zzzzz = pbuffer.data(idx_gh + 293);

    auto ts_zzzz_xxxxx = pbuffer.data(idx_gh + 294);

    auto ts_zzzz_xxxxy = pbuffer.data(idx_gh + 295);

    auto ts_zzzz_xxxxz = pbuffer.data(idx_gh + 296);

    auto ts_zzzz_xxxyy = pbuffer.data(idx_gh + 297);

    auto ts_zzzz_xxxyz = pbuffer.data(idx_gh + 298);

    auto ts_zzzz_xxxzz = pbuffer.data(idx_gh + 299);

    auto ts_zzzz_xxyyy = pbuffer.data(idx_gh + 300);

    auto ts_zzzz_xxyyz = pbuffer.data(idx_gh + 301);

    auto ts_zzzz_xxyzz = pbuffer.data(idx_gh + 302);

    auto ts_zzzz_xxzzz = pbuffer.data(idx_gh + 303);

    auto ts_zzzz_xyyyy = pbuffer.data(idx_gh + 304);

    auto ts_zzzz_xyyyz = pbuffer.data(idx_gh + 305);

    auto ts_zzzz_xyyzz = pbuffer.data(idx_gh + 306);

    auto ts_zzzz_xyzzz = pbuffer.data(idx_gh + 307);

    auto ts_zzzz_xzzzz = pbuffer.data(idx_gh + 308);

    auto ts_zzzz_yyyyy = pbuffer.data(idx_gh + 309);

    auto ts_zzzz_yyyyz = pbuffer.data(idx_gh + 310);

    auto ts_zzzz_yyyzz = pbuffer.data(idx_gh + 311);

    auto ts_zzzz_yyzzz = pbuffer.data(idx_gh + 312);

    auto ts_zzzz_yzzzz = pbuffer.data(idx_gh + 313);

    auto ts_zzzz_zzzzz = pbuffer.data(idx_gh + 314);

    // Set up components of auxiliary buffer : HG

    auto ts_xxxxx_xxxx = pbuffer.data(idx_hg);

    auto ts_xxxxx_xxxy = pbuffer.data(idx_hg + 1);

    auto ts_xxxxx_xxxz = pbuffer.data(idx_hg + 2);

    auto ts_xxxxx_xxyy = pbuffer.data(idx_hg + 3);

    auto ts_xxxxx_xxyz = pbuffer.data(idx_hg + 4);

    auto ts_xxxxx_xxzz = pbuffer.data(idx_hg + 5);

    auto ts_xxxxx_xyyy = pbuffer.data(idx_hg + 6);

    auto ts_xxxxx_xyyz = pbuffer.data(idx_hg + 7);

    auto ts_xxxxx_xyzz = pbuffer.data(idx_hg + 8);

    auto ts_xxxxx_xzzz = pbuffer.data(idx_hg + 9);

    auto ts_xxxxx_yyyy = pbuffer.data(idx_hg + 10);

    auto ts_xxxxx_yyyz = pbuffer.data(idx_hg + 11);

    auto ts_xxxxx_yyzz = pbuffer.data(idx_hg + 12);

    auto ts_xxxxx_yzzz = pbuffer.data(idx_hg + 13);

    auto ts_xxxxx_zzzz = pbuffer.data(idx_hg + 14);

    auto ts_xxxxy_xxxx = pbuffer.data(idx_hg + 15);

    auto ts_xxxxy_xxxy = pbuffer.data(idx_hg + 16);

    auto ts_xxxxy_xxxz = pbuffer.data(idx_hg + 17);

    auto ts_xxxxy_xxyy = pbuffer.data(idx_hg + 18);

    auto ts_xxxxy_xxyz = pbuffer.data(idx_hg + 19);

    auto ts_xxxxy_xxzz = pbuffer.data(idx_hg + 20);

    auto ts_xxxxy_xyyy = pbuffer.data(idx_hg + 21);

    auto ts_xxxxy_xyyz = pbuffer.data(idx_hg + 22);

    auto ts_xxxxy_xyzz = pbuffer.data(idx_hg + 23);

    auto ts_xxxxy_xzzz = pbuffer.data(idx_hg + 24);

    auto ts_xxxxy_yyyy = pbuffer.data(idx_hg + 25);

    auto ts_xxxxy_yyyz = pbuffer.data(idx_hg + 26);

    auto ts_xxxxy_yyzz = pbuffer.data(idx_hg + 27);

    auto ts_xxxxy_yzzz = pbuffer.data(idx_hg + 28);

    auto ts_xxxxy_zzzz = pbuffer.data(idx_hg + 29);

    auto ts_xxxxz_xxxx = pbuffer.data(idx_hg + 30);

    auto ts_xxxxz_xxxy = pbuffer.data(idx_hg + 31);

    auto ts_xxxxz_xxxz = pbuffer.data(idx_hg + 32);

    auto ts_xxxxz_xxyy = pbuffer.data(idx_hg + 33);

    auto ts_xxxxz_xxyz = pbuffer.data(idx_hg + 34);

    auto ts_xxxxz_xxzz = pbuffer.data(idx_hg + 35);

    auto ts_xxxxz_xyyy = pbuffer.data(idx_hg + 36);

    auto ts_xxxxz_xyyz = pbuffer.data(idx_hg + 37);

    auto ts_xxxxz_xyzz = pbuffer.data(idx_hg + 38);

    auto ts_xxxxz_xzzz = pbuffer.data(idx_hg + 39);

    auto ts_xxxxz_yyyy = pbuffer.data(idx_hg + 40);

    auto ts_xxxxz_yyyz = pbuffer.data(idx_hg + 41);

    auto ts_xxxxz_yyzz = pbuffer.data(idx_hg + 42);

    auto ts_xxxxz_yzzz = pbuffer.data(idx_hg + 43);

    auto ts_xxxxz_zzzz = pbuffer.data(idx_hg + 44);

    auto ts_xxxyy_xxxx = pbuffer.data(idx_hg + 45);

    auto ts_xxxyy_xxxy = pbuffer.data(idx_hg + 46);

    auto ts_xxxyy_xxxz = pbuffer.data(idx_hg + 47);

    auto ts_xxxyy_xxyy = pbuffer.data(idx_hg + 48);

    auto ts_xxxyy_xxyz = pbuffer.data(idx_hg + 49);

    auto ts_xxxyy_xxzz = pbuffer.data(idx_hg + 50);

    auto ts_xxxyy_xyyy = pbuffer.data(idx_hg + 51);

    auto ts_xxxyy_xyyz = pbuffer.data(idx_hg + 52);

    auto ts_xxxyy_xyzz = pbuffer.data(idx_hg + 53);

    auto ts_xxxyy_xzzz = pbuffer.data(idx_hg + 54);

    auto ts_xxxyy_yyyy = pbuffer.data(idx_hg + 55);

    auto ts_xxxyy_yyyz = pbuffer.data(idx_hg + 56);

    auto ts_xxxyy_yyzz = pbuffer.data(idx_hg + 57);

    auto ts_xxxyy_yzzz = pbuffer.data(idx_hg + 58);

    auto ts_xxxyy_zzzz = pbuffer.data(idx_hg + 59);

    auto ts_xxxyz_xxxx = pbuffer.data(idx_hg + 60);

    auto ts_xxxyz_xxxy = pbuffer.data(idx_hg + 61);

    auto ts_xxxyz_xxxz = pbuffer.data(idx_hg + 62);

    auto ts_xxxyz_xxyy = pbuffer.data(idx_hg + 63);

    auto ts_xxxyz_xxyz = pbuffer.data(idx_hg + 64);

    auto ts_xxxyz_xxzz = pbuffer.data(idx_hg + 65);

    auto ts_xxxyz_xyyy = pbuffer.data(idx_hg + 66);

    auto ts_xxxyz_xyyz = pbuffer.data(idx_hg + 67);

    auto ts_xxxyz_xyzz = pbuffer.data(idx_hg + 68);

    auto ts_xxxyz_xzzz = pbuffer.data(idx_hg + 69);

    auto ts_xxxyz_yyyy = pbuffer.data(idx_hg + 70);

    auto ts_xxxyz_yyyz = pbuffer.data(idx_hg + 71);

    auto ts_xxxyz_yyzz = pbuffer.data(idx_hg + 72);

    auto ts_xxxyz_yzzz = pbuffer.data(idx_hg + 73);

    auto ts_xxxyz_zzzz = pbuffer.data(idx_hg + 74);

    auto ts_xxxzz_xxxx = pbuffer.data(idx_hg + 75);

    auto ts_xxxzz_xxxy = pbuffer.data(idx_hg + 76);

    auto ts_xxxzz_xxxz = pbuffer.data(idx_hg + 77);

    auto ts_xxxzz_xxyy = pbuffer.data(idx_hg + 78);

    auto ts_xxxzz_xxyz = pbuffer.data(idx_hg + 79);

    auto ts_xxxzz_xxzz = pbuffer.data(idx_hg + 80);

    auto ts_xxxzz_xyyy = pbuffer.data(idx_hg + 81);

    auto ts_xxxzz_xyyz = pbuffer.data(idx_hg + 82);

    auto ts_xxxzz_xyzz = pbuffer.data(idx_hg + 83);

    auto ts_xxxzz_xzzz = pbuffer.data(idx_hg + 84);

    auto ts_xxxzz_yyyy = pbuffer.data(idx_hg + 85);

    auto ts_xxxzz_yyyz = pbuffer.data(idx_hg + 86);

    auto ts_xxxzz_yyzz = pbuffer.data(idx_hg + 87);

    auto ts_xxxzz_yzzz = pbuffer.data(idx_hg + 88);

    auto ts_xxxzz_zzzz = pbuffer.data(idx_hg + 89);

    auto ts_xxyyy_xxxx = pbuffer.data(idx_hg + 90);

    auto ts_xxyyy_xxxy = pbuffer.data(idx_hg + 91);

    auto ts_xxyyy_xxxz = pbuffer.data(idx_hg + 92);

    auto ts_xxyyy_xxyy = pbuffer.data(idx_hg + 93);

    auto ts_xxyyy_xxyz = pbuffer.data(idx_hg + 94);

    auto ts_xxyyy_xxzz = pbuffer.data(idx_hg + 95);

    auto ts_xxyyy_xyyy = pbuffer.data(idx_hg + 96);

    auto ts_xxyyy_xyyz = pbuffer.data(idx_hg + 97);

    auto ts_xxyyy_xyzz = pbuffer.data(idx_hg + 98);

    auto ts_xxyyy_xzzz = pbuffer.data(idx_hg + 99);

    auto ts_xxyyy_yyyy = pbuffer.data(idx_hg + 100);

    auto ts_xxyyy_yyyz = pbuffer.data(idx_hg + 101);

    auto ts_xxyyy_yyzz = pbuffer.data(idx_hg + 102);

    auto ts_xxyyy_yzzz = pbuffer.data(idx_hg + 103);

    auto ts_xxyyy_zzzz = pbuffer.data(idx_hg + 104);

    auto ts_xxyyz_xxxx = pbuffer.data(idx_hg + 105);

    auto ts_xxyyz_xxxy = pbuffer.data(idx_hg + 106);

    auto ts_xxyyz_xxxz = pbuffer.data(idx_hg + 107);

    auto ts_xxyyz_xxyy = pbuffer.data(idx_hg + 108);

    auto ts_xxyyz_xxyz = pbuffer.data(idx_hg + 109);

    auto ts_xxyyz_xxzz = pbuffer.data(idx_hg + 110);

    auto ts_xxyyz_xyyy = pbuffer.data(idx_hg + 111);

    auto ts_xxyyz_xyyz = pbuffer.data(idx_hg + 112);

    auto ts_xxyyz_xyzz = pbuffer.data(idx_hg + 113);

    auto ts_xxyyz_xzzz = pbuffer.data(idx_hg + 114);

    auto ts_xxyyz_yyyy = pbuffer.data(idx_hg + 115);

    auto ts_xxyyz_yyyz = pbuffer.data(idx_hg + 116);

    auto ts_xxyyz_yyzz = pbuffer.data(idx_hg + 117);

    auto ts_xxyyz_yzzz = pbuffer.data(idx_hg + 118);

    auto ts_xxyyz_zzzz = pbuffer.data(idx_hg + 119);

    auto ts_xxyzz_xxxx = pbuffer.data(idx_hg + 120);

    auto ts_xxyzz_xxxy = pbuffer.data(idx_hg + 121);

    auto ts_xxyzz_xxxz = pbuffer.data(idx_hg + 122);

    auto ts_xxyzz_xxyy = pbuffer.data(idx_hg + 123);

    auto ts_xxyzz_xxyz = pbuffer.data(idx_hg + 124);

    auto ts_xxyzz_xxzz = pbuffer.data(idx_hg + 125);

    auto ts_xxyzz_xyyy = pbuffer.data(idx_hg + 126);

    auto ts_xxyzz_xyyz = pbuffer.data(idx_hg + 127);

    auto ts_xxyzz_xyzz = pbuffer.data(idx_hg + 128);

    auto ts_xxyzz_xzzz = pbuffer.data(idx_hg + 129);

    auto ts_xxyzz_yyyy = pbuffer.data(idx_hg + 130);

    auto ts_xxyzz_yyyz = pbuffer.data(idx_hg + 131);

    auto ts_xxyzz_yyzz = pbuffer.data(idx_hg + 132);

    auto ts_xxyzz_yzzz = pbuffer.data(idx_hg + 133);

    auto ts_xxyzz_zzzz = pbuffer.data(idx_hg + 134);

    auto ts_xxzzz_xxxx = pbuffer.data(idx_hg + 135);

    auto ts_xxzzz_xxxy = pbuffer.data(idx_hg + 136);

    auto ts_xxzzz_xxxz = pbuffer.data(idx_hg + 137);

    auto ts_xxzzz_xxyy = pbuffer.data(idx_hg + 138);

    auto ts_xxzzz_xxyz = pbuffer.data(idx_hg + 139);

    auto ts_xxzzz_xxzz = pbuffer.data(idx_hg + 140);

    auto ts_xxzzz_xyyy = pbuffer.data(idx_hg + 141);

    auto ts_xxzzz_xyyz = pbuffer.data(idx_hg + 142);

    auto ts_xxzzz_xyzz = pbuffer.data(idx_hg + 143);

    auto ts_xxzzz_xzzz = pbuffer.data(idx_hg + 144);

    auto ts_xxzzz_yyyy = pbuffer.data(idx_hg + 145);

    auto ts_xxzzz_yyyz = pbuffer.data(idx_hg + 146);

    auto ts_xxzzz_yyzz = pbuffer.data(idx_hg + 147);

    auto ts_xxzzz_yzzz = pbuffer.data(idx_hg + 148);

    auto ts_xxzzz_zzzz = pbuffer.data(idx_hg + 149);

    auto ts_xyyyy_xxxx = pbuffer.data(idx_hg + 150);

    auto ts_xyyyy_xxxy = pbuffer.data(idx_hg + 151);

    auto ts_xyyyy_xxxz = pbuffer.data(idx_hg + 152);

    auto ts_xyyyy_xxyy = pbuffer.data(idx_hg + 153);

    auto ts_xyyyy_xxyz = pbuffer.data(idx_hg + 154);

    auto ts_xyyyy_xxzz = pbuffer.data(idx_hg + 155);

    auto ts_xyyyy_xyyy = pbuffer.data(idx_hg + 156);

    auto ts_xyyyy_xyyz = pbuffer.data(idx_hg + 157);

    auto ts_xyyyy_xyzz = pbuffer.data(idx_hg + 158);

    auto ts_xyyyy_xzzz = pbuffer.data(idx_hg + 159);

    auto ts_xyyyy_yyyy = pbuffer.data(idx_hg + 160);

    auto ts_xyyyy_yyyz = pbuffer.data(idx_hg + 161);

    auto ts_xyyyy_yyzz = pbuffer.data(idx_hg + 162);

    auto ts_xyyyy_yzzz = pbuffer.data(idx_hg + 163);

    auto ts_xyyyy_zzzz = pbuffer.data(idx_hg + 164);

    auto ts_xyyyz_xxxx = pbuffer.data(idx_hg + 165);

    auto ts_xyyyz_xxxy = pbuffer.data(idx_hg + 166);

    auto ts_xyyyz_xxxz = pbuffer.data(idx_hg + 167);

    auto ts_xyyyz_xxyy = pbuffer.data(idx_hg + 168);

    auto ts_xyyyz_xxyz = pbuffer.data(idx_hg + 169);

    auto ts_xyyyz_xxzz = pbuffer.data(idx_hg + 170);

    auto ts_xyyyz_xyyy = pbuffer.data(idx_hg + 171);

    auto ts_xyyyz_xyyz = pbuffer.data(idx_hg + 172);

    auto ts_xyyyz_xyzz = pbuffer.data(idx_hg + 173);

    auto ts_xyyyz_xzzz = pbuffer.data(idx_hg + 174);

    auto ts_xyyyz_yyyy = pbuffer.data(idx_hg + 175);

    auto ts_xyyyz_yyyz = pbuffer.data(idx_hg + 176);

    auto ts_xyyyz_yyzz = pbuffer.data(idx_hg + 177);

    auto ts_xyyyz_yzzz = pbuffer.data(idx_hg + 178);

    auto ts_xyyyz_zzzz = pbuffer.data(idx_hg + 179);

    auto ts_xyyzz_xxxx = pbuffer.data(idx_hg + 180);

    auto ts_xyyzz_xxxy = pbuffer.data(idx_hg + 181);

    auto ts_xyyzz_xxxz = pbuffer.data(idx_hg + 182);

    auto ts_xyyzz_xxyy = pbuffer.data(idx_hg + 183);

    auto ts_xyyzz_xxyz = pbuffer.data(idx_hg + 184);

    auto ts_xyyzz_xxzz = pbuffer.data(idx_hg + 185);

    auto ts_xyyzz_xyyy = pbuffer.data(idx_hg + 186);

    auto ts_xyyzz_xyyz = pbuffer.data(idx_hg + 187);

    auto ts_xyyzz_xyzz = pbuffer.data(idx_hg + 188);

    auto ts_xyyzz_xzzz = pbuffer.data(idx_hg + 189);

    auto ts_xyyzz_yyyy = pbuffer.data(idx_hg + 190);

    auto ts_xyyzz_yyyz = pbuffer.data(idx_hg + 191);

    auto ts_xyyzz_yyzz = pbuffer.data(idx_hg + 192);

    auto ts_xyyzz_yzzz = pbuffer.data(idx_hg + 193);

    auto ts_xyyzz_zzzz = pbuffer.data(idx_hg + 194);

    auto ts_xyzzz_xxxx = pbuffer.data(idx_hg + 195);

    auto ts_xyzzz_xxxy = pbuffer.data(idx_hg + 196);

    auto ts_xyzzz_xxxz = pbuffer.data(idx_hg + 197);

    auto ts_xyzzz_xxyy = pbuffer.data(idx_hg + 198);

    auto ts_xyzzz_xxyz = pbuffer.data(idx_hg + 199);

    auto ts_xyzzz_xxzz = pbuffer.data(idx_hg + 200);

    auto ts_xyzzz_xyyy = pbuffer.data(idx_hg + 201);

    auto ts_xyzzz_xyyz = pbuffer.data(idx_hg + 202);

    auto ts_xyzzz_xyzz = pbuffer.data(idx_hg + 203);

    auto ts_xyzzz_xzzz = pbuffer.data(idx_hg + 204);

    auto ts_xyzzz_yyyy = pbuffer.data(idx_hg + 205);

    auto ts_xyzzz_yyyz = pbuffer.data(idx_hg + 206);

    auto ts_xyzzz_yyzz = pbuffer.data(idx_hg + 207);

    auto ts_xyzzz_yzzz = pbuffer.data(idx_hg + 208);

    auto ts_xyzzz_zzzz = pbuffer.data(idx_hg + 209);

    auto ts_xzzzz_xxxx = pbuffer.data(idx_hg + 210);

    auto ts_xzzzz_xxxy = pbuffer.data(idx_hg + 211);

    auto ts_xzzzz_xxxz = pbuffer.data(idx_hg + 212);

    auto ts_xzzzz_xxyy = pbuffer.data(idx_hg + 213);

    auto ts_xzzzz_xxyz = pbuffer.data(idx_hg + 214);

    auto ts_xzzzz_xxzz = pbuffer.data(idx_hg + 215);

    auto ts_xzzzz_xyyy = pbuffer.data(idx_hg + 216);

    auto ts_xzzzz_xyyz = pbuffer.data(idx_hg + 217);

    auto ts_xzzzz_xyzz = pbuffer.data(idx_hg + 218);

    auto ts_xzzzz_xzzz = pbuffer.data(idx_hg + 219);

    auto ts_xzzzz_yyyy = pbuffer.data(idx_hg + 220);

    auto ts_xzzzz_yyyz = pbuffer.data(idx_hg + 221);

    auto ts_xzzzz_yyzz = pbuffer.data(idx_hg + 222);

    auto ts_xzzzz_yzzz = pbuffer.data(idx_hg + 223);

    auto ts_xzzzz_zzzz = pbuffer.data(idx_hg + 224);

    auto ts_yyyyy_xxxx = pbuffer.data(idx_hg + 225);

    auto ts_yyyyy_xxxy = pbuffer.data(idx_hg + 226);

    auto ts_yyyyy_xxxz = pbuffer.data(idx_hg + 227);

    auto ts_yyyyy_xxyy = pbuffer.data(idx_hg + 228);

    auto ts_yyyyy_xxyz = pbuffer.data(idx_hg + 229);

    auto ts_yyyyy_xxzz = pbuffer.data(idx_hg + 230);

    auto ts_yyyyy_xyyy = pbuffer.data(idx_hg + 231);

    auto ts_yyyyy_xyyz = pbuffer.data(idx_hg + 232);

    auto ts_yyyyy_xyzz = pbuffer.data(idx_hg + 233);

    auto ts_yyyyy_xzzz = pbuffer.data(idx_hg + 234);

    auto ts_yyyyy_yyyy = pbuffer.data(idx_hg + 235);

    auto ts_yyyyy_yyyz = pbuffer.data(idx_hg + 236);

    auto ts_yyyyy_yyzz = pbuffer.data(idx_hg + 237);

    auto ts_yyyyy_yzzz = pbuffer.data(idx_hg + 238);

    auto ts_yyyyy_zzzz = pbuffer.data(idx_hg + 239);

    auto ts_yyyyz_xxxx = pbuffer.data(idx_hg + 240);

    auto ts_yyyyz_xxxy = pbuffer.data(idx_hg + 241);

    auto ts_yyyyz_xxxz = pbuffer.data(idx_hg + 242);

    auto ts_yyyyz_xxyy = pbuffer.data(idx_hg + 243);

    auto ts_yyyyz_xxyz = pbuffer.data(idx_hg + 244);

    auto ts_yyyyz_xxzz = pbuffer.data(idx_hg + 245);

    auto ts_yyyyz_xyyy = pbuffer.data(idx_hg + 246);

    auto ts_yyyyz_xyyz = pbuffer.data(idx_hg + 247);

    auto ts_yyyyz_xyzz = pbuffer.data(idx_hg + 248);

    auto ts_yyyyz_xzzz = pbuffer.data(idx_hg + 249);

    auto ts_yyyyz_yyyy = pbuffer.data(idx_hg + 250);

    auto ts_yyyyz_yyyz = pbuffer.data(idx_hg + 251);

    auto ts_yyyyz_yyzz = pbuffer.data(idx_hg + 252);

    auto ts_yyyyz_yzzz = pbuffer.data(idx_hg + 253);

    auto ts_yyyyz_zzzz = pbuffer.data(idx_hg + 254);

    auto ts_yyyzz_xxxx = pbuffer.data(idx_hg + 255);

    auto ts_yyyzz_xxxy = pbuffer.data(idx_hg + 256);

    auto ts_yyyzz_xxxz = pbuffer.data(idx_hg + 257);

    auto ts_yyyzz_xxyy = pbuffer.data(idx_hg + 258);

    auto ts_yyyzz_xxyz = pbuffer.data(idx_hg + 259);

    auto ts_yyyzz_xxzz = pbuffer.data(idx_hg + 260);

    auto ts_yyyzz_xyyy = pbuffer.data(idx_hg + 261);

    auto ts_yyyzz_xyyz = pbuffer.data(idx_hg + 262);

    auto ts_yyyzz_xyzz = pbuffer.data(idx_hg + 263);

    auto ts_yyyzz_xzzz = pbuffer.data(idx_hg + 264);

    auto ts_yyyzz_yyyy = pbuffer.data(idx_hg + 265);

    auto ts_yyyzz_yyyz = pbuffer.data(idx_hg + 266);

    auto ts_yyyzz_yyzz = pbuffer.data(idx_hg + 267);

    auto ts_yyyzz_yzzz = pbuffer.data(idx_hg + 268);

    auto ts_yyyzz_zzzz = pbuffer.data(idx_hg + 269);

    auto ts_yyzzz_xxxx = pbuffer.data(idx_hg + 270);

    auto ts_yyzzz_xxxy = pbuffer.data(idx_hg + 271);

    auto ts_yyzzz_xxxz = pbuffer.data(idx_hg + 272);

    auto ts_yyzzz_xxyy = pbuffer.data(idx_hg + 273);

    auto ts_yyzzz_xxyz = pbuffer.data(idx_hg + 274);

    auto ts_yyzzz_xxzz = pbuffer.data(idx_hg + 275);

    auto ts_yyzzz_xyyy = pbuffer.data(idx_hg + 276);

    auto ts_yyzzz_xyyz = pbuffer.data(idx_hg + 277);

    auto ts_yyzzz_xyzz = pbuffer.data(idx_hg + 278);

    auto ts_yyzzz_xzzz = pbuffer.data(idx_hg + 279);

    auto ts_yyzzz_yyyy = pbuffer.data(idx_hg + 280);

    auto ts_yyzzz_yyyz = pbuffer.data(idx_hg + 281);

    auto ts_yyzzz_yyzz = pbuffer.data(idx_hg + 282);

    auto ts_yyzzz_yzzz = pbuffer.data(idx_hg + 283);

    auto ts_yyzzz_zzzz = pbuffer.data(idx_hg + 284);

    auto ts_yzzzz_xxxx = pbuffer.data(idx_hg + 285);

    auto ts_yzzzz_xxxy = pbuffer.data(idx_hg + 286);

    auto ts_yzzzz_xxxz = pbuffer.data(idx_hg + 287);

    auto ts_yzzzz_xxyy = pbuffer.data(idx_hg + 288);

    auto ts_yzzzz_xxyz = pbuffer.data(idx_hg + 289);

    auto ts_yzzzz_xxzz = pbuffer.data(idx_hg + 290);

    auto ts_yzzzz_xyyy = pbuffer.data(idx_hg + 291);

    auto ts_yzzzz_xyyz = pbuffer.data(idx_hg + 292);

    auto ts_yzzzz_xyzz = pbuffer.data(idx_hg + 293);

    auto ts_yzzzz_xzzz = pbuffer.data(idx_hg + 294);

    auto ts_yzzzz_yyyy = pbuffer.data(idx_hg + 295);

    auto ts_yzzzz_yyyz = pbuffer.data(idx_hg + 296);

    auto ts_yzzzz_yyzz = pbuffer.data(idx_hg + 297);

    auto ts_yzzzz_yzzz = pbuffer.data(idx_hg + 298);

    auto ts_yzzzz_zzzz = pbuffer.data(idx_hg + 299);

    auto ts_zzzzz_xxxx = pbuffer.data(idx_hg + 300);

    auto ts_zzzzz_xxxy = pbuffer.data(idx_hg + 301);

    auto ts_zzzzz_xxxz = pbuffer.data(idx_hg + 302);

    auto ts_zzzzz_xxyy = pbuffer.data(idx_hg + 303);

    auto ts_zzzzz_xxyz = pbuffer.data(idx_hg + 304);

    auto ts_zzzzz_xxzz = pbuffer.data(idx_hg + 305);

    auto ts_zzzzz_xyyy = pbuffer.data(idx_hg + 306);

    auto ts_zzzzz_xyyz = pbuffer.data(idx_hg + 307);

    auto ts_zzzzz_xyzz = pbuffer.data(idx_hg + 308);

    auto ts_zzzzz_xzzz = pbuffer.data(idx_hg + 309);

    auto ts_zzzzz_yyyy = pbuffer.data(idx_hg + 310);

    auto ts_zzzzz_yyyz = pbuffer.data(idx_hg + 311);

    auto ts_zzzzz_yyzz = pbuffer.data(idx_hg + 312);

    auto ts_zzzzz_yzzz = pbuffer.data(idx_hg + 313);

    auto ts_zzzzz_zzzz = pbuffer.data(idx_hg + 314);

    // Set up components of auxiliary buffer : HH

    auto ts_xxxxx_xxxxx = pbuffer.data(idx_hh);

    auto ts_xxxxx_xxxxy = pbuffer.data(idx_hh + 1);

    auto ts_xxxxx_xxxxz = pbuffer.data(idx_hh + 2);

    auto ts_xxxxx_xxxyy = pbuffer.data(idx_hh + 3);

    auto ts_xxxxx_xxxyz = pbuffer.data(idx_hh + 4);

    auto ts_xxxxx_xxxzz = pbuffer.data(idx_hh + 5);

    auto ts_xxxxx_xxyyy = pbuffer.data(idx_hh + 6);

    auto ts_xxxxx_xxyyz = pbuffer.data(idx_hh + 7);

    auto ts_xxxxx_xxyzz = pbuffer.data(idx_hh + 8);

    auto ts_xxxxx_xxzzz = pbuffer.data(idx_hh + 9);

    auto ts_xxxxx_xyyyy = pbuffer.data(idx_hh + 10);

    auto ts_xxxxx_xyyyz = pbuffer.data(idx_hh + 11);

    auto ts_xxxxx_xyyzz = pbuffer.data(idx_hh + 12);

    auto ts_xxxxx_xyzzz = pbuffer.data(idx_hh + 13);

    auto ts_xxxxx_xzzzz = pbuffer.data(idx_hh + 14);

    auto ts_xxxxx_yyyyy = pbuffer.data(idx_hh + 15);

    auto ts_xxxxx_yyyyz = pbuffer.data(idx_hh + 16);

    auto ts_xxxxx_yyyzz = pbuffer.data(idx_hh + 17);

    auto ts_xxxxx_yyzzz = pbuffer.data(idx_hh + 18);

    auto ts_xxxxx_yzzzz = pbuffer.data(idx_hh + 19);

    auto ts_xxxxx_zzzzz = pbuffer.data(idx_hh + 20);

    auto ts_xxxxy_xxxxx = pbuffer.data(idx_hh + 21);

    auto ts_xxxxy_xxxxy = pbuffer.data(idx_hh + 22);

    auto ts_xxxxy_xxxxz = pbuffer.data(idx_hh + 23);

    auto ts_xxxxy_xxxyy = pbuffer.data(idx_hh + 24);

    auto ts_xxxxy_xxxyz = pbuffer.data(idx_hh + 25);

    auto ts_xxxxy_xxxzz = pbuffer.data(idx_hh + 26);

    auto ts_xxxxy_xxyyy = pbuffer.data(idx_hh + 27);

    auto ts_xxxxy_xxyyz = pbuffer.data(idx_hh + 28);

    auto ts_xxxxy_xxyzz = pbuffer.data(idx_hh + 29);

    auto ts_xxxxy_xxzzz = pbuffer.data(idx_hh + 30);

    auto ts_xxxxy_xyyyy = pbuffer.data(idx_hh + 31);

    auto ts_xxxxy_xyyyz = pbuffer.data(idx_hh + 32);

    auto ts_xxxxy_xyyzz = pbuffer.data(idx_hh + 33);

    auto ts_xxxxy_xyzzz = pbuffer.data(idx_hh + 34);

    auto ts_xxxxy_xzzzz = pbuffer.data(idx_hh + 35);

    auto ts_xxxxy_yyyyy = pbuffer.data(idx_hh + 36);

    auto ts_xxxxy_yyyyz = pbuffer.data(idx_hh + 37);

    auto ts_xxxxy_yyyzz = pbuffer.data(idx_hh + 38);

    auto ts_xxxxy_yyzzz = pbuffer.data(idx_hh + 39);

    auto ts_xxxxy_yzzzz = pbuffer.data(idx_hh + 40);

    auto ts_xxxxy_zzzzz = pbuffer.data(idx_hh + 41);

    auto ts_xxxxz_xxxxx = pbuffer.data(idx_hh + 42);

    auto ts_xxxxz_xxxxy = pbuffer.data(idx_hh + 43);

    auto ts_xxxxz_xxxxz = pbuffer.data(idx_hh + 44);

    auto ts_xxxxz_xxxyy = pbuffer.data(idx_hh + 45);

    auto ts_xxxxz_xxxyz = pbuffer.data(idx_hh + 46);

    auto ts_xxxxz_xxxzz = pbuffer.data(idx_hh + 47);

    auto ts_xxxxz_xxyyy = pbuffer.data(idx_hh + 48);

    auto ts_xxxxz_xxyyz = pbuffer.data(idx_hh + 49);

    auto ts_xxxxz_xxyzz = pbuffer.data(idx_hh + 50);

    auto ts_xxxxz_xxzzz = pbuffer.data(idx_hh + 51);

    auto ts_xxxxz_xyyyy = pbuffer.data(idx_hh + 52);

    auto ts_xxxxz_xyyyz = pbuffer.data(idx_hh + 53);

    auto ts_xxxxz_xyyzz = pbuffer.data(idx_hh + 54);

    auto ts_xxxxz_xyzzz = pbuffer.data(idx_hh + 55);

    auto ts_xxxxz_xzzzz = pbuffer.data(idx_hh + 56);

    auto ts_xxxxz_yyyyy = pbuffer.data(idx_hh + 57);

    auto ts_xxxxz_yyyyz = pbuffer.data(idx_hh + 58);

    auto ts_xxxxz_yyyzz = pbuffer.data(idx_hh + 59);

    auto ts_xxxxz_yyzzz = pbuffer.data(idx_hh + 60);

    auto ts_xxxxz_yzzzz = pbuffer.data(idx_hh + 61);

    auto ts_xxxxz_zzzzz = pbuffer.data(idx_hh + 62);

    auto ts_xxxyy_xxxxx = pbuffer.data(idx_hh + 63);

    auto ts_xxxyy_xxxxy = pbuffer.data(idx_hh + 64);

    auto ts_xxxyy_xxxxz = pbuffer.data(idx_hh + 65);

    auto ts_xxxyy_xxxyy = pbuffer.data(idx_hh + 66);

    auto ts_xxxyy_xxxyz = pbuffer.data(idx_hh + 67);

    auto ts_xxxyy_xxxzz = pbuffer.data(idx_hh + 68);

    auto ts_xxxyy_xxyyy = pbuffer.data(idx_hh + 69);

    auto ts_xxxyy_xxyyz = pbuffer.data(idx_hh + 70);

    auto ts_xxxyy_xxyzz = pbuffer.data(idx_hh + 71);

    auto ts_xxxyy_xxzzz = pbuffer.data(idx_hh + 72);

    auto ts_xxxyy_xyyyy = pbuffer.data(idx_hh + 73);

    auto ts_xxxyy_xyyyz = pbuffer.data(idx_hh + 74);

    auto ts_xxxyy_xyyzz = pbuffer.data(idx_hh + 75);

    auto ts_xxxyy_xyzzz = pbuffer.data(idx_hh + 76);

    auto ts_xxxyy_xzzzz = pbuffer.data(idx_hh + 77);

    auto ts_xxxyy_yyyyy = pbuffer.data(idx_hh + 78);

    auto ts_xxxyy_yyyyz = pbuffer.data(idx_hh + 79);

    auto ts_xxxyy_yyyzz = pbuffer.data(idx_hh + 80);

    auto ts_xxxyy_yyzzz = pbuffer.data(idx_hh + 81);

    auto ts_xxxyy_yzzzz = pbuffer.data(idx_hh + 82);

    auto ts_xxxyy_zzzzz = pbuffer.data(idx_hh + 83);

    auto ts_xxxyz_xxxxx = pbuffer.data(idx_hh + 84);

    auto ts_xxxyz_xxxxy = pbuffer.data(idx_hh + 85);

    auto ts_xxxyz_xxxxz = pbuffer.data(idx_hh + 86);

    auto ts_xxxyz_xxxyy = pbuffer.data(idx_hh + 87);

    auto ts_xxxyz_xxxyz = pbuffer.data(idx_hh + 88);

    auto ts_xxxyz_xxxzz = pbuffer.data(idx_hh + 89);

    auto ts_xxxyz_xxyyy = pbuffer.data(idx_hh + 90);

    auto ts_xxxyz_xxyyz = pbuffer.data(idx_hh + 91);

    auto ts_xxxyz_xxyzz = pbuffer.data(idx_hh + 92);

    auto ts_xxxyz_xxzzz = pbuffer.data(idx_hh + 93);

    auto ts_xxxyz_xyyyy = pbuffer.data(idx_hh + 94);

    auto ts_xxxyz_xyyyz = pbuffer.data(idx_hh + 95);

    auto ts_xxxyz_xyyzz = pbuffer.data(idx_hh + 96);

    auto ts_xxxyz_xyzzz = pbuffer.data(idx_hh + 97);

    auto ts_xxxyz_xzzzz = pbuffer.data(idx_hh + 98);

    auto ts_xxxyz_yyyyy = pbuffer.data(idx_hh + 99);

    auto ts_xxxyz_yyyyz = pbuffer.data(idx_hh + 100);

    auto ts_xxxyz_yyyzz = pbuffer.data(idx_hh + 101);

    auto ts_xxxyz_yyzzz = pbuffer.data(idx_hh + 102);

    auto ts_xxxyz_yzzzz = pbuffer.data(idx_hh + 103);

    auto ts_xxxyz_zzzzz = pbuffer.data(idx_hh + 104);

    auto ts_xxxzz_xxxxx = pbuffer.data(idx_hh + 105);

    auto ts_xxxzz_xxxxy = pbuffer.data(idx_hh + 106);

    auto ts_xxxzz_xxxxz = pbuffer.data(idx_hh + 107);

    auto ts_xxxzz_xxxyy = pbuffer.data(idx_hh + 108);

    auto ts_xxxzz_xxxyz = pbuffer.data(idx_hh + 109);

    auto ts_xxxzz_xxxzz = pbuffer.data(idx_hh + 110);

    auto ts_xxxzz_xxyyy = pbuffer.data(idx_hh + 111);

    auto ts_xxxzz_xxyyz = pbuffer.data(idx_hh + 112);

    auto ts_xxxzz_xxyzz = pbuffer.data(idx_hh + 113);

    auto ts_xxxzz_xxzzz = pbuffer.data(idx_hh + 114);

    auto ts_xxxzz_xyyyy = pbuffer.data(idx_hh + 115);

    auto ts_xxxzz_xyyyz = pbuffer.data(idx_hh + 116);

    auto ts_xxxzz_xyyzz = pbuffer.data(idx_hh + 117);

    auto ts_xxxzz_xyzzz = pbuffer.data(idx_hh + 118);

    auto ts_xxxzz_xzzzz = pbuffer.data(idx_hh + 119);

    auto ts_xxxzz_yyyyy = pbuffer.data(idx_hh + 120);

    auto ts_xxxzz_yyyyz = pbuffer.data(idx_hh + 121);

    auto ts_xxxzz_yyyzz = pbuffer.data(idx_hh + 122);

    auto ts_xxxzz_yyzzz = pbuffer.data(idx_hh + 123);

    auto ts_xxxzz_yzzzz = pbuffer.data(idx_hh + 124);

    auto ts_xxxzz_zzzzz = pbuffer.data(idx_hh + 125);

    auto ts_xxyyy_xxxxx = pbuffer.data(idx_hh + 126);

    auto ts_xxyyy_xxxxy = pbuffer.data(idx_hh + 127);

    auto ts_xxyyy_xxxxz = pbuffer.data(idx_hh + 128);

    auto ts_xxyyy_xxxyy = pbuffer.data(idx_hh + 129);

    auto ts_xxyyy_xxxyz = pbuffer.data(idx_hh + 130);

    auto ts_xxyyy_xxxzz = pbuffer.data(idx_hh + 131);

    auto ts_xxyyy_xxyyy = pbuffer.data(idx_hh + 132);

    auto ts_xxyyy_xxyyz = pbuffer.data(idx_hh + 133);

    auto ts_xxyyy_xxyzz = pbuffer.data(idx_hh + 134);

    auto ts_xxyyy_xxzzz = pbuffer.data(idx_hh + 135);

    auto ts_xxyyy_xyyyy = pbuffer.data(idx_hh + 136);

    auto ts_xxyyy_xyyyz = pbuffer.data(idx_hh + 137);

    auto ts_xxyyy_xyyzz = pbuffer.data(idx_hh + 138);

    auto ts_xxyyy_xyzzz = pbuffer.data(idx_hh + 139);

    auto ts_xxyyy_xzzzz = pbuffer.data(idx_hh + 140);

    auto ts_xxyyy_yyyyy = pbuffer.data(idx_hh + 141);

    auto ts_xxyyy_yyyyz = pbuffer.data(idx_hh + 142);

    auto ts_xxyyy_yyyzz = pbuffer.data(idx_hh + 143);

    auto ts_xxyyy_yyzzz = pbuffer.data(idx_hh + 144);

    auto ts_xxyyy_yzzzz = pbuffer.data(idx_hh + 145);

    auto ts_xxyyy_zzzzz = pbuffer.data(idx_hh + 146);

    auto ts_xxyyz_xxxxx = pbuffer.data(idx_hh + 147);

    auto ts_xxyyz_xxxxy = pbuffer.data(idx_hh + 148);

    auto ts_xxyyz_xxxxz = pbuffer.data(idx_hh + 149);

    auto ts_xxyyz_xxxyy = pbuffer.data(idx_hh + 150);

    auto ts_xxyyz_xxxyz = pbuffer.data(idx_hh + 151);

    auto ts_xxyyz_xxxzz = pbuffer.data(idx_hh + 152);

    auto ts_xxyyz_xxyyy = pbuffer.data(idx_hh + 153);

    auto ts_xxyyz_xxyyz = pbuffer.data(idx_hh + 154);

    auto ts_xxyyz_xxyzz = pbuffer.data(idx_hh + 155);

    auto ts_xxyyz_xxzzz = pbuffer.data(idx_hh + 156);

    auto ts_xxyyz_xyyyy = pbuffer.data(idx_hh + 157);

    auto ts_xxyyz_xyyyz = pbuffer.data(idx_hh + 158);

    auto ts_xxyyz_xyyzz = pbuffer.data(idx_hh + 159);

    auto ts_xxyyz_xyzzz = pbuffer.data(idx_hh + 160);

    auto ts_xxyyz_xzzzz = pbuffer.data(idx_hh + 161);

    auto ts_xxyyz_yyyyy = pbuffer.data(idx_hh + 162);

    auto ts_xxyyz_yyyyz = pbuffer.data(idx_hh + 163);

    auto ts_xxyyz_yyyzz = pbuffer.data(idx_hh + 164);

    auto ts_xxyyz_yyzzz = pbuffer.data(idx_hh + 165);

    auto ts_xxyyz_yzzzz = pbuffer.data(idx_hh + 166);

    auto ts_xxyyz_zzzzz = pbuffer.data(idx_hh + 167);

    auto ts_xxyzz_xxxxx = pbuffer.data(idx_hh + 168);

    auto ts_xxyzz_xxxxy = pbuffer.data(idx_hh + 169);

    auto ts_xxyzz_xxxxz = pbuffer.data(idx_hh + 170);

    auto ts_xxyzz_xxxyy = pbuffer.data(idx_hh + 171);

    auto ts_xxyzz_xxxyz = pbuffer.data(idx_hh + 172);

    auto ts_xxyzz_xxxzz = pbuffer.data(idx_hh + 173);

    auto ts_xxyzz_xxyyy = pbuffer.data(idx_hh + 174);

    auto ts_xxyzz_xxyyz = pbuffer.data(idx_hh + 175);

    auto ts_xxyzz_xxyzz = pbuffer.data(idx_hh + 176);

    auto ts_xxyzz_xxzzz = pbuffer.data(idx_hh + 177);

    auto ts_xxyzz_xyyyy = pbuffer.data(idx_hh + 178);

    auto ts_xxyzz_xyyyz = pbuffer.data(idx_hh + 179);

    auto ts_xxyzz_xyyzz = pbuffer.data(idx_hh + 180);

    auto ts_xxyzz_xyzzz = pbuffer.data(idx_hh + 181);

    auto ts_xxyzz_xzzzz = pbuffer.data(idx_hh + 182);

    auto ts_xxyzz_yyyyy = pbuffer.data(idx_hh + 183);

    auto ts_xxyzz_yyyyz = pbuffer.data(idx_hh + 184);

    auto ts_xxyzz_yyyzz = pbuffer.data(idx_hh + 185);

    auto ts_xxyzz_yyzzz = pbuffer.data(idx_hh + 186);

    auto ts_xxyzz_yzzzz = pbuffer.data(idx_hh + 187);

    auto ts_xxyzz_zzzzz = pbuffer.data(idx_hh + 188);

    auto ts_xxzzz_xxxxx = pbuffer.data(idx_hh + 189);

    auto ts_xxzzz_xxxxy = pbuffer.data(idx_hh + 190);

    auto ts_xxzzz_xxxxz = pbuffer.data(idx_hh + 191);

    auto ts_xxzzz_xxxyy = pbuffer.data(idx_hh + 192);

    auto ts_xxzzz_xxxyz = pbuffer.data(idx_hh + 193);

    auto ts_xxzzz_xxxzz = pbuffer.data(idx_hh + 194);

    auto ts_xxzzz_xxyyy = pbuffer.data(idx_hh + 195);

    auto ts_xxzzz_xxyyz = pbuffer.data(idx_hh + 196);

    auto ts_xxzzz_xxyzz = pbuffer.data(idx_hh + 197);

    auto ts_xxzzz_xxzzz = pbuffer.data(idx_hh + 198);

    auto ts_xxzzz_xyyyy = pbuffer.data(idx_hh + 199);

    auto ts_xxzzz_xyyyz = pbuffer.data(idx_hh + 200);

    auto ts_xxzzz_xyyzz = pbuffer.data(idx_hh + 201);

    auto ts_xxzzz_xyzzz = pbuffer.data(idx_hh + 202);

    auto ts_xxzzz_xzzzz = pbuffer.data(idx_hh + 203);

    auto ts_xxzzz_yyyyy = pbuffer.data(idx_hh + 204);

    auto ts_xxzzz_yyyyz = pbuffer.data(idx_hh + 205);

    auto ts_xxzzz_yyyzz = pbuffer.data(idx_hh + 206);

    auto ts_xxzzz_yyzzz = pbuffer.data(idx_hh + 207);

    auto ts_xxzzz_yzzzz = pbuffer.data(idx_hh + 208);

    auto ts_xxzzz_zzzzz = pbuffer.data(idx_hh + 209);

    auto ts_xyyyy_xxxxx = pbuffer.data(idx_hh + 210);

    auto ts_xyyyy_xxxxy = pbuffer.data(idx_hh + 211);

    auto ts_xyyyy_xxxxz = pbuffer.data(idx_hh + 212);

    auto ts_xyyyy_xxxyy = pbuffer.data(idx_hh + 213);

    auto ts_xyyyy_xxxyz = pbuffer.data(idx_hh + 214);

    auto ts_xyyyy_xxxzz = pbuffer.data(idx_hh + 215);

    auto ts_xyyyy_xxyyy = pbuffer.data(idx_hh + 216);

    auto ts_xyyyy_xxyyz = pbuffer.data(idx_hh + 217);

    auto ts_xyyyy_xxyzz = pbuffer.data(idx_hh + 218);

    auto ts_xyyyy_xxzzz = pbuffer.data(idx_hh + 219);

    auto ts_xyyyy_xyyyy = pbuffer.data(idx_hh + 220);

    auto ts_xyyyy_xyyyz = pbuffer.data(idx_hh + 221);

    auto ts_xyyyy_xyyzz = pbuffer.data(idx_hh + 222);

    auto ts_xyyyy_xyzzz = pbuffer.data(idx_hh + 223);

    auto ts_xyyyy_xzzzz = pbuffer.data(idx_hh + 224);

    auto ts_xyyyy_yyyyy = pbuffer.data(idx_hh + 225);

    auto ts_xyyyy_yyyyz = pbuffer.data(idx_hh + 226);

    auto ts_xyyyy_yyyzz = pbuffer.data(idx_hh + 227);

    auto ts_xyyyy_yyzzz = pbuffer.data(idx_hh + 228);

    auto ts_xyyyy_yzzzz = pbuffer.data(idx_hh + 229);

    auto ts_xyyyy_zzzzz = pbuffer.data(idx_hh + 230);

    auto ts_xyyyz_xxxxx = pbuffer.data(idx_hh + 231);

    auto ts_xyyyz_xxxxy = pbuffer.data(idx_hh + 232);

    auto ts_xyyyz_xxxxz = pbuffer.data(idx_hh + 233);

    auto ts_xyyyz_xxxyy = pbuffer.data(idx_hh + 234);

    auto ts_xyyyz_xxxyz = pbuffer.data(idx_hh + 235);

    auto ts_xyyyz_xxxzz = pbuffer.data(idx_hh + 236);

    auto ts_xyyyz_xxyyy = pbuffer.data(idx_hh + 237);

    auto ts_xyyyz_xxyyz = pbuffer.data(idx_hh + 238);

    auto ts_xyyyz_xxyzz = pbuffer.data(idx_hh + 239);

    auto ts_xyyyz_xxzzz = pbuffer.data(idx_hh + 240);

    auto ts_xyyyz_xyyyy = pbuffer.data(idx_hh + 241);

    auto ts_xyyyz_xyyyz = pbuffer.data(idx_hh + 242);

    auto ts_xyyyz_xyyzz = pbuffer.data(idx_hh + 243);

    auto ts_xyyyz_xyzzz = pbuffer.data(idx_hh + 244);

    auto ts_xyyyz_xzzzz = pbuffer.data(idx_hh + 245);

    auto ts_xyyyz_yyyyy = pbuffer.data(idx_hh + 246);

    auto ts_xyyyz_yyyyz = pbuffer.data(idx_hh + 247);

    auto ts_xyyyz_yyyzz = pbuffer.data(idx_hh + 248);

    auto ts_xyyyz_yyzzz = pbuffer.data(idx_hh + 249);

    auto ts_xyyyz_yzzzz = pbuffer.data(idx_hh + 250);

    auto ts_xyyyz_zzzzz = pbuffer.data(idx_hh + 251);

    auto ts_xyyzz_xxxxx = pbuffer.data(idx_hh + 252);

    auto ts_xyyzz_xxxxy = pbuffer.data(idx_hh + 253);

    auto ts_xyyzz_xxxxz = pbuffer.data(idx_hh + 254);

    auto ts_xyyzz_xxxyy = pbuffer.data(idx_hh + 255);

    auto ts_xyyzz_xxxyz = pbuffer.data(idx_hh + 256);

    auto ts_xyyzz_xxxzz = pbuffer.data(idx_hh + 257);

    auto ts_xyyzz_xxyyy = pbuffer.data(idx_hh + 258);

    auto ts_xyyzz_xxyyz = pbuffer.data(idx_hh + 259);

    auto ts_xyyzz_xxyzz = pbuffer.data(idx_hh + 260);

    auto ts_xyyzz_xxzzz = pbuffer.data(idx_hh + 261);

    auto ts_xyyzz_xyyyy = pbuffer.data(idx_hh + 262);

    auto ts_xyyzz_xyyyz = pbuffer.data(idx_hh + 263);

    auto ts_xyyzz_xyyzz = pbuffer.data(idx_hh + 264);

    auto ts_xyyzz_xyzzz = pbuffer.data(idx_hh + 265);

    auto ts_xyyzz_xzzzz = pbuffer.data(idx_hh + 266);

    auto ts_xyyzz_yyyyy = pbuffer.data(idx_hh + 267);

    auto ts_xyyzz_yyyyz = pbuffer.data(idx_hh + 268);

    auto ts_xyyzz_yyyzz = pbuffer.data(idx_hh + 269);

    auto ts_xyyzz_yyzzz = pbuffer.data(idx_hh + 270);

    auto ts_xyyzz_yzzzz = pbuffer.data(idx_hh + 271);

    auto ts_xyyzz_zzzzz = pbuffer.data(idx_hh + 272);

    auto ts_xyzzz_xxxxx = pbuffer.data(idx_hh + 273);

    auto ts_xyzzz_xxxxy = pbuffer.data(idx_hh + 274);

    auto ts_xyzzz_xxxxz = pbuffer.data(idx_hh + 275);

    auto ts_xyzzz_xxxyy = pbuffer.data(idx_hh + 276);

    auto ts_xyzzz_xxxyz = pbuffer.data(idx_hh + 277);

    auto ts_xyzzz_xxxzz = pbuffer.data(idx_hh + 278);

    auto ts_xyzzz_xxyyy = pbuffer.data(idx_hh + 279);

    auto ts_xyzzz_xxyyz = pbuffer.data(idx_hh + 280);

    auto ts_xyzzz_xxyzz = pbuffer.data(idx_hh + 281);

    auto ts_xyzzz_xxzzz = pbuffer.data(idx_hh + 282);

    auto ts_xyzzz_xyyyy = pbuffer.data(idx_hh + 283);

    auto ts_xyzzz_xyyyz = pbuffer.data(idx_hh + 284);

    auto ts_xyzzz_xyyzz = pbuffer.data(idx_hh + 285);

    auto ts_xyzzz_xyzzz = pbuffer.data(idx_hh + 286);

    auto ts_xyzzz_xzzzz = pbuffer.data(idx_hh + 287);

    auto ts_xyzzz_yyyyy = pbuffer.data(idx_hh + 288);

    auto ts_xyzzz_yyyyz = pbuffer.data(idx_hh + 289);

    auto ts_xyzzz_yyyzz = pbuffer.data(idx_hh + 290);

    auto ts_xyzzz_yyzzz = pbuffer.data(idx_hh + 291);

    auto ts_xyzzz_yzzzz = pbuffer.data(idx_hh + 292);

    auto ts_xyzzz_zzzzz = pbuffer.data(idx_hh + 293);

    auto ts_xzzzz_xxxxx = pbuffer.data(idx_hh + 294);

    auto ts_xzzzz_xxxxy = pbuffer.data(idx_hh + 295);

    auto ts_xzzzz_xxxxz = pbuffer.data(idx_hh + 296);

    auto ts_xzzzz_xxxyy = pbuffer.data(idx_hh + 297);

    auto ts_xzzzz_xxxyz = pbuffer.data(idx_hh + 298);

    auto ts_xzzzz_xxxzz = pbuffer.data(idx_hh + 299);

    auto ts_xzzzz_xxyyy = pbuffer.data(idx_hh + 300);

    auto ts_xzzzz_xxyyz = pbuffer.data(idx_hh + 301);

    auto ts_xzzzz_xxyzz = pbuffer.data(idx_hh + 302);

    auto ts_xzzzz_xxzzz = pbuffer.data(idx_hh + 303);

    auto ts_xzzzz_xyyyy = pbuffer.data(idx_hh + 304);

    auto ts_xzzzz_xyyyz = pbuffer.data(idx_hh + 305);

    auto ts_xzzzz_xyyzz = pbuffer.data(idx_hh + 306);

    auto ts_xzzzz_xyzzz = pbuffer.data(idx_hh + 307);

    auto ts_xzzzz_xzzzz = pbuffer.data(idx_hh + 308);

    auto ts_xzzzz_yyyyy = pbuffer.data(idx_hh + 309);

    auto ts_xzzzz_yyyyz = pbuffer.data(idx_hh + 310);

    auto ts_xzzzz_yyyzz = pbuffer.data(idx_hh + 311);

    auto ts_xzzzz_yyzzz = pbuffer.data(idx_hh + 312);

    auto ts_xzzzz_yzzzz = pbuffer.data(idx_hh + 313);

    auto ts_xzzzz_zzzzz = pbuffer.data(idx_hh + 314);

    auto ts_yyyyy_xxxxx = pbuffer.data(idx_hh + 315);

    auto ts_yyyyy_xxxxy = pbuffer.data(idx_hh + 316);

    auto ts_yyyyy_xxxxz = pbuffer.data(idx_hh + 317);

    auto ts_yyyyy_xxxyy = pbuffer.data(idx_hh + 318);

    auto ts_yyyyy_xxxyz = pbuffer.data(idx_hh + 319);

    auto ts_yyyyy_xxxzz = pbuffer.data(idx_hh + 320);

    auto ts_yyyyy_xxyyy = pbuffer.data(idx_hh + 321);

    auto ts_yyyyy_xxyyz = pbuffer.data(idx_hh + 322);

    auto ts_yyyyy_xxyzz = pbuffer.data(idx_hh + 323);

    auto ts_yyyyy_xxzzz = pbuffer.data(idx_hh + 324);

    auto ts_yyyyy_xyyyy = pbuffer.data(idx_hh + 325);

    auto ts_yyyyy_xyyyz = pbuffer.data(idx_hh + 326);

    auto ts_yyyyy_xyyzz = pbuffer.data(idx_hh + 327);

    auto ts_yyyyy_xyzzz = pbuffer.data(idx_hh + 328);

    auto ts_yyyyy_xzzzz = pbuffer.data(idx_hh + 329);

    auto ts_yyyyy_yyyyy = pbuffer.data(idx_hh + 330);

    auto ts_yyyyy_yyyyz = pbuffer.data(idx_hh + 331);

    auto ts_yyyyy_yyyzz = pbuffer.data(idx_hh + 332);

    auto ts_yyyyy_yyzzz = pbuffer.data(idx_hh + 333);

    auto ts_yyyyy_yzzzz = pbuffer.data(idx_hh + 334);

    auto ts_yyyyy_zzzzz = pbuffer.data(idx_hh + 335);

    auto ts_yyyyz_xxxxx = pbuffer.data(idx_hh + 336);

    auto ts_yyyyz_xxxxy = pbuffer.data(idx_hh + 337);

    auto ts_yyyyz_xxxxz = pbuffer.data(idx_hh + 338);

    auto ts_yyyyz_xxxyy = pbuffer.data(idx_hh + 339);

    auto ts_yyyyz_xxxyz = pbuffer.data(idx_hh + 340);

    auto ts_yyyyz_xxxzz = pbuffer.data(idx_hh + 341);

    auto ts_yyyyz_xxyyy = pbuffer.data(idx_hh + 342);

    auto ts_yyyyz_xxyyz = pbuffer.data(idx_hh + 343);

    auto ts_yyyyz_xxyzz = pbuffer.data(idx_hh + 344);

    auto ts_yyyyz_xxzzz = pbuffer.data(idx_hh + 345);

    auto ts_yyyyz_xyyyy = pbuffer.data(idx_hh + 346);

    auto ts_yyyyz_xyyyz = pbuffer.data(idx_hh + 347);

    auto ts_yyyyz_xyyzz = pbuffer.data(idx_hh + 348);

    auto ts_yyyyz_xyzzz = pbuffer.data(idx_hh + 349);

    auto ts_yyyyz_xzzzz = pbuffer.data(idx_hh + 350);

    auto ts_yyyyz_yyyyy = pbuffer.data(idx_hh + 351);

    auto ts_yyyyz_yyyyz = pbuffer.data(idx_hh + 352);

    auto ts_yyyyz_yyyzz = pbuffer.data(idx_hh + 353);

    auto ts_yyyyz_yyzzz = pbuffer.data(idx_hh + 354);

    auto ts_yyyyz_yzzzz = pbuffer.data(idx_hh + 355);

    auto ts_yyyyz_zzzzz = pbuffer.data(idx_hh + 356);

    auto ts_yyyzz_xxxxx = pbuffer.data(idx_hh + 357);

    auto ts_yyyzz_xxxxy = pbuffer.data(idx_hh + 358);

    auto ts_yyyzz_xxxxz = pbuffer.data(idx_hh + 359);

    auto ts_yyyzz_xxxyy = pbuffer.data(idx_hh + 360);

    auto ts_yyyzz_xxxyz = pbuffer.data(idx_hh + 361);

    auto ts_yyyzz_xxxzz = pbuffer.data(idx_hh + 362);

    auto ts_yyyzz_xxyyy = pbuffer.data(idx_hh + 363);

    auto ts_yyyzz_xxyyz = pbuffer.data(idx_hh + 364);

    auto ts_yyyzz_xxyzz = pbuffer.data(idx_hh + 365);

    auto ts_yyyzz_xxzzz = pbuffer.data(idx_hh + 366);

    auto ts_yyyzz_xyyyy = pbuffer.data(idx_hh + 367);

    auto ts_yyyzz_xyyyz = pbuffer.data(idx_hh + 368);

    auto ts_yyyzz_xyyzz = pbuffer.data(idx_hh + 369);

    auto ts_yyyzz_xyzzz = pbuffer.data(idx_hh + 370);

    auto ts_yyyzz_xzzzz = pbuffer.data(idx_hh + 371);

    auto ts_yyyzz_yyyyy = pbuffer.data(idx_hh + 372);

    auto ts_yyyzz_yyyyz = pbuffer.data(idx_hh + 373);

    auto ts_yyyzz_yyyzz = pbuffer.data(idx_hh + 374);

    auto ts_yyyzz_yyzzz = pbuffer.data(idx_hh + 375);

    auto ts_yyyzz_yzzzz = pbuffer.data(idx_hh + 376);

    auto ts_yyyzz_zzzzz = pbuffer.data(idx_hh + 377);

    auto ts_yyzzz_xxxxx = pbuffer.data(idx_hh + 378);

    auto ts_yyzzz_xxxxy = pbuffer.data(idx_hh + 379);

    auto ts_yyzzz_xxxxz = pbuffer.data(idx_hh + 380);

    auto ts_yyzzz_xxxyy = pbuffer.data(idx_hh + 381);

    auto ts_yyzzz_xxxyz = pbuffer.data(idx_hh + 382);

    auto ts_yyzzz_xxxzz = pbuffer.data(idx_hh + 383);

    auto ts_yyzzz_xxyyy = pbuffer.data(idx_hh + 384);

    auto ts_yyzzz_xxyyz = pbuffer.data(idx_hh + 385);

    auto ts_yyzzz_xxyzz = pbuffer.data(idx_hh + 386);

    auto ts_yyzzz_xxzzz = pbuffer.data(idx_hh + 387);

    auto ts_yyzzz_xyyyy = pbuffer.data(idx_hh + 388);

    auto ts_yyzzz_xyyyz = pbuffer.data(idx_hh + 389);

    auto ts_yyzzz_xyyzz = pbuffer.data(idx_hh + 390);

    auto ts_yyzzz_xyzzz = pbuffer.data(idx_hh + 391);

    auto ts_yyzzz_xzzzz = pbuffer.data(idx_hh + 392);

    auto ts_yyzzz_yyyyy = pbuffer.data(idx_hh + 393);

    auto ts_yyzzz_yyyyz = pbuffer.data(idx_hh + 394);

    auto ts_yyzzz_yyyzz = pbuffer.data(idx_hh + 395);

    auto ts_yyzzz_yyzzz = pbuffer.data(idx_hh + 396);

    auto ts_yyzzz_yzzzz = pbuffer.data(idx_hh + 397);

    auto ts_yyzzz_zzzzz = pbuffer.data(idx_hh + 398);

    auto ts_yzzzz_xxxxx = pbuffer.data(idx_hh + 399);

    auto ts_yzzzz_xxxxy = pbuffer.data(idx_hh + 400);

    auto ts_yzzzz_xxxxz = pbuffer.data(idx_hh + 401);

    auto ts_yzzzz_xxxyy = pbuffer.data(idx_hh + 402);

    auto ts_yzzzz_xxxyz = pbuffer.data(idx_hh + 403);

    auto ts_yzzzz_xxxzz = pbuffer.data(idx_hh + 404);

    auto ts_yzzzz_xxyyy = pbuffer.data(idx_hh + 405);

    auto ts_yzzzz_xxyyz = pbuffer.data(idx_hh + 406);

    auto ts_yzzzz_xxyzz = pbuffer.data(idx_hh + 407);

    auto ts_yzzzz_xxzzz = pbuffer.data(idx_hh + 408);

    auto ts_yzzzz_xyyyy = pbuffer.data(idx_hh + 409);

    auto ts_yzzzz_xyyyz = pbuffer.data(idx_hh + 410);

    auto ts_yzzzz_xyyzz = pbuffer.data(idx_hh + 411);

    auto ts_yzzzz_xyzzz = pbuffer.data(idx_hh + 412);

    auto ts_yzzzz_xzzzz = pbuffer.data(idx_hh + 413);

    auto ts_yzzzz_yyyyy = pbuffer.data(idx_hh + 414);

    auto ts_yzzzz_yyyyz = pbuffer.data(idx_hh + 415);

    auto ts_yzzzz_yyyzz = pbuffer.data(idx_hh + 416);

    auto ts_yzzzz_yyzzz = pbuffer.data(idx_hh + 417);

    auto ts_yzzzz_yzzzz = pbuffer.data(idx_hh + 418);

    auto ts_yzzzz_zzzzz = pbuffer.data(idx_hh + 419);

    auto ts_zzzzz_xxxxx = pbuffer.data(idx_hh + 420);

    auto ts_zzzzz_xxxxy = pbuffer.data(idx_hh + 421);

    auto ts_zzzzz_xxxxz = pbuffer.data(idx_hh + 422);

    auto ts_zzzzz_xxxyy = pbuffer.data(idx_hh + 423);

    auto ts_zzzzz_xxxyz = pbuffer.data(idx_hh + 424);

    auto ts_zzzzz_xxxzz = pbuffer.data(idx_hh + 425);

    auto ts_zzzzz_xxyyy = pbuffer.data(idx_hh + 426);

    auto ts_zzzzz_xxyyz = pbuffer.data(idx_hh + 427);

    auto ts_zzzzz_xxyzz = pbuffer.data(idx_hh + 428);

    auto ts_zzzzz_xxzzz = pbuffer.data(idx_hh + 429);

    auto ts_zzzzz_xyyyy = pbuffer.data(idx_hh + 430);

    auto ts_zzzzz_xyyyz = pbuffer.data(idx_hh + 431);

    auto ts_zzzzz_xyyzz = pbuffer.data(idx_hh + 432);

    auto ts_zzzzz_xyzzz = pbuffer.data(idx_hh + 433);

    auto ts_zzzzz_xzzzz = pbuffer.data(idx_hh + 434);

    auto ts_zzzzz_yyyyy = pbuffer.data(idx_hh + 435);

    auto ts_zzzzz_yyyyz = pbuffer.data(idx_hh + 436);

    auto ts_zzzzz_yyyzz = pbuffer.data(idx_hh + 437);

    auto ts_zzzzz_yyzzz = pbuffer.data(idx_hh + 438);

    auto ts_zzzzz_yzzzz = pbuffer.data(idx_hh + 439);

    auto ts_zzzzz_zzzzz = pbuffer.data(idx_hh + 440);

    // Set up 0-21 components of targeted buffer : HH

    auto gs_x_xxxxx_xxxxx = pbuffer.data(idx_g_hh);

    auto gs_x_xxxxx_xxxxy = pbuffer.data(idx_g_hh + 1);

    auto gs_x_xxxxx_xxxxz = pbuffer.data(idx_g_hh + 2);

    auto gs_x_xxxxx_xxxyy = pbuffer.data(idx_g_hh + 3);

    auto gs_x_xxxxx_xxxyz = pbuffer.data(idx_g_hh + 4);

    auto gs_x_xxxxx_xxxzz = pbuffer.data(idx_g_hh + 5);

    auto gs_x_xxxxx_xxyyy = pbuffer.data(idx_g_hh + 6);

    auto gs_x_xxxxx_xxyyz = pbuffer.data(idx_g_hh + 7);

    auto gs_x_xxxxx_xxyzz = pbuffer.data(idx_g_hh + 8);

    auto gs_x_xxxxx_xxzzz = pbuffer.data(idx_g_hh + 9);

    auto gs_x_xxxxx_xyyyy = pbuffer.data(idx_g_hh + 10);

    auto gs_x_xxxxx_xyyyz = pbuffer.data(idx_g_hh + 11);

    auto gs_x_xxxxx_xyyzz = pbuffer.data(idx_g_hh + 12);

    auto gs_x_xxxxx_xyzzz = pbuffer.data(idx_g_hh + 13);

    auto gs_x_xxxxx_xzzzz = pbuffer.data(idx_g_hh + 14);

    auto gs_x_xxxxx_yyyyy = pbuffer.data(idx_g_hh + 15);

    auto gs_x_xxxxx_yyyyz = pbuffer.data(idx_g_hh + 16);

    auto gs_x_xxxxx_yyyzz = pbuffer.data(idx_g_hh + 17);

    auto gs_x_xxxxx_yyzzz = pbuffer.data(idx_g_hh + 18);

    auto gs_x_xxxxx_yzzzz = pbuffer.data(idx_g_hh + 19);

    auto gs_x_xxxxx_zzzzz = pbuffer.data(idx_g_hh + 20);

    #pragma omp simd aligned(gc_x, gs_x_xxxxx_xxxxx, gs_x_xxxxx_xxxxy, gs_x_xxxxx_xxxxz, gs_x_xxxxx_xxxyy, gs_x_xxxxx_xxxyz, gs_x_xxxxx_xxxzz, gs_x_xxxxx_xxyyy, gs_x_xxxxx_xxyyz, gs_x_xxxxx_xxyzz, gs_x_xxxxx_xxzzz, gs_x_xxxxx_xyyyy, gs_x_xxxxx_xyyyz, gs_x_xxxxx_xyyzz, gs_x_xxxxx_xyzzz, gs_x_xxxxx_xzzzz, gs_x_xxxxx_yyyyy, gs_x_xxxxx_yyyyz, gs_x_xxxxx_yyyzz, gs_x_xxxxx_yyzzz, gs_x_xxxxx_yzzzz, gs_x_xxxxx_zzzzz, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxzz, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyzz, ts_xxxx_xxzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyzz, ts_xxxx_xyzzz, ts_xxxx_xzzzz, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyzz, ts_xxxx_yyzzz, ts_xxxx_yzzzz, ts_xxxx_zzzzz, ts_xxxxx_xxxx, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxxz, ts_xxxxx_xxxy, ts_xxxxx_xxxyy, ts_xxxxx_xxxyz, ts_xxxxx_xxxz, ts_xxxxx_xxxzz, ts_xxxxx_xxyy, ts_xxxxx_xxyyy, ts_xxxxx_xxyyz, ts_xxxxx_xxyz, ts_xxxxx_xxyzz, ts_xxxxx_xxzz, ts_xxxxx_xxzzz, ts_xxxxx_xyyy, ts_xxxxx_xyyyy, ts_xxxxx_xyyyz, ts_xxxxx_xyyz, ts_xxxxx_xyyzz, ts_xxxxx_xyzz, ts_xxxxx_xyzzz, ts_xxxxx_xzzz, ts_xxxxx_xzzzz, ts_xxxxx_yyyy, ts_xxxxx_yyyyy, ts_xxxxx_yyyyz, ts_xxxxx_yyyz, ts_xxxxx_yyyzz, ts_xxxxx_yyzz, ts_xxxxx_yyzzz, ts_xxxxx_yzzz, ts_xxxxx_yzzzz, ts_xxxxx_zzzz, ts_xxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxx_xxxxx[i] = 10.0 * ts_xxxx_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxxxy[i] = 10.0 * ts_xxxx_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxxxz[i] = 10.0 * ts_xxxx_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxxyy[i] = 10.0 * ts_xxxx_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxxyz[i] = 10.0 * ts_xxxx_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxxzz[i] = 10.0 * ts_xxxx_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxyyy[i] = 10.0 * ts_xxxx_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxyyz[i] = 10.0 * ts_xxxx_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxyzz[i] = 10.0 * ts_xxxx_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xxzzz[i] = 10.0 * ts_xxxx_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyyyy[i] = 10.0 * ts_xxxx_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyyyz[i] = 10.0 * ts_xxxx_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyyzz[i] = 10.0 * ts_xxxx_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xyzzz[i] = 10.0 * ts_xxxx_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xzzzz[i] = 10.0 * ts_xxxx_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyyyy[i] = 10.0 * ts_xxxx_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyyyz[i] = 10.0 * ts_xxxx_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyyzz[i] = 10.0 * ts_xxxx_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yyzzz[i] = 10.0 * ts_xxxx_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yzzzz[i] = 10.0 * ts_xxxx_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_zzzzz[i] = 10.0 * ts_xxxx_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 21-42 components of targeted buffer : HH

    auto gs_x_xxxxy_xxxxx = pbuffer.data(idx_g_hh + 21);

    auto gs_x_xxxxy_xxxxy = pbuffer.data(idx_g_hh + 22);

    auto gs_x_xxxxy_xxxxz = pbuffer.data(idx_g_hh + 23);

    auto gs_x_xxxxy_xxxyy = pbuffer.data(idx_g_hh + 24);

    auto gs_x_xxxxy_xxxyz = pbuffer.data(idx_g_hh + 25);

    auto gs_x_xxxxy_xxxzz = pbuffer.data(idx_g_hh + 26);

    auto gs_x_xxxxy_xxyyy = pbuffer.data(idx_g_hh + 27);

    auto gs_x_xxxxy_xxyyz = pbuffer.data(idx_g_hh + 28);

    auto gs_x_xxxxy_xxyzz = pbuffer.data(idx_g_hh + 29);

    auto gs_x_xxxxy_xxzzz = pbuffer.data(idx_g_hh + 30);

    auto gs_x_xxxxy_xyyyy = pbuffer.data(idx_g_hh + 31);

    auto gs_x_xxxxy_xyyyz = pbuffer.data(idx_g_hh + 32);

    auto gs_x_xxxxy_xyyzz = pbuffer.data(idx_g_hh + 33);

    auto gs_x_xxxxy_xyzzz = pbuffer.data(idx_g_hh + 34);

    auto gs_x_xxxxy_xzzzz = pbuffer.data(idx_g_hh + 35);

    auto gs_x_xxxxy_yyyyy = pbuffer.data(idx_g_hh + 36);

    auto gs_x_xxxxy_yyyyz = pbuffer.data(idx_g_hh + 37);

    auto gs_x_xxxxy_yyyzz = pbuffer.data(idx_g_hh + 38);

    auto gs_x_xxxxy_yyzzz = pbuffer.data(idx_g_hh + 39);

    auto gs_x_xxxxy_yzzzz = pbuffer.data(idx_g_hh + 40);

    auto gs_x_xxxxy_zzzzz = pbuffer.data(idx_g_hh + 41);

    #pragma omp simd aligned(gc_x, gs_x_xxxxy_xxxxx, gs_x_xxxxy_xxxxy, gs_x_xxxxy_xxxxz, gs_x_xxxxy_xxxyy, gs_x_xxxxy_xxxyz, gs_x_xxxxy_xxxzz, gs_x_xxxxy_xxyyy, gs_x_xxxxy_xxyyz, gs_x_xxxxy_xxyzz, gs_x_xxxxy_xxzzz, gs_x_xxxxy_xyyyy, gs_x_xxxxy_xyyyz, gs_x_xxxxy_xyyzz, gs_x_xxxxy_xyzzz, gs_x_xxxxy_xzzzz, gs_x_xxxxy_yyyyy, gs_x_xxxxy_yyyyz, gs_x_xxxxy_yyyzz, gs_x_xxxxy_yyzzz, gs_x_xxxxy_yzzzz, gs_x_xxxxy_zzzzz, ts_xxxxy_xxxx, ts_xxxxy_xxxxx, ts_xxxxy_xxxxy, ts_xxxxy_xxxxz, ts_xxxxy_xxxy, ts_xxxxy_xxxyy, ts_xxxxy_xxxyz, ts_xxxxy_xxxz, ts_xxxxy_xxxzz, ts_xxxxy_xxyy, ts_xxxxy_xxyyy, ts_xxxxy_xxyyz, ts_xxxxy_xxyz, ts_xxxxy_xxyzz, ts_xxxxy_xxzz, ts_xxxxy_xxzzz, ts_xxxxy_xyyy, ts_xxxxy_xyyyy, ts_xxxxy_xyyyz, ts_xxxxy_xyyz, ts_xxxxy_xyyzz, ts_xxxxy_xyzz, ts_xxxxy_xyzzz, ts_xxxxy_xzzz, ts_xxxxy_xzzzz, ts_xxxxy_yyyy, ts_xxxxy_yyyyy, ts_xxxxy_yyyyz, ts_xxxxy_yyyz, ts_xxxxy_yyyzz, ts_xxxxy_yyzz, ts_xxxxy_yyzzz, ts_xxxxy_yzzz, ts_xxxxy_yzzzz, ts_xxxxy_zzzz, ts_xxxxy_zzzzz, ts_xxxy_xxxxx, ts_xxxy_xxxxy, ts_xxxy_xxxxz, ts_xxxy_xxxyy, ts_xxxy_xxxyz, ts_xxxy_xxxzz, ts_xxxy_xxyyy, ts_xxxy_xxyyz, ts_xxxy_xxyzz, ts_xxxy_xxzzz, ts_xxxy_xyyyy, ts_xxxy_xyyyz, ts_xxxy_xyyzz, ts_xxxy_xyzzz, ts_xxxy_xzzzz, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyzz, ts_xxxy_yyzzz, ts_xxxy_yzzzz, ts_xxxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxy_xxxxx[i] = 8.0 * ts_xxxy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxxxy[i] = 8.0 * ts_xxxy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxxxz[i] = 8.0 * ts_xxxy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxxyy[i] = 8.0 * ts_xxxy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxxyz[i] = 8.0 * ts_xxxy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxxzz[i] = 8.0 * ts_xxxy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxyyy[i] = 8.0 * ts_xxxy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxyyz[i] = 8.0 * ts_xxxy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxyzz[i] = 8.0 * ts_xxxy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xxzzz[i] = 8.0 * ts_xxxy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyyyy[i] = 8.0 * ts_xxxy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyyyz[i] = 8.0 * ts_xxxy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyyzz[i] = 8.0 * ts_xxxy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xyzzz[i] = 8.0 * ts_xxxy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xzzzz[i] = 8.0 * ts_xxxy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyyyy[i] = 8.0 * ts_xxxy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyyyz[i] = 8.0 * ts_xxxy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyyzz[i] = 8.0 * ts_xxxy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yyzzz[i] = 8.0 * ts_xxxy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yzzzz[i] = 8.0 * ts_xxxy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_zzzzz[i] = 8.0 * ts_xxxy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 42-63 components of targeted buffer : HH

    auto gs_x_xxxxz_xxxxx = pbuffer.data(idx_g_hh + 42);

    auto gs_x_xxxxz_xxxxy = pbuffer.data(idx_g_hh + 43);

    auto gs_x_xxxxz_xxxxz = pbuffer.data(idx_g_hh + 44);

    auto gs_x_xxxxz_xxxyy = pbuffer.data(idx_g_hh + 45);

    auto gs_x_xxxxz_xxxyz = pbuffer.data(idx_g_hh + 46);

    auto gs_x_xxxxz_xxxzz = pbuffer.data(idx_g_hh + 47);

    auto gs_x_xxxxz_xxyyy = pbuffer.data(idx_g_hh + 48);

    auto gs_x_xxxxz_xxyyz = pbuffer.data(idx_g_hh + 49);

    auto gs_x_xxxxz_xxyzz = pbuffer.data(idx_g_hh + 50);

    auto gs_x_xxxxz_xxzzz = pbuffer.data(idx_g_hh + 51);

    auto gs_x_xxxxz_xyyyy = pbuffer.data(idx_g_hh + 52);

    auto gs_x_xxxxz_xyyyz = pbuffer.data(idx_g_hh + 53);

    auto gs_x_xxxxz_xyyzz = pbuffer.data(idx_g_hh + 54);

    auto gs_x_xxxxz_xyzzz = pbuffer.data(idx_g_hh + 55);

    auto gs_x_xxxxz_xzzzz = pbuffer.data(idx_g_hh + 56);

    auto gs_x_xxxxz_yyyyy = pbuffer.data(idx_g_hh + 57);

    auto gs_x_xxxxz_yyyyz = pbuffer.data(idx_g_hh + 58);

    auto gs_x_xxxxz_yyyzz = pbuffer.data(idx_g_hh + 59);

    auto gs_x_xxxxz_yyzzz = pbuffer.data(idx_g_hh + 60);

    auto gs_x_xxxxz_yzzzz = pbuffer.data(idx_g_hh + 61);

    auto gs_x_xxxxz_zzzzz = pbuffer.data(idx_g_hh + 62);

    #pragma omp simd aligned(gc_x, gs_x_xxxxz_xxxxx, gs_x_xxxxz_xxxxy, gs_x_xxxxz_xxxxz, gs_x_xxxxz_xxxyy, gs_x_xxxxz_xxxyz, gs_x_xxxxz_xxxzz, gs_x_xxxxz_xxyyy, gs_x_xxxxz_xxyyz, gs_x_xxxxz_xxyzz, gs_x_xxxxz_xxzzz, gs_x_xxxxz_xyyyy, gs_x_xxxxz_xyyyz, gs_x_xxxxz_xyyzz, gs_x_xxxxz_xyzzz, gs_x_xxxxz_xzzzz, gs_x_xxxxz_yyyyy, gs_x_xxxxz_yyyyz, gs_x_xxxxz_yyyzz, gs_x_xxxxz_yyzzz, gs_x_xxxxz_yzzzz, gs_x_xxxxz_zzzzz, ts_xxxxz_xxxx, ts_xxxxz_xxxxx, ts_xxxxz_xxxxy, ts_xxxxz_xxxxz, ts_xxxxz_xxxy, ts_xxxxz_xxxyy, ts_xxxxz_xxxyz, ts_xxxxz_xxxz, ts_xxxxz_xxxzz, ts_xxxxz_xxyy, ts_xxxxz_xxyyy, ts_xxxxz_xxyyz, ts_xxxxz_xxyz, ts_xxxxz_xxyzz, ts_xxxxz_xxzz, ts_xxxxz_xxzzz, ts_xxxxz_xyyy, ts_xxxxz_xyyyy, ts_xxxxz_xyyyz, ts_xxxxz_xyyz, ts_xxxxz_xyyzz, ts_xxxxz_xyzz, ts_xxxxz_xyzzz, ts_xxxxz_xzzz, ts_xxxxz_xzzzz, ts_xxxxz_yyyy, ts_xxxxz_yyyyy, ts_xxxxz_yyyyz, ts_xxxxz_yyyz, ts_xxxxz_yyyzz, ts_xxxxz_yyzz, ts_xxxxz_yyzzz, ts_xxxxz_yzzz, ts_xxxxz_yzzzz, ts_xxxxz_zzzz, ts_xxxxz_zzzzz, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxxz, ts_xxxz_xxxyy, ts_xxxz_xxxyz, ts_xxxz_xxxzz, ts_xxxz_xxyyy, ts_xxxz_xxyyz, ts_xxxz_xxyzz, ts_xxxz_xxzzz, ts_xxxz_xyyyy, ts_xxxz_xyyyz, ts_xxxz_xyyzz, ts_xxxz_xyzzz, ts_xxxz_xzzzz, ts_xxxz_yyyyy, ts_xxxz_yyyyz, ts_xxxz_yyyzz, ts_xxxz_yyzzz, ts_xxxz_yzzzz, ts_xxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxz_xxxxx[i] = 8.0 * ts_xxxz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxxxy[i] = 8.0 * ts_xxxz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxxxz[i] = 8.0 * ts_xxxz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxxyy[i] = 8.0 * ts_xxxz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxxyz[i] = 8.0 * ts_xxxz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxxzz[i] = 8.0 * ts_xxxz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxyyy[i] = 8.0 * ts_xxxz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxyyz[i] = 8.0 * ts_xxxz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxyzz[i] = 8.0 * ts_xxxz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xxzzz[i] = 8.0 * ts_xxxz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyyyy[i] = 8.0 * ts_xxxz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyyyz[i] = 8.0 * ts_xxxz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyyzz[i] = 8.0 * ts_xxxz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xyzzz[i] = 8.0 * ts_xxxz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xzzzz[i] = 8.0 * ts_xxxz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyyyy[i] = 8.0 * ts_xxxz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyyyz[i] = 8.0 * ts_xxxz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyyzz[i] = 8.0 * ts_xxxz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yyzzz[i] = 8.0 * ts_xxxz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yzzzz[i] = 8.0 * ts_xxxz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_zzzzz[i] = 8.0 * ts_xxxz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 63-84 components of targeted buffer : HH

    auto gs_x_xxxyy_xxxxx = pbuffer.data(idx_g_hh + 63);

    auto gs_x_xxxyy_xxxxy = pbuffer.data(idx_g_hh + 64);

    auto gs_x_xxxyy_xxxxz = pbuffer.data(idx_g_hh + 65);

    auto gs_x_xxxyy_xxxyy = pbuffer.data(idx_g_hh + 66);

    auto gs_x_xxxyy_xxxyz = pbuffer.data(idx_g_hh + 67);

    auto gs_x_xxxyy_xxxzz = pbuffer.data(idx_g_hh + 68);

    auto gs_x_xxxyy_xxyyy = pbuffer.data(idx_g_hh + 69);

    auto gs_x_xxxyy_xxyyz = pbuffer.data(idx_g_hh + 70);

    auto gs_x_xxxyy_xxyzz = pbuffer.data(idx_g_hh + 71);

    auto gs_x_xxxyy_xxzzz = pbuffer.data(idx_g_hh + 72);

    auto gs_x_xxxyy_xyyyy = pbuffer.data(idx_g_hh + 73);

    auto gs_x_xxxyy_xyyyz = pbuffer.data(idx_g_hh + 74);

    auto gs_x_xxxyy_xyyzz = pbuffer.data(idx_g_hh + 75);

    auto gs_x_xxxyy_xyzzz = pbuffer.data(idx_g_hh + 76);

    auto gs_x_xxxyy_xzzzz = pbuffer.data(idx_g_hh + 77);

    auto gs_x_xxxyy_yyyyy = pbuffer.data(idx_g_hh + 78);

    auto gs_x_xxxyy_yyyyz = pbuffer.data(idx_g_hh + 79);

    auto gs_x_xxxyy_yyyzz = pbuffer.data(idx_g_hh + 80);

    auto gs_x_xxxyy_yyzzz = pbuffer.data(idx_g_hh + 81);

    auto gs_x_xxxyy_yzzzz = pbuffer.data(idx_g_hh + 82);

    auto gs_x_xxxyy_zzzzz = pbuffer.data(idx_g_hh + 83);

    #pragma omp simd aligned(gc_x, gs_x_xxxyy_xxxxx, gs_x_xxxyy_xxxxy, gs_x_xxxyy_xxxxz, gs_x_xxxyy_xxxyy, gs_x_xxxyy_xxxyz, gs_x_xxxyy_xxxzz, gs_x_xxxyy_xxyyy, gs_x_xxxyy_xxyyz, gs_x_xxxyy_xxyzz, gs_x_xxxyy_xxzzz, gs_x_xxxyy_xyyyy, gs_x_xxxyy_xyyyz, gs_x_xxxyy_xyyzz, gs_x_xxxyy_xyzzz, gs_x_xxxyy_xzzzz, gs_x_xxxyy_yyyyy, gs_x_xxxyy_yyyyz, gs_x_xxxyy_yyyzz, gs_x_xxxyy_yyzzz, gs_x_xxxyy_yzzzz, gs_x_xxxyy_zzzzz, ts_xxxyy_xxxx, ts_xxxyy_xxxxx, ts_xxxyy_xxxxy, ts_xxxyy_xxxxz, ts_xxxyy_xxxy, ts_xxxyy_xxxyy, ts_xxxyy_xxxyz, ts_xxxyy_xxxz, ts_xxxyy_xxxzz, ts_xxxyy_xxyy, ts_xxxyy_xxyyy, ts_xxxyy_xxyyz, ts_xxxyy_xxyz, ts_xxxyy_xxyzz, ts_xxxyy_xxzz, ts_xxxyy_xxzzz, ts_xxxyy_xyyy, ts_xxxyy_xyyyy, ts_xxxyy_xyyyz, ts_xxxyy_xyyz, ts_xxxyy_xyyzz, ts_xxxyy_xyzz, ts_xxxyy_xyzzz, ts_xxxyy_xzzz, ts_xxxyy_xzzzz, ts_xxxyy_yyyy, ts_xxxyy_yyyyy, ts_xxxyy_yyyyz, ts_xxxyy_yyyz, ts_xxxyy_yyyzz, ts_xxxyy_yyzz, ts_xxxyy_yyzzz, ts_xxxyy_yzzz, ts_xxxyy_yzzzz, ts_xxxyy_zzzz, ts_xxxyy_zzzzz, ts_xxyy_xxxxx, ts_xxyy_xxxxy, ts_xxyy_xxxxz, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxxzz, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyzz, ts_xxyy_xxzzz, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyzz, ts_xxyy_xyzzz, ts_xxyy_xzzzz, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyzz, ts_xxyy_yyzzz, ts_xxyy_yzzzz, ts_xxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyy_xxxxx[i] = 6.0 * ts_xxyy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxxxy[i] = 6.0 * ts_xxyy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxxxz[i] = 6.0 * ts_xxyy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxxyy[i] = 6.0 * ts_xxyy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxxyz[i] = 6.0 * ts_xxyy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxxzz[i] = 6.0 * ts_xxyy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxyyy[i] = 6.0 * ts_xxyy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxyyz[i] = 6.0 * ts_xxyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxyzz[i] = 6.0 * ts_xxyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xxzzz[i] = 6.0 * ts_xxyy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyyyy[i] = 6.0 * ts_xxyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyyyz[i] = 6.0 * ts_xxyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyyzz[i] = 6.0 * ts_xxyy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xyzzz[i] = 6.0 * ts_xxyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xzzzz[i] = 6.0 * ts_xxyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyyyy[i] = 6.0 * ts_xxyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyyyz[i] = 6.0 * ts_xxyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyyzz[i] = 6.0 * ts_xxyy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yyzzz[i] = 6.0 * ts_xxyy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yzzzz[i] = 6.0 * ts_xxyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_zzzzz[i] = 6.0 * ts_xxyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 84-105 components of targeted buffer : HH

    auto gs_x_xxxyz_xxxxx = pbuffer.data(idx_g_hh + 84);

    auto gs_x_xxxyz_xxxxy = pbuffer.data(idx_g_hh + 85);

    auto gs_x_xxxyz_xxxxz = pbuffer.data(idx_g_hh + 86);

    auto gs_x_xxxyz_xxxyy = pbuffer.data(idx_g_hh + 87);

    auto gs_x_xxxyz_xxxyz = pbuffer.data(idx_g_hh + 88);

    auto gs_x_xxxyz_xxxzz = pbuffer.data(idx_g_hh + 89);

    auto gs_x_xxxyz_xxyyy = pbuffer.data(idx_g_hh + 90);

    auto gs_x_xxxyz_xxyyz = pbuffer.data(idx_g_hh + 91);

    auto gs_x_xxxyz_xxyzz = pbuffer.data(idx_g_hh + 92);

    auto gs_x_xxxyz_xxzzz = pbuffer.data(idx_g_hh + 93);

    auto gs_x_xxxyz_xyyyy = pbuffer.data(idx_g_hh + 94);

    auto gs_x_xxxyz_xyyyz = pbuffer.data(idx_g_hh + 95);

    auto gs_x_xxxyz_xyyzz = pbuffer.data(idx_g_hh + 96);

    auto gs_x_xxxyz_xyzzz = pbuffer.data(idx_g_hh + 97);

    auto gs_x_xxxyz_xzzzz = pbuffer.data(idx_g_hh + 98);

    auto gs_x_xxxyz_yyyyy = pbuffer.data(idx_g_hh + 99);

    auto gs_x_xxxyz_yyyyz = pbuffer.data(idx_g_hh + 100);

    auto gs_x_xxxyz_yyyzz = pbuffer.data(idx_g_hh + 101);

    auto gs_x_xxxyz_yyzzz = pbuffer.data(idx_g_hh + 102);

    auto gs_x_xxxyz_yzzzz = pbuffer.data(idx_g_hh + 103);

    auto gs_x_xxxyz_zzzzz = pbuffer.data(idx_g_hh + 104);

    #pragma omp simd aligned(gc_x, gs_x_xxxyz_xxxxx, gs_x_xxxyz_xxxxy, gs_x_xxxyz_xxxxz, gs_x_xxxyz_xxxyy, gs_x_xxxyz_xxxyz, gs_x_xxxyz_xxxzz, gs_x_xxxyz_xxyyy, gs_x_xxxyz_xxyyz, gs_x_xxxyz_xxyzz, gs_x_xxxyz_xxzzz, gs_x_xxxyz_xyyyy, gs_x_xxxyz_xyyyz, gs_x_xxxyz_xyyzz, gs_x_xxxyz_xyzzz, gs_x_xxxyz_xzzzz, gs_x_xxxyz_yyyyy, gs_x_xxxyz_yyyyz, gs_x_xxxyz_yyyzz, gs_x_xxxyz_yyzzz, gs_x_xxxyz_yzzzz, gs_x_xxxyz_zzzzz, ts_xxxyz_xxxx, ts_xxxyz_xxxxx, ts_xxxyz_xxxxy, ts_xxxyz_xxxxz, ts_xxxyz_xxxy, ts_xxxyz_xxxyy, ts_xxxyz_xxxyz, ts_xxxyz_xxxz, ts_xxxyz_xxxzz, ts_xxxyz_xxyy, ts_xxxyz_xxyyy, ts_xxxyz_xxyyz, ts_xxxyz_xxyz, ts_xxxyz_xxyzz, ts_xxxyz_xxzz, ts_xxxyz_xxzzz, ts_xxxyz_xyyy, ts_xxxyz_xyyyy, ts_xxxyz_xyyyz, ts_xxxyz_xyyz, ts_xxxyz_xyyzz, ts_xxxyz_xyzz, ts_xxxyz_xyzzz, ts_xxxyz_xzzz, ts_xxxyz_xzzzz, ts_xxxyz_yyyy, ts_xxxyz_yyyyy, ts_xxxyz_yyyyz, ts_xxxyz_yyyz, ts_xxxyz_yyyzz, ts_xxxyz_yyzz, ts_xxxyz_yyzzz, ts_xxxyz_yzzz, ts_xxxyz_yzzzz, ts_xxxyz_zzzz, ts_xxxyz_zzzzz, ts_xxyz_xxxxx, ts_xxyz_xxxxy, ts_xxyz_xxxxz, ts_xxyz_xxxyy, ts_xxyz_xxxyz, ts_xxyz_xxxzz, ts_xxyz_xxyyy, ts_xxyz_xxyyz, ts_xxyz_xxyzz, ts_xxyz_xxzzz, ts_xxyz_xyyyy, ts_xxyz_xyyyz, ts_xxyz_xyyzz, ts_xxyz_xyzzz, ts_xxyz_xzzzz, ts_xxyz_yyyyy, ts_xxyz_yyyyz, ts_xxyz_yyyzz, ts_xxyz_yyzzz, ts_xxyz_yzzzz, ts_xxyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyz_xxxxx[i] = 6.0 * ts_xxyz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxxxy[i] = 6.0 * ts_xxyz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxxxz[i] = 6.0 * ts_xxyz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxxyy[i] = 6.0 * ts_xxyz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxxyz[i] = 6.0 * ts_xxyz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxxzz[i] = 6.0 * ts_xxyz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxyyy[i] = 6.0 * ts_xxyz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxyyz[i] = 6.0 * ts_xxyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxyzz[i] = 6.0 * ts_xxyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xxzzz[i] = 6.0 * ts_xxyz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyyyy[i] = 6.0 * ts_xxyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyyyz[i] = 6.0 * ts_xxyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyyzz[i] = 6.0 * ts_xxyz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xyzzz[i] = 6.0 * ts_xxyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xzzzz[i] = 6.0 * ts_xxyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyyyy[i] = 6.0 * ts_xxyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyyyz[i] = 6.0 * ts_xxyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyyzz[i] = 6.0 * ts_xxyz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yyzzz[i] = 6.0 * ts_xxyz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yzzzz[i] = 6.0 * ts_xxyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_zzzzz[i] = 6.0 * ts_xxyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 105-126 components of targeted buffer : HH

    auto gs_x_xxxzz_xxxxx = pbuffer.data(idx_g_hh + 105);

    auto gs_x_xxxzz_xxxxy = pbuffer.data(idx_g_hh + 106);

    auto gs_x_xxxzz_xxxxz = pbuffer.data(idx_g_hh + 107);

    auto gs_x_xxxzz_xxxyy = pbuffer.data(idx_g_hh + 108);

    auto gs_x_xxxzz_xxxyz = pbuffer.data(idx_g_hh + 109);

    auto gs_x_xxxzz_xxxzz = pbuffer.data(idx_g_hh + 110);

    auto gs_x_xxxzz_xxyyy = pbuffer.data(idx_g_hh + 111);

    auto gs_x_xxxzz_xxyyz = pbuffer.data(idx_g_hh + 112);

    auto gs_x_xxxzz_xxyzz = pbuffer.data(idx_g_hh + 113);

    auto gs_x_xxxzz_xxzzz = pbuffer.data(idx_g_hh + 114);

    auto gs_x_xxxzz_xyyyy = pbuffer.data(idx_g_hh + 115);

    auto gs_x_xxxzz_xyyyz = pbuffer.data(idx_g_hh + 116);

    auto gs_x_xxxzz_xyyzz = pbuffer.data(idx_g_hh + 117);

    auto gs_x_xxxzz_xyzzz = pbuffer.data(idx_g_hh + 118);

    auto gs_x_xxxzz_xzzzz = pbuffer.data(idx_g_hh + 119);

    auto gs_x_xxxzz_yyyyy = pbuffer.data(idx_g_hh + 120);

    auto gs_x_xxxzz_yyyyz = pbuffer.data(idx_g_hh + 121);

    auto gs_x_xxxzz_yyyzz = pbuffer.data(idx_g_hh + 122);

    auto gs_x_xxxzz_yyzzz = pbuffer.data(idx_g_hh + 123);

    auto gs_x_xxxzz_yzzzz = pbuffer.data(idx_g_hh + 124);

    auto gs_x_xxxzz_zzzzz = pbuffer.data(idx_g_hh + 125);

    #pragma omp simd aligned(gc_x, gs_x_xxxzz_xxxxx, gs_x_xxxzz_xxxxy, gs_x_xxxzz_xxxxz, gs_x_xxxzz_xxxyy, gs_x_xxxzz_xxxyz, gs_x_xxxzz_xxxzz, gs_x_xxxzz_xxyyy, gs_x_xxxzz_xxyyz, gs_x_xxxzz_xxyzz, gs_x_xxxzz_xxzzz, gs_x_xxxzz_xyyyy, gs_x_xxxzz_xyyyz, gs_x_xxxzz_xyyzz, gs_x_xxxzz_xyzzz, gs_x_xxxzz_xzzzz, gs_x_xxxzz_yyyyy, gs_x_xxxzz_yyyyz, gs_x_xxxzz_yyyzz, gs_x_xxxzz_yyzzz, gs_x_xxxzz_yzzzz, gs_x_xxxzz_zzzzz, ts_xxxzz_xxxx, ts_xxxzz_xxxxx, ts_xxxzz_xxxxy, ts_xxxzz_xxxxz, ts_xxxzz_xxxy, ts_xxxzz_xxxyy, ts_xxxzz_xxxyz, ts_xxxzz_xxxz, ts_xxxzz_xxxzz, ts_xxxzz_xxyy, ts_xxxzz_xxyyy, ts_xxxzz_xxyyz, ts_xxxzz_xxyz, ts_xxxzz_xxyzz, ts_xxxzz_xxzz, ts_xxxzz_xxzzz, ts_xxxzz_xyyy, ts_xxxzz_xyyyy, ts_xxxzz_xyyyz, ts_xxxzz_xyyz, ts_xxxzz_xyyzz, ts_xxxzz_xyzz, ts_xxxzz_xyzzz, ts_xxxzz_xzzz, ts_xxxzz_xzzzz, ts_xxxzz_yyyy, ts_xxxzz_yyyyy, ts_xxxzz_yyyyz, ts_xxxzz_yyyz, ts_xxxzz_yyyzz, ts_xxxzz_yyzz, ts_xxxzz_yyzzz, ts_xxxzz_yzzz, ts_xxxzz_yzzzz, ts_xxxzz_zzzz, ts_xxxzz_zzzzz, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxxz, ts_xxzz_xxxyy, ts_xxzz_xxxyz, ts_xxzz_xxxzz, ts_xxzz_xxyyy, ts_xxzz_xxyyz, ts_xxzz_xxyzz, ts_xxzz_xxzzz, ts_xxzz_xyyyy, ts_xxzz_xyyyz, ts_xxzz_xyyzz, ts_xxzz_xyzzz, ts_xxzz_xzzzz, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyzz, ts_xxzz_yyzzz, ts_xxzz_yzzzz, ts_xxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxzz_xxxxx[i] = 6.0 * ts_xxzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxxzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxxxy[i] = 6.0 * ts_xxzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxxxz[i] = 6.0 * ts_xxzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxxyy[i] = 6.0 * ts_xxzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxxyz[i] = 6.0 * ts_xxzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxxzz[i] = 6.0 * ts_xxzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxyyy[i] = 6.0 * ts_xxzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxyyz[i] = 6.0 * ts_xxzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxyzz[i] = 6.0 * ts_xxzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xxzzz[i] = 6.0 * ts_xxzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyyyy[i] = 6.0 * ts_xxzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyyyz[i] = 6.0 * ts_xxzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyyzz[i] = 6.0 * ts_xxzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xyzzz[i] = 6.0 * ts_xxzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xzzzz[i] = 6.0 * ts_xxzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyyyy[i] = 6.0 * ts_xxzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyyyz[i] = 6.0 * ts_xxzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyyzz[i] = 6.0 * ts_xxzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yyzzz[i] = 6.0 * ts_xxzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yzzzz[i] = 6.0 * ts_xxzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_zzzzz[i] = 6.0 * ts_xxzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 126-147 components of targeted buffer : HH

    auto gs_x_xxyyy_xxxxx = pbuffer.data(idx_g_hh + 126);

    auto gs_x_xxyyy_xxxxy = pbuffer.data(idx_g_hh + 127);

    auto gs_x_xxyyy_xxxxz = pbuffer.data(idx_g_hh + 128);

    auto gs_x_xxyyy_xxxyy = pbuffer.data(idx_g_hh + 129);

    auto gs_x_xxyyy_xxxyz = pbuffer.data(idx_g_hh + 130);

    auto gs_x_xxyyy_xxxzz = pbuffer.data(idx_g_hh + 131);

    auto gs_x_xxyyy_xxyyy = pbuffer.data(idx_g_hh + 132);

    auto gs_x_xxyyy_xxyyz = pbuffer.data(idx_g_hh + 133);

    auto gs_x_xxyyy_xxyzz = pbuffer.data(idx_g_hh + 134);

    auto gs_x_xxyyy_xxzzz = pbuffer.data(idx_g_hh + 135);

    auto gs_x_xxyyy_xyyyy = pbuffer.data(idx_g_hh + 136);

    auto gs_x_xxyyy_xyyyz = pbuffer.data(idx_g_hh + 137);

    auto gs_x_xxyyy_xyyzz = pbuffer.data(idx_g_hh + 138);

    auto gs_x_xxyyy_xyzzz = pbuffer.data(idx_g_hh + 139);

    auto gs_x_xxyyy_xzzzz = pbuffer.data(idx_g_hh + 140);

    auto gs_x_xxyyy_yyyyy = pbuffer.data(idx_g_hh + 141);

    auto gs_x_xxyyy_yyyyz = pbuffer.data(idx_g_hh + 142);

    auto gs_x_xxyyy_yyyzz = pbuffer.data(idx_g_hh + 143);

    auto gs_x_xxyyy_yyzzz = pbuffer.data(idx_g_hh + 144);

    auto gs_x_xxyyy_yzzzz = pbuffer.data(idx_g_hh + 145);

    auto gs_x_xxyyy_zzzzz = pbuffer.data(idx_g_hh + 146);

    #pragma omp simd aligned(gc_x, gs_x_xxyyy_xxxxx, gs_x_xxyyy_xxxxy, gs_x_xxyyy_xxxxz, gs_x_xxyyy_xxxyy, gs_x_xxyyy_xxxyz, gs_x_xxyyy_xxxzz, gs_x_xxyyy_xxyyy, gs_x_xxyyy_xxyyz, gs_x_xxyyy_xxyzz, gs_x_xxyyy_xxzzz, gs_x_xxyyy_xyyyy, gs_x_xxyyy_xyyyz, gs_x_xxyyy_xyyzz, gs_x_xxyyy_xyzzz, gs_x_xxyyy_xzzzz, gs_x_xxyyy_yyyyy, gs_x_xxyyy_yyyyz, gs_x_xxyyy_yyyzz, gs_x_xxyyy_yyzzz, gs_x_xxyyy_yzzzz, gs_x_xxyyy_zzzzz, ts_xxyyy_xxxx, ts_xxyyy_xxxxx, ts_xxyyy_xxxxy, ts_xxyyy_xxxxz, ts_xxyyy_xxxy, ts_xxyyy_xxxyy, ts_xxyyy_xxxyz, ts_xxyyy_xxxz, ts_xxyyy_xxxzz, ts_xxyyy_xxyy, ts_xxyyy_xxyyy, ts_xxyyy_xxyyz, ts_xxyyy_xxyz, ts_xxyyy_xxyzz, ts_xxyyy_xxzz, ts_xxyyy_xxzzz, ts_xxyyy_xyyy, ts_xxyyy_xyyyy, ts_xxyyy_xyyyz, ts_xxyyy_xyyz, ts_xxyyy_xyyzz, ts_xxyyy_xyzz, ts_xxyyy_xyzzz, ts_xxyyy_xzzz, ts_xxyyy_xzzzz, ts_xxyyy_yyyy, ts_xxyyy_yyyyy, ts_xxyyy_yyyyz, ts_xxyyy_yyyz, ts_xxyyy_yyyzz, ts_xxyyy_yyzz, ts_xxyyy_yyzzz, ts_xxyyy_yzzz, ts_xxyyy_yzzzz, ts_xxyyy_zzzz, ts_xxyyy_zzzzz, ts_xyyy_xxxxx, ts_xyyy_xxxxy, ts_xyyy_xxxxz, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxxzz, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyzz, ts_xyyy_xxzzz, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyzz, ts_xyyy_xyzzz, ts_xyyy_xzzzz, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyzz, ts_xyyy_yyzzz, ts_xyyy_yzzzz, ts_xyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyy_xxxxx[i] = 4.0 * ts_xyyy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxxxy[i] = 4.0 * ts_xyyy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxxxz[i] = 4.0 * ts_xyyy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxxyy[i] = 4.0 * ts_xyyy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxxyz[i] = 4.0 * ts_xyyy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxxzz[i] = 4.0 * ts_xyyy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxyyy[i] = 4.0 * ts_xyyy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxyyz[i] = 4.0 * ts_xyyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxyzz[i] = 4.0 * ts_xyyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xxzzz[i] = 4.0 * ts_xyyy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyyyy[i] = 4.0 * ts_xyyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyyyz[i] = 4.0 * ts_xyyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyyzz[i] = 4.0 * ts_xyyy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xyzzz[i] = 4.0 * ts_xyyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xzzzz[i] = 4.0 * ts_xyyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyyyy[i] = 4.0 * ts_xyyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyyyz[i] = 4.0 * ts_xyyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyyzz[i] = 4.0 * ts_xyyy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yyzzz[i] = 4.0 * ts_xyyy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yzzzz[i] = 4.0 * ts_xyyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_zzzzz[i] = 4.0 * ts_xyyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 147-168 components of targeted buffer : HH

    auto gs_x_xxyyz_xxxxx = pbuffer.data(idx_g_hh + 147);

    auto gs_x_xxyyz_xxxxy = pbuffer.data(idx_g_hh + 148);

    auto gs_x_xxyyz_xxxxz = pbuffer.data(idx_g_hh + 149);

    auto gs_x_xxyyz_xxxyy = pbuffer.data(idx_g_hh + 150);

    auto gs_x_xxyyz_xxxyz = pbuffer.data(idx_g_hh + 151);

    auto gs_x_xxyyz_xxxzz = pbuffer.data(idx_g_hh + 152);

    auto gs_x_xxyyz_xxyyy = pbuffer.data(idx_g_hh + 153);

    auto gs_x_xxyyz_xxyyz = pbuffer.data(idx_g_hh + 154);

    auto gs_x_xxyyz_xxyzz = pbuffer.data(idx_g_hh + 155);

    auto gs_x_xxyyz_xxzzz = pbuffer.data(idx_g_hh + 156);

    auto gs_x_xxyyz_xyyyy = pbuffer.data(idx_g_hh + 157);

    auto gs_x_xxyyz_xyyyz = pbuffer.data(idx_g_hh + 158);

    auto gs_x_xxyyz_xyyzz = pbuffer.data(idx_g_hh + 159);

    auto gs_x_xxyyz_xyzzz = pbuffer.data(idx_g_hh + 160);

    auto gs_x_xxyyz_xzzzz = pbuffer.data(idx_g_hh + 161);

    auto gs_x_xxyyz_yyyyy = pbuffer.data(idx_g_hh + 162);

    auto gs_x_xxyyz_yyyyz = pbuffer.data(idx_g_hh + 163);

    auto gs_x_xxyyz_yyyzz = pbuffer.data(idx_g_hh + 164);

    auto gs_x_xxyyz_yyzzz = pbuffer.data(idx_g_hh + 165);

    auto gs_x_xxyyz_yzzzz = pbuffer.data(idx_g_hh + 166);

    auto gs_x_xxyyz_zzzzz = pbuffer.data(idx_g_hh + 167);

    #pragma omp simd aligned(gc_x, gs_x_xxyyz_xxxxx, gs_x_xxyyz_xxxxy, gs_x_xxyyz_xxxxz, gs_x_xxyyz_xxxyy, gs_x_xxyyz_xxxyz, gs_x_xxyyz_xxxzz, gs_x_xxyyz_xxyyy, gs_x_xxyyz_xxyyz, gs_x_xxyyz_xxyzz, gs_x_xxyyz_xxzzz, gs_x_xxyyz_xyyyy, gs_x_xxyyz_xyyyz, gs_x_xxyyz_xyyzz, gs_x_xxyyz_xyzzz, gs_x_xxyyz_xzzzz, gs_x_xxyyz_yyyyy, gs_x_xxyyz_yyyyz, gs_x_xxyyz_yyyzz, gs_x_xxyyz_yyzzz, gs_x_xxyyz_yzzzz, gs_x_xxyyz_zzzzz, ts_xxyyz_xxxx, ts_xxyyz_xxxxx, ts_xxyyz_xxxxy, ts_xxyyz_xxxxz, ts_xxyyz_xxxy, ts_xxyyz_xxxyy, ts_xxyyz_xxxyz, ts_xxyyz_xxxz, ts_xxyyz_xxxzz, ts_xxyyz_xxyy, ts_xxyyz_xxyyy, ts_xxyyz_xxyyz, ts_xxyyz_xxyz, ts_xxyyz_xxyzz, ts_xxyyz_xxzz, ts_xxyyz_xxzzz, ts_xxyyz_xyyy, ts_xxyyz_xyyyy, ts_xxyyz_xyyyz, ts_xxyyz_xyyz, ts_xxyyz_xyyzz, ts_xxyyz_xyzz, ts_xxyyz_xyzzz, ts_xxyyz_xzzz, ts_xxyyz_xzzzz, ts_xxyyz_yyyy, ts_xxyyz_yyyyy, ts_xxyyz_yyyyz, ts_xxyyz_yyyz, ts_xxyyz_yyyzz, ts_xxyyz_yyzz, ts_xxyyz_yyzzz, ts_xxyyz_yzzz, ts_xxyyz_yzzzz, ts_xxyyz_zzzz, ts_xxyyz_zzzzz, ts_xyyz_xxxxx, ts_xyyz_xxxxy, ts_xyyz_xxxxz, ts_xyyz_xxxyy, ts_xyyz_xxxyz, ts_xyyz_xxxzz, ts_xyyz_xxyyy, ts_xyyz_xxyyz, ts_xyyz_xxyzz, ts_xyyz_xxzzz, ts_xyyz_xyyyy, ts_xyyz_xyyyz, ts_xyyz_xyyzz, ts_xyyz_xyzzz, ts_xyyz_xzzzz, ts_xyyz_yyyyy, ts_xyyz_yyyyz, ts_xyyz_yyyzz, ts_xyyz_yyzzz, ts_xyyz_yzzzz, ts_xyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyz_xxxxx[i] = 4.0 * ts_xyyz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxxxy[i] = 4.0 * ts_xyyz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxxxz[i] = 4.0 * ts_xyyz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxxyy[i] = 4.0 * ts_xyyz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxxyz[i] = 4.0 * ts_xyyz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxxzz[i] = 4.0 * ts_xyyz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxyyy[i] = 4.0 * ts_xyyz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxyyz[i] = 4.0 * ts_xyyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxyzz[i] = 4.0 * ts_xyyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xxzzz[i] = 4.0 * ts_xyyz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyyyy[i] = 4.0 * ts_xyyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyyyz[i] = 4.0 * ts_xyyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyyzz[i] = 4.0 * ts_xyyz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xyzzz[i] = 4.0 * ts_xyyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xzzzz[i] = 4.0 * ts_xyyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyyyy[i] = 4.0 * ts_xyyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyyyz[i] = 4.0 * ts_xyyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyyzz[i] = 4.0 * ts_xyyz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yyzzz[i] = 4.0 * ts_xyyz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yzzzz[i] = 4.0 * ts_xyyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_zzzzz[i] = 4.0 * ts_xyyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 168-189 components of targeted buffer : HH

    auto gs_x_xxyzz_xxxxx = pbuffer.data(idx_g_hh + 168);

    auto gs_x_xxyzz_xxxxy = pbuffer.data(idx_g_hh + 169);

    auto gs_x_xxyzz_xxxxz = pbuffer.data(idx_g_hh + 170);

    auto gs_x_xxyzz_xxxyy = pbuffer.data(idx_g_hh + 171);

    auto gs_x_xxyzz_xxxyz = pbuffer.data(idx_g_hh + 172);

    auto gs_x_xxyzz_xxxzz = pbuffer.data(idx_g_hh + 173);

    auto gs_x_xxyzz_xxyyy = pbuffer.data(idx_g_hh + 174);

    auto gs_x_xxyzz_xxyyz = pbuffer.data(idx_g_hh + 175);

    auto gs_x_xxyzz_xxyzz = pbuffer.data(idx_g_hh + 176);

    auto gs_x_xxyzz_xxzzz = pbuffer.data(idx_g_hh + 177);

    auto gs_x_xxyzz_xyyyy = pbuffer.data(idx_g_hh + 178);

    auto gs_x_xxyzz_xyyyz = pbuffer.data(idx_g_hh + 179);

    auto gs_x_xxyzz_xyyzz = pbuffer.data(idx_g_hh + 180);

    auto gs_x_xxyzz_xyzzz = pbuffer.data(idx_g_hh + 181);

    auto gs_x_xxyzz_xzzzz = pbuffer.data(idx_g_hh + 182);

    auto gs_x_xxyzz_yyyyy = pbuffer.data(idx_g_hh + 183);

    auto gs_x_xxyzz_yyyyz = pbuffer.data(idx_g_hh + 184);

    auto gs_x_xxyzz_yyyzz = pbuffer.data(idx_g_hh + 185);

    auto gs_x_xxyzz_yyzzz = pbuffer.data(idx_g_hh + 186);

    auto gs_x_xxyzz_yzzzz = pbuffer.data(idx_g_hh + 187);

    auto gs_x_xxyzz_zzzzz = pbuffer.data(idx_g_hh + 188);

    #pragma omp simd aligned(gc_x, gs_x_xxyzz_xxxxx, gs_x_xxyzz_xxxxy, gs_x_xxyzz_xxxxz, gs_x_xxyzz_xxxyy, gs_x_xxyzz_xxxyz, gs_x_xxyzz_xxxzz, gs_x_xxyzz_xxyyy, gs_x_xxyzz_xxyyz, gs_x_xxyzz_xxyzz, gs_x_xxyzz_xxzzz, gs_x_xxyzz_xyyyy, gs_x_xxyzz_xyyyz, gs_x_xxyzz_xyyzz, gs_x_xxyzz_xyzzz, gs_x_xxyzz_xzzzz, gs_x_xxyzz_yyyyy, gs_x_xxyzz_yyyyz, gs_x_xxyzz_yyyzz, gs_x_xxyzz_yyzzz, gs_x_xxyzz_yzzzz, gs_x_xxyzz_zzzzz, ts_xxyzz_xxxx, ts_xxyzz_xxxxx, ts_xxyzz_xxxxy, ts_xxyzz_xxxxz, ts_xxyzz_xxxy, ts_xxyzz_xxxyy, ts_xxyzz_xxxyz, ts_xxyzz_xxxz, ts_xxyzz_xxxzz, ts_xxyzz_xxyy, ts_xxyzz_xxyyy, ts_xxyzz_xxyyz, ts_xxyzz_xxyz, ts_xxyzz_xxyzz, ts_xxyzz_xxzz, ts_xxyzz_xxzzz, ts_xxyzz_xyyy, ts_xxyzz_xyyyy, ts_xxyzz_xyyyz, ts_xxyzz_xyyz, ts_xxyzz_xyyzz, ts_xxyzz_xyzz, ts_xxyzz_xyzzz, ts_xxyzz_xzzz, ts_xxyzz_xzzzz, ts_xxyzz_yyyy, ts_xxyzz_yyyyy, ts_xxyzz_yyyyz, ts_xxyzz_yyyz, ts_xxyzz_yyyzz, ts_xxyzz_yyzz, ts_xxyzz_yyzzz, ts_xxyzz_yzzz, ts_xxyzz_yzzzz, ts_xxyzz_zzzz, ts_xxyzz_zzzzz, ts_xyzz_xxxxx, ts_xyzz_xxxxy, ts_xyzz_xxxxz, ts_xyzz_xxxyy, ts_xyzz_xxxyz, ts_xyzz_xxxzz, ts_xyzz_xxyyy, ts_xyzz_xxyyz, ts_xyzz_xxyzz, ts_xyzz_xxzzz, ts_xyzz_xyyyy, ts_xyzz_xyyyz, ts_xyzz_xyyzz, ts_xyzz_xyzzz, ts_xyzz_xzzzz, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyzz, ts_xyzz_yyzzz, ts_xyzz_yzzzz, ts_xyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyzz_xxxxx[i] = 4.0 * ts_xyzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxxxy[i] = 4.0 * ts_xyzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxxxz[i] = 4.0 * ts_xyzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxxyy[i] = 4.0 * ts_xyzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxxyz[i] = 4.0 * ts_xyzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxxzz[i] = 4.0 * ts_xyzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxyyy[i] = 4.0 * ts_xyzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxyyz[i] = 4.0 * ts_xyzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxyzz[i] = 4.0 * ts_xyzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xxzzz[i] = 4.0 * ts_xyzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyyyy[i] = 4.0 * ts_xyzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyyyz[i] = 4.0 * ts_xyzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyyzz[i] = 4.0 * ts_xyzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xyzzz[i] = 4.0 * ts_xyzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xzzzz[i] = 4.0 * ts_xyzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyyyy[i] = 4.0 * ts_xyzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyyyz[i] = 4.0 * ts_xyzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyyzz[i] = 4.0 * ts_xyzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yyzzz[i] = 4.0 * ts_xyzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yzzzz[i] = 4.0 * ts_xyzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_zzzzz[i] = 4.0 * ts_xyzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 189-210 components of targeted buffer : HH

    auto gs_x_xxzzz_xxxxx = pbuffer.data(idx_g_hh + 189);

    auto gs_x_xxzzz_xxxxy = pbuffer.data(idx_g_hh + 190);

    auto gs_x_xxzzz_xxxxz = pbuffer.data(idx_g_hh + 191);

    auto gs_x_xxzzz_xxxyy = pbuffer.data(idx_g_hh + 192);

    auto gs_x_xxzzz_xxxyz = pbuffer.data(idx_g_hh + 193);

    auto gs_x_xxzzz_xxxzz = pbuffer.data(idx_g_hh + 194);

    auto gs_x_xxzzz_xxyyy = pbuffer.data(idx_g_hh + 195);

    auto gs_x_xxzzz_xxyyz = pbuffer.data(idx_g_hh + 196);

    auto gs_x_xxzzz_xxyzz = pbuffer.data(idx_g_hh + 197);

    auto gs_x_xxzzz_xxzzz = pbuffer.data(idx_g_hh + 198);

    auto gs_x_xxzzz_xyyyy = pbuffer.data(idx_g_hh + 199);

    auto gs_x_xxzzz_xyyyz = pbuffer.data(idx_g_hh + 200);

    auto gs_x_xxzzz_xyyzz = pbuffer.data(idx_g_hh + 201);

    auto gs_x_xxzzz_xyzzz = pbuffer.data(idx_g_hh + 202);

    auto gs_x_xxzzz_xzzzz = pbuffer.data(idx_g_hh + 203);

    auto gs_x_xxzzz_yyyyy = pbuffer.data(idx_g_hh + 204);

    auto gs_x_xxzzz_yyyyz = pbuffer.data(idx_g_hh + 205);

    auto gs_x_xxzzz_yyyzz = pbuffer.data(idx_g_hh + 206);

    auto gs_x_xxzzz_yyzzz = pbuffer.data(idx_g_hh + 207);

    auto gs_x_xxzzz_yzzzz = pbuffer.data(idx_g_hh + 208);

    auto gs_x_xxzzz_zzzzz = pbuffer.data(idx_g_hh + 209);

    #pragma omp simd aligned(gc_x, gs_x_xxzzz_xxxxx, gs_x_xxzzz_xxxxy, gs_x_xxzzz_xxxxz, gs_x_xxzzz_xxxyy, gs_x_xxzzz_xxxyz, gs_x_xxzzz_xxxzz, gs_x_xxzzz_xxyyy, gs_x_xxzzz_xxyyz, gs_x_xxzzz_xxyzz, gs_x_xxzzz_xxzzz, gs_x_xxzzz_xyyyy, gs_x_xxzzz_xyyyz, gs_x_xxzzz_xyyzz, gs_x_xxzzz_xyzzz, gs_x_xxzzz_xzzzz, gs_x_xxzzz_yyyyy, gs_x_xxzzz_yyyyz, gs_x_xxzzz_yyyzz, gs_x_xxzzz_yyzzz, gs_x_xxzzz_yzzzz, gs_x_xxzzz_zzzzz, ts_xxzzz_xxxx, ts_xxzzz_xxxxx, ts_xxzzz_xxxxy, ts_xxzzz_xxxxz, ts_xxzzz_xxxy, ts_xxzzz_xxxyy, ts_xxzzz_xxxyz, ts_xxzzz_xxxz, ts_xxzzz_xxxzz, ts_xxzzz_xxyy, ts_xxzzz_xxyyy, ts_xxzzz_xxyyz, ts_xxzzz_xxyz, ts_xxzzz_xxyzz, ts_xxzzz_xxzz, ts_xxzzz_xxzzz, ts_xxzzz_xyyy, ts_xxzzz_xyyyy, ts_xxzzz_xyyyz, ts_xxzzz_xyyz, ts_xxzzz_xyyzz, ts_xxzzz_xyzz, ts_xxzzz_xyzzz, ts_xxzzz_xzzz, ts_xxzzz_xzzzz, ts_xxzzz_yyyy, ts_xxzzz_yyyyy, ts_xxzzz_yyyyz, ts_xxzzz_yyyz, ts_xxzzz_yyyzz, ts_xxzzz_yyzz, ts_xxzzz_yyzzz, ts_xxzzz_yzzz, ts_xxzzz_yzzzz, ts_xxzzz_zzzz, ts_xxzzz_zzzzz, ts_xzzz_xxxxx, ts_xzzz_xxxxy, ts_xzzz_xxxxz, ts_xzzz_xxxyy, ts_xzzz_xxxyz, ts_xzzz_xxxzz, ts_xzzz_xxyyy, ts_xzzz_xxyyz, ts_xzzz_xxyzz, ts_xzzz_xxzzz, ts_xzzz_xyyyy, ts_xzzz_xyyyz, ts_xzzz_xyyzz, ts_xzzz_xyzzz, ts_xzzz_xzzzz, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyzz, ts_xzzz_yyzzz, ts_xzzz_yzzzz, ts_xzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzzz_xxxxx[i] = 4.0 * ts_xzzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxxxy[i] = 4.0 * ts_xzzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxxxz[i] = 4.0 * ts_xzzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxxyy[i] = 4.0 * ts_xzzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxxyz[i] = 4.0 * ts_xzzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxxzz[i] = 4.0 * ts_xzzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxyyy[i] = 4.0 * ts_xzzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxyyz[i] = 4.0 * ts_xzzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxyzz[i] = 4.0 * ts_xzzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xxzzz[i] = 4.0 * ts_xzzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyyyy[i] = 4.0 * ts_xzzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyyyz[i] = 4.0 * ts_xzzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyyzz[i] = 4.0 * ts_xzzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xyzzz[i] = 4.0 * ts_xzzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xzzzz[i] = 4.0 * ts_xzzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyyyy[i] = 4.0 * ts_xzzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyyyz[i] = 4.0 * ts_xzzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyyzz[i] = 4.0 * ts_xzzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yyzzz[i] = 4.0 * ts_xzzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yzzzz[i] = 4.0 * ts_xzzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_zzzzz[i] = 4.0 * ts_xzzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 210-231 components of targeted buffer : HH

    auto gs_x_xyyyy_xxxxx = pbuffer.data(idx_g_hh + 210);

    auto gs_x_xyyyy_xxxxy = pbuffer.data(idx_g_hh + 211);

    auto gs_x_xyyyy_xxxxz = pbuffer.data(idx_g_hh + 212);

    auto gs_x_xyyyy_xxxyy = pbuffer.data(idx_g_hh + 213);

    auto gs_x_xyyyy_xxxyz = pbuffer.data(idx_g_hh + 214);

    auto gs_x_xyyyy_xxxzz = pbuffer.data(idx_g_hh + 215);

    auto gs_x_xyyyy_xxyyy = pbuffer.data(idx_g_hh + 216);

    auto gs_x_xyyyy_xxyyz = pbuffer.data(idx_g_hh + 217);

    auto gs_x_xyyyy_xxyzz = pbuffer.data(idx_g_hh + 218);

    auto gs_x_xyyyy_xxzzz = pbuffer.data(idx_g_hh + 219);

    auto gs_x_xyyyy_xyyyy = pbuffer.data(idx_g_hh + 220);

    auto gs_x_xyyyy_xyyyz = pbuffer.data(idx_g_hh + 221);

    auto gs_x_xyyyy_xyyzz = pbuffer.data(idx_g_hh + 222);

    auto gs_x_xyyyy_xyzzz = pbuffer.data(idx_g_hh + 223);

    auto gs_x_xyyyy_xzzzz = pbuffer.data(idx_g_hh + 224);

    auto gs_x_xyyyy_yyyyy = pbuffer.data(idx_g_hh + 225);

    auto gs_x_xyyyy_yyyyz = pbuffer.data(idx_g_hh + 226);

    auto gs_x_xyyyy_yyyzz = pbuffer.data(idx_g_hh + 227);

    auto gs_x_xyyyy_yyzzz = pbuffer.data(idx_g_hh + 228);

    auto gs_x_xyyyy_yzzzz = pbuffer.data(idx_g_hh + 229);

    auto gs_x_xyyyy_zzzzz = pbuffer.data(idx_g_hh + 230);

    #pragma omp simd aligned(gc_x, gs_x_xyyyy_xxxxx, gs_x_xyyyy_xxxxy, gs_x_xyyyy_xxxxz, gs_x_xyyyy_xxxyy, gs_x_xyyyy_xxxyz, gs_x_xyyyy_xxxzz, gs_x_xyyyy_xxyyy, gs_x_xyyyy_xxyyz, gs_x_xyyyy_xxyzz, gs_x_xyyyy_xxzzz, gs_x_xyyyy_xyyyy, gs_x_xyyyy_xyyyz, gs_x_xyyyy_xyyzz, gs_x_xyyyy_xyzzz, gs_x_xyyyy_xzzzz, gs_x_xyyyy_yyyyy, gs_x_xyyyy_yyyyz, gs_x_xyyyy_yyyzz, gs_x_xyyyy_yyzzz, gs_x_xyyyy_yzzzz, gs_x_xyyyy_zzzzz, ts_xyyyy_xxxx, ts_xyyyy_xxxxx, ts_xyyyy_xxxxy, ts_xyyyy_xxxxz, ts_xyyyy_xxxy, ts_xyyyy_xxxyy, ts_xyyyy_xxxyz, ts_xyyyy_xxxz, ts_xyyyy_xxxzz, ts_xyyyy_xxyy, ts_xyyyy_xxyyy, ts_xyyyy_xxyyz, ts_xyyyy_xxyz, ts_xyyyy_xxyzz, ts_xyyyy_xxzz, ts_xyyyy_xxzzz, ts_xyyyy_xyyy, ts_xyyyy_xyyyy, ts_xyyyy_xyyyz, ts_xyyyy_xyyz, ts_xyyyy_xyyzz, ts_xyyyy_xyzz, ts_xyyyy_xyzzz, ts_xyyyy_xzzz, ts_xyyyy_xzzzz, ts_xyyyy_yyyy, ts_xyyyy_yyyyy, ts_xyyyy_yyyyz, ts_xyyyy_yyyz, ts_xyyyy_yyyzz, ts_xyyyy_yyzz, ts_xyyyy_yyzzz, ts_xyyyy_yzzz, ts_xyyyy_yzzzz, ts_xyyyy_zzzz, ts_xyyyy_zzzzz, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxzz, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyzz, ts_yyyy_xxzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyzz, ts_yyyy_xyzzz, ts_yyyy_xzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyzz, ts_yyyy_yyzzz, ts_yyyy_yzzzz, ts_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyy_xxxxx[i] = 2.0 * ts_yyyy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxxxy[i] = 2.0 * ts_yyyy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxxxz[i] = 2.0 * ts_yyyy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxxyy[i] = 2.0 * ts_yyyy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxxyz[i] = 2.0 * ts_yyyy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxxzz[i] = 2.0 * ts_yyyy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxyyy[i] = 2.0 * ts_yyyy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxyyz[i] = 2.0 * ts_yyyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxyzz[i] = 2.0 * ts_yyyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xxzzz[i] = 2.0 * ts_yyyy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyyyy[i] = 2.0 * ts_yyyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyyyz[i] = 2.0 * ts_yyyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyyzz[i] = 2.0 * ts_yyyy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xyzzz[i] = 2.0 * ts_yyyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xzzzz[i] = 2.0 * ts_yyyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyyyy[i] = 2.0 * ts_yyyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyyyz[i] = 2.0 * ts_yyyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyyzz[i] = 2.0 * ts_yyyy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yyzzz[i] = 2.0 * ts_yyyy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yzzzz[i] = 2.0 * ts_yyyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_zzzzz[i] = 2.0 * ts_yyyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 231-252 components of targeted buffer : HH

    auto gs_x_xyyyz_xxxxx = pbuffer.data(idx_g_hh + 231);

    auto gs_x_xyyyz_xxxxy = pbuffer.data(idx_g_hh + 232);

    auto gs_x_xyyyz_xxxxz = pbuffer.data(idx_g_hh + 233);

    auto gs_x_xyyyz_xxxyy = pbuffer.data(idx_g_hh + 234);

    auto gs_x_xyyyz_xxxyz = pbuffer.data(idx_g_hh + 235);

    auto gs_x_xyyyz_xxxzz = pbuffer.data(idx_g_hh + 236);

    auto gs_x_xyyyz_xxyyy = pbuffer.data(idx_g_hh + 237);

    auto gs_x_xyyyz_xxyyz = pbuffer.data(idx_g_hh + 238);

    auto gs_x_xyyyz_xxyzz = pbuffer.data(idx_g_hh + 239);

    auto gs_x_xyyyz_xxzzz = pbuffer.data(idx_g_hh + 240);

    auto gs_x_xyyyz_xyyyy = pbuffer.data(idx_g_hh + 241);

    auto gs_x_xyyyz_xyyyz = pbuffer.data(idx_g_hh + 242);

    auto gs_x_xyyyz_xyyzz = pbuffer.data(idx_g_hh + 243);

    auto gs_x_xyyyz_xyzzz = pbuffer.data(idx_g_hh + 244);

    auto gs_x_xyyyz_xzzzz = pbuffer.data(idx_g_hh + 245);

    auto gs_x_xyyyz_yyyyy = pbuffer.data(idx_g_hh + 246);

    auto gs_x_xyyyz_yyyyz = pbuffer.data(idx_g_hh + 247);

    auto gs_x_xyyyz_yyyzz = pbuffer.data(idx_g_hh + 248);

    auto gs_x_xyyyz_yyzzz = pbuffer.data(idx_g_hh + 249);

    auto gs_x_xyyyz_yzzzz = pbuffer.data(idx_g_hh + 250);

    auto gs_x_xyyyz_zzzzz = pbuffer.data(idx_g_hh + 251);

    #pragma omp simd aligned(gc_x, gs_x_xyyyz_xxxxx, gs_x_xyyyz_xxxxy, gs_x_xyyyz_xxxxz, gs_x_xyyyz_xxxyy, gs_x_xyyyz_xxxyz, gs_x_xyyyz_xxxzz, gs_x_xyyyz_xxyyy, gs_x_xyyyz_xxyyz, gs_x_xyyyz_xxyzz, gs_x_xyyyz_xxzzz, gs_x_xyyyz_xyyyy, gs_x_xyyyz_xyyyz, gs_x_xyyyz_xyyzz, gs_x_xyyyz_xyzzz, gs_x_xyyyz_xzzzz, gs_x_xyyyz_yyyyy, gs_x_xyyyz_yyyyz, gs_x_xyyyz_yyyzz, gs_x_xyyyz_yyzzz, gs_x_xyyyz_yzzzz, gs_x_xyyyz_zzzzz, ts_xyyyz_xxxx, ts_xyyyz_xxxxx, ts_xyyyz_xxxxy, ts_xyyyz_xxxxz, ts_xyyyz_xxxy, ts_xyyyz_xxxyy, ts_xyyyz_xxxyz, ts_xyyyz_xxxz, ts_xyyyz_xxxzz, ts_xyyyz_xxyy, ts_xyyyz_xxyyy, ts_xyyyz_xxyyz, ts_xyyyz_xxyz, ts_xyyyz_xxyzz, ts_xyyyz_xxzz, ts_xyyyz_xxzzz, ts_xyyyz_xyyy, ts_xyyyz_xyyyy, ts_xyyyz_xyyyz, ts_xyyyz_xyyz, ts_xyyyz_xyyzz, ts_xyyyz_xyzz, ts_xyyyz_xyzzz, ts_xyyyz_xzzz, ts_xyyyz_xzzzz, ts_xyyyz_yyyy, ts_xyyyz_yyyyy, ts_xyyyz_yyyyz, ts_xyyyz_yyyz, ts_xyyyz_yyyzz, ts_xyyyz_yyzz, ts_xyyyz_yyzzz, ts_xyyyz_yzzz, ts_xyyyz_yzzzz, ts_xyyyz_zzzz, ts_xyyyz_zzzzz, ts_yyyz_xxxxx, ts_yyyz_xxxxy, ts_yyyz_xxxxz, ts_yyyz_xxxyy, ts_yyyz_xxxyz, ts_yyyz_xxxzz, ts_yyyz_xxyyy, ts_yyyz_xxyyz, ts_yyyz_xxyzz, ts_yyyz_xxzzz, ts_yyyz_xyyyy, ts_yyyz_xyyyz, ts_yyyz_xyyzz, ts_yyyz_xyzzz, ts_yyyz_xzzzz, ts_yyyz_yyyyy, ts_yyyz_yyyyz, ts_yyyz_yyyzz, ts_yyyz_yyzzz, ts_yyyz_yzzzz, ts_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyz_xxxxx[i] = 2.0 * ts_yyyz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxxxy[i] = 2.0 * ts_yyyz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxxxz[i] = 2.0 * ts_yyyz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxxyy[i] = 2.0 * ts_yyyz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxxyz[i] = 2.0 * ts_yyyz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxxzz[i] = 2.0 * ts_yyyz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxyyy[i] = 2.0 * ts_yyyz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxyyz[i] = 2.0 * ts_yyyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxyzz[i] = 2.0 * ts_yyyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xxzzz[i] = 2.0 * ts_yyyz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyyyy[i] = 2.0 * ts_yyyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyyyz[i] = 2.0 * ts_yyyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyyzz[i] = 2.0 * ts_yyyz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xyzzz[i] = 2.0 * ts_yyyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xzzzz[i] = 2.0 * ts_yyyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyyyy[i] = 2.0 * ts_yyyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyyyz[i] = 2.0 * ts_yyyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyyzz[i] = 2.0 * ts_yyyz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yyzzz[i] = 2.0 * ts_yyyz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yzzzz[i] = 2.0 * ts_yyyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_zzzzz[i] = 2.0 * ts_yyyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 252-273 components of targeted buffer : HH

    auto gs_x_xyyzz_xxxxx = pbuffer.data(idx_g_hh + 252);

    auto gs_x_xyyzz_xxxxy = pbuffer.data(idx_g_hh + 253);

    auto gs_x_xyyzz_xxxxz = pbuffer.data(idx_g_hh + 254);

    auto gs_x_xyyzz_xxxyy = pbuffer.data(idx_g_hh + 255);

    auto gs_x_xyyzz_xxxyz = pbuffer.data(idx_g_hh + 256);

    auto gs_x_xyyzz_xxxzz = pbuffer.data(idx_g_hh + 257);

    auto gs_x_xyyzz_xxyyy = pbuffer.data(idx_g_hh + 258);

    auto gs_x_xyyzz_xxyyz = pbuffer.data(idx_g_hh + 259);

    auto gs_x_xyyzz_xxyzz = pbuffer.data(idx_g_hh + 260);

    auto gs_x_xyyzz_xxzzz = pbuffer.data(idx_g_hh + 261);

    auto gs_x_xyyzz_xyyyy = pbuffer.data(idx_g_hh + 262);

    auto gs_x_xyyzz_xyyyz = pbuffer.data(idx_g_hh + 263);

    auto gs_x_xyyzz_xyyzz = pbuffer.data(idx_g_hh + 264);

    auto gs_x_xyyzz_xyzzz = pbuffer.data(idx_g_hh + 265);

    auto gs_x_xyyzz_xzzzz = pbuffer.data(idx_g_hh + 266);

    auto gs_x_xyyzz_yyyyy = pbuffer.data(idx_g_hh + 267);

    auto gs_x_xyyzz_yyyyz = pbuffer.data(idx_g_hh + 268);

    auto gs_x_xyyzz_yyyzz = pbuffer.data(idx_g_hh + 269);

    auto gs_x_xyyzz_yyzzz = pbuffer.data(idx_g_hh + 270);

    auto gs_x_xyyzz_yzzzz = pbuffer.data(idx_g_hh + 271);

    auto gs_x_xyyzz_zzzzz = pbuffer.data(idx_g_hh + 272);

    #pragma omp simd aligned(gc_x, gs_x_xyyzz_xxxxx, gs_x_xyyzz_xxxxy, gs_x_xyyzz_xxxxz, gs_x_xyyzz_xxxyy, gs_x_xyyzz_xxxyz, gs_x_xyyzz_xxxzz, gs_x_xyyzz_xxyyy, gs_x_xyyzz_xxyyz, gs_x_xyyzz_xxyzz, gs_x_xyyzz_xxzzz, gs_x_xyyzz_xyyyy, gs_x_xyyzz_xyyyz, gs_x_xyyzz_xyyzz, gs_x_xyyzz_xyzzz, gs_x_xyyzz_xzzzz, gs_x_xyyzz_yyyyy, gs_x_xyyzz_yyyyz, gs_x_xyyzz_yyyzz, gs_x_xyyzz_yyzzz, gs_x_xyyzz_yzzzz, gs_x_xyyzz_zzzzz, ts_xyyzz_xxxx, ts_xyyzz_xxxxx, ts_xyyzz_xxxxy, ts_xyyzz_xxxxz, ts_xyyzz_xxxy, ts_xyyzz_xxxyy, ts_xyyzz_xxxyz, ts_xyyzz_xxxz, ts_xyyzz_xxxzz, ts_xyyzz_xxyy, ts_xyyzz_xxyyy, ts_xyyzz_xxyyz, ts_xyyzz_xxyz, ts_xyyzz_xxyzz, ts_xyyzz_xxzz, ts_xyyzz_xxzzz, ts_xyyzz_xyyy, ts_xyyzz_xyyyy, ts_xyyzz_xyyyz, ts_xyyzz_xyyz, ts_xyyzz_xyyzz, ts_xyyzz_xyzz, ts_xyyzz_xyzzz, ts_xyyzz_xzzz, ts_xyyzz_xzzzz, ts_xyyzz_yyyy, ts_xyyzz_yyyyy, ts_xyyzz_yyyyz, ts_xyyzz_yyyz, ts_xyyzz_yyyzz, ts_xyyzz_yyzz, ts_xyyzz_yyzzz, ts_xyyzz_yzzz, ts_xyyzz_yzzzz, ts_xyyzz_zzzz, ts_xyyzz_zzzzz, ts_yyzz_xxxxx, ts_yyzz_xxxxy, ts_yyzz_xxxxz, ts_yyzz_xxxyy, ts_yyzz_xxxyz, ts_yyzz_xxxzz, ts_yyzz_xxyyy, ts_yyzz_xxyyz, ts_yyzz_xxyzz, ts_yyzz_xxzzz, ts_yyzz_xyyyy, ts_yyzz_xyyyz, ts_yyzz_xyyzz, ts_yyzz_xyzzz, ts_yyzz_xzzzz, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyzz, ts_yyzz_yyzzz, ts_yyzz_yzzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyzz_xxxxx[i] = 2.0 * ts_yyzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxxxy[i] = 2.0 * ts_yyzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxxxz[i] = 2.0 * ts_yyzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxxyy[i] = 2.0 * ts_yyzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxxyz[i] = 2.0 * ts_yyzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxxzz[i] = 2.0 * ts_yyzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxyyy[i] = 2.0 * ts_yyzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxyyz[i] = 2.0 * ts_yyzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxyzz[i] = 2.0 * ts_yyzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xxzzz[i] = 2.0 * ts_yyzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyyyy[i] = 2.0 * ts_yyzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyyyz[i] = 2.0 * ts_yyzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyyzz[i] = 2.0 * ts_yyzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xyzzz[i] = 2.0 * ts_yyzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xzzzz[i] = 2.0 * ts_yyzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyyyy[i] = 2.0 * ts_yyzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyyyz[i] = 2.0 * ts_yyzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyyzz[i] = 2.0 * ts_yyzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yyzzz[i] = 2.0 * ts_yyzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yzzzz[i] = 2.0 * ts_yyzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_zzzzz[i] = 2.0 * ts_yyzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 273-294 components of targeted buffer : HH

    auto gs_x_xyzzz_xxxxx = pbuffer.data(idx_g_hh + 273);

    auto gs_x_xyzzz_xxxxy = pbuffer.data(idx_g_hh + 274);

    auto gs_x_xyzzz_xxxxz = pbuffer.data(idx_g_hh + 275);

    auto gs_x_xyzzz_xxxyy = pbuffer.data(idx_g_hh + 276);

    auto gs_x_xyzzz_xxxyz = pbuffer.data(idx_g_hh + 277);

    auto gs_x_xyzzz_xxxzz = pbuffer.data(idx_g_hh + 278);

    auto gs_x_xyzzz_xxyyy = pbuffer.data(idx_g_hh + 279);

    auto gs_x_xyzzz_xxyyz = pbuffer.data(idx_g_hh + 280);

    auto gs_x_xyzzz_xxyzz = pbuffer.data(idx_g_hh + 281);

    auto gs_x_xyzzz_xxzzz = pbuffer.data(idx_g_hh + 282);

    auto gs_x_xyzzz_xyyyy = pbuffer.data(idx_g_hh + 283);

    auto gs_x_xyzzz_xyyyz = pbuffer.data(idx_g_hh + 284);

    auto gs_x_xyzzz_xyyzz = pbuffer.data(idx_g_hh + 285);

    auto gs_x_xyzzz_xyzzz = pbuffer.data(idx_g_hh + 286);

    auto gs_x_xyzzz_xzzzz = pbuffer.data(idx_g_hh + 287);

    auto gs_x_xyzzz_yyyyy = pbuffer.data(idx_g_hh + 288);

    auto gs_x_xyzzz_yyyyz = pbuffer.data(idx_g_hh + 289);

    auto gs_x_xyzzz_yyyzz = pbuffer.data(idx_g_hh + 290);

    auto gs_x_xyzzz_yyzzz = pbuffer.data(idx_g_hh + 291);

    auto gs_x_xyzzz_yzzzz = pbuffer.data(idx_g_hh + 292);

    auto gs_x_xyzzz_zzzzz = pbuffer.data(idx_g_hh + 293);

    #pragma omp simd aligned(gc_x, gs_x_xyzzz_xxxxx, gs_x_xyzzz_xxxxy, gs_x_xyzzz_xxxxz, gs_x_xyzzz_xxxyy, gs_x_xyzzz_xxxyz, gs_x_xyzzz_xxxzz, gs_x_xyzzz_xxyyy, gs_x_xyzzz_xxyyz, gs_x_xyzzz_xxyzz, gs_x_xyzzz_xxzzz, gs_x_xyzzz_xyyyy, gs_x_xyzzz_xyyyz, gs_x_xyzzz_xyyzz, gs_x_xyzzz_xyzzz, gs_x_xyzzz_xzzzz, gs_x_xyzzz_yyyyy, gs_x_xyzzz_yyyyz, gs_x_xyzzz_yyyzz, gs_x_xyzzz_yyzzz, gs_x_xyzzz_yzzzz, gs_x_xyzzz_zzzzz, ts_xyzzz_xxxx, ts_xyzzz_xxxxx, ts_xyzzz_xxxxy, ts_xyzzz_xxxxz, ts_xyzzz_xxxy, ts_xyzzz_xxxyy, ts_xyzzz_xxxyz, ts_xyzzz_xxxz, ts_xyzzz_xxxzz, ts_xyzzz_xxyy, ts_xyzzz_xxyyy, ts_xyzzz_xxyyz, ts_xyzzz_xxyz, ts_xyzzz_xxyzz, ts_xyzzz_xxzz, ts_xyzzz_xxzzz, ts_xyzzz_xyyy, ts_xyzzz_xyyyy, ts_xyzzz_xyyyz, ts_xyzzz_xyyz, ts_xyzzz_xyyzz, ts_xyzzz_xyzz, ts_xyzzz_xyzzz, ts_xyzzz_xzzz, ts_xyzzz_xzzzz, ts_xyzzz_yyyy, ts_xyzzz_yyyyy, ts_xyzzz_yyyyz, ts_xyzzz_yyyz, ts_xyzzz_yyyzz, ts_xyzzz_yyzz, ts_xyzzz_yyzzz, ts_xyzzz_yzzz, ts_xyzzz_yzzzz, ts_xyzzz_zzzz, ts_xyzzz_zzzzz, ts_yzzz_xxxxx, ts_yzzz_xxxxy, ts_yzzz_xxxxz, ts_yzzz_xxxyy, ts_yzzz_xxxyz, ts_yzzz_xxxzz, ts_yzzz_xxyyy, ts_yzzz_xxyyz, ts_yzzz_xxyzz, ts_yzzz_xxzzz, ts_yzzz_xyyyy, ts_yzzz_xyyyz, ts_yzzz_xyyzz, ts_yzzz_xyzzz, ts_yzzz_xzzzz, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyzz, ts_yzzz_yyzzz, ts_yzzz_yzzzz, ts_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzzz_xxxxx[i] = 2.0 * ts_yzzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxxxy[i] = 2.0 * ts_yzzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxxxz[i] = 2.0 * ts_yzzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxxyy[i] = 2.0 * ts_yzzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxxyz[i] = 2.0 * ts_yzzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxxzz[i] = 2.0 * ts_yzzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxyyy[i] = 2.0 * ts_yzzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxyyz[i] = 2.0 * ts_yzzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxyzz[i] = 2.0 * ts_yzzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xxzzz[i] = 2.0 * ts_yzzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyyyy[i] = 2.0 * ts_yzzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyyyz[i] = 2.0 * ts_yzzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyyzz[i] = 2.0 * ts_yzzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xyzzz[i] = 2.0 * ts_yzzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xzzzz[i] = 2.0 * ts_yzzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyyyy[i] = 2.0 * ts_yzzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyyyz[i] = 2.0 * ts_yzzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyyzz[i] = 2.0 * ts_yzzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yyzzz[i] = 2.0 * ts_yzzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yzzzz[i] = 2.0 * ts_yzzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_zzzzz[i] = 2.0 * ts_yzzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 294-315 components of targeted buffer : HH

    auto gs_x_xzzzz_xxxxx = pbuffer.data(idx_g_hh + 294);

    auto gs_x_xzzzz_xxxxy = pbuffer.data(idx_g_hh + 295);

    auto gs_x_xzzzz_xxxxz = pbuffer.data(idx_g_hh + 296);

    auto gs_x_xzzzz_xxxyy = pbuffer.data(idx_g_hh + 297);

    auto gs_x_xzzzz_xxxyz = pbuffer.data(idx_g_hh + 298);

    auto gs_x_xzzzz_xxxzz = pbuffer.data(idx_g_hh + 299);

    auto gs_x_xzzzz_xxyyy = pbuffer.data(idx_g_hh + 300);

    auto gs_x_xzzzz_xxyyz = pbuffer.data(idx_g_hh + 301);

    auto gs_x_xzzzz_xxyzz = pbuffer.data(idx_g_hh + 302);

    auto gs_x_xzzzz_xxzzz = pbuffer.data(idx_g_hh + 303);

    auto gs_x_xzzzz_xyyyy = pbuffer.data(idx_g_hh + 304);

    auto gs_x_xzzzz_xyyyz = pbuffer.data(idx_g_hh + 305);

    auto gs_x_xzzzz_xyyzz = pbuffer.data(idx_g_hh + 306);

    auto gs_x_xzzzz_xyzzz = pbuffer.data(idx_g_hh + 307);

    auto gs_x_xzzzz_xzzzz = pbuffer.data(idx_g_hh + 308);

    auto gs_x_xzzzz_yyyyy = pbuffer.data(idx_g_hh + 309);

    auto gs_x_xzzzz_yyyyz = pbuffer.data(idx_g_hh + 310);

    auto gs_x_xzzzz_yyyzz = pbuffer.data(idx_g_hh + 311);

    auto gs_x_xzzzz_yyzzz = pbuffer.data(idx_g_hh + 312);

    auto gs_x_xzzzz_yzzzz = pbuffer.data(idx_g_hh + 313);

    auto gs_x_xzzzz_zzzzz = pbuffer.data(idx_g_hh + 314);

    #pragma omp simd aligned(gc_x, gs_x_xzzzz_xxxxx, gs_x_xzzzz_xxxxy, gs_x_xzzzz_xxxxz, gs_x_xzzzz_xxxyy, gs_x_xzzzz_xxxyz, gs_x_xzzzz_xxxzz, gs_x_xzzzz_xxyyy, gs_x_xzzzz_xxyyz, gs_x_xzzzz_xxyzz, gs_x_xzzzz_xxzzz, gs_x_xzzzz_xyyyy, gs_x_xzzzz_xyyyz, gs_x_xzzzz_xyyzz, gs_x_xzzzz_xyzzz, gs_x_xzzzz_xzzzz, gs_x_xzzzz_yyyyy, gs_x_xzzzz_yyyyz, gs_x_xzzzz_yyyzz, gs_x_xzzzz_yyzzz, gs_x_xzzzz_yzzzz, gs_x_xzzzz_zzzzz, ts_xzzzz_xxxx, ts_xzzzz_xxxxx, ts_xzzzz_xxxxy, ts_xzzzz_xxxxz, ts_xzzzz_xxxy, ts_xzzzz_xxxyy, ts_xzzzz_xxxyz, ts_xzzzz_xxxz, ts_xzzzz_xxxzz, ts_xzzzz_xxyy, ts_xzzzz_xxyyy, ts_xzzzz_xxyyz, ts_xzzzz_xxyz, ts_xzzzz_xxyzz, ts_xzzzz_xxzz, ts_xzzzz_xxzzz, ts_xzzzz_xyyy, ts_xzzzz_xyyyy, ts_xzzzz_xyyyz, ts_xzzzz_xyyz, ts_xzzzz_xyyzz, ts_xzzzz_xyzz, ts_xzzzz_xyzzz, ts_xzzzz_xzzz, ts_xzzzz_xzzzz, ts_xzzzz_yyyy, ts_xzzzz_yyyyy, ts_xzzzz_yyyyz, ts_xzzzz_yyyz, ts_xzzzz_yyyzz, ts_xzzzz_yyzz, ts_xzzzz_yyzzz, ts_xzzzz_yzzz, ts_xzzzz_yzzzz, ts_xzzzz_zzzz, ts_xzzzz_zzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzzz_xxxxx[i] = 2.0 * ts_zzzz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxxxy[i] = 2.0 * ts_zzzz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxxxz[i] = 2.0 * ts_zzzz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxxyy[i] = 2.0 * ts_zzzz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxxyz[i] = 2.0 * ts_zzzz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxxzz[i] = 2.0 * ts_zzzz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxyyy[i] = 2.0 * ts_zzzz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxyyz[i] = 2.0 * ts_zzzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxyzz[i] = 2.0 * ts_zzzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xxzzz[i] = 2.0 * ts_zzzz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyyyy[i] = 2.0 * ts_zzzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyyyz[i] = 2.0 * ts_zzzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyyzz[i] = 2.0 * ts_zzzz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xyzzz[i] = 2.0 * ts_zzzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xzzzz[i] = 2.0 * ts_zzzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyyyy[i] = 2.0 * ts_zzzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyyyz[i] = 2.0 * ts_zzzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyyzz[i] = 2.0 * ts_zzzz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yyzzz[i] = 2.0 * ts_zzzz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yzzzz[i] = 2.0 * ts_zzzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_zzzzz[i] = 2.0 * ts_zzzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 315-336 components of targeted buffer : HH

    auto gs_x_yyyyy_xxxxx = pbuffer.data(idx_g_hh + 315);

    auto gs_x_yyyyy_xxxxy = pbuffer.data(idx_g_hh + 316);

    auto gs_x_yyyyy_xxxxz = pbuffer.data(idx_g_hh + 317);

    auto gs_x_yyyyy_xxxyy = pbuffer.data(idx_g_hh + 318);

    auto gs_x_yyyyy_xxxyz = pbuffer.data(idx_g_hh + 319);

    auto gs_x_yyyyy_xxxzz = pbuffer.data(idx_g_hh + 320);

    auto gs_x_yyyyy_xxyyy = pbuffer.data(idx_g_hh + 321);

    auto gs_x_yyyyy_xxyyz = pbuffer.data(idx_g_hh + 322);

    auto gs_x_yyyyy_xxyzz = pbuffer.data(idx_g_hh + 323);

    auto gs_x_yyyyy_xxzzz = pbuffer.data(idx_g_hh + 324);

    auto gs_x_yyyyy_xyyyy = pbuffer.data(idx_g_hh + 325);

    auto gs_x_yyyyy_xyyyz = pbuffer.data(idx_g_hh + 326);

    auto gs_x_yyyyy_xyyzz = pbuffer.data(idx_g_hh + 327);

    auto gs_x_yyyyy_xyzzz = pbuffer.data(idx_g_hh + 328);

    auto gs_x_yyyyy_xzzzz = pbuffer.data(idx_g_hh + 329);

    auto gs_x_yyyyy_yyyyy = pbuffer.data(idx_g_hh + 330);

    auto gs_x_yyyyy_yyyyz = pbuffer.data(idx_g_hh + 331);

    auto gs_x_yyyyy_yyyzz = pbuffer.data(idx_g_hh + 332);

    auto gs_x_yyyyy_yyzzz = pbuffer.data(idx_g_hh + 333);

    auto gs_x_yyyyy_yzzzz = pbuffer.data(idx_g_hh + 334);

    auto gs_x_yyyyy_zzzzz = pbuffer.data(idx_g_hh + 335);

    #pragma omp simd aligned(gc_x, gs_x_yyyyy_xxxxx, gs_x_yyyyy_xxxxy, gs_x_yyyyy_xxxxz, gs_x_yyyyy_xxxyy, gs_x_yyyyy_xxxyz, gs_x_yyyyy_xxxzz, gs_x_yyyyy_xxyyy, gs_x_yyyyy_xxyyz, gs_x_yyyyy_xxyzz, gs_x_yyyyy_xxzzz, gs_x_yyyyy_xyyyy, gs_x_yyyyy_xyyyz, gs_x_yyyyy_xyyzz, gs_x_yyyyy_xyzzz, gs_x_yyyyy_xzzzz, gs_x_yyyyy_yyyyy, gs_x_yyyyy_yyyyz, gs_x_yyyyy_yyyzz, gs_x_yyyyy_yyzzz, gs_x_yyyyy_yzzzz, gs_x_yyyyy_zzzzz, ts_yyyyy_xxxx, ts_yyyyy_xxxxx, ts_yyyyy_xxxxy, ts_yyyyy_xxxxz, ts_yyyyy_xxxy, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxxz, ts_yyyyy_xxxzz, ts_yyyyy_xxyy, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyz, ts_yyyyy_xxyzz, ts_yyyyy_xxzz, ts_yyyyy_xxzzz, ts_yyyyy_xyyy, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzz, ts_yyyyy_xyzzz, ts_yyyyy_xzzz, ts_yyyyy_xzzzz, ts_yyyyy_yyyy, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzz, ts_yyyyy_yzzzz, ts_yyyyy_zzzz, ts_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyy_xxxxx[i] = 10.0 * ts_yyyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxxxy[i] = 8.0 * ts_yyyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxxxz[i] = 8.0 * ts_yyyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxxyy[i] = 6.0 * ts_yyyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxxyz[i] = 6.0 * ts_yyyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxxzz[i] = 6.0 * ts_yyyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxyyy[i] = 4.0 * ts_yyyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxyyz[i] = 4.0 * ts_yyyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxyzz[i] = 4.0 * ts_yyyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xxzzz[i] = 4.0 * ts_yyyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyyyy[i] = 2.0 * ts_yyyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyyyz[i] = 2.0 * ts_yyyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyyzz[i] = 2.0 * ts_yyyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xyzzz[i] = 2.0 * ts_yyyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xzzzz[i] = 2.0 * ts_yyyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyyyy[i] = 2.0 * ts_yyyyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyyyz[i] = 2.0 * ts_yyyyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyyzz[i] = 2.0 * ts_yyyyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yyzzz[i] = 2.0 * ts_yyyyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yzzzz[i] = 2.0 * ts_yyyyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_zzzzz[i] = 2.0 * ts_yyyyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 336-357 components of targeted buffer : HH

    auto gs_x_yyyyz_xxxxx = pbuffer.data(idx_g_hh + 336);

    auto gs_x_yyyyz_xxxxy = pbuffer.data(idx_g_hh + 337);

    auto gs_x_yyyyz_xxxxz = pbuffer.data(idx_g_hh + 338);

    auto gs_x_yyyyz_xxxyy = pbuffer.data(idx_g_hh + 339);

    auto gs_x_yyyyz_xxxyz = pbuffer.data(idx_g_hh + 340);

    auto gs_x_yyyyz_xxxzz = pbuffer.data(idx_g_hh + 341);

    auto gs_x_yyyyz_xxyyy = pbuffer.data(idx_g_hh + 342);

    auto gs_x_yyyyz_xxyyz = pbuffer.data(idx_g_hh + 343);

    auto gs_x_yyyyz_xxyzz = pbuffer.data(idx_g_hh + 344);

    auto gs_x_yyyyz_xxzzz = pbuffer.data(idx_g_hh + 345);

    auto gs_x_yyyyz_xyyyy = pbuffer.data(idx_g_hh + 346);

    auto gs_x_yyyyz_xyyyz = pbuffer.data(idx_g_hh + 347);

    auto gs_x_yyyyz_xyyzz = pbuffer.data(idx_g_hh + 348);

    auto gs_x_yyyyz_xyzzz = pbuffer.data(idx_g_hh + 349);

    auto gs_x_yyyyz_xzzzz = pbuffer.data(idx_g_hh + 350);

    auto gs_x_yyyyz_yyyyy = pbuffer.data(idx_g_hh + 351);

    auto gs_x_yyyyz_yyyyz = pbuffer.data(idx_g_hh + 352);

    auto gs_x_yyyyz_yyyzz = pbuffer.data(idx_g_hh + 353);

    auto gs_x_yyyyz_yyzzz = pbuffer.data(idx_g_hh + 354);

    auto gs_x_yyyyz_yzzzz = pbuffer.data(idx_g_hh + 355);

    auto gs_x_yyyyz_zzzzz = pbuffer.data(idx_g_hh + 356);

    #pragma omp simd aligned(gc_x, gs_x_yyyyz_xxxxx, gs_x_yyyyz_xxxxy, gs_x_yyyyz_xxxxz, gs_x_yyyyz_xxxyy, gs_x_yyyyz_xxxyz, gs_x_yyyyz_xxxzz, gs_x_yyyyz_xxyyy, gs_x_yyyyz_xxyyz, gs_x_yyyyz_xxyzz, gs_x_yyyyz_xxzzz, gs_x_yyyyz_xyyyy, gs_x_yyyyz_xyyyz, gs_x_yyyyz_xyyzz, gs_x_yyyyz_xyzzz, gs_x_yyyyz_xzzzz, gs_x_yyyyz_yyyyy, gs_x_yyyyz_yyyyz, gs_x_yyyyz_yyyzz, gs_x_yyyyz_yyzzz, gs_x_yyyyz_yzzzz, gs_x_yyyyz_zzzzz, ts_yyyyz_xxxx, ts_yyyyz_xxxxx, ts_yyyyz_xxxxy, ts_yyyyz_xxxxz, ts_yyyyz_xxxy, ts_yyyyz_xxxyy, ts_yyyyz_xxxyz, ts_yyyyz_xxxz, ts_yyyyz_xxxzz, ts_yyyyz_xxyy, ts_yyyyz_xxyyy, ts_yyyyz_xxyyz, ts_yyyyz_xxyz, ts_yyyyz_xxyzz, ts_yyyyz_xxzz, ts_yyyyz_xxzzz, ts_yyyyz_xyyy, ts_yyyyz_xyyyy, ts_yyyyz_xyyyz, ts_yyyyz_xyyz, ts_yyyyz_xyyzz, ts_yyyyz_xyzz, ts_yyyyz_xyzzz, ts_yyyyz_xzzz, ts_yyyyz_xzzzz, ts_yyyyz_yyyy, ts_yyyyz_yyyyy, ts_yyyyz_yyyyz, ts_yyyyz_yyyz, ts_yyyyz_yyyzz, ts_yyyyz_yyzz, ts_yyyyz_yyzzz, ts_yyyyz_yzzz, ts_yyyyz_yzzzz, ts_yyyyz_zzzz, ts_yyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyz_xxxxx[i] = 10.0 * ts_yyyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxxxy[i] = 8.0 * ts_yyyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxxxz[i] = 8.0 * ts_yyyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxxyy[i] = 6.0 * ts_yyyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxxyz[i] = 6.0 * ts_yyyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxxzz[i] = 6.0 * ts_yyyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxyyy[i] = 4.0 * ts_yyyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxyyz[i] = 4.0 * ts_yyyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxyzz[i] = 4.0 * ts_yyyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xxzzz[i] = 4.0 * ts_yyyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyyyy[i] = 2.0 * ts_yyyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyyyz[i] = 2.0 * ts_yyyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyyzz[i] = 2.0 * ts_yyyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xyzzz[i] = 2.0 * ts_yyyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xzzzz[i] = 2.0 * ts_yyyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyyyy[i] = 2.0 * ts_yyyyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyyyz[i] = 2.0 * ts_yyyyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyyzz[i] = 2.0 * ts_yyyyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yyzzz[i] = 2.0 * ts_yyyyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yzzzz[i] = 2.0 * ts_yyyyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_zzzzz[i] = 2.0 * ts_yyyyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 357-378 components of targeted buffer : HH

    auto gs_x_yyyzz_xxxxx = pbuffer.data(idx_g_hh + 357);

    auto gs_x_yyyzz_xxxxy = pbuffer.data(idx_g_hh + 358);

    auto gs_x_yyyzz_xxxxz = pbuffer.data(idx_g_hh + 359);

    auto gs_x_yyyzz_xxxyy = pbuffer.data(idx_g_hh + 360);

    auto gs_x_yyyzz_xxxyz = pbuffer.data(idx_g_hh + 361);

    auto gs_x_yyyzz_xxxzz = pbuffer.data(idx_g_hh + 362);

    auto gs_x_yyyzz_xxyyy = pbuffer.data(idx_g_hh + 363);

    auto gs_x_yyyzz_xxyyz = pbuffer.data(idx_g_hh + 364);

    auto gs_x_yyyzz_xxyzz = pbuffer.data(idx_g_hh + 365);

    auto gs_x_yyyzz_xxzzz = pbuffer.data(idx_g_hh + 366);

    auto gs_x_yyyzz_xyyyy = pbuffer.data(idx_g_hh + 367);

    auto gs_x_yyyzz_xyyyz = pbuffer.data(idx_g_hh + 368);

    auto gs_x_yyyzz_xyyzz = pbuffer.data(idx_g_hh + 369);

    auto gs_x_yyyzz_xyzzz = pbuffer.data(idx_g_hh + 370);

    auto gs_x_yyyzz_xzzzz = pbuffer.data(idx_g_hh + 371);

    auto gs_x_yyyzz_yyyyy = pbuffer.data(idx_g_hh + 372);

    auto gs_x_yyyzz_yyyyz = pbuffer.data(idx_g_hh + 373);

    auto gs_x_yyyzz_yyyzz = pbuffer.data(idx_g_hh + 374);

    auto gs_x_yyyzz_yyzzz = pbuffer.data(idx_g_hh + 375);

    auto gs_x_yyyzz_yzzzz = pbuffer.data(idx_g_hh + 376);

    auto gs_x_yyyzz_zzzzz = pbuffer.data(idx_g_hh + 377);

    #pragma omp simd aligned(gc_x, gs_x_yyyzz_xxxxx, gs_x_yyyzz_xxxxy, gs_x_yyyzz_xxxxz, gs_x_yyyzz_xxxyy, gs_x_yyyzz_xxxyz, gs_x_yyyzz_xxxzz, gs_x_yyyzz_xxyyy, gs_x_yyyzz_xxyyz, gs_x_yyyzz_xxyzz, gs_x_yyyzz_xxzzz, gs_x_yyyzz_xyyyy, gs_x_yyyzz_xyyyz, gs_x_yyyzz_xyyzz, gs_x_yyyzz_xyzzz, gs_x_yyyzz_xzzzz, gs_x_yyyzz_yyyyy, gs_x_yyyzz_yyyyz, gs_x_yyyzz_yyyzz, gs_x_yyyzz_yyzzz, gs_x_yyyzz_yzzzz, gs_x_yyyzz_zzzzz, ts_yyyzz_xxxx, ts_yyyzz_xxxxx, ts_yyyzz_xxxxy, ts_yyyzz_xxxxz, ts_yyyzz_xxxy, ts_yyyzz_xxxyy, ts_yyyzz_xxxyz, ts_yyyzz_xxxz, ts_yyyzz_xxxzz, ts_yyyzz_xxyy, ts_yyyzz_xxyyy, ts_yyyzz_xxyyz, ts_yyyzz_xxyz, ts_yyyzz_xxyzz, ts_yyyzz_xxzz, ts_yyyzz_xxzzz, ts_yyyzz_xyyy, ts_yyyzz_xyyyy, ts_yyyzz_xyyyz, ts_yyyzz_xyyz, ts_yyyzz_xyyzz, ts_yyyzz_xyzz, ts_yyyzz_xyzzz, ts_yyyzz_xzzz, ts_yyyzz_xzzzz, ts_yyyzz_yyyy, ts_yyyzz_yyyyy, ts_yyyzz_yyyyz, ts_yyyzz_yyyz, ts_yyyzz_yyyzz, ts_yyyzz_yyzz, ts_yyyzz_yyzzz, ts_yyyzz_yzzz, ts_yyyzz_yzzzz, ts_yyyzz_zzzz, ts_yyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyzz_xxxxx[i] = 10.0 * ts_yyyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxxxy[i] = 8.0 * ts_yyyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxxxz[i] = 8.0 * ts_yyyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxxyy[i] = 6.0 * ts_yyyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxxyz[i] = 6.0 * ts_yyyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxxzz[i] = 6.0 * ts_yyyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxyyy[i] = 4.0 * ts_yyyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxyyz[i] = 4.0 * ts_yyyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxyzz[i] = 4.0 * ts_yyyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xxzzz[i] = 4.0 * ts_yyyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyyyy[i] = 2.0 * ts_yyyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyyyz[i] = 2.0 * ts_yyyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyyzz[i] = 2.0 * ts_yyyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xyzzz[i] = 2.0 * ts_yyyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xzzzz[i] = 2.0 * ts_yyyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyyyy[i] = 2.0 * ts_yyyzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyyyz[i] = 2.0 * ts_yyyzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyyzz[i] = 2.0 * ts_yyyzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yyzzz[i] = 2.0 * ts_yyyzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yzzzz[i] = 2.0 * ts_yyyzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_zzzzz[i] = 2.0 * ts_yyyzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 378-399 components of targeted buffer : HH

    auto gs_x_yyzzz_xxxxx = pbuffer.data(idx_g_hh + 378);

    auto gs_x_yyzzz_xxxxy = pbuffer.data(idx_g_hh + 379);

    auto gs_x_yyzzz_xxxxz = pbuffer.data(idx_g_hh + 380);

    auto gs_x_yyzzz_xxxyy = pbuffer.data(idx_g_hh + 381);

    auto gs_x_yyzzz_xxxyz = pbuffer.data(idx_g_hh + 382);

    auto gs_x_yyzzz_xxxzz = pbuffer.data(idx_g_hh + 383);

    auto gs_x_yyzzz_xxyyy = pbuffer.data(idx_g_hh + 384);

    auto gs_x_yyzzz_xxyyz = pbuffer.data(idx_g_hh + 385);

    auto gs_x_yyzzz_xxyzz = pbuffer.data(idx_g_hh + 386);

    auto gs_x_yyzzz_xxzzz = pbuffer.data(idx_g_hh + 387);

    auto gs_x_yyzzz_xyyyy = pbuffer.data(idx_g_hh + 388);

    auto gs_x_yyzzz_xyyyz = pbuffer.data(idx_g_hh + 389);

    auto gs_x_yyzzz_xyyzz = pbuffer.data(idx_g_hh + 390);

    auto gs_x_yyzzz_xyzzz = pbuffer.data(idx_g_hh + 391);

    auto gs_x_yyzzz_xzzzz = pbuffer.data(idx_g_hh + 392);

    auto gs_x_yyzzz_yyyyy = pbuffer.data(idx_g_hh + 393);

    auto gs_x_yyzzz_yyyyz = pbuffer.data(idx_g_hh + 394);

    auto gs_x_yyzzz_yyyzz = pbuffer.data(idx_g_hh + 395);

    auto gs_x_yyzzz_yyzzz = pbuffer.data(idx_g_hh + 396);

    auto gs_x_yyzzz_yzzzz = pbuffer.data(idx_g_hh + 397);

    auto gs_x_yyzzz_zzzzz = pbuffer.data(idx_g_hh + 398);

    #pragma omp simd aligned(gc_x, gs_x_yyzzz_xxxxx, gs_x_yyzzz_xxxxy, gs_x_yyzzz_xxxxz, gs_x_yyzzz_xxxyy, gs_x_yyzzz_xxxyz, gs_x_yyzzz_xxxzz, gs_x_yyzzz_xxyyy, gs_x_yyzzz_xxyyz, gs_x_yyzzz_xxyzz, gs_x_yyzzz_xxzzz, gs_x_yyzzz_xyyyy, gs_x_yyzzz_xyyyz, gs_x_yyzzz_xyyzz, gs_x_yyzzz_xyzzz, gs_x_yyzzz_xzzzz, gs_x_yyzzz_yyyyy, gs_x_yyzzz_yyyyz, gs_x_yyzzz_yyyzz, gs_x_yyzzz_yyzzz, gs_x_yyzzz_yzzzz, gs_x_yyzzz_zzzzz, ts_yyzzz_xxxx, ts_yyzzz_xxxxx, ts_yyzzz_xxxxy, ts_yyzzz_xxxxz, ts_yyzzz_xxxy, ts_yyzzz_xxxyy, ts_yyzzz_xxxyz, ts_yyzzz_xxxz, ts_yyzzz_xxxzz, ts_yyzzz_xxyy, ts_yyzzz_xxyyy, ts_yyzzz_xxyyz, ts_yyzzz_xxyz, ts_yyzzz_xxyzz, ts_yyzzz_xxzz, ts_yyzzz_xxzzz, ts_yyzzz_xyyy, ts_yyzzz_xyyyy, ts_yyzzz_xyyyz, ts_yyzzz_xyyz, ts_yyzzz_xyyzz, ts_yyzzz_xyzz, ts_yyzzz_xyzzz, ts_yyzzz_xzzz, ts_yyzzz_xzzzz, ts_yyzzz_yyyy, ts_yyzzz_yyyyy, ts_yyzzz_yyyyz, ts_yyzzz_yyyz, ts_yyzzz_yyyzz, ts_yyzzz_yyzz, ts_yyzzz_yyzzz, ts_yyzzz_yzzz, ts_yyzzz_yzzzz, ts_yyzzz_zzzz, ts_yyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzzz_xxxxx[i] = 10.0 * ts_yyzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxxxy[i] = 8.0 * ts_yyzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxxxz[i] = 8.0 * ts_yyzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxxyy[i] = 6.0 * ts_yyzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxxyz[i] = 6.0 * ts_yyzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxxzz[i] = 6.0 * ts_yyzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxyyy[i] = 4.0 * ts_yyzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxyyz[i] = 4.0 * ts_yyzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxyzz[i] = 4.0 * ts_yyzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xxzzz[i] = 4.0 * ts_yyzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyyyy[i] = 2.0 * ts_yyzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyyyz[i] = 2.0 * ts_yyzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyyzz[i] = 2.0 * ts_yyzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xyzzz[i] = 2.0 * ts_yyzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xzzzz[i] = 2.0 * ts_yyzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyyyy[i] = 2.0 * ts_yyzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyyyz[i] = 2.0 * ts_yyzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyyzz[i] = 2.0 * ts_yyzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yyzzz[i] = 2.0 * ts_yyzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yzzzz[i] = 2.0 * ts_yyzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_zzzzz[i] = 2.0 * ts_yyzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 399-420 components of targeted buffer : HH

    auto gs_x_yzzzz_xxxxx = pbuffer.data(idx_g_hh + 399);

    auto gs_x_yzzzz_xxxxy = pbuffer.data(idx_g_hh + 400);

    auto gs_x_yzzzz_xxxxz = pbuffer.data(idx_g_hh + 401);

    auto gs_x_yzzzz_xxxyy = pbuffer.data(idx_g_hh + 402);

    auto gs_x_yzzzz_xxxyz = pbuffer.data(idx_g_hh + 403);

    auto gs_x_yzzzz_xxxzz = pbuffer.data(idx_g_hh + 404);

    auto gs_x_yzzzz_xxyyy = pbuffer.data(idx_g_hh + 405);

    auto gs_x_yzzzz_xxyyz = pbuffer.data(idx_g_hh + 406);

    auto gs_x_yzzzz_xxyzz = pbuffer.data(idx_g_hh + 407);

    auto gs_x_yzzzz_xxzzz = pbuffer.data(idx_g_hh + 408);

    auto gs_x_yzzzz_xyyyy = pbuffer.data(idx_g_hh + 409);

    auto gs_x_yzzzz_xyyyz = pbuffer.data(idx_g_hh + 410);

    auto gs_x_yzzzz_xyyzz = pbuffer.data(idx_g_hh + 411);

    auto gs_x_yzzzz_xyzzz = pbuffer.data(idx_g_hh + 412);

    auto gs_x_yzzzz_xzzzz = pbuffer.data(idx_g_hh + 413);

    auto gs_x_yzzzz_yyyyy = pbuffer.data(idx_g_hh + 414);

    auto gs_x_yzzzz_yyyyz = pbuffer.data(idx_g_hh + 415);

    auto gs_x_yzzzz_yyyzz = pbuffer.data(idx_g_hh + 416);

    auto gs_x_yzzzz_yyzzz = pbuffer.data(idx_g_hh + 417);

    auto gs_x_yzzzz_yzzzz = pbuffer.data(idx_g_hh + 418);

    auto gs_x_yzzzz_zzzzz = pbuffer.data(idx_g_hh + 419);

    #pragma omp simd aligned(gc_x, gs_x_yzzzz_xxxxx, gs_x_yzzzz_xxxxy, gs_x_yzzzz_xxxxz, gs_x_yzzzz_xxxyy, gs_x_yzzzz_xxxyz, gs_x_yzzzz_xxxzz, gs_x_yzzzz_xxyyy, gs_x_yzzzz_xxyyz, gs_x_yzzzz_xxyzz, gs_x_yzzzz_xxzzz, gs_x_yzzzz_xyyyy, gs_x_yzzzz_xyyyz, gs_x_yzzzz_xyyzz, gs_x_yzzzz_xyzzz, gs_x_yzzzz_xzzzz, gs_x_yzzzz_yyyyy, gs_x_yzzzz_yyyyz, gs_x_yzzzz_yyyzz, gs_x_yzzzz_yyzzz, gs_x_yzzzz_yzzzz, gs_x_yzzzz_zzzzz, ts_yzzzz_xxxx, ts_yzzzz_xxxxx, ts_yzzzz_xxxxy, ts_yzzzz_xxxxz, ts_yzzzz_xxxy, ts_yzzzz_xxxyy, ts_yzzzz_xxxyz, ts_yzzzz_xxxz, ts_yzzzz_xxxzz, ts_yzzzz_xxyy, ts_yzzzz_xxyyy, ts_yzzzz_xxyyz, ts_yzzzz_xxyz, ts_yzzzz_xxyzz, ts_yzzzz_xxzz, ts_yzzzz_xxzzz, ts_yzzzz_xyyy, ts_yzzzz_xyyyy, ts_yzzzz_xyyyz, ts_yzzzz_xyyz, ts_yzzzz_xyyzz, ts_yzzzz_xyzz, ts_yzzzz_xyzzz, ts_yzzzz_xzzz, ts_yzzzz_xzzzz, ts_yzzzz_yyyy, ts_yzzzz_yyyyy, ts_yzzzz_yyyyz, ts_yzzzz_yyyz, ts_yzzzz_yyyzz, ts_yzzzz_yyzz, ts_yzzzz_yyzzz, ts_yzzzz_yzzz, ts_yzzzz_yzzzz, ts_yzzzz_zzzz, ts_yzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzzz_xxxxx[i] = 10.0 * ts_yzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxxxy[i] = 8.0 * ts_yzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxxxz[i] = 8.0 * ts_yzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxxyy[i] = 6.0 * ts_yzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxxyz[i] = 6.0 * ts_yzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxxzz[i] = 6.0 * ts_yzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxyyy[i] = 4.0 * ts_yzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxyyz[i] = 4.0 * ts_yzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxyzz[i] = 4.0 * ts_yzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xxzzz[i] = 4.0 * ts_yzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyyyy[i] = 2.0 * ts_yzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyyyz[i] = 2.0 * ts_yzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyyzz[i] = 2.0 * ts_yzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xyzzz[i] = 2.0 * ts_yzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xzzzz[i] = 2.0 * ts_yzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyyyy[i] = 2.0 * ts_yzzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyyyz[i] = 2.0 * ts_yzzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyyzz[i] = 2.0 * ts_yzzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yyzzz[i] = 2.0 * ts_yzzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yzzzz[i] = 2.0 * ts_yzzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_zzzzz[i] = 2.0 * ts_yzzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 420-441 components of targeted buffer : HH

    auto gs_x_zzzzz_xxxxx = pbuffer.data(idx_g_hh + 420);

    auto gs_x_zzzzz_xxxxy = pbuffer.data(idx_g_hh + 421);

    auto gs_x_zzzzz_xxxxz = pbuffer.data(idx_g_hh + 422);

    auto gs_x_zzzzz_xxxyy = pbuffer.data(idx_g_hh + 423);

    auto gs_x_zzzzz_xxxyz = pbuffer.data(idx_g_hh + 424);

    auto gs_x_zzzzz_xxxzz = pbuffer.data(idx_g_hh + 425);

    auto gs_x_zzzzz_xxyyy = pbuffer.data(idx_g_hh + 426);

    auto gs_x_zzzzz_xxyyz = pbuffer.data(idx_g_hh + 427);

    auto gs_x_zzzzz_xxyzz = pbuffer.data(idx_g_hh + 428);

    auto gs_x_zzzzz_xxzzz = pbuffer.data(idx_g_hh + 429);

    auto gs_x_zzzzz_xyyyy = pbuffer.data(idx_g_hh + 430);

    auto gs_x_zzzzz_xyyyz = pbuffer.data(idx_g_hh + 431);

    auto gs_x_zzzzz_xyyzz = pbuffer.data(idx_g_hh + 432);

    auto gs_x_zzzzz_xyzzz = pbuffer.data(idx_g_hh + 433);

    auto gs_x_zzzzz_xzzzz = pbuffer.data(idx_g_hh + 434);

    auto gs_x_zzzzz_yyyyy = pbuffer.data(idx_g_hh + 435);

    auto gs_x_zzzzz_yyyyz = pbuffer.data(idx_g_hh + 436);

    auto gs_x_zzzzz_yyyzz = pbuffer.data(idx_g_hh + 437);

    auto gs_x_zzzzz_yyzzz = pbuffer.data(idx_g_hh + 438);

    auto gs_x_zzzzz_yzzzz = pbuffer.data(idx_g_hh + 439);

    auto gs_x_zzzzz_zzzzz = pbuffer.data(idx_g_hh + 440);

    #pragma omp simd aligned(gc_x, gs_x_zzzzz_xxxxx, gs_x_zzzzz_xxxxy, gs_x_zzzzz_xxxxz, gs_x_zzzzz_xxxyy, gs_x_zzzzz_xxxyz, gs_x_zzzzz_xxxzz, gs_x_zzzzz_xxyyy, gs_x_zzzzz_xxyyz, gs_x_zzzzz_xxyzz, gs_x_zzzzz_xxzzz, gs_x_zzzzz_xyyyy, gs_x_zzzzz_xyyyz, gs_x_zzzzz_xyyzz, gs_x_zzzzz_xyzzz, gs_x_zzzzz_xzzzz, gs_x_zzzzz_yyyyy, gs_x_zzzzz_yyyyz, gs_x_zzzzz_yyyzz, gs_x_zzzzz_yyzzz, gs_x_zzzzz_yzzzz, gs_x_zzzzz_zzzzz, ts_zzzzz_xxxx, ts_zzzzz_xxxxx, ts_zzzzz_xxxxy, ts_zzzzz_xxxxz, ts_zzzzz_xxxy, ts_zzzzz_xxxyy, ts_zzzzz_xxxyz, ts_zzzzz_xxxz, ts_zzzzz_xxxzz, ts_zzzzz_xxyy, ts_zzzzz_xxyyy, ts_zzzzz_xxyyz, ts_zzzzz_xxyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyy, ts_zzzzz_xyyyy, ts_zzzzz_xyyyz, ts_zzzzz_xyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyy, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzzz_xxxxx[i] = 10.0 * ts_zzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxxxy[i] = 8.0 * ts_zzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxxxz[i] = 8.0 * ts_zzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxxyy[i] = 6.0 * ts_zzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxxyz[i] = 6.0 * ts_zzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxxzz[i] = 6.0 * ts_zzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxyyy[i] = 4.0 * ts_zzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxyyz[i] = 4.0 * ts_zzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxyzz[i] = 4.0 * ts_zzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xxzzz[i] = 4.0 * ts_zzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyyyy[i] = 2.0 * ts_zzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyyyz[i] = 2.0 * ts_zzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyyzz[i] = 2.0 * ts_zzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xyzzz[i] = 2.0 * ts_zzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xzzzz[i] = 2.0 * ts_zzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyyyy[i] = 2.0 * ts_zzzzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyyyz[i] = 2.0 * ts_zzzzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyyzz[i] = 2.0 * ts_zzzzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yyzzz[i] = 2.0 * ts_zzzzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yzzzz[i] = 2.0 * ts_zzzzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_zzzzz[i] = 2.0 * ts_zzzzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 441-462 components of targeted buffer : HH

    auto gs_y_xxxxx_xxxxx = pbuffer.data(idx_g_hh + 441);

    auto gs_y_xxxxx_xxxxy = pbuffer.data(idx_g_hh + 442);

    auto gs_y_xxxxx_xxxxz = pbuffer.data(idx_g_hh + 443);

    auto gs_y_xxxxx_xxxyy = pbuffer.data(idx_g_hh + 444);

    auto gs_y_xxxxx_xxxyz = pbuffer.data(idx_g_hh + 445);

    auto gs_y_xxxxx_xxxzz = pbuffer.data(idx_g_hh + 446);

    auto gs_y_xxxxx_xxyyy = pbuffer.data(idx_g_hh + 447);

    auto gs_y_xxxxx_xxyyz = pbuffer.data(idx_g_hh + 448);

    auto gs_y_xxxxx_xxyzz = pbuffer.data(idx_g_hh + 449);

    auto gs_y_xxxxx_xxzzz = pbuffer.data(idx_g_hh + 450);

    auto gs_y_xxxxx_xyyyy = pbuffer.data(idx_g_hh + 451);

    auto gs_y_xxxxx_xyyyz = pbuffer.data(idx_g_hh + 452);

    auto gs_y_xxxxx_xyyzz = pbuffer.data(idx_g_hh + 453);

    auto gs_y_xxxxx_xyzzz = pbuffer.data(idx_g_hh + 454);

    auto gs_y_xxxxx_xzzzz = pbuffer.data(idx_g_hh + 455);

    auto gs_y_xxxxx_yyyyy = pbuffer.data(idx_g_hh + 456);

    auto gs_y_xxxxx_yyyyz = pbuffer.data(idx_g_hh + 457);

    auto gs_y_xxxxx_yyyzz = pbuffer.data(idx_g_hh + 458);

    auto gs_y_xxxxx_yyzzz = pbuffer.data(idx_g_hh + 459);

    auto gs_y_xxxxx_yzzzz = pbuffer.data(idx_g_hh + 460);

    auto gs_y_xxxxx_zzzzz = pbuffer.data(idx_g_hh + 461);

    #pragma omp simd aligned(gc_y, gs_y_xxxxx_xxxxx, gs_y_xxxxx_xxxxy, gs_y_xxxxx_xxxxz, gs_y_xxxxx_xxxyy, gs_y_xxxxx_xxxyz, gs_y_xxxxx_xxxzz, gs_y_xxxxx_xxyyy, gs_y_xxxxx_xxyyz, gs_y_xxxxx_xxyzz, gs_y_xxxxx_xxzzz, gs_y_xxxxx_xyyyy, gs_y_xxxxx_xyyyz, gs_y_xxxxx_xyyzz, gs_y_xxxxx_xyzzz, gs_y_xxxxx_xzzzz, gs_y_xxxxx_yyyyy, gs_y_xxxxx_yyyyz, gs_y_xxxxx_yyyzz, gs_y_xxxxx_yyzzz, gs_y_xxxxx_yzzzz, gs_y_xxxxx_zzzzz, ts_xxxxx_xxxx, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxxz, ts_xxxxx_xxxy, ts_xxxxx_xxxyy, ts_xxxxx_xxxyz, ts_xxxxx_xxxz, ts_xxxxx_xxxzz, ts_xxxxx_xxyy, ts_xxxxx_xxyyy, ts_xxxxx_xxyyz, ts_xxxxx_xxyz, ts_xxxxx_xxyzz, ts_xxxxx_xxzz, ts_xxxxx_xxzzz, ts_xxxxx_xyyy, ts_xxxxx_xyyyy, ts_xxxxx_xyyyz, ts_xxxxx_xyyz, ts_xxxxx_xyyzz, ts_xxxxx_xyzz, ts_xxxxx_xyzzz, ts_xxxxx_xzzz, ts_xxxxx_xzzzz, ts_xxxxx_yyyy, ts_xxxxx_yyyyy, ts_xxxxx_yyyyz, ts_xxxxx_yyyz, ts_xxxxx_yyyzz, ts_xxxxx_yyzz, ts_xxxxx_yyzzz, ts_xxxxx_yzzz, ts_xxxxx_yzzzz, ts_xxxxx_zzzz, ts_xxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxx_xxxxx[i] = 2.0 * ts_xxxxx_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxxxy[i] = 2.0 * ts_xxxxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxxxz[i] = 2.0 * ts_xxxxx_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxxyy[i] = 4.0 * ts_xxxxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxxyz[i] = 2.0 * ts_xxxxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxxzz[i] = 2.0 * ts_xxxxx_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxyyy[i] = 6.0 * ts_xxxxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxyyz[i] = 4.0 * ts_xxxxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxyzz[i] = 2.0 * ts_xxxxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xxzzz[i] = 2.0 * ts_xxxxx_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyyyy[i] = 8.0 * ts_xxxxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyyyz[i] = 6.0 * ts_xxxxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyyzz[i] = 4.0 * ts_xxxxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xyzzz[i] = 2.0 * ts_xxxxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xzzzz[i] = 2.0 * ts_xxxxx_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyyyy[i] = 10.0 * ts_xxxxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyyyz[i] = 8.0 * ts_xxxxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyyzz[i] = 6.0 * ts_xxxxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yyzzz[i] = 4.0 * ts_xxxxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yzzzz[i] = 2.0 * ts_xxxxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_zzzzz[i] = 2.0 * ts_xxxxx_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 462-483 components of targeted buffer : HH

    auto gs_y_xxxxy_xxxxx = pbuffer.data(idx_g_hh + 462);

    auto gs_y_xxxxy_xxxxy = pbuffer.data(idx_g_hh + 463);

    auto gs_y_xxxxy_xxxxz = pbuffer.data(idx_g_hh + 464);

    auto gs_y_xxxxy_xxxyy = pbuffer.data(idx_g_hh + 465);

    auto gs_y_xxxxy_xxxyz = pbuffer.data(idx_g_hh + 466);

    auto gs_y_xxxxy_xxxzz = pbuffer.data(idx_g_hh + 467);

    auto gs_y_xxxxy_xxyyy = pbuffer.data(idx_g_hh + 468);

    auto gs_y_xxxxy_xxyyz = pbuffer.data(idx_g_hh + 469);

    auto gs_y_xxxxy_xxyzz = pbuffer.data(idx_g_hh + 470);

    auto gs_y_xxxxy_xxzzz = pbuffer.data(idx_g_hh + 471);

    auto gs_y_xxxxy_xyyyy = pbuffer.data(idx_g_hh + 472);

    auto gs_y_xxxxy_xyyyz = pbuffer.data(idx_g_hh + 473);

    auto gs_y_xxxxy_xyyzz = pbuffer.data(idx_g_hh + 474);

    auto gs_y_xxxxy_xyzzz = pbuffer.data(idx_g_hh + 475);

    auto gs_y_xxxxy_xzzzz = pbuffer.data(idx_g_hh + 476);

    auto gs_y_xxxxy_yyyyy = pbuffer.data(idx_g_hh + 477);

    auto gs_y_xxxxy_yyyyz = pbuffer.data(idx_g_hh + 478);

    auto gs_y_xxxxy_yyyzz = pbuffer.data(idx_g_hh + 479);

    auto gs_y_xxxxy_yyzzz = pbuffer.data(idx_g_hh + 480);

    auto gs_y_xxxxy_yzzzz = pbuffer.data(idx_g_hh + 481);

    auto gs_y_xxxxy_zzzzz = pbuffer.data(idx_g_hh + 482);

    #pragma omp simd aligned(gc_y, gs_y_xxxxy_xxxxx, gs_y_xxxxy_xxxxy, gs_y_xxxxy_xxxxz, gs_y_xxxxy_xxxyy, gs_y_xxxxy_xxxyz, gs_y_xxxxy_xxxzz, gs_y_xxxxy_xxyyy, gs_y_xxxxy_xxyyz, gs_y_xxxxy_xxyzz, gs_y_xxxxy_xxzzz, gs_y_xxxxy_xyyyy, gs_y_xxxxy_xyyyz, gs_y_xxxxy_xyyzz, gs_y_xxxxy_xyzzz, gs_y_xxxxy_xzzzz, gs_y_xxxxy_yyyyy, gs_y_xxxxy_yyyyz, gs_y_xxxxy_yyyzz, gs_y_xxxxy_yyzzz, gs_y_xxxxy_yzzzz, gs_y_xxxxy_zzzzz, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxzz, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyzz, ts_xxxx_xxzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyzz, ts_xxxx_xyzzz, ts_xxxx_xzzzz, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyzz, ts_xxxx_yyzzz, ts_xxxx_yzzzz, ts_xxxx_zzzzz, ts_xxxxy_xxxx, ts_xxxxy_xxxxx, ts_xxxxy_xxxxy, ts_xxxxy_xxxxz, ts_xxxxy_xxxy, ts_xxxxy_xxxyy, ts_xxxxy_xxxyz, ts_xxxxy_xxxz, ts_xxxxy_xxxzz, ts_xxxxy_xxyy, ts_xxxxy_xxyyy, ts_xxxxy_xxyyz, ts_xxxxy_xxyz, ts_xxxxy_xxyzz, ts_xxxxy_xxzz, ts_xxxxy_xxzzz, ts_xxxxy_xyyy, ts_xxxxy_xyyyy, ts_xxxxy_xyyyz, ts_xxxxy_xyyz, ts_xxxxy_xyyzz, ts_xxxxy_xyzz, ts_xxxxy_xyzzz, ts_xxxxy_xzzz, ts_xxxxy_xzzzz, ts_xxxxy_yyyy, ts_xxxxy_yyyyy, ts_xxxxy_yyyyz, ts_xxxxy_yyyz, ts_xxxxy_yyyzz, ts_xxxxy_yyzz, ts_xxxxy_yyzzz, ts_xxxxy_yzzz, ts_xxxxy_yzzzz, ts_xxxxy_zzzz, ts_xxxxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxy_xxxxx[i] = 2.0 * ts_xxxx_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxxxy[i] = 2.0 * ts_xxxx_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxxxz[i] = 2.0 * ts_xxxx_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxxyy[i] = 2.0 * ts_xxxx_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxxyz[i] = 2.0 * ts_xxxx_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxxzz[i] = 2.0 * ts_xxxx_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxyyy[i] = 2.0 * ts_xxxx_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxyyz[i] = 2.0 * ts_xxxx_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxyzz[i] = 2.0 * ts_xxxx_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xxzzz[i] = 2.0 * ts_xxxx_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyyyy[i] = 2.0 * ts_xxxx_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyyyz[i] = 2.0 * ts_xxxx_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyyzz[i] = 2.0 * ts_xxxx_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xyzzz[i] = 2.0 * ts_xxxx_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xzzzz[i] = 2.0 * ts_xxxx_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyyyy[i] = 2.0 * ts_xxxx_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxxxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyyyz[i] = 2.0 * ts_xxxx_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyyzz[i] = 2.0 * ts_xxxx_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yyzzz[i] = 2.0 * ts_xxxx_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yzzzz[i] = 2.0 * ts_xxxx_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_zzzzz[i] = 2.0 * ts_xxxx_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 483-504 components of targeted buffer : HH

    auto gs_y_xxxxz_xxxxx = pbuffer.data(idx_g_hh + 483);

    auto gs_y_xxxxz_xxxxy = pbuffer.data(idx_g_hh + 484);

    auto gs_y_xxxxz_xxxxz = pbuffer.data(idx_g_hh + 485);

    auto gs_y_xxxxz_xxxyy = pbuffer.data(idx_g_hh + 486);

    auto gs_y_xxxxz_xxxyz = pbuffer.data(idx_g_hh + 487);

    auto gs_y_xxxxz_xxxzz = pbuffer.data(idx_g_hh + 488);

    auto gs_y_xxxxz_xxyyy = pbuffer.data(idx_g_hh + 489);

    auto gs_y_xxxxz_xxyyz = pbuffer.data(idx_g_hh + 490);

    auto gs_y_xxxxz_xxyzz = pbuffer.data(idx_g_hh + 491);

    auto gs_y_xxxxz_xxzzz = pbuffer.data(idx_g_hh + 492);

    auto gs_y_xxxxz_xyyyy = pbuffer.data(idx_g_hh + 493);

    auto gs_y_xxxxz_xyyyz = pbuffer.data(idx_g_hh + 494);

    auto gs_y_xxxxz_xyyzz = pbuffer.data(idx_g_hh + 495);

    auto gs_y_xxxxz_xyzzz = pbuffer.data(idx_g_hh + 496);

    auto gs_y_xxxxz_xzzzz = pbuffer.data(idx_g_hh + 497);

    auto gs_y_xxxxz_yyyyy = pbuffer.data(idx_g_hh + 498);

    auto gs_y_xxxxz_yyyyz = pbuffer.data(idx_g_hh + 499);

    auto gs_y_xxxxz_yyyzz = pbuffer.data(idx_g_hh + 500);

    auto gs_y_xxxxz_yyzzz = pbuffer.data(idx_g_hh + 501);

    auto gs_y_xxxxz_yzzzz = pbuffer.data(idx_g_hh + 502);

    auto gs_y_xxxxz_zzzzz = pbuffer.data(idx_g_hh + 503);

    #pragma omp simd aligned(gc_y, gs_y_xxxxz_xxxxx, gs_y_xxxxz_xxxxy, gs_y_xxxxz_xxxxz, gs_y_xxxxz_xxxyy, gs_y_xxxxz_xxxyz, gs_y_xxxxz_xxxzz, gs_y_xxxxz_xxyyy, gs_y_xxxxz_xxyyz, gs_y_xxxxz_xxyzz, gs_y_xxxxz_xxzzz, gs_y_xxxxz_xyyyy, gs_y_xxxxz_xyyyz, gs_y_xxxxz_xyyzz, gs_y_xxxxz_xyzzz, gs_y_xxxxz_xzzzz, gs_y_xxxxz_yyyyy, gs_y_xxxxz_yyyyz, gs_y_xxxxz_yyyzz, gs_y_xxxxz_yyzzz, gs_y_xxxxz_yzzzz, gs_y_xxxxz_zzzzz, ts_xxxxz_xxxx, ts_xxxxz_xxxxx, ts_xxxxz_xxxxy, ts_xxxxz_xxxxz, ts_xxxxz_xxxy, ts_xxxxz_xxxyy, ts_xxxxz_xxxyz, ts_xxxxz_xxxz, ts_xxxxz_xxxzz, ts_xxxxz_xxyy, ts_xxxxz_xxyyy, ts_xxxxz_xxyyz, ts_xxxxz_xxyz, ts_xxxxz_xxyzz, ts_xxxxz_xxzz, ts_xxxxz_xxzzz, ts_xxxxz_xyyy, ts_xxxxz_xyyyy, ts_xxxxz_xyyyz, ts_xxxxz_xyyz, ts_xxxxz_xyyzz, ts_xxxxz_xyzz, ts_xxxxz_xyzzz, ts_xxxxz_xzzz, ts_xxxxz_xzzzz, ts_xxxxz_yyyy, ts_xxxxz_yyyyy, ts_xxxxz_yyyyz, ts_xxxxz_yyyz, ts_xxxxz_yyyzz, ts_xxxxz_yyzz, ts_xxxxz_yyzzz, ts_xxxxz_yzzz, ts_xxxxz_yzzzz, ts_xxxxz_zzzz, ts_xxxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxz_xxxxx[i] = 2.0 * ts_xxxxz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxxxy[i] = 2.0 * ts_xxxxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxxxz[i] = 2.0 * ts_xxxxz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxxyy[i] = 4.0 * ts_xxxxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxxyz[i] = 2.0 * ts_xxxxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxxzz[i] = 2.0 * ts_xxxxz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxyyy[i] = 6.0 * ts_xxxxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxyyz[i] = 4.0 * ts_xxxxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxyzz[i] = 2.0 * ts_xxxxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xxzzz[i] = 2.0 * ts_xxxxz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyyyy[i] = 8.0 * ts_xxxxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyyyz[i] = 6.0 * ts_xxxxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyyzz[i] = 4.0 * ts_xxxxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xyzzz[i] = 2.0 * ts_xxxxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xzzzz[i] = 2.0 * ts_xxxxz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyyyy[i] = 10.0 * ts_xxxxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyyyz[i] = 8.0 * ts_xxxxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyyzz[i] = 6.0 * ts_xxxxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yyzzz[i] = 4.0 * ts_xxxxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yzzzz[i] = 2.0 * ts_xxxxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_zzzzz[i] = 2.0 * ts_xxxxz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 504-525 components of targeted buffer : HH

    auto gs_y_xxxyy_xxxxx = pbuffer.data(idx_g_hh + 504);

    auto gs_y_xxxyy_xxxxy = pbuffer.data(idx_g_hh + 505);

    auto gs_y_xxxyy_xxxxz = pbuffer.data(idx_g_hh + 506);

    auto gs_y_xxxyy_xxxyy = pbuffer.data(idx_g_hh + 507);

    auto gs_y_xxxyy_xxxyz = pbuffer.data(idx_g_hh + 508);

    auto gs_y_xxxyy_xxxzz = pbuffer.data(idx_g_hh + 509);

    auto gs_y_xxxyy_xxyyy = pbuffer.data(idx_g_hh + 510);

    auto gs_y_xxxyy_xxyyz = pbuffer.data(idx_g_hh + 511);

    auto gs_y_xxxyy_xxyzz = pbuffer.data(idx_g_hh + 512);

    auto gs_y_xxxyy_xxzzz = pbuffer.data(idx_g_hh + 513);

    auto gs_y_xxxyy_xyyyy = pbuffer.data(idx_g_hh + 514);

    auto gs_y_xxxyy_xyyyz = pbuffer.data(idx_g_hh + 515);

    auto gs_y_xxxyy_xyyzz = pbuffer.data(idx_g_hh + 516);

    auto gs_y_xxxyy_xyzzz = pbuffer.data(idx_g_hh + 517);

    auto gs_y_xxxyy_xzzzz = pbuffer.data(idx_g_hh + 518);

    auto gs_y_xxxyy_yyyyy = pbuffer.data(idx_g_hh + 519);

    auto gs_y_xxxyy_yyyyz = pbuffer.data(idx_g_hh + 520);

    auto gs_y_xxxyy_yyyzz = pbuffer.data(idx_g_hh + 521);

    auto gs_y_xxxyy_yyzzz = pbuffer.data(idx_g_hh + 522);

    auto gs_y_xxxyy_yzzzz = pbuffer.data(idx_g_hh + 523);

    auto gs_y_xxxyy_zzzzz = pbuffer.data(idx_g_hh + 524);

    #pragma omp simd aligned(gc_y, gs_y_xxxyy_xxxxx, gs_y_xxxyy_xxxxy, gs_y_xxxyy_xxxxz, gs_y_xxxyy_xxxyy, gs_y_xxxyy_xxxyz, gs_y_xxxyy_xxxzz, gs_y_xxxyy_xxyyy, gs_y_xxxyy_xxyyz, gs_y_xxxyy_xxyzz, gs_y_xxxyy_xxzzz, gs_y_xxxyy_xyyyy, gs_y_xxxyy_xyyyz, gs_y_xxxyy_xyyzz, gs_y_xxxyy_xyzzz, gs_y_xxxyy_xzzzz, gs_y_xxxyy_yyyyy, gs_y_xxxyy_yyyyz, gs_y_xxxyy_yyyzz, gs_y_xxxyy_yyzzz, gs_y_xxxyy_yzzzz, gs_y_xxxyy_zzzzz, ts_xxxy_xxxxx, ts_xxxy_xxxxy, ts_xxxy_xxxxz, ts_xxxy_xxxyy, ts_xxxy_xxxyz, ts_xxxy_xxxzz, ts_xxxy_xxyyy, ts_xxxy_xxyyz, ts_xxxy_xxyzz, ts_xxxy_xxzzz, ts_xxxy_xyyyy, ts_xxxy_xyyyz, ts_xxxy_xyyzz, ts_xxxy_xyzzz, ts_xxxy_xzzzz, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyzz, ts_xxxy_yyzzz, ts_xxxy_yzzzz, ts_xxxy_zzzzz, ts_xxxyy_xxxx, ts_xxxyy_xxxxx, ts_xxxyy_xxxxy, ts_xxxyy_xxxxz, ts_xxxyy_xxxy, ts_xxxyy_xxxyy, ts_xxxyy_xxxyz, ts_xxxyy_xxxz, ts_xxxyy_xxxzz, ts_xxxyy_xxyy, ts_xxxyy_xxyyy, ts_xxxyy_xxyyz, ts_xxxyy_xxyz, ts_xxxyy_xxyzz, ts_xxxyy_xxzz, ts_xxxyy_xxzzz, ts_xxxyy_xyyy, ts_xxxyy_xyyyy, ts_xxxyy_xyyyz, ts_xxxyy_xyyz, ts_xxxyy_xyyzz, ts_xxxyy_xyzz, ts_xxxyy_xyzzz, ts_xxxyy_xzzz, ts_xxxyy_xzzzz, ts_xxxyy_yyyy, ts_xxxyy_yyyyy, ts_xxxyy_yyyyz, ts_xxxyy_yyyz, ts_xxxyy_yyyzz, ts_xxxyy_yyzz, ts_xxxyy_yyzzz, ts_xxxyy_yzzz, ts_xxxyy_yzzzz, ts_xxxyy_zzzz, ts_xxxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyy_xxxxx[i] = 4.0 * ts_xxxy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxxxy[i] = 4.0 * ts_xxxy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxxxz[i] = 4.0 * ts_xxxy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxxyy[i] = 4.0 * ts_xxxy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxxyz[i] = 4.0 * ts_xxxy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxxzz[i] = 4.0 * ts_xxxy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxyyy[i] = 4.0 * ts_xxxy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxyyz[i] = 4.0 * ts_xxxy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxyzz[i] = 4.0 * ts_xxxy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xxzzz[i] = 4.0 * ts_xxxy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyyyy[i] = 4.0 * ts_xxxy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyyyz[i] = 4.0 * ts_xxxy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyyzz[i] = 4.0 * ts_xxxy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xyzzz[i] = 4.0 * ts_xxxy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xzzzz[i] = 4.0 * ts_xxxy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyyyy[i] = 4.0 * ts_xxxy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxxyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyyyz[i] = 4.0 * ts_xxxy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyyzz[i] = 4.0 * ts_xxxy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yyzzz[i] = 4.0 * ts_xxxy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yzzzz[i] = 4.0 * ts_xxxy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_zzzzz[i] = 4.0 * ts_xxxy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 525-546 components of targeted buffer : HH

    auto gs_y_xxxyz_xxxxx = pbuffer.data(idx_g_hh + 525);

    auto gs_y_xxxyz_xxxxy = pbuffer.data(idx_g_hh + 526);

    auto gs_y_xxxyz_xxxxz = pbuffer.data(idx_g_hh + 527);

    auto gs_y_xxxyz_xxxyy = pbuffer.data(idx_g_hh + 528);

    auto gs_y_xxxyz_xxxyz = pbuffer.data(idx_g_hh + 529);

    auto gs_y_xxxyz_xxxzz = pbuffer.data(idx_g_hh + 530);

    auto gs_y_xxxyz_xxyyy = pbuffer.data(idx_g_hh + 531);

    auto gs_y_xxxyz_xxyyz = pbuffer.data(idx_g_hh + 532);

    auto gs_y_xxxyz_xxyzz = pbuffer.data(idx_g_hh + 533);

    auto gs_y_xxxyz_xxzzz = pbuffer.data(idx_g_hh + 534);

    auto gs_y_xxxyz_xyyyy = pbuffer.data(idx_g_hh + 535);

    auto gs_y_xxxyz_xyyyz = pbuffer.data(idx_g_hh + 536);

    auto gs_y_xxxyz_xyyzz = pbuffer.data(idx_g_hh + 537);

    auto gs_y_xxxyz_xyzzz = pbuffer.data(idx_g_hh + 538);

    auto gs_y_xxxyz_xzzzz = pbuffer.data(idx_g_hh + 539);

    auto gs_y_xxxyz_yyyyy = pbuffer.data(idx_g_hh + 540);

    auto gs_y_xxxyz_yyyyz = pbuffer.data(idx_g_hh + 541);

    auto gs_y_xxxyz_yyyzz = pbuffer.data(idx_g_hh + 542);

    auto gs_y_xxxyz_yyzzz = pbuffer.data(idx_g_hh + 543);

    auto gs_y_xxxyz_yzzzz = pbuffer.data(idx_g_hh + 544);

    auto gs_y_xxxyz_zzzzz = pbuffer.data(idx_g_hh + 545);

    #pragma omp simd aligned(gc_y, gs_y_xxxyz_xxxxx, gs_y_xxxyz_xxxxy, gs_y_xxxyz_xxxxz, gs_y_xxxyz_xxxyy, gs_y_xxxyz_xxxyz, gs_y_xxxyz_xxxzz, gs_y_xxxyz_xxyyy, gs_y_xxxyz_xxyyz, gs_y_xxxyz_xxyzz, gs_y_xxxyz_xxzzz, gs_y_xxxyz_xyyyy, gs_y_xxxyz_xyyyz, gs_y_xxxyz_xyyzz, gs_y_xxxyz_xyzzz, gs_y_xxxyz_xzzzz, gs_y_xxxyz_yyyyy, gs_y_xxxyz_yyyyz, gs_y_xxxyz_yyyzz, gs_y_xxxyz_yyzzz, gs_y_xxxyz_yzzzz, gs_y_xxxyz_zzzzz, ts_xxxyz_xxxx, ts_xxxyz_xxxxx, ts_xxxyz_xxxxy, ts_xxxyz_xxxxz, ts_xxxyz_xxxy, ts_xxxyz_xxxyy, ts_xxxyz_xxxyz, ts_xxxyz_xxxz, ts_xxxyz_xxxzz, ts_xxxyz_xxyy, ts_xxxyz_xxyyy, ts_xxxyz_xxyyz, ts_xxxyz_xxyz, ts_xxxyz_xxyzz, ts_xxxyz_xxzz, ts_xxxyz_xxzzz, ts_xxxyz_xyyy, ts_xxxyz_xyyyy, ts_xxxyz_xyyyz, ts_xxxyz_xyyz, ts_xxxyz_xyyzz, ts_xxxyz_xyzz, ts_xxxyz_xyzzz, ts_xxxyz_xzzz, ts_xxxyz_xzzzz, ts_xxxyz_yyyy, ts_xxxyz_yyyyy, ts_xxxyz_yyyyz, ts_xxxyz_yyyz, ts_xxxyz_yyyzz, ts_xxxyz_yyzz, ts_xxxyz_yyzzz, ts_xxxyz_yzzz, ts_xxxyz_yzzzz, ts_xxxyz_zzzz, ts_xxxyz_zzzzz, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxxz, ts_xxxz_xxxyy, ts_xxxz_xxxyz, ts_xxxz_xxxzz, ts_xxxz_xxyyy, ts_xxxz_xxyyz, ts_xxxz_xxyzz, ts_xxxz_xxzzz, ts_xxxz_xyyyy, ts_xxxz_xyyyz, ts_xxxz_xyyzz, ts_xxxz_xyzzz, ts_xxxz_xzzzz, ts_xxxz_yyyyy, ts_xxxz_yyyyz, ts_xxxz_yyyzz, ts_xxxz_yyzzz, ts_xxxz_yzzzz, ts_xxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyz_xxxxx[i] = 2.0 * ts_xxxz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxxxy[i] = 2.0 * ts_xxxz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxxxz[i] = 2.0 * ts_xxxz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxxyy[i] = 2.0 * ts_xxxz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxxyz[i] = 2.0 * ts_xxxz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxxzz[i] = 2.0 * ts_xxxz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxyyy[i] = 2.0 * ts_xxxz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxyyz[i] = 2.0 * ts_xxxz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxyzz[i] = 2.0 * ts_xxxz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xxzzz[i] = 2.0 * ts_xxxz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyyyy[i] = 2.0 * ts_xxxz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyyyz[i] = 2.0 * ts_xxxz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyyzz[i] = 2.0 * ts_xxxz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xyzzz[i] = 2.0 * ts_xxxz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xzzzz[i] = 2.0 * ts_xxxz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyyyy[i] = 2.0 * ts_xxxz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxxyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyyyz[i] = 2.0 * ts_xxxz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyyzz[i] = 2.0 * ts_xxxz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yyzzz[i] = 2.0 * ts_xxxz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yzzzz[i] = 2.0 * ts_xxxz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_zzzzz[i] = 2.0 * ts_xxxz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 546-567 components of targeted buffer : HH

    auto gs_y_xxxzz_xxxxx = pbuffer.data(idx_g_hh + 546);

    auto gs_y_xxxzz_xxxxy = pbuffer.data(idx_g_hh + 547);

    auto gs_y_xxxzz_xxxxz = pbuffer.data(idx_g_hh + 548);

    auto gs_y_xxxzz_xxxyy = pbuffer.data(idx_g_hh + 549);

    auto gs_y_xxxzz_xxxyz = pbuffer.data(idx_g_hh + 550);

    auto gs_y_xxxzz_xxxzz = pbuffer.data(idx_g_hh + 551);

    auto gs_y_xxxzz_xxyyy = pbuffer.data(idx_g_hh + 552);

    auto gs_y_xxxzz_xxyyz = pbuffer.data(idx_g_hh + 553);

    auto gs_y_xxxzz_xxyzz = pbuffer.data(idx_g_hh + 554);

    auto gs_y_xxxzz_xxzzz = pbuffer.data(idx_g_hh + 555);

    auto gs_y_xxxzz_xyyyy = pbuffer.data(idx_g_hh + 556);

    auto gs_y_xxxzz_xyyyz = pbuffer.data(idx_g_hh + 557);

    auto gs_y_xxxzz_xyyzz = pbuffer.data(idx_g_hh + 558);

    auto gs_y_xxxzz_xyzzz = pbuffer.data(idx_g_hh + 559);

    auto gs_y_xxxzz_xzzzz = pbuffer.data(idx_g_hh + 560);

    auto gs_y_xxxzz_yyyyy = pbuffer.data(idx_g_hh + 561);

    auto gs_y_xxxzz_yyyyz = pbuffer.data(idx_g_hh + 562);

    auto gs_y_xxxzz_yyyzz = pbuffer.data(idx_g_hh + 563);

    auto gs_y_xxxzz_yyzzz = pbuffer.data(idx_g_hh + 564);

    auto gs_y_xxxzz_yzzzz = pbuffer.data(idx_g_hh + 565);

    auto gs_y_xxxzz_zzzzz = pbuffer.data(idx_g_hh + 566);

    #pragma omp simd aligned(gc_y, gs_y_xxxzz_xxxxx, gs_y_xxxzz_xxxxy, gs_y_xxxzz_xxxxz, gs_y_xxxzz_xxxyy, gs_y_xxxzz_xxxyz, gs_y_xxxzz_xxxzz, gs_y_xxxzz_xxyyy, gs_y_xxxzz_xxyyz, gs_y_xxxzz_xxyzz, gs_y_xxxzz_xxzzz, gs_y_xxxzz_xyyyy, gs_y_xxxzz_xyyyz, gs_y_xxxzz_xyyzz, gs_y_xxxzz_xyzzz, gs_y_xxxzz_xzzzz, gs_y_xxxzz_yyyyy, gs_y_xxxzz_yyyyz, gs_y_xxxzz_yyyzz, gs_y_xxxzz_yyzzz, gs_y_xxxzz_yzzzz, gs_y_xxxzz_zzzzz, ts_xxxzz_xxxx, ts_xxxzz_xxxxx, ts_xxxzz_xxxxy, ts_xxxzz_xxxxz, ts_xxxzz_xxxy, ts_xxxzz_xxxyy, ts_xxxzz_xxxyz, ts_xxxzz_xxxz, ts_xxxzz_xxxzz, ts_xxxzz_xxyy, ts_xxxzz_xxyyy, ts_xxxzz_xxyyz, ts_xxxzz_xxyz, ts_xxxzz_xxyzz, ts_xxxzz_xxzz, ts_xxxzz_xxzzz, ts_xxxzz_xyyy, ts_xxxzz_xyyyy, ts_xxxzz_xyyyz, ts_xxxzz_xyyz, ts_xxxzz_xyyzz, ts_xxxzz_xyzz, ts_xxxzz_xyzzz, ts_xxxzz_xzzz, ts_xxxzz_xzzzz, ts_xxxzz_yyyy, ts_xxxzz_yyyyy, ts_xxxzz_yyyyz, ts_xxxzz_yyyz, ts_xxxzz_yyyzz, ts_xxxzz_yyzz, ts_xxxzz_yyzzz, ts_xxxzz_yzzz, ts_xxxzz_yzzzz, ts_xxxzz_zzzz, ts_xxxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxzz_xxxxx[i] = 2.0 * ts_xxxzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxxxy[i] = 2.0 * ts_xxxzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxxxz[i] = 2.0 * ts_xxxzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxxyy[i] = 4.0 * ts_xxxzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxxyz[i] = 2.0 * ts_xxxzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxxzz[i] = 2.0 * ts_xxxzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxyyy[i] = 6.0 * ts_xxxzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxyyz[i] = 4.0 * ts_xxxzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxyzz[i] = 2.0 * ts_xxxzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xxzzz[i] = 2.0 * ts_xxxzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyyyy[i] = 8.0 * ts_xxxzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyyyz[i] = 6.0 * ts_xxxzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyyzz[i] = 4.0 * ts_xxxzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xyzzz[i] = 2.0 * ts_xxxzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xzzzz[i] = 2.0 * ts_xxxzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyyyy[i] = 10.0 * ts_xxxzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyyyz[i] = 8.0 * ts_xxxzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyyzz[i] = 6.0 * ts_xxxzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yyzzz[i] = 4.0 * ts_xxxzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yzzzz[i] = 2.0 * ts_xxxzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_zzzzz[i] = 2.0 * ts_xxxzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 567-588 components of targeted buffer : HH

    auto gs_y_xxyyy_xxxxx = pbuffer.data(idx_g_hh + 567);

    auto gs_y_xxyyy_xxxxy = pbuffer.data(idx_g_hh + 568);

    auto gs_y_xxyyy_xxxxz = pbuffer.data(idx_g_hh + 569);

    auto gs_y_xxyyy_xxxyy = pbuffer.data(idx_g_hh + 570);

    auto gs_y_xxyyy_xxxyz = pbuffer.data(idx_g_hh + 571);

    auto gs_y_xxyyy_xxxzz = pbuffer.data(idx_g_hh + 572);

    auto gs_y_xxyyy_xxyyy = pbuffer.data(idx_g_hh + 573);

    auto gs_y_xxyyy_xxyyz = pbuffer.data(idx_g_hh + 574);

    auto gs_y_xxyyy_xxyzz = pbuffer.data(idx_g_hh + 575);

    auto gs_y_xxyyy_xxzzz = pbuffer.data(idx_g_hh + 576);

    auto gs_y_xxyyy_xyyyy = pbuffer.data(idx_g_hh + 577);

    auto gs_y_xxyyy_xyyyz = pbuffer.data(idx_g_hh + 578);

    auto gs_y_xxyyy_xyyzz = pbuffer.data(idx_g_hh + 579);

    auto gs_y_xxyyy_xyzzz = pbuffer.data(idx_g_hh + 580);

    auto gs_y_xxyyy_xzzzz = pbuffer.data(idx_g_hh + 581);

    auto gs_y_xxyyy_yyyyy = pbuffer.data(idx_g_hh + 582);

    auto gs_y_xxyyy_yyyyz = pbuffer.data(idx_g_hh + 583);

    auto gs_y_xxyyy_yyyzz = pbuffer.data(idx_g_hh + 584);

    auto gs_y_xxyyy_yyzzz = pbuffer.data(idx_g_hh + 585);

    auto gs_y_xxyyy_yzzzz = pbuffer.data(idx_g_hh + 586);

    auto gs_y_xxyyy_zzzzz = pbuffer.data(idx_g_hh + 587);

    #pragma omp simd aligned(gc_y, gs_y_xxyyy_xxxxx, gs_y_xxyyy_xxxxy, gs_y_xxyyy_xxxxz, gs_y_xxyyy_xxxyy, gs_y_xxyyy_xxxyz, gs_y_xxyyy_xxxzz, gs_y_xxyyy_xxyyy, gs_y_xxyyy_xxyyz, gs_y_xxyyy_xxyzz, gs_y_xxyyy_xxzzz, gs_y_xxyyy_xyyyy, gs_y_xxyyy_xyyyz, gs_y_xxyyy_xyyzz, gs_y_xxyyy_xyzzz, gs_y_xxyyy_xzzzz, gs_y_xxyyy_yyyyy, gs_y_xxyyy_yyyyz, gs_y_xxyyy_yyyzz, gs_y_xxyyy_yyzzz, gs_y_xxyyy_yzzzz, gs_y_xxyyy_zzzzz, ts_xxyy_xxxxx, ts_xxyy_xxxxy, ts_xxyy_xxxxz, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxxzz, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyzz, ts_xxyy_xxzzz, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyzz, ts_xxyy_xyzzz, ts_xxyy_xzzzz, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyzz, ts_xxyy_yyzzz, ts_xxyy_yzzzz, ts_xxyy_zzzzz, ts_xxyyy_xxxx, ts_xxyyy_xxxxx, ts_xxyyy_xxxxy, ts_xxyyy_xxxxz, ts_xxyyy_xxxy, ts_xxyyy_xxxyy, ts_xxyyy_xxxyz, ts_xxyyy_xxxz, ts_xxyyy_xxxzz, ts_xxyyy_xxyy, ts_xxyyy_xxyyy, ts_xxyyy_xxyyz, ts_xxyyy_xxyz, ts_xxyyy_xxyzz, ts_xxyyy_xxzz, ts_xxyyy_xxzzz, ts_xxyyy_xyyy, ts_xxyyy_xyyyy, ts_xxyyy_xyyyz, ts_xxyyy_xyyz, ts_xxyyy_xyyzz, ts_xxyyy_xyzz, ts_xxyyy_xyzzz, ts_xxyyy_xzzz, ts_xxyyy_xzzzz, ts_xxyyy_yyyy, ts_xxyyy_yyyyy, ts_xxyyy_yyyyz, ts_xxyyy_yyyz, ts_xxyyy_yyyzz, ts_xxyyy_yyzz, ts_xxyyy_yyzzz, ts_xxyyy_yzzz, ts_xxyyy_yzzzz, ts_xxyyy_zzzz, ts_xxyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyy_xxxxx[i] = 6.0 * ts_xxyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxxxy[i] = 6.0 * ts_xxyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxxxz[i] = 6.0 * ts_xxyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxxyy[i] = 6.0 * ts_xxyy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxxyz[i] = 6.0 * ts_xxyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxxzz[i] = 6.0 * ts_xxyy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxyyy[i] = 6.0 * ts_xxyy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxyyz[i] = 6.0 * ts_xxyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxyzz[i] = 6.0 * ts_xxyy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xxzzz[i] = 6.0 * ts_xxyy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyyyy[i] = 6.0 * ts_xxyy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyyyz[i] = 6.0 * ts_xxyy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyyzz[i] = 6.0 * ts_xxyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xyzzz[i] = 6.0 * ts_xxyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xzzzz[i] = 6.0 * ts_xxyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyyyy[i] = 6.0 * ts_xxyy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyyyz[i] = 6.0 * ts_xxyy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyyzz[i] = 6.0 * ts_xxyy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yyzzz[i] = 6.0 * ts_xxyy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yzzzz[i] = 6.0 * ts_xxyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_zzzzz[i] = 6.0 * ts_xxyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 588-609 components of targeted buffer : HH

    auto gs_y_xxyyz_xxxxx = pbuffer.data(idx_g_hh + 588);

    auto gs_y_xxyyz_xxxxy = pbuffer.data(idx_g_hh + 589);

    auto gs_y_xxyyz_xxxxz = pbuffer.data(idx_g_hh + 590);

    auto gs_y_xxyyz_xxxyy = pbuffer.data(idx_g_hh + 591);

    auto gs_y_xxyyz_xxxyz = pbuffer.data(idx_g_hh + 592);

    auto gs_y_xxyyz_xxxzz = pbuffer.data(idx_g_hh + 593);

    auto gs_y_xxyyz_xxyyy = pbuffer.data(idx_g_hh + 594);

    auto gs_y_xxyyz_xxyyz = pbuffer.data(idx_g_hh + 595);

    auto gs_y_xxyyz_xxyzz = pbuffer.data(idx_g_hh + 596);

    auto gs_y_xxyyz_xxzzz = pbuffer.data(idx_g_hh + 597);

    auto gs_y_xxyyz_xyyyy = pbuffer.data(idx_g_hh + 598);

    auto gs_y_xxyyz_xyyyz = pbuffer.data(idx_g_hh + 599);

    auto gs_y_xxyyz_xyyzz = pbuffer.data(idx_g_hh + 600);

    auto gs_y_xxyyz_xyzzz = pbuffer.data(idx_g_hh + 601);

    auto gs_y_xxyyz_xzzzz = pbuffer.data(idx_g_hh + 602);

    auto gs_y_xxyyz_yyyyy = pbuffer.data(idx_g_hh + 603);

    auto gs_y_xxyyz_yyyyz = pbuffer.data(idx_g_hh + 604);

    auto gs_y_xxyyz_yyyzz = pbuffer.data(idx_g_hh + 605);

    auto gs_y_xxyyz_yyzzz = pbuffer.data(idx_g_hh + 606);

    auto gs_y_xxyyz_yzzzz = pbuffer.data(idx_g_hh + 607);

    auto gs_y_xxyyz_zzzzz = pbuffer.data(idx_g_hh + 608);

    #pragma omp simd aligned(gc_y, gs_y_xxyyz_xxxxx, gs_y_xxyyz_xxxxy, gs_y_xxyyz_xxxxz, gs_y_xxyyz_xxxyy, gs_y_xxyyz_xxxyz, gs_y_xxyyz_xxxzz, gs_y_xxyyz_xxyyy, gs_y_xxyyz_xxyyz, gs_y_xxyyz_xxyzz, gs_y_xxyyz_xxzzz, gs_y_xxyyz_xyyyy, gs_y_xxyyz_xyyyz, gs_y_xxyyz_xyyzz, gs_y_xxyyz_xyzzz, gs_y_xxyyz_xzzzz, gs_y_xxyyz_yyyyy, gs_y_xxyyz_yyyyz, gs_y_xxyyz_yyyzz, gs_y_xxyyz_yyzzz, gs_y_xxyyz_yzzzz, gs_y_xxyyz_zzzzz, ts_xxyyz_xxxx, ts_xxyyz_xxxxx, ts_xxyyz_xxxxy, ts_xxyyz_xxxxz, ts_xxyyz_xxxy, ts_xxyyz_xxxyy, ts_xxyyz_xxxyz, ts_xxyyz_xxxz, ts_xxyyz_xxxzz, ts_xxyyz_xxyy, ts_xxyyz_xxyyy, ts_xxyyz_xxyyz, ts_xxyyz_xxyz, ts_xxyyz_xxyzz, ts_xxyyz_xxzz, ts_xxyyz_xxzzz, ts_xxyyz_xyyy, ts_xxyyz_xyyyy, ts_xxyyz_xyyyz, ts_xxyyz_xyyz, ts_xxyyz_xyyzz, ts_xxyyz_xyzz, ts_xxyyz_xyzzz, ts_xxyyz_xzzz, ts_xxyyz_xzzzz, ts_xxyyz_yyyy, ts_xxyyz_yyyyy, ts_xxyyz_yyyyz, ts_xxyyz_yyyz, ts_xxyyz_yyyzz, ts_xxyyz_yyzz, ts_xxyyz_yyzzz, ts_xxyyz_yzzz, ts_xxyyz_yzzzz, ts_xxyyz_zzzz, ts_xxyyz_zzzzz, ts_xxyz_xxxxx, ts_xxyz_xxxxy, ts_xxyz_xxxxz, ts_xxyz_xxxyy, ts_xxyz_xxxyz, ts_xxyz_xxxzz, ts_xxyz_xxyyy, ts_xxyz_xxyyz, ts_xxyz_xxyzz, ts_xxyz_xxzzz, ts_xxyz_xyyyy, ts_xxyz_xyyyz, ts_xxyz_xyyzz, ts_xxyz_xyzzz, ts_xxyz_xzzzz, ts_xxyz_yyyyy, ts_xxyz_yyyyz, ts_xxyz_yyyzz, ts_xxyz_yyzzz, ts_xxyz_yzzzz, ts_xxyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyz_xxxxx[i] = 4.0 * ts_xxyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxxxy[i] = 4.0 * ts_xxyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxxxz[i] = 4.0 * ts_xxyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxxyy[i] = 4.0 * ts_xxyz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxxyz[i] = 4.0 * ts_xxyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxxzz[i] = 4.0 * ts_xxyz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxyyy[i] = 4.0 * ts_xxyz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxyyz[i] = 4.0 * ts_xxyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxyzz[i] = 4.0 * ts_xxyz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xxzzz[i] = 4.0 * ts_xxyz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyyyy[i] = 4.0 * ts_xxyz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyyyz[i] = 4.0 * ts_xxyz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyyzz[i] = 4.0 * ts_xxyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xyzzz[i] = 4.0 * ts_xxyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xzzzz[i] = 4.0 * ts_xxyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyyyy[i] = 4.0 * ts_xxyz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyyyz[i] = 4.0 * ts_xxyz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyyzz[i] = 4.0 * ts_xxyz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yyzzz[i] = 4.0 * ts_xxyz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yzzzz[i] = 4.0 * ts_xxyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_zzzzz[i] = 4.0 * ts_xxyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 609-630 components of targeted buffer : HH

    auto gs_y_xxyzz_xxxxx = pbuffer.data(idx_g_hh + 609);

    auto gs_y_xxyzz_xxxxy = pbuffer.data(idx_g_hh + 610);

    auto gs_y_xxyzz_xxxxz = pbuffer.data(idx_g_hh + 611);

    auto gs_y_xxyzz_xxxyy = pbuffer.data(idx_g_hh + 612);

    auto gs_y_xxyzz_xxxyz = pbuffer.data(idx_g_hh + 613);

    auto gs_y_xxyzz_xxxzz = pbuffer.data(idx_g_hh + 614);

    auto gs_y_xxyzz_xxyyy = pbuffer.data(idx_g_hh + 615);

    auto gs_y_xxyzz_xxyyz = pbuffer.data(idx_g_hh + 616);

    auto gs_y_xxyzz_xxyzz = pbuffer.data(idx_g_hh + 617);

    auto gs_y_xxyzz_xxzzz = pbuffer.data(idx_g_hh + 618);

    auto gs_y_xxyzz_xyyyy = pbuffer.data(idx_g_hh + 619);

    auto gs_y_xxyzz_xyyyz = pbuffer.data(idx_g_hh + 620);

    auto gs_y_xxyzz_xyyzz = pbuffer.data(idx_g_hh + 621);

    auto gs_y_xxyzz_xyzzz = pbuffer.data(idx_g_hh + 622);

    auto gs_y_xxyzz_xzzzz = pbuffer.data(idx_g_hh + 623);

    auto gs_y_xxyzz_yyyyy = pbuffer.data(idx_g_hh + 624);

    auto gs_y_xxyzz_yyyyz = pbuffer.data(idx_g_hh + 625);

    auto gs_y_xxyzz_yyyzz = pbuffer.data(idx_g_hh + 626);

    auto gs_y_xxyzz_yyzzz = pbuffer.data(idx_g_hh + 627);

    auto gs_y_xxyzz_yzzzz = pbuffer.data(idx_g_hh + 628);

    auto gs_y_xxyzz_zzzzz = pbuffer.data(idx_g_hh + 629);

    #pragma omp simd aligned(gc_y, gs_y_xxyzz_xxxxx, gs_y_xxyzz_xxxxy, gs_y_xxyzz_xxxxz, gs_y_xxyzz_xxxyy, gs_y_xxyzz_xxxyz, gs_y_xxyzz_xxxzz, gs_y_xxyzz_xxyyy, gs_y_xxyzz_xxyyz, gs_y_xxyzz_xxyzz, gs_y_xxyzz_xxzzz, gs_y_xxyzz_xyyyy, gs_y_xxyzz_xyyyz, gs_y_xxyzz_xyyzz, gs_y_xxyzz_xyzzz, gs_y_xxyzz_xzzzz, gs_y_xxyzz_yyyyy, gs_y_xxyzz_yyyyz, gs_y_xxyzz_yyyzz, gs_y_xxyzz_yyzzz, gs_y_xxyzz_yzzzz, gs_y_xxyzz_zzzzz, ts_xxyzz_xxxx, ts_xxyzz_xxxxx, ts_xxyzz_xxxxy, ts_xxyzz_xxxxz, ts_xxyzz_xxxy, ts_xxyzz_xxxyy, ts_xxyzz_xxxyz, ts_xxyzz_xxxz, ts_xxyzz_xxxzz, ts_xxyzz_xxyy, ts_xxyzz_xxyyy, ts_xxyzz_xxyyz, ts_xxyzz_xxyz, ts_xxyzz_xxyzz, ts_xxyzz_xxzz, ts_xxyzz_xxzzz, ts_xxyzz_xyyy, ts_xxyzz_xyyyy, ts_xxyzz_xyyyz, ts_xxyzz_xyyz, ts_xxyzz_xyyzz, ts_xxyzz_xyzz, ts_xxyzz_xyzzz, ts_xxyzz_xzzz, ts_xxyzz_xzzzz, ts_xxyzz_yyyy, ts_xxyzz_yyyyy, ts_xxyzz_yyyyz, ts_xxyzz_yyyz, ts_xxyzz_yyyzz, ts_xxyzz_yyzz, ts_xxyzz_yyzzz, ts_xxyzz_yzzz, ts_xxyzz_yzzzz, ts_xxyzz_zzzz, ts_xxyzz_zzzzz, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxxz, ts_xxzz_xxxyy, ts_xxzz_xxxyz, ts_xxzz_xxxzz, ts_xxzz_xxyyy, ts_xxzz_xxyyz, ts_xxzz_xxyzz, ts_xxzz_xxzzz, ts_xxzz_xyyyy, ts_xxzz_xyyyz, ts_xxzz_xyyzz, ts_xxzz_xyzzz, ts_xxzz_xzzzz, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyzz, ts_xxzz_yyzzz, ts_xxzz_yzzzz, ts_xxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyzz_xxxxx[i] = 2.0 * ts_xxzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxxxy[i] = 2.0 * ts_xxzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxxxz[i] = 2.0 * ts_xxzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxxyy[i] = 2.0 * ts_xxzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxxyz[i] = 2.0 * ts_xxzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxxzz[i] = 2.0 * ts_xxzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxyyy[i] = 2.0 * ts_xxzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxyyz[i] = 2.0 * ts_xxzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxyzz[i] = 2.0 * ts_xxzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xxzzz[i] = 2.0 * ts_xxzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyyyy[i] = 2.0 * ts_xxzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyyyz[i] = 2.0 * ts_xxzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyyzz[i] = 2.0 * ts_xxzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xyzzz[i] = 2.0 * ts_xxzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xzzzz[i] = 2.0 * ts_xxzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyyyy[i] = 2.0 * ts_xxzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyyyz[i] = 2.0 * ts_xxzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyyzz[i] = 2.0 * ts_xxzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yyzzz[i] = 2.0 * ts_xxzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yzzzz[i] = 2.0 * ts_xxzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_zzzzz[i] = 2.0 * ts_xxzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 630-651 components of targeted buffer : HH

    auto gs_y_xxzzz_xxxxx = pbuffer.data(idx_g_hh + 630);

    auto gs_y_xxzzz_xxxxy = pbuffer.data(idx_g_hh + 631);

    auto gs_y_xxzzz_xxxxz = pbuffer.data(idx_g_hh + 632);

    auto gs_y_xxzzz_xxxyy = pbuffer.data(idx_g_hh + 633);

    auto gs_y_xxzzz_xxxyz = pbuffer.data(idx_g_hh + 634);

    auto gs_y_xxzzz_xxxzz = pbuffer.data(idx_g_hh + 635);

    auto gs_y_xxzzz_xxyyy = pbuffer.data(idx_g_hh + 636);

    auto gs_y_xxzzz_xxyyz = pbuffer.data(idx_g_hh + 637);

    auto gs_y_xxzzz_xxyzz = pbuffer.data(idx_g_hh + 638);

    auto gs_y_xxzzz_xxzzz = pbuffer.data(idx_g_hh + 639);

    auto gs_y_xxzzz_xyyyy = pbuffer.data(idx_g_hh + 640);

    auto gs_y_xxzzz_xyyyz = pbuffer.data(idx_g_hh + 641);

    auto gs_y_xxzzz_xyyzz = pbuffer.data(idx_g_hh + 642);

    auto gs_y_xxzzz_xyzzz = pbuffer.data(idx_g_hh + 643);

    auto gs_y_xxzzz_xzzzz = pbuffer.data(idx_g_hh + 644);

    auto gs_y_xxzzz_yyyyy = pbuffer.data(idx_g_hh + 645);

    auto gs_y_xxzzz_yyyyz = pbuffer.data(idx_g_hh + 646);

    auto gs_y_xxzzz_yyyzz = pbuffer.data(idx_g_hh + 647);

    auto gs_y_xxzzz_yyzzz = pbuffer.data(idx_g_hh + 648);

    auto gs_y_xxzzz_yzzzz = pbuffer.data(idx_g_hh + 649);

    auto gs_y_xxzzz_zzzzz = pbuffer.data(idx_g_hh + 650);

    #pragma omp simd aligned(gc_y, gs_y_xxzzz_xxxxx, gs_y_xxzzz_xxxxy, gs_y_xxzzz_xxxxz, gs_y_xxzzz_xxxyy, gs_y_xxzzz_xxxyz, gs_y_xxzzz_xxxzz, gs_y_xxzzz_xxyyy, gs_y_xxzzz_xxyyz, gs_y_xxzzz_xxyzz, gs_y_xxzzz_xxzzz, gs_y_xxzzz_xyyyy, gs_y_xxzzz_xyyyz, gs_y_xxzzz_xyyzz, gs_y_xxzzz_xyzzz, gs_y_xxzzz_xzzzz, gs_y_xxzzz_yyyyy, gs_y_xxzzz_yyyyz, gs_y_xxzzz_yyyzz, gs_y_xxzzz_yyzzz, gs_y_xxzzz_yzzzz, gs_y_xxzzz_zzzzz, ts_xxzzz_xxxx, ts_xxzzz_xxxxx, ts_xxzzz_xxxxy, ts_xxzzz_xxxxz, ts_xxzzz_xxxy, ts_xxzzz_xxxyy, ts_xxzzz_xxxyz, ts_xxzzz_xxxz, ts_xxzzz_xxxzz, ts_xxzzz_xxyy, ts_xxzzz_xxyyy, ts_xxzzz_xxyyz, ts_xxzzz_xxyz, ts_xxzzz_xxyzz, ts_xxzzz_xxzz, ts_xxzzz_xxzzz, ts_xxzzz_xyyy, ts_xxzzz_xyyyy, ts_xxzzz_xyyyz, ts_xxzzz_xyyz, ts_xxzzz_xyyzz, ts_xxzzz_xyzz, ts_xxzzz_xyzzz, ts_xxzzz_xzzz, ts_xxzzz_xzzzz, ts_xxzzz_yyyy, ts_xxzzz_yyyyy, ts_xxzzz_yyyyz, ts_xxzzz_yyyz, ts_xxzzz_yyyzz, ts_xxzzz_yyzz, ts_xxzzz_yyzzz, ts_xxzzz_yzzz, ts_xxzzz_yzzzz, ts_xxzzz_zzzz, ts_xxzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzzz_xxxxx[i] = 2.0 * ts_xxzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxxxy[i] = 2.0 * ts_xxzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxxxz[i] = 2.0 * ts_xxzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxxyy[i] = 4.0 * ts_xxzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxxyz[i] = 2.0 * ts_xxzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxxzz[i] = 2.0 * ts_xxzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxyyy[i] = 6.0 * ts_xxzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxyyz[i] = 4.0 * ts_xxzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxyzz[i] = 2.0 * ts_xxzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xxzzz[i] = 2.0 * ts_xxzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyyyy[i] = 8.0 * ts_xxzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyyyz[i] = 6.0 * ts_xxzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyyzz[i] = 4.0 * ts_xxzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xyzzz[i] = 2.0 * ts_xxzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xzzzz[i] = 2.0 * ts_xxzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyyyy[i] = 10.0 * ts_xxzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyyyz[i] = 8.0 * ts_xxzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyyzz[i] = 6.0 * ts_xxzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yyzzz[i] = 4.0 * ts_xxzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yzzzz[i] = 2.0 * ts_xxzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_zzzzz[i] = 2.0 * ts_xxzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 651-672 components of targeted buffer : HH

    auto gs_y_xyyyy_xxxxx = pbuffer.data(idx_g_hh + 651);

    auto gs_y_xyyyy_xxxxy = pbuffer.data(idx_g_hh + 652);

    auto gs_y_xyyyy_xxxxz = pbuffer.data(idx_g_hh + 653);

    auto gs_y_xyyyy_xxxyy = pbuffer.data(idx_g_hh + 654);

    auto gs_y_xyyyy_xxxyz = pbuffer.data(idx_g_hh + 655);

    auto gs_y_xyyyy_xxxzz = pbuffer.data(idx_g_hh + 656);

    auto gs_y_xyyyy_xxyyy = pbuffer.data(idx_g_hh + 657);

    auto gs_y_xyyyy_xxyyz = pbuffer.data(idx_g_hh + 658);

    auto gs_y_xyyyy_xxyzz = pbuffer.data(idx_g_hh + 659);

    auto gs_y_xyyyy_xxzzz = pbuffer.data(idx_g_hh + 660);

    auto gs_y_xyyyy_xyyyy = pbuffer.data(idx_g_hh + 661);

    auto gs_y_xyyyy_xyyyz = pbuffer.data(idx_g_hh + 662);

    auto gs_y_xyyyy_xyyzz = pbuffer.data(idx_g_hh + 663);

    auto gs_y_xyyyy_xyzzz = pbuffer.data(idx_g_hh + 664);

    auto gs_y_xyyyy_xzzzz = pbuffer.data(idx_g_hh + 665);

    auto gs_y_xyyyy_yyyyy = pbuffer.data(idx_g_hh + 666);

    auto gs_y_xyyyy_yyyyz = pbuffer.data(idx_g_hh + 667);

    auto gs_y_xyyyy_yyyzz = pbuffer.data(idx_g_hh + 668);

    auto gs_y_xyyyy_yyzzz = pbuffer.data(idx_g_hh + 669);

    auto gs_y_xyyyy_yzzzz = pbuffer.data(idx_g_hh + 670);

    auto gs_y_xyyyy_zzzzz = pbuffer.data(idx_g_hh + 671);

    #pragma omp simd aligned(gc_y, gs_y_xyyyy_xxxxx, gs_y_xyyyy_xxxxy, gs_y_xyyyy_xxxxz, gs_y_xyyyy_xxxyy, gs_y_xyyyy_xxxyz, gs_y_xyyyy_xxxzz, gs_y_xyyyy_xxyyy, gs_y_xyyyy_xxyyz, gs_y_xyyyy_xxyzz, gs_y_xyyyy_xxzzz, gs_y_xyyyy_xyyyy, gs_y_xyyyy_xyyyz, gs_y_xyyyy_xyyzz, gs_y_xyyyy_xyzzz, gs_y_xyyyy_xzzzz, gs_y_xyyyy_yyyyy, gs_y_xyyyy_yyyyz, gs_y_xyyyy_yyyzz, gs_y_xyyyy_yyzzz, gs_y_xyyyy_yzzzz, gs_y_xyyyy_zzzzz, ts_xyyy_xxxxx, ts_xyyy_xxxxy, ts_xyyy_xxxxz, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxxzz, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyzz, ts_xyyy_xxzzz, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyzz, ts_xyyy_xyzzz, ts_xyyy_xzzzz, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyzz, ts_xyyy_yyzzz, ts_xyyy_yzzzz, ts_xyyy_zzzzz, ts_xyyyy_xxxx, ts_xyyyy_xxxxx, ts_xyyyy_xxxxy, ts_xyyyy_xxxxz, ts_xyyyy_xxxy, ts_xyyyy_xxxyy, ts_xyyyy_xxxyz, ts_xyyyy_xxxz, ts_xyyyy_xxxzz, ts_xyyyy_xxyy, ts_xyyyy_xxyyy, ts_xyyyy_xxyyz, ts_xyyyy_xxyz, ts_xyyyy_xxyzz, ts_xyyyy_xxzz, ts_xyyyy_xxzzz, ts_xyyyy_xyyy, ts_xyyyy_xyyyy, ts_xyyyy_xyyyz, ts_xyyyy_xyyz, ts_xyyyy_xyyzz, ts_xyyyy_xyzz, ts_xyyyy_xyzzz, ts_xyyyy_xzzz, ts_xyyyy_xzzzz, ts_xyyyy_yyyy, ts_xyyyy_yyyyy, ts_xyyyy_yyyyz, ts_xyyyy_yyyz, ts_xyyyy_yyyzz, ts_xyyyy_yyzz, ts_xyyyy_yyzzz, ts_xyyyy_yzzz, ts_xyyyy_yzzzz, ts_xyyyy_zzzz, ts_xyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyy_xxxxx[i] = 8.0 * ts_xyyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxxxy[i] = 8.0 * ts_xyyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxxxz[i] = 8.0 * ts_xyyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxxyy[i] = 8.0 * ts_xyyy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxxyz[i] = 8.0 * ts_xyyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxxzz[i] = 8.0 * ts_xyyy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxyyy[i] = 8.0 * ts_xyyy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxyyz[i] = 8.0 * ts_xyyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxyzz[i] = 8.0 * ts_xyyy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xxzzz[i] = 8.0 * ts_xyyy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyyyy[i] = 8.0 * ts_xyyy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyyyz[i] = 8.0 * ts_xyyy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyyzz[i] = 8.0 * ts_xyyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xyzzz[i] = 8.0 * ts_xyyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xzzzz[i] = 8.0 * ts_xyyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyyyy[i] = 8.0 * ts_xyyy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyyyz[i] = 8.0 * ts_xyyy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyyzz[i] = 8.0 * ts_xyyy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yyzzz[i] = 8.0 * ts_xyyy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yzzzz[i] = 8.0 * ts_xyyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_zzzzz[i] = 8.0 * ts_xyyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 672-693 components of targeted buffer : HH

    auto gs_y_xyyyz_xxxxx = pbuffer.data(idx_g_hh + 672);

    auto gs_y_xyyyz_xxxxy = pbuffer.data(idx_g_hh + 673);

    auto gs_y_xyyyz_xxxxz = pbuffer.data(idx_g_hh + 674);

    auto gs_y_xyyyz_xxxyy = pbuffer.data(idx_g_hh + 675);

    auto gs_y_xyyyz_xxxyz = pbuffer.data(idx_g_hh + 676);

    auto gs_y_xyyyz_xxxzz = pbuffer.data(idx_g_hh + 677);

    auto gs_y_xyyyz_xxyyy = pbuffer.data(idx_g_hh + 678);

    auto gs_y_xyyyz_xxyyz = pbuffer.data(idx_g_hh + 679);

    auto gs_y_xyyyz_xxyzz = pbuffer.data(idx_g_hh + 680);

    auto gs_y_xyyyz_xxzzz = pbuffer.data(idx_g_hh + 681);

    auto gs_y_xyyyz_xyyyy = pbuffer.data(idx_g_hh + 682);

    auto gs_y_xyyyz_xyyyz = pbuffer.data(idx_g_hh + 683);

    auto gs_y_xyyyz_xyyzz = pbuffer.data(idx_g_hh + 684);

    auto gs_y_xyyyz_xyzzz = pbuffer.data(idx_g_hh + 685);

    auto gs_y_xyyyz_xzzzz = pbuffer.data(idx_g_hh + 686);

    auto gs_y_xyyyz_yyyyy = pbuffer.data(idx_g_hh + 687);

    auto gs_y_xyyyz_yyyyz = pbuffer.data(idx_g_hh + 688);

    auto gs_y_xyyyz_yyyzz = pbuffer.data(idx_g_hh + 689);

    auto gs_y_xyyyz_yyzzz = pbuffer.data(idx_g_hh + 690);

    auto gs_y_xyyyz_yzzzz = pbuffer.data(idx_g_hh + 691);

    auto gs_y_xyyyz_zzzzz = pbuffer.data(idx_g_hh + 692);

    #pragma omp simd aligned(gc_y, gs_y_xyyyz_xxxxx, gs_y_xyyyz_xxxxy, gs_y_xyyyz_xxxxz, gs_y_xyyyz_xxxyy, gs_y_xyyyz_xxxyz, gs_y_xyyyz_xxxzz, gs_y_xyyyz_xxyyy, gs_y_xyyyz_xxyyz, gs_y_xyyyz_xxyzz, gs_y_xyyyz_xxzzz, gs_y_xyyyz_xyyyy, gs_y_xyyyz_xyyyz, gs_y_xyyyz_xyyzz, gs_y_xyyyz_xyzzz, gs_y_xyyyz_xzzzz, gs_y_xyyyz_yyyyy, gs_y_xyyyz_yyyyz, gs_y_xyyyz_yyyzz, gs_y_xyyyz_yyzzz, gs_y_xyyyz_yzzzz, gs_y_xyyyz_zzzzz, ts_xyyyz_xxxx, ts_xyyyz_xxxxx, ts_xyyyz_xxxxy, ts_xyyyz_xxxxz, ts_xyyyz_xxxy, ts_xyyyz_xxxyy, ts_xyyyz_xxxyz, ts_xyyyz_xxxz, ts_xyyyz_xxxzz, ts_xyyyz_xxyy, ts_xyyyz_xxyyy, ts_xyyyz_xxyyz, ts_xyyyz_xxyz, ts_xyyyz_xxyzz, ts_xyyyz_xxzz, ts_xyyyz_xxzzz, ts_xyyyz_xyyy, ts_xyyyz_xyyyy, ts_xyyyz_xyyyz, ts_xyyyz_xyyz, ts_xyyyz_xyyzz, ts_xyyyz_xyzz, ts_xyyyz_xyzzz, ts_xyyyz_xzzz, ts_xyyyz_xzzzz, ts_xyyyz_yyyy, ts_xyyyz_yyyyy, ts_xyyyz_yyyyz, ts_xyyyz_yyyz, ts_xyyyz_yyyzz, ts_xyyyz_yyzz, ts_xyyyz_yyzzz, ts_xyyyz_yzzz, ts_xyyyz_yzzzz, ts_xyyyz_zzzz, ts_xyyyz_zzzzz, ts_xyyz_xxxxx, ts_xyyz_xxxxy, ts_xyyz_xxxxz, ts_xyyz_xxxyy, ts_xyyz_xxxyz, ts_xyyz_xxxzz, ts_xyyz_xxyyy, ts_xyyz_xxyyz, ts_xyyz_xxyzz, ts_xyyz_xxzzz, ts_xyyz_xyyyy, ts_xyyz_xyyyz, ts_xyyz_xyyzz, ts_xyyz_xyzzz, ts_xyyz_xzzzz, ts_xyyz_yyyyy, ts_xyyz_yyyyz, ts_xyyz_yyyzz, ts_xyyz_yyzzz, ts_xyyz_yzzzz, ts_xyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyz_xxxxx[i] = 6.0 * ts_xyyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxxxy[i] = 6.0 * ts_xyyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxxxz[i] = 6.0 * ts_xyyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxxyy[i] = 6.0 * ts_xyyz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxxyz[i] = 6.0 * ts_xyyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxxzz[i] = 6.0 * ts_xyyz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxyyy[i] = 6.0 * ts_xyyz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxyyz[i] = 6.0 * ts_xyyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxyzz[i] = 6.0 * ts_xyyz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xxzzz[i] = 6.0 * ts_xyyz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyyyy[i] = 6.0 * ts_xyyz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyyyz[i] = 6.0 * ts_xyyz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyyzz[i] = 6.0 * ts_xyyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xyzzz[i] = 6.0 * ts_xyyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xzzzz[i] = 6.0 * ts_xyyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyyyy[i] = 6.0 * ts_xyyz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyyyz[i] = 6.0 * ts_xyyz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyyzz[i] = 6.0 * ts_xyyz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yyzzz[i] = 6.0 * ts_xyyz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yzzzz[i] = 6.0 * ts_xyyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_zzzzz[i] = 6.0 * ts_xyyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 693-714 components of targeted buffer : HH

    auto gs_y_xyyzz_xxxxx = pbuffer.data(idx_g_hh + 693);

    auto gs_y_xyyzz_xxxxy = pbuffer.data(idx_g_hh + 694);

    auto gs_y_xyyzz_xxxxz = pbuffer.data(idx_g_hh + 695);

    auto gs_y_xyyzz_xxxyy = pbuffer.data(idx_g_hh + 696);

    auto gs_y_xyyzz_xxxyz = pbuffer.data(idx_g_hh + 697);

    auto gs_y_xyyzz_xxxzz = pbuffer.data(idx_g_hh + 698);

    auto gs_y_xyyzz_xxyyy = pbuffer.data(idx_g_hh + 699);

    auto gs_y_xyyzz_xxyyz = pbuffer.data(idx_g_hh + 700);

    auto gs_y_xyyzz_xxyzz = pbuffer.data(idx_g_hh + 701);

    auto gs_y_xyyzz_xxzzz = pbuffer.data(idx_g_hh + 702);

    auto gs_y_xyyzz_xyyyy = pbuffer.data(idx_g_hh + 703);

    auto gs_y_xyyzz_xyyyz = pbuffer.data(idx_g_hh + 704);

    auto gs_y_xyyzz_xyyzz = pbuffer.data(idx_g_hh + 705);

    auto gs_y_xyyzz_xyzzz = pbuffer.data(idx_g_hh + 706);

    auto gs_y_xyyzz_xzzzz = pbuffer.data(idx_g_hh + 707);

    auto gs_y_xyyzz_yyyyy = pbuffer.data(idx_g_hh + 708);

    auto gs_y_xyyzz_yyyyz = pbuffer.data(idx_g_hh + 709);

    auto gs_y_xyyzz_yyyzz = pbuffer.data(idx_g_hh + 710);

    auto gs_y_xyyzz_yyzzz = pbuffer.data(idx_g_hh + 711);

    auto gs_y_xyyzz_yzzzz = pbuffer.data(idx_g_hh + 712);

    auto gs_y_xyyzz_zzzzz = pbuffer.data(idx_g_hh + 713);

    #pragma omp simd aligned(gc_y, gs_y_xyyzz_xxxxx, gs_y_xyyzz_xxxxy, gs_y_xyyzz_xxxxz, gs_y_xyyzz_xxxyy, gs_y_xyyzz_xxxyz, gs_y_xyyzz_xxxzz, gs_y_xyyzz_xxyyy, gs_y_xyyzz_xxyyz, gs_y_xyyzz_xxyzz, gs_y_xyyzz_xxzzz, gs_y_xyyzz_xyyyy, gs_y_xyyzz_xyyyz, gs_y_xyyzz_xyyzz, gs_y_xyyzz_xyzzz, gs_y_xyyzz_xzzzz, gs_y_xyyzz_yyyyy, gs_y_xyyzz_yyyyz, gs_y_xyyzz_yyyzz, gs_y_xyyzz_yyzzz, gs_y_xyyzz_yzzzz, gs_y_xyyzz_zzzzz, ts_xyyzz_xxxx, ts_xyyzz_xxxxx, ts_xyyzz_xxxxy, ts_xyyzz_xxxxz, ts_xyyzz_xxxy, ts_xyyzz_xxxyy, ts_xyyzz_xxxyz, ts_xyyzz_xxxz, ts_xyyzz_xxxzz, ts_xyyzz_xxyy, ts_xyyzz_xxyyy, ts_xyyzz_xxyyz, ts_xyyzz_xxyz, ts_xyyzz_xxyzz, ts_xyyzz_xxzz, ts_xyyzz_xxzzz, ts_xyyzz_xyyy, ts_xyyzz_xyyyy, ts_xyyzz_xyyyz, ts_xyyzz_xyyz, ts_xyyzz_xyyzz, ts_xyyzz_xyzz, ts_xyyzz_xyzzz, ts_xyyzz_xzzz, ts_xyyzz_xzzzz, ts_xyyzz_yyyy, ts_xyyzz_yyyyy, ts_xyyzz_yyyyz, ts_xyyzz_yyyz, ts_xyyzz_yyyzz, ts_xyyzz_yyzz, ts_xyyzz_yyzzz, ts_xyyzz_yzzz, ts_xyyzz_yzzzz, ts_xyyzz_zzzz, ts_xyyzz_zzzzz, ts_xyzz_xxxxx, ts_xyzz_xxxxy, ts_xyzz_xxxxz, ts_xyzz_xxxyy, ts_xyzz_xxxyz, ts_xyzz_xxxzz, ts_xyzz_xxyyy, ts_xyzz_xxyyz, ts_xyzz_xxyzz, ts_xyzz_xxzzz, ts_xyzz_xyyyy, ts_xyzz_xyyyz, ts_xyzz_xyyzz, ts_xyzz_xyzzz, ts_xyzz_xzzzz, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyzz, ts_xyzz_yyzzz, ts_xyzz_yzzzz, ts_xyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyzz_xxxxx[i] = 4.0 * ts_xyzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxxxy[i] = 4.0 * ts_xyzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxxxz[i] = 4.0 * ts_xyzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxxyy[i] = 4.0 * ts_xyzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxxyz[i] = 4.0 * ts_xyzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxxzz[i] = 4.0 * ts_xyzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxyyy[i] = 4.0 * ts_xyzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxyyz[i] = 4.0 * ts_xyzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxyzz[i] = 4.0 * ts_xyzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xxzzz[i] = 4.0 * ts_xyzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyyyy[i] = 4.0 * ts_xyzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyyyz[i] = 4.0 * ts_xyzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyyzz[i] = 4.0 * ts_xyzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xyzzz[i] = 4.0 * ts_xyzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xzzzz[i] = 4.0 * ts_xyzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyyyy[i] = 4.0 * ts_xyzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyyyz[i] = 4.0 * ts_xyzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyyzz[i] = 4.0 * ts_xyzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yyzzz[i] = 4.0 * ts_xyzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yzzzz[i] = 4.0 * ts_xyzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_zzzzz[i] = 4.0 * ts_xyzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 714-735 components of targeted buffer : HH

    auto gs_y_xyzzz_xxxxx = pbuffer.data(idx_g_hh + 714);

    auto gs_y_xyzzz_xxxxy = pbuffer.data(idx_g_hh + 715);

    auto gs_y_xyzzz_xxxxz = pbuffer.data(idx_g_hh + 716);

    auto gs_y_xyzzz_xxxyy = pbuffer.data(idx_g_hh + 717);

    auto gs_y_xyzzz_xxxyz = pbuffer.data(idx_g_hh + 718);

    auto gs_y_xyzzz_xxxzz = pbuffer.data(idx_g_hh + 719);

    auto gs_y_xyzzz_xxyyy = pbuffer.data(idx_g_hh + 720);

    auto gs_y_xyzzz_xxyyz = pbuffer.data(idx_g_hh + 721);

    auto gs_y_xyzzz_xxyzz = pbuffer.data(idx_g_hh + 722);

    auto gs_y_xyzzz_xxzzz = pbuffer.data(idx_g_hh + 723);

    auto gs_y_xyzzz_xyyyy = pbuffer.data(idx_g_hh + 724);

    auto gs_y_xyzzz_xyyyz = pbuffer.data(idx_g_hh + 725);

    auto gs_y_xyzzz_xyyzz = pbuffer.data(idx_g_hh + 726);

    auto gs_y_xyzzz_xyzzz = pbuffer.data(idx_g_hh + 727);

    auto gs_y_xyzzz_xzzzz = pbuffer.data(idx_g_hh + 728);

    auto gs_y_xyzzz_yyyyy = pbuffer.data(idx_g_hh + 729);

    auto gs_y_xyzzz_yyyyz = pbuffer.data(idx_g_hh + 730);

    auto gs_y_xyzzz_yyyzz = pbuffer.data(idx_g_hh + 731);

    auto gs_y_xyzzz_yyzzz = pbuffer.data(idx_g_hh + 732);

    auto gs_y_xyzzz_yzzzz = pbuffer.data(idx_g_hh + 733);

    auto gs_y_xyzzz_zzzzz = pbuffer.data(idx_g_hh + 734);

    #pragma omp simd aligned(gc_y, gs_y_xyzzz_xxxxx, gs_y_xyzzz_xxxxy, gs_y_xyzzz_xxxxz, gs_y_xyzzz_xxxyy, gs_y_xyzzz_xxxyz, gs_y_xyzzz_xxxzz, gs_y_xyzzz_xxyyy, gs_y_xyzzz_xxyyz, gs_y_xyzzz_xxyzz, gs_y_xyzzz_xxzzz, gs_y_xyzzz_xyyyy, gs_y_xyzzz_xyyyz, gs_y_xyzzz_xyyzz, gs_y_xyzzz_xyzzz, gs_y_xyzzz_xzzzz, gs_y_xyzzz_yyyyy, gs_y_xyzzz_yyyyz, gs_y_xyzzz_yyyzz, gs_y_xyzzz_yyzzz, gs_y_xyzzz_yzzzz, gs_y_xyzzz_zzzzz, ts_xyzzz_xxxx, ts_xyzzz_xxxxx, ts_xyzzz_xxxxy, ts_xyzzz_xxxxz, ts_xyzzz_xxxy, ts_xyzzz_xxxyy, ts_xyzzz_xxxyz, ts_xyzzz_xxxz, ts_xyzzz_xxxzz, ts_xyzzz_xxyy, ts_xyzzz_xxyyy, ts_xyzzz_xxyyz, ts_xyzzz_xxyz, ts_xyzzz_xxyzz, ts_xyzzz_xxzz, ts_xyzzz_xxzzz, ts_xyzzz_xyyy, ts_xyzzz_xyyyy, ts_xyzzz_xyyyz, ts_xyzzz_xyyz, ts_xyzzz_xyyzz, ts_xyzzz_xyzz, ts_xyzzz_xyzzz, ts_xyzzz_xzzz, ts_xyzzz_xzzzz, ts_xyzzz_yyyy, ts_xyzzz_yyyyy, ts_xyzzz_yyyyz, ts_xyzzz_yyyz, ts_xyzzz_yyyzz, ts_xyzzz_yyzz, ts_xyzzz_yyzzz, ts_xyzzz_yzzz, ts_xyzzz_yzzzz, ts_xyzzz_zzzz, ts_xyzzz_zzzzz, ts_xzzz_xxxxx, ts_xzzz_xxxxy, ts_xzzz_xxxxz, ts_xzzz_xxxyy, ts_xzzz_xxxyz, ts_xzzz_xxxzz, ts_xzzz_xxyyy, ts_xzzz_xxyyz, ts_xzzz_xxyzz, ts_xzzz_xxzzz, ts_xzzz_xyyyy, ts_xzzz_xyyyz, ts_xzzz_xyyzz, ts_xzzz_xyzzz, ts_xzzz_xzzzz, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyzz, ts_xzzz_yyzzz, ts_xzzz_yzzzz, ts_xzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzzz_xxxxx[i] = 2.0 * ts_xzzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxxxy[i] = 2.0 * ts_xzzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxxxz[i] = 2.0 * ts_xzzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxxyy[i] = 2.0 * ts_xzzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxxyz[i] = 2.0 * ts_xzzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxxzz[i] = 2.0 * ts_xzzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxyyy[i] = 2.0 * ts_xzzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxyyz[i] = 2.0 * ts_xzzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxyzz[i] = 2.0 * ts_xzzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xxzzz[i] = 2.0 * ts_xzzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyyyy[i] = 2.0 * ts_xzzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyyyz[i] = 2.0 * ts_xzzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyyzz[i] = 2.0 * ts_xzzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xyzzz[i] = 2.0 * ts_xzzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xzzzz[i] = 2.0 * ts_xzzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyyyy[i] = 2.0 * ts_xzzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyyyz[i] = 2.0 * ts_xzzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyyzz[i] = 2.0 * ts_xzzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yyzzz[i] = 2.0 * ts_xzzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yzzzz[i] = 2.0 * ts_xzzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_zzzzz[i] = 2.0 * ts_xzzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 735-756 components of targeted buffer : HH

    auto gs_y_xzzzz_xxxxx = pbuffer.data(idx_g_hh + 735);

    auto gs_y_xzzzz_xxxxy = pbuffer.data(idx_g_hh + 736);

    auto gs_y_xzzzz_xxxxz = pbuffer.data(idx_g_hh + 737);

    auto gs_y_xzzzz_xxxyy = pbuffer.data(idx_g_hh + 738);

    auto gs_y_xzzzz_xxxyz = pbuffer.data(idx_g_hh + 739);

    auto gs_y_xzzzz_xxxzz = pbuffer.data(idx_g_hh + 740);

    auto gs_y_xzzzz_xxyyy = pbuffer.data(idx_g_hh + 741);

    auto gs_y_xzzzz_xxyyz = pbuffer.data(idx_g_hh + 742);

    auto gs_y_xzzzz_xxyzz = pbuffer.data(idx_g_hh + 743);

    auto gs_y_xzzzz_xxzzz = pbuffer.data(idx_g_hh + 744);

    auto gs_y_xzzzz_xyyyy = pbuffer.data(idx_g_hh + 745);

    auto gs_y_xzzzz_xyyyz = pbuffer.data(idx_g_hh + 746);

    auto gs_y_xzzzz_xyyzz = pbuffer.data(idx_g_hh + 747);

    auto gs_y_xzzzz_xyzzz = pbuffer.data(idx_g_hh + 748);

    auto gs_y_xzzzz_xzzzz = pbuffer.data(idx_g_hh + 749);

    auto gs_y_xzzzz_yyyyy = pbuffer.data(idx_g_hh + 750);

    auto gs_y_xzzzz_yyyyz = pbuffer.data(idx_g_hh + 751);

    auto gs_y_xzzzz_yyyzz = pbuffer.data(idx_g_hh + 752);

    auto gs_y_xzzzz_yyzzz = pbuffer.data(idx_g_hh + 753);

    auto gs_y_xzzzz_yzzzz = pbuffer.data(idx_g_hh + 754);

    auto gs_y_xzzzz_zzzzz = pbuffer.data(idx_g_hh + 755);

    #pragma omp simd aligned(gc_y, gs_y_xzzzz_xxxxx, gs_y_xzzzz_xxxxy, gs_y_xzzzz_xxxxz, gs_y_xzzzz_xxxyy, gs_y_xzzzz_xxxyz, gs_y_xzzzz_xxxzz, gs_y_xzzzz_xxyyy, gs_y_xzzzz_xxyyz, gs_y_xzzzz_xxyzz, gs_y_xzzzz_xxzzz, gs_y_xzzzz_xyyyy, gs_y_xzzzz_xyyyz, gs_y_xzzzz_xyyzz, gs_y_xzzzz_xyzzz, gs_y_xzzzz_xzzzz, gs_y_xzzzz_yyyyy, gs_y_xzzzz_yyyyz, gs_y_xzzzz_yyyzz, gs_y_xzzzz_yyzzz, gs_y_xzzzz_yzzzz, gs_y_xzzzz_zzzzz, ts_xzzzz_xxxx, ts_xzzzz_xxxxx, ts_xzzzz_xxxxy, ts_xzzzz_xxxxz, ts_xzzzz_xxxy, ts_xzzzz_xxxyy, ts_xzzzz_xxxyz, ts_xzzzz_xxxz, ts_xzzzz_xxxzz, ts_xzzzz_xxyy, ts_xzzzz_xxyyy, ts_xzzzz_xxyyz, ts_xzzzz_xxyz, ts_xzzzz_xxyzz, ts_xzzzz_xxzz, ts_xzzzz_xxzzz, ts_xzzzz_xyyy, ts_xzzzz_xyyyy, ts_xzzzz_xyyyz, ts_xzzzz_xyyz, ts_xzzzz_xyyzz, ts_xzzzz_xyzz, ts_xzzzz_xyzzz, ts_xzzzz_xzzz, ts_xzzzz_xzzzz, ts_xzzzz_yyyy, ts_xzzzz_yyyyy, ts_xzzzz_yyyyz, ts_xzzzz_yyyz, ts_xzzzz_yyyzz, ts_xzzzz_yyzz, ts_xzzzz_yyzzz, ts_xzzzz_yzzz, ts_xzzzz_yzzzz, ts_xzzzz_zzzz, ts_xzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzzz_xxxxx[i] = 2.0 * ts_xzzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxxxy[i] = 2.0 * ts_xzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxxxz[i] = 2.0 * ts_xzzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxxyy[i] = 4.0 * ts_xzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxxyz[i] = 2.0 * ts_xzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxxzz[i] = 2.0 * ts_xzzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxyyy[i] = 6.0 * ts_xzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxyyz[i] = 4.0 * ts_xzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxyzz[i] = 2.0 * ts_xzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xxzzz[i] = 2.0 * ts_xzzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyyyy[i] = 8.0 * ts_xzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyyyz[i] = 6.0 * ts_xzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyyzz[i] = 4.0 * ts_xzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xyzzz[i] = 2.0 * ts_xzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xzzzz[i] = 2.0 * ts_xzzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyyyy[i] = 10.0 * ts_xzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyyyz[i] = 8.0 * ts_xzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyyzz[i] = 6.0 * ts_xzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yyzzz[i] = 4.0 * ts_xzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yzzzz[i] = 2.0 * ts_xzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_zzzzz[i] = 2.0 * ts_xzzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 756-777 components of targeted buffer : HH

    auto gs_y_yyyyy_xxxxx = pbuffer.data(idx_g_hh + 756);

    auto gs_y_yyyyy_xxxxy = pbuffer.data(idx_g_hh + 757);

    auto gs_y_yyyyy_xxxxz = pbuffer.data(idx_g_hh + 758);

    auto gs_y_yyyyy_xxxyy = pbuffer.data(idx_g_hh + 759);

    auto gs_y_yyyyy_xxxyz = pbuffer.data(idx_g_hh + 760);

    auto gs_y_yyyyy_xxxzz = pbuffer.data(idx_g_hh + 761);

    auto gs_y_yyyyy_xxyyy = pbuffer.data(idx_g_hh + 762);

    auto gs_y_yyyyy_xxyyz = pbuffer.data(idx_g_hh + 763);

    auto gs_y_yyyyy_xxyzz = pbuffer.data(idx_g_hh + 764);

    auto gs_y_yyyyy_xxzzz = pbuffer.data(idx_g_hh + 765);

    auto gs_y_yyyyy_xyyyy = pbuffer.data(idx_g_hh + 766);

    auto gs_y_yyyyy_xyyyz = pbuffer.data(idx_g_hh + 767);

    auto gs_y_yyyyy_xyyzz = pbuffer.data(idx_g_hh + 768);

    auto gs_y_yyyyy_xyzzz = pbuffer.data(idx_g_hh + 769);

    auto gs_y_yyyyy_xzzzz = pbuffer.data(idx_g_hh + 770);

    auto gs_y_yyyyy_yyyyy = pbuffer.data(idx_g_hh + 771);

    auto gs_y_yyyyy_yyyyz = pbuffer.data(idx_g_hh + 772);

    auto gs_y_yyyyy_yyyzz = pbuffer.data(idx_g_hh + 773);

    auto gs_y_yyyyy_yyzzz = pbuffer.data(idx_g_hh + 774);

    auto gs_y_yyyyy_yzzzz = pbuffer.data(idx_g_hh + 775);

    auto gs_y_yyyyy_zzzzz = pbuffer.data(idx_g_hh + 776);

    #pragma omp simd aligned(gc_y, gs_y_yyyyy_xxxxx, gs_y_yyyyy_xxxxy, gs_y_yyyyy_xxxxz, gs_y_yyyyy_xxxyy, gs_y_yyyyy_xxxyz, gs_y_yyyyy_xxxzz, gs_y_yyyyy_xxyyy, gs_y_yyyyy_xxyyz, gs_y_yyyyy_xxyzz, gs_y_yyyyy_xxzzz, gs_y_yyyyy_xyyyy, gs_y_yyyyy_xyyyz, gs_y_yyyyy_xyyzz, gs_y_yyyyy_xyzzz, gs_y_yyyyy_xzzzz, gs_y_yyyyy_yyyyy, gs_y_yyyyy_yyyyz, gs_y_yyyyy_yyyzz, gs_y_yyyyy_yyzzz, gs_y_yyyyy_yzzzz, gs_y_yyyyy_zzzzz, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxzz, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyzz, ts_yyyy_xxzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyzz, ts_yyyy_xyzzz, ts_yyyy_xzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyzz, ts_yyyy_yyzzz, ts_yyyy_yzzzz, ts_yyyy_zzzzz, ts_yyyyy_xxxx, ts_yyyyy_xxxxx, ts_yyyyy_xxxxy, ts_yyyyy_xxxxz, ts_yyyyy_xxxy, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxxz, ts_yyyyy_xxxzz, ts_yyyyy_xxyy, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyz, ts_yyyyy_xxyzz, ts_yyyyy_xxzz, ts_yyyyy_xxzzz, ts_yyyyy_xyyy, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzz, ts_yyyyy_xyzzz, ts_yyyyy_xzzz, ts_yyyyy_xzzzz, ts_yyyyy_yyyy, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzz, ts_yyyyy_yzzzz, ts_yyyyy_zzzz, ts_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyy_xxxxx[i] = 10.0 * ts_yyyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxxxy[i] = 10.0 * ts_yyyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxxxz[i] = 10.0 * ts_yyyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxxyy[i] = 10.0 * ts_yyyy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxxyz[i] = 10.0 * ts_yyyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxxzz[i] = 10.0 * ts_yyyy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxyyy[i] = 10.0 * ts_yyyy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxyyz[i] = 10.0 * ts_yyyy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxyzz[i] = 10.0 * ts_yyyy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xxzzz[i] = 10.0 * ts_yyyy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyyyy[i] = 10.0 * ts_yyyy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyyyz[i] = 10.0 * ts_yyyy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyyzz[i] = 10.0 * ts_yyyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xyzzz[i] = 10.0 * ts_yyyy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xzzzz[i] = 10.0 * ts_yyyy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyyyy[i] = 10.0 * ts_yyyy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyyyz[i] = 10.0 * ts_yyyy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyyzz[i] = 10.0 * ts_yyyy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yyzzz[i] = 10.0 * ts_yyyy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yzzzz[i] = 10.0 * ts_yyyy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_zzzzz[i] = 10.0 * ts_yyyy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 777-798 components of targeted buffer : HH

    auto gs_y_yyyyz_xxxxx = pbuffer.data(idx_g_hh + 777);

    auto gs_y_yyyyz_xxxxy = pbuffer.data(idx_g_hh + 778);

    auto gs_y_yyyyz_xxxxz = pbuffer.data(idx_g_hh + 779);

    auto gs_y_yyyyz_xxxyy = pbuffer.data(idx_g_hh + 780);

    auto gs_y_yyyyz_xxxyz = pbuffer.data(idx_g_hh + 781);

    auto gs_y_yyyyz_xxxzz = pbuffer.data(idx_g_hh + 782);

    auto gs_y_yyyyz_xxyyy = pbuffer.data(idx_g_hh + 783);

    auto gs_y_yyyyz_xxyyz = pbuffer.data(idx_g_hh + 784);

    auto gs_y_yyyyz_xxyzz = pbuffer.data(idx_g_hh + 785);

    auto gs_y_yyyyz_xxzzz = pbuffer.data(idx_g_hh + 786);

    auto gs_y_yyyyz_xyyyy = pbuffer.data(idx_g_hh + 787);

    auto gs_y_yyyyz_xyyyz = pbuffer.data(idx_g_hh + 788);

    auto gs_y_yyyyz_xyyzz = pbuffer.data(idx_g_hh + 789);

    auto gs_y_yyyyz_xyzzz = pbuffer.data(idx_g_hh + 790);

    auto gs_y_yyyyz_xzzzz = pbuffer.data(idx_g_hh + 791);

    auto gs_y_yyyyz_yyyyy = pbuffer.data(idx_g_hh + 792);

    auto gs_y_yyyyz_yyyyz = pbuffer.data(idx_g_hh + 793);

    auto gs_y_yyyyz_yyyzz = pbuffer.data(idx_g_hh + 794);

    auto gs_y_yyyyz_yyzzz = pbuffer.data(idx_g_hh + 795);

    auto gs_y_yyyyz_yzzzz = pbuffer.data(idx_g_hh + 796);

    auto gs_y_yyyyz_zzzzz = pbuffer.data(idx_g_hh + 797);

    #pragma omp simd aligned(gc_y, gs_y_yyyyz_xxxxx, gs_y_yyyyz_xxxxy, gs_y_yyyyz_xxxxz, gs_y_yyyyz_xxxyy, gs_y_yyyyz_xxxyz, gs_y_yyyyz_xxxzz, gs_y_yyyyz_xxyyy, gs_y_yyyyz_xxyyz, gs_y_yyyyz_xxyzz, gs_y_yyyyz_xxzzz, gs_y_yyyyz_xyyyy, gs_y_yyyyz_xyyyz, gs_y_yyyyz_xyyzz, gs_y_yyyyz_xyzzz, gs_y_yyyyz_xzzzz, gs_y_yyyyz_yyyyy, gs_y_yyyyz_yyyyz, gs_y_yyyyz_yyyzz, gs_y_yyyyz_yyzzz, gs_y_yyyyz_yzzzz, gs_y_yyyyz_zzzzz, ts_yyyyz_xxxx, ts_yyyyz_xxxxx, ts_yyyyz_xxxxy, ts_yyyyz_xxxxz, ts_yyyyz_xxxy, ts_yyyyz_xxxyy, ts_yyyyz_xxxyz, ts_yyyyz_xxxz, ts_yyyyz_xxxzz, ts_yyyyz_xxyy, ts_yyyyz_xxyyy, ts_yyyyz_xxyyz, ts_yyyyz_xxyz, ts_yyyyz_xxyzz, ts_yyyyz_xxzz, ts_yyyyz_xxzzz, ts_yyyyz_xyyy, ts_yyyyz_xyyyy, ts_yyyyz_xyyyz, ts_yyyyz_xyyz, ts_yyyyz_xyyzz, ts_yyyyz_xyzz, ts_yyyyz_xyzzz, ts_yyyyz_xzzz, ts_yyyyz_xzzzz, ts_yyyyz_yyyy, ts_yyyyz_yyyyy, ts_yyyyz_yyyyz, ts_yyyyz_yyyz, ts_yyyyz_yyyzz, ts_yyyyz_yyzz, ts_yyyyz_yyzzz, ts_yyyyz_yzzz, ts_yyyyz_yzzzz, ts_yyyyz_zzzz, ts_yyyyz_zzzzz, ts_yyyz_xxxxx, ts_yyyz_xxxxy, ts_yyyz_xxxxz, ts_yyyz_xxxyy, ts_yyyz_xxxyz, ts_yyyz_xxxzz, ts_yyyz_xxyyy, ts_yyyz_xxyyz, ts_yyyz_xxyzz, ts_yyyz_xxzzz, ts_yyyz_xyyyy, ts_yyyz_xyyyz, ts_yyyz_xyyzz, ts_yyyz_xyzzz, ts_yyyz_xzzzz, ts_yyyz_yyyyy, ts_yyyz_yyyyz, ts_yyyz_yyyzz, ts_yyyz_yyzzz, ts_yyyz_yzzzz, ts_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyz_xxxxx[i] = 8.0 * ts_yyyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxxxy[i] = 8.0 * ts_yyyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxxxz[i] = 8.0 * ts_yyyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxxyy[i] = 8.0 * ts_yyyz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxxyz[i] = 8.0 * ts_yyyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxxzz[i] = 8.0 * ts_yyyz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxyyy[i] = 8.0 * ts_yyyz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxyyz[i] = 8.0 * ts_yyyz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxyzz[i] = 8.0 * ts_yyyz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xxzzz[i] = 8.0 * ts_yyyz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyyyy[i] = 8.0 * ts_yyyz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyyyz[i] = 8.0 * ts_yyyz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyyzz[i] = 8.0 * ts_yyyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xyzzz[i] = 8.0 * ts_yyyz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xzzzz[i] = 8.0 * ts_yyyz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyyyy[i] = 8.0 * ts_yyyz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyyyz[i] = 8.0 * ts_yyyz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyyzz[i] = 8.0 * ts_yyyz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yyzzz[i] = 8.0 * ts_yyyz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yzzzz[i] = 8.0 * ts_yyyz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_zzzzz[i] = 8.0 * ts_yyyz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 798-819 components of targeted buffer : HH

    auto gs_y_yyyzz_xxxxx = pbuffer.data(idx_g_hh + 798);

    auto gs_y_yyyzz_xxxxy = pbuffer.data(idx_g_hh + 799);

    auto gs_y_yyyzz_xxxxz = pbuffer.data(idx_g_hh + 800);

    auto gs_y_yyyzz_xxxyy = pbuffer.data(idx_g_hh + 801);

    auto gs_y_yyyzz_xxxyz = pbuffer.data(idx_g_hh + 802);

    auto gs_y_yyyzz_xxxzz = pbuffer.data(idx_g_hh + 803);

    auto gs_y_yyyzz_xxyyy = pbuffer.data(idx_g_hh + 804);

    auto gs_y_yyyzz_xxyyz = pbuffer.data(idx_g_hh + 805);

    auto gs_y_yyyzz_xxyzz = pbuffer.data(idx_g_hh + 806);

    auto gs_y_yyyzz_xxzzz = pbuffer.data(idx_g_hh + 807);

    auto gs_y_yyyzz_xyyyy = pbuffer.data(idx_g_hh + 808);

    auto gs_y_yyyzz_xyyyz = pbuffer.data(idx_g_hh + 809);

    auto gs_y_yyyzz_xyyzz = pbuffer.data(idx_g_hh + 810);

    auto gs_y_yyyzz_xyzzz = pbuffer.data(idx_g_hh + 811);

    auto gs_y_yyyzz_xzzzz = pbuffer.data(idx_g_hh + 812);

    auto gs_y_yyyzz_yyyyy = pbuffer.data(idx_g_hh + 813);

    auto gs_y_yyyzz_yyyyz = pbuffer.data(idx_g_hh + 814);

    auto gs_y_yyyzz_yyyzz = pbuffer.data(idx_g_hh + 815);

    auto gs_y_yyyzz_yyzzz = pbuffer.data(idx_g_hh + 816);

    auto gs_y_yyyzz_yzzzz = pbuffer.data(idx_g_hh + 817);

    auto gs_y_yyyzz_zzzzz = pbuffer.data(idx_g_hh + 818);

    #pragma omp simd aligned(gc_y, gs_y_yyyzz_xxxxx, gs_y_yyyzz_xxxxy, gs_y_yyyzz_xxxxz, gs_y_yyyzz_xxxyy, gs_y_yyyzz_xxxyz, gs_y_yyyzz_xxxzz, gs_y_yyyzz_xxyyy, gs_y_yyyzz_xxyyz, gs_y_yyyzz_xxyzz, gs_y_yyyzz_xxzzz, gs_y_yyyzz_xyyyy, gs_y_yyyzz_xyyyz, gs_y_yyyzz_xyyzz, gs_y_yyyzz_xyzzz, gs_y_yyyzz_xzzzz, gs_y_yyyzz_yyyyy, gs_y_yyyzz_yyyyz, gs_y_yyyzz_yyyzz, gs_y_yyyzz_yyzzz, gs_y_yyyzz_yzzzz, gs_y_yyyzz_zzzzz, ts_yyyzz_xxxx, ts_yyyzz_xxxxx, ts_yyyzz_xxxxy, ts_yyyzz_xxxxz, ts_yyyzz_xxxy, ts_yyyzz_xxxyy, ts_yyyzz_xxxyz, ts_yyyzz_xxxz, ts_yyyzz_xxxzz, ts_yyyzz_xxyy, ts_yyyzz_xxyyy, ts_yyyzz_xxyyz, ts_yyyzz_xxyz, ts_yyyzz_xxyzz, ts_yyyzz_xxzz, ts_yyyzz_xxzzz, ts_yyyzz_xyyy, ts_yyyzz_xyyyy, ts_yyyzz_xyyyz, ts_yyyzz_xyyz, ts_yyyzz_xyyzz, ts_yyyzz_xyzz, ts_yyyzz_xyzzz, ts_yyyzz_xzzz, ts_yyyzz_xzzzz, ts_yyyzz_yyyy, ts_yyyzz_yyyyy, ts_yyyzz_yyyyz, ts_yyyzz_yyyz, ts_yyyzz_yyyzz, ts_yyyzz_yyzz, ts_yyyzz_yyzzz, ts_yyyzz_yzzz, ts_yyyzz_yzzzz, ts_yyyzz_zzzz, ts_yyyzz_zzzzz, ts_yyzz_xxxxx, ts_yyzz_xxxxy, ts_yyzz_xxxxz, ts_yyzz_xxxyy, ts_yyzz_xxxyz, ts_yyzz_xxxzz, ts_yyzz_xxyyy, ts_yyzz_xxyyz, ts_yyzz_xxyzz, ts_yyzz_xxzzz, ts_yyzz_xyyyy, ts_yyzz_xyyyz, ts_yyzz_xyyzz, ts_yyzz_xyzzz, ts_yyzz_xzzzz, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyzz, ts_yyzz_yyzzz, ts_yyzz_yzzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyzz_xxxxx[i] = 6.0 * ts_yyzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxxxy[i] = 6.0 * ts_yyzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxxxz[i] = 6.0 * ts_yyzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxxyy[i] = 6.0 * ts_yyzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxxyz[i] = 6.0 * ts_yyzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxxzz[i] = 6.0 * ts_yyzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxyyy[i] = 6.0 * ts_yyzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxyyz[i] = 6.0 * ts_yyzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxyzz[i] = 6.0 * ts_yyzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xxzzz[i] = 6.0 * ts_yyzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyyyy[i] = 6.0 * ts_yyzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyyyz[i] = 6.0 * ts_yyzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyyzz[i] = 6.0 * ts_yyzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xyzzz[i] = 6.0 * ts_yyzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xzzzz[i] = 6.0 * ts_yyzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyyyy[i] = 6.0 * ts_yyzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyyyz[i] = 6.0 * ts_yyzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyyzz[i] = 6.0 * ts_yyzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yyzzz[i] = 6.0 * ts_yyzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yzzzz[i] = 6.0 * ts_yyzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_zzzzz[i] = 6.0 * ts_yyzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 819-840 components of targeted buffer : HH

    auto gs_y_yyzzz_xxxxx = pbuffer.data(idx_g_hh + 819);

    auto gs_y_yyzzz_xxxxy = pbuffer.data(idx_g_hh + 820);

    auto gs_y_yyzzz_xxxxz = pbuffer.data(idx_g_hh + 821);

    auto gs_y_yyzzz_xxxyy = pbuffer.data(idx_g_hh + 822);

    auto gs_y_yyzzz_xxxyz = pbuffer.data(idx_g_hh + 823);

    auto gs_y_yyzzz_xxxzz = pbuffer.data(idx_g_hh + 824);

    auto gs_y_yyzzz_xxyyy = pbuffer.data(idx_g_hh + 825);

    auto gs_y_yyzzz_xxyyz = pbuffer.data(idx_g_hh + 826);

    auto gs_y_yyzzz_xxyzz = pbuffer.data(idx_g_hh + 827);

    auto gs_y_yyzzz_xxzzz = pbuffer.data(idx_g_hh + 828);

    auto gs_y_yyzzz_xyyyy = pbuffer.data(idx_g_hh + 829);

    auto gs_y_yyzzz_xyyyz = pbuffer.data(idx_g_hh + 830);

    auto gs_y_yyzzz_xyyzz = pbuffer.data(idx_g_hh + 831);

    auto gs_y_yyzzz_xyzzz = pbuffer.data(idx_g_hh + 832);

    auto gs_y_yyzzz_xzzzz = pbuffer.data(idx_g_hh + 833);

    auto gs_y_yyzzz_yyyyy = pbuffer.data(idx_g_hh + 834);

    auto gs_y_yyzzz_yyyyz = pbuffer.data(idx_g_hh + 835);

    auto gs_y_yyzzz_yyyzz = pbuffer.data(idx_g_hh + 836);

    auto gs_y_yyzzz_yyzzz = pbuffer.data(idx_g_hh + 837);

    auto gs_y_yyzzz_yzzzz = pbuffer.data(idx_g_hh + 838);

    auto gs_y_yyzzz_zzzzz = pbuffer.data(idx_g_hh + 839);

    #pragma omp simd aligned(gc_y, gs_y_yyzzz_xxxxx, gs_y_yyzzz_xxxxy, gs_y_yyzzz_xxxxz, gs_y_yyzzz_xxxyy, gs_y_yyzzz_xxxyz, gs_y_yyzzz_xxxzz, gs_y_yyzzz_xxyyy, gs_y_yyzzz_xxyyz, gs_y_yyzzz_xxyzz, gs_y_yyzzz_xxzzz, gs_y_yyzzz_xyyyy, gs_y_yyzzz_xyyyz, gs_y_yyzzz_xyyzz, gs_y_yyzzz_xyzzz, gs_y_yyzzz_xzzzz, gs_y_yyzzz_yyyyy, gs_y_yyzzz_yyyyz, gs_y_yyzzz_yyyzz, gs_y_yyzzz_yyzzz, gs_y_yyzzz_yzzzz, gs_y_yyzzz_zzzzz, ts_yyzzz_xxxx, ts_yyzzz_xxxxx, ts_yyzzz_xxxxy, ts_yyzzz_xxxxz, ts_yyzzz_xxxy, ts_yyzzz_xxxyy, ts_yyzzz_xxxyz, ts_yyzzz_xxxz, ts_yyzzz_xxxzz, ts_yyzzz_xxyy, ts_yyzzz_xxyyy, ts_yyzzz_xxyyz, ts_yyzzz_xxyz, ts_yyzzz_xxyzz, ts_yyzzz_xxzz, ts_yyzzz_xxzzz, ts_yyzzz_xyyy, ts_yyzzz_xyyyy, ts_yyzzz_xyyyz, ts_yyzzz_xyyz, ts_yyzzz_xyyzz, ts_yyzzz_xyzz, ts_yyzzz_xyzzz, ts_yyzzz_xzzz, ts_yyzzz_xzzzz, ts_yyzzz_yyyy, ts_yyzzz_yyyyy, ts_yyzzz_yyyyz, ts_yyzzz_yyyz, ts_yyzzz_yyyzz, ts_yyzzz_yyzz, ts_yyzzz_yyzzz, ts_yyzzz_yzzz, ts_yyzzz_yzzzz, ts_yyzzz_zzzz, ts_yyzzz_zzzzz, ts_yzzz_xxxxx, ts_yzzz_xxxxy, ts_yzzz_xxxxz, ts_yzzz_xxxyy, ts_yzzz_xxxyz, ts_yzzz_xxxzz, ts_yzzz_xxyyy, ts_yzzz_xxyyz, ts_yzzz_xxyzz, ts_yzzz_xxzzz, ts_yzzz_xyyyy, ts_yzzz_xyyyz, ts_yzzz_xyyzz, ts_yzzz_xyzzz, ts_yzzz_xzzzz, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyzz, ts_yzzz_yyzzz, ts_yzzz_yzzzz, ts_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzzz_xxxxx[i] = 4.0 * ts_yzzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxxxy[i] = 4.0 * ts_yzzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxxxz[i] = 4.0 * ts_yzzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxxyy[i] = 4.0 * ts_yzzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxxyz[i] = 4.0 * ts_yzzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxxzz[i] = 4.0 * ts_yzzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxyyy[i] = 4.0 * ts_yzzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxyyz[i] = 4.0 * ts_yzzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxyzz[i] = 4.0 * ts_yzzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xxzzz[i] = 4.0 * ts_yzzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyyyy[i] = 4.0 * ts_yzzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyyyz[i] = 4.0 * ts_yzzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyyzz[i] = 4.0 * ts_yzzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xyzzz[i] = 4.0 * ts_yzzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xzzzz[i] = 4.0 * ts_yzzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyyyy[i] = 4.0 * ts_yzzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyyyz[i] = 4.0 * ts_yzzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyyzz[i] = 4.0 * ts_yzzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yyzzz[i] = 4.0 * ts_yzzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yzzzz[i] = 4.0 * ts_yzzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_zzzzz[i] = 4.0 * ts_yzzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 840-861 components of targeted buffer : HH

    auto gs_y_yzzzz_xxxxx = pbuffer.data(idx_g_hh + 840);

    auto gs_y_yzzzz_xxxxy = pbuffer.data(idx_g_hh + 841);

    auto gs_y_yzzzz_xxxxz = pbuffer.data(idx_g_hh + 842);

    auto gs_y_yzzzz_xxxyy = pbuffer.data(idx_g_hh + 843);

    auto gs_y_yzzzz_xxxyz = pbuffer.data(idx_g_hh + 844);

    auto gs_y_yzzzz_xxxzz = pbuffer.data(idx_g_hh + 845);

    auto gs_y_yzzzz_xxyyy = pbuffer.data(idx_g_hh + 846);

    auto gs_y_yzzzz_xxyyz = pbuffer.data(idx_g_hh + 847);

    auto gs_y_yzzzz_xxyzz = pbuffer.data(idx_g_hh + 848);

    auto gs_y_yzzzz_xxzzz = pbuffer.data(idx_g_hh + 849);

    auto gs_y_yzzzz_xyyyy = pbuffer.data(idx_g_hh + 850);

    auto gs_y_yzzzz_xyyyz = pbuffer.data(idx_g_hh + 851);

    auto gs_y_yzzzz_xyyzz = pbuffer.data(idx_g_hh + 852);

    auto gs_y_yzzzz_xyzzz = pbuffer.data(idx_g_hh + 853);

    auto gs_y_yzzzz_xzzzz = pbuffer.data(idx_g_hh + 854);

    auto gs_y_yzzzz_yyyyy = pbuffer.data(idx_g_hh + 855);

    auto gs_y_yzzzz_yyyyz = pbuffer.data(idx_g_hh + 856);

    auto gs_y_yzzzz_yyyzz = pbuffer.data(idx_g_hh + 857);

    auto gs_y_yzzzz_yyzzz = pbuffer.data(idx_g_hh + 858);

    auto gs_y_yzzzz_yzzzz = pbuffer.data(idx_g_hh + 859);

    auto gs_y_yzzzz_zzzzz = pbuffer.data(idx_g_hh + 860);

    #pragma omp simd aligned(gc_y, gs_y_yzzzz_xxxxx, gs_y_yzzzz_xxxxy, gs_y_yzzzz_xxxxz, gs_y_yzzzz_xxxyy, gs_y_yzzzz_xxxyz, gs_y_yzzzz_xxxzz, gs_y_yzzzz_xxyyy, gs_y_yzzzz_xxyyz, gs_y_yzzzz_xxyzz, gs_y_yzzzz_xxzzz, gs_y_yzzzz_xyyyy, gs_y_yzzzz_xyyyz, gs_y_yzzzz_xyyzz, gs_y_yzzzz_xyzzz, gs_y_yzzzz_xzzzz, gs_y_yzzzz_yyyyy, gs_y_yzzzz_yyyyz, gs_y_yzzzz_yyyzz, gs_y_yzzzz_yyzzz, gs_y_yzzzz_yzzzz, gs_y_yzzzz_zzzzz, ts_yzzzz_xxxx, ts_yzzzz_xxxxx, ts_yzzzz_xxxxy, ts_yzzzz_xxxxz, ts_yzzzz_xxxy, ts_yzzzz_xxxyy, ts_yzzzz_xxxyz, ts_yzzzz_xxxz, ts_yzzzz_xxxzz, ts_yzzzz_xxyy, ts_yzzzz_xxyyy, ts_yzzzz_xxyyz, ts_yzzzz_xxyz, ts_yzzzz_xxyzz, ts_yzzzz_xxzz, ts_yzzzz_xxzzz, ts_yzzzz_xyyy, ts_yzzzz_xyyyy, ts_yzzzz_xyyyz, ts_yzzzz_xyyz, ts_yzzzz_xyyzz, ts_yzzzz_xyzz, ts_yzzzz_xyzzz, ts_yzzzz_xzzz, ts_yzzzz_xzzzz, ts_yzzzz_yyyy, ts_yzzzz_yyyyy, ts_yzzzz_yyyyz, ts_yzzzz_yyyz, ts_yzzzz_yyyzz, ts_yzzzz_yyzz, ts_yzzzz_yyzzz, ts_yzzzz_yzzz, ts_yzzzz_yzzzz, ts_yzzzz_zzzz, ts_yzzzz_zzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzzz_xxxxx[i] = 2.0 * ts_zzzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxxxy[i] = 2.0 * ts_zzzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxxxz[i] = 2.0 * ts_zzzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxxyy[i] = 2.0 * ts_zzzz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxxyz[i] = 2.0 * ts_zzzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxxzz[i] = 2.0 * ts_zzzz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxyyy[i] = 2.0 * ts_zzzz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxyyz[i] = 2.0 * ts_zzzz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxyzz[i] = 2.0 * ts_zzzz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xxzzz[i] = 2.0 * ts_zzzz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyyyy[i] = 2.0 * ts_zzzz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyyyz[i] = 2.0 * ts_zzzz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyyzz[i] = 2.0 * ts_zzzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xyzzz[i] = 2.0 * ts_zzzz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xzzzz[i] = 2.0 * ts_zzzz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyyyy[i] = 2.0 * ts_zzzz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyyyz[i] = 2.0 * ts_zzzz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyyzz[i] = 2.0 * ts_zzzz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yyzzz[i] = 2.0 * ts_zzzz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yzzzz[i] = 2.0 * ts_zzzz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_zzzzz[i] = 2.0 * ts_zzzz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 861-882 components of targeted buffer : HH

    auto gs_y_zzzzz_xxxxx = pbuffer.data(idx_g_hh + 861);

    auto gs_y_zzzzz_xxxxy = pbuffer.data(idx_g_hh + 862);

    auto gs_y_zzzzz_xxxxz = pbuffer.data(idx_g_hh + 863);

    auto gs_y_zzzzz_xxxyy = pbuffer.data(idx_g_hh + 864);

    auto gs_y_zzzzz_xxxyz = pbuffer.data(idx_g_hh + 865);

    auto gs_y_zzzzz_xxxzz = pbuffer.data(idx_g_hh + 866);

    auto gs_y_zzzzz_xxyyy = pbuffer.data(idx_g_hh + 867);

    auto gs_y_zzzzz_xxyyz = pbuffer.data(idx_g_hh + 868);

    auto gs_y_zzzzz_xxyzz = pbuffer.data(idx_g_hh + 869);

    auto gs_y_zzzzz_xxzzz = pbuffer.data(idx_g_hh + 870);

    auto gs_y_zzzzz_xyyyy = pbuffer.data(idx_g_hh + 871);

    auto gs_y_zzzzz_xyyyz = pbuffer.data(idx_g_hh + 872);

    auto gs_y_zzzzz_xyyzz = pbuffer.data(idx_g_hh + 873);

    auto gs_y_zzzzz_xyzzz = pbuffer.data(idx_g_hh + 874);

    auto gs_y_zzzzz_xzzzz = pbuffer.data(idx_g_hh + 875);

    auto gs_y_zzzzz_yyyyy = pbuffer.data(idx_g_hh + 876);

    auto gs_y_zzzzz_yyyyz = pbuffer.data(idx_g_hh + 877);

    auto gs_y_zzzzz_yyyzz = pbuffer.data(idx_g_hh + 878);

    auto gs_y_zzzzz_yyzzz = pbuffer.data(idx_g_hh + 879);

    auto gs_y_zzzzz_yzzzz = pbuffer.data(idx_g_hh + 880);

    auto gs_y_zzzzz_zzzzz = pbuffer.data(idx_g_hh + 881);

    #pragma omp simd aligned(gc_y, gs_y_zzzzz_xxxxx, gs_y_zzzzz_xxxxy, gs_y_zzzzz_xxxxz, gs_y_zzzzz_xxxyy, gs_y_zzzzz_xxxyz, gs_y_zzzzz_xxxzz, gs_y_zzzzz_xxyyy, gs_y_zzzzz_xxyyz, gs_y_zzzzz_xxyzz, gs_y_zzzzz_xxzzz, gs_y_zzzzz_xyyyy, gs_y_zzzzz_xyyyz, gs_y_zzzzz_xyyzz, gs_y_zzzzz_xyzzz, gs_y_zzzzz_xzzzz, gs_y_zzzzz_yyyyy, gs_y_zzzzz_yyyyz, gs_y_zzzzz_yyyzz, gs_y_zzzzz_yyzzz, gs_y_zzzzz_yzzzz, gs_y_zzzzz_zzzzz, ts_zzzzz_xxxx, ts_zzzzz_xxxxx, ts_zzzzz_xxxxy, ts_zzzzz_xxxxz, ts_zzzzz_xxxy, ts_zzzzz_xxxyy, ts_zzzzz_xxxyz, ts_zzzzz_xxxz, ts_zzzzz_xxxzz, ts_zzzzz_xxyy, ts_zzzzz_xxyyy, ts_zzzzz_xxyyz, ts_zzzzz_xxyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyy, ts_zzzzz_xyyyy, ts_zzzzz_xyyyz, ts_zzzzz_xyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyy, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzzz_xxxxx[i] = 2.0 * ts_zzzzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxxxy[i] = 2.0 * ts_zzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxxxz[i] = 2.0 * ts_zzzzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxxyy[i] = 4.0 * ts_zzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxxyz[i] = 2.0 * ts_zzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxxzz[i] = 2.0 * ts_zzzzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxyyy[i] = 6.0 * ts_zzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxyyz[i] = 4.0 * ts_zzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxyzz[i] = 2.0 * ts_zzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xxzzz[i] = 2.0 * ts_zzzzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyyyy[i] = 8.0 * ts_zzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyyyz[i] = 6.0 * ts_zzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyyzz[i] = 4.0 * ts_zzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xyzzz[i] = 2.0 * ts_zzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xzzzz[i] = 2.0 * ts_zzzzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyyyy[i] = 10.0 * ts_zzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyyyz[i] = 8.0 * ts_zzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyyzz[i] = 6.0 * ts_zzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yyzzz[i] = 4.0 * ts_zzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yzzzz[i] = 2.0 * ts_zzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_zzzzz[i] = 2.0 * ts_zzzzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 882-903 components of targeted buffer : HH

    auto gs_z_xxxxx_xxxxx = pbuffer.data(idx_g_hh + 882);

    auto gs_z_xxxxx_xxxxy = pbuffer.data(idx_g_hh + 883);

    auto gs_z_xxxxx_xxxxz = pbuffer.data(idx_g_hh + 884);

    auto gs_z_xxxxx_xxxyy = pbuffer.data(idx_g_hh + 885);

    auto gs_z_xxxxx_xxxyz = pbuffer.data(idx_g_hh + 886);

    auto gs_z_xxxxx_xxxzz = pbuffer.data(idx_g_hh + 887);

    auto gs_z_xxxxx_xxyyy = pbuffer.data(idx_g_hh + 888);

    auto gs_z_xxxxx_xxyyz = pbuffer.data(idx_g_hh + 889);

    auto gs_z_xxxxx_xxyzz = pbuffer.data(idx_g_hh + 890);

    auto gs_z_xxxxx_xxzzz = pbuffer.data(idx_g_hh + 891);

    auto gs_z_xxxxx_xyyyy = pbuffer.data(idx_g_hh + 892);

    auto gs_z_xxxxx_xyyyz = pbuffer.data(idx_g_hh + 893);

    auto gs_z_xxxxx_xyyzz = pbuffer.data(idx_g_hh + 894);

    auto gs_z_xxxxx_xyzzz = pbuffer.data(idx_g_hh + 895);

    auto gs_z_xxxxx_xzzzz = pbuffer.data(idx_g_hh + 896);

    auto gs_z_xxxxx_yyyyy = pbuffer.data(idx_g_hh + 897);

    auto gs_z_xxxxx_yyyyz = pbuffer.data(idx_g_hh + 898);

    auto gs_z_xxxxx_yyyzz = pbuffer.data(idx_g_hh + 899);

    auto gs_z_xxxxx_yyzzz = pbuffer.data(idx_g_hh + 900);

    auto gs_z_xxxxx_yzzzz = pbuffer.data(idx_g_hh + 901);

    auto gs_z_xxxxx_zzzzz = pbuffer.data(idx_g_hh + 902);

    #pragma omp simd aligned(gc_z, gs_z_xxxxx_xxxxx, gs_z_xxxxx_xxxxy, gs_z_xxxxx_xxxxz, gs_z_xxxxx_xxxyy, gs_z_xxxxx_xxxyz, gs_z_xxxxx_xxxzz, gs_z_xxxxx_xxyyy, gs_z_xxxxx_xxyyz, gs_z_xxxxx_xxyzz, gs_z_xxxxx_xxzzz, gs_z_xxxxx_xyyyy, gs_z_xxxxx_xyyyz, gs_z_xxxxx_xyyzz, gs_z_xxxxx_xyzzz, gs_z_xxxxx_xzzzz, gs_z_xxxxx_yyyyy, gs_z_xxxxx_yyyyz, gs_z_xxxxx_yyyzz, gs_z_xxxxx_yyzzz, gs_z_xxxxx_yzzzz, gs_z_xxxxx_zzzzz, ts_xxxxx_xxxx, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxxz, ts_xxxxx_xxxy, ts_xxxxx_xxxyy, ts_xxxxx_xxxyz, ts_xxxxx_xxxz, ts_xxxxx_xxxzz, ts_xxxxx_xxyy, ts_xxxxx_xxyyy, ts_xxxxx_xxyyz, ts_xxxxx_xxyz, ts_xxxxx_xxyzz, ts_xxxxx_xxzz, ts_xxxxx_xxzzz, ts_xxxxx_xyyy, ts_xxxxx_xyyyy, ts_xxxxx_xyyyz, ts_xxxxx_xyyz, ts_xxxxx_xyyzz, ts_xxxxx_xyzz, ts_xxxxx_xyzzz, ts_xxxxx_xzzz, ts_xxxxx_xzzzz, ts_xxxxx_yyyy, ts_xxxxx_yyyyy, ts_xxxxx_yyyyz, ts_xxxxx_yyyz, ts_xxxxx_yyyzz, ts_xxxxx_yyzz, ts_xxxxx_yyzzz, ts_xxxxx_yzzz, ts_xxxxx_yzzzz, ts_xxxxx_zzzz, ts_xxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxx_xxxxx[i] = 2.0 * ts_xxxxx_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxxxy[i] = 2.0 * ts_xxxxx_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxxxz[i] = 2.0 * ts_xxxxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxxyy[i] = 2.0 * ts_xxxxx_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxxyz[i] = 2.0 * ts_xxxxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxxzz[i] = 4.0 * ts_xxxxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxyyy[i] = 2.0 * ts_xxxxx_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxyyz[i] = 2.0 * ts_xxxxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxyzz[i] = 4.0 * ts_xxxxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xxzzz[i] = 6.0 * ts_xxxxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyyyy[i] = 2.0 * ts_xxxxx_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyyyz[i] = 2.0 * ts_xxxxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyyzz[i] = 4.0 * ts_xxxxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xyzzz[i] = 6.0 * ts_xxxxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xzzzz[i] = 8.0 * ts_xxxxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyyyy[i] = 2.0 * ts_xxxxx_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyyyz[i] = 2.0 * ts_xxxxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyyzz[i] = 4.0 * ts_xxxxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yyzzz[i] = 6.0 * ts_xxxxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yzzzz[i] = 8.0 * ts_xxxxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_zzzzz[i] = 10.0 * ts_xxxxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 903-924 components of targeted buffer : HH

    auto gs_z_xxxxy_xxxxx = pbuffer.data(idx_g_hh + 903);

    auto gs_z_xxxxy_xxxxy = pbuffer.data(idx_g_hh + 904);

    auto gs_z_xxxxy_xxxxz = pbuffer.data(idx_g_hh + 905);

    auto gs_z_xxxxy_xxxyy = pbuffer.data(idx_g_hh + 906);

    auto gs_z_xxxxy_xxxyz = pbuffer.data(idx_g_hh + 907);

    auto gs_z_xxxxy_xxxzz = pbuffer.data(idx_g_hh + 908);

    auto gs_z_xxxxy_xxyyy = pbuffer.data(idx_g_hh + 909);

    auto gs_z_xxxxy_xxyyz = pbuffer.data(idx_g_hh + 910);

    auto gs_z_xxxxy_xxyzz = pbuffer.data(idx_g_hh + 911);

    auto gs_z_xxxxy_xxzzz = pbuffer.data(idx_g_hh + 912);

    auto gs_z_xxxxy_xyyyy = pbuffer.data(idx_g_hh + 913);

    auto gs_z_xxxxy_xyyyz = pbuffer.data(idx_g_hh + 914);

    auto gs_z_xxxxy_xyyzz = pbuffer.data(idx_g_hh + 915);

    auto gs_z_xxxxy_xyzzz = pbuffer.data(idx_g_hh + 916);

    auto gs_z_xxxxy_xzzzz = pbuffer.data(idx_g_hh + 917);

    auto gs_z_xxxxy_yyyyy = pbuffer.data(idx_g_hh + 918);

    auto gs_z_xxxxy_yyyyz = pbuffer.data(idx_g_hh + 919);

    auto gs_z_xxxxy_yyyzz = pbuffer.data(idx_g_hh + 920);

    auto gs_z_xxxxy_yyzzz = pbuffer.data(idx_g_hh + 921);

    auto gs_z_xxxxy_yzzzz = pbuffer.data(idx_g_hh + 922);

    auto gs_z_xxxxy_zzzzz = pbuffer.data(idx_g_hh + 923);

    #pragma omp simd aligned(gc_z, gs_z_xxxxy_xxxxx, gs_z_xxxxy_xxxxy, gs_z_xxxxy_xxxxz, gs_z_xxxxy_xxxyy, gs_z_xxxxy_xxxyz, gs_z_xxxxy_xxxzz, gs_z_xxxxy_xxyyy, gs_z_xxxxy_xxyyz, gs_z_xxxxy_xxyzz, gs_z_xxxxy_xxzzz, gs_z_xxxxy_xyyyy, gs_z_xxxxy_xyyyz, gs_z_xxxxy_xyyzz, gs_z_xxxxy_xyzzz, gs_z_xxxxy_xzzzz, gs_z_xxxxy_yyyyy, gs_z_xxxxy_yyyyz, gs_z_xxxxy_yyyzz, gs_z_xxxxy_yyzzz, gs_z_xxxxy_yzzzz, gs_z_xxxxy_zzzzz, ts_xxxxy_xxxx, ts_xxxxy_xxxxx, ts_xxxxy_xxxxy, ts_xxxxy_xxxxz, ts_xxxxy_xxxy, ts_xxxxy_xxxyy, ts_xxxxy_xxxyz, ts_xxxxy_xxxz, ts_xxxxy_xxxzz, ts_xxxxy_xxyy, ts_xxxxy_xxyyy, ts_xxxxy_xxyyz, ts_xxxxy_xxyz, ts_xxxxy_xxyzz, ts_xxxxy_xxzz, ts_xxxxy_xxzzz, ts_xxxxy_xyyy, ts_xxxxy_xyyyy, ts_xxxxy_xyyyz, ts_xxxxy_xyyz, ts_xxxxy_xyyzz, ts_xxxxy_xyzz, ts_xxxxy_xyzzz, ts_xxxxy_xzzz, ts_xxxxy_xzzzz, ts_xxxxy_yyyy, ts_xxxxy_yyyyy, ts_xxxxy_yyyyz, ts_xxxxy_yyyz, ts_xxxxy_yyyzz, ts_xxxxy_yyzz, ts_xxxxy_yyzzz, ts_xxxxy_yzzz, ts_xxxxy_yzzzz, ts_xxxxy_zzzz, ts_xxxxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxy_xxxxx[i] = 2.0 * ts_xxxxy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxxxy[i] = 2.0 * ts_xxxxy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxxxz[i] = 2.0 * ts_xxxxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxxyy[i] = 2.0 * ts_xxxxy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxxyz[i] = 2.0 * ts_xxxxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxxzz[i] = 4.0 * ts_xxxxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxyyy[i] = 2.0 * ts_xxxxy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxyyz[i] = 2.0 * ts_xxxxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxyzz[i] = 4.0 * ts_xxxxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xxzzz[i] = 6.0 * ts_xxxxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyyyy[i] = 2.0 * ts_xxxxy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyyyz[i] = 2.0 * ts_xxxxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyyzz[i] = 4.0 * ts_xxxxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xyzzz[i] = 6.0 * ts_xxxxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xzzzz[i] = 8.0 * ts_xxxxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyyyy[i] = 2.0 * ts_xxxxy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyyyz[i] = 2.0 * ts_xxxxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyyzz[i] = 4.0 * ts_xxxxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yyzzz[i] = 6.0 * ts_xxxxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yzzzz[i] = 8.0 * ts_xxxxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_zzzzz[i] = 10.0 * ts_xxxxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 924-945 components of targeted buffer : HH

    auto gs_z_xxxxz_xxxxx = pbuffer.data(idx_g_hh + 924);

    auto gs_z_xxxxz_xxxxy = pbuffer.data(idx_g_hh + 925);

    auto gs_z_xxxxz_xxxxz = pbuffer.data(idx_g_hh + 926);

    auto gs_z_xxxxz_xxxyy = pbuffer.data(idx_g_hh + 927);

    auto gs_z_xxxxz_xxxyz = pbuffer.data(idx_g_hh + 928);

    auto gs_z_xxxxz_xxxzz = pbuffer.data(idx_g_hh + 929);

    auto gs_z_xxxxz_xxyyy = pbuffer.data(idx_g_hh + 930);

    auto gs_z_xxxxz_xxyyz = pbuffer.data(idx_g_hh + 931);

    auto gs_z_xxxxz_xxyzz = pbuffer.data(idx_g_hh + 932);

    auto gs_z_xxxxz_xxzzz = pbuffer.data(idx_g_hh + 933);

    auto gs_z_xxxxz_xyyyy = pbuffer.data(idx_g_hh + 934);

    auto gs_z_xxxxz_xyyyz = pbuffer.data(idx_g_hh + 935);

    auto gs_z_xxxxz_xyyzz = pbuffer.data(idx_g_hh + 936);

    auto gs_z_xxxxz_xyzzz = pbuffer.data(idx_g_hh + 937);

    auto gs_z_xxxxz_xzzzz = pbuffer.data(idx_g_hh + 938);

    auto gs_z_xxxxz_yyyyy = pbuffer.data(idx_g_hh + 939);

    auto gs_z_xxxxz_yyyyz = pbuffer.data(idx_g_hh + 940);

    auto gs_z_xxxxz_yyyzz = pbuffer.data(idx_g_hh + 941);

    auto gs_z_xxxxz_yyzzz = pbuffer.data(idx_g_hh + 942);

    auto gs_z_xxxxz_yzzzz = pbuffer.data(idx_g_hh + 943);

    auto gs_z_xxxxz_zzzzz = pbuffer.data(idx_g_hh + 944);

    #pragma omp simd aligned(gc_z, gs_z_xxxxz_xxxxx, gs_z_xxxxz_xxxxy, gs_z_xxxxz_xxxxz, gs_z_xxxxz_xxxyy, gs_z_xxxxz_xxxyz, gs_z_xxxxz_xxxzz, gs_z_xxxxz_xxyyy, gs_z_xxxxz_xxyyz, gs_z_xxxxz_xxyzz, gs_z_xxxxz_xxzzz, gs_z_xxxxz_xyyyy, gs_z_xxxxz_xyyyz, gs_z_xxxxz_xyyzz, gs_z_xxxxz_xyzzz, gs_z_xxxxz_xzzzz, gs_z_xxxxz_yyyyy, gs_z_xxxxz_yyyyz, gs_z_xxxxz_yyyzz, gs_z_xxxxz_yyzzz, gs_z_xxxxz_yzzzz, gs_z_xxxxz_zzzzz, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxzz, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyzz, ts_xxxx_xxzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyzz, ts_xxxx_xyzzz, ts_xxxx_xzzzz, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyzz, ts_xxxx_yyzzz, ts_xxxx_yzzzz, ts_xxxx_zzzzz, ts_xxxxz_xxxx, ts_xxxxz_xxxxx, ts_xxxxz_xxxxy, ts_xxxxz_xxxxz, ts_xxxxz_xxxy, ts_xxxxz_xxxyy, ts_xxxxz_xxxyz, ts_xxxxz_xxxz, ts_xxxxz_xxxzz, ts_xxxxz_xxyy, ts_xxxxz_xxyyy, ts_xxxxz_xxyyz, ts_xxxxz_xxyz, ts_xxxxz_xxyzz, ts_xxxxz_xxzz, ts_xxxxz_xxzzz, ts_xxxxz_xyyy, ts_xxxxz_xyyyy, ts_xxxxz_xyyyz, ts_xxxxz_xyyz, ts_xxxxz_xyyzz, ts_xxxxz_xyzz, ts_xxxxz_xyzzz, ts_xxxxz_xzzz, ts_xxxxz_xzzzz, ts_xxxxz_yyyy, ts_xxxxz_yyyyy, ts_xxxxz_yyyyz, ts_xxxxz_yyyz, ts_xxxxz_yyyzz, ts_xxxxz_yyzz, ts_xxxxz_yyzzz, ts_xxxxz_yzzz, ts_xxxxz_yzzzz, ts_xxxxz_zzzz, ts_xxxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxz_xxxxx[i] = 2.0 * ts_xxxx_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxxxy[i] = 2.0 * ts_xxxx_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxxxz[i] = 2.0 * ts_xxxx_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxxyy[i] = 2.0 * ts_xxxx_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxxyz[i] = 2.0 * ts_xxxx_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxxzz[i] = 2.0 * ts_xxxx_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxyyy[i] = 2.0 * ts_xxxx_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxyyz[i] = 2.0 * ts_xxxx_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxyzz[i] = 2.0 * ts_xxxx_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xxzzz[i] = 2.0 * ts_xxxx_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyyyy[i] = 2.0 * ts_xxxx_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyyyz[i] = 2.0 * ts_xxxx_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyyzz[i] = 2.0 * ts_xxxx_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xyzzz[i] = 2.0 * ts_xxxx_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xzzzz[i] = 2.0 * ts_xxxx_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyyyy[i] = 2.0 * ts_xxxx_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyyyz[i] = 2.0 * ts_xxxx_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyyzz[i] = 2.0 * ts_xxxx_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yyzzz[i] = 2.0 * ts_xxxx_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yzzzz[i] = 2.0 * ts_xxxx_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_zzzzz[i] = 2.0 * ts_xxxx_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxxxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 945-966 components of targeted buffer : HH

    auto gs_z_xxxyy_xxxxx = pbuffer.data(idx_g_hh + 945);

    auto gs_z_xxxyy_xxxxy = pbuffer.data(idx_g_hh + 946);

    auto gs_z_xxxyy_xxxxz = pbuffer.data(idx_g_hh + 947);

    auto gs_z_xxxyy_xxxyy = pbuffer.data(idx_g_hh + 948);

    auto gs_z_xxxyy_xxxyz = pbuffer.data(idx_g_hh + 949);

    auto gs_z_xxxyy_xxxzz = pbuffer.data(idx_g_hh + 950);

    auto gs_z_xxxyy_xxyyy = pbuffer.data(idx_g_hh + 951);

    auto gs_z_xxxyy_xxyyz = pbuffer.data(idx_g_hh + 952);

    auto gs_z_xxxyy_xxyzz = pbuffer.data(idx_g_hh + 953);

    auto gs_z_xxxyy_xxzzz = pbuffer.data(idx_g_hh + 954);

    auto gs_z_xxxyy_xyyyy = pbuffer.data(idx_g_hh + 955);

    auto gs_z_xxxyy_xyyyz = pbuffer.data(idx_g_hh + 956);

    auto gs_z_xxxyy_xyyzz = pbuffer.data(idx_g_hh + 957);

    auto gs_z_xxxyy_xyzzz = pbuffer.data(idx_g_hh + 958);

    auto gs_z_xxxyy_xzzzz = pbuffer.data(idx_g_hh + 959);

    auto gs_z_xxxyy_yyyyy = pbuffer.data(idx_g_hh + 960);

    auto gs_z_xxxyy_yyyyz = pbuffer.data(idx_g_hh + 961);

    auto gs_z_xxxyy_yyyzz = pbuffer.data(idx_g_hh + 962);

    auto gs_z_xxxyy_yyzzz = pbuffer.data(idx_g_hh + 963);

    auto gs_z_xxxyy_yzzzz = pbuffer.data(idx_g_hh + 964);

    auto gs_z_xxxyy_zzzzz = pbuffer.data(idx_g_hh + 965);

    #pragma omp simd aligned(gc_z, gs_z_xxxyy_xxxxx, gs_z_xxxyy_xxxxy, gs_z_xxxyy_xxxxz, gs_z_xxxyy_xxxyy, gs_z_xxxyy_xxxyz, gs_z_xxxyy_xxxzz, gs_z_xxxyy_xxyyy, gs_z_xxxyy_xxyyz, gs_z_xxxyy_xxyzz, gs_z_xxxyy_xxzzz, gs_z_xxxyy_xyyyy, gs_z_xxxyy_xyyyz, gs_z_xxxyy_xyyzz, gs_z_xxxyy_xyzzz, gs_z_xxxyy_xzzzz, gs_z_xxxyy_yyyyy, gs_z_xxxyy_yyyyz, gs_z_xxxyy_yyyzz, gs_z_xxxyy_yyzzz, gs_z_xxxyy_yzzzz, gs_z_xxxyy_zzzzz, ts_xxxyy_xxxx, ts_xxxyy_xxxxx, ts_xxxyy_xxxxy, ts_xxxyy_xxxxz, ts_xxxyy_xxxy, ts_xxxyy_xxxyy, ts_xxxyy_xxxyz, ts_xxxyy_xxxz, ts_xxxyy_xxxzz, ts_xxxyy_xxyy, ts_xxxyy_xxyyy, ts_xxxyy_xxyyz, ts_xxxyy_xxyz, ts_xxxyy_xxyzz, ts_xxxyy_xxzz, ts_xxxyy_xxzzz, ts_xxxyy_xyyy, ts_xxxyy_xyyyy, ts_xxxyy_xyyyz, ts_xxxyy_xyyz, ts_xxxyy_xyyzz, ts_xxxyy_xyzz, ts_xxxyy_xyzzz, ts_xxxyy_xzzz, ts_xxxyy_xzzzz, ts_xxxyy_yyyy, ts_xxxyy_yyyyy, ts_xxxyy_yyyyz, ts_xxxyy_yyyz, ts_xxxyy_yyyzz, ts_xxxyy_yyzz, ts_xxxyy_yyzzz, ts_xxxyy_yzzz, ts_xxxyy_yzzzz, ts_xxxyy_zzzz, ts_xxxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyy_xxxxx[i] = 2.0 * ts_xxxyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxxxy[i] = 2.0 * ts_xxxyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxxxz[i] = 2.0 * ts_xxxyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxxyy[i] = 2.0 * ts_xxxyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxxyz[i] = 2.0 * ts_xxxyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxxzz[i] = 4.0 * ts_xxxyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxyyy[i] = 2.0 * ts_xxxyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxyyz[i] = 2.0 * ts_xxxyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxyzz[i] = 4.0 * ts_xxxyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xxzzz[i] = 6.0 * ts_xxxyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyyyy[i] = 2.0 * ts_xxxyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyyyz[i] = 2.0 * ts_xxxyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyyzz[i] = 4.0 * ts_xxxyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xyzzz[i] = 6.0 * ts_xxxyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xzzzz[i] = 8.0 * ts_xxxyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyyyy[i] = 2.0 * ts_xxxyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyyyz[i] = 2.0 * ts_xxxyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyyzz[i] = 4.0 * ts_xxxyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yyzzz[i] = 6.0 * ts_xxxyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yzzzz[i] = 8.0 * ts_xxxyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_zzzzz[i] = 10.0 * ts_xxxyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 966-987 components of targeted buffer : HH

    auto gs_z_xxxyz_xxxxx = pbuffer.data(idx_g_hh + 966);

    auto gs_z_xxxyz_xxxxy = pbuffer.data(idx_g_hh + 967);

    auto gs_z_xxxyz_xxxxz = pbuffer.data(idx_g_hh + 968);

    auto gs_z_xxxyz_xxxyy = pbuffer.data(idx_g_hh + 969);

    auto gs_z_xxxyz_xxxyz = pbuffer.data(idx_g_hh + 970);

    auto gs_z_xxxyz_xxxzz = pbuffer.data(idx_g_hh + 971);

    auto gs_z_xxxyz_xxyyy = pbuffer.data(idx_g_hh + 972);

    auto gs_z_xxxyz_xxyyz = pbuffer.data(idx_g_hh + 973);

    auto gs_z_xxxyz_xxyzz = pbuffer.data(idx_g_hh + 974);

    auto gs_z_xxxyz_xxzzz = pbuffer.data(idx_g_hh + 975);

    auto gs_z_xxxyz_xyyyy = pbuffer.data(idx_g_hh + 976);

    auto gs_z_xxxyz_xyyyz = pbuffer.data(idx_g_hh + 977);

    auto gs_z_xxxyz_xyyzz = pbuffer.data(idx_g_hh + 978);

    auto gs_z_xxxyz_xyzzz = pbuffer.data(idx_g_hh + 979);

    auto gs_z_xxxyz_xzzzz = pbuffer.data(idx_g_hh + 980);

    auto gs_z_xxxyz_yyyyy = pbuffer.data(idx_g_hh + 981);

    auto gs_z_xxxyz_yyyyz = pbuffer.data(idx_g_hh + 982);

    auto gs_z_xxxyz_yyyzz = pbuffer.data(idx_g_hh + 983);

    auto gs_z_xxxyz_yyzzz = pbuffer.data(idx_g_hh + 984);

    auto gs_z_xxxyz_yzzzz = pbuffer.data(idx_g_hh + 985);

    auto gs_z_xxxyz_zzzzz = pbuffer.data(idx_g_hh + 986);

    #pragma omp simd aligned(gc_z, gs_z_xxxyz_xxxxx, gs_z_xxxyz_xxxxy, gs_z_xxxyz_xxxxz, gs_z_xxxyz_xxxyy, gs_z_xxxyz_xxxyz, gs_z_xxxyz_xxxzz, gs_z_xxxyz_xxyyy, gs_z_xxxyz_xxyyz, gs_z_xxxyz_xxyzz, gs_z_xxxyz_xxzzz, gs_z_xxxyz_xyyyy, gs_z_xxxyz_xyyyz, gs_z_xxxyz_xyyzz, gs_z_xxxyz_xyzzz, gs_z_xxxyz_xzzzz, gs_z_xxxyz_yyyyy, gs_z_xxxyz_yyyyz, gs_z_xxxyz_yyyzz, gs_z_xxxyz_yyzzz, gs_z_xxxyz_yzzzz, gs_z_xxxyz_zzzzz, ts_xxxy_xxxxx, ts_xxxy_xxxxy, ts_xxxy_xxxxz, ts_xxxy_xxxyy, ts_xxxy_xxxyz, ts_xxxy_xxxzz, ts_xxxy_xxyyy, ts_xxxy_xxyyz, ts_xxxy_xxyzz, ts_xxxy_xxzzz, ts_xxxy_xyyyy, ts_xxxy_xyyyz, ts_xxxy_xyyzz, ts_xxxy_xyzzz, ts_xxxy_xzzzz, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyzz, ts_xxxy_yyzzz, ts_xxxy_yzzzz, ts_xxxy_zzzzz, ts_xxxyz_xxxx, ts_xxxyz_xxxxx, ts_xxxyz_xxxxy, ts_xxxyz_xxxxz, ts_xxxyz_xxxy, ts_xxxyz_xxxyy, ts_xxxyz_xxxyz, ts_xxxyz_xxxz, ts_xxxyz_xxxzz, ts_xxxyz_xxyy, ts_xxxyz_xxyyy, ts_xxxyz_xxyyz, ts_xxxyz_xxyz, ts_xxxyz_xxyzz, ts_xxxyz_xxzz, ts_xxxyz_xxzzz, ts_xxxyz_xyyy, ts_xxxyz_xyyyy, ts_xxxyz_xyyyz, ts_xxxyz_xyyz, ts_xxxyz_xyyzz, ts_xxxyz_xyzz, ts_xxxyz_xyzzz, ts_xxxyz_xzzz, ts_xxxyz_xzzzz, ts_xxxyz_yyyy, ts_xxxyz_yyyyy, ts_xxxyz_yyyyz, ts_xxxyz_yyyz, ts_xxxyz_yyyzz, ts_xxxyz_yyzz, ts_xxxyz_yyzzz, ts_xxxyz_yzzz, ts_xxxyz_yzzzz, ts_xxxyz_zzzz, ts_xxxyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyz_xxxxx[i] = 2.0 * ts_xxxy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxxxy[i] = 2.0 * ts_xxxy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxxxz[i] = 2.0 * ts_xxxy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxxyy[i] = 2.0 * ts_xxxy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxxyz[i] = 2.0 * ts_xxxy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxxzz[i] = 2.0 * ts_xxxy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxyyy[i] = 2.0 * ts_xxxy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxyyz[i] = 2.0 * ts_xxxy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxyzz[i] = 2.0 * ts_xxxy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xxzzz[i] = 2.0 * ts_xxxy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyyyy[i] = 2.0 * ts_xxxy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyyyz[i] = 2.0 * ts_xxxy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyyzz[i] = 2.0 * ts_xxxy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xyzzz[i] = 2.0 * ts_xxxy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xzzzz[i] = 2.0 * ts_xxxy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyyyy[i] = 2.0 * ts_xxxy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyyyz[i] = 2.0 * ts_xxxy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyyzz[i] = 2.0 * ts_xxxy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yyzzz[i] = 2.0 * ts_xxxy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yzzzz[i] = 2.0 * ts_xxxy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_zzzzz[i] = 2.0 * ts_xxxy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxxyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 987-1008 components of targeted buffer : HH

    auto gs_z_xxxzz_xxxxx = pbuffer.data(idx_g_hh + 987);

    auto gs_z_xxxzz_xxxxy = pbuffer.data(idx_g_hh + 988);

    auto gs_z_xxxzz_xxxxz = pbuffer.data(idx_g_hh + 989);

    auto gs_z_xxxzz_xxxyy = pbuffer.data(idx_g_hh + 990);

    auto gs_z_xxxzz_xxxyz = pbuffer.data(idx_g_hh + 991);

    auto gs_z_xxxzz_xxxzz = pbuffer.data(idx_g_hh + 992);

    auto gs_z_xxxzz_xxyyy = pbuffer.data(idx_g_hh + 993);

    auto gs_z_xxxzz_xxyyz = pbuffer.data(idx_g_hh + 994);

    auto gs_z_xxxzz_xxyzz = pbuffer.data(idx_g_hh + 995);

    auto gs_z_xxxzz_xxzzz = pbuffer.data(idx_g_hh + 996);

    auto gs_z_xxxzz_xyyyy = pbuffer.data(idx_g_hh + 997);

    auto gs_z_xxxzz_xyyyz = pbuffer.data(idx_g_hh + 998);

    auto gs_z_xxxzz_xyyzz = pbuffer.data(idx_g_hh + 999);

    auto gs_z_xxxzz_xyzzz = pbuffer.data(idx_g_hh + 1000);

    auto gs_z_xxxzz_xzzzz = pbuffer.data(idx_g_hh + 1001);

    auto gs_z_xxxzz_yyyyy = pbuffer.data(idx_g_hh + 1002);

    auto gs_z_xxxzz_yyyyz = pbuffer.data(idx_g_hh + 1003);

    auto gs_z_xxxzz_yyyzz = pbuffer.data(idx_g_hh + 1004);

    auto gs_z_xxxzz_yyzzz = pbuffer.data(idx_g_hh + 1005);

    auto gs_z_xxxzz_yzzzz = pbuffer.data(idx_g_hh + 1006);

    auto gs_z_xxxzz_zzzzz = pbuffer.data(idx_g_hh + 1007);

    #pragma omp simd aligned(gc_z, gs_z_xxxzz_xxxxx, gs_z_xxxzz_xxxxy, gs_z_xxxzz_xxxxz, gs_z_xxxzz_xxxyy, gs_z_xxxzz_xxxyz, gs_z_xxxzz_xxxzz, gs_z_xxxzz_xxyyy, gs_z_xxxzz_xxyyz, gs_z_xxxzz_xxyzz, gs_z_xxxzz_xxzzz, gs_z_xxxzz_xyyyy, gs_z_xxxzz_xyyyz, gs_z_xxxzz_xyyzz, gs_z_xxxzz_xyzzz, gs_z_xxxzz_xzzzz, gs_z_xxxzz_yyyyy, gs_z_xxxzz_yyyyz, gs_z_xxxzz_yyyzz, gs_z_xxxzz_yyzzz, gs_z_xxxzz_yzzzz, gs_z_xxxzz_zzzzz, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxxz, ts_xxxz_xxxyy, ts_xxxz_xxxyz, ts_xxxz_xxxzz, ts_xxxz_xxyyy, ts_xxxz_xxyyz, ts_xxxz_xxyzz, ts_xxxz_xxzzz, ts_xxxz_xyyyy, ts_xxxz_xyyyz, ts_xxxz_xyyzz, ts_xxxz_xyzzz, ts_xxxz_xzzzz, ts_xxxz_yyyyy, ts_xxxz_yyyyz, ts_xxxz_yyyzz, ts_xxxz_yyzzz, ts_xxxz_yzzzz, ts_xxxz_zzzzz, ts_xxxzz_xxxx, ts_xxxzz_xxxxx, ts_xxxzz_xxxxy, ts_xxxzz_xxxxz, ts_xxxzz_xxxy, ts_xxxzz_xxxyy, ts_xxxzz_xxxyz, ts_xxxzz_xxxz, ts_xxxzz_xxxzz, ts_xxxzz_xxyy, ts_xxxzz_xxyyy, ts_xxxzz_xxyyz, ts_xxxzz_xxyz, ts_xxxzz_xxyzz, ts_xxxzz_xxzz, ts_xxxzz_xxzzz, ts_xxxzz_xyyy, ts_xxxzz_xyyyy, ts_xxxzz_xyyyz, ts_xxxzz_xyyz, ts_xxxzz_xyyzz, ts_xxxzz_xyzz, ts_xxxzz_xyzzz, ts_xxxzz_xzzz, ts_xxxzz_xzzzz, ts_xxxzz_yyyy, ts_xxxzz_yyyyy, ts_xxxzz_yyyyz, ts_xxxzz_yyyz, ts_xxxzz_yyyzz, ts_xxxzz_yyzz, ts_xxxzz_yyzzz, ts_xxxzz_yzzz, ts_xxxzz_yzzzz, ts_xxxzz_zzzz, ts_xxxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxzz_xxxxx[i] = 4.0 * ts_xxxz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxxxy[i] = 4.0 * ts_xxxz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxxxz[i] = 4.0 * ts_xxxz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxxyy[i] = 4.0 * ts_xxxz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxxyz[i] = 4.0 * ts_xxxz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxxzz[i] = 4.0 * ts_xxxz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxyyy[i] = 4.0 * ts_xxxz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxyyz[i] = 4.0 * ts_xxxz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxyzz[i] = 4.0 * ts_xxxz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xxzzz[i] = 4.0 * ts_xxxz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyyyy[i] = 4.0 * ts_xxxz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyyyz[i] = 4.0 * ts_xxxz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyyzz[i] = 4.0 * ts_xxxz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xyzzz[i] = 4.0 * ts_xxxz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xzzzz[i] = 4.0 * ts_xxxz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyyyy[i] = 4.0 * ts_xxxz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyyyz[i] = 4.0 * ts_xxxz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyyzz[i] = 4.0 * ts_xxxz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yyzzz[i] = 4.0 * ts_xxxz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxxzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yzzzz[i] = 4.0 * ts_xxxz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxxzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_zzzzz[i] = 4.0 * ts_xxxz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxxzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1008-1029 components of targeted buffer : HH

    auto gs_z_xxyyy_xxxxx = pbuffer.data(idx_g_hh + 1008);

    auto gs_z_xxyyy_xxxxy = pbuffer.data(idx_g_hh + 1009);

    auto gs_z_xxyyy_xxxxz = pbuffer.data(idx_g_hh + 1010);

    auto gs_z_xxyyy_xxxyy = pbuffer.data(idx_g_hh + 1011);

    auto gs_z_xxyyy_xxxyz = pbuffer.data(idx_g_hh + 1012);

    auto gs_z_xxyyy_xxxzz = pbuffer.data(idx_g_hh + 1013);

    auto gs_z_xxyyy_xxyyy = pbuffer.data(idx_g_hh + 1014);

    auto gs_z_xxyyy_xxyyz = pbuffer.data(idx_g_hh + 1015);

    auto gs_z_xxyyy_xxyzz = pbuffer.data(idx_g_hh + 1016);

    auto gs_z_xxyyy_xxzzz = pbuffer.data(idx_g_hh + 1017);

    auto gs_z_xxyyy_xyyyy = pbuffer.data(idx_g_hh + 1018);

    auto gs_z_xxyyy_xyyyz = pbuffer.data(idx_g_hh + 1019);

    auto gs_z_xxyyy_xyyzz = pbuffer.data(idx_g_hh + 1020);

    auto gs_z_xxyyy_xyzzz = pbuffer.data(idx_g_hh + 1021);

    auto gs_z_xxyyy_xzzzz = pbuffer.data(idx_g_hh + 1022);

    auto gs_z_xxyyy_yyyyy = pbuffer.data(idx_g_hh + 1023);

    auto gs_z_xxyyy_yyyyz = pbuffer.data(idx_g_hh + 1024);

    auto gs_z_xxyyy_yyyzz = pbuffer.data(idx_g_hh + 1025);

    auto gs_z_xxyyy_yyzzz = pbuffer.data(idx_g_hh + 1026);

    auto gs_z_xxyyy_yzzzz = pbuffer.data(idx_g_hh + 1027);

    auto gs_z_xxyyy_zzzzz = pbuffer.data(idx_g_hh + 1028);

    #pragma omp simd aligned(gc_z, gs_z_xxyyy_xxxxx, gs_z_xxyyy_xxxxy, gs_z_xxyyy_xxxxz, gs_z_xxyyy_xxxyy, gs_z_xxyyy_xxxyz, gs_z_xxyyy_xxxzz, gs_z_xxyyy_xxyyy, gs_z_xxyyy_xxyyz, gs_z_xxyyy_xxyzz, gs_z_xxyyy_xxzzz, gs_z_xxyyy_xyyyy, gs_z_xxyyy_xyyyz, gs_z_xxyyy_xyyzz, gs_z_xxyyy_xyzzz, gs_z_xxyyy_xzzzz, gs_z_xxyyy_yyyyy, gs_z_xxyyy_yyyyz, gs_z_xxyyy_yyyzz, gs_z_xxyyy_yyzzz, gs_z_xxyyy_yzzzz, gs_z_xxyyy_zzzzz, ts_xxyyy_xxxx, ts_xxyyy_xxxxx, ts_xxyyy_xxxxy, ts_xxyyy_xxxxz, ts_xxyyy_xxxy, ts_xxyyy_xxxyy, ts_xxyyy_xxxyz, ts_xxyyy_xxxz, ts_xxyyy_xxxzz, ts_xxyyy_xxyy, ts_xxyyy_xxyyy, ts_xxyyy_xxyyz, ts_xxyyy_xxyz, ts_xxyyy_xxyzz, ts_xxyyy_xxzz, ts_xxyyy_xxzzz, ts_xxyyy_xyyy, ts_xxyyy_xyyyy, ts_xxyyy_xyyyz, ts_xxyyy_xyyz, ts_xxyyy_xyyzz, ts_xxyyy_xyzz, ts_xxyyy_xyzzz, ts_xxyyy_xzzz, ts_xxyyy_xzzzz, ts_xxyyy_yyyy, ts_xxyyy_yyyyy, ts_xxyyy_yyyyz, ts_xxyyy_yyyz, ts_xxyyy_yyyzz, ts_xxyyy_yyzz, ts_xxyyy_yyzzz, ts_xxyyy_yzzz, ts_xxyyy_yzzzz, ts_xxyyy_zzzz, ts_xxyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyy_xxxxx[i] = 2.0 * ts_xxyyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxxxy[i] = 2.0 * ts_xxyyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxxxz[i] = 2.0 * ts_xxyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxxyy[i] = 2.0 * ts_xxyyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxxyz[i] = 2.0 * ts_xxyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxxzz[i] = 4.0 * ts_xxyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxyyy[i] = 2.0 * ts_xxyyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxyyz[i] = 2.0 * ts_xxyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxyzz[i] = 4.0 * ts_xxyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xxzzz[i] = 6.0 * ts_xxyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyyyy[i] = 2.0 * ts_xxyyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyyyz[i] = 2.0 * ts_xxyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyyzz[i] = 4.0 * ts_xxyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xyzzz[i] = 6.0 * ts_xxyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xzzzz[i] = 8.0 * ts_xxyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyyyy[i] = 2.0 * ts_xxyyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyyyz[i] = 2.0 * ts_xxyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyyzz[i] = 4.0 * ts_xxyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yyzzz[i] = 6.0 * ts_xxyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yzzzz[i] = 8.0 * ts_xxyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_zzzzz[i] = 10.0 * ts_xxyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1029-1050 components of targeted buffer : HH

    auto gs_z_xxyyz_xxxxx = pbuffer.data(idx_g_hh + 1029);

    auto gs_z_xxyyz_xxxxy = pbuffer.data(idx_g_hh + 1030);

    auto gs_z_xxyyz_xxxxz = pbuffer.data(idx_g_hh + 1031);

    auto gs_z_xxyyz_xxxyy = pbuffer.data(idx_g_hh + 1032);

    auto gs_z_xxyyz_xxxyz = pbuffer.data(idx_g_hh + 1033);

    auto gs_z_xxyyz_xxxzz = pbuffer.data(idx_g_hh + 1034);

    auto gs_z_xxyyz_xxyyy = pbuffer.data(idx_g_hh + 1035);

    auto gs_z_xxyyz_xxyyz = pbuffer.data(idx_g_hh + 1036);

    auto gs_z_xxyyz_xxyzz = pbuffer.data(idx_g_hh + 1037);

    auto gs_z_xxyyz_xxzzz = pbuffer.data(idx_g_hh + 1038);

    auto gs_z_xxyyz_xyyyy = pbuffer.data(idx_g_hh + 1039);

    auto gs_z_xxyyz_xyyyz = pbuffer.data(idx_g_hh + 1040);

    auto gs_z_xxyyz_xyyzz = pbuffer.data(idx_g_hh + 1041);

    auto gs_z_xxyyz_xyzzz = pbuffer.data(idx_g_hh + 1042);

    auto gs_z_xxyyz_xzzzz = pbuffer.data(idx_g_hh + 1043);

    auto gs_z_xxyyz_yyyyy = pbuffer.data(idx_g_hh + 1044);

    auto gs_z_xxyyz_yyyyz = pbuffer.data(idx_g_hh + 1045);

    auto gs_z_xxyyz_yyyzz = pbuffer.data(idx_g_hh + 1046);

    auto gs_z_xxyyz_yyzzz = pbuffer.data(idx_g_hh + 1047);

    auto gs_z_xxyyz_yzzzz = pbuffer.data(idx_g_hh + 1048);

    auto gs_z_xxyyz_zzzzz = pbuffer.data(idx_g_hh + 1049);

    #pragma omp simd aligned(gc_z, gs_z_xxyyz_xxxxx, gs_z_xxyyz_xxxxy, gs_z_xxyyz_xxxxz, gs_z_xxyyz_xxxyy, gs_z_xxyyz_xxxyz, gs_z_xxyyz_xxxzz, gs_z_xxyyz_xxyyy, gs_z_xxyyz_xxyyz, gs_z_xxyyz_xxyzz, gs_z_xxyyz_xxzzz, gs_z_xxyyz_xyyyy, gs_z_xxyyz_xyyyz, gs_z_xxyyz_xyyzz, gs_z_xxyyz_xyzzz, gs_z_xxyyz_xzzzz, gs_z_xxyyz_yyyyy, gs_z_xxyyz_yyyyz, gs_z_xxyyz_yyyzz, gs_z_xxyyz_yyzzz, gs_z_xxyyz_yzzzz, gs_z_xxyyz_zzzzz, ts_xxyy_xxxxx, ts_xxyy_xxxxy, ts_xxyy_xxxxz, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxxzz, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyzz, ts_xxyy_xxzzz, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyzz, ts_xxyy_xyzzz, ts_xxyy_xzzzz, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyzz, ts_xxyy_yyzzz, ts_xxyy_yzzzz, ts_xxyy_zzzzz, ts_xxyyz_xxxx, ts_xxyyz_xxxxx, ts_xxyyz_xxxxy, ts_xxyyz_xxxxz, ts_xxyyz_xxxy, ts_xxyyz_xxxyy, ts_xxyyz_xxxyz, ts_xxyyz_xxxz, ts_xxyyz_xxxzz, ts_xxyyz_xxyy, ts_xxyyz_xxyyy, ts_xxyyz_xxyyz, ts_xxyyz_xxyz, ts_xxyyz_xxyzz, ts_xxyyz_xxzz, ts_xxyyz_xxzzz, ts_xxyyz_xyyy, ts_xxyyz_xyyyy, ts_xxyyz_xyyyz, ts_xxyyz_xyyz, ts_xxyyz_xyyzz, ts_xxyyz_xyzz, ts_xxyyz_xyzzz, ts_xxyyz_xzzz, ts_xxyyz_xzzzz, ts_xxyyz_yyyy, ts_xxyyz_yyyyy, ts_xxyyz_yyyyz, ts_xxyyz_yyyz, ts_xxyyz_yyyzz, ts_xxyyz_yyzz, ts_xxyyz_yyzzz, ts_xxyyz_yzzz, ts_xxyyz_yzzzz, ts_xxyyz_zzzz, ts_xxyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyz_xxxxx[i] = 2.0 * ts_xxyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxxxy[i] = 2.0 * ts_xxyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxxxz[i] = 2.0 * ts_xxyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxxyy[i] = 2.0 * ts_xxyy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxxyz[i] = 2.0 * ts_xxyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxxzz[i] = 2.0 * ts_xxyy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxyyy[i] = 2.0 * ts_xxyy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxyyz[i] = 2.0 * ts_xxyy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxyzz[i] = 2.0 * ts_xxyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xxzzz[i] = 2.0 * ts_xxyy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyyyy[i] = 2.0 * ts_xxyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyyyz[i] = 2.0 * ts_xxyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyyzz[i] = 2.0 * ts_xxyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xyzzz[i] = 2.0 * ts_xxyy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xzzzz[i] = 2.0 * ts_xxyy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyyyy[i] = 2.0 * ts_xxyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyyyz[i] = 2.0 * ts_xxyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyyzz[i] = 2.0 * ts_xxyy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yyzzz[i] = 2.0 * ts_xxyy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yzzzz[i] = 2.0 * ts_xxyy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_zzzzz[i] = 2.0 * ts_xxyy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1050-1071 components of targeted buffer : HH

    auto gs_z_xxyzz_xxxxx = pbuffer.data(idx_g_hh + 1050);

    auto gs_z_xxyzz_xxxxy = pbuffer.data(idx_g_hh + 1051);

    auto gs_z_xxyzz_xxxxz = pbuffer.data(idx_g_hh + 1052);

    auto gs_z_xxyzz_xxxyy = pbuffer.data(idx_g_hh + 1053);

    auto gs_z_xxyzz_xxxyz = pbuffer.data(idx_g_hh + 1054);

    auto gs_z_xxyzz_xxxzz = pbuffer.data(idx_g_hh + 1055);

    auto gs_z_xxyzz_xxyyy = pbuffer.data(idx_g_hh + 1056);

    auto gs_z_xxyzz_xxyyz = pbuffer.data(idx_g_hh + 1057);

    auto gs_z_xxyzz_xxyzz = pbuffer.data(idx_g_hh + 1058);

    auto gs_z_xxyzz_xxzzz = pbuffer.data(idx_g_hh + 1059);

    auto gs_z_xxyzz_xyyyy = pbuffer.data(idx_g_hh + 1060);

    auto gs_z_xxyzz_xyyyz = pbuffer.data(idx_g_hh + 1061);

    auto gs_z_xxyzz_xyyzz = pbuffer.data(idx_g_hh + 1062);

    auto gs_z_xxyzz_xyzzz = pbuffer.data(idx_g_hh + 1063);

    auto gs_z_xxyzz_xzzzz = pbuffer.data(idx_g_hh + 1064);

    auto gs_z_xxyzz_yyyyy = pbuffer.data(idx_g_hh + 1065);

    auto gs_z_xxyzz_yyyyz = pbuffer.data(idx_g_hh + 1066);

    auto gs_z_xxyzz_yyyzz = pbuffer.data(idx_g_hh + 1067);

    auto gs_z_xxyzz_yyzzz = pbuffer.data(idx_g_hh + 1068);

    auto gs_z_xxyzz_yzzzz = pbuffer.data(idx_g_hh + 1069);

    auto gs_z_xxyzz_zzzzz = pbuffer.data(idx_g_hh + 1070);

    #pragma omp simd aligned(gc_z, gs_z_xxyzz_xxxxx, gs_z_xxyzz_xxxxy, gs_z_xxyzz_xxxxz, gs_z_xxyzz_xxxyy, gs_z_xxyzz_xxxyz, gs_z_xxyzz_xxxzz, gs_z_xxyzz_xxyyy, gs_z_xxyzz_xxyyz, gs_z_xxyzz_xxyzz, gs_z_xxyzz_xxzzz, gs_z_xxyzz_xyyyy, gs_z_xxyzz_xyyyz, gs_z_xxyzz_xyyzz, gs_z_xxyzz_xyzzz, gs_z_xxyzz_xzzzz, gs_z_xxyzz_yyyyy, gs_z_xxyzz_yyyyz, gs_z_xxyzz_yyyzz, gs_z_xxyzz_yyzzz, gs_z_xxyzz_yzzzz, gs_z_xxyzz_zzzzz, ts_xxyz_xxxxx, ts_xxyz_xxxxy, ts_xxyz_xxxxz, ts_xxyz_xxxyy, ts_xxyz_xxxyz, ts_xxyz_xxxzz, ts_xxyz_xxyyy, ts_xxyz_xxyyz, ts_xxyz_xxyzz, ts_xxyz_xxzzz, ts_xxyz_xyyyy, ts_xxyz_xyyyz, ts_xxyz_xyyzz, ts_xxyz_xyzzz, ts_xxyz_xzzzz, ts_xxyz_yyyyy, ts_xxyz_yyyyz, ts_xxyz_yyyzz, ts_xxyz_yyzzz, ts_xxyz_yzzzz, ts_xxyz_zzzzz, ts_xxyzz_xxxx, ts_xxyzz_xxxxx, ts_xxyzz_xxxxy, ts_xxyzz_xxxxz, ts_xxyzz_xxxy, ts_xxyzz_xxxyy, ts_xxyzz_xxxyz, ts_xxyzz_xxxz, ts_xxyzz_xxxzz, ts_xxyzz_xxyy, ts_xxyzz_xxyyy, ts_xxyzz_xxyyz, ts_xxyzz_xxyz, ts_xxyzz_xxyzz, ts_xxyzz_xxzz, ts_xxyzz_xxzzz, ts_xxyzz_xyyy, ts_xxyzz_xyyyy, ts_xxyzz_xyyyz, ts_xxyzz_xyyz, ts_xxyzz_xyyzz, ts_xxyzz_xyzz, ts_xxyzz_xyzzz, ts_xxyzz_xzzz, ts_xxyzz_xzzzz, ts_xxyzz_yyyy, ts_xxyzz_yyyyy, ts_xxyzz_yyyyz, ts_xxyzz_yyyz, ts_xxyzz_yyyzz, ts_xxyzz_yyzz, ts_xxyzz_yyzzz, ts_xxyzz_yzzz, ts_xxyzz_yzzzz, ts_xxyzz_zzzz, ts_xxyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyzz_xxxxx[i] = 4.0 * ts_xxyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxxxy[i] = 4.0 * ts_xxyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxxxz[i] = 4.0 * ts_xxyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxxyy[i] = 4.0 * ts_xxyz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxxyz[i] = 4.0 * ts_xxyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxxzz[i] = 4.0 * ts_xxyz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxyyy[i] = 4.0 * ts_xxyz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxyyz[i] = 4.0 * ts_xxyz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxyzz[i] = 4.0 * ts_xxyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xxzzz[i] = 4.0 * ts_xxyz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyyyy[i] = 4.0 * ts_xxyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyyyz[i] = 4.0 * ts_xxyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyyzz[i] = 4.0 * ts_xxyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xyzzz[i] = 4.0 * ts_xxyz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xzzzz[i] = 4.0 * ts_xxyz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyyyy[i] = 4.0 * ts_xxyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyyyz[i] = 4.0 * ts_xxyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyyzz[i] = 4.0 * ts_xxyz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yyzzz[i] = 4.0 * ts_xxyz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yzzzz[i] = 4.0 * ts_xxyz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_zzzzz[i] = 4.0 * ts_xxyz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1071-1092 components of targeted buffer : HH

    auto gs_z_xxzzz_xxxxx = pbuffer.data(idx_g_hh + 1071);

    auto gs_z_xxzzz_xxxxy = pbuffer.data(idx_g_hh + 1072);

    auto gs_z_xxzzz_xxxxz = pbuffer.data(idx_g_hh + 1073);

    auto gs_z_xxzzz_xxxyy = pbuffer.data(idx_g_hh + 1074);

    auto gs_z_xxzzz_xxxyz = pbuffer.data(idx_g_hh + 1075);

    auto gs_z_xxzzz_xxxzz = pbuffer.data(idx_g_hh + 1076);

    auto gs_z_xxzzz_xxyyy = pbuffer.data(idx_g_hh + 1077);

    auto gs_z_xxzzz_xxyyz = pbuffer.data(idx_g_hh + 1078);

    auto gs_z_xxzzz_xxyzz = pbuffer.data(idx_g_hh + 1079);

    auto gs_z_xxzzz_xxzzz = pbuffer.data(idx_g_hh + 1080);

    auto gs_z_xxzzz_xyyyy = pbuffer.data(idx_g_hh + 1081);

    auto gs_z_xxzzz_xyyyz = pbuffer.data(idx_g_hh + 1082);

    auto gs_z_xxzzz_xyyzz = pbuffer.data(idx_g_hh + 1083);

    auto gs_z_xxzzz_xyzzz = pbuffer.data(idx_g_hh + 1084);

    auto gs_z_xxzzz_xzzzz = pbuffer.data(idx_g_hh + 1085);

    auto gs_z_xxzzz_yyyyy = pbuffer.data(idx_g_hh + 1086);

    auto gs_z_xxzzz_yyyyz = pbuffer.data(idx_g_hh + 1087);

    auto gs_z_xxzzz_yyyzz = pbuffer.data(idx_g_hh + 1088);

    auto gs_z_xxzzz_yyzzz = pbuffer.data(idx_g_hh + 1089);

    auto gs_z_xxzzz_yzzzz = pbuffer.data(idx_g_hh + 1090);

    auto gs_z_xxzzz_zzzzz = pbuffer.data(idx_g_hh + 1091);

    #pragma omp simd aligned(gc_z, gs_z_xxzzz_xxxxx, gs_z_xxzzz_xxxxy, gs_z_xxzzz_xxxxz, gs_z_xxzzz_xxxyy, gs_z_xxzzz_xxxyz, gs_z_xxzzz_xxxzz, gs_z_xxzzz_xxyyy, gs_z_xxzzz_xxyyz, gs_z_xxzzz_xxyzz, gs_z_xxzzz_xxzzz, gs_z_xxzzz_xyyyy, gs_z_xxzzz_xyyyz, gs_z_xxzzz_xyyzz, gs_z_xxzzz_xyzzz, gs_z_xxzzz_xzzzz, gs_z_xxzzz_yyyyy, gs_z_xxzzz_yyyyz, gs_z_xxzzz_yyyzz, gs_z_xxzzz_yyzzz, gs_z_xxzzz_yzzzz, gs_z_xxzzz_zzzzz, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxxz, ts_xxzz_xxxyy, ts_xxzz_xxxyz, ts_xxzz_xxxzz, ts_xxzz_xxyyy, ts_xxzz_xxyyz, ts_xxzz_xxyzz, ts_xxzz_xxzzz, ts_xxzz_xyyyy, ts_xxzz_xyyyz, ts_xxzz_xyyzz, ts_xxzz_xyzzz, ts_xxzz_xzzzz, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyzz, ts_xxzz_yyzzz, ts_xxzz_yzzzz, ts_xxzz_zzzzz, ts_xxzzz_xxxx, ts_xxzzz_xxxxx, ts_xxzzz_xxxxy, ts_xxzzz_xxxxz, ts_xxzzz_xxxy, ts_xxzzz_xxxyy, ts_xxzzz_xxxyz, ts_xxzzz_xxxz, ts_xxzzz_xxxzz, ts_xxzzz_xxyy, ts_xxzzz_xxyyy, ts_xxzzz_xxyyz, ts_xxzzz_xxyz, ts_xxzzz_xxyzz, ts_xxzzz_xxzz, ts_xxzzz_xxzzz, ts_xxzzz_xyyy, ts_xxzzz_xyyyy, ts_xxzzz_xyyyz, ts_xxzzz_xyyz, ts_xxzzz_xyyzz, ts_xxzzz_xyzz, ts_xxzzz_xyzzz, ts_xxzzz_xzzz, ts_xxzzz_xzzzz, ts_xxzzz_yyyy, ts_xxzzz_yyyyy, ts_xxzzz_yyyyz, ts_xxzzz_yyyz, ts_xxzzz_yyyzz, ts_xxzzz_yyzz, ts_xxzzz_yyzzz, ts_xxzzz_yzzz, ts_xxzzz_yzzzz, ts_xxzzz_zzzz, ts_xxzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzzz_xxxxx[i] = 6.0 * ts_xxzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxxxy[i] = 6.0 * ts_xxzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxxxz[i] = 6.0 * ts_xxzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxxyy[i] = 6.0 * ts_xxzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxxyz[i] = 6.0 * ts_xxzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxxzz[i] = 6.0 * ts_xxzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxyyy[i] = 6.0 * ts_xxzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxyyz[i] = 6.0 * ts_xxzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxyzz[i] = 6.0 * ts_xxzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xxzzz[i] = 6.0 * ts_xxzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyyyy[i] = 6.0 * ts_xxzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyyyz[i] = 6.0 * ts_xxzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyyzz[i] = 6.0 * ts_xxzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xyzzz[i] = 6.0 * ts_xxzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xzzzz[i] = 6.0 * ts_xxzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyyyy[i] = 6.0 * ts_xxzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyyyz[i] = 6.0 * ts_xxzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyyzz[i] = 6.0 * ts_xxzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yyzzz[i] = 6.0 * ts_xxzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yzzzz[i] = 6.0 * ts_xxzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_zzzzz[i] = 6.0 * ts_xxzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1092-1113 components of targeted buffer : HH

    auto gs_z_xyyyy_xxxxx = pbuffer.data(idx_g_hh + 1092);

    auto gs_z_xyyyy_xxxxy = pbuffer.data(idx_g_hh + 1093);

    auto gs_z_xyyyy_xxxxz = pbuffer.data(idx_g_hh + 1094);

    auto gs_z_xyyyy_xxxyy = pbuffer.data(idx_g_hh + 1095);

    auto gs_z_xyyyy_xxxyz = pbuffer.data(idx_g_hh + 1096);

    auto gs_z_xyyyy_xxxzz = pbuffer.data(idx_g_hh + 1097);

    auto gs_z_xyyyy_xxyyy = pbuffer.data(idx_g_hh + 1098);

    auto gs_z_xyyyy_xxyyz = pbuffer.data(idx_g_hh + 1099);

    auto gs_z_xyyyy_xxyzz = pbuffer.data(idx_g_hh + 1100);

    auto gs_z_xyyyy_xxzzz = pbuffer.data(idx_g_hh + 1101);

    auto gs_z_xyyyy_xyyyy = pbuffer.data(idx_g_hh + 1102);

    auto gs_z_xyyyy_xyyyz = pbuffer.data(idx_g_hh + 1103);

    auto gs_z_xyyyy_xyyzz = pbuffer.data(idx_g_hh + 1104);

    auto gs_z_xyyyy_xyzzz = pbuffer.data(idx_g_hh + 1105);

    auto gs_z_xyyyy_xzzzz = pbuffer.data(idx_g_hh + 1106);

    auto gs_z_xyyyy_yyyyy = pbuffer.data(idx_g_hh + 1107);

    auto gs_z_xyyyy_yyyyz = pbuffer.data(idx_g_hh + 1108);

    auto gs_z_xyyyy_yyyzz = pbuffer.data(idx_g_hh + 1109);

    auto gs_z_xyyyy_yyzzz = pbuffer.data(idx_g_hh + 1110);

    auto gs_z_xyyyy_yzzzz = pbuffer.data(idx_g_hh + 1111);

    auto gs_z_xyyyy_zzzzz = pbuffer.data(idx_g_hh + 1112);

    #pragma omp simd aligned(gc_z, gs_z_xyyyy_xxxxx, gs_z_xyyyy_xxxxy, gs_z_xyyyy_xxxxz, gs_z_xyyyy_xxxyy, gs_z_xyyyy_xxxyz, gs_z_xyyyy_xxxzz, gs_z_xyyyy_xxyyy, gs_z_xyyyy_xxyyz, gs_z_xyyyy_xxyzz, gs_z_xyyyy_xxzzz, gs_z_xyyyy_xyyyy, gs_z_xyyyy_xyyyz, gs_z_xyyyy_xyyzz, gs_z_xyyyy_xyzzz, gs_z_xyyyy_xzzzz, gs_z_xyyyy_yyyyy, gs_z_xyyyy_yyyyz, gs_z_xyyyy_yyyzz, gs_z_xyyyy_yyzzz, gs_z_xyyyy_yzzzz, gs_z_xyyyy_zzzzz, ts_xyyyy_xxxx, ts_xyyyy_xxxxx, ts_xyyyy_xxxxy, ts_xyyyy_xxxxz, ts_xyyyy_xxxy, ts_xyyyy_xxxyy, ts_xyyyy_xxxyz, ts_xyyyy_xxxz, ts_xyyyy_xxxzz, ts_xyyyy_xxyy, ts_xyyyy_xxyyy, ts_xyyyy_xxyyz, ts_xyyyy_xxyz, ts_xyyyy_xxyzz, ts_xyyyy_xxzz, ts_xyyyy_xxzzz, ts_xyyyy_xyyy, ts_xyyyy_xyyyy, ts_xyyyy_xyyyz, ts_xyyyy_xyyz, ts_xyyyy_xyyzz, ts_xyyyy_xyzz, ts_xyyyy_xyzzz, ts_xyyyy_xzzz, ts_xyyyy_xzzzz, ts_xyyyy_yyyy, ts_xyyyy_yyyyy, ts_xyyyy_yyyyz, ts_xyyyy_yyyz, ts_xyyyy_yyyzz, ts_xyyyy_yyzz, ts_xyyyy_yyzzz, ts_xyyyy_yzzz, ts_xyyyy_yzzzz, ts_xyyyy_zzzz, ts_xyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyy_xxxxx[i] = 2.0 * ts_xyyyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxxxy[i] = 2.0 * ts_xyyyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxxxz[i] = 2.0 * ts_xyyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxxyy[i] = 2.0 * ts_xyyyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxxyz[i] = 2.0 * ts_xyyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxxzz[i] = 4.0 * ts_xyyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxyyy[i] = 2.0 * ts_xyyyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxyyz[i] = 2.0 * ts_xyyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxyzz[i] = 4.0 * ts_xyyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xxzzz[i] = 6.0 * ts_xyyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyyyy[i] = 2.0 * ts_xyyyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyyyz[i] = 2.0 * ts_xyyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyyzz[i] = 4.0 * ts_xyyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xyzzz[i] = 6.0 * ts_xyyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xzzzz[i] = 8.0 * ts_xyyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyyyy[i] = 2.0 * ts_xyyyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyyyz[i] = 2.0 * ts_xyyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyyzz[i] = 4.0 * ts_xyyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yyzzz[i] = 6.0 * ts_xyyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yzzzz[i] = 8.0 * ts_xyyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_zzzzz[i] = 10.0 * ts_xyyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1113-1134 components of targeted buffer : HH

    auto gs_z_xyyyz_xxxxx = pbuffer.data(idx_g_hh + 1113);

    auto gs_z_xyyyz_xxxxy = pbuffer.data(idx_g_hh + 1114);

    auto gs_z_xyyyz_xxxxz = pbuffer.data(idx_g_hh + 1115);

    auto gs_z_xyyyz_xxxyy = pbuffer.data(idx_g_hh + 1116);

    auto gs_z_xyyyz_xxxyz = pbuffer.data(idx_g_hh + 1117);

    auto gs_z_xyyyz_xxxzz = pbuffer.data(idx_g_hh + 1118);

    auto gs_z_xyyyz_xxyyy = pbuffer.data(idx_g_hh + 1119);

    auto gs_z_xyyyz_xxyyz = pbuffer.data(idx_g_hh + 1120);

    auto gs_z_xyyyz_xxyzz = pbuffer.data(idx_g_hh + 1121);

    auto gs_z_xyyyz_xxzzz = pbuffer.data(idx_g_hh + 1122);

    auto gs_z_xyyyz_xyyyy = pbuffer.data(idx_g_hh + 1123);

    auto gs_z_xyyyz_xyyyz = pbuffer.data(idx_g_hh + 1124);

    auto gs_z_xyyyz_xyyzz = pbuffer.data(idx_g_hh + 1125);

    auto gs_z_xyyyz_xyzzz = pbuffer.data(idx_g_hh + 1126);

    auto gs_z_xyyyz_xzzzz = pbuffer.data(idx_g_hh + 1127);

    auto gs_z_xyyyz_yyyyy = pbuffer.data(idx_g_hh + 1128);

    auto gs_z_xyyyz_yyyyz = pbuffer.data(idx_g_hh + 1129);

    auto gs_z_xyyyz_yyyzz = pbuffer.data(idx_g_hh + 1130);

    auto gs_z_xyyyz_yyzzz = pbuffer.data(idx_g_hh + 1131);

    auto gs_z_xyyyz_yzzzz = pbuffer.data(idx_g_hh + 1132);

    auto gs_z_xyyyz_zzzzz = pbuffer.data(idx_g_hh + 1133);

    #pragma omp simd aligned(gc_z, gs_z_xyyyz_xxxxx, gs_z_xyyyz_xxxxy, gs_z_xyyyz_xxxxz, gs_z_xyyyz_xxxyy, gs_z_xyyyz_xxxyz, gs_z_xyyyz_xxxzz, gs_z_xyyyz_xxyyy, gs_z_xyyyz_xxyyz, gs_z_xyyyz_xxyzz, gs_z_xyyyz_xxzzz, gs_z_xyyyz_xyyyy, gs_z_xyyyz_xyyyz, gs_z_xyyyz_xyyzz, gs_z_xyyyz_xyzzz, gs_z_xyyyz_xzzzz, gs_z_xyyyz_yyyyy, gs_z_xyyyz_yyyyz, gs_z_xyyyz_yyyzz, gs_z_xyyyz_yyzzz, gs_z_xyyyz_yzzzz, gs_z_xyyyz_zzzzz, ts_xyyy_xxxxx, ts_xyyy_xxxxy, ts_xyyy_xxxxz, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxxzz, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyzz, ts_xyyy_xxzzz, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyzz, ts_xyyy_xyzzz, ts_xyyy_xzzzz, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyzz, ts_xyyy_yyzzz, ts_xyyy_yzzzz, ts_xyyy_zzzzz, ts_xyyyz_xxxx, ts_xyyyz_xxxxx, ts_xyyyz_xxxxy, ts_xyyyz_xxxxz, ts_xyyyz_xxxy, ts_xyyyz_xxxyy, ts_xyyyz_xxxyz, ts_xyyyz_xxxz, ts_xyyyz_xxxzz, ts_xyyyz_xxyy, ts_xyyyz_xxyyy, ts_xyyyz_xxyyz, ts_xyyyz_xxyz, ts_xyyyz_xxyzz, ts_xyyyz_xxzz, ts_xyyyz_xxzzz, ts_xyyyz_xyyy, ts_xyyyz_xyyyy, ts_xyyyz_xyyyz, ts_xyyyz_xyyz, ts_xyyyz_xyyzz, ts_xyyyz_xyzz, ts_xyyyz_xyzzz, ts_xyyyz_xzzz, ts_xyyyz_xzzzz, ts_xyyyz_yyyy, ts_xyyyz_yyyyy, ts_xyyyz_yyyyz, ts_xyyyz_yyyz, ts_xyyyz_yyyzz, ts_xyyyz_yyzz, ts_xyyyz_yyzzz, ts_xyyyz_yzzz, ts_xyyyz_yzzzz, ts_xyyyz_zzzz, ts_xyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyz_xxxxx[i] = 2.0 * ts_xyyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxxxy[i] = 2.0 * ts_xyyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxxxz[i] = 2.0 * ts_xyyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxxyy[i] = 2.0 * ts_xyyy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxxyz[i] = 2.0 * ts_xyyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxxzz[i] = 2.0 * ts_xyyy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxyyy[i] = 2.0 * ts_xyyy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxyyz[i] = 2.0 * ts_xyyy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxyzz[i] = 2.0 * ts_xyyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xxzzz[i] = 2.0 * ts_xyyy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyyyy[i] = 2.0 * ts_xyyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyyyz[i] = 2.0 * ts_xyyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyyzz[i] = 2.0 * ts_xyyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xyzzz[i] = 2.0 * ts_xyyy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xzzzz[i] = 2.0 * ts_xyyy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyyyy[i] = 2.0 * ts_xyyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyyyz[i] = 2.0 * ts_xyyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyyzz[i] = 2.0 * ts_xyyy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yyzzz[i] = 2.0 * ts_xyyy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yzzzz[i] = 2.0 * ts_xyyy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_zzzzz[i] = 2.0 * ts_xyyy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xyyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1134-1155 components of targeted buffer : HH

    auto gs_z_xyyzz_xxxxx = pbuffer.data(idx_g_hh + 1134);

    auto gs_z_xyyzz_xxxxy = pbuffer.data(idx_g_hh + 1135);

    auto gs_z_xyyzz_xxxxz = pbuffer.data(idx_g_hh + 1136);

    auto gs_z_xyyzz_xxxyy = pbuffer.data(idx_g_hh + 1137);

    auto gs_z_xyyzz_xxxyz = pbuffer.data(idx_g_hh + 1138);

    auto gs_z_xyyzz_xxxzz = pbuffer.data(idx_g_hh + 1139);

    auto gs_z_xyyzz_xxyyy = pbuffer.data(idx_g_hh + 1140);

    auto gs_z_xyyzz_xxyyz = pbuffer.data(idx_g_hh + 1141);

    auto gs_z_xyyzz_xxyzz = pbuffer.data(idx_g_hh + 1142);

    auto gs_z_xyyzz_xxzzz = pbuffer.data(idx_g_hh + 1143);

    auto gs_z_xyyzz_xyyyy = pbuffer.data(idx_g_hh + 1144);

    auto gs_z_xyyzz_xyyyz = pbuffer.data(idx_g_hh + 1145);

    auto gs_z_xyyzz_xyyzz = pbuffer.data(idx_g_hh + 1146);

    auto gs_z_xyyzz_xyzzz = pbuffer.data(idx_g_hh + 1147);

    auto gs_z_xyyzz_xzzzz = pbuffer.data(idx_g_hh + 1148);

    auto gs_z_xyyzz_yyyyy = pbuffer.data(idx_g_hh + 1149);

    auto gs_z_xyyzz_yyyyz = pbuffer.data(idx_g_hh + 1150);

    auto gs_z_xyyzz_yyyzz = pbuffer.data(idx_g_hh + 1151);

    auto gs_z_xyyzz_yyzzz = pbuffer.data(idx_g_hh + 1152);

    auto gs_z_xyyzz_yzzzz = pbuffer.data(idx_g_hh + 1153);

    auto gs_z_xyyzz_zzzzz = pbuffer.data(idx_g_hh + 1154);

    #pragma omp simd aligned(gc_z, gs_z_xyyzz_xxxxx, gs_z_xyyzz_xxxxy, gs_z_xyyzz_xxxxz, gs_z_xyyzz_xxxyy, gs_z_xyyzz_xxxyz, gs_z_xyyzz_xxxzz, gs_z_xyyzz_xxyyy, gs_z_xyyzz_xxyyz, gs_z_xyyzz_xxyzz, gs_z_xyyzz_xxzzz, gs_z_xyyzz_xyyyy, gs_z_xyyzz_xyyyz, gs_z_xyyzz_xyyzz, gs_z_xyyzz_xyzzz, gs_z_xyyzz_xzzzz, gs_z_xyyzz_yyyyy, gs_z_xyyzz_yyyyz, gs_z_xyyzz_yyyzz, gs_z_xyyzz_yyzzz, gs_z_xyyzz_yzzzz, gs_z_xyyzz_zzzzz, ts_xyyz_xxxxx, ts_xyyz_xxxxy, ts_xyyz_xxxxz, ts_xyyz_xxxyy, ts_xyyz_xxxyz, ts_xyyz_xxxzz, ts_xyyz_xxyyy, ts_xyyz_xxyyz, ts_xyyz_xxyzz, ts_xyyz_xxzzz, ts_xyyz_xyyyy, ts_xyyz_xyyyz, ts_xyyz_xyyzz, ts_xyyz_xyzzz, ts_xyyz_xzzzz, ts_xyyz_yyyyy, ts_xyyz_yyyyz, ts_xyyz_yyyzz, ts_xyyz_yyzzz, ts_xyyz_yzzzz, ts_xyyz_zzzzz, ts_xyyzz_xxxx, ts_xyyzz_xxxxx, ts_xyyzz_xxxxy, ts_xyyzz_xxxxz, ts_xyyzz_xxxy, ts_xyyzz_xxxyy, ts_xyyzz_xxxyz, ts_xyyzz_xxxz, ts_xyyzz_xxxzz, ts_xyyzz_xxyy, ts_xyyzz_xxyyy, ts_xyyzz_xxyyz, ts_xyyzz_xxyz, ts_xyyzz_xxyzz, ts_xyyzz_xxzz, ts_xyyzz_xxzzz, ts_xyyzz_xyyy, ts_xyyzz_xyyyy, ts_xyyzz_xyyyz, ts_xyyzz_xyyz, ts_xyyzz_xyyzz, ts_xyyzz_xyzz, ts_xyyzz_xyzzz, ts_xyyzz_xzzz, ts_xyyzz_xzzzz, ts_xyyzz_yyyy, ts_xyyzz_yyyyy, ts_xyyzz_yyyyz, ts_xyyzz_yyyz, ts_xyyzz_yyyzz, ts_xyyzz_yyzz, ts_xyyzz_yyzzz, ts_xyyzz_yzzz, ts_xyyzz_yzzzz, ts_xyyzz_zzzz, ts_xyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyzz_xxxxx[i] = 4.0 * ts_xyyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxxxy[i] = 4.0 * ts_xyyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxxxz[i] = 4.0 * ts_xyyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxxyy[i] = 4.0 * ts_xyyz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxxyz[i] = 4.0 * ts_xyyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxxzz[i] = 4.0 * ts_xyyz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxyyy[i] = 4.0 * ts_xyyz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxyyz[i] = 4.0 * ts_xyyz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxyzz[i] = 4.0 * ts_xyyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xxzzz[i] = 4.0 * ts_xyyz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyyyy[i] = 4.0 * ts_xyyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyyyz[i] = 4.0 * ts_xyyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyyzz[i] = 4.0 * ts_xyyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xyzzz[i] = 4.0 * ts_xyyz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xzzzz[i] = 4.0 * ts_xyyz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyyyy[i] = 4.0 * ts_xyyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyyyz[i] = 4.0 * ts_xyyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyyzz[i] = 4.0 * ts_xyyz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yyzzz[i] = 4.0 * ts_xyyz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yzzzz[i] = 4.0 * ts_xyyz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_zzzzz[i] = 4.0 * ts_xyyz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xyyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1155-1176 components of targeted buffer : HH

    auto gs_z_xyzzz_xxxxx = pbuffer.data(idx_g_hh + 1155);

    auto gs_z_xyzzz_xxxxy = pbuffer.data(idx_g_hh + 1156);

    auto gs_z_xyzzz_xxxxz = pbuffer.data(idx_g_hh + 1157);

    auto gs_z_xyzzz_xxxyy = pbuffer.data(idx_g_hh + 1158);

    auto gs_z_xyzzz_xxxyz = pbuffer.data(idx_g_hh + 1159);

    auto gs_z_xyzzz_xxxzz = pbuffer.data(idx_g_hh + 1160);

    auto gs_z_xyzzz_xxyyy = pbuffer.data(idx_g_hh + 1161);

    auto gs_z_xyzzz_xxyyz = pbuffer.data(idx_g_hh + 1162);

    auto gs_z_xyzzz_xxyzz = pbuffer.data(idx_g_hh + 1163);

    auto gs_z_xyzzz_xxzzz = pbuffer.data(idx_g_hh + 1164);

    auto gs_z_xyzzz_xyyyy = pbuffer.data(idx_g_hh + 1165);

    auto gs_z_xyzzz_xyyyz = pbuffer.data(idx_g_hh + 1166);

    auto gs_z_xyzzz_xyyzz = pbuffer.data(idx_g_hh + 1167);

    auto gs_z_xyzzz_xyzzz = pbuffer.data(idx_g_hh + 1168);

    auto gs_z_xyzzz_xzzzz = pbuffer.data(idx_g_hh + 1169);

    auto gs_z_xyzzz_yyyyy = pbuffer.data(idx_g_hh + 1170);

    auto gs_z_xyzzz_yyyyz = pbuffer.data(idx_g_hh + 1171);

    auto gs_z_xyzzz_yyyzz = pbuffer.data(idx_g_hh + 1172);

    auto gs_z_xyzzz_yyzzz = pbuffer.data(idx_g_hh + 1173);

    auto gs_z_xyzzz_yzzzz = pbuffer.data(idx_g_hh + 1174);

    auto gs_z_xyzzz_zzzzz = pbuffer.data(idx_g_hh + 1175);

    #pragma omp simd aligned(gc_z, gs_z_xyzzz_xxxxx, gs_z_xyzzz_xxxxy, gs_z_xyzzz_xxxxz, gs_z_xyzzz_xxxyy, gs_z_xyzzz_xxxyz, gs_z_xyzzz_xxxzz, gs_z_xyzzz_xxyyy, gs_z_xyzzz_xxyyz, gs_z_xyzzz_xxyzz, gs_z_xyzzz_xxzzz, gs_z_xyzzz_xyyyy, gs_z_xyzzz_xyyyz, gs_z_xyzzz_xyyzz, gs_z_xyzzz_xyzzz, gs_z_xyzzz_xzzzz, gs_z_xyzzz_yyyyy, gs_z_xyzzz_yyyyz, gs_z_xyzzz_yyyzz, gs_z_xyzzz_yyzzz, gs_z_xyzzz_yzzzz, gs_z_xyzzz_zzzzz, ts_xyzz_xxxxx, ts_xyzz_xxxxy, ts_xyzz_xxxxz, ts_xyzz_xxxyy, ts_xyzz_xxxyz, ts_xyzz_xxxzz, ts_xyzz_xxyyy, ts_xyzz_xxyyz, ts_xyzz_xxyzz, ts_xyzz_xxzzz, ts_xyzz_xyyyy, ts_xyzz_xyyyz, ts_xyzz_xyyzz, ts_xyzz_xyzzz, ts_xyzz_xzzzz, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyzz, ts_xyzz_yyzzz, ts_xyzz_yzzzz, ts_xyzz_zzzzz, ts_xyzzz_xxxx, ts_xyzzz_xxxxx, ts_xyzzz_xxxxy, ts_xyzzz_xxxxz, ts_xyzzz_xxxy, ts_xyzzz_xxxyy, ts_xyzzz_xxxyz, ts_xyzzz_xxxz, ts_xyzzz_xxxzz, ts_xyzzz_xxyy, ts_xyzzz_xxyyy, ts_xyzzz_xxyyz, ts_xyzzz_xxyz, ts_xyzzz_xxyzz, ts_xyzzz_xxzz, ts_xyzzz_xxzzz, ts_xyzzz_xyyy, ts_xyzzz_xyyyy, ts_xyzzz_xyyyz, ts_xyzzz_xyyz, ts_xyzzz_xyyzz, ts_xyzzz_xyzz, ts_xyzzz_xyzzz, ts_xyzzz_xzzz, ts_xyzzz_xzzzz, ts_xyzzz_yyyy, ts_xyzzz_yyyyy, ts_xyzzz_yyyyz, ts_xyzzz_yyyz, ts_xyzzz_yyyzz, ts_xyzzz_yyzz, ts_xyzzz_yyzzz, ts_xyzzz_yzzz, ts_xyzzz_yzzzz, ts_xyzzz_zzzz, ts_xyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzzz_xxxxx[i] = 6.0 * ts_xyzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxxxy[i] = 6.0 * ts_xyzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxxxz[i] = 6.0 * ts_xyzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxxyy[i] = 6.0 * ts_xyzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxxyz[i] = 6.0 * ts_xyzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxxzz[i] = 6.0 * ts_xyzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxyyy[i] = 6.0 * ts_xyzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxyyz[i] = 6.0 * ts_xyzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxyzz[i] = 6.0 * ts_xyzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xxzzz[i] = 6.0 * ts_xyzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyyyy[i] = 6.0 * ts_xyzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyyyz[i] = 6.0 * ts_xyzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyyzz[i] = 6.0 * ts_xyzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xyzzz[i] = 6.0 * ts_xyzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xzzzz[i] = 6.0 * ts_xyzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyyyy[i] = 6.0 * ts_xyzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyyyz[i] = 6.0 * ts_xyzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyyzz[i] = 6.0 * ts_xyzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yyzzz[i] = 6.0 * ts_xyzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yzzzz[i] = 6.0 * ts_xyzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_zzzzz[i] = 6.0 * ts_xyzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xyzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1176-1197 components of targeted buffer : HH

    auto gs_z_xzzzz_xxxxx = pbuffer.data(idx_g_hh + 1176);

    auto gs_z_xzzzz_xxxxy = pbuffer.data(idx_g_hh + 1177);

    auto gs_z_xzzzz_xxxxz = pbuffer.data(idx_g_hh + 1178);

    auto gs_z_xzzzz_xxxyy = pbuffer.data(idx_g_hh + 1179);

    auto gs_z_xzzzz_xxxyz = pbuffer.data(idx_g_hh + 1180);

    auto gs_z_xzzzz_xxxzz = pbuffer.data(idx_g_hh + 1181);

    auto gs_z_xzzzz_xxyyy = pbuffer.data(idx_g_hh + 1182);

    auto gs_z_xzzzz_xxyyz = pbuffer.data(idx_g_hh + 1183);

    auto gs_z_xzzzz_xxyzz = pbuffer.data(idx_g_hh + 1184);

    auto gs_z_xzzzz_xxzzz = pbuffer.data(idx_g_hh + 1185);

    auto gs_z_xzzzz_xyyyy = pbuffer.data(idx_g_hh + 1186);

    auto gs_z_xzzzz_xyyyz = pbuffer.data(idx_g_hh + 1187);

    auto gs_z_xzzzz_xyyzz = pbuffer.data(idx_g_hh + 1188);

    auto gs_z_xzzzz_xyzzz = pbuffer.data(idx_g_hh + 1189);

    auto gs_z_xzzzz_xzzzz = pbuffer.data(idx_g_hh + 1190);

    auto gs_z_xzzzz_yyyyy = pbuffer.data(idx_g_hh + 1191);

    auto gs_z_xzzzz_yyyyz = pbuffer.data(idx_g_hh + 1192);

    auto gs_z_xzzzz_yyyzz = pbuffer.data(idx_g_hh + 1193);

    auto gs_z_xzzzz_yyzzz = pbuffer.data(idx_g_hh + 1194);

    auto gs_z_xzzzz_yzzzz = pbuffer.data(idx_g_hh + 1195);

    auto gs_z_xzzzz_zzzzz = pbuffer.data(idx_g_hh + 1196);

    #pragma omp simd aligned(gc_z, gs_z_xzzzz_xxxxx, gs_z_xzzzz_xxxxy, gs_z_xzzzz_xxxxz, gs_z_xzzzz_xxxyy, gs_z_xzzzz_xxxyz, gs_z_xzzzz_xxxzz, gs_z_xzzzz_xxyyy, gs_z_xzzzz_xxyyz, gs_z_xzzzz_xxyzz, gs_z_xzzzz_xxzzz, gs_z_xzzzz_xyyyy, gs_z_xzzzz_xyyyz, gs_z_xzzzz_xyyzz, gs_z_xzzzz_xyzzz, gs_z_xzzzz_xzzzz, gs_z_xzzzz_yyyyy, gs_z_xzzzz_yyyyz, gs_z_xzzzz_yyyzz, gs_z_xzzzz_yyzzz, gs_z_xzzzz_yzzzz, gs_z_xzzzz_zzzzz, ts_xzzz_xxxxx, ts_xzzz_xxxxy, ts_xzzz_xxxxz, ts_xzzz_xxxyy, ts_xzzz_xxxyz, ts_xzzz_xxxzz, ts_xzzz_xxyyy, ts_xzzz_xxyyz, ts_xzzz_xxyzz, ts_xzzz_xxzzz, ts_xzzz_xyyyy, ts_xzzz_xyyyz, ts_xzzz_xyyzz, ts_xzzz_xyzzz, ts_xzzz_xzzzz, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyzz, ts_xzzz_yyzzz, ts_xzzz_yzzzz, ts_xzzz_zzzzz, ts_xzzzz_xxxx, ts_xzzzz_xxxxx, ts_xzzzz_xxxxy, ts_xzzzz_xxxxz, ts_xzzzz_xxxy, ts_xzzzz_xxxyy, ts_xzzzz_xxxyz, ts_xzzzz_xxxz, ts_xzzzz_xxxzz, ts_xzzzz_xxyy, ts_xzzzz_xxyyy, ts_xzzzz_xxyyz, ts_xzzzz_xxyz, ts_xzzzz_xxyzz, ts_xzzzz_xxzz, ts_xzzzz_xxzzz, ts_xzzzz_xyyy, ts_xzzzz_xyyyy, ts_xzzzz_xyyyz, ts_xzzzz_xyyz, ts_xzzzz_xyyzz, ts_xzzzz_xyzz, ts_xzzzz_xyzzz, ts_xzzzz_xzzz, ts_xzzzz_xzzzz, ts_xzzzz_yyyy, ts_xzzzz_yyyyy, ts_xzzzz_yyyyz, ts_xzzzz_yyyz, ts_xzzzz_yyyzz, ts_xzzzz_yyzz, ts_xzzzz_yyzzz, ts_xzzzz_yzzz, ts_xzzzz_yzzzz, ts_xzzzz_zzzz, ts_xzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzzz_xxxxx[i] = 8.0 * ts_xzzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxxxy[i] = 8.0 * ts_xzzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxxxz[i] = 8.0 * ts_xzzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxxyy[i] = 8.0 * ts_xzzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxxyz[i] = 8.0 * ts_xzzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxxzz[i] = 8.0 * ts_xzzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxyyy[i] = 8.0 * ts_xzzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxyyz[i] = 8.0 * ts_xzzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxyzz[i] = 8.0 * ts_xzzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xxzzz[i] = 8.0 * ts_xzzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyyyy[i] = 8.0 * ts_xzzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyyyz[i] = 8.0 * ts_xzzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyyzz[i] = 8.0 * ts_xzzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xyzzz[i] = 8.0 * ts_xzzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xzzzz[i] = 8.0 * ts_xzzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyyyy[i] = 8.0 * ts_xzzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyyyz[i] = 8.0 * ts_xzzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyyzz[i] = 8.0 * ts_xzzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yyzzz[i] = 8.0 * ts_xzzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yzzzz[i] = 8.0 * ts_xzzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_zzzzz[i] = 8.0 * ts_xzzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1197-1218 components of targeted buffer : HH

    auto gs_z_yyyyy_xxxxx = pbuffer.data(idx_g_hh + 1197);

    auto gs_z_yyyyy_xxxxy = pbuffer.data(idx_g_hh + 1198);

    auto gs_z_yyyyy_xxxxz = pbuffer.data(idx_g_hh + 1199);

    auto gs_z_yyyyy_xxxyy = pbuffer.data(idx_g_hh + 1200);

    auto gs_z_yyyyy_xxxyz = pbuffer.data(idx_g_hh + 1201);

    auto gs_z_yyyyy_xxxzz = pbuffer.data(idx_g_hh + 1202);

    auto gs_z_yyyyy_xxyyy = pbuffer.data(idx_g_hh + 1203);

    auto gs_z_yyyyy_xxyyz = pbuffer.data(idx_g_hh + 1204);

    auto gs_z_yyyyy_xxyzz = pbuffer.data(idx_g_hh + 1205);

    auto gs_z_yyyyy_xxzzz = pbuffer.data(idx_g_hh + 1206);

    auto gs_z_yyyyy_xyyyy = pbuffer.data(idx_g_hh + 1207);

    auto gs_z_yyyyy_xyyyz = pbuffer.data(idx_g_hh + 1208);

    auto gs_z_yyyyy_xyyzz = pbuffer.data(idx_g_hh + 1209);

    auto gs_z_yyyyy_xyzzz = pbuffer.data(idx_g_hh + 1210);

    auto gs_z_yyyyy_xzzzz = pbuffer.data(idx_g_hh + 1211);

    auto gs_z_yyyyy_yyyyy = pbuffer.data(idx_g_hh + 1212);

    auto gs_z_yyyyy_yyyyz = pbuffer.data(idx_g_hh + 1213);

    auto gs_z_yyyyy_yyyzz = pbuffer.data(idx_g_hh + 1214);

    auto gs_z_yyyyy_yyzzz = pbuffer.data(idx_g_hh + 1215);

    auto gs_z_yyyyy_yzzzz = pbuffer.data(idx_g_hh + 1216);

    auto gs_z_yyyyy_zzzzz = pbuffer.data(idx_g_hh + 1217);

    #pragma omp simd aligned(gc_z, gs_z_yyyyy_xxxxx, gs_z_yyyyy_xxxxy, gs_z_yyyyy_xxxxz, gs_z_yyyyy_xxxyy, gs_z_yyyyy_xxxyz, gs_z_yyyyy_xxxzz, gs_z_yyyyy_xxyyy, gs_z_yyyyy_xxyyz, gs_z_yyyyy_xxyzz, gs_z_yyyyy_xxzzz, gs_z_yyyyy_xyyyy, gs_z_yyyyy_xyyyz, gs_z_yyyyy_xyyzz, gs_z_yyyyy_xyzzz, gs_z_yyyyy_xzzzz, gs_z_yyyyy_yyyyy, gs_z_yyyyy_yyyyz, gs_z_yyyyy_yyyzz, gs_z_yyyyy_yyzzz, gs_z_yyyyy_yzzzz, gs_z_yyyyy_zzzzz, ts_yyyyy_xxxx, ts_yyyyy_xxxxx, ts_yyyyy_xxxxy, ts_yyyyy_xxxxz, ts_yyyyy_xxxy, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxxz, ts_yyyyy_xxxzz, ts_yyyyy_xxyy, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyz, ts_yyyyy_xxyzz, ts_yyyyy_xxzz, ts_yyyyy_xxzzz, ts_yyyyy_xyyy, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzz, ts_yyyyy_xyzzz, ts_yyyyy_xzzz, ts_yyyyy_xzzzz, ts_yyyyy_yyyy, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzz, ts_yyyyy_yzzzz, ts_yyyyy_zzzz, ts_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyy_xxxxx[i] = 2.0 * ts_yyyyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxxxy[i] = 2.0 * ts_yyyyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxxxz[i] = 2.0 * ts_yyyyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxxyy[i] = 2.0 * ts_yyyyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxxyz[i] = 2.0 * ts_yyyyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxxzz[i] = 4.0 * ts_yyyyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxyyy[i] = 2.0 * ts_yyyyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxyyz[i] = 2.0 * ts_yyyyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxyzz[i] = 4.0 * ts_yyyyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xxzzz[i] = 6.0 * ts_yyyyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyyyy[i] = 2.0 * ts_yyyyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyyyz[i] = 2.0 * ts_yyyyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyyzz[i] = 4.0 * ts_yyyyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xyzzz[i] = 6.0 * ts_yyyyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xzzzz[i] = 8.0 * ts_yyyyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyyyy[i] = 2.0 * ts_yyyyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyyyz[i] = 2.0 * ts_yyyyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyyzz[i] = 4.0 * ts_yyyyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yyzzz[i] = 6.0 * ts_yyyyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yzzzz[i] = 8.0 * ts_yyyyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_zzzzz[i] = 10.0 * ts_yyyyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1218-1239 components of targeted buffer : HH

    auto gs_z_yyyyz_xxxxx = pbuffer.data(idx_g_hh + 1218);

    auto gs_z_yyyyz_xxxxy = pbuffer.data(idx_g_hh + 1219);

    auto gs_z_yyyyz_xxxxz = pbuffer.data(idx_g_hh + 1220);

    auto gs_z_yyyyz_xxxyy = pbuffer.data(idx_g_hh + 1221);

    auto gs_z_yyyyz_xxxyz = pbuffer.data(idx_g_hh + 1222);

    auto gs_z_yyyyz_xxxzz = pbuffer.data(idx_g_hh + 1223);

    auto gs_z_yyyyz_xxyyy = pbuffer.data(idx_g_hh + 1224);

    auto gs_z_yyyyz_xxyyz = pbuffer.data(idx_g_hh + 1225);

    auto gs_z_yyyyz_xxyzz = pbuffer.data(idx_g_hh + 1226);

    auto gs_z_yyyyz_xxzzz = pbuffer.data(idx_g_hh + 1227);

    auto gs_z_yyyyz_xyyyy = pbuffer.data(idx_g_hh + 1228);

    auto gs_z_yyyyz_xyyyz = pbuffer.data(idx_g_hh + 1229);

    auto gs_z_yyyyz_xyyzz = pbuffer.data(idx_g_hh + 1230);

    auto gs_z_yyyyz_xyzzz = pbuffer.data(idx_g_hh + 1231);

    auto gs_z_yyyyz_xzzzz = pbuffer.data(idx_g_hh + 1232);

    auto gs_z_yyyyz_yyyyy = pbuffer.data(idx_g_hh + 1233);

    auto gs_z_yyyyz_yyyyz = pbuffer.data(idx_g_hh + 1234);

    auto gs_z_yyyyz_yyyzz = pbuffer.data(idx_g_hh + 1235);

    auto gs_z_yyyyz_yyzzz = pbuffer.data(idx_g_hh + 1236);

    auto gs_z_yyyyz_yzzzz = pbuffer.data(idx_g_hh + 1237);

    auto gs_z_yyyyz_zzzzz = pbuffer.data(idx_g_hh + 1238);

    #pragma omp simd aligned(gc_z, gs_z_yyyyz_xxxxx, gs_z_yyyyz_xxxxy, gs_z_yyyyz_xxxxz, gs_z_yyyyz_xxxyy, gs_z_yyyyz_xxxyz, gs_z_yyyyz_xxxzz, gs_z_yyyyz_xxyyy, gs_z_yyyyz_xxyyz, gs_z_yyyyz_xxyzz, gs_z_yyyyz_xxzzz, gs_z_yyyyz_xyyyy, gs_z_yyyyz_xyyyz, gs_z_yyyyz_xyyzz, gs_z_yyyyz_xyzzz, gs_z_yyyyz_xzzzz, gs_z_yyyyz_yyyyy, gs_z_yyyyz_yyyyz, gs_z_yyyyz_yyyzz, gs_z_yyyyz_yyzzz, gs_z_yyyyz_yzzzz, gs_z_yyyyz_zzzzz, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxzz, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyzz, ts_yyyy_xxzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyzz, ts_yyyy_xyzzz, ts_yyyy_xzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyzz, ts_yyyy_yyzzz, ts_yyyy_yzzzz, ts_yyyy_zzzzz, ts_yyyyz_xxxx, ts_yyyyz_xxxxx, ts_yyyyz_xxxxy, ts_yyyyz_xxxxz, ts_yyyyz_xxxy, ts_yyyyz_xxxyy, ts_yyyyz_xxxyz, ts_yyyyz_xxxz, ts_yyyyz_xxxzz, ts_yyyyz_xxyy, ts_yyyyz_xxyyy, ts_yyyyz_xxyyz, ts_yyyyz_xxyz, ts_yyyyz_xxyzz, ts_yyyyz_xxzz, ts_yyyyz_xxzzz, ts_yyyyz_xyyy, ts_yyyyz_xyyyy, ts_yyyyz_xyyyz, ts_yyyyz_xyyz, ts_yyyyz_xyyzz, ts_yyyyz_xyzz, ts_yyyyz_xyzzz, ts_yyyyz_xzzz, ts_yyyyz_xzzzz, ts_yyyyz_yyyy, ts_yyyyz_yyyyy, ts_yyyyz_yyyyz, ts_yyyyz_yyyz, ts_yyyyz_yyyzz, ts_yyyyz_yyzz, ts_yyyyz_yyzzz, ts_yyyyz_yzzz, ts_yyyyz_yzzzz, ts_yyyyz_zzzz, ts_yyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyz_xxxxx[i] = 2.0 * ts_yyyy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxxxy[i] = 2.0 * ts_yyyy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxxxz[i] = 2.0 * ts_yyyy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxxyy[i] = 2.0 * ts_yyyy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxxyz[i] = 2.0 * ts_yyyy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxxzz[i] = 2.0 * ts_yyyy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxyyy[i] = 2.0 * ts_yyyy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxyyz[i] = 2.0 * ts_yyyy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxyzz[i] = 2.0 * ts_yyyy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xxzzz[i] = 2.0 * ts_yyyy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyyyy[i] = 2.0 * ts_yyyy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyyyz[i] = 2.0 * ts_yyyy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyyzz[i] = 2.0 * ts_yyyy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xyzzz[i] = 2.0 * ts_yyyy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xzzzz[i] = 2.0 * ts_yyyy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyyyy[i] = 2.0 * ts_yyyy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyyyz[i] = 2.0 * ts_yyyy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyyzz[i] = 2.0 * ts_yyyy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yyzzz[i] = 2.0 * ts_yyyy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yzzzz[i] = 2.0 * ts_yyyy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_zzzzz[i] = 2.0 * ts_yyyy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yyyyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1239-1260 components of targeted buffer : HH

    auto gs_z_yyyzz_xxxxx = pbuffer.data(idx_g_hh + 1239);

    auto gs_z_yyyzz_xxxxy = pbuffer.data(idx_g_hh + 1240);

    auto gs_z_yyyzz_xxxxz = pbuffer.data(idx_g_hh + 1241);

    auto gs_z_yyyzz_xxxyy = pbuffer.data(idx_g_hh + 1242);

    auto gs_z_yyyzz_xxxyz = pbuffer.data(idx_g_hh + 1243);

    auto gs_z_yyyzz_xxxzz = pbuffer.data(idx_g_hh + 1244);

    auto gs_z_yyyzz_xxyyy = pbuffer.data(idx_g_hh + 1245);

    auto gs_z_yyyzz_xxyyz = pbuffer.data(idx_g_hh + 1246);

    auto gs_z_yyyzz_xxyzz = pbuffer.data(idx_g_hh + 1247);

    auto gs_z_yyyzz_xxzzz = pbuffer.data(idx_g_hh + 1248);

    auto gs_z_yyyzz_xyyyy = pbuffer.data(idx_g_hh + 1249);

    auto gs_z_yyyzz_xyyyz = pbuffer.data(idx_g_hh + 1250);

    auto gs_z_yyyzz_xyyzz = pbuffer.data(idx_g_hh + 1251);

    auto gs_z_yyyzz_xyzzz = pbuffer.data(idx_g_hh + 1252);

    auto gs_z_yyyzz_xzzzz = pbuffer.data(idx_g_hh + 1253);

    auto gs_z_yyyzz_yyyyy = pbuffer.data(idx_g_hh + 1254);

    auto gs_z_yyyzz_yyyyz = pbuffer.data(idx_g_hh + 1255);

    auto gs_z_yyyzz_yyyzz = pbuffer.data(idx_g_hh + 1256);

    auto gs_z_yyyzz_yyzzz = pbuffer.data(idx_g_hh + 1257);

    auto gs_z_yyyzz_yzzzz = pbuffer.data(idx_g_hh + 1258);

    auto gs_z_yyyzz_zzzzz = pbuffer.data(idx_g_hh + 1259);

    #pragma omp simd aligned(gc_z, gs_z_yyyzz_xxxxx, gs_z_yyyzz_xxxxy, gs_z_yyyzz_xxxxz, gs_z_yyyzz_xxxyy, gs_z_yyyzz_xxxyz, gs_z_yyyzz_xxxzz, gs_z_yyyzz_xxyyy, gs_z_yyyzz_xxyyz, gs_z_yyyzz_xxyzz, gs_z_yyyzz_xxzzz, gs_z_yyyzz_xyyyy, gs_z_yyyzz_xyyyz, gs_z_yyyzz_xyyzz, gs_z_yyyzz_xyzzz, gs_z_yyyzz_xzzzz, gs_z_yyyzz_yyyyy, gs_z_yyyzz_yyyyz, gs_z_yyyzz_yyyzz, gs_z_yyyzz_yyzzz, gs_z_yyyzz_yzzzz, gs_z_yyyzz_zzzzz, ts_yyyz_xxxxx, ts_yyyz_xxxxy, ts_yyyz_xxxxz, ts_yyyz_xxxyy, ts_yyyz_xxxyz, ts_yyyz_xxxzz, ts_yyyz_xxyyy, ts_yyyz_xxyyz, ts_yyyz_xxyzz, ts_yyyz_xxzzz, ts_yyyz_xyyyy, ts_yyyz_xyyyz, ts_yyyz_xyyzz, ts_yyyz_xyzzz, ts_yyyz_xzzzz, ts_yyyz_yyyyy, ts_yyyz_yyyyz, ts_yyyz_yyyzz, ts_yyyz_yyzzz, ts_yyyz_yzzzz, ts_yyyz_zzzzz, ts_yyyzz_xxxx, ts_yyyzz_xxxxx, ts_yyyzz_xxxxy, ts_yyyzz_xxxxz, ts_yyyzz_xxxy, ts_yyyzz_xxxyy, ts_yyyzz_xxxyz, ts_yyyzz_xxxz, ts_yyyzz_xxxzz, ts_yyyzz_xxyy, ts_yyyzz_xxyyy, ts_yyyzz_xxyyz, ts_yyyzz_xxyz, ts_yyyzz_xxyzz, ts_yyyzz_xxzz, ts_yyyzz_xxzzz, ts_yyyzz_xyyy, ts_yyyzz_xyyyy, ts_yyyzz_xyyyz, ts_yyyzz_xyyz, ts_yyyzz_xyyzz, ts_yyyzz_xyzz, ts_yyyzz_xyzzz, ts_yyyzz_xzzz, ts_yyyzz_xzzzz, ts_yyyzz_yyyy, ts_yyyzz_yyyyy, ts_yyyzz_yyyyz, ts_yyyzz_yyyz, ts_yyyzz_yyyzz, ts_yyyzz_yyzz, ts_yyyzz_yyzzz, ts_yyyzz_yzzz, ts_yyyzz_yzzzz, ts_yyyzz_zzzz, ts_yyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyzz_xxxxx[i] = 4.0 * ts_yyyz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxxxy[i] = 4.0 * ts_yyyz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxxxz[i] = 4.0 * ts_yyyz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxxyy[i] = 4.0 * ts_yyyz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxxyz[i] = 4.0 * ts_yyyz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxxzz[i] = 4.0 * ts_yyyz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxyyy[i] = 4.0 * ts_yyyz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxyyz[i] = 4.0 * ts_yyyz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxyzz[i] = 4.0 * ts_yyyz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xxzzz[i] = 4.0 * ts_yyyz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyyyy[i] = 4.0 * ts_yyyz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyyyz[i] = 4.0 * ts_yyyz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyyzz[i] = 4.0 * ts_yyyz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xyzzz[i] = 4.0 * ts_yyyz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xzzzz[i] = 4.0 * ts_yyyz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyyyy[i] = 4.0 * ts_yyyz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyyyz[i] = 4.0 * ts_yyyz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyyzz[i] = 4.0 * ts_yyyz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yyzzz[i] = 4.0 * ts_yyyz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyyzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yzzzz[i] = 4.0 * ts_yyyz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyyzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_zzzzz[i] = 4.0 * ts_yyyz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yyyzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1260-1281 components of targeted buffer : HH

    auto gs_z_yyzzz_xxxxx = pbuffer.data(idx_g_hh + 1260);

    auto gs_z_yyzzz_xxxxy = pbuffer.data(idx_g_hh + 1261);

    auto gs_z_yyzzz_xxxxz = pbuffer.data(idx_g_hh + 1262);

    auto gs_z_yyzzz_xxxyy = pbuffer.data(idx_g_hh + 1263);

    auto gs_z_yyzzz_xxxyz = pbuffer.data(idx_g_hh + 1264);

    auto gs_z_yyzzz_xxxzz = pbuffer.data(idx_g_hh + 1265);

    auto gs_z_yyzzz_xxyyy = pbuffer.data(idx_g_hh + 1266);

    auto gs_z_yyzzz_xxyyz = pbuffer.data(idx_g_hh + 1267);

    auto gs_z_yyzzz_xxyzz = pbuffer.data(idx_g_hh + 1268);

    auto gs_z_yyzzz_xxzzz = pbuffer.data(idx_g_hh + 1269);

    auto gs_z_yyzzz_xyyyy = pbuffer.data(idx_g_hh + 1270);

    auto gs_z_yyzzz_xyyyz = pbuffer.data(idx_g_hh + 1271);

    auto gs_z_yyzzz_xyyzz = pbuffer.data(idx_g_hh + 1272);

    auto gs_z_yyzzz_xyzzz = pbuffer.data(idx_g_hh + 1273);

    auto gs_z_yyzzz_xzzzz = pbuffer.data(idx_g_hh + 1274);

    auto gs_z_yyzzz_yyyyy = pbuffer.data(idx_g_hh + 1275);

    auto gs_z_yyzzz_yyyyz = pbuffer.data(idx_g_hh + 1276);

    auto gs_z_yyzzz_yyyzz = pbuffer.data(idx_g_hh + 1277);

    auto gs_z_yyzzz_yyzzz = pbuffer.data(idx_g_hh + 1278);

    auto gs_z_yyzzz_yzzzz = pbuffer.data(idx_g_hh + 1279);

    auto gs_z_yyzzz_zzzzz = pbuffer.data(idx_g_hh + 1280);

    #pragma omp simd aligned(gc_z, gs_z_yyzzz_xxxxx, gs_z_yyzzz_xxxxy, gs_z_yyzzz_xxxxz, gs_z_yyzzz_xxxyy, gs_z_yyzzz_xxxyz, gs_z_yyzzz_xxxzz, gs_z_yyzzz_xxyyy, gs_z_yyzzz_xxyyz, gs_z_yyzzz_xxyzz, gs_z_yyzzz_xxzzz, gs_z_yyzzz_xyyyy, gs_z_yyzzz_xyyyz, gs_z_yyzzz_xyyzz, gs_z_yyzzz_xyzzz, gs_z_yyzzz_xzzzz, gs_z_yyzzz_yyyyy, gs_z_yyzzz_yyyyz, gs_z_yyzzz_yyyzz, gs_z_yyzzz_yyzzz, gs_z_yyzzz_yzzzz, gs_z_yyzzz_zzzzz, ts_yyzz_xxxxx, ts_yyzz_xxxxy, ts_yyzz_xxxxz, ts_yyzz_xxxyy, ts_yyzz_xxxyz, ts_yyzz_xxxzz, ts_yyzz_xxyyy, ts_yyzz_xxyyz, ts_yyzz_xxyzz, ts_yyzz_xxzzz, ts_yyzz_xyyyy, ts_yyzz_xyyyz, ts_yyzz_xyyzz, ts_yyzz_xyzzz, ts_yyzz_xzzzz, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyzz, ts_yyzz_yyzzz, ts_yyzz_yzzzz, ts_yyzz_zzzzz, ts_yyzzz_xxxx, ts_yyzzz_xxxxx, ts_yyzzz_xxxxy, ts_yyzzz_xxxxz, ts_yyzzz_xxxy, ts_yyzzz_xxxyy, ts_yyzzz_xxxyz, ts_yyzzz_xxxz, ts_yyzzz_xxxzz, ts_yyzzz_xxyy, ts_yyzzz_xxyyy, ts_yyzzz_xxyyz, ts_yyzzz_xxyz, ts_yyzzz_xxyzz, ts_yyzzz_xxzz, ts_yyzzz_xxzzz, ts_yyzzz_xyyy, ts_yyzzz_xyyyy, ts_yyzzz_xyyyz, ts_yyzzz_xyyz, ts_yyzzz_xyyzz, ts_yyzzz_xyzz, ts_yyzzz_xyzzz, ts_yyzzz_xzzz, ts_yyzzz_xzzzz, ts_yyzzz_yyyy, ts_yyzzz_yyyyy, ts_yyzzz_yyyyz, ts_yyzzz_yyyz, ts_yyzzz_yyyzz, ts_yyzzz_yyzz, ts_yyzzz_yyzzz, ts_yyzzz_yzzz, ts_yyzzz_yzzzz, ts_yyzzz_zzzz, ts_yyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzzz_xxxxx[i] = 6.0 * ts_yyzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxxxy[i] = 6.0 * ts_yyzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxxxz[i] = 6.0 * ts_yyzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxxyy[i] = 6.0 * ts_yyzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxxyz[i] = 6.0 * ts_yyzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxxzz[i] = 6.0 * ts_yyzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxyyy[i] = 6.0 * ts_yyzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxyyz[i] = 6.0 * ts_yyzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxyzz[i] = 6.0 * ts_yyzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xxzzz[i] = 6.0 * ts_yyzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyyyy[i] = 6.0 * ts_yyzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyyyz[i] = 6.0 * ts_yyzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyyzz[i] = 6.0 * ts_yyzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xyzzz[i] = 6.0 * ts_yyzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xzzzz[i] = 6.0 * ts_yyzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyyyy[i] = 6.0 * ts_yyzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyyyz[i] = 6.0 * ts_yyzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyyzz[i] = 6.0 * ts_yyzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yyzzz[i] = 6.0 * ts_yyzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yzzzz[i] = 6.0 * ts_yyzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_zzzzz[i] = 6.0 * ts_yyzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yyzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1281-1302 components of targeted buffer : HH

    auto gs_z_yzzzz_xxxxx = pbuffer.data(idx_g_hh + 1281);

    auto gs_z_yzzzz_xxxxy = pbuffer.data(idx_g_hh + 1282);

    auto gs_z_yzzzz_xxxxz = pbuffer.data(idx_g_hh + 1283);

    auto gs_z_yzzzz_xxxyy = pbuffer.data(idx_g_hh + 1284);

    auto gs_z_yzzzz_xxxyz = pbuffer.data(idx_g_hh + 1285);

    auto gs_z_yzzzz_xxxzz = pbuffer.data(idx_g_hh + 1286);

    auto gs_z_yzzzz_xxyyy = pbuffer.data(idx_g_hh + 1287);

    auto gs_z_yzzzz_xxyyz = pbuffer.data(idx_g_hh + 1288);

    auto gs_z_yzzzz_xxyzz = pbuffer.data(idx_g_hh + 1289);

    auto gs_z_yzzzz_xxzzz = pbuffer.data(idx_g_hh + 1290);

    auto gs_z_yzzzz_xyyyy = pbuffer.data(idx_g_hh + 1291);

    auto gs_z_yzzzz_xyyyz = pbuffer.data(idx_g_hh + 1292);

    auto gs_z_yzzzz_xyyzz = pbuffer.data(idx_g_hh + 1293);

    auto gs_z_yzzzz_xyzzz = pbuffer.data(idx_g_hh + 1294);

    auto gs_z_yzzzz_xzzzz = pbuffer.data(idx_g_hh + 1295);

    auto gs_z_yzzzz_yyyyy = pbuffer.data(idx_g_hh + 1296);

    auto gs_z_yzzzz_yyyyz = pbuffer.data(idx_g_hh + 1297);

    auto gs_z_yzzzz_yyyzz = pbuffer.data(idx_g_hh + 1298);

    auto gs_z_yzzzz_yyzzz = pbuffer.data(idx_g_hh + 1299);

    auto gs_z_yzzzz_yzzzz = pbuffer.data(idx_g_hh + 1300);

    auto gs_z_yzzzz_zzzzz = pbuffer.data(idx_g_hh + 1301);

    #pragma omp simd aligned(gc_z, gs_z_yzzzz_xxxxx, gs_z_yzzzz_xxxxy, gs_z_yzzzz_xxxxz, gs_z_yzzzz_xxxyy, gs_z_yzzzz_xxxyz, gs_z_yzzzz_xxxzz, gs_z_yzzzz_xxyyy, gs_z_yzzzz_xxyyz, gs_z_yzzzz_xxyzz, gs_z_yzzzz_xxzzz, gs_z_yzzzz_xyyyy, gs_z_yzzzz_xyyyz, gs_z_yzzzz_xyyzz, gs_z_yzzzz_xyzzz, gs_z_yzzzz_xzzzz, gs_z_yzzzz_yyyyy, gs_z_yzzzz_yyyyz, gs_z_yzzzz_yyyzz, gs_z_yzzzz_yyzzz, gs_z_yzzzz_yzzzz, gs_z_yzzzz_zzzzz, ts_yzzz_xxxxx, ts_yzzz_xxxxy, ts_yzzz_xxxxz, ts_yzzz_xxxyy, ts_yzzz_xxxyz, ts_yzzz_xxxzz, ts_yzzz_xxyyy, ts_yzzz_xxyyz, ts_yzzz_xxyzz, ts_yzzz_xxzzz, ts_yzzz_xyyyy, ts_yzzz_xyyyz, ts_yzzz_xyyzz, ts_yzzz_xyzzz, ts_yzzz_xzzzz, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyzz, ts_yzzz_yyzzz, ts_yzzz_yzzzz, ts_yzzz_zzzzz, ts_yzzzz_xxxx, ts_yzzzz_xxxxx, ts_yzzzz_xxxxy, ts_yzzzz_xxxxz, ts_yzzzz_xxxy, ts_yzzzz_xxxyy, ts_yzzzz_xxxyz, ts_yzzzz_xxxz, ts_yzzzz_xxxzz, ts_yzzzz_xxyy, ts_yzzzz_xxyyy, ts_yzzzz_xxyyz, ts_yzzzz_xxyz, ts_yzzzz_xxyzz, ts_yzzzz_xxzz, ts_yzzzz_xxzzz, ts_yzzzz_xyyy, ts_yzzzz_xyyyy, ts_yzzzz_xyyyz, ts_yzzzz_xyyz, ts_yzzzz_xyyzz, ts_yzzzz_xyzz, ts_yzzzz_xyzzz, ts_yzzzz_xzzz, ts_yzzzz_xzzzz, ts_yzzzz_yyyy, ts_yzzzz_yyyyy, ts_yzzzz_yyyyz, ts_yzzzz_yyyz, ts_yzzzz_yyyzz, ts_yzzzz_yyzz, ts_yzzzz_yyzzz, ts_yzzzz_yzzz, ts_yzzzz_yzzzz, ts_yzzzz_zzzz, ts_yzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzzz_xxxxx[i] = 8.0 * ts_yzzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxxxy[i] = 8.0 * ts_yzzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxxxz[i] = 8.0 * ts_yzzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxxyy[i] = 8.0 * ts_yzzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxxyz[i] = 8.0 * ts_yzzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxxzz[i] = 8.0 * ts_yzzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxyyy[i] = 8.0 * ts_yzzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxyyz[i] = 8.0 * ts_yzzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxyzz[i] = 8.0 * ts_yzzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xxzzz[i] = 8.0 * ts_yzzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyyyy[i] = 8.0 * ts_yzzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyyyz[i] = 8.0 * ts_yzzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyyzz[i] = 8.0 * ts_yzzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xyzzz[i] = 8.0 * ts_yzzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xzzzz[i] = 8.0 * ts_yzzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyyyy[i] = 8.0 * ts_yzzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyyyz[i] = 8.0 * ts_yzzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyyzz[i] = 8.0 * ts_yzzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yyzzz[i] = 8.0 * ts_yzzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yzzzz[i] = 8.0 * ts_yzzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_zzzzz[i] = 8.0 * ts_yzzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 1302-1323 components of targeted buffer : HH

    auto gs_z_zzzzz_xxxxx = pbuffer.data(idx_g_hh + 1302);

    auto gs_z_zzzzz_xxxxy = pbuffer.data(idx_g_hh + 1303);

    auto gs_z_zzzzz_xxxxz = pbuffer.data(idx_g_hh + 1304);

    auto gs_z_zzzzz_xxxyy = pbuffer.data(idx_g_hh + 1305);

    auto gs_z_zzzzz_xxxyz = pbuffer.data(idx_g_hh + 1306);

    auto gs_z_zzzzz_xxxzz = pbuffer.data(idx_g_hh + 1307);

    auto gs_z_zzzzz_xxyyy = pbuffer.data(idx_g_hh + 1308);

    auto gs_z_zzzzz_xxyyz = pbuffer.data(idx_g_hh + 1309);

    auto gs_z_zzzzz_xxyzz = pbuffer.data(idx_g_hh + 1310);

    auto gs_z_zzzzz_xxzzz = pbuffer.data(idx_g_hh + 1311);

    auto gs_z_zzzzz_xyyyy = pbuffer.data(idx_g_hh + 1312);

    auto gs_z_zzzzz_xyyyz = pbuffer.data(idx_g_hh + 1313);

    auto gs_z_zzzzz_xyyzz = pbuffer.data(idx_g_hh + 1314);

    auto gs_z_zzzzz_xyzzz = pbuffer.data(idx_g_hh + 1315);

    auto gs_z_zzzzz_xzzzz = pbuffer.data(idx_g_hh + 1316);

    auto gs_z_zzzzz_yyyyy = pbuffer.data(idx_g_hh + 1317);

    auto gs_z_zzzzz_yyyyz = pbuffer.data(idx_g_hh + 1318);

    auto gs_z_zzzzz_yyyzz = pbuffer.data(idx_g_hh + 1319);

    auto gs_z_zzzzz_yyzzz = pbuffer.data(idx_g_hh + 1320);

    auto gs_z_zzzzz_yzzzz = pbuffer.data(idx_g_hh + 1321);

    auto gs_z_zzzzz_zzzzz = pbuffer.data(idx_g_hh + 1322);

    #pragma omp simd aligned(gc_z, gs_z_zzzzz_xxxxx, gs_z_zzzzz_xxxxy, gs_z_zzzzz_xxxxz, gs_z_zzzzz_xxxyy, gs_z_zzzzz_xxxyz, gs_z_zzzzz_xxxzz, gs_z_zzzzz_xxyyy, gs_z_zzzzz_xxyyz, gs_z_zzzzz_xxyzz, gs_z_zzzzz_xxzzz, gs_z_zzzzz_xyyyy, gs_z_zzzzz_xyyyz, gs_z_zzzzz_xyyzz, gs_z_zzzzz_xyzzz, gs_z_zzzzz_xzzzz, gs_z_zzzzz_yyyyy, gs_z_zzzzz_yyyyz, gs_z_zzzzz_yyyzz, gs_z_zzzzz_yyzzz, gs_z_zzzzz_yzzzz, gs_z_zzzzz_zzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, ts_zzzzz_xxxx, ts_zzzzz_xxxxx, ts_zzzzz_xxxxy, ts_zzzzz_xxxxz, ts_zzzzz_xxxy, ts_zzzzz_xxxyy, ts_zzzzz_xxxyz, ts_zzzzz_xxxz, ts_zzzzz_xxxzz, ts_zzzzz_xxyy, ts_zzzzz_xxyyy, ts_zzzzz_xxyyz, ts_zzzzz_xxyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyy, ts_zzzzz_xyyyy, ts_zzzzz_xyyyz, ts_zzzzz_xyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyy, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzzz_xxxxx[i] = 10.0 * ts_zzzz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxxxy[i] = 10.0 * ts_zzzz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxxxz[i] = 10.0 * ts_zzzz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxxyy[i] = 10.0 * ts_zzzz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxxyz[i] = 10.0 * ts_zzzz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxxzz[i] = 10.0 * ts_zzzz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxyyy[i] = 10.0 * ts_zzzz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxyyz[i] = 10.0 * ts_zzzz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxyzz[i] = 10.0 * ts_zzzz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xxzzz[i] = 10.0 * ts_zzzz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyyyy[i] = 10.0 * ts_zzzz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyyyz[i] = 10.0 * ts_zzzz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyyzz[i] = 10.0 * ts_zzzz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xyzzz[i] = 10.0 * ts_zzzz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xzzzz[i] = 10.0 * ts_zzzz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzzzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyyyy[i] = 10.0 * ts_zzzz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyyyz[i] = 10.0 * ts_zzzz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyyzz[i] = 10.0 * ts_zzzz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yyzzz[i] = 10.0 * ts_zzzz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzzzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yzzzz[i] = 10.0 * ts_zzzz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzzzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_zzzzz[i] = 10.0 * ts_zzzz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_zzzzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_zzzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

