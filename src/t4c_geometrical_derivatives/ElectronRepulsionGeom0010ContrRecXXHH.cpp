#include "ElectronRepulsionGeom0010ContrRecXXHH.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxhh(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxhh,
                                            const size_t idx_xxgh,
                                            const size_t idx_geom_10_xxgh,
                                            const size_t idx_geom_10_xxgi,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSGH

            const auto gh_off = idx_xxgh + (i * bcomps + j) * 315;

            auto g_xxxx_xxxxx = cbuffer.data(gh_off + 0);

            auto g_xxxx_xxxxy = cbuffer.data(gh_off + 1);

            auto g_xxxx_xxxxz = cbuffer.data(gh_off + 2);

            auto g_xxxx_xxxyy = cbuffer.data(gh_off + 3);

            auto g_xxxx_xxxyz = cbuffer.data(gh_off + 4);

            auto g_xxxx_xxxzz = cbuffer.data(gh_off + 5);

            auto g_xxxx_xxyyy = cbuffer.data(gh_off + 6);

            auto g_xxxx_xxyyz = cbuffer.data(gh_off + 7);

            auto g_xxxx_xxyzz = cbuffer.data(gh_off + 8);

            auto g_xxxx_xxzzz = cbuffer.data(gh_off + 9);

            auto g_xxxx_xyyyy = cbuffer.data(gh_off + 10);

            auto g_xxxx_xyyyz = cbuffer.data(gh_off + 11);

            auto g_xxxx_xyyzz = cbuffer.data(gh_off + 12);

            auto g_xxxx_xyzzz = cbuffer.data(gh_off + 13);

            auto g_xxxx_xzzzz = cbuffer.data(gh_off + 14);

            auto g_xxxx_yyyyy = cbuffer.data(gh_off + 15);

            auto g_xxxx_yyyyz = cbuffer.data(gh_off + 16);

            auto g_xxxx_yyyzz = cbuffer.data(gh_off + 17);

            auto g_xxxx_yyzzz = cbuffer.data(gh_off + 18);

            auto g_xxxx_yzzzz = cbuffer.data(gh_off + 19);

            auto g_xxxx_zzzzz = cbuffer.data(gh_off + 20);

            auto g_xxxy_xxxxx = cbuffer.data(gh_off + 21);

            auto g_xxxy_xxxxy = cbuffer.data(gh_off + 22);

            auto g_xxxy_xxxxz = cbuffer.data(gh_off + 23);

            auto g_xxxy_xxxyy = cbuffer.data(gh_off + 24);

            auto g_xxxy_xxxyz = cbuffer.data(gh_off + 25);

            auto g_xxxy_xxxzz = cbuffer.data(gh_off + 26);

            auto g_xxxy_xxyyy = cbuffer.data(gh_off + 27);

            auto g_xxxy_xxyyz = cbuffer.data(gh_off + 28);

            auto g_xxxy_xxyzz = cbuffer.data(gh_off + 29);

            auto g_xxxy_xxzzz = cbuffer.data(gh_off + 30);

            auto g_xxxy_xyyyy = cbuffer.data(gh_off + 31);

            auto g_xxxy_xyyyz = cbuffer.data(gh_off + 32);

            auto g_xxxy_xyyzz = cbuffer.data(gh_off + 33);

            auto g_xxxy_xyzzz = cbuffer.data(gh_off + 34);

            auto g_xxxy_xzzzz = cbuffer.data(gh_off + 35);

            auto g_xxxy_yyyyy = cbuffer.data(gh_off + 36);

            auto g_xxxy_yyyyz = cbuffer.data(gh_off + 37);

            auto g_xxxy_yyyzz = cbuffer.data(gh_off + 38);

            auto g_xxxy_yyzzz = cbuffer.data(gh_off + 39);

            auto g_xxxy_yzzzz = cbuffer.data(gh_off + 40);

            auto g_xxxy_zzzzz = cbuffer.data(gh_off + 41);

            auto g_xxxz_xxxxx = cbuffer.data(gh_off + 42);

            auto g_xxxz_xxxxy = cbuffer.data(gh_off + 43);

            auto g_xxxz_xxxxz = cbuffer.data(gh_off + 44);

            auto g_xxxz_xxxyy = cbuffer.data(gh_off + 45);

            auto g_xxxz_xxxyz = cbuffer.data(gh_off + 46);

            auto g_xxxz_xxxzz = cbuffer.data(gh_off + 47);

            auto g_xxxz_xxyyy = cbuffer.data(gh_off + 48);

            auto g_xxxz_xxyyz = cbuffer.data(gh_off + 49);

            auto g_xxxz_xxyzz = cbuffer.data(gh_off + 50);

            auto g_xxxz_xxzzz = cbuffer.data(gh_off + 51);

            auto g_xxxz_xyyyy = cbuffer.data(gh_off + 52);

            auto g_xxxz_xyyyz = cbuffer.data(gh_off + 53);

            auto g_xxxz_xyyzz = cbuffer.data(gh_off + 54);

            auto g_xxxz_xyzzz = cbuffer.data(gh_off + 55);

            auto g_xxxz_xzzzz = cbuffer.data(gh_off + 56);

            auto g_xxxz_yyyyy = cbuffer.data(gh_off + 57);

            auto g_xxxz_yyyyz = cbuffer.data(gh_off + 58);

            auto g_xxxz_yyyzz = cbuffer.data(gh_off + 59);

            auto g_xxxz_yyzzz = cbuffer.data(gh_off + 60);

            auto g_xxxz_yzzzz = cbuffer.data(gh_off + 61);

            auto g_xxxz_zzzzz = cbuffer.data(gh_off + 62);

            auto g_xxyy_xxxxx = cbuffer.data(gh_off + 63);

            auto g_xxyy_xxxxy = cbuffer.data(gh_off + 64);

            auto g_xxyy_xxxxz = cbuffer.data(gh_off + 65);

            auto g_xxyy_xxxyy = cbuffer.data(gh_off + 66);

            auto g_xxyy_xxxyz = cbuffer.data(gh_off + 67);

            auto g_xxyy_xxxzz = cbuffer.data(gh_off + 68);

            auto g_xxyy_xxyyy = cbuffer.data(gh_off + 69);

            auto g_xxyy_xxyyz = cbuffer.data(gh_off + 70);

            auto g_xxyy_xxyzz = cbuffer.data(gh_off + 71);

            auto g_xxyy_xxzzz = cbuffer.data(gh_off + 72);

            auto g_xxyy_xyyyy = cbuffer.data(gh_off + 73);

            auto g_xxyy_xyyyz = cbuffer.data(gh_off + 74);

            auto g_xxyy_xyyzz = cbuffer.data(gh_off + 75);

            auto g_xxyy_xyzzz = cbuffer.data(gh_off + 76);

            auto g_xxyy_xzzzz = cbuffer.data(gh_off + 77);

            auto g_xxyy_yyyyy = cbuffer.data(gh_off + 78);

            auto g_xxyy_yyyyz = cbuffer.data(gh_off + 79);

            auto g_xxyy_yyyzz = cbuffer.data(gh_off + 80);

            auto g_xxyy_yyzzz = cbuffer.data(gh_off + 81);

            auto g_xxyy_yzzzz = cbuffer.data(gh_off + 82);

            auto g_xxyy_zzzzz = cbuffer.data(gh_off + 83);

            auto g_xxyz_xxxxx = cbuffer.data(gh_off + 84);

            auto g_xxyz_xxxxy = cbuffer.data(gh_off + 85);

            auto g_xxyz_xxxxz = cbuffer.data(gh_off + 86);

            auto g_xxyz_xxxyy = cbuffer.data(gh_off + 87);

            auto g_xxyz_xxxyz = cbuffer.data(gh_off + 88);

            auto g_xxyz_xxxzz = cbuffer.data(gh_off + 89);

            auto g_xxyz_xxyyy = cbuffer.data(gh_off + 90);

            auto g_xxyz_xxyyz = cbuffer.data(gh_off + 91);

            auto g_xxyz_xxyzz = cbuffer.data(gh_off + 92);

            auto g_xxyz_xxzzz = cbuffer.data(gh_off + 93);

            auto g_xxyz_xyyyy = cbuffer.data(gh_off + 94);

            auto g_xxyz_xyyyz = cbuffer.data(gh_off + 95);

            auto g_xxyz_xyyzz = cbuffer.data(gh_off + 96);

            auto g_xxyz_xyzzz = cbuffer.data(gh_off + 97);

            auto g_xxyz_xzzzz = cbuffer.data(gh_off + 98);

            auto g_xxyz_yyyyy = cbuffer.data(gh_off + 99);

            auto g_xxyz_yyyyz = cbuffer.data(gh_off + 100);

            auto g_xxyz_yyyzz = cbuffer.data(gh_off + 101);

            auto g_xxyz_yyzzz = cbuffer.data(gh_off + 102);

            auto g_xxyz_yzzzz = cbuffer.data(gh_off + 103);

            auto g_xxyz_zzzzz = cbuffer.data(gh_off + 104);

            auto g_xxzz_xxxxx = cbuffer.data(gh_off + 105);

            auto g_xxzz_xxxxy = cbuffer.data(gh_off + 106);

            auto g_xxzz_xxxxz = cbuffer.data(gh_off + 107);

            auto g_xxzz_xxxyy = cbuffer.data(gh_off + 108);

            auto g_xxzz_xxxyz = cbuffer.data(gh_off + 109);

            auto g_xxzz_xxxzz = cbuffer.data(gh_off + 110);

            auto g_xxzz_xxyyy = cbuffer.data(gh_off + 111);

            auto g_xxzz_xxyyz = cbuffer.data(gh_off + 112);

            auto g_xxzz_xxyzz = cbuffer.data(gh_off + 113);

            auto g_xxzz_xxzzz = cbuffer.data(gh_off + 114);

            auto g_xxzz_xyyyy = cbuffer.data(gh_off + 115);

            auto g_xxzz_xyyyz = cbuffer.data(gh_off + 116);

            auto g_xxzz_xyyzz = cbuffer.data(gh_off + 117);

            auto g_xxzz_xyzzz = cbuffer.data(gh_off + 118);

            auto g_xxzz_xzzzz = cbuffer.data(gh_off + 119);

            auto g_xxzz_yyyyy = cbuffer.data(gh_off + 120);

            auto g_xxzz_yyyyz = cbuffer.data(gh_off + 121);

            auto g_xxzz_yyyzz = cbuffer.data(gh_off + 122);

            auto g_xxzz_yyzzz = cbuffer.data(gh_off + 123);

            auto g_xxzz_yzzzz = cbuffer.data(gh_off + 124);

            auto g_xxzz_zzzzz = cbuffer.data(gh_off + 125);

            auto g_xyyy_xxxxx = cbuffer.data(gh_off + 126);

            auto g_xyyy_xxxxy = cbuffer.data(gh_off + 127);

            auto g_xyyy_xxxxz = cbuffer.data(gh_off + 128);

            auto g_xyyy_xxxyy = cbuffer.data(gh_off + 129);

            auto g_xyyy_xxxyz = cbuffer.data(gh_off + 130);

            auto g_xyyy_xxxzz = cbuffer.data(gh_off + 131);

            auto g_xyyy_xxyyy = cbuffer.data(gh_off + 132);

            auto g_xyyy_xxyyz = cbuffer.data(gh_off + 133);

            auto g_xyyy_xxyzz = cbuffer.data(gh_off + 134);

            auto g_xyyy_xxzzz = cbuffer.data(gh_off + 135);

            auto g_xyyy_xyyyy = cbuffer.data(gh_off + 136);

            auto g_xyyy_xyyyz = cbuffer.data(gh_off + 137);

            auto g_xyyy_xyyzz = cbuffer.data(gh_off + 138);

            auto g_xyyy_xyzzz = cbuffer.data(gh_off + 139);

            auto g_xyyy_xzzzz = cbuffer.data(gh_off + 140);

            auto g_xyyy_yyyyy = cbuffer.data(gh_off + 141);

            auto g_xyyy_yyyyz = cbuffer.data(gh_off + 142);

            auto g_xyyy_yyyzz = cbuffer.data(gh_off + 143);

            auto g_xyyy_yyzzz = cbuffer.data(gh_off + 144);

            auto g_xyyy_yzzzz = cbuffer.data(gh_off + 145);

            auto g_xyyy_zzzzz = cbuffer.data(gh_off + 146);

            auto g_xyyz_xxxxx = cbuffer.data(gh_off + 147);

            auto g_xyyz_xxxxy = cbuffer.data(gh_off + 148);

            auto g_xyyz_xxxxz = cbuffer.data(gh_off + 149);

            auto g_xyyz_xxxyy = cbuffer.data(gh_off + 150);

            auto g_xyyz_xxxyz = cbuffer.data(gh_off + 151);

            auto g_xyyz_xxxzz = cbuffer.data(gh_off + 152);

            auto g_xyyz_xxyyy = cbuffer.data(gh_off + 153);

            auto g_xyyz_xxyyz = cbuffer.data(gh_off + 154);

            auto g_xyyz_xxyzz = cbuffer.data(gh_off + 155);

            auto g_xyyz_xxzzz = cbuffer.data(gh_off + 156);

            auto g_xyyz_xyyyy = cbuffer.data(gh_off + 157);

            auto g_xyyz_xyyyz = cbuffer.data(gh_off + 158);

            auto g_xyyz_xyyzz = cbuffer.data(gh_off + 159);

            auto g_xyyz_xyzzz = cbuffer.data(gh_off + 160);

            auto g_xyyz_xzzzz = cbuffer.data(gh_off + 161);

            auto g_xyyz_yyyyy = cbuffer.data(gh_off + 162);

            auto g_xyyz_yyyyz = cbuffer.data(gh_off + 163);

            auto g_xyyz_yyyzz = cbuffer.data(gh_off + 164);

            auto g_xyyz_yyzzz = cbuffer.data(gh_off + 165);

            auto g_xyyz_yzzzz = cbuffer.data(gh_off + 166);

            auto g_xyyz_zzzzz = cbuffer.data(gh_off + 167);

            auto g_xyzz_xxxxx = cbuffer.data(gh_off + 168);

            auto g_xyzz_xxxxy = cbuffer.data(gh_off + 169);

            auto g_xyzz_xxxxz = cbuffer.data(gh_off + 170);

            auto g_xyzz_xxxyy = cbuffer.data(gh_off + 171);

            auto g_xyzz_xxxyz = cbuffer.data(gh_off + 172);

            auto g_xyzz_xxxzz = cbuffer.data(gh_off + 173);

            auto g_xyzz_xxyyy = cbuffer.data(gh_off + 174);

            auto g_xyzz_xxyyz = cbuffer.data(gh_off + 175);

            auto g_xyzz_xxyzz = cbuffer.data(gh_off + 176);

            auto g_xyzz_xxzzz = cbuffer.data(gh_off + 177);

            auto g_xyzz_xyyyy = cbuffer.data(gh_off + 178);

            auto g_xyzz_xyyyz = cbuffer.data(gh_off + 179);

            auto g_xyzz_xyyzz = cbuffer.data(gh_off + 180);

            auto g_xyzz_xyzzz = cbuffer.data(gh_off + 181);

            auto g_xyzz_xzzzz = cbuffer.data(gh_off + 182);

            auto g_xyzz_yyyyy = cbuffer.data(gh_off + 183);

            auto g_xyzz_yyyyz = cbuffer.data(gh_off + 184);

            auto g_xyzz_yyyzz = cbuffer.data(gh_off + 185);

            auto g_xyzz_yyzzz = cbuffer.data(gh_off + 186);

            auto g_xyzz_yzzzz = cbuffer.data(gh_off + 187);

            auto g_xyzz_zzzzz = cbuffer.data(gh_off + 188);

            auto g_xzzz_xxxxx = cbuffer.data(gh_off + 189);

            auto g_xzzz_xxxxy = cbuffer.data(gh_off + 190);

            auto g_xzzz_xxxxz = cbuffer.data(gh_off + 191);

            auto g_xzzz_xxxyy = cbuffer.data(gh_off + 192);

            auto g_xzzz_xxxyz = cbuffer.data(gh_off + 193);

            auto g_xzzz_xxxzz = cbuffer.data(gh_off + 194);

            auto g_xzzz_xxyyy = cbuffer.data(gh_off + 195);

            auto g_xzzz_xxyyz = cbuffer.data(gh_off + 196);

            auto g_xzzz_xxyzz = cbuffer.data(gh_off + 197);

            auto g_xzzz_xxzzz = cbuffer.data(gh_off + 198);

            auto g_xzzz_xyyyy = cbuffer.data(gh_off + 199);

            auto g_xzzz_xyyyz = cbuffer.data(gh_off + 200);

            auto g_xzzz_xyyzz = cbuffer.data(gh_off + 201);

            auto g_xzzz_xyzzz = cbuffer.data(gh_off + 202);

            auto g_xzzz_xzzzz = cbuffer.data(gh_off + 203);

            auto g_xzzz_yyyyy = cbuffer.data(gh_off + 204);

            auto g_xzzz_yyyyz = cbuffer.data(gh_off + 205);

            auto g_xzzz_yyyzz = cbuffer.data(gh_off + 206);

            auto g_xzzz_yyzzz = cbuffer.data(gh_off + 207);

            auto g_xzzz_yzzzz = cbuffer.data(gh_off + 208);

            auto g_xzzz_zzzzz = cbuffer.data(gh_off + 209);

            auto g_yyyy_xxxxx = cbuffer.data(gh_off + 210);

            auto g_yyyy_xxxxy = cbuffer.data(gh_off + 211);

            auto g_yyyy_xxxxz = cbuffer.data(gh_off + 212);

            auto g_yyyy_xxxyy = cbuffer.data(gh_off + 213);

            auto g_yyyy_xxxyz = cbuffer.data(gh_off + 214);

            auto g_yyyy_xxxzz = cbuffer.data(gh_off + 215);

            auto g_yyyy_xxyyy = cbuffer.data(gh_off + 216);

            auto g_yyyy_xxyyz = cbuffer.data(gh_off + 217);

            auto g_yyyy_xxyzz = cbuffer.data(gh_off + 218);

            auto g_yyyy_xxzzz = cbuffer.data(gh_off + 219);

            auto g_yyyy_xyyyy = cbuffer.data(gh_off + 220);

            auto g_yyyy_xyyyz = cbuffer.data(gh_off + 221);

            auto g_yyyy_xyyzz = cbuffer.data(gh_off + 222);

            auto g_yyyy_xyzzz = cbuffer.data(gh_off + 223);

            auto g_yyyy_xzzzz = cbuffer.data(gh_off + 224);

            auto g_yyyy_yyyyy = cbuffer.data(gh_off + 225);

            auto g_yyyy_yyyyz = cbuffer.data(gh_off + 226);

            auto g_yyyy_yyyzz = cbuffer.data(gh_off + 227);

            auto g_yyyy_yyzzz = cbuffer.data(gh_off + 228);

            auto g_yyyy_yzzzz = cbuffer.data(gh_off + 229);

            auto g_yyyy_zzzzz = cbuffer.data(gh_off + 230);

            auto g_yyyz_xxxxx = cbuffer.data(gh_off + 231);

            auto g_yyyz_xxxxy = cbuffer.data(gh_off + 232);

            auto g_yyyz_xxxxz = cbuffer.data(gh_off + 233);

            auto g_yyyz_xxxyy = cbuffer.data(gh_off + 234);

            auto g_yyyz_xxxyz = cbuffer.data(gh_off + 235);

            auto g_yyyz_xxxzz = cbuffer.data(gh_off + 236);

            auto g_yyyz_xxyyy = cbuffer.data(gh_off + 237);

            auto g_yyyz_xxyyz = cbuffer.data(gh_off + 238);

            auto g_yyyz_xxyzz = cbuffer.data(gh_off + 239);

            auto g_yyyz_xxzzz = cbuffer.data(gh_off + 240);

            auto g_yyyz_xyyyy = cbuffer.data(gh_off + 241);

            auto g_yyyz_xyyyz = cbuffer.data(gh_off + 242);

            auto g_yyyz_xyyzz = cbuffer.data(gh_off + 243);

            auto g_yyyz_xyzzz = cbuffer.data(gh_off + 244);

            auto g_yyyz_xzzzz = cbuffer.data(gh_off + 245);

            auto g_yyyz_yyyyy = cbuffer.data(gh_off + 246);

            auto g_yyyz_yyyyz = cbuffer.data(gh_off + 247);

            auto g_yyyz_yyyzz = cbuffer.data(gh_off + 248);

            auto g_yyyz_yyzzz = cbuffer.data(gh_off + 249);

            auto g_yyyz_yzzzz = cbuffer.data(gh_off + 250);

            auto g_yyyz_zzzzz = cbuffer.data(gh_off + 251);

            auto g_yyzz_xxxxx = cbuffer.data(gh_off + 252);

            auto g_yyzz_xxxxy = cbuffer.data(gh_off + 253);

            auto g_yyzz_xxxxz = cbuffer.data(gh_off + 254);

            auto g_yyzz_xxxyy = cbuffer.data(gh_off + 255);

            auto g_yyzz_xxxyz = cbuffer.data(gh_off + 256);

            auto g_yyzz_xxxzz = cbuffer.data(gh_off + 257);

            auto g_yyzz_xxyyy = cbuffer.data(gh_off + 258);

            auto g_yyzz_xxyyz = cbuffer.data(gh_off + 259);

            auto g_yyzz_xxyzz = cbuffer.data(gh_off + 260);

            auto g_yyzz_xxzzz = cbuffer.data(gh_off + 261);

            auto g_yyzz_xyyyy = cbuffer.data(gh_off + 262);

            auto g_yyzz_xyyyz = cbuffer.data(gh_off + 263);

            auto g_yyzz_xyyzz = cbuffer.data(gh_off + 264);

            auto g_yyzz_xyzzz = cbuffer.data(gh_off + 265);

            auto g_yyzz_xzzzz = cbuffer.data(gh_off + 266);

            auto g_yyzz_yyyyy = cbuffer.data(gh_off + 267);

            auto g_yyzz_yyyyz = cbuffer.data(gh_off + 268);

            auto g_yyzz_yyyzz = cbuffer.data(gh_off + 269);

            auto g_yyzz_yyzzz = cbuffer.data(gh_off + 270);

            auto g_yyzz_yzzzz = cbuffer.data(gh_off + 271);

            auto g_yyzz_zzzzz = cbuffer.data(gh_off + 272);

            auto g_yzzz_xxxxx = cbuffer.data(gh_off + 273);

            auto g_yzzz_xxxxy = cbuffer.data(gh_off + 274);

            auto g_yzzz_xxxxz = cbuffer.data(gh_off + 275);

            auto g_yzzz_xxxyy = cbuffer.data(gh_off + 276);

            auto g_yzzz_xxxyz = cbuffer.data(gh_off + 277);

            auto g_yzzz_xxxzz = cbuffer.data(gh_off + 278);

            auto g_yzzz_xxyyy = cbuffer.data(gh_off + 279);

            auto g_yzzz_xxyyz = cbuffer.data(gh_off + 280);

            auto g_yzzz_xxyzz = cbuffer.data(gh_off + 281);

            auto g_yzzz_xxzzz = cbuffer.data(gh_off + 282);

            auto g_yzzz_xyyyy = cbuffer.data(gh_off + 283);

            auto g_yzzz_xyyyz = cbuffer.data(gh_off + 284);

            auto g_yzzz_xyyzz = cbuffer.data(gh_off + 285);

            auto g_yzzz_xyzzz = cbuffer.data(gh_off + 286);

            auto g_yzzz_xzzzz = cbuffer.data(gh_off + 287);

            auto g_yzzz_yyyyy = cbuffer.data(gh_off + 288);

            auto g_yzzz_yyyyz = cbuffer.data(gh_off + 289);

            auto g_yzzz_yyyzz = cbuffer.data(gh_off + 290);

            auto g_yzzz_yyzzz = cbuffer.data(gh_off + 291);

            auto g_yzzz_yzzzz = cbuffer.data(gh_off + 292);

            auto g_yzzz_zzzzz = cbuffer.data(gh_off + 293);

            auto g_zzzz_xxxxx = cbuffer.data(gh_off + 294);

            auto g_zzzz_xxxxy = cbuffer.data(gh_off + 295);

            auto g_zzzz_xxxxz = cbuffer.data(gh_off + 296);

            auto g_zzzz_xxxyy = cbuffer.data(gh_off + 297);

            auto g_zzzz_xxxyz = cbuffer.data(gh_off + 298);

            auto g_zzzz_xxxzz = cbuffer.data(gh_off + 299);

            auto g_zzzz_xxyyy = cbuffer.data(gh_off + 300);

            auto g_zzzz_xxyyz = cbuffer.data(gh_off + 301);

            auto g_zzzz_xxyzz = cbuffer.data(gh_off + 302);

            auto g_zzzz_xxzzz = cbuffer.data(gh_off + 303);

            auto g_zzzz_xyyyy = cbuffer.data(gh_off + 304);

            auto g_zzzz_xyyyz = cbuffer.data(gh_off + 305);

            auto g_zzzz_xyyzz = cbuffer.data(gh_off + 306);

            auto g_zzzz_xyzzz = cbuffer.data(gh_off + 307);

            auto g_zzzz_xzzzz = cbuffer.data(gh_off + 308);

            auto g_zzzz_yyyyy = cbuffer.data(gh_off + 309);

            auto g_zzzz_yyyyz = cbuffer.data(gh_off + 310);

            auto g_zzzz_yyyzz = cbuffer.data(gh_off + 311);

            auto g_zzzz_yyzzz = cbuffer.data(gh_off + 312);

            auto g_zzzz_yzzzz = cbuffer.data(gh_off + 313);

            auto g_zzzz_zzzzz = cbuffer.data(gh_off + 314);

            /// Set up components of auxilary buffer : SSGH

            const auto gh_geom_10_off = idx_geom_10_xxgh + (i * bcomps + j) * 315;

            auto g_x_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_y_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 5);

            auto g_y_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 6);

            auto g_y_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 7);

            auto g_y_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 8);

            auto g_y_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 9);

            auto g_y_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 10);

            auto g_y_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 11);

            auto g_y_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 12);

            auto g_y_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 13);

            auto g_y_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 14);

            auto g_y_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 15);

            auto g_y_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 16);

            auto g_y_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 17);

            auto g_y_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 18);

            auto g_y_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 19);

            auto g_y_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 20);

            auto g_y_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 21);

            auto g_y_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 22);

            auto g_y_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 23);

            auto g_y_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 24);

            auto g_y_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 25);

            auto g_y_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 26);

            auto g_y_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 27);

            auto g_y_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 28);

            auto g_y_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 29);

            auto g_y_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 30);

            auto g_y_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 31);

            auto g_y_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 32);

            auto g_y_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 33);

            auto g_y_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 34);

            auto g_y_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 35);

            auto g_y_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 36);

            auto g_y_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 37);

            auto g_y_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 38);

            auto g_y_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 39);

            auto g_y_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 40);

            auto g_y_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 41);

            auto g_y_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 42);

            auto g_y_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 43);

            auto g_y_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 44);

            auto g_y_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 45);

            auto g_y_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 46);

            auto g_y_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 47);

            auto g_y_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 48);

            auto g_y_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 49);

            auto g_y_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 50);

            auto g_y_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 51);

            auto g_y_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 52);

            auto g_y_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 53);

            auto g_y_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 54);

            auto g_y_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 55);

            auto g_y_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 56);

            auto g_y_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 57);

            auto g_y_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 58);

            auto g_y_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 59);

            auto g_y_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 60);

            auto g_y_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 61);

            auto g_y_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 62);

            auto g_y_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 63);

            auto g_y_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 64);

            auto g_y_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 65);

            auto g_y_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 66);

            auto g_y_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 67);

            auto g_y_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 68);

            auto g_y_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 69);

            auto g_y_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 70);

            auto g_y_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 71);

            auto g_y_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 72);

            auto g_y_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 73);

            auto g_y_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 74);

            auto g_y_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 75);

            auto g_y_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 76);

            auto g_y_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 77);

            auto g_y_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 78);

            auto g_y_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 79);

            auto g_y_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 80);

            auto g_y_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 81);

            auto g_y_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 82);

            auto g_y_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 83);

            auto g_y_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 84);

            auto g_y_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 85);

            auto g_y_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 86);

            auto g_y_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 87);

            auto g_y_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 88);

            auto g_y_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 89);

            auto g_y_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 90);

            auto g_y_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 91);

            auto g_y_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 92);

            auto g_y_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 93);

            auto g_y_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 94);

            auto g_y_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 95);

            auto g_y_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 96);

            auto g_y_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 97);

            auto g_y_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 98);

            auto g_y_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 99);

            auto g_y_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 100);

            auto g_y_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 101);

            auto g_y_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 102);

            auto g_y_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 103);

            auto g_y_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 104);

            auto g_y_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 105);

            auto g_y_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 106);

            auto g_y_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 107);

            auto g_y_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 108);

            auto g_y_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 109);

            auto g_y_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 110);

            auto g_y_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 111);

            auto g_y_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 112);

            auto g_y_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 113);

            auto g_y_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 114);

            auto g_y_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 115);

            auto g_y_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 116);

            auto g_y_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 117);

            auto g_y_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 118);

            auto g_y_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 119);

            auto g_y_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 120);

            auto g_y_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 121);

            auto g_y_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 122);

            auto g_y_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 123);

            auto g_y_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 124);

            auto g_y_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 125);

            auto g_y_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 126);

            auto g_y_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 127);

            auto g_y_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 128);

            auto g_y_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 129);

            auto g_y_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 130);

            auto g_y_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 131);

            auto g_y_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 132);

            auto g_y_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 133);

            auto g_y_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 134);

            auto g_y_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 135);

            auto g_y_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 136);

            auto g_y_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 137);

            auto g_y_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 138);

            auto g_y_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 139);

            auto g_y_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 140);

            auto g_y_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 141);

            auto g_y_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 142);

            auto g_y_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 143);

            auto g_y_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 144);

            auto g_y_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 145);

            auto g_y_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 146);

            auto g_y_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 147);

            auto g_y_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 148);

            auto g_y_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 149);

            auto g_y_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 150);

            auto g_y_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 151);

            auto g_y_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 152);

            auto g_y_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 153);

            auto g_y_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 154);

            auto g_y_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 155);

            auto g_y_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 156);

            auto g_y_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 157);

            auto g_y_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 158);

            auto g_y_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 159);

            auto g_y_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 160);

            auto g_y_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 161);

            auto g_y_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 162);

            auto g_y_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 163);

            auto g_y_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 164);

            auto g_y_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 165);

            auto g_y_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 166);

            auto g_y_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 167);

            auto g_y_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 168);

            auto g_y_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 169);

            auto g_y_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 170);

            auto g_y_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 171);

            auto g_y_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 172);

            auto g_y_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 173);

            auto g_y_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 174);

            auto g_y_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 175);

            auto g_y_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 176);

            auto g_y_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 177);

            auto g_y_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 178);

            auto g_y_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 179);

            auto g_y_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 180);

            auto g_y_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 181);

            auto g_y_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 182);

            auto g_y_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 183);

            auto g_y_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 184);

            auto g_y_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 185);

            auto g_y_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 186);

            auto g_y_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 187);

            auto g_y_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 188);

            auto g_y_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 189);

            auto g_y_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 190);

            auto g_y_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 191);

            auto g_y_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 192);

            auto g_y_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 193);

            auto g_y_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 194);

            auto g_y_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 195);

            auto g_y_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 196);

            auto g_y_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 197);

            auto g_y_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 198);

            auto g_y_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 199);

            auto g_y_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 200);

            auto g_y_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 201);

            auto g_y_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 202);

            auto g_y_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 203);

            auto g_y_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 204);

            auto g_y_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 205);

            auto g_y_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 206);

            auto g_y_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 207);

            auto g_y_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 208);

            auto g_y_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 209);

            auto g_y_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 210);

            auto g_y_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 211);

            auto g_y_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 212);

            auto g_y_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 213);

            auto g_y_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 214);

            auto g_y_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 215);

            auto g_y_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 216);

            auto g_y_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 217);

            auto g_y_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 218);

            auto g_y_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 219);

            auto g_y_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 220);

            auto g_y_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 221);

            auto g_y_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 222);

            auto g_y_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 223);

            auto g_y_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 224);

            auto g_y_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 225);

            auto g_y_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 226);

            auto g_y_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 227);

            auto g_y_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 228);

            auto g_y_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 229);

            auto g_y_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 230);

            auto g_y_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 231);

            auto g_y_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 232);

            auto g_y_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 233);

            auto g_y_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 234);

            auto g_y_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 235);

            auto g_y_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 236);

            auto g_y_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 237);

            auto g_y_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 238);

            auto g_y_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 239);

            auto g_y_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 240);

            auto g_y_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 241);

            auto g_y_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 242);

            auto g_y_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 243);

            auto g_y_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 244);

            auto g_y_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 245);

            auto g_y_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 246);

            auto g_y_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 247);

            auto g_y_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 248);

            auto g_y_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 249);

            auto g_y_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 250);

            auto g_y_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 251);

            auto g_y_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 252);

            auto g_y_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 253);

            auto g_y_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 254);

            auto g_y_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 255);

            auto g_y_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 256);

            auto g_y_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 257);

            auto g_y_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 258);

            auto g_y_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 259);

            auto g_y_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 260);

            auto g_y_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 261);

            auto g_y_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 262);

            auto g_y_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 263);

            auto g_y_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 264);

            auto g_y_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 265);

            auto g_y_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 266);

            auto g_y_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 267);

            auto g_y_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 268);

            auto g_y_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 269);

            auto g_y_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 270);

            auto g_y_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 271);

            auto g_y_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 272);

            auto g_y_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 273);

            auto g_y_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 274);

            auto g_y_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 275);

            auto g_y_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 276);

            auto g_y_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 277);

            auto g_y_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 278);

            auto g_y_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 279);

            auto g_y_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 280);

            auto g_y_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 281);

            auto g_y_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 282);

            auto g_y_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 283);

            auto g_y_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 284);

            auto g_y_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 285);

            auto g_y_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 286);

            auto g_y_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 287);

            auto g_y_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 288);

            auto g_y_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 289);

            auto g_y_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 290);

            auto g_y_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 291);

            auto g_y_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 292);

            auto g_y_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 293);

            auto g_y_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 294);

            auto g_y_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 295);

            auto g_y_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 296);

            auto g_y_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 297);

            auto g_y_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 298);

            auto g_y_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 299);

            auto g_y_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 300);

            auto g_y_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 301);

            auto g_y_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 302);

            auto g_y_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 303);

            auto g_y_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 304);

            auto g_y_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 305);

            auto g_y_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 306);

            auto g_y_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 307);

            auto g_y_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 308);

            auto g_y_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 309);

            auto g_y_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 310);

            auto g_y_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 311);

            auto g_y_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 312);

            auto g_y_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 313);

            auto g_y_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 314);

            auto g_z_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 5);

            auto g_z_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 6);

            auto g_z_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 7);

            auto g_z_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 8);

            auto g_z_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 9);

            auto g_z_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 10);

            auto g_z_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 11);

            auto g_z_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 12);

            auto g_z_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 13);

            auto g_z_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 14);

            auto g_z_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 15);

            auto g_z_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 16);

            auto g_z_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 17);

            auto g_z_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 18);

            auto g_z_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 19);

            auto g_z_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 20);

            auto g_z_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 21);

            auto g_z_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 22);

            auto g_z_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 23);

            auto g_z_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 24);

            auto g_z_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 25);

            auto g_z_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 26);

            auto g_z_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 27);

            auto g_z_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 28);

            auto g_z_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 29);

            auto g_z_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 30);

            auto g_z_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 31);

            auto g_z_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 32);

            auto g_z_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 33);

            auto g_z_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 34);

            auto g_z_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 35);

            auto g_z_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 36);

            auto g_z_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 37);

            auto g_z_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 38);

            auto g_z_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 39);

            auto g_z_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 40);

            auto g_z_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 41);

            auto g_z_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 42);

            auto g_z_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 43);

            auto g_z_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 44);

            auto g_z_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 45);

            auto g_z_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 46);

            auto g_z_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 47);

            auto g_z_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 48);

            auto g_z_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 49);

            auto g_z_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 50);

            auto g_z_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 51);

            auto g_z_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 52);

            auto g_z_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 53);

            auto g_z_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 54);

            auto g_z_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 55);

            auto g_z_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 56);

            auto g_z_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 57);

            auto g_z_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 58);

            auto g_z_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 59);

            auto g_z_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 60);

            auto g_z_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 61);

            auto g_z_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 62);

            auto g_z_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 63);

            auto g_z_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 64);

            auto g_z_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 65);

            auto g_z_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 66);

            auto g_z_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 67);

            auto g_z_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 68);

            auto g_z_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 69);

            auto g_z_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 70);

            auto g_z_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 71);

            auto g_z_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 72);

            auto g_z_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 73);

            auto g_z_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 74);

            auto g_z_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 75);

            auto g_z_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 76);

            auto g_z_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 77);

            auto g_z_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 78);

            auto g_z_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 79);

            auto g_z_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 80);

            auto g_z_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 81);

            auto g_z_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 82);

            auto g_z_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 83);

            auto g_z_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 84);

            auto g_z_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 85);

            auto g_z_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 86);

            auto g_z_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 87);

            auto g_z_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 88);

            auto g_z_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 89);

            auto g_z_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 90);

            auto g_z_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 91);

            auto g_z_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 92);

            auto g_z_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 93);

            auto g_z_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 94);

            auto g_z_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 95);

            auto g_z_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 96);

            auto g_z_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 97);

            auto g_z_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 98);

            auto g_z_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 99);

            auto g_z_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 100);

            auto g_z_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 101);

            auto g_z_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 102);

            auto g_z_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 103);

            auto g_z_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 104);

            auto g_z_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 105);

            auto g_z_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 106);

            auto g_z_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 107);

            auto g_z_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 108);

            auto g_z_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 109);

            auto g_z_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 110);

            auto g_z_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 111);

            auto g_z_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 112);

            auto g_z_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 113);

            auto g_z_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 114);

            auto g_z_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 115);

            auto g_z_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 116);

            auto g_z_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 117);

            auto g_z_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 118);

            auto g_z_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 119);

            auto g_z_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 120);

            auto g_z_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 121);

            auto g_z_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 122);

            auto g_z_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 123);

            auto g_z_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 124);

            auto g_z_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 125);

            auto g_z_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 126);

            auto g_z_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 127);

            auto g_z_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 128);

            auto g_z_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 129);

            auto g_z_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 130);

            auto g_z_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 131);

            auto g_z_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 132);

            auto g_z_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 133);

            auto g_z_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 134);

            auto g_z_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 135);

            auto g_z_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 136);

            auto g_z_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 137);

            auto g_z_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 138);

            auto g_z_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 139);

            auto g_z_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 140);

            auto g_z_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 141);

            auto g_z_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 142);

            auto g_z_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 143);

            auto g_z_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 144);

            auto g_z_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 145);

            auto g_z_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 146);

            auto g_z_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 147);

            auto g_z_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 148);

            auto g_z_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 149);

            auto g_z_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 150);

            auto g_z_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 151);

            auto g_z_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 152);

            auto g_z_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 153);

            auto g_z_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 154);

            auto g_z_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 155);

            auto g_z_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 156);

            auto g_z_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 157);

            auto g_z_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 158);

            auto g_z_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 159);

            auto g_z_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 160);

            auto g_z_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 161);

            auto g_z_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 162);

            auto g_z_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 163);

            auto g_z_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 164);

            auto g_z_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 165);

            auto g_z_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 166);

            auto g_z_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 167);

            auto g_z_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 168);

            auto g_z_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 169);

            auto g_z_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 170);

            auto g_z_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 171);

            auto g_z_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 172);

            auto g_z_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 173);

            auto g_z_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 174);

            auto g_z_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 175);

            auto g_z_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 176);

            auto g_z_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 177);

            auto g_z_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 178);

            auto g_z_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 179);

            auto g_z_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 180);

            auto g_z_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 181);

            auto g_z_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 182);

            auto g_z_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 183);

            auto g_z_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 184);

            auto g_z_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 185);

            auto g_z_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 186);

            auto g_z_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 187);

            auto g_z_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 188);

            auto g_z_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 189);

            auto g_z_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 190);

            auto g_z_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 191);

            auto g_z_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 192);

            auto g_z_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 193);

            auto g_z_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 194);

            auto g_z_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 195);

            auto g_z_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 196);

            auto g_z_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 197);

            auto g_z_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 198);

            auto g_z_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 199);

            auto g_z_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 200);

            auto g_z_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 201);

            auto g_z_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 202);

            auto g_z_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 203);

            auto g_z_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 204);

            auto g_z_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 205);

            auto g_z_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 206);

            auto g_z_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 207);

            auto g_z_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 208);

            auto g_z_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 209);

            auto g_z_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 210);

            auto g_z_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 211);

            auto g_z_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 212);

            auto g_z_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 213);

            auto g_z_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 214);

            auto g_z_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 215);

            auto g_z_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 216);

            auto g_z_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 217);

            auto g_z_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 218);

            auto g_z_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 219);

            auto g_z_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 220);

            auto g_z_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 221);

            auto g_z_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 222);

            auto g_z_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 223);

            auto g_z_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 224);

            auto g_z_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 225);

            auto g_z_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 226);

            auto g_z_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 227);

            auto g_z_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 228);

            auto g_z_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 229);

            auto g_z_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 230);

            auto g_z_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 231);

            auto g_z_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 232);

            auto g_z_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 233);

            auto g_z_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 234);

            auto g_z_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 235);

            auto g_z_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 236);

            auto g_z_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 237);

            auto g_z_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 238);

            auto g_z_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 239);

            auto g_z_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 240);

            auto g_z_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 241);

            auto g_z_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 242);

            auto g_z_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 243);

            auto g_z_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 244);

            auto g_z_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 245);

            auto g_z_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 246);

            auto g_z_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 247);

            auto g_z_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 248);

            auto g_z_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 249);

            auto g_z_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 250);

            auto g_z_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 251);

            auto g_z_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 252);

            auto g_z_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 253);

            auto g_z_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 254);

            auto g_z_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 255);

            auto g_z_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 256);

            auto g_z_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 257);

            auto g_z_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 258);

            auto g_z_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 259);

            auto g_z_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 260);

            auto g_z_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 261);

            auto g_z_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 262);

            auto g_z_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 263);

            auto g_z_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 264);

            auto g_z_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 265);

            auto g_z_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 266);

            auto g_z_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 267);

            auto g_z_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 268);

            auto g_z_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 269);

            auto g_z_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 270);

            auto g_z_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 271);

            auto g_z_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 272);

            auto g_z_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 273);

            auto g_z_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 274);

            auto g_z_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 275);

            auto g_z_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 276);

            auto g_z_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 277);

            auto g_z_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 278);

            auto g_z_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 279);

            auto g_z_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 280);

            auto g_z_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 281);

            auto g_z_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 282);

            auto g_z_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 283);

            auto g_z_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 284);

            auto g_z_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 285);

            auto g_z_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 286);

            auto g_z_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 287);

            auto g_z_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 288);

            auto g_z_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 289);

            auto g_z_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 290);

            auto g_z_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 291);

            auto g_z_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 292);

            auto g_z_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 293);

            auto g_z_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 294);

            auto g_z_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 295);

            auto g_z_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 296);

            auto g_z_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 297);

            auto g_z_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 298);

            auto g_z_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 299);

            auto g_z_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 300);

            auto g_z_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 301);

            auto g_z_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 302);

            auto g_z_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 303);

            auto g_z_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 304);

            auto g_z_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 305);

            auto g_z_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 306);

            auto g_z_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 307);

            auto g_z_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 308);

            auto g_z_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 309);

            auto g_z_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 310);

            auto g_z_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 311);

            auto g_z_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 312);

            auto g_z_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 313);

            auto g_z_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 314);

            /// Set up components of auxilary buffer : SSGI

            const auto gi_geom_10_off = idx_geom_10_xxgi + (i * bcomps + j) * 420;

            auto g_x_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_x_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 335);

            auto g_x_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 363);

            auto g_x_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 377);

            auto g_x_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 391);

            auto g_x_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 398);

            auto g_x_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 419);

            auto g_y_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 5);

            auto g_y_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 6);

            auto g_y_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 7);

            auto g_y_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 8);

            auto g_y_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 9);

            auto g_y_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 10);

            auto g_y_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 11);

            auto g_y_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 12);

            auto g_y_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 13);

            auto g_y_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 14);

            auto g_y_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 15);

            auto g_y_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 16);

            auto g_y_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 17);

            auto g_y_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 18);

            auto g_y_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 19);

            auto g_y_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 20);

            auto g_y_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 21);

            auto g_y_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 22);

            auto g_y_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 23);

            auto g_y_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 24);

            auto g_y_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 25);

            auto g_y_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 26);

            auto g_y_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 27);

            auto g_y_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 28);

            auto g_y_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 29);

            auto g_y_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 30);

            auto g_y_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 31);

            auto g_y_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 32);

            auto g_y_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 33);

            auto g_y_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 34);

            auto g_y_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 35);

            auto g_y_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 36);

            auto g_y_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 37);

            auto g_y_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 38);

            auto g_y_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 39);

            auto g_y_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 40);

            auto g_y_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 41);

            auto g_y_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 42);

            auto g_y_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 43);

            auto g_y_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 44);

            auto g_y_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 45);

            auto g_y_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 46);

            auto g_y_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 47);

            auto g_y_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 48);

            auto g_y_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 49);

            auto g_y_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 50);

            auto g_y_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 51);

            auto g_y_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 52);

            auto g_y_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 53);

            auto g_y_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 54);

            auto g_y_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 55);

            auto g_y_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 56);

            auto g_y_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 57);

            auto g_y_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 58);

            auto g_y_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 59);

            auto g_y_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 60);

            auto g_y_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 61);

            auto g_y_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 62);

            auto g_y_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 63);

            auto g_y_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 64);

            auto g_y_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 65);

            auto g_y_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 66);

            auto g_y_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 67);

            auto g_y_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 68);

            auto g_y_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 69);

            auto g_y_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 70);

            auto g_y_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 71);

            auto g_y_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 72);

            auto g_y_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 73);

            auto g_y_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 74);

            auto g_y_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 75);

            auto g_y_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 76);

            auto g_y_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 77);

            auto g_y_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 78);

            auto g_y_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 79);

            auto g_y_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 80);

            auto g_y_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 81);

            auto g_y_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 82);

            auto g_y_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 83);

            auto g_y_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 84);

            auto g_y_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 85);

            auto g_y_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 86);

            auto g_y_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 87);

            auto g_y_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 88);

            auto g_y_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 89);

            auto g_y_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 90);

            auto g_y_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 91);

            auto g_y_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 92);

            auto g_y_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 93);

            auto g_y_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 94);

            auto g_y_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 95);

            auto g_y_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 96);

            auto g_y_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 97);

            auto g_y_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 98);

            auto g_y_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 99);

            auto g_y_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 100);

            auto g_y_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 101);

            auto g_y_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 102);

            auto g_y_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 103);

            auto g_y_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 104);

            auto g_y_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 105);

            auto g_y_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 106);

            auto g_y_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 107);

            auto g_y_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 108);

            auto g_y_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 109);

            auto g_y_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 110);

            auto g_y_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 111);

            auto g_y_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 112);

            auto g_y_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 113);

            auto g_y_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 114);

            auto g_y_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 115);

            auto g_y_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 116);

            auto g_y_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 117);

            auto g_y_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 118);

            auto g_y_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 119);

            auto g_y_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 120);

            auto g_y_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 121);

            auto g_y_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 122);

            auto g_y_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 123);

            auto g_y_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 124);

            auto g_y_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 125);

            auto g_y_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 126);

            auto g_y_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 127);

            auto g_y_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 128);

            auto g_y_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 129);

            auto g_y_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 130);

            auto g_y_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 131);

            auto g_y_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 132);

            auto g_y_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 133);

            auto g_y_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 134);

            auto g_y_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 135);

            auto g_y_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 136);

            auto g_y_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 137);

            auto g_y_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 138);

            auto g_y_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 139);

            auto g_y_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 140);

            auto g_y_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 141);

            auto g_y_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 142);

            auto g_y_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 143);

            auto g_y_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 144);

            auto g_y_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 145);

            auto g_y_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 146);

            auto g_y_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 147);

            auto g_y_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 148);

            auto g_y_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 149);

            auto g_y_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 150);

            auto g_y_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 151);

            auto g_y_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 152);

            auto g_y_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 153);

            auto g_y_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 154);

            auto g_y_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 155);

            auto g_y_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 156);

            auto g_y_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 157);

            auto g_y_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 158);

            auto g_y_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 159);

            auto g_y_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 160);

            auto g_y_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 161);

            auto g_y_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 162);

            auto g_y_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 163);

            auto g_y_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 164);

            auto g_y_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 165);

            auto g_y_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 166);

            auto g_y_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 167);

            auto g_y_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 168);

            auto g_y_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 169);

            auto g_y_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 170);

            auto g_y_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 171);

            auto g_y_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 172);

            auto g_y_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 173);

            auto g_y_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 174);

            auto g_y_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 175);

            auto g_y_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 176);

            auto g_y_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 177);

            auto g_y_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 178);

            auto g_y_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 179);

            auto g_y_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 180);

            auto g_y_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 181);

            auto g_y_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 182);

            auto g_y_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 183);

            auto g_y_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 184);

            auto g_y_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 185);

            auto g_y_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 186);

            auto g_y_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 187);

            auto g_y_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 188);

            auto g_y_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 189);

            auto g_y_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 190);

            auto g_y_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 191);

            auto g_y_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 192);

            auto g_y_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 193);

            auto g_y_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 194);

            auto g_y_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 195);

            auto g_y_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 196);

            auto g_y_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 197);

            auto g_y_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 198);

            auto g_y_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 199);

            auto g_y_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 200);

            auto g_y_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 201);

            auto g_y_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 202);

            auto g_y_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 203);

            auto g_y_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 204);

            auto g_y_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 205);

            auto g_y_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 206);

            auto g_y_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 207);

            auto g_y_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 208);

            auto g_y_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 209);

            auto g_y_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 210);

            auto g_y_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 211);

            auto g_y_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 212);

            auto g_y_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 213);

            auto g_y_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 214);

            auto g_y_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 215);

            auto g_y_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 216);

            auto g_y_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 217);

            auto g_y_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 218);

            auto g_y_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 219);

            auto g_y_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 220);

            auto g_y_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 221);

            auto g_y_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 222);

            auto g_y_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 223);

            auto g_y_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 224);

            auto g_y_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 225);

            auto g_y_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 226);

            auto g_y_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 227);

            auto g_y_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 228);

            auto g_y_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 229);

            auto g_y_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 230);

            auto g_y_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 231);

            auto g_y_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 232);

            auto g_y_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 233);

            auto g_y_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 234);

            auto g_y_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 235);

            auto g_y_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 236);

            auto g_y_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 237);

            auto g_y_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 238);

            auto g_y_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 239);

            auto g_y_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 240);

            auto g_y_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 241);

            auto g_y_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 242);

            auto g_y_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 243);

            auto g_y_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 244);

            auto g_y_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 245);

            auto g_y_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 246);

            auto g_y_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 247);

            auto g_y_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 248);

            auto g_y_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 249);

            auto g_y_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 250);

            auto g_y_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 251);

            auto g_y_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 252);

            auto g_y_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 253);

            auto g_y_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 254);

            auto g_y_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 255);

            auto g_y_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 256);

            auto g_y_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 257);

            auto g_y_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 258);

            auto g_y_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 259);

            auto g_y_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 260);

            auto g_y_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 261);

            auto g_y_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 262);

            auto g_y_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 263);

            auto g_y_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 264);

            auto g_y_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 265);

            auto g_y_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 266);

            auto g_y_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 267);

            auto g_y_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 268);

            auto g_y_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 269);

            auto g_y_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 270);

            auto g_y_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 271);

            auto g_y_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 272);

            auto g_y_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 273);

            auto g_y_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 274);

            auto g_y_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 275);

            auto g_y_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 276);

            auto g_y_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 277);

            auto g_y_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 278);

            auto g_y_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 279);

            auto g_y_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 280);

            auto g_y_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 281);

            auto g_y_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 282);

            auto g_y_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 283);

            auto g_y_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 284);

            auto g_y_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 285);

            auto g_y_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 286);

            auto g_y_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 287);

            auto g_y_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 288);

            auto g_y_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 289);

            auto g_y_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 290);

            auto g_y_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 291);

            auto g_y_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 292);

            auto g_y_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 293);

            auto g_y_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 294);

            auto g_y_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 295);

            auto g_y_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 296);

            auto g_y_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 297);

            auto g_y_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 298);

            auto g_y_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 299);

            auto g_y_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 300);

            auto g_y_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 301);

            auto g_y_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 302);

            auto g_y_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 303);

            auto g_y_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 304);

            auto g_y_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 305);

            auto g_y_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 306);

            auto g_y_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 307);

            auto g_y_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 308);

            auto g_y_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 309);

            auto g_y_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 310);

            auto g_y_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 311);

            auto g_y_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 312);

            auto g_y_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 313);

            auto g_y_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 314);

            auto g_y_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 315);

            auto g_y_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 316);

            auto g_y_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 317);

            auto g_y_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 318);

            auto g_y_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 319);

            auto g_y_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 320);

            auto g_y_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 321);

            auto g_y_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 322);

            auto g_y_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 323);

            auto g_y_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 324);

            auto g_y_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 325);

            auto g_y_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 326);

            auto g_y_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 327);

            auto g_y_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 328);

            auto g_y_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 329);

            auto g_y_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 330);

            auto g_y_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 331);

            auto g_y_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 332);

            auto g_y_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 333);

            auto g_y_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 334);

            auto g_y_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 335);

            auto g_y_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 336);

            auto g_y_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 337);

            auto g_y_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 338);

            auto g_y_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 339);

            auto g_y_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 340);

            auto g_y_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 341);

            auto g_y_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 342);

            auto g_y_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 343);

            auto g_y_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 344);

            auto g_y_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 345);

            auto g_y_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 346);

            auto g_y_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 347);

            auto g_y_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 348);

            auto g_y_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 349);

            auto g_y_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 350);

            auto g_y_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 351);

            auto g_y_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 352);

            auto g_y_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 353);

            auto g_y_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 354);

            auto g_y_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 355);

            auto g_y_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 356);

            auto g_y_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 357);

            auto g_y_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 358);

            auto g_y_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 359);

            auto g_y_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 360);

            auto g_y_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 361);

            auto g_y_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 362);

            auto g_y_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 363);

            auto g_y_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 364);

            auto g_y_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 365);

            auto g_y_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 366);

            auto g_y_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 367);

            auto g_y_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 368);

            auto g_y_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 369);

            auto g_y_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 370);

            auto g_y_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 371);

            auto g_y_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 372);

            auto g_y_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 373);

            auto g_y_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 374);

            auto g_y_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 375);

            auto g_y_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 376);

            auto g_y_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 377);

            auto g_y_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 378);

            auto g_y_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 379);

            auto g_y_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 380);

            auto g_y_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 381);

            auto g_y_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 382);

            auto g_y_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 383);

            auto g_y_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 384);

            auto g_y_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 385);

            auto g_y_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 386);

            auto g_y_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 387);

            auto g_y_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 388);

            auto g_y_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 389);

            auto g_y_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 390);

            auto g_y_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 391);

            auto g_y_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 392);

            auto g_y_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 393);

            auto g_y_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 394);

            auto g_y_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 395);

            auto g_y_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 396);

            auto g_y_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 397);

            auto g_y_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 398);

            auto g_y_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 399);

            auto g_y_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 400);

            auto g_y_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 401);

            auto g_y_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 402);

            auto g_y_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 403);

            auto g_y_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 404);

            auto g_y_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 405);

            auto g_y_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 406);

            auto g_y_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 407);

            auto g_y_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 408);

            auto g_y_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 409);

            auto g_y_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 410);

            auto g_y_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 411);

            auto g_y_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 412);

            auto g_y_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 413);

            auto g_y_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 414);

            auto g_y_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 415);

            auto g_y_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 416);

            auto g_y_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 417);

            auto g_y_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 418);

            auto g_y_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 419);

            auto g_z_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 5);

            auto g_z_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 6);

            auto g_z_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 7);

            auto g_z_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 8);

            auto g_z_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 9);

            auto g_z_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 10);

            auto g_z_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 11);

            auto g_z_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 12);

            auto g_z_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 13);

            auto g_z_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 14);

            auto g_z_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 15);

            auto g_z_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 16);

            auto g_z_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 17);

            auto g_z_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 18);

            auto g_z_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 19);

            auto g_z_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 20);

            auto g_z_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 21);

            auto g_z_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 22);

            auto g_z_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 23);

            auto g_z_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 24);

            auto g_z_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 25);

            auto g_z_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 26);

            auto g_z_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 27);

            auto g_z_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 28);

            auto g_z_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 29);

            auto g_z_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 30);

            auto g_z_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 31);

            auto g_z_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 32);

            auto g_z_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 33);

            auto g_z_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 34);

            auto g_z_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 35);

            auto g_z_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 36);

            auto g_z_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 37);

            auto g_z_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 38);

            auto g_z_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 39);

            auto g_z_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 40);

            auto g_z_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 41);

            auto g_z_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 42);

            auto g_z_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 43);

            auto g_z_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 44);

            auto g_z_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 45);

            auto g_z_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 46);

            auto g_z_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 47);

            auto g_z_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 48);

            auto g_z_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 49);

            auto g_z_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 50);

            auto g_z_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 51);

            auto g_z_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 52);

            auto g_z_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 53);

            auto g_z_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 54);

            auto g_z_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 55);

            auto g_z_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 56);

            auto g_z_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 57);

            auto g_z_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 58);

            auto g_z_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 59);

            auto g_z_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 60);

            auto g_z_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 61);

            auto g_z_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 62);

            auto g_z_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 63);

            auto g_z_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 64);

            auto g_z_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 65);

            auto g_z_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 66);

            auto g_z_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 67);

            auto g_z_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 68);

            auto g_z_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 69);

            auto g_z_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 70);

            auto g_z_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 71);

            auto g_z_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 72);

            auto g_z_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 73);

            auto g_z_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 74);

            auto g_z_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 75);

            auto g_z_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 76);

            auto g_z_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 77);

            auto g_z_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 78);

            auto g_z_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 79);

            auto g_z_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 80);

            auto g_z_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 81);

            auto g_z_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 82);

            auto g_z_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 83);

            auto g_z_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 84);

            auto g_z_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 85);

            auto g_z_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 86);

            auto g_z_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 87);

            auto g_z_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 88);

            auto g_z_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 89);

            auto g_z_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 90);

            auto g_z_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 91);

            auto g_z_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 92);

            auto g_z_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 93);

            auto g_z_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 94);

            auto g_z_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 95);

            auto g_z_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 96);

            auto g_z_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 97);

            auto g_z_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 98);

            auto g_z_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 99);

            auto g_z_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 100);

            auto g_z_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 101);

            auto g_z_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 102);

            auto g_z_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 103);

            auto g_z_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 104);

            auto g_z_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 105);

            auto g_z_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 106);

            auto g_z_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 107);

            auto g_z_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 108);

            auto g_z_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 109);

            auto g_z_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 110);

            auto g_z_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 111);

            auto g_z_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 112);

            auto g_z_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 113);

            auto g_z_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 114);

            auto g_z_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 115);

            auto g_z_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 116);

            auto g_z_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 117);

            auto g_z_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 118);

            auto g_z_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 119);

            auto g_z_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 120);

            auto g_z_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 121);

            auto g_z_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 122);

            auto g_z_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 123);

            auto g_z_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 124);

            auto g_z_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 125);

            auto g_z_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 126);

            auto g_z_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 127);

            auto g_z_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 128);

            auto g_z_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 129);

            auto g_z_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 130);

            auto g_z_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 131);

            auto g_z_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 132);

            auto g_z_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 133);

            auto g_z_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 134);

            auto g_z_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 135);

            auto g_z_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 136);

            auto g_z_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 137);

            auto g_z_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 138);

            auto g_z_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 139);

            auto g_z_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 140);

            auto g_z_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 141);

            auto g_z_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 142);

            auto g_z_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 143);

            auto g_z_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 144);

            auto g_z_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 145);

            auto g_z_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 146);

            auto g_z_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 147);

            auto g_z_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 148);

            auto g_z_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 149);

            auto g_z_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 150);

            auto g_z_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 151);

            auto g_z_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 152);

            auto g_z_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 153);

            auto g_z_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 154);

            auto g_z_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 155);

            auto g_z_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 156);

            auto g_z_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 157);

            auto g_z_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 158);

            auto g_z_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 159);

            auto g_z_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 160);

            auto g_z_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 161);

            auto g_z_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 162);

            auto g_z_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 163);

            auto g_z_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 164);

            auto g_z_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 165);

            auto g_z_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 166);

            auto g_z_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 167);

            auto g_z_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 168);

            auto g_z_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 169);

            auto g_z_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 170);

            auto g_z_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 171);

            auto g_z_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 172);

            auto g_z_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 173);

            auto g_z_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 174);

            auto g_z_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 175);

            auto g_z_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 176);

            auto g_z_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 177);

            auto g_z_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 178);

            auto g_z_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 179);

            auto g_z_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 180);

            auto g_z_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 181);

            auto g_z_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 182);

            auto g_z_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 183);

            auto g_z_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 184);

            auto g_z_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 185);

            auto g_z_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 186);

            auto g_z_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 187);

            auto g_z_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 188);

            auto g_z_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 189);

            auto g_z_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 190);

            auto g_z_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 191);

            auto g_z_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 192);

            auto g_z_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 193);

            auto g_z_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 194);

            auto g_z_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 195);

            auto g_z_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 196);

            auto g_z_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 197);

            auto g_z_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 198);

            auto g_z_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 199);

            auto g_z_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 200);

            auto g_z_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 201);

            auto g_z_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 202);

            auto g_z_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 203);

            auto g_z_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 204);

            auto g_z_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 205);

            auto g_z_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 206);

            auto g_z_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 207);

            auto g_z_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 208);

            auto g_z_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 209);

            auto g_z_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 210);

            auto g_z_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 211);

            auto g_z_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 212);

            auto g_z_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 213);

            auto g_z_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 214);

            auto g_z_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 215);

            auto g_z_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 216);

            auto g_z_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 217);

            auto g_z_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 218);

            auto g_z_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 219);

            auto g_z_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 220);

            auto g_z_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 221);

            auto g_z_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 222);

            auto g_z_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 223);

            auto g_z_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 224);

            auto g_z_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 225);

            auto g_z_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 226);

            auto g_z_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 227);

            auto g_z_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 228);

            auto g_z_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 229);

            auto g_z_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 230);

            auto g_z_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 231);

            auto g_z_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 232);

            auto g_z_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 233);

            auto g_z_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 234);

            auto g_z_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 235);

            auto g_z_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 236);

            auto g_z_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 237);

            auto g_z_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 238);

            auto g_z_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 239);

            auto g_z_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 240);

            auto g_z_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 241);

            auto g_z_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 242);

            auto g_z_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 243);

            auto g_z_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 244);

            auto g_z_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 245);

            auto g_z_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 246);

            auto g_z_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 247);

            auto g_z_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 248);

            auto g_z_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 249);

            auto g_z_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 250);

            auto g_z_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 251);

            auto g_z_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 252);

            auto g_z_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 253);

            auto g_z_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 254);

            auto g_z_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 255);

            auto g_z_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 256);

            auto g_z_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 257);

            auto g_z_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 258);

            auto g_z_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 259);

            auto g_z_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 260);

            auto g_z_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 261);

            auto g_z_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 262);

            auto g_z_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 263);

            auto g_z_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 264);

            auto g_z_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 265);

            auto g_z_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 266);

            auto g_z_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 267);

            auto g_z_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 268);

            auto g_z_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 269);

            auto g_z_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 270);

            auto g_z_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 271);

            auto g_z_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 272);

            auto g_z_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 273);

            auto g_z_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 274);

            auto g_z_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 275);

            auto g_z_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 276);

            auto g_z_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 277);

            auto g_z_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 278);

            auto g_z_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 279);

            auto g_z_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 280);

            auto g_z_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 281);

            auto g_z_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 282);

            auto g_z_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 283);

            auto g_z_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 284);

            auto g_z_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 285);

            auto g_z_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 286);

            auto g_z_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 287);

            auto g_z_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 288);

            auto g_z_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 289);

            auto g_z_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 290);

            auto g_z_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 291);

            auto g_z_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 292);

            auto g_z_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 293);

            auto g_z_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 294);

            auto g_z_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 295);

            auto g_z_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 296);

            auto g_z_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 297);

            auto g_z_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 298);

            auto g_z_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 299);

            auto g_z_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 300);

            auto g_z_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 301);

            auto g_z_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 302);

            auto g_z_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 303);

            auto g_z_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 304);

            auto g_z_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 305);

            auto g_z_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 306);

            auto g_z_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 307);

            auto g_z_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 308);

            auto g_z_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 309);

            auto g_z_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 310);

            auto g_z_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 311);

            auto g_z_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 312);

            auto g_z_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 313);

            auto g_z_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 314);

            auto g_z_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 315);

            auto g_z_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 316);

            auto g_z_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 317);

            auto g_z_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 318);

            auto g_z_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 319);

            auto g_z_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 320);

            auto g_z_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 321);

            auto g_z_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 322);

            auto g_z_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 323);

            auto g_z_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 324);

            auto g_z_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 325);

            auto g_z_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 326);

            auto g_z_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 327);

            auto g_z_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 328);

            auto g_z_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 329);

            auto g_z_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 330);

            auto g_z_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 331);

            auto g_z_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 332);

            auto g_z_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 333);

            auto g_z_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 334);

            auto g_z_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 335);

            auto g_z_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 336);

            auto g_z_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 337);

            auto g_z_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 338);

            auto g_z_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 339);

            auto g_z_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 340);

            auto g_z_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 341);

            auto g_z_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 342);

            auto g_z_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 343);

            auto g_z_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 344);

            auto g_z_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 345);

            auto g_z_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 346);

            auto g_z_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 347);

            auto g_z_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 348);

            auto g_z_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 349);

            auto g_z_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 350);

            auto g_z_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 351);

            auto g_z_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 352);

            auto g_z_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 353);

            auto g_z_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 354);

            auto g_z_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 355);

            auto g_z_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 356);

            auto g_z_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 357);

            auto g_z_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 358);

            auto g_z_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 359);

            auto g_z_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 360);

            auto g_z_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 361);

            auto g_z_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 362);

            auto g_z_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 363);

            auto g_z_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 364);

            auto g_z_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 365);

            auto g_z_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 366);

            auto g_z_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 367);

            auto g_z_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 368);

            auto g_z_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 369);

            auto g_z_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 370);

            auto g_z_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 371);

            auto g_z_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 372);

            auto g_z_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 373);

            auto g_z_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 374);

            auto g_z_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 375);

            auto g_z_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 376);

            auto g_z_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 377);

            auto g_z_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 378);

            auto g_z_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 379);

            auto g_z_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 380);

            auto g_z_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 381);

            auto g_z_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 382);

            auto g_z_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 383);

            auto g_z_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 384);

            auto g_z_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 385);

            auto g_z_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 386);

            auto g_z_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 387);

            auto g_z_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 388);

            auto g_z_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 389);

            auto g_z_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 390);

            auto g_z_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 391);

            auto g_z_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 392);

            auto g_z_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 393);

            auto g_z_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 394);

            auto g_z_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 395);

            auto g_z_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 396);

            auto g_z_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 397);

            auto g_z_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 398);

            auto g_z_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 399);

            auto g_z_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 400);

            auto g_z_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 401);

            auto g_z_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 402);

            auto g_z_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 403);

            auto g_z_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 404);

            auto g_z_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 405);

            auto g_z_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 406);

            auto g_z_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 407);

            auto g_z_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 408);

            auto g_z_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 409);

            auto g_z_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 410);

            auto g_z_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 411);

            auto g_z_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 412);

            auto g_z_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 413);

            auto g_z_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 414);

            auto g_z_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 415);

            auto g_z_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 416);

            auto g_z_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 417);

            auto g_z_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 418);

            auto g_z_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 419);

            /// set up bra offset for contr_buffer_xxhh

            const auto hh_geom_10_off = idx_geom_10_xxhh + (i * bcomps + j) * 441;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_zzzzz, g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_zzzzz, g_xxxx_xxxxx, g_xxxx_xxxxy, g_xxxx_xxxxz, g_xxxx_xxxyy, g_xxxx_xxxyz, g_xxxx_xxxzz, g_xxxx_xxyyy, g_xxxx_xxyyz, g_xxxx_xxyzz, g_xxxx_xxzzz, g_xxxx_xyyyy, g_xxxx_xyyyz, g_xxxx_xyyzz, g_xxxx_xyzzz, g_xxxx_xzzzz, g_xxxx_yyyyy, g_xxxx_yyyyz, g_xxxx_yyyzz, g_xxxx_yyzzz, g_xxxx_yzzzz, g_xxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xxxxx[k] = -g_xxxx_xxxxx[k] - g_x_0_xxxx_xxxxx[k] * cd_x[k] + g_x_0_xxxx_xxxxxx[k];

                g_x_0_xxxxx_xxxxy[k] = -g_xxxx_xxxxy[k] - g_x_0_xxxx_xxxxy[k] * cd_x[k] + g_x_0_xxxx_xxxxxy[k];

                g_x_0_xxxxx_xxxxz[k] = -g_xxxx_xxxxz[k] - g_x_0_xxxx_xxxxz[k] * cd_x[k] + g_x_0_xxxx_xxxxxz[k];

                g_x_0_xxxxx_xxxyy[k] = -g_xxxx_xxxyy[k] - g_x_0_xxxx_xxxyy[k] * cd_x[k] + g_x_0_xxxx_xxxxyy[k];

                g_x_0_xxxxx_xxxyz[k] = -g_xxxx_xxxyz[k] - g_x_0_xxxx_xxxyz[k] * cd_x[k] + g_x_0_xxxx_xxxxyz[k];

                g_x_0_xxxxx_xxxzz[k] = -g_xxxx_xxxzz[k] - g_x_0_xxxx_xxxzz[k] * cd_x[k] + g_x_0_xxxx_xxxxzz[k];

                g_x_0_xxxxx_xxyyy[k] = -g_xxxx_xxyyy[k] - g_x_0_xxxx_xxyyy[k] * cd_x[k] + g_x_0_xxxx_xxxyyy[k];

                g_x_0_xxxxx_xxyyz[k] = -g_xxxx_xxyyz[k] - g_x_0_xxxx_xxyyz[k] * cd_x[k] + g_x_0_xxxx_xxxyyz[k];

                g_x_0_xxxxx_xxyzz[k] = -g_xxxx_xxyzz[k] - g_x_0_xxxx_xxyzz[k] * cd_x[k] + g_x_0_xxxx_xxxyzz[k];

                g_x_0_xxxxx_xxzzz[k] = -g_xxxx_xxzzz[k] - g_x_0_xxxx_xxzzz[k] * cd_x[k] + g_x_0_xxxx_xxxzzz[k];

                g_x_0_xxxxx_xyyyy[k] = -g_xxxx_xyyyy[k] - g_x_0_xxxx_xyyyy[k] * cd_x[k] + g_x_0_xxxx_xxyyyy[k];

                g_x_0_xxxxx_xyyyz[k] = -g_xxxx_xyyyz[k] - g_x_0_xxxx_xyyyz[k] * cd_x[k] + g_x_0_xxxx_xxyyyz[k];

                g_x_0_xxxxx_xyyzz[k] = -g_xxxx_xyyzz[k] - g_x_0_xxxx_xyyzz[k] * cd_x[k] + g_x_0_xxxx_xxyyzz[k];

                g_x_0_xxxxx_xyzzz[k] = -g_xxxx_xyzzz[k] - g_x_0_xxxx_xyzzz[k] * cd_x[k] + g_x_0_xxxx_xxyzzz[k];

                g_x_0_xxxxx_xzzzz[k] = -g_xxxx_xzzzz[k] - g_x_0_xxxx_xzzzz[k] * cd_x[k] + g_x_0_xxxx_xxzzzz[k];

                g_x_0_xxxxx_yyyyy[k] = -g_xxxx_yyyyy[k] - g_x_0_xxxx_yyyyy[k] * cd_x[k] + g_x_0_xxxx_xyyyyy[k];

                g_x_0_xxxxx_yyyyz[k] = -g_xxxx_yyyyz[k] - g_x_0_xxxx_yyyyz[k] * cd_x[k] + g_x_0_xxxx_xyyyyz[k];

                g_x_0_xxxxx_yyyzz[k] = -g_xxxx_yyyzz[k] - g_x_0_xxxx_yyyzz[k] * cd_x[k] + g_x_0_xxxx_xyyyzz[k];

                g_x_0_xxxxx_yyzzz[k] = -g_xxxx_yyzzz[k] - g_x_0_xxxx_yyzzz[k] * cd_x[k] + g_x_0_xxxx_xyyzzz[k];

                g_x_0_xxxxx_yzzzz[k] = -g_xxxx_yzzzz[k] - g_x_0_xxxx_yzzzz[k] * cd_x[k] + g_x_0_xxxx_xyzzzz[k];

                g_x_0_xxxxx_zzzzz[k] = -g_xxxx_zzzzz[k] - g_x_0_xxxx_zzzzz[k] * cd_x[k] + g_x_0_xxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_zzzzz, g_x_0_xxxxy_xxxxx, g_x_0_xxxxy_xxxxy, g_x_0_xxxxy_xxxxz, g_x_0_xxxxy_xxxyy, g_x_0_xxxxy_xxxyz, g_x_0_xxxxy_xxxzz, g_x_0_xxxxy_xxyyy, g_x_0_xxxxy_xxyyz, g_x_0_xxxxy_xxyzz, g_x_0_xxxxy_xxzzz, g_x_0_xxxxy_xyyyy, g_x_0_xxxxy_xyyyz, g_x_0_xxxxy_xyyzz, g_x_0_xxxxy_xyzzz, g_x_0_xxxxy_xzzzz, g_x_0_xxxxy_yyyyy, g_x_0_xxxxy_yyyyz, g_x_0_xxxxy_yyyzz, g_x_0_xxxxy_yyzzz, g_x_0_xxxxy_yzzzz, g_x_0_xxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xxxxx[k] = -g_x_0_xxxx_xxxxx[k] * cd_y[k] + g_x_0_xxxx_xxxxxy[k];

                g_x_0_xxxxy_xxxxy[k] = -g_x_0_xxxx_xxxxy[k] * cd_y[k] + g_x_0_xxxx_xxxxyy[k];

                g_x_0_xxxxy_xxxxz[k] = -g_x_0_xxxx_xxxxz[k] * cd_y[k] + g_x_0_xxxx_xxxxyz[k];

                g_x_0_xxxxy_xxxyy[k] = -g_x_0_xxxx_xxxyy[k] * cd_y[k] + g_x_0_xxxx_xxxyyy[k];

                g_x_0_xxxxy_xxxyz[k] = -g_x_0_xxxx_xxxyz[k] * cd_y[k] + g_x_0_xxxx_xxxyyz[k];

                g_x_0_xxxxy_xxxzz[k] = -g_x_0_xxxx_xxxzz[k] * cd_y[k] + g_x_0_xxxx_xxxyzz[k];

                g_x_0_xxxxy_xxyyy[k] = -g_x_0_xxxx_xxyyy[k] * cd_y[k] + g_x_0_xxxx_xxyyyy[k];

                g_x_0_xxxxy_xxyyz[k] = -g_x_0_xxxx_xxyyz[k] * cd_y[k] + g_x_0_xxxx_xxyyyz[k];

                g_x_0_xxxxy_xxyzz[k] = -g_x_0_xxxx_xxyzz[k] * cd_y[k] + g_x_0_xxxx_xxyyzz[k];

                g_x_0_xxxxy_xxzzz[k] = -g_x_0_xxxx_xxzzz[k] * cd_y[k] + g_x_0_xxxx_xxyzzz[k];

                g_x_0_xxxxy_xyyyy[k] = -g_x_0_xxxx_xyyyy[k] * cd_y[k] + g_x_0_xxxx_xyyyyy[k];

                g_x_0_xxxxy_xyyyz[k] = -g_x_0_xxxx_xyyyz[k] * cd_y[k] + g_x_0_xxxx_xyyyyz[k];

                g_x_0_xxxxy_xyyzz[k] = -g_x_0_xxxx_xyyzz[k] * cd_y[k] + g_x_0_xxxx_xyyyzz[k];

                g_x_0_xxxxy_xyzzz[k] = -g_x_0_xxxx_xyzzz[k] * cd_y[k] + g_x_0_xxxx_xyyzzz[k];

                g_x_0_xxxxy_xzzzz[k] = -g_x_0_xxxx_xzzzz[k] * cd_y[k] + g_x_0_xxxx_xyzzzz[k];

                g_x_0_xxxxy_yyyyy[k] = -g_x_0_xxxx_yyyyy[k] * cd_y[k] + g_x_0_xxxx_yyyyyy[k];

                g_x_0_xxxxy_yyyyz[k] = -g_x_0_xxxx_yyyyz[k] * cd_y[k] + g_x_0_xxxx_yyyyyz[k];

                g_x_0_xxxxy_yyyzz[k] = -g_x_0_xxxx_yyyzz[k] * cd_y[k] + g_x_0_xxxx_yyyyzz[k];

                g_x_0_xxxxy_yyzzz[k] = -g_x_0_xxxx_yyzzz[k] * cd_y[k] + g_x_0_xxxx_yyyzzz[k];

                g_x_0_xxxxy_yzzzz[k] = -g_x_0_xxxx_yzzzz[k] * cd_y[k] + g_x_0_xxxx_yyzzzz[k];

                g_x_0_xxxxy_zzzzz[k] = -g_x_0_xxxx_zzzzz[k] * cd_y[k] + g_x_0_xxxx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_z, g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_zzzzz, g_x_0_xxxx_zzzzzz, g_x_0_xxxxz_xxxxx, g_x_0_xxxxz_xxxxy, g_x_0_xxxxz_xxxxz, g_x_0_xxxxz_xxxyy, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxzz, g_x_0_xxxxz_xxyyy, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxzzz, g_x_0_xxxxz_xyyyy, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xzzzz, g_x_0_xxxxz_yyyyy, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xxxxx[k] = -g_x_0_xxxx_xxxxx[k] * cd_z[k] + g_x_0_xxxx_xxxxxz[k];

                g_x_0_xxxxz_xxxxy[k] = -g_x_0_xxxx_xxxxy[k] * cd_z[k] + g_x_0_xxxx_xxxxyz[k];

                g_x_0_xxxxz_xxxxz[k] = -g_x_0_xxxx_xxxxz[k] * cd_z[k] + g_x_0_xxxx_xxxxzz[k];

                g_x_0_xxxxz_xxxyy[k] = -g_x_0_xxxx_xxxyy[k] * cd_z[k] + g_x_0_xxxx_xxxyyz[k];

                g_x_0_xxxxz_xxxyz[k] = -g_x_0_xxxx_xxxyz[k] * cd_z[k] + g_x_0_xxxx_xxxyzz[k];

                g_x_0_xxxxz_xxxzz[k] = -g_x_0_xxxx_xxxzz[k] * cd_z[k] + g_x_0_xxxx_xxxzzz[k];

                g_x_0_xxxxz_xxyyy[k] = -g_x_0_xxxx_xxyyy[k] * cd_z[k] + g_x_0_xxxx_xxyyyz[k];

                g_x_0_xxxxz_xxyyz[k] = -g_x_0_xxxx_xxyyz[k] * cd_z[k] + g_x_0_xxxx_xxyyzz[k];

                g_x_0_xxxxz_xxyzz[k] = -g_x_0_xxxx_xxyzz[k] * cd_z[k] + g_x_0_xxxx_xxyzzz[k];

                g_x_0_xxxxz_xxzzz[k] = -g_x_0_xxxx_xxzzz[k] * cd_z[k] + g_x_0_xxxx_xxzzzz[k];

                g_x_0_xxxxz_xyyyy[k] = -g_x_0_xxxx_xyyyy[k] * cd_z[k] + g_x_0_xxxx_xyyyyz[k];

                g_x_0_xxxxz_xyyyz[k] = -g_x_0_xxxx_xyyyz[k] * cd_z[k] + g_x_0_xxxx_xyyyzz[k];

                g_x_0_xxxxz_xyyzz[k] = -g_x_0_xxxx_xyyzz[k] * cd_z[k] + g_x_0_xxxx_xyyzzz[k];

                g_x_0_xxxxz_xyzzz[k] = -g_x_0_xxxx_xyzzz[k] * cd_z[k] + g_x_0_xxxx_xyzzzz[k];

                g_x_0_xxxxz_xzzzz[k] = -g_x_0_xxxx_xzzzz[k] * cd_z[k] + g_x_0_xxxx_xzzzzz[k];

                g_x_0_xxxxz_yyyyy[k] = -g_x_0_xxxx_yyyyy[k] * cd_z[k] + g_x_0_xxxx_yyyyyz[k];

                g_x_0_xxxxz_yyyyz[k] = -g_x_0_xxxx_yyyyz[k] * cd_z[k] + g_x_0_xxxx_yyyyzz[k];

                g_x_0_xxxxz_yyyzz[k] = -g_x_0_xxxx_yyyzz[k] * cd_z[k] + g_x_0_xxxx_yyyzzz[k];

                g_x_0_xxxxz_yyzzz[k] = -g_x_0_xxxx_yyzzz[k] * cd_z[k] + g_x_0_xxxx_yyzzzz[k];

                g_x_0_xxxxz_yzzzz[k] = -g_x_0_xxxx_yzzzz[k] * cd_z[k] + g_x_0_xxxx_yzzzzz[k];

                g_x_0_xxxxz_zzzzz[k] = -g_x_0_xxxx_zzzzz[k] * cd_z[k] + g_x_0_xxxx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_y, g_x_0_xxxy_xxxxx, g_x_0_xxxy_xxxxxy, g_x_0_xxxy_xxxxy, g_x_0_xxxy_xxxxyy, g_x_0_xxxy_xxxxyz, g_x_0_xxxy_xxxxz, g_x_0_xxxy_xxxyy, g_x_0_xxxy_xxxyyy, g_x_0_xxxy_xxxyyz, g_x_0_xxxy_xxxyz, g_x_0_xxxy_xxxyzz, g_x_0_xxxy_xxxzz, g_x_0_xxxy_xxyyy, g_x_0_xxxy_xxyyyy, g_x_0_xxxy_xxyyyz, g_x_0_xxxy_xxyyz, g_x_0_xxxy_xxyyzz, g_x_0_xxxy_xxyzz, g_x_0_xxxy_xxyzzz, g_x_0_xxxy_xxzzz, g_x_0_xxxy_xyyyy, g_x_0_xxxy_xyyyyy, g_x_0_xxxy_xyyyyz, g_x_0_xxxy_xyyyz, g_x_0_xxxy_xyyyzz, g_x_0_xxxy_xyyzz, g_x_0_xxxy_xyyzzz, g_x_0_xxxy_xyzzz, g_x_0_xxxy_xyzzzz, g_x_0_xxxy_xzzzz, g_x_0_xxxy_yyyyy, g_x_0_xxxy_yyyyyy, g_x_0_xxxy_yyyyyz, g_x_0_xxxy_yyyyz, g_x_0_xxxy_yyyyzz, g_x_0_xxxy_yyyzz, g_x_0_xxxy_yyyzzz, g_x_0_xxxy_yyzzz, g_x_0_xxxy_yyzzzz, g_x_0_xxxy_yzzzz, g_x_0_xxxy_yzzzzz, g_x_0_xxxy_zzzzz, g_x_0_xxxyy_xxxxx, g_x_0_xxxyy_xxxxy, g_x_0_xxxyy_xxxxz, g_x_0_xxxyy_xxxyy, g_x_0_xxxyy_xxxyz, g_x_0_xxxyy_xxxzz, g_x_0_xxxyy_xxyyy, g_x_0_xxxyy_xxyyz, g_x_0_xxxyy_xxyzz, g_x_0_xxxyy_xxzzz, g_x_0_xxxyy_xyyyy, g_x_0_xxxyy_xyyyz, g_x_0_xxxyy_xyyzz, g_x_0_xxxyy_xyzzz, g_x_0_xxxyy_xzzzz, g_x_0_xxxyy_yyyyy, g_x_0_xxxyy_yyyyz, g_x_0_xxxyy_yyyzz, g_x_0_xxxyy_yyzzz, g_x_0_xxxyy_yzzzz, g_x_0_xxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xxxxx[k] = -g_x_0_xxxy_xxxxx[k] * cd_y[k] + g_x_0_xxxy_xxxxxy[k];

                g_x_0_xxxyy_xxxxy[k] = -g_x_0_xxxy_xxxxy[k] * cd_y[k] + g_x_0_xxxy_xxxxyy[k];

                g_x_0_xxxyy_xxxxz[k] = -g_x_0_xxxy_xxxxz[k] * cd_y[k] + g_x_0_xxxy_xxxxyz[k];

                g_x_0_xxxyy_xxxyy[k] = -g_x_0_xxxy_xxxyy[k] * cd_y[k] + g_x_0_xxxy_xxxyyy[k];

                g_x_0_xxxyy_xxxyz[k] = -g_x_0_xxxy_xxxyz[k] * cd_y[k] + g_x_0_xxxy_xxxyyz[k];

                g_x_0_xxxyy_xxxzz[k] = -g_x_0_xxxy_xxxzz[k] * cd_y[k] + g_x_0_xxxy_xxxyzz[k];

                g_x_0_xxxyy_xxyyy[k] = -g_x_0_xxxy_xxyyy[k] * cd_y[k] + g_x_0_xxxy_xxyyyy[k];

                g_x_0_xxxyy_xxyyz[k] = -g_x_0_xxxy_xxyyz[k] * cd_y[k] + g_x_0_xxxy_xxyyyz[k];

                g_x_0_xxxyy_xxyzz[k] = -g_x_0_xxxy_xxyzz[k] * cd_y[k] + g_x_0_xxxy_xxyyzz[k];

                g_x_0_xxxyy_xxzzz[k] = -g_x_0_xxxy_xxzzz[k] * cd_y[k] + g_x_0_xxxy_xxyzzz[k];

                g_x_0_xxxyy_xyyyy[k] = -g_x_0_xxxy_xyyyy[k] * cd_y[k] + g_x_0_xxxy_xyyyyy[k];

                g_x_0_xxxyy_xyyyz[k] = -g_x_0_xxxy_xyyyz[k] * cd_y[k] + g_x_0_xxxy_xyyyyz[k];

                g_x_0_xxxyy_xyyzz[k] = -g_x_0_xxxy_xyyzz[k] * cd_y[k] + g_x_0_xxxy_xyyyzz[k];

                g_x_0_xxxyy_xyzzz[k] = -g_x_0_xxxy_xyzzz[k] * cd_y[k] + g_x_0_xxxy_xyyzzz[k];

                g_x_0_xxxyy_xzzzz[k] = -g_x_0_xxxy_xzzzz[k] * cd_y[k] + g_x_0_xxxy_xyzzzz[k];

                g_x_0_xxxyy_yyyyy[k] = -g_x_0_xxxy_yyyyy[k] * cd_y[k] + g_x_0_xxxy_yyyyyy[k];

                g_x_0_xxxyy_yyyyz[k] = -g_x_0_xxxy_yyyyz[k] * cd_y[k] + g_x_0_xxxy_yyyyyz[k];

                g_x_0_xxxyy_yyyzz[k] = -g_x_0_xxxy_yyyzz[k] * cd_y[k] + g_x_0_xxxy_yyyyzz[k];

                g_x_0_xxxyy_yyzzz[k] = -g_x_0_xxxy_yyzzz[k] * cd_y[k] + g_x_0_xxxy_yyyzzz[k];

                g_x_0_xxxyy_yzzzz[k] = -g_x_0_xxxy_yzzzz[k] * cd_y[k] + g_x_0_xxxy_yyzzzz[k];

                g_x_0_xxxyy_zzzzz[k] = -g_x_0_xxxy_zzzzz[k] * cd_y[k] + g_x_0_xxxy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyz_xxxxx, g_x_0_xxxyz_xxxxy, g_x_0_xxxyz_xxxxz, g_x_0_xxxyz_xxxyy, g_x_0_xxxyz_xxxyz, g_x_0_xxxyz_xxxzz, g_x_0_xxxyz_xxyyy, g_x_0_xxxyz_xxyyz, g_x_0_xxxyz_xxyzz, g_x_0_xxxyz_xxzzz, g_x_0_xxxyz_xyyyy, g_x_0_xxxyz_xyyyz, g_x_0_xxxyz_xyyzz, g_x_0_xxxyz_xyzzz, g_x_0_xxxyz_xzzzz, g_x_0_xxxyz_yyyyy, g_x_0_xxxyz_yyyyz, g_x_0_xxxyz_yyyzz, g_x_0_xxxyz_yyzzz, g_x_0_xxxyz_yzzzz, g_x_0_xxxyz_zzzzz, g_x_0_xxxz_xxxxx, g_x_0_xxxz_xxxxxy, g_x_0_xxxz_xxxxy, g_x_0_xxxz_xxxxyy, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxz, g_x_0_xxxz_xxxyy, g_x_0_xxxz_xxxyyy, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxzz, g_x_0_xxxz_xxyyy, g_x_0_xxxz_xxyyyy, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxzzz, g_x_0_xxxz_xyyyy, g_x_0_xxxz_xyyyyy, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xzzzz, g_x_0_xxxz_yyyyy, g_x_0_xxxz_yyyyyy, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xxxxx[k] = -g_x_0_xxxz_xxxxx[k] * cd_y[k] + g_x_0_xxxz_xxxxxy[k];

                g_x_0_xxxyz_xxxxy[k] = -g_x_0_xxxz_xxxxy[k] * cd_y[k] + g_x_0_xxxz_xxxxyy[k];

                g_x_0_xxxyz_xxxxz[k] = -g_x_0_xxxz_xxxxz[k] * cd_y[k] + g_x_0_xxxz_xxxxyz[k];

                g_x_0_xxxyz_xxxyy[k] = -g_x_0_xxxz_xxxyy[k] * cd_y[k] + g_x_0_xxxz_xxxyyy[k];

                g_x_0_xxxyz_xxxyz[k] = -g_x_0_xxxz_xxxyz[k] * cd_y[k] + g_x_0_xxxz_xxxyyz[k];

                g_x_0_xxxyz_xxxzz[k] = -g_x_0_xxxz_xxxzz[k] * cd_y[k] + g_x_0_xxxz_xxxyzz[k];

                g_x_0_xxxyz_xxyyy[k] = -g_x_0_xxxz_xxyyy[k] * cd_y[k] + g_x_0_xxxz_xxyyyy[k];

                g_x_0_xxxyz_xxyyz[k] = -g_x_0_xxxz_xxyyz[k] * cd_y[k] + g_x_0_xxxz_xxyyyz[k];

                g_x_0_xxxyz_xxyzz[k] = -g_x_0_xxxz_xxyzz[k] * cd_y[k] + g_x_0_xxxz_xxyyzz[k];

                g_x_0_xxxyz_xxzzz[k] = -g_x_0_xxxz_xxzzz[k] * cd_y[k] + g_x_0_xxxz_xxyzzz[k];

                g_x_0_xxxyz_xyyyy[k] = -g_x_0_xxxz_xyyyy[k] * cd_y[k] + g_x_0_xxxz_xyyyyy[k];

                g_x_0_xxxyz_xyyyz[k] = -g_x_0_xxxz_xyyyz[k] * cd_y[k] + g_x_0_xxxz_xyyyyz[k];

                g_x_0_xxxyz_xyyzz[k] = -g_x_0_xxxz_xyyzz[k] * cd_y[k] + g_x_0_xxxz_xyyyzz[k];

                g_x_0_xxxyz_xyzzz[k] = -g_x_0_xxxz_xyzzz[k] * cd_y[k] + g_x_0_xxxz_xyyzzz[k];

                g_x_0_xxxyz_xzzzz[k] = -g_x_0_xxxz_xzzzz[k] * cd_y[k] + g_x_0_xxxz_xyzzzz[k];

                g_x_0_xxxyz_yyyyy[k] = -g_x_0_xxxz_yyyyy[k] * cd_y[k] + g_x_0_xxxz_yyyyyy[k];

                g_x_0_xxxyz_yyyyz[k] = -g_x_0_xxxz_yyyyz[k] * cd_y[k] + g_x_0_xxxz_yyyyyz[k];

                g_x_0_xxxyz_yyyzz[k] = -g_x_0_xxxz_yyyzz[k] * cd_y[k] + g_x_0_xxxz_yyyyzz[k];

                g_x_0_xxxyz_yyzzz[k] = -g_x_0_xxxz_yyzzz[k] * cd_y[k] + g_x_0_xxxz_yyyzzz[k];

                g_x_0_xxxyz_yzzzz[k] = -g_x_0_xxxz_yzzzz[k] * cd_y[k] + g_x_0_xxxz_yyzzzz[k];

                g_x_0_xxxyz_zzzzz[k] = -g_x_0_xxxz_zzzzz[k] * cd_y[k] + g_x_0_xxxz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_z, g_x_0_xxxz_xxxxx, g_x_0_xxxz_xxxxxz, g_x_0_xxxz_xxxxy, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxz, g_x_0_xxxz_xxxxzz, g_x_0_xxxz_xxxyy, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxzz, g_x_0_xxxz_xxxzzz, g_x_0_xxxz_xxyyy, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxzzz, g_x_0_xxxz_xxzzzz, g_x_0_xxxz_xyyyy, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xzzzz, g_x_0_xxxz_xzzzzz, g_x_0_xxxz_yyyyy, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_zzzzz, g_x_0_xxxz_zzzzzz, g_x_0_xxxzz_xxxxx, g_x_0_xxxzz_xxxxy, g_x_0_xxxzz_xxxxz, g_x_0_xxxzz_xxxyy, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxzz, g_x_0_xxxzz_xxyyy, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxzzz, g_x_0_xxxzz_xyyyy, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xzzzz, g_x_0_xxxzz_yyyyy, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xxxxx[k] = -g_x_0_xxxz_xxxxx[k] * cd_z[k] + g_x_0_xxxz_xxxxxz[k];

                g_x_0_xxxzz_xxxxy[k] = -g_x_0_xxxz_xxxxy[k] * cd_z[k] + g_x_0_xxxz_xxxxyz[k];

                g_x_0_xxxzz_xxxxz[k] = -g_x_0_xxxz_xxxxz[k] * cd_z[k] + g_x_0_xxxz_xxxxzz[k];

                g_x_0_xxxzz_xxxyy[k] = -g_x_0_xxxz_xxxyy[k] * cd_z[k] + g_x_0_xxxz_xxxyyz[k];

                g_x_0_xxxzz_xxxyz[k] = -g_x_0_xxxz_xxxyz[k] * cd_z[k] + g_x_0_xxxz_xxxyzz[k];

                g_x_0_xxxzz_xxxzz[k] = -g_x_0_xxxz_xxxzz[k] * cd_z[k] + g_x_0_xxxz_xxxzzz[k];

                g_x_0_xxxzz_xxyyy[k] = -g_x_0_xxxz_xxyyy[k] * cd_z[k] + g_x_0_xxxz_xxyyyz[k];

                g_x_0_xxxzz_xxyyz[k] = -g_x_0_xxxz_xxyyz[k] * cd_z[k] + g_x_0_xxxz_xxyyzz[k];

                g_x_0_xxxzz_xxyzz[k] = -g_x_0_xxxz_xxyzz[k] * cd_z[k] + g_x_0_xxxz_xxyzzz[k];

                g_x_0_xxxzz_xxzzz[k] = -g_x_0_xxxz_xxzzz[k] * cd_z[k] + g_x_0_xxxz_xxzzzz[k];

                g_x_0_xxxzz_xyyyy[k] = -g_x_0_xxxz_xyyyy[k] * cd_z[k] + g_x_0_xxxz_xyyyyz[k];

                g_x_0_xxxzz_xyyyz[k] = -g_x_0_xxxz_xyyyz[k] * cd_z[k] + g_x_0_xxxz_xyyyzz[k];

                g_x_0_xxxzz_xyyzz[k] = -g_x_0_xxxz_xyyzz[k] * cd_z[k] + g_x_0_xxxz_xyyzzz[k];

                g_x_0_xxxzz_xyzzz[k] = -g_x_0_xxxz_xyzzz[k] * cd_z[k] + g_x_0_xxxz_xyzzzz[k];

                g_x_0_xxxzz_xzzzz[k] = -g_x_0_xxxz_xzzzz[k] * cd_z[k] + g_x_0_xxxz_xzzzzz[k];

                g_x_0_xxxzz_yyyyy[k] = -g_x_0_xxxz_yyyyy[k] * cd_z[k] + g_x_0_xxxz_yyyyyz[k];

                g_x_0_xxxzz_yyyyz[k] = -g_x_0_xxxz_yyyyz[k] * cd_z[k] + g_x_0_xxxz_yyyyzz[k];

                g_x_0_xxxzz_yyyzz[k] = -g_x_0_xxxz_yyyzz[k] * cd_z[k] + g_x_0_xxxz_yyyzzz[k];

                g_x_0_xxxzz_yyzzz[k] = -g_x_0_xxxz_yyzzz[k] * cd_z[k] + g_x_0_xxxz_yyzzzz[k];

                g_x_0_xxxzz_yzzzz[k] = -g_x_0_xxxz_yzzzz[k] * cd_z[k] + g_x_0_xxxz_yzzzzz[k];

                g_x_0_xxxzz_zzzzz[k] = -g_x_0_xxxz_zzzzz[k] * cd_z[k] + g_x_0_xxxz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_y, g_x_0_xxyy_xxxxx, g_x_0_xxyy_xxxxxy, g_x_0_xxyy_xxxxy, g_x_0_xxyy_xxxxyy, g_x_0_xxyy_xxxxyz, g_x_0_xxyy_xxxxz, g_x_0_xxyy_xxxyy, g_x_0_xxyy_xxxyyy, g_x_0_xxyy_xxxyyz, g_x_0_xxyy_xxxyz, g_x_0_xxyy_xxxyzz, g_x_0_xxyy_xxxzz, g_x_0_xxyy_xxyyy, g_x_0_xxyy_xxyyyy, g_x_0_xxyy_xxyyyz, g_x_0_xxyy_xxyyz, g_x_0_xxyy_xxyyzz, g_x_0_xxyy_xxyzz, g_x_0_xxyy_xxyzzz, g_x_0_xxyy_xxzzz, g_x_0_xxyy_xyyyy, g_x_0_xxyy_xyyyyy, g_x_0_xxyy_xyyyyz, g_x_0_xxyy_xyyyz, g_x_0_xxyy_xyyyzz, g_x_0_xxyy_xyyzz, g_x_0_xxyy_xyyzzz, g_x_0_xxyy_xyzzz, g_x_0_xxyy_xyzzzz, g_x_0_xxyy_xzzzz, g_x_0_xxyy_yyyyy, g_x_0_xxyy_yyyyyy, g_x_0_xxyy_yyyyyz, g_x_0_xxyy_yyyyz, g_x_0_xxyy_yyyyzz, g_x_0_xxyy_yyyzz, g_x_0_xxyy_yyyzzz, g_x_0_xxyy_yyzzz, g_x_0_xxyy_yyzzzz, g_x_0_xxyy_yzzzz, g_x_0_xxyy_yzzzzz, g_x_0_xxyy_zzzzz, g_x_0_xxyyy_xxxxx, g_x_0_xxyyy_xxxxy, g_x_0_xxyyy_xxxxz, g_x_0_xxyyy_xxxyy, g_x_0_xxyyy_xxxyz, g_x_0_xxyyy_xxxzz, g_x_0_xxyyy_xxyyy, g_x_0_xxyyy_xxyyz, g_x_0_xxyyy_xxyzz, g_x_0_xxyyy_xxzzz, g_x_0_xxyyy_xyyyy, g_x_0_xxyyy_xyyyz, g_x_0_xxyyy_xyyzz, g_x_0_xxyyy_xyzzz, g_x_0_xxyyy_xzzzz, g_x_0_xxyyy_yyyyy, g_x_0_xxyyy_yyyyz, g_x_0_xxyyy_yyyzz, g_x_0_xxyyy_yyzzz, g_x_0_xxyyy_yzzzz, g_x_0_xxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xxxxx[k] = -g_x_0_xxyy_xxxxx[k] * cd_y[k] + g_x_0_xxyy_xxxxxy[k];

                g_x_0_xxyyy_xxxxy[k] = -g_x_0_xxyy_xxxxy[k] * cd_y[k] + g_x_0_xxyy_xxxxyy[k];

                g_x_0_xxyyy_xxxxz[k] = -g_x_0_xxyy_xxxxz[k] * cd_y[k] + g_x_0_xxyy_xxxxyz[k];

                g_x_0_xxyyy_xxxyy[k] = -g_x_0_xxyy_xxxyy[k] * cd_y[k] + g_x_0_xxyy_xxxyyy[k];

                g_x_0_xxyyy_xxxyz[k] = -g_x_0_xxyy_xxxyz[k] * cd_y[k] + g_x_0_xxyy_xxxyyz[k];

                g_x_0_xxyyy_xxxzz[k] = -g_x_0_xxyy_xxxzz[k] * cd_y[k] + g_x_0_xxyy_xxxyzz[k];

                g_x_0_xxyyy_xxyyy[k] = -g_x_0_xxyy_xxyyy[k] * cd_y[k] + g_x_0_xxyy_xxyyyy[k];

                g_x_0_xxyyy_xxyyz[k] = -g_x_0_xxyy_xxyyz[k] * cd_y[k] + g_x_0_xxyy_xxyyyz[k];

                g_x_0_xxyyy_xxyzz[k] = -g_x_0_xxyy_xxyzz[k] * cd_y[k] + g_x_0_xxyy_xxyyzz[k];

                g_x_0_xxyyy_xxzzz[k] = -g_x_0_xxyy_xxzzz[k] * cd_y[k] + g_x_0_xxyy_xxyzzz[k];

                g_x_0_xxyyy_xyyyy[k] = -g_x_0_xxyy_xyyyy[k] * cd_y[k] + g_x_0_xxyy_xyyyyy[k];

                g_x_0_xxyyy_xyyyz[k] = -g_x_0_xxyy_xyyyz[k] * cd_y[k] + g_x_0_xxyy_xyyyyz[k];

                g_x_0_xxyyy_xyyzz[k] = -g_x_0_xxyy_xyyzz[k] * cd_y[k] + g_x_0_xxyy_xyyyzz[k];

                g_x_0_xxyyy_xyzzz[k] = -g_x_0_xxyy_xyzzz[k] * cd_y[k] + g_x_0_xxyy_xyyzzz[k];

                g_x_0_xxyyy_xzzzz[k] = -g_x_0_xxyy_xzzzz[k] * cd_y[k] + g_x_0_xxyy_xyzzzz[k];

                g_x_0_xxyyy_yyyyy[k] = -g_x_0_xxyy_yyyyy[k] * cd_y[k] + g_x_0_xxyy_yyyyyy[k];

                g_x_0_xxyyy_yyyyz[k] = -g_x_0_xxyy_yyyyz[k] * cd_y[k] + g_x_0_xxyy_yyyyyz[k];

                g_x_0_xxyyy_yyyzz[k] = -g_x_0_xxyy_yyyzz[k] * cd_y[k] + g_x_0_xxyy_yyyyzz[k];

                g_x_0_xxyyy_yyzzz[k] = -g_x_0_xxyy_yyzzz[k] * cd_y[k] + g_x_0_xxyy_yyyzzz[k];

                g_x_0_xxyyy_yzzzz[k] = -g_x_0_xxyy_yzzzz[k] * cd_y[k] + g_x_0_xxyy_yyzzzz[k];

                g_x_0_xxyyy_zzzzz[k] = -g_x_0_xxyy_zzzzz[k] * cd_y[k] + g_x_0_xxyy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyz_xxxxx, g_x_0_xxyyz_xxxxy, g_x_0_xxyyz_xxxxz, g_x_0_xxyyz_xxxyy, g_x_0_xxyyz_xxxyz, g_x_0_xxyyz_xxxzz, g_x_0_xxyyz_xxyyy, g_x_0_xxyyz_xxyyz, g_x_0_xxyyz_xxyzz, g_x_0_xxyyz_xxzzz, g_x_0_xxyyz_xyyyy, g_x_0_xxyyz_xyyyz, g_x_0_xxyyz_xyyzz, g_x_0_xxyyz_xyzzz, g_x_0_xxyyz_xzzzz, g_x_0_xxyyz_yyyyy, g_x_0_xxyyz_yyyyz, g_x_0_xxyyz_yyyzz, g_x_0_xxyyz_yyzzz, g_x_0_xxyyz_yzzzz, g_x_0_xxyyz_zzzzz, g_x_0_xxyz_xxxxx, g_x_0_xxyz_xxxxxy, g_x_0_xxyz_xxxxy, g_x_0_xxyz_xxxxyy, g_x_0_xxyz_xxxxyz, g_x_0_xxyz_xxxxz, g_x_0_xxyz_xxxyy, g_x_0_xxyz_xxxyyy, g_x_0_xxyz_xxxyyz, g_x_0_xxyz_xxxyz, g_x_0_xxyz_xxxyzz, g_x_0_xxyz_xxxzz, g_x_0_xxyz_xxyyy, g_x_0_xxyz_xxyyyy, g_x_0_xxyz_xxyyyz, g_x_0_xxyz_xxyyz, g_x_0_xxyz_xxyyzz, g_x_0_xxyz_xxyzz, g_x_0_xxyz_xxyzzz, g_x_0_xxyz_xxzzz, g_x_0_xxyz_xyyyy, g_x_0_xxyz_xyyyyy, g_x_0_xxyz_xyyyyz, g_x_0_xxyz_xyyyz, g_x_0_xxyz_xyyyzz, g_x_0_xxyz_xyyzz, g_x_0_xxyz_xyyzzz, g_x_0_xxyz_xyzzz, g_x_0_xxyz_xyzzzz, g_x_0_xxyz_xzzzz, g_x_0_xxyz_yyyyy, g_x_0_xxyz_yyyyyy, g_x_0_xxyz_yyyyyz, g_x_0_xxyz_yyyyz, g_x_0_xxyz_yyyyzz, g_x_0_xxyz_yyyzz, g_x_0_xxyz_yyyzzz, g_x_0_xxyz_yyzzz, g_x_0_xxyz_yyzzzz, g_x_0_xxyz_yzzzz, g_x_0_xxyz_yzzzzz, g_x_0_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xxxxx[k] = -g_x_0_xxyz_xxxxx[k] * cd_y[k] + g_x_0_xxyz_xxxxxy[k];

                g_x_0_xxyyz_xxxxy[k] = -g_x_0_xxyz_xxxxy[k] * cd_y[k] + g_x_0_xxyz_xxxxyy[k];

                g_x_0_xxyyz_xxxxz[k] = -g_x_0_xxyz_xxxxz[k] * cd_y[k] + g_x_0_xxyz_xxxxyz[k];

                g_x_0_xxyyz_xxxyy[k] = -g_x_0_xxyz_xxxyy[k] * cd_y[k] + g_x_0_xxyz_xxxyyy[k];

                g_x_0_xxyyz_xxxyz[k] = -g_x_0_xxyz_xxxyz[k] * cd_y[k] + g_x_0_xxyz_xxxyyz[k];

                g_x_0_xxyyz_xxxzz[k] = -g_x_0_xxyz_xxxzz[k] * cd_y[k] + g_x_0_xxyz_xxxyzz[k];

                g_x_0_xxyyz_xxyyy[k] = -g_x_0_xxyz_xxyyy[k] * cd_y[k] + g_x_0_xxyz_xxyyyy[k];

                g_x_0_xxyyz_xxyyz[k] = -g_x_0_xxyz_xxyyz[k] * cd_y[k] + g_x_0_xxyz_xxyyyz[k];

                g_x_0_xxyyz_xxyzz[k] = -g_x_0_xxyz_xxyzz[k] * cd_y[k] + g_x_0_xxyz_xxyyzz[k];

                g_x_0_xxyyz_xxzzz[k] = -g_x_0_xxyz_xxzzz[k] * cd_y[k] + g_x_0_xxyz_xxyzzz[k];

                g_x_0_xxyyz_xyyyy[k] = -g_x_0_xxyz_xyyyy[k] * cd_y[k] + g_x_0_xxyz_xyyyyy[k];

                g_x_0_xxyyz_xyyyz[k] = -g_x_0_xxyz_xyyyz[k] * cd_y[k] + g_x_0_xxyz_xyyyyz[k];

                g_x_0_xxyyz_xyyzz[k] = -g_x_0_xxyz_xyyzz[k] * cd_y[k] + g_x_0_xxyz_xyyyzz[k];

                g_x_0_xxyyz_xyzzz[k] = -g_x_0_xxyz_xyzzz[k] * cd_y[k] + g_x_0_xxyz_xyyzzz[k];

                g_x_0_xxyyz_xzzzz[k] = -g_x_0_xxyz_xzzzz[k] * cd_y[k] + g_x_0_xxyz_xyzzzz[k];

                g_x_0_xxyyz_yyyyy[k] = -g_x_0_xxyz_yyyyy[k] * cd_y[k] + g_x_0_xxyz_yyyyyy[k];

                g_x_0_xxyyz_yyyyz[k] = -g_x_0_xxyz_yyyyz[k] * cd_y[k] + g_x_0_xxyz_yyyyyz[k];

                g_x_0_xxyyz_yyyzz[k] = -g_x_0_xxyz_yyyzz[k] * cd_y[k] + g_x_0_xxyz_yyyyzz[k];

                g_x_0_xxyyz_yyzzz[k] = -g_x_0_xxyz_yyzzz[k] * cd_y[k] + g_x_0_xxyz_yyyzzz[k];

                g_x_0_xxyyz_yzzzz[k] = -g_x_0_xxyz_yzzzz[k] * cd_y[k] + g_x_0_xxyz_yyzzzz[k];

                g_x_0_xxyyz_zzzzz[k] = -g_x_0_xxyz_zzzzz[k] * cd_y[k] + g_x_0_xxyz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzz_xxxxx, g_x_0_xxyzz_xxxxy, g_x_0_xxyzz_xxxxz, g_x_0_xxyzz_xxxyy, g_x_0_xxyzz_xxxyz, g_x_0_xxyzz_xxxzz, g_x_0_xxyzz_xxyyy, g_x_0_xxyzz_xxyyz, g_x_0_xxyzz_xxyzz, g_x_0_xxyzz_xxzzz, g_x_0_xxyzz_xyyyy, g_x_0_xxyzz_xyyyz, g_x_0_xxyzz_xyyzz, g_x_0_xxyzz_xyzzz, g_x_0_xxyzz_xzzzz, g_x_0_xxyzz_yyyyy, g_x_0_xxyzz_yyyyz, g_x_0_xxyzz_yyyzz, g_x_0_xxyzz_yyzzz, g_x_0_xxyzz_yzzzz, g_x_0_xxyzz_zzzzz, g_x_0_xxzz_xxxxx, g_x_0_xxzz_xxxxxy, g_x_0_xxzz_xxxxy, g_x_0_xxzz_xxxxyy, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxz, g_x_0_xxzz_xxxyy, g_x_0_xxzz_xxxyyy, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxzz, g_x_0_xxzz_xxyyy, g_x_0_xxzz_xxyyyy, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxzzz, g_x_0_xxzz_xyyyy, g_x_0_xxzz_xyyyyy, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xzzzz, g_x_0_xxzz_yyyyy, g_x_0_xxzz_yyyyyy, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xxxxx[k] = -g_x_0_xxzz_xxxxx[k] * cd_y[k] + g_x_0_xxzz_xxxxxy[k];

                g_x_0_xxyzz_xxxxy[k] = -g_x_0_xxzz_xxxxy[k] * cd_y[k] + g_x_0_xxzz_xxxxyy[k];

                g_x_0_xxyzz_xxxxz[k] = -g_x_0_xxzz_xxxxz[k] * cd_y[k] + g_x_0_xxzz_xxxxyz[k];

                g_x_0_xxyzz_xxxyy[k] = -g_x_0_xxzz_xxxyy[k] * cd_y[k] + g_x_0_xxzz_xxxyyy[k];

                g_x_0_xxyzz_xxxyz[k] = -g_x_0_xxzz_xxxyz[k] * cd_y[k] + g_x_0_xxzz_xxxyyz[k];

                g_x_0_xxyzz_xxxzz[k] = -g_x_0_xxzz_xxxzz[k] * cd_y[k] + g_x_0_xxzz_xxxyzz[k];

                g_x_0_xxyzz_xxyyy[k] = -g_x_0_xxzz_xxyyy[k] * cd_y[k] + g_x_0_xxzz_xxyyyy[k];

                g_x_0_xxyzz_xxyyz[k] = -g_x_0_xxzz_xxyyz[k] * cd_y[k] + g_x_0_xxzz_xxyyyz[k];

                g_x_0_xxyzz_xxyzz[k] = -g_x_0_xxzz_xxyzz[k] * cd_y[k] + g_x_0_xxzz_xxyyzz[k];

                g_x_0_xxyzz_xxzzz[k] = -g_x_0_xxzz_xxzzz[k] * cd_y[k] + g_x_0_xxzz_xxyzzz[k];

                g_x_0_xxyzz_xyyyy[k] = -g_x_0_xxzz_xyyyy[k] * cd_y[k] + g_x_0_xxzz_xyyyyy[k];

                g_x_0_xxyzz_xyyyz[k] = -g_x_0_xxzz_xyyyz[k] * cd_y[k] + g_x_0_xxzz_xyyyyz[k];

                g_x_0_xxyzz_xyyzz[k] = -g_x_0_xxzz_xyyzz[k] * cd_y[k] + g_x_0_xxzz_xyyyzz[k];

                g_x_0_xxyzz_xyzzz[k] = -g_x_0_xxzz_xyzzz[k] * cd_y[k] + g_x_0_xxzz_xyyzzz[k];

                g_x_0_xxyzz_xzzzz[k] = -g_x_0_xxzz_xzzzz[k] * cd_y[k] + g_x_0_xxzz_xyzzzz[k];

                g_x_0_xxyzz_yyyyy[k] = -g_x_0_xxzz_yyyyy[k] * cd_y[k] + g_x_0_xxzz_yyyyyy[k];

                g_x_0_xxyzz_yyyyz[k] = -g_x_0_xxzz_yyyyz[k] * cd_y[k] + g_x_0_xxzz_yyyyyz[k];

                g_x_0_xxyzz_yyyzz[k] = -g_x_0_xxzz_yyyzz[k] * cd_y[k] + g_x_0_xxzz_yyyyzz[k];

                g_x_0_xxyzz_yyzzz[k] = -g_x_0_xxzz_yyzzz[k] * cd_y[k] + g_x_0_xxzz_yyyzzz[k];

                g_x_0_xxyzz_yzzzz[k] = -g_x_0_xxzz_yzzzz[k] * cd_y[k] + g_x_0_xxzz_yyzzzz[k];

                g_x_0_xxyzz_zzzzz[k] = -g_x_0_xxzz_zzzzz[k] * cd_y[k] + g_x_0_xxzz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_z, g_x_0_xxzz_xxxxx, g_x_0_xxzz_xxxxxz, g_x_0_xxzz_xxxxy, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxz, g_x_0_xxzz_xxxxzz, g_x_0_xxzz_xxxyy, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxzz, g_x_0_xxzz_xxxzzz, g_x_0_xxzz_xxyyy, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxzzz, g_x_0_xxzz_xxzzzz, g_x_0_xxzz_xyyyy, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xzzzz, g_x_0_xxzz_xzzzzz, g_x_0_xxzz_yyyyy, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_zzzzz, g_x_0_xxzz_zzzzzz, g_x_0_xxzzz_xxxxx, g_x_0_xxzzz_xxxxy, g_x_0_xxzzz_xxxxz, g_x_0_xxzzz_xxxyy, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxzz, g_x_0_xxzzz_xxyyy, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxzzz, g_x_0_xxzzz_xyyyy, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xzzzz, g_x_0_xxzzz_yyyyy, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xxxxx[k] = -g_x_0_xxzz_xxxxx[k] * cd_z[k] + g_x_0_xxzz_xxxxxz[k];

                g_x_0_xxzzz_xxxxy[k] = -g_x_0_xxzz_xxxxy[k] * cd_z[k] + g_x_0_xxzz_xxxxyz[k];

                g_x_0_xxzzz_xxxxz[k] = -g_x_0_xxzz_xxxxz[k] * cd_z[k] + g_x_0_xxzz_xxxxzz[k];

                g_x_0_xxzzz_xxxyy[k] = -g_x_0_xxzz_xxxyy[k] * cd_z[k] + g_x_0_xxzz_xxxyyz[k];

                g_x_0_xxzzz_xxxyz[k] = -g_x_0_xxzz_xxxyz[k] * cd_z[k] + g_x_0_xxzz_xxxyzz[k];

                g_x_0_xxzzz_xxxzz[k] = -g_x_0_xxzz_xxxzz[k] * cd_z[k] + g_x_0_xxzz_xxxzzz[k];

                g_x_0_xxzzz_xxyyy[k] = -g_x_0_xxzz_xxyyy[k] * cd_z[k] + g_x_0_xxzz_xxyyyz[k];

                g_x_0_xxzzz_xxyyz[k] = -g_x_0_xxzz_xxyyz[k] * cd_z[k] + g_x_0_xxzz_xxyyzz[k];

                g_x_0_xxzzz_xxyzz[k] = -g_x_0_xxzz_xxyzz[k] * cd_z[k] + g_x_0_xxzz_xxyzzz[k];

                g_x_0_xxzzz_xxzzz[k] = -g_x_0_xxzz_xxzzz[k] * cd_z[k] + g_x_0_xxzz_xxzzzz[k];

                g_x_0_xxzzz_xyyyy[k] = -g_x_0_xxzz_xyyyy[k] * cd_z[k] + g_x_0_xxzz_xyyyyz[k];

                g_x_0_xxzzz_xyyyz[k] = -g_x_0_xxzz_xyyyz[k] * cd_z[k] + g_x_0_xxzz_xyyyzz[k];

                g_x_0_xxzzz_xyyzz[k] = -g_x_0_xxzz_xyyzz[k] * cd_z[k] + g_x_0_xxzz_xyyzzz[k];

                g_x_0_xxzzz_xyzzz[k] = -g_x_0_xxzz_xyzzz[k] * cd_z[k] + g_x_0_xxzz_xyzzzz[k];

                g_x_0_xxzzz_xzzzz[k] = -g_x_0_xxzz_xzzzz[k] * cd_z[k] + g_x_0_xxzz_xzzzzz[k];

                g_x_0_xxzzz_yyyyy[k] = -g_x_0_xxzz_yyyyy[k] * cd_z[k] + g_x_0_xxzz_yyyyyz[k];

                g_x_0_xxzzz_yyyyz[k] = -g_x_0_xxzz_yyyyz[k] * cd_z[k] + g_x_0_xxzz_yyyyzz[k];

                g_x_0_xxzzz_yyyzz[k] = -g_x_0_xxzz_yyyzz[k] * cd_z[k] + g_x_0_xxzz_yyyzzz[k];

                g_x_0_xxzzz_yyzzz[k] = -g_x_0_xxzz_yyzzz[k] * cd_z[k] + g_x_0_xxzz_yyzzzz[k];

                g_x_0_xxzzz_yzzzz[k] = -g_x_0_xxzz_yzzzz[k] * cd_z[k] + g_x_0_xxzz_yzzzzz[k];

                g_x_0_xxzzz_zzzzz[k] = -g_x_0_xxzz_zzzzz[k] * cd_z[k] + g_x_0_xxzz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_y, g_x_0_xyyy_xxxxx, g_x_0_xyyy_xxxxxy, g_x_0_xyyy_xxxxy, g_x_0_xyyy_xxxxyy, g_x_0_xyyy_xxxxyz, g_x_0_xyyy_xxxxz, g_x_0_xyyy_xxxyy, g_x_0_xyyy_xxxyyy, g_x_0_xyyy_xxxyyz, g_x_0_xyyy_xxxyz, g_x_0_xyyy_xxxyzz, g_x_0_xyyy_xxxzz, g_x_0_xyyy_xxyyy, g_x_0_xyyy_xxyyyy, g_x_0_xyyy_xxyyyz, g_x_0_xyyy_xxyyz, g_x_0_xyyy_xxyyzz, g_x_0_xyyy_xxyzz, g_x_0_xyyy_xxyzzz, g_x_0_xyyy_xxzzz, g_x_0_xyyy_xyyyy, g_x_0_xyyy_xyyyyy, g_x_0_xyyy_xyyyyz, g_x_0_xyyy_xyyyz, g_x_0_xyyy_xyyyzz, g_x_0_xyyy_xyyzz, g_x_0_xyyy_xyyzzz, g_x_0_xyyy_xyzzz, g_x_0_xyyy_xyzzzz, g_x_0_xyyy_xzzzz, g_x_0_xyyy_yyyyy, g_x_0_xyyy_yyyyyy, g_x_0_xyyy_yyyyyz, g_x_0_xyyy_yyyyz, g_x_0_xyyy_yyyyzz, g_x_0_xyyy_yyyzz, g_x_0_xyyy_yyyzzz, g_x_0_xyyy_yyzzz, g_x_0_xyyy_yyzzzz, g_x_0_xyyy_yzzzz, g_x_0_xyyy_yzzzzz, g_x_0_xyyy_zzzzz, g_x_0_xyyyy_xxxxx, g_x_0_xyyyy_xxxxy, g_x_0_xyyyy_xxxxz, g_x_0_xyyyy_xxxyy, g_x_0_xyyyy_xxxyz, g_x_0_xyyyy_xxxzz, g_x_0_xyyyy_xxyyy, g_x_0_xyyyy_xxyyz, g_x_0_xyyyy_xxyzz, g_x_0_xyyyy_xxzzz, g_x_0_xyyyy_xyyyy, g_x_0_xyyyy_xyyyz, g_x_0_xyyyy_xyyzz, g_x_0_xyyyy_xyzzz, g_x_0_xyyyy_xzzzz, g_x_0_xyyyy_yyyyy, g_x_0_xyyyy_yyyyz, g_x_0_xyyyy_yyyzz, g_x_0_xyyyy_yyzzz, g_x_0_xyyyy_yzzzz, g_x_0_xyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xxxxx[k] = -g_x_0_xyyy_xxxxx[k] * cd_y[k] + g_x_0_xyyy_xxxxxy[k];

                g_x_0_xyyyy_xxxxy[k] = -g_x_0_xyyy_xxxxy[k] * cd_y[k] + g_x_0_xyyy_xxxxyy[k];

                g_x_0_xyyyy_xxxxz[k] = -g_x_0_xyyy_xxxxz[k] * cd_y[k] + g_x_0_xyyy_xxxxyz[k];

                g_x_0_xyyyy_xxxyy[k] = -g_x_0_xyyy_xxxyy[k] * cd_y[k] + g_x_0_xyyy_xxxyyy[k];

                g_x_0_xyyyy_xxxyz[k] = -g_x_0_xyyy_xxxyz[k] * cd_y[k] + g_x_0_xyyy_xxxyyz[k];

                g_x_0_xyyyy_xxxzz[k] = -g_x_0_xyyy_xxxzz[k] * cd_y[k] + g_x_0_xyyy_xxxyzz[k];

                g_x_0_xyyyy_xxyyy[k] = -g_x_0_xyyy_xxyyy[k] * cd_y[k] + g_x_0_xyyy_xxyyyy[k];

                g_x_0_xyyyy_xxyyz[k] = -g_x_0_xyyy_xxyyz[k] * cd_y[k] + g_x_0_xyyy_xxyyyz[k];

                g_x_0_xyyyy_xxyzz[k] = -g_x_0_xyyy_xxyzz[k] * cd_y[k] + g_x_0_xyyy_xxyyzz[k];

                g_x_0_xyyyy_xxzzz[k] = -g_x_0_xyyy_xxzzz[k] * cd_y[k] + g_x_0_xyyy_xxyzzz[k];

                g_x_0_xyyyy_xyyyy[k] = -g_x_0_xyyy_xyyyy[k] * cd_y[k] + g_x_0_xyyy_xyyyyy[k];

                g_x_0_xyyyy_xyyyz[k] = -g_x_0_xyyy_xyyyz[k] * cd_y[k] + g_x_0_xyyy_xyyyyz[k];

                g_x_0_xyyyy_xyyzz[k] = -g_x_0_xyyy_xyyzz[k] * cd_y[k] + g_x_0_xyyy_xyyyzz[k];

                g_x_0_xyyyy_xyzzz[k] = -g_x_0_xyyy_xyzzz[k] * cd_y[k] + g_x_0_xyyy_xyyzzz[k];

                g_x_0_xyyyy_xzzzz[k] = -g_x_0_xyyy_xzzzz[k] * cd_y[k] + g_x_0_xyyy_xyzzzz[k];

                g_x_0_xyyyy_yyyyy[k] = -g_x_0_xyyy_yyyyy[k] * cd_y[k] + g_x_0_xyyy_yyyyyy[k];

                g_x_0_xyyyy_yyyyz[k] = -g_x_0_xyyy_yyyyz[k] * cd_y[k] + g_x_0_xyyy_yyyyyz[k];

                g_x_0_xyyyy_yyyzz[k] = -g_x_0_xyyy_yyyzz[k] * cd_y[k] + g_x_0_xyyy_yyyyzz[k];

                g_x_0_xyyyy_yyzzz[k] = -g_x_0_xyyy_yyzzz[k] * cd_y[k] + g_x_0_xyyy_yyyzzz[k];

                g_x_0_xyyyy_yzzzz[k] = -g_x_0_xyyy_yzzzz[k] * cd_y[k] + g_x_0_xyyy_yyzzzz[k];

                g_x_0_xyyyy_zzzzz[k] = -g_x_0_xyyy_zzzzz[k] * cd_y[k] + g_x_0_xyyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyz_xxxxx, g_x_0_xyyyz_xxxxy, g_x_0_xyyyz_xxxxz, g_x_0_xyyyz_xxxyy, g_x_0_xyyyz_xxxyz, g_x_0_xyyyz_xxxzz, g_x_0_xyyyz_xxyyy, g_x_0_xyyyz_xxyyz, g_x_0_xyyyz_xxyzz, g_x_0_xyyyz_xxzzz, g_x_0_xyyyz_xyyyy, g_x_0_xyyyz_xyyyz, g_x_0_xyyyz_xyyzz, g_x_0_xyyyz_xyzzz, g_x_0_xyyyz_xzzzz, g_x_0_xyyyz_yyyyy, g_x_0_xyyyz_yyyyz, g_x_0_xyyyz_yyyzz, g_x_0_xyyyz_yyzzz, g_x_0_xyyyz_yzzzz, g_x_0_xyyyz_zzzzz, g_x_0_xyyz_xxxxx, g_x_0_xyyz_xxxxxy, g_x_0_xyyz_xxxxy, g_x_0_xyyz_xxxxyy, g_x_0_xyyz_xxxxyz, g_x_0_xyyz_xxxxz, g_x_0_xyyz_xxxyy, g_x_0_xyyz_xxxyyy, g_x_0_xyyz_xxxyyz, g_x_0_xyyz_xxxyz, g_x_0_xyyz_xxxyzz, g_x_0_xyyz_xxxzz, g_x_0_xyyz_xxyyy, g_x_0_xyyz_xxyyyy, g_x_0_xyyz_xxyyyz, g_x_0_xyyz_xxyyz, g_x_0_xyyz_xxyyzz, g_x_0_xyyz_xxyzz, g_x_0_xyyz_xxyzzz, g_x_0_xyyz_xxzzz, g_x_0_xyyz_xyyyy, g_x_0_xyyz_xyyyyy, g_x_0_xyyz_xyyyyz, g_x_0_xyyz_xyyyz, g_x_0_xyyz_xyyyzz, g_x_0_xyyz_xyyzz, g_x_0_xyyz_xyyzzz, g_x_0_xyyz_xyzzz, g_x_0_xyyz_xyzzzz, g_x_0_xyyz_xzzzz, g_x_0_xyyz_yyyyy, g_x_0_xyyz_yyyyyy, g_x_0_xyyz_yyyyyz, g_x_0_xyyz_yyyyz, g_x_0_xyyz_yyyyzz, g_x_0_xyyz_yyyzz, g_x_0_xyyz_yyyzzz, g_x_0_xyyz_yyzzz, g_x_0_xyyz_yyzzzz, g_x_0_xyyz_yzzzz, g_x_0_xyyz_yzzzzz, g_x_0_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xxxxx[k] = -g_x_0_xyyz_xxxxx[k] * cd_y[k] + g_x_0_xyyz_xxxxxy[k];

                g_x_0_xyyyz_xxxxy[k] = -g_x_0_xyyz_xxxxy[k] * cd_y[k] + g_x_0_xyyz_xxxxyy[k];

                g_x_0_xyyyz_xxxxz[k] = -g_x_0_xyyz_xxxxz[k] * cd_y[k] + g_x_0_xyyz_xxxxyz[k];

                g_x_0_xyyyz_xxxyy[k] = -g_x_0_xyyz_xxxyy[k] * cd_y[k] + g_x_0_xyyz_xxxyyy[k];

                g_x_0_xyyyz_xxxyz[k] = -g_x_0_xyyz_xxxyz[k] * cd_y[k] + g_x_0_xyyz_xxxyyz[k];

                g_x_0_xyyyz_xxxzz[k] = -g_x_0_xyyz_xxxzz[k] * cd_y[k] + g_x_0_xyyz_xxxyzz[k];

                g_x_0_xyyyz_xxyyy[k] = -g_x_0_xyyz_xxyyy[k] * cd_y[k] + g_x_0_xyyz_xxyyyy[k];

                g_x_0_xyyyz_xxyyz[k] = -g_x_0_xyyz_xxyyz[k] * cd_y[k] + g_x_0_xyyz_xxyyyz[k];

                g_x_0_xyyyz_xxyzz[k] = -g_x_0_xyyz_xxyzz[k] * cd_y[k] + g_x_0_xyyz_xxyyzz[k];

                g_x_0_xyyyz_xxzzz[k] = -g_x_0_xyyz_xxzzz[k] * cd_y[k] + g_x_0_xyyz_xxyzzz[k];

                g_x_0_xyyyz_xyyyy[k] = -g_x_0_xyyz_xyyyy[k] * cd_y[k] + g_x_0_xyyz_xyyyyy[k];

                g_x_0_xyyyz_xyyyz[k] = -g_x_0_xyyz_xyyyz[k] * cd_y[k] + g_x_0_xyyz_xyyyyz[k];

                g_x_0_xyyyz_xyyzz[k] = -g_x_0_xyyz_xyyzz[k] * cd_y[k] + g_x_0_xyyz_xyyyzz[k];

                g_x_0_xyyyz_xyzzz[k] = -g_x_0_xyyz_xyzzz[k] * cd_y[k] + g_x_0_xyyz_xyyzzz[k];

                g_x_0_xyyyz_xzzzz[k] = -g_x_0_xyyz_xzzzz[k] * cd_y[k] + g_x_0_xyyz_xyzzzz[k];

                g_x_0_xyyyz_yyyyy[k] = -g_x_0_xyyz_yyyyy[k] * cd_y[k] + g_x_0_xyyz_yyyyyy[k];

                g_x_0_xyyyz_yyyyz[k] = -g_x_0_xyyz_yyyyz[k] * cd_y[k] + g_x_0_xyyz_yyyyyz[k];

                g_x_0_xyyyz_yyyzz[k] = -g_x_0_xyyz_yyyzz[k] * cd_y[k] + g_x_0_xyyz_yyyyzz[k];

                g_x_0_xyyyz_yyzzz[k] = -g_x_0_xyyz_yyzzz[k] * cd_y[k] + g_x_0_xyyz_yyyzzz[k];

                g_x_0_xyyyz_yzzzz[k] = -g_x_0_xyyz_yzzzz[k] * cd_y[k] + g_x_0_xyyz_yyzzzz[k];

                g_x_0_xyyyz_zzzzz[k] = -g_x_0_xyyz_zzzzz[k] * cd_y[k] + g_x_0_xyyz_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzz_xxxxx, g_x_0_xyyzz_xxxxy, g_x_0_xyyzz_xxxxz, g_x_0_xyyzz_xxxyy, g_x_0_xyyzz_xxxyz, g_x_0_xyyzz_xxxzz, g_x_0_xyyzz_xxyyy, g_x_0_xyyzz_xxyyz, g_x_0_xyyzz_xxyzz, g_x_0_xyyzz_xxzzz, g_x_0_xyyzz_xyyyy, g_x_0_xyyzz_xyyyz, g_x_0_xyyzz_xyyzz, g_x_0_xyyzz_xyzzz, g_x_0_xyyzz_xzzzz, g_x_0_xyyzz_yyyyy, g_x_0_xyyzz_yyyyz, g_x_0_xyyzz_yyyzz, g_x_0_xyyzz_yyzzz, g_x_0_xyyzz_yzzzz, g_x_0_xyyzz_zzzzz, g_x_0_xyzz_xxxxx, g_x_0_xyzz_xxxxxy, g_x_0_xyzz_xxxxy, g_x_0_xyzz_xxxxyy, g_x_0_xyzz_xxxxyz, g_x_0_xyzz_xxxxz, g_x_0_xyzz_xxxyy, g_x_0_xyzz_xxxyyy, g_x_0_xyzz_xxxyyz, g_x_0_xyzz_xxxyz, g_x_0_xyzz_xxxyzz, g_x_0_xyzz_xxxzz, g_x_0_xyzz_xxyyy, g_x_0_xyzz_xxyyyy, g_x_0_xyzz_xxyyyz, g_x_0_xyzz_xxyyz, g_x_0_xyzz_xxyyzz, g_x_0_xyzz_xxyzz, g_x_0_xyzz_xxyzzz, g_x_0_xyzz_xxzzz, g_x_0_xyzz_xyyyy, g_x_0_xyzz_xyyyyy, g_x_0_xyzz_xyyyyz, g_x_0_xyzz_xyyyz, g_x_0_xyzz_xyyyzz, g_x_0_xyzz_xyyzz, g_x_0_xyzz_xyyzzz, g_x_0_xyzz_xyzzz, g_x_0_xyzz_xyzzzz, g_x_0_xyzz_xzzzz, g_x_0_xyzz_yyyyy, g_x_0_xyzz_yyyyyy, g_x_0_xyzz_yyyyyz, g_x_0_xyzz_yyyyz, g_x_0_xyzz_yyyyzz, g_x_0_xyzz_yyyzz, g_x_0_xyzz_yyyzzz, g_x_0_xyzz_yyzzz, g_x_0_xyzz_yyzzzz, g_x_0_xyzz_yzzzz, g_x_0_xyzz_yzzzzz, g_x_0_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xxxxx[k] = -g_x_0_xyzz_xxxxx[k] * cd_y[k] + g_x_0_xyzz_xxxxxy[k];

                g_x_0_xyyzz_xxxxy[k] = -g_x_0_xyzz_xxxxy[k] * cd_y[k] + g_x_0_xyzz_xxxxyy[k];

                g_x_0_xyyzz_xxxxz[k] = -g_x_0_xyzz_xxxxz[k] * cd_y[k] + g_x_0_xyzz_xxxxyz[k];

                g_x_0_xyyzz_xxxyy[k] = -g_x_0_xyzz_xxxyy[k] * cd_y[k] + g_x_0_xyzz_xxxyyy[k];

                g_x_0_xyyzz_xxxyz[k] = -g_x_0_xyzz_xxxyz[k] * cd_y[k] + g_x_0_xyzz_xxxyyz[k];

                g_x_0_xyyzz_xxxzz[k] = -g_x_0_xyzz_xxxzz[k] * cd_y[k] + g_x_0_xyzz_xxxyzz[k];

                g_x_0_xyyzz_xxyyy[k] = -g_x_0_xyzz_xxyyy[k] * cd_y[k] + g_x_0_xyzz_xxyyyy[k];

                g_x_0_xyyzz_xxyyz[k] = -g_x_0_xyzz_xxyyz[k] * cd_y[k] + g_x_0_xyzz_xxyyyz[k];

                g_x_0_xyyzz_xxyzz[k] = -g_x_0_xyzz_xxyzz[k] * cd_y[k] + g_x_0_xyzz_xxyyzz[k];

                g_x_0_xyyzz_xxzzz[k] = -g_x_0_xyzz_xxzzz[k] * cd_y[k] + g_x_0_xyzz_xxyzzz[k];

                g_x_0_xyyzz_xyyyy[k] = -g_x_0_xyzz_xyyyy[k] * cd_y[k] + g_x_0_xyzz_xyyyyy[k];

                g_x_0_xyyzz_xyyyz[k] = -g_x_0_xyzz_xyyyz[k] * cd_y[k] + g_x_0_xyzz_xyyyyz[k];

                g_x_0_xyyzz_xyyzz[k] = -g_x_0_xyzz_xyyzz[k] * cd_y[k] + g_x_0_xyzz_xyyyzz[k];

                g_x_0_xyyzz_xyzzz[k] = -g_x_0_xyzz_xyzzz[k] * cd_y[k] + g_x_0_xyzz_xyyzzz[k];

                g_x_0_xyyzz_xzzzz[k] = -g_x_0_xyzz_xzzzz[k] * cd_y[k] + g_x_0_xyzz_xyzzzz[k];

                g_x_0_xyyzz_yyyyy[k] = -g_x_0_xyzz_yyyyy[k] * cd_y[k] + g_x_0_xyzz_yyyyyy[k];

                g_x_0_xyyzz_yyyyz[k] = -g_x_0_xyzz_yyyyz[k] * cd_y[k] + g_x_0_xyzz_yyyyyz[k];

                g_x_0_xyyzz_yyyzz[k] = -g_x_0_xyzz_yyyzz[k] * cd_y[k] + g_x_0_xyzz_yyyyzz[k];

                g_x_0_xyyzz_yyzzz[k] = -g_x_0_xyzz_yyzzz[k] * cd_y[k] + g_x_0_xyzz_yyyzzz[k];

                g_x_0_xyyzz_yzzzz[k] = -g_x_0_xyzz_yzzzz[k] * cd_y[k] + g_x_0_xyzz_yyzzzz[k];

                g_x_0_xyyzz_zzzzz[k] = -g_x_0_xyzz_zzzzz[k] * cd_y[k] + g_x_0_xyzz_yzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzz_xxxxx, g_x_0_xyzzz_xxxxy, g_x_0_xyzzz_xxxxz, g_x_0_xyzzz_xxxyy, g_x_0_xyzzz_xxxyz, g_x_0_xyzzz_xxxzz, g_x_0_xyzzz_xxyyy, g_x_0_xyzzz_xxyyz, g_x_0_xyzzz_xxyzz, g_x_0_xyzzz_xxzzz, g_x_0_xyzzz_xyyyy, g_x_0_xyzzz_xyyyz, g_x_0_xyzzz_xyyzz, g_x_0_xyzzz_xyzzz, g_x_0_xyzzz_xzzzz, g_x_0_xyzzz_yyyyy, g_x_0_xyzzz_yyyyz, g_x_0_xyzzz_yyyzz, g_x_0_xyzzz_yyzzz, g_x_0_xyzzz_yzzzz, g_x_0_xyzzz_zzzzz, g_x_0_xzzz_xxxxx, g_x_0_xzzz_xxxxxy, g_x_0_xzzz_xxxxy, g_x_0_xzzz_xxxxyy, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxz, g_x_0_xzzz_xxxyy, g_x_0_xzzz_xxxyyy, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxzz, g_x_0_xzzz_xxyyy, g_x_0_xzzz_xxyyyy, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxzzz, g_x_0_xzzz_xyyyy, g_x_0_xzzz_xyyyyy, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xzzzz, g_x_0_xzzz_yyyyy, g_x_0_xzzz_yyyyyy, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xxxxx[k] = -g_x_0_xzzz_xxxxx[k] * cd_y[k] + g_x_0_xzzz_xxxxxy[k];

                g_x_0_xyzzz_xxxxy[k] = -g_x_0_xzzz_xxxxy[k] * cd_y[k] + g_x_0_xzzz_xxxxyy[k];

                g_x_0_xyzzz_xxxxz[k] = -g_x_0_xzzz_xxxxz[k] * cd_y[k] + g_x_0_xzzz_xxxxyz[k];

                g_x_0_xyzzz_xxxyy[k] = -g_x_0_xzzz_xxxyy[k] * cd_y[k] + g_x_0_xzzz_xxxyyy[k];

                g_x_0_xyzzz_xxxyz[k] = -g_x_0_xzzz_xxxyz[k] * cd_y[k] + g_x_0_xzzz_xxxyyz[k];

                g_x_0_xyzzz_xxxzz[k] = -g_x_0_xzzz_xxxzz[k] * cd_y[k] + g_x_0_xzzz_xxxyzz[k];

                g_x_0_xyzzz_xxyyy[k] = -g_x_0_xzzz_xxyyy[k] * cd_y[k] + g_x_0_xzzz_xxyyyy[k];

                g_x_0_xyzzz_xxyyz[k] = -g_x_0_xzzz_xxyyz[k] * cd_y[k] + g_x_0_xzzz_xxyyyz[k];

                g_x_0_xyzzz_xxyzz[k] = -g_x_0_xzzz_xxyzz[k] * cd_y[k] + g_x_0_xzzz_xxyyzz[k];

                g_x_0_xyzzz_xxzzz[k] = -g_x_0_xzzz_xxzzz[k] * cd_y[k] + g_x_0_xzzz_xxyzzz[k];

                g_x_0_xyzzz_xyyyy[k] = -g_x_0_xzzz_xyyyy[k] * cd_y[k] + g_x_0_xzzz_xyyyyy[k];

                g_x_0_xyzzz_xyyyz[k] = -g_x_0_xzzz_xyyyz[k] * cd_y[k] + g_x_0_xzzz_xyyyyz[k];

                g_x_0_xyzzz_xyyzz[k] = -g_x_0_xzzz_xyyzz[k] * cd_y[k] + g_x_0_xzzz_xyyyzz[k];

                g_x_0_xyzzz_xyzzz[k] = -g_x_0_xzzz_xyzzz[k] * cd_y[k] + g_x_0_xzzz_xyyzzz[k];

                g_x_0_xyzzz_xzzzz[k] = -g_x_0_xzzz_xzzzz[k] * cd_y[k] + g_x_0_xzzz_xyzzzz[k];

                g_x_0_xyzzz_yyyyy[k] = -g_x_0_xzzz_yyyyy[k] * cd_y[k] + g_x_0_xzzz_yyyyyy[k];

                g_x_0_xyzzz_yyyyz[k] = -g_x_0_xzzz_yyyyz[k] * cd_y[k] + g_x_0_xzzz_yyyyyz[k];

                g_x_0_xyzzz_yyyzz[k] = -g_x_0_xzzz_yyyzz[k] * cd_y[k] + g_x_0_xzzz_yyyyzz[k];

                g_x_0_xyzzz_yyzzz[k] = -g_x_0_xzzz_yyzzz[k] * cd_y[k] + g_x_0_xzzz_yyyzzz[k];

                g_x_0_xyzzz_yzzzz[k] = -g_x_0_xzzz_yzzzz[k] * cd_y[k] + g_x_0_xzzz_yyzzzz[k];

                g_x_0_xyzzz_zzzzz[k] = -g_x_0_xzzz_zzzzz[k] * cd_y[k] + g_x_0_xzzz_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_x_0_xzzz_xxxxx, g_x_0_xzzz_xxxxxz, g_x_0_xzzz_xxxxy, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxz, g_x_0_xzzz_xxxxzz, g_x_0_xzzz_xxxyy, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxzz, g_x_0_xzzz_xxxzzz, g_x_0_xzzz_xxyyy, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxzzz, g_x_0_xzzz_xxzzzz, g_x_0_xzzz_xyyyy, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xzzzz, g_x_0_xzzz_xzzzzz, g_x_0_xzzz_yyyyy, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_zzzzz, g_x_0_xzzz_zzzzzz, g_x_0_xzzzz_xxxxx, g_x_0_xzzzz_xxxxy, g_x_0_xzzzz_xxxxz, g_x_0_xzzzz_xxxyy, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxzz, g_x_0_xzzzz_xxyyy, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxzzz, g_x_0_xzzzz_xyyyy, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xzzzz, g_x_0_xzzzz_yyyyy, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xxxxx[k] = -g_x_0_xzzz_xxxxx[k] * cd_z[k] + g_x_0_xzzz_xxxxxz[k];

                g_x_0_xzzzz_xxxxy[k] = -g_x_0_xzzz_xxxxy[k] * cd_z[k] + g_x_0_xzzz_xxxxyz[k];

                g_x_0_xzzzz_xxxxz[k] = -g_x_0_xzzz_xxxxz[k] * cd_z[k] + g_x_0_xzzz_xxxxzz[k];

                g_x_0_xzzzz_xxxyy[k] = -g_x_0_xzzz_xxxyy[k] * cd_z[k] + g_x_0_xzzz_xxxyyz[k];

                g_x_0_xzzzz_xxxyz[k] = -g_x_0_xzzz_xxxyz[k] * cd_z[k] + g_x_0_xzzz_xxxyzz[k];

                g_x_0_xzzzz_xxxzz[k] = -g_x_0_xzzz_xxxzz[k] * cd_z[k] + g_x_0_xzzz_xxxzzz[k];

                g_x_0_xzzzz_xxyyy[k] = -g_x_0_xzzz_xxyyy[k] * cd_z[k] + g_x_0_xzzz_xxyyyz[k];

                g_x_0_xzzzz_xxyyz[k] = -g_x_0_xzzz_xxyyz[k] * cd_z[k] + g_x_0_xzzz_xxyyzz[k];

                g_x_0_xzzzz_xxyzz[k] = -g_x_0_xzzz_xxyzz[k] * cd_z[k] + g_x_0_xzzz_xxyzzz[k];

                g_x_0_xzzzz_xxzzz[k] = -g_x_0_xzzz_xxzzz[k] * cd_z[k] + g_x_0_xzzz_xxzzzz[k];

                g_x_0_xzzzz_xyyyy[k] = -g_x_0_xzzz_xyyyy[k] * cd_z[k] + g_x_0_xzzz_xyyyyz[k];

                g_x_0_xzzzz_xyyyz[k] = -g_x_0_xzzz_xyyyz[k] * cd_z[k] + g_x_0_xzzz_xyyyzz[k];

                g_x_0_xzzzz_xyyzz[k] = -g_x_0_xzzz_xyyzz[k] * cd_z[k] + g_x_0_xzzz_xyyzzz[k];

                g_x_0_xzzzz_xyzzz[k] = -g_x_0_xzzz_xyzzz[k] * cd_z[k] + g_x_0_xzzz_xyzzzz[k];

                g_x_0_xzzzz_xzzzz[k] = -g_x_0_xzzz_xzzzz[k] * cd_z[k] + g_x_0_xzzz_xzzzzz[k];

                g_x_0_xzzzz_yyyyy[k] = -g_x_0_xzzz_yyyyy[k] * cd_z[k] + g_x_0_xzzz_yyyyyz[k];

                g_x_0_xzzzz_yyyyz[k] = -g_x_0_xzzz_yyyyz[k] * cd_z[k] + g_x_0_xzzz_yyyyzz[k];

                g_x_0_xzzzz_yyyzz[k] = -g_x_0_xzzz_yyyzz[k] * cd_z[k] + g_x_0_xzzz_yyyzzz[k];

                g_x_0_xzzzz_yyzzz[k] = -g_x_0_xzzz_yyzzz[k] * cd_z[k] + g_x_0_xzzz_yyzzzz[k];

                g_x_0_xzzzz_yzzzz[k] = -g_x_0_xzzz_yzzzz[k] * cd_z[k] + g_x_0_xzzz_yzzzzz[k];

                g_x_0_xzzzz_zzzzz[k] = -g_x_0_xzzz_zzzzz[k] * cd_z[k] + g_x_0_xzzz_zzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_y, g_x_0_yyyy_xxxxx, g_x_0_yyyy_xxxxxy, g_x_0_yyyy_xxxxy, g_x_0_yyyy_xxxxyy, g_x_0_yyyy_xxxxyz, g_x_0_yyyy_xxxxz, g_x_0_yyyy_xxxyy, g_x_0_yyyy_xxxyyy, g_x_0_yyyy_xxxyyz, g_x_0_yyyy_xxxyz, g_x_0_yyyy_xxxyzz, g_x_0_yyyy_xxxzz, g_x_0_yyyy_xxyyy, g_x_0_yyyy_xxyyyy, g_x_0_yyyy_xxyyyz, g_x_0_yyyy_xxyyz, g_x_0_yyyy_xxyyzz, g_x_0_yyyy_xxyzz, g_x_0_yyyy_xxyzzz, g_x_0_yyyy_xxzzz, g_x_0_yyyy_xyyyy, g_x_0_yyyy_xyyyyy, g_x_0_yyyy_xyyyyz, g_x_0_yyyy_xyyyz, g_x_0_yyyy_xyyyzz, g_x_0_yyyy_xyyzz, g_x_0_yyyy_xyyzzz, g_x_0_yyyy_xyzzz, g_x_0_yyyy_xyzzzz, g_x_0_yyyy_xzzzz, g_x_0_yyyy_yyyyy, g_x_0_yyyy_yyyyyy, g_x_0_yyyy_yyyyyz, g_x_0_yyyy_yyyyz, g_x_0_yyyy_yyyyzz, g_x_0_yyyy_yyyzz, g_x_0_yyyy_yyyzzz, g_x_0_yyyy_yyzzz, g_x_0_yyyy_yyzzzz, g_x_0_yyyy_yzzzz, g_x_0_yyyy_yzzzzz, g_x_0_yyyy_zzzzz, g_x_0_yyyyy_xxxxx, g_x_0_yyyyy_xxxxy, g_x_0_yyyyy_xxxxz, g_x_0_yyyyy_xxxyy, g_x_0_yyyyy_xxxyz, g_x_0_yyyyy_xxxzz, g_x_0_yyyyy_xxyyy, g_x_0_yyyyy_xxyyz, g_x_0_yyyyy_xxyzz, g_x_0_yyyyy_xxzzz, g_x_0_yyyyy_xyyyy, g_x_0_yyyyy_xyyyz, g_x_0_yyyyy_xyyzz, g_x_0_yyyyy_xyzzz, g_x_0_yyyyy_xzzzz, g_x_0_yyyyy_yyyyy, g_x_0_yyyyy_yyyyz, g_x_0_yyyyy_yyyzz, g_x_0_yyyyy_yyzzz, g_x_0_yyyyy_yzzzz, g_x_0_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xxxxx[k] = -g_x_0_yyyy_xxxxx[k] * cd_y[k] + g_x_0_yyyy_xxxxxy[k];

                g_x_0_yyyyy_xxxxy[k] = -g_x_0_yyyy_xxxxy[k] * cd_y[k] + g_x_0_yyyy_xxxxyy[k];

                g_x_0_yyyyy_xxxxz[k] = -g_x_0_yyyy_xxxxz[k] * cd_y[k] + g_x_0_yyyy_xxxxyz[k];

                g_x_0_yyyyy_xxxyy[k] = -g_x_0_yyyy_xxxyy[k] * cd_y[k] + g_x_0_yyyy_xxxyyy[k];

                g_x_0_yyyyy_xxxyz[k] = -g_x_0_yyyy_xxxyz[k] * cd_y[k] + g_x_0_yyyy_xxxyyz[k];

                g_x_0_yyyyy_xxxzz[k] = -g_x_0_yyyy_xxxzz[k] * cd_y[k] + g_x_0_yyyy_xxxyzz[k];

                g_x_0_yyyyy_xxyyy[k] = -g_x_0_yyyy_xxyyy[k] * cd_y[k] + g_x_0_yyyy_xxyyyy[k];

                g_x_0_yyyyy_xxyyz[k] = -g_x_0_yyyy_xxyyz[k] * cd_y[k] + g_x_0_yyyy_xxyyyz[k];

                g_x_0_yyyyy_xxyzz[k] = -g_x_0_yyyy_xxyzz[k] * cd_y[k] + g_x_0_yyyy_xxyyzz[k];

                g_x_0_yyyyy_xxzzz[k] = -g_x_0_yyyy_xxzzz[k] * cd_y[k] + g_x_0_yyyy_xxyzzz[k];

                g_x_0_yyyyy_xyyyy[k] = -g_x_0_yyyy_xyyyy[k] * cd_y[k] + g_x_0_yyyy_xyyyyy[k];

                g_x_0_yyyyy_xyyyz[k] = -g_x_0_yyyy_xyyyz[k] * cd_y[k] + g_x_0_yyyy_xyyyyz[k];

                g_x_0_yyyyy_xyyzz[k] = -g_x_0_yyyy_xyyzz[k] * cd_y[k] + g_x_0_yyyy_xyyyzz[k];

                g_x_0_yyyyy_xyzzz[k] = -g_x_0_yyyy_xyzzz[k] * cd_y[k] + g_x_0_yyyy_xyyzzz[k];

                g_x_0_yyyyy_xzzzz[k] = -g_x_0_yyyy_xzzzz[k] * cd_y[k] + g_x_0_yyyy_xyzzzz[k];

                g_x_0_yyyyy_yyyyy[k] = -g_x_0_yyyy_yyyyy[k] * cd_y[k] + g_x_0_yyyy_yyyyyy[k];

                g_x_0_yyyyy_yyyyz[k] = -g_x_0_yyyy_yyyyz[k] * cd_y[k] + g_x_0_yyyy_yyyyyz[k];

                g_x_0_yyyyy_yyyzz[k] = -g_x_0_yyyy_yyyzz[k] * cd_y[k] + g_x_0_yyyy_yyyyzz[k];

                g_x_0_yyyyy_yyzzz[k] = -g_x_0_yyyy_yyzzz[k] * cd_y[k] + g_x_0_yyyy_yyyzzz[k];

                g_x_0_yyyyy_yzzzz[k] = -g_x_0_yyyy_yzzzz[k] * cd_y[k] + g_x_0_yyyy_yyzzzz[k];

                g_x_0_yyyyy_zzzzz[k] = -g_x_0_yyyy_zzzzz[k] * cd_y[k] + g_x_0_yyyy_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 356);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyz_xxxxx, g_x_0_yyyyz_xxxxy, g_x_0_yyyyz_xxxxz, g_x_0_yyyyz_xxxyy, g_x_0_yyyyz_xxxyz, g_x_0_yyyyz_xxxzz, g_x_0_yyyyz_xxyyy, g_x_0_yyyyz_xxyyz, g_x_0_yyyyz_xxyzz, g_x_0_yyyyz_xxzzz, g_x_0_yyyyz_xyyyy, g_x_0_yyyyz_xyyyz, g_x_0_yyyyz_xyyzz, g_x_0_yyyyz_xyzzz, g_x_0_yyyyz_xzzzz, g_x_0_yyyyz_yyyyy, g_x_0_yyyyz_yyyyz, g_x_0_yyyyz_yyyzz, g_x_0_yyyyz_yyzzz, g_x_0_yyyyz_yzzzz, g_x_0_yyyyz_zzzzz, g_x_0_yyyz_xxxxx, g_x_0_yyyz_xxxxxy, g_x_0_yyyz_xxxxy, g_x_0_yyyz_xxxxyy, g_x_0_yyyz_xxxxyz, g_x_0_yyyz_xxxxz, g_x_0_yyyz_xxxyy, g_x_0_yyyz_xxxyyy, g_x_0_yyyz_xxxyyz, g_x_0_yyyz_xxxyz, g_x_0_yyyz_xxxyzz, g_x_0_yyyz_xxxzz, g_x_0_yyyz_xxyyy, g_x_0_yyyz_xxyyyy, g_x_0_yyyz_xxyyyz, g_x_0_yyyz_xxyyz, g_x_0_yyyz_xxyyzz, g_x_0_yyyz_xxyzz, g_x_0_yyyz_xxyzzz, g_x_0_yyyz_xxzzz, g_x_0_yyyz_xyyyy, g_x_0_yyyz_xyyyyy, g_x_0_yyyz_xyyyyz, g_x_0_yyyz_xyyyz, g_x_0_yyyz_xyyyzz, g_x_0_yyyz_xyyzz, g_x_0_yyyz_xyyzzz, g_x_0_yyyz_xyzzz, g_x_0_yyyz_xyzzzz, g_x_0_yyyz_xzzzz, g_x_0_yyyz_yyyyy, g_x_0_yyyz_yyyyyy, g_x_0_yyyz_yyyyyz, g_x_0_yyyz_yyyyz, g_x_0_yyyz_yyyyzz, g_x_0_yyyz_yyyzz, g_x_0_yyyz_yyyzzz, g_x_0_yyyz_yyzzz, g_x_0_yyyz_yyzzzz, g_x_0_yyyz_yzzzz, g_x_0_yyyz_yzzzzz, g_x_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xxxxx[k] = -g_x_0_yyyz_xxxxx[k] * cd_y[k] + g_x_0_yyyz_xxxxxy[k];

                g_x_0_yyyyz_xxxxy[k] = -g_x_0_yyyz_xxxxy[k] * cd_y[k] + g_x_0_yyyz_xxxxyy[k];

                g_x_0_yyyyz_xxxxz[k] = -g_x_0_yyyz_xxxxz[k] * cd_y[k] + g_x_0_yyyz_xxxxyz[k];

                g_x_0_yyyyz_xxxyy[k] = -g_x_0_yyyz_xxxyy[k] * cd_y[k] + g_x_0_yyyz_xxxyyy[k];

                g_x_0_yyyyz_xxxyz[k] = -g_x_0_yyyz_xxxyz[k] * cd_y[k] + g_x_0_yyyz_xxxyyz[k];

                g_x_0_yyyyz_xxxzz[k] = -g_x_0_yyyz_xxxzz[k] * cd_y[k] + g_x_0_yyyz_xxxyzz[k];

                g_x_0_yyyyz_xxyyy[k] = -g_x_0_yyyz_xxyyy[k] * cd_y[k] + g_x_0_yyyz_xxyyyy[k];

                g_x_0_yyyyz_xxyyz[k] = -g_x_0_yyyz_xxyyz[k] * cd_y[k] + g_x_0_yyyz_xxyyyz[k];

                g_x_0_yyyyz_xxyzz[k] = -g_x_0_yyyz_xxyzz[k] * cd_y[k] + g_x_0_yyyz_xxyyzz[k];

                g_x_0_yyyyz_xxzzz[k] = -g_x_0_yyyz_xxzzz[k] * cd_y[k] + g_x_0_yyyz_xxyzzz[k];

                g_x_0_yyyyz_xyyyy[k] = -g_x_0_yyyz_xyyyy[k] * cd_y[k] + g_x_0_yyyz_xyyyyy[k];

                g_x_0_yyyyz_xyyyz[k] = -g_x_0_yyyz_xyyyz[k] * cd_y[k] + g_x_0_yyyz_xyyyyz[k];

                g_x_0_yyyyz_xyyzz[k] = -g_x_0_yyyz_xyyzz[k] * cd_y[k] + g_x_0_yyyz_xyyyzz[k];

                g_x_0_yyyyz_xyzzz[k] = -g_x_0_yyyz_xyzzz[k] * cd_y[k] + g_x_0_yyyz_xyyzzz[k];

                g_x_0_yyyyz_xzzzz[k] = -g_x_0_yyyz_xzzzz[k] * cd_y[k] + g_x_0_yyyz_xyzzzz[k];

                g_x_0_yyyyz_yyyyy[k] = -g_x_0_yyyz_yyyyy[k] * cd_y[k] + g_x_0_yyyz_yyyyyy[k];

                g_x_0_yyyyz_yyyyz[k] = -g_x_0_yyyz_yyyyz[k] * cd_y[k] + g_x_0_yyyz_yyyyyz[k];

                g_x_0_yyyyz_yyyzz[k] = -g_x_0_yyyz_yyyzz[k] * cd_y[k] + g_x_0_yyyz_yyyyzz[k];

                g_x_0_yyyyz_yyzzz[k] = -g_x_0_yyyz_yyzzz[k] * cd_y[k] + g_x_0_yyyz_yyyzzz[k];

                g_x_0_yyyyz_yzzzz[k] = -g_x_0_yyyz_yzzzz[k] * cd_y[k] + g_x_0_yyyz_yyzzzz[k];

                g_x_0_yyyyz_zzzzz[k] = -g_x_0_yyyz_zzzzz[k] * cd_y[k] + g_x_0_yyyz_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 363);

            auto g_x_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 377);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzz_xxxxx, g_x_0_yyyzz_xxxxy, g_x_0_yyyzz_xxxxz, g_x_0_yyyzz_xxxyy, g_x_0_yyyzz_xxxyz, g_x_0_yyyzz_xxxzz, g_x_0_yyyzz_xxyyy, g_x_0_yyyzz_xxyyz, g_x_0_yyyzz_xxyzz, g_x_0_yyyzz_xxzzz, g_x_0_yyyzz_xyyyy, g_x_0_yyyzz_xyyyz, g_x_0_yyyzz_xyyzz, g_x_0_yyyzz_xyzzz, g_x_0_yyyzz_xzzzz, g_x_0_yyyzz_yyyyy, g_x_0_yyyzz_yyyyz, g_x_0_yyyzz_yyyzz, g_x_0_yyyzz_yyzzz, g_x_0_yyyzz_yzzzz, g_x_0_yyyzz_zzzzz, g_x_0_yyzz_xxxxx, g_x_0_yyzz_xxxxxy, g_x_0_yyzz_xxxxy, g_x_0_yyzz_xxxxyy, g_x_0_yyzz_xxxxyz, g_x_0_yyzz_xxxxz, g_x_0_yyzz_xxxyy, g_x_0_yyzz_xxxyyy, g_x_0_yyzz_xxxyyz, g_x_0_yyzz_xxxyz, g_x_0_yyzz_xxxyzz, g_x_0_yyzz_xxxzz, g_x_0_yyzz_xxyyy, g_x_0_yyzz_xxyyyy, g_x_0_yyzz_xxyyyz, g_x_0_yyzz_xxyyz, g_x_0_yyzz_xxyyzz, g_x_0_yyzz_xxyzz, g_x_0_yyzz_xxyzzz, g_x_0_yyzz_xxzzz, g_x_0_yyzz_xyyyy, g_x_0_yyzz_xyyyyy, g_x_0_yyzz_xyyyyz, g_x_0_yyzz_xyyyz, g_x_0_yyzz_xyyyzz, g_x_0_yyzz_xyyzz, g_x_0_yyzz_xyyzzz, g_x_0_yyzz_xyzzz, g_x_0_yyzz_xyzzzz, g_x_0_yyzz_xzzzz, g_x_0_yyzz_yyyyy, g_x_0_yyzz_yyyyyy, g_x_0_yyzz_yyyyyz, g_x_0_yyzz_yyyyz, g_x_0_yyzz_yyyyzz, g_x_0_yyzz_yyyzz, g_x_0_yyzz_yyyzzz, g_x_0_yyzz_yyzzz, g_x_0_yyzz_yyzzzz, g_x_0_yyzz_yzzzz, g_x_0_yyzz_yzzzzz, g_x_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xxxxx[k] = -g_x_0_yyzz_xxxxx[k] * cd_y[k] + g_x_0_yyzz_xxxxxy[k];

                g_x_0_yyyzz_xxxxy[k] = -g_x_0_yyzz_xxxxy[k] * cd_y[k] + g_x_0_yyzz_xxxxyy[k];

                g_x_0_yyyzz_xxxxz[k] = -g_x_0_yyzz_xxxxz[k] * cd_y[k] + g_x_0_yyzz_xxxxyz[k];

                g_x_0_yyyzz_xxxyy[k] = -g_x_0_yyzz_xxxyy[k] * cd_y[k] + g_x_0_yyzz_xxxyyy[k];

                g_x_0_yyyzz_xxxyz[k] = -g_x_0_yyzz_xxxyz[k] * cd_y[k] + g_x_0_yyzz_xxxyyz[k];

                g_x_0_yyyzz_xxxzz[k] = -g_x_0_yyzz_xxxzz[k] * cd_y[k] + g_x_0_yyzz_xxxyzz[k];

                g_x_0_yyyzz_xxyyy[k] = -g_x_0_yyzz_xxyyy[k] * cd_y[k] + g_x_0_yyzz_xxyyyy[k];

                g_x_0_yyyzz_xxyyz[k] = -g_x_0_yyzz_xxyyz[k] * cd_y[k] + g_x_0_yyzz_xxyyyz[k];

                g_x_0_yyyzz_xxyzz[k] = -g_x_0_yyzz_xxyzz[k] * cd_y[k] + g_x_0_yyzz_xxyyzz[k];

                g_x_0_yyyzz_xxzzz[k] = -g_x_0_yyzz_xxzzz[k] * cd_y[k] + g_x_0_yyzz_xxyzzz[k];

                g_x_0_yyyzz_xyyyy[k] = -g_x_0_yyzz_xyyyy[k] * cd_y[k] + g_x_0_yyzz_xyyyyy[k];

                g_x_0_yyyzz_xyyyz[k] = -g_x_0_yyzz_xyyyz[k] * cd_y[k] + g_x_0_yyzz_xyyyyz[k];

                g_x_0_yyyzz_xyyzz[k] = -g_x_0_yyzz_xyyzz[k] * cd_y[k] + g_x_0_yyzz_xyyyzz[k];

                g_x_0_yyyzz_xyzzz[k] = -g_x_0_yyzz_xyzzz[k] * cd_y[k] + g_x_0_yyzz_xyyzzz[k];

                g_x_0_yyyzz_xzzzz[k] = -g_x_0_yyzz_xzzzz[k] * cd_y[k] + g_x_0_yyzz_xyzzzz[k];

                g_x_0_yyyzz_yyyyy[k] = -g_x_0_yyzz_yyyyy[k] * cd_y[k] + g_x_0_yyzz_yyyyyy[k];

                g_x_0_yyyzz_yyyyz[k] = -g_x_0_yyzz_yyyyz[k] * cd_y[k] + g_x_0_yyzz_yyyyyz[k];

                g_x_0_yyyzz_yyyzz[k] = -g_x_0_yyzz_yyyzz[k] * cd_y[k] + g_x_0_yyzz_yyyyzz[k];

                g_x_0_yyyzz_yyzzz[k] = -g_x_0_yyzz_yyzzz[k] * cd_y[k] + g_x_0_yyzz_yyyzzz[k];

                g_x_0_yyyzz_yzzzz[k] = -g_x_0_yyzz_yzzzz[k] * cd_y[k] + g_x_0_yyzz_yyzzzz[k];

                g_x_0_yyyzz_zzzzz[k] = -g_x_0_yyzz_zzzzz[k] * cd_y[k] + g_x_0_yyzz_yzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 391);

            auto g_x_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 398);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzz_xxxxx, g_x_0_yyzzz_xxxxy, g_x_0_yyzzz_xxxxz, g_x_0_yyzzz_xxxyy, g_x_0_yyzzz_xxxyz, g_x_0_yyzzz_xxxzz, g_x_0_yyzzz_xxyyy, g_x_0_yyzzz_xxyyz, g_x_0_yyzzz_xxyzz, g_x_0_yyzzz_xxzzz, g_x_0_yyzzz_xyyyy, g_x_0_yyzzz_xyyyz, g_x_0_yyzzz_xyyzz, g_x_0_yyzzz_xyzzz, g_x_0_yyzzz_xzzzz, g_x_0_yyzzz_yyyyy, g_x_0_yyzzz_yyyyz, g_x_0_yyzzz_yyyzz, g_x_0_yyzzz_yyzzz, g_x_0_yyzzz_yzzzz, g_x_0_yyzzz_zzzzz, g_x_0_yzzz_xxxxx, g_x_0_yzzz_xxxxxy, g_x_0_yzzz_xxxxy, g_x_0_yzzz_xxxxyy, g_x_0_yzzz_xxxxyz, g_x_0_yzzz_xxxxz, g_x_0_yzzz_xxxyy, g_x_0_yzzz_xxxyyy, g_x_0_yzzz_xxxyyz, g_x_0_yzzz_xxxyz, g_x_0_yzzz_xxxyzz, g_x_0_yzzz_xxxzz, g_x_0_yzzz_xxyyy, g_x_0_yzzz_xxyyyy, g_x_0_yzzz_xxyyyz, g_x_0_yzzz_xxyyz, g_x_0_yzzz_xxyyzz, g_x_0_yzzz_xxyzz, g_x_0_yzzz_xxyzzz, g_x_0_yzzz_xxzzz, g_x_0_yzzz_xyyyy, g_x_0_yzzz_xyyyyy, g_x_0_yzzz_xyyyyz, g_x_0_yzzz_xyyyz, g_x_0_yzzz_xyyyzz, g_x_0_yzzz_xyyzz, g_x_0_yzzz_xyyzzz, g_x_0_yzzz_xyzzz, g_x_0_yzzz_xyzzzz, g_x_0_yzzz_xzzzz, g_x_0_yzzz_yyyyy, g_x_0_yzzz_yyyyyy, g_x_0_yzzz_yyyyyz, g_x_0_yzzz_yyyyz, g_x_0_yzzz_yyyyzz, g_x_0_yzzz_yyyzz, g_x_0_yzzz_yyyzzz, g_x_0_yzzz_yyzzz, g_x_0_yzzz_yyzzzz, g_x_0_yzzz_yzzzz, g_x_0_yzzz_yzzzzz, g_x_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xxxxx[k] = -g_x_0_yzzz_xxxxx[k] * cd_y[k] + g_x_0_yzzz_xxxxxy[k];

                g_x_0_yyzzz_xxxxy[k] = -g_x_0_yzzz_xxxxy[k] * cd_y[k] + g_x_0_yzzz_xxxxyy[k];

                g_x_0_yyzzz_xxxxz[k] = -g_x_0_yzzz_xxxxz[k] * cd_y[k] + g_x_0_yzzz_xxxxyz[k];

                g_x_0_yyzzz_xxxyy[k] = -g_x_0_yzzz_xxxyy[k] * cd_y[k] + g_x_0_yzzz_xxxyyy[k];

                g_x_0_yyzzz_xxxyz[k] = -g_x_0_yzzz_xxxyz[k] * cd_y[k] + g_x_0_yzzz_xxxyyz[k];

                g_x_0_yyzzz_xxxzz[k] = -g_x_0_yzzz_xxxzz[k] * cd_y[k] + g_x_0_yzzz_xxxyzz[k];

                g_x_0_yyzzz_xxyyy[k] = -g_x_0_yzzz_xxyyy[k] * cd_y[k] + g_x_0_yzzz_xxyyyy[k];

                g_x_0_yyzzz_xxyyz[k] = -g_x_0_yzzz_xxyyz[k] * cd_y[k] + g_x_0_yzzz_xxyyyz[k];

                g_x_0_yyzzz_xxyzz[k] = -g_x_0_yzzz_xxyzz[k] * cd_y[k] + g_x_0_yzzz_xxyyzz[k];

                g_x_0_yyzzz_xxzzz[k] = -g_x_0_yzzz_xxzzz[k] * cd_y[k] + g_x_0_yzzz_xxyzzz[k];

                g_x_0_yyzzz_xyyyy[k] = -g_x_0_yzzz_xyyyy[k] * cd_y[k] + g_x_0_yzzz_xyyyyy[k];

                g_x_0_yyzzz_xyyyz[k] = -g_x_0_yzzz_xyyyz[k] * cd_y[k] + g_x_0_yzzz_xyyyyz[k];

                g_x_0_yyzzz_xyyzz[k] = -g_x_0_yzzz_xyyzz[k] * cd_y[k] + g_x_0_yzzz_xyyyzz[k];

                g_x_0_yyzzz_xyzzz[k] = -g_x_0_yzzz_xyzzz[k] * cd_y[k] + g_x_0_yzzz_xyyzzz[k];

                g_x_0_yyzzz_xzzzz[k] = -g_x_0_yzzz_xzzzz[k] * cd_y[k] + g_x_0_yzzz_xyzzzz[k];

                g_x_0_yyzzz_yyyyy[k] = -g_x_0_yzzz_yyyyy[k] * cd_y[k] + g_x_0_yzzz_yyyyyy[k];

                g_x_0_yyzzz_yyyyz[k] = -g_x_0_yzzz_yyyyz[k] * cd_y[k] + g_x_0_yzzz_yyyyyz[k];

                g_x_0_yyzzz_yyyzz[k] = -g_x_0_yzzz_yyyzz[k] * cd_y[k] + g_x_0_yzzz_yyyyzz[k];

                g_x_0_yyzzz_yyzzz[k] = -g_x_0_yzzz_yyzzz[k] * cd_y[k] + g_x_0_yzzz_yyyzzz[k];

                g_x_0_yyzzz_yzzzz[k] = -g_x_0_yzzz_yzzzz[k] * cd_y[k] + g_x_0_yzzz_yyzzzz[k];

                g_x_0_yyzzz_zzzzz[k] = -g_x_0_yzzz_zzzzz[k] * cd_y[k] + g_x_0_yzzz_yzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzz_xxxxx, g_x_0_yzzzz_xxxxy, g_x_0_yzzzz_xxxxz, g_x_0_yzzzz_xxxyy, g_x_0_yzzzz_xxxyz, g_x_0_yzzzz_xxxzz, g_x_0_yzzzz_xxyyy, g_x_0_yzzzz_xxyyz, g_x_0_yzzzz_xxyzz, g_x_0_yzzzz_xxzzz, g_x_0_yzzzz_xyyyy, g_x_0_yzzzz_xyyyz, g_x_0_yzzzz_xyyzz, g_x_0_yzzzz_xyzzz, g_x_0_yzzzz_xzzzz, g_x_0_yzzzz_yyyyy, g_x_0_yzzzz_yyyyz, g_x_0_yzzzz_yyyzz, g_x_0_yzzzz_yyzzz, g_x_0_yzzzz_yzzzz, g_x_0_yzzzz_zzzzz, g_x_0_zzzz_xxxxx, g_x_0_zzzz_xxxxxy, g_x_0_zzzz_xxxxy, g_x_0_zzzz_xxxxyy, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxz, g_x_0_zzzz_xxxyy, g_x_0_zzzz_xxxyyy, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxzz, g_x_0_zzzz_xxyyy, g_x_0_zzzz_xxyyyy, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxzzz, g_x_0_zzzz_xyyyy, g_x_0_zzzz_xyyyyy, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xzzzz, g_x_0_zzzz_yyyyy, g_x_0_zzzz_yyyyyy, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xxxxx[k] = -g_x_0_zzzz_xxxxx[k] * cd_y[k] + g_x_0_zzzz_xxxxxy[k];

                g_x_0_yzzzz_xxxxy[k] = -g_x_0_zzzz_xxxxy[k] * cd_y[k] + g_x_0_zzzz_xxxxyy[k];

                g_x_0_yzzzz_xxxxz[k] = -g_x_0_zzzz_xxxxz[k] * cd_y[k] + g_x_0_zzzz_xxxxyz[k];

                g_x_0_yzzzz_xxxyy[k] = -g_x_0_zzzz_xxxyy[k] * cd_y[k] + g_x_0_zzzz_xxxyyy[k];

                g_x_0_yzzzz_xxxyz[k] = -g_x_0_zzzz_xxxyz[k] * cd_y[k] + g_x_0_zzzz_xxxyyz[k];

                g_x_0_yzzzz_xxxzz[k] = -g_x_0_zzzz_xxxzz[k] * cd_y[k] + g_x_0_zzzz_xxxyzz[k];

                g_x_0_yzzzz_xxyyy[k] = -g_x_0_zzzz_xxyyy[k] * cd_y[k] + g_x_0_zzzz_xxyyyy[k];

                g_x_0_yzzzz_xxyyz[k] = -g_x_0_zzzz_xxyyz[k] * cd_y[k] + g_x_0_zzzz_xxyyyz[k];

                g_x_0_yzzzz_xxyzz[k] = -g_x_0_zzzz_xxyzz[k] * cd_y[k] + g_x_0_zzzz_xxyyzz[k];

                g_x_0_yzzzz_xxzzz[k] = -g_x_0_zzzz_xxzzz[k] * cd_y[k] + g_x_0_zzzz_xxyzzz[k];

                g_x_0_yzzzz_xyyyy[k] = -g_x_0_zzzz_xyyyy[k] * cd_y[k] + g_x_0_zzzz_xyyyyy[k];

                g_x_0_yzzzz_xyyyz[k] = -g_x_0_zzzz_xyyyz[k] * cd_y[k] + g_x_0_zzzz_xyyyyz[k];

                g_x_0_yzzzz_xyyzz[k] = -g_x_0_zzzz_xyyzz[k] * cd_y[k] + g_x_0_zzzz_xyyyzz[k];

                g_x_0_yzzzz_xyzzz[k] = -g_x_0_zzzz_xyzzz[k] * cd_y[k] + g_x_0_zzzz_xyyzzz[k];

                g_x_0_yzzzz_xzzzz[k] = -g_x_0_zzzz_xzzzz[k] * cd_y[k] + g_x_0_zzzz_xyzzzz[k];

                g_x_0_yzzzz_yyyyy[k] = -g_x_0_zzzz_yyyyy[k] * cd_y[k] + g_x_0_zzzz_yyyyyy[k];

                g_x_0_yzzzz_yyyyz[k] = -g_x_0_zzzz_yyyyz[k] * cd_y[k] + g_x_0_zzzz_yyyyyz[k];

                g_x_0_yzzzz_yyyzz[k] = -g_x_0_zzzz_yyyzz[k] * cd_y[k] + g_x_0_zzzz_yyyyzz[k];

                g_x_0_yzzzz_yyzzz[k] = -g_x_0_zzzz_yyzzz[k] * cd_y[k] + g_x_0_zzzz_yyyzzz[k];

                g_x_0_yzzzz_yzzzz[k] = -g_x_0_zzzz_yzzzz[k] * cd_y[k] + g_x_0_zzzz_yyzzzz[k];

                g_x_0_yzzzz_zzzzz[k] = -g_x_0_zzzz_zzzzz[k] * cd_y[k] + g_x_0_zzzz_yzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 420);

            auto g_x_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 421);

            auto g_x_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 422);

            auto g_x_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 423);

            auto g_x_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 424);

            auto g_x_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 425);

            auto g_x_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 426);

            auto g_x_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 427);

            auto g_x_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 428);

            auto g_x_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 429);

            auto g_x_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 430);

            auto g_x_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 431);

            auto g_x_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 432);

            auto g_x_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 433);

            auto g_x_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 434);

            auto g_x_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 435);

            auto g_x_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 436);

            auto g_x_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 437);

            auto g_x_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 438);

            auto g_x_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 439);

            auto g_x_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 0 * acomps * bcomps + 440);

            #pragma omp simd aligned(cd_z, g_x_0_zzzz_xxxxx, g_x_0_zzzz_xxxxxz, g_x_0_zzzz_xxxxy, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxz, g_x_0_zzzz_xxxxzz, g_x_0_zzzz_xxxyy, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxzz, g_x_0_zzzz_xxxzzz, g_x_0_zzzz_xxyyy, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxzzz, g_x_0_zzzz_xxzzzz, g_x_0_zzzz_xyyyy, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xzzzz, g_x_0_zzzz_xzzzzz, g_x_0_zzzz_yyyyy, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_zzzzz, g_x_0_zzzz_zzzzzz, g_x_0_zzzzz_xxxxx, g_x_0_zzzzz_xxxxy, g_x_0_zzzzz_xxxxz, g_x_0_zzzzz_xxxyy, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxzz, g_x_0_zzzzz_xxyyy, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxzzz, g_x_0_zzzzz_xyyyy, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xzzzz, g_x_0_zzzzz_yyyyy, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xxxxx[k] = -g_x_0_zzzz_xxxxx[k] * cd_z[k] + g_x_0_zzzz_xxxxxz[k];

                g_x_0_zzzzz_xxxxy[k] = -g_x_0_zzzz_xxxxy[k] * cd_z[k] + g_x_0_zzzz_xxxxyz[k];

                g_x_0_zzzzz_xxxxz[k] = -g_x_0_zzzz_xxxxz[k] * cd_z[k] + g_x_0_zzzz_xxxxzz[k];

                g_x_0_zzzzz_xxxyy[k] = -g_x_0_zzzz_xxxyy[k] * cd_z[k] + g_x_0_zzzz_xxxyyz[k];

                g_x_0_zzzzz_xxxyz[k] = -g_x_0_zzzz_xxxyz[k] * cd_z[k] + g_x_0_zzzz_xxxyzz[k];

                g_x_0_zzzzz_xxxzz[k] = -g_x_0_zzzz_xxxzz[k] * cd_z[k] + g_x_0_zzzz_xxxzzz[k];

                g_x_0_zzzzz_xxyyy[k] = -g_x_0_zzzz_xxyyy[k] * cd_z[k] + g_x_0_zzzz_xxyyyz[k];

                g_x_0_zzzzz_xxyyz[k] = -g_x_0_zzzz_xxyyz[k] * cd_z[k] + g_x_0_zzzz_xxyyzz[k];

                g_x_0_zzzzz_xxyzz[k] = -g_x_0_zzzz_xxyzz[k] * cd_z[k] + g_x_0_zzzz_xxyzzz[k];

                g_x_0_zzzzz_xxzzz[k] = -g_x_0_zzzz_xxzzz[k] * cd_z[k] + g_x_0_zzzz_xxzzzz[k];

                g_x_0_zzzzz_xyyyy[k] = -g_x_0_zzzz_xyyyy[k] * cd_z[k] + g_x_0_zzzz_xyyyyz[k];

                g_x_0_zzzzz_xyyyz[k] = -g_x_0_zzzz_xyyyz[k] * cd_z[k] + g_x_0_zzzz_xyyyzz[k];

                g_x_0_zzzzz_xyyzz[k] = -g_x_0_zzzz_xyyzz[k] * cd_z[k] + g_x_0_zzzz_xyyzzz[k];

                g_x_0_zzzzz_xyzzz[k] = -g_x_0_zzzz_xyzzz[k] * cd_z[k] + g_x_0_zzzz_xyzzzz[k];

                g_x_0_zzzzz_xzzzz[k] = -g_x_0_zzzz_xzzzz[k] * cd_z[k] + g_x_0_zzzz_xzzzzz[k];

                g_x_0_zzzzz_yyyyy[k] = -g_x_0_zzzz_yyyyy[k] * cd_z[k] + g_x_0_zzzz_yyyyyz[k];

                g_x_0_zzzzz_yyyyz[k] = -g_x_0_zzzz_yyyyz[k] * cd_z[k] + g_x_0_zzzz_yyyyzz[k];

                g_x_0_zzzzz_yyyzz[k] = -g_x_0_zzzz_yyyzz[k] * cd_z[k] + g_x_0_zzzz_yyyzzz[k];

                g_x_0_zzzzz_yyzzz[k] = -g_x_0_zzzz_yyzzz[k] * cd_z[k] + g_x_0_zzzz_yyzzzz[k];

                g_x_0_zzzzz_yzzzz[k] = -g_x_0_zzzz_yzzzz[k] * cd_z[k] + g_x_0_zzzz_yzzzzz[k];

                g_x_0_zzzzz_zzzzz[k] = -g_x_0_zzzz_zzzzz[k] * cd_z[k] + g_x_0_zzzz_zzzzzz[k];
            }
            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 5);

            auto g_y_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 6);

            auto g_y_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 7);

            auto g_y_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 8);

            auto g_y_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 9);

            auto g_y_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 10);

            auto g_y_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 11);

            auto g_y_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 12);

            auto g_y_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 13);

            auto g_y_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 14);

            auto g_y_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 15);

            auto g_y_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 16);

            auto g_y_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 17);

            auto g_y_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 18);

            auto g_y_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 19);

            auto g_y_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_y_0_xxxx_xxxxx, g_y_0_xxxx_xxxxxx, g_y_0_xxxx_xxxxxy, g_y_0_xxxx_xxxxxz, g_y_0_xxxx_xxxxy, g_y_0_xxxx_xxxxyy, g_y_0_xxxx_xxxxyz, g_y_0_xxxx_xxxxz, g_y_0_xxxx_xxxxzz, g_y_0_xxxx_xxxyy, g_y_0_xxxx_xxxyyy, g_y_0_xxxx_xxxyyz, g_y_0_xxxx_xxxyz, g_y_0_xxxx_xxxyzz, g_y_0_xxxx_xxxzz, g_y_0_xxxx_xxxzzz, g_y_0_xxxx_xxyyy, g_y_0_xxxx_xxyyyy, g_y_0_xxxx_xxyyyz, g_y_0_xxxx_xxyyz, g_y_0_xxxx_xxyyzz, g_y_0_xxxx_xxyzz, g_y_0_xxxx_xxyzzz, g_y_0_xxxx_xxzzz, g_y_0_xxxx_xxzzzz, g_y_0_xxxx_xyyyy, g_y_0_xxxx_xyyyyy, g_y_0_xxxx_xyyyyz, g_y_0_xxxx_xyyyz, g_y_0_xxxx_xyyyzz, g_y_0_xxxx_xyyzz, g_y_0_xxxx_xyyzzz, g_y_0_xxxx_xyzzz, g_y_0_xxxx_xyzzzz, g_y_0_xxxx_xzzzz, g_y_0_xxxx_xzzzzz, g_y_0_xxxx_yyyyy, g_y_0_xxxx_yyyyz, g_y_0_xxxx_yyyzz, g_y_0_xxxx_yyzzz, g_y_0_xxxx_yzzzz, g_y_0_xxxx_zzzzz, g_y_0_xxxxx_xxxxx, g_y_0_xxxxx_xxxxy, g_y_0_xxxxx_xxxxz, g_y_0_xxxxx_xxxyy, g_y_0_xxxxx_xxxyz, g_y_0_xxxxx_xxxzz, g_y_0_xxxxx_xxyyy, g_y_0_xxxxx_xxyyz, g_y_0_xxxxx_xxyzz, g_y_0_xxxxx_xxzzz, g_y_0_xxxxx_xyyyy, g_y_0_xxxxx_xyyyz, g_y_0_xxxxx_xyyzz, g_y_0_xxxxx_xyzzz, g_y_0_xxxxx_xzzzz, g_y_0_xxxxx_yyyyy, g_y_0_xxxxx_yyyyz, g_y_0_xxxxx_yyyzz, g_y_0_xxxxx_yyzzz, g_y_0_xxxxx_yzzzz, g_y_0_xxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xxxxx[k] = -g_y_0_xxxx_xxxxx[k] * cd_x[k] + g_y_0_xxxx_xxxxxx[k];

                g_y_0_xxxxx_xxxxy[k] = -g_y_0_xxxx_xxxxy[k] * cd_x[k] + g_y_0_xxxx_xxxxxy[k];

                g_y_0_xxxxx_xxxxz[k] = -g_y_0_xxxx_xxxxz[k] * cd_x[k] + g_y_0_xxxx_xxxxxz[k];

                g_y_0_xxxxx_xxxyy[k] = -g_y_0_xxxx_xxxyy[k] * cd_x[k] + g_y_0_xxxx_xxxxyy[k];

                g_y_0_xxxxx_xxxyz[k] = -g_y_0_xxxx_xxxyz[k] * cd_x[k] + g_y_0_xxxx_xxxxyz[k];

                g_y_0_xxxxx_xxxzz[k] = -g_y_0_xxxx_xxxzz[k] * cd_x[k] + g_y_0_xxxx_xxxxzz[k];

                g_y_0_xxxxx_xxyyy[k] = -g_y_0_xxxx_xxyyy[k] * cd_x[k] + g_y_0_xxxx_xxxyyy[k];

                g_y_0_xxxxx_xxyyz[k] = -g_y_0_xxxx_xxyyz[k] * cd_x[k] + g_y_0_xxxx_xxxyyz[k];

                g_y_0_xxxxx_xxyzz[k] = -g_y_0_xxxx_xxyzz[k] * cd_x[k] + g_y_0_xxxx_xxxyzz[k];

                g_y_0_xxxxx_xxzzz[k] = -g_y_0_xxxx_xxzzz[k] * cd_x[k] + g_y_0_xxxx_xxxzzz[k];

                g_y_0_xxxxx_xyyyy[k] = -g_y_0_xxxx_xyyyy[k] * cd_x[k] + g_y_0_xxxx_xxyyyy[k];

                g_y_0_xxxxx_xyyyz[k] = -g_y_0_xxxx_xyyyz[k] * cd_x[k] + g_y_0_xxxx_xxyyyz[k];

                g_y_0_xxxxx_xyyzz[k] = -g_y_0_xxxx_xyyzz[k] * cd_x[k] + g_y_0_xxxx_xxyyzz[k];

                g_y_0_xxxxx_xyzzz[k] = -g_y_0_xxxx_xyzzz[k] * cd_x[k] + g_y_0_xxxx_xxyzzz[k];

                g_y_0_xxxxx_xzzzz[k] = -g_y_0_xxxx_xzzzz[k] * cd_x[k] + g_y_0_xxxx_xxzzzz[k];

                g_y_0_xxxxx_yyyyy[k] = -g_y_0_xxxx_yyyyy[k] * cd_x[k] + g_y_0_xxxx_xyyyyy[k];

                g_y_0_xxxxx_yyyyz[k] = -g_y_0_xxxx_yyyyz[k] * cd_x[k] + g_y_0_xxxx_xyyyyz[k];

                g_y_0_xxxxx_yyyzz[k] = -g_y_0_xxxx_yyyzz[k] * cd_x[k] + g_y_0_xxxx_xyyyzz[k];

                g_y_0_xxxxx_yyzzz[k] = -g_y_0_xxxx_yyzzz[k] * cd_x[k] + g_y_0_xxxx_xyyzzz[k];

                g_y_0_xxxxx_yzzzz[k] = -g_y_0_xxxx_yzzzz[k] * cd_x[k] + g_y_0_xxxx_xyzzzz[k];

                g_y_0_xxxxx_zzzzz[k] = -g_y_0_xxxx_zzzzz[k] * cd_x[k] + g_y_0_xxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 21);

            auto g_y_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 22);

            auto g_y_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 23);

            auto g_y_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 24);

            auto g_y_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 25);

            auto g_y_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 26);

            auto g_y_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 27);

            auto g_y_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 28);

            auto g_y_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 29);

            auto g_y_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 30);

            auto g_y_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 31);

            auto g_y_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 32);

            auto g_y_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 33);

            auto g_y_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 34);

            auto g_y_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 35);

            auto g_y_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 36);

            auto g_y_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 37);

            auto g_y_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 38);

            auto g_y_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 39);

            auto g_y_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 40);

            auto g_y_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxy_xxxxx, g_y_0_xxxxy_xxxxy, g_y_0_xxxxy_xxxxz, g_y_0_xxxxy_xxxyy, g_y_0_xxxxy_xxxyz, g_y_0_xxxxy_xxxzz, g_y_0_xxxxy_xxyyy, g_y_0_xxxxy_xxyyz, g_y_0_xxxxy_xxyzz, g_y_0_xxxxy_xxzzz, g_y_0_xxxxy_xyyyy, g_y_0_xxxxy_xyyyz, g_y_0_xxxxy_xyyzz, g_y_0_xxxxy_xyzzz, g_y_0_xxxxy_xzzzz, g_y_0_xxxxy_yyyyy, g_y_0_xxxxy_yyyyz, g_y_0_xxxxy_yyyzz, g_y_0_xxxxy_yyzzz, g_y_0_xxxxy_yzzzz, g_y_0_xxxxy_zzzzz, g_y_0_xxxy_xxxxx, g_y_0_xxxy_xxxxxx, g_y_0_xxxy_xxxxxy, g_y_0_xxxy_xxxxxz, g_y_0_xxxy_xxxxy, g_y_0_xxxy_xxxxyy, g_y_0_xxxy_xxxxyz, g_y_0_xxxy_xxxxz, g_y_0_xxxy_xxxxzz, g_y_0_xxxy_xxxyy, g_y_0_xxxy_xxxyyy, g_y_0_xxxy_xxxyyz, g_y_0_xxxy_xxxyz, g_y_0_xxxy_xxxyzz, g_y_0_xxxy_xxxzz, g_y_0_xxxy_xxxzzz, g_y_0_xxxy_xxyyy, g_y_0_xxxy_xxyyyy, g_y_0_xxxy_xxyyyz, g_y_0_xxxy_xxyyz, g_y_0_xxxy_xxyyzz, g_y_0_xxxy_xxyzz, g_y_0_xxxy_xxyzzz, g_y_0_xxxy_xxzzz, g_y_0_xxxy_xxzzzz, g_y_0_xxxy_xyyyy, g_y_0_xxxy_xyyyyy, g_y_0_xxxy_xyyyyz, g_y_0_xxxy_xyyyz, g_y_0_xxxy_xyyyzz, g_y_0_xxxy_xyyzz, g_y_0_xxxy_xyyzzz, g_y_0_xxxy_xyzzz, g_y_0_xxxy_xyzzzz, g_y_0_xxxy_xzzzz, g_y_0_xxxy_xzzzzz, g_y_0_xxxy_yyyyy, g_y_0_xxxy_yyyyz, g_y_0_xxxy_yyyzz, g_y_0_xxxy_yyzzz, g_y_0_xxxy_yzzzz, g_y_0_xxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xxxxx[k] = -g_y_0_xxxy_xxxxx[k] * cd_x[k] + g_y_0_xxxy_xxxxxx[k];

                g_y_0_xxxxy_xxxxy[k] = -g_y_0_xxxy_xxxxy[k] * cd_x[k] + g_y_0_xxxy_xxxxxy[k];

                g_y_0_xxxxy_xxxxz[k] = -g_y_0_xxxy_xxxxz[k] * cd_x[k] + g_y_0_xxxy_xxxxxz[k];

                g_y_0_xxxxy_xxxyy[k] = -g_y_0_xxxy_xxxyy[k] * cd_x[k] + g_y_0_xxxy_xxxxyy[k];

                g_y_0_xxxxy_xxxyz[k] = -g_y_0_xxxy_xxxyz[k] * cd_x[k] + g_y_0_xxxy_xxxxyz[k];

                g_y_0_xxxxy_xxxzz[k] = -g_y_0_xxxy_xxxzz[k] * cd_x[k] + g_y_0_xxxy_xxxxzz[k];

                g_y_0_xxxxy_xxyyy[k] = -g_y_0_xxxy_xxyyy[k] * cd_x[k] + g_y_0_xxxy_xxxyyy[k];

                g_y_0_xxxxy_xxyyz[k] = -g_y_0_xxxy_xxyyz[k] * cd_x[k] + g_y_0_xxxy_xxxyyz[k];

                g_y_0_xxxxy_xxyzz[k] = -g_y_0_xxxy_xxyzz[k] * cd_x[k] + g_y_0_xxxy_xxxyzz[k];

                g_y_0_xxxxy_xxzzz[k] = -g_y_0_xxxy_xxzzz[k] * cd_x[k] + g_y_0_xxxy_xxxzzz[k];

                g_y_0_xxxxy_xyyyy[k] = -g_y_0_xxxy_xyyyy[k] * cd_x[k] + g_y_0_xxxy_xxyyyy[k];

                g_y_0_xxxxy_xyyyz[k] = -g_y_0_xxxy_xyyyz[k] * cd_x[k] + g_y_0_xxxy_xxyyyz[k];

                g_y_0_xxxxy_xyyzz[k] = -g_y_0_xxxy_xyyzz[k] * cd_x[k] + g_y_0_xxxy_xxyyzz[k];

                g_y_0_xxxxy_xyzzz[k] = -g_y_0_xxxy_xyzzz[k] * cd_x[k] + g_y_0_xxxy_xxyzzz[k];

                g_y_0_xxxxy_xzzzz[k] = -g_y_0_xxxy_xzzzz[k] * cd_x[k] + g_y_0_xxxy_xxzzzz[k];

                g_y_0_xxxxy_yyyyy[k] = -g_y_0_xxxy_yyyyy[k] * cd_x[k] + g_y_0_xxxy_xyyyyy[k];

                g_y_0_xxxxy_yyyyz[k] = -g_y_0_xxxy_yyyyz[k] * cd_x[k] + g_y_0_xxxy_xyyyyz[k];

                g_y_0_xxxxy_yyyzz[k] = -g_y_0_xxxy_yyyzz[k] * cd_x[k] + g_y_0_xxxy_xyyyzz[k];

                g_y_0_xxxxy_yyzzz[k] = -g_y_0_xxxy_yyzzz[k] * cd_x[k] + g_y_0_xxxy_xyyzzz[k];

                g_y_0_xxxxy_yzzzz[k] = -g_y_0_xxxy_yzzzz[k] * cd_x[k] + g_y_0_xxxy_xyzzzz[k];

                g_y_0_xxxxy_zzzzz[k] = -g_y_0_xxxy_zzzzz[k] * cd_x[k] + g_y_0_xxxy_xzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 42);

            auto g_y_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 43);

            auto g_y_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 44);

            auto g_y_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 45);

            auto g_y_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 46);

            auto g_y_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 47);

            auto g_y_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 48);

            auto g_y_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 49);

            auto g_y_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 50);

            auto g_y_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 51);

            auto g_y_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 52);

            auto g_y_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 53);

            auto g_y_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 54);

            auto g_y_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 55);

            auto g_y_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 56);

            auto g_y_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 57);

            auto g_y_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 58);

            auto g_y_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 59);

            auto g_y_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 60);

            auto g_y_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 61);

            auto g_y_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxz_xxxxx, g_y_0_xxxxz_xxxxy, g_y_0_xxxxz_xxxxz, g_y_0_xxxxz_xxxyy, g_y_0_xxxxz_xxxyz, g_y_0_xxxxz_xxxzz, g_y_0_xxxxz_xxyyy, g_y_0_xxxxz_xxyyz, g_y_0_xxxxz_xxyzz, g_y_0_xxxxz_xxzzz, g_y_0_xxxxz_xyyyy, g_y_0_xxxxz_xyyyz, g_y_0_xxxxz_xyyzz, g_y_0_xxxxz_xyzzz, g_y_0_xxxxz_xzzzz, g_y_0_xxxxz_yyyyy, g_y_0_xxxxz_yyyyz, g_y_0_xxxxz_yyyzz, g_y_0_xxxxz_yyzzz, g_y_0_xxxxz_yzzzz, g_y_0_xxxxz_zzzzz, g_y_0_xxxz_xxxxx, g_y_0_xxxz_xxxxxx, g_y_0_xxxz_xxxxxy, g_y_0_xxxz_xxxxxz, g_y_0_xxxz_xxxxy, g_y_0_xxxz_xxxxyy, g_y_0_xxxz_xxxxyz, g_y_0_xxxz_xxxxz, g_y_0_xxxz_xxxxzz, g_y_0_xxxz_xxxyy, g_y_0_xxxz_xxxyyy, g_y_0_xxxz_xxxyyz, g_y_0_xxxz_xxxyz, g_y_0_xxxz_xxxyzz, g_y_0_xxxz_xxxzz, g_y_0_xxxz_xxxzzz, g_y_0_xxxz_xxyyy, g_y_0_xxxz_xxyyyy, g_y_0_xxxz_xxyyyz, g_y_0_xxxz_xxyyz, g_y_0_xxxz_xxyyzz, g_y_0_xxxz_xxyzz, g_y_0_xxxz_xxyzzz, g_y_0_xxxz_xxzzz, g_y_0_xxxz_xxzzzz, g_y_0_xxxz_xyyyy, g_y_0_xxxz_xyyyyy, g_y_0_xxxz_xyyyyz, g_y_0_xxxz_xyyyz, g_y_0_xxxz_xyyyzz, g_y_0_xxxz_xyyzz, g_y_0_xxxz_xyyzzz, g_y_0_xxxz_xyzzz, g_y_0_xxxz_xyzzzz, g_y_0_xxxz_xzzzz, g_y_0_xxxz_xzzzzz, g_y_0_xxxz_yyyyy, g_y_0_xxxz_yyyyz, g_y_0_xxxz_yyyzz, g_y_0_xxxz_yyzzz, g_y_0_xxxz_yzzzz, g_y_0_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xxxxx[k] = -g_y_0_xxxz_xxxxx[k] * cd_x[k] + g_y_0_xxxz_xxxxxx[k];

                g_y_0_xxxxz_xxxxy[k] = -g_y_0_xxxz_xxxxy[k] * cd_x[k] + g_y_0_xxxz_xxxxxy[k];

                g_y_0_xxxxz_xxxxz[k] = -g_y_0_xxxz_xxxxz[k] * cd_x[k] + g_y_0_xxxz_xxxxxz[k];

                g_y_0_xxxxz_xxxyy[k] = -g_y_0_xxxz_xxxyy[k] * cd_x[k] + g_y_0_xxxz_xxxxyy[k];

                g_y_0_xxxxz_xxxyz[k] = -g_y_0_xxxz_xxxyz[k] * cd_x[k] + g_y_0_xxxz_xxxxyz[k];

                g_y_0_xxxxz_xxxzz[k] = -g_y_0_xxxz_xxxzz[k] * cd_x[k] + g_y_0_xxxz_xxxxzz[k];

                g_y_0_xxxxz_xxyyy[k] = -g_y_0_xxxz_xxyyy[k] * cd_x[k] + g_y_0_xxxz_xxxyyy[k];

                g_y_0_xxxxz_xxyyz[k] = -g_y_0_xxxz_xxyyz[k] * cd_x[k] + g_y_0_xxxz_xxxyyz[k];

                g_y_0_xxxxz_xxyzz[k] = -g_y_0_xxxz_xxyzz[k] * cd_x[k] + g_y_0_xxxz_xxxyzz[k];

                g_y_0_xxxxz_xxzzz[k] = -g_y_0_xxxz_xxzzz[k] * cd_x[k] + g_y_0_xxxz_xxxzzz[k];

                g_y_0_xxxxz_xyyyy[k] = -g_y_0_xxxz_xyyyy[k] * cd_x[k] + g_y_0_xxxz_xxyyyy[k];

                g_y_0_xxxxz_xyyyz[k] = -g_y_0_xxxz_xyyyz[k] * cd_x[k] + g_y_0_xxxz_xxyyyz[k];

                g_y_0_xxxxz_xyyzz[k] = -g_y_0_xxxz_xyyzz[k] * cd_x[k] + g_y_0_xxxz_xxyyzz[k];

                g_y_0_xxxxz_xyzzz[k] = -g_y_0_xxxz_xyzzz[k] * cd_x[k] + g_y_0_xxxz_xxyzzz[k];

                g_y_0_xxxxz_xzzzz[k] = -g_y_0_xxxz_xzzzz[k] * cd_x[k] + g_y_0_xxxz_xxzzzz[k];

                g_y_0_xxxxz_yyyyy[k] = -g_y_0_xxxz_yyyyy[k] * cd_x[k] + g_y_0_xxxz_xyyyyy[k];

                g_y_0_xxxxz_yyyyz[k] = -g_y_0_xxxz_yyyyz[k] * cd_x[k] + g_y_0_xxxz_xyyyyz[k];

                g_y_0_xxxxz_yyyzz[k] = -g_y_0_xxxz_yyyzz[k] * cd_x[k] + g_y_0_xxxz_xyyyzz[k];

                g_y_0_xxxxz_yyzzz[k] = -g_y_0_xxxz_yyzzz[k] * cd_x[k] + g_y_0_xxxz_xyyzzz[k];

                g_y_0_xxxxz_yzzzz[k] = -g_y_0_xxxz_yzzzz[k] * cd_x[k] + g_y_0_xxxz_xyzzzz[k];

                g_y_0_xxxxz_zzzzz[k] = -g_y_0_xxxz_zzzzz[k] * cd_x[k] + g_y_0_xxxz_xzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 63);

            auto g_y_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 64);

            auto g_y_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 65);

            auto g_y_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 66);

            auto g_y_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 67);

            auto g_y_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 68);

            auto g_y_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 69);

            auto g_y_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 70);

            auto g_y_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 71);

            auto g_y_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 72);

            auto g_y_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 73);

            auto g_y_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 74);

            auto g_y_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 75);

            auto g_y_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 76);

            auto g_y_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 77);

            auto g_y_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 78);

            auto g_y_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 79);

            auto g_y_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 80);

            auto g_y_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 81);

            auto g_y_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 82);

            auto g_y_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyy_xxxxx, g_y_0_xxxyy_xxxxy, g_y_0_xxxyy_xxxxz, g_y_0_xxxyy_xxxyy, g_y_0_xxxyy_xxxyz, g_y_0_xxxyy_xxxzz, g_y_0_xxxyy_xxyyy, g_y_0_xxxyy_xxyyz, g_y_0_xxxyy_xxyzz, g_y_0_xxxyy_xxzzz, g_y_0_xxxyy_xyyyy, g_y_0_xxxyy_xyyyz, g_y_0_xxxyy_xyyzz, g_y_0_xxxyy_xyzzz, g_y_0_xxxyy_xzzzz, g_y_0_xxxyy_yyyyy, g_y_0_xxxyy_yyyyz, g_y_0_xxxyy_yyyzz, g_y_0_xxxyy_yyzzz, g_y_0_xxxyy_yzzzz, g_y_0_xxxyy_zzzzz, g_y_0_xxyy_xxxxx, g_y_0_xxyy_xxxxxx, g_y_0_xxyy_xxxxxy, g_y_0_xxyy_xxxxxz, g_y_0_xxyy_xxxxy, g_y_0_xxyy_xxxxyy, g_y_0_xxyy_xxxxyz, g_y_0_xxyy_xxxxz, g_y_0_xxyy_xxxxzz, g_y_0_xxyy_xxxyy, g_y_0_xxyy_xxxyyy, g_y_0_xxyy_xxxyyz, g_y_0_xxyy_xxxyz, g_y_0_xxyy_xxxyzz, g_y_0_xxyy_xxxzz, g_y_0_xxyy_xxxzzz, g_y_0_xxyy_xxyyy, g_y_0_xxyy_xxyyyy, g_y_0_xxyy_xxyyyz, g_y_0_xxyy_xxyyz, g_y_0_xxyy_xxyyzz, g_y_0_xxyy_xxyzz, g_y_0_xxyy_xxyzzz, g_y_0_xxyy_xxzzz, g_y_0_xxyy_xxzzzz, g_y_0_xxyy_xyyyy, g_y_0_xxyy_xyyyyy, g_y_0_xxyy_xyyyyz, g_y_0_xxyy_xyyyz, g_y_0_xxyy_xyyyzz, g_y_0_xxyy_xyyzz, g_y_0_xxyy_xyyzzz, g_y_0_xxyy_xyzzz, g_y_0_xxyy_xyzzzz, g_y_0_xxyy_xzzzz, g_y_0_xxyy_xzzzzz, g_y_0_xxyy_yyyyy, g_y_0_xxyy_yyyyz, g_y_0_xxyy_yyyzz, g_y_0_xxyy_yyzzz, g_y_0_xxyy_yzzzz, g_y_0_xxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xxxxx[k] = -g_y_0_xxyy_xxxxx[k] * cd_x[k] + g_y_0_xxyy_xxxxxx[k];

                g_y_0_xxxyy_xxxxy[k] = -g_y_0_xxyy_xxxxy[k] * cd_x[k] + g_y_0_xxyy_xxxxxy[k];

                g_y_0_xxxyy_xxxxz[k] = -g_y_0_xxyy_xxxxz[k] * cd_x[k] + g_y_0_xxyy_xxxxxz[k];

                g_y_0_xxxyy_xxxyy[k] = -g_y_0_xxyy_xxxyy[k] * cd_x[k] + g_y_0_xxyy_xxxxyy[k];

                g_y_0_xxxyy_xxxyz[k] = -g_y_0_xxyy_xxxyz[k] * cd_x[k] + g_y_0_xxyy_xxxxyz[k];

                g_y_0_xxxyy_xxxzz[k] = -g_y_0_xxyy_xxxzz[k] * cd_x[k] + g_y_0_xxyy_xxxxzz[k];

                g_y_0_xxxyy_xxyyy[k] = -g_y_0_xxyy_xxyyy[k] * cd_x[k] + g_y_0_xxyy_xxxyyy[k];

                g_y_0_xxxyy_xxyyz[k] = -g_y_0_xxyy_xxyyz[k] * cd_x[k] + g_y_0_xxyy_xxxyyz[k];

                g_y_0_xxxyy_xxyzz[k] = -g_y_0_xxyy_xxyzz[k] * cd_x[k] + g_y_0_xxyy_xxxyzz[k];

                g_y_0_xxxyy_xxzzz[k] = -g_y_0_xxyy_xxzzz[k] * cd_x[k] + g_y_0_xxyy_xxxzzz[k];

                g_y_0_xxxyy_xyyyy[k] = -g_y_0_xxyy_xyyyy[k] * cd_x[k] + g_y_0_xxyy_xxyyyy[k];

                g_y_0_xxxyy_xyyyz[k] = -g_y_0_xxyy_xyyyz[k] * cd_x[k] + g_y_0_xxyy_xxyyyz[k];

                g_y_0_xxxyy_xyyzz[k] = -g_y_0_xxyy_xyyzz[k] * cd_x[k] + g_y_0_xxyy_xxyyzz[k];

                g_y_0_xxxyy_xyzzz[k] = -g_y_0_xxyy_xyzzz[k] * cd_x[k] + g_y_0_xxyy_xxyzzz[k];

                g_y_0_xxxyy_xzzzz[k] = -g_y_0_xxyy_xzzzz[k] * cd_x[k] + g_y_0_xxyy_xxzzzz[k];

                g_y_0_xxxyy_yyyyy[k] = -g_y_0_xxyy_yyyyy[k] * cd_x[k] + g_y_0_xxyy_xyyyyy[k];

                g_y_0_xxxyy_yyyyz[k] = -g_y_0_xxyy_yyyyz[k] * cd_x[k] + g_y_0_xxyy_xyyyyz[k];

                g_y_0_xxxyy_yyyzz[k] = -g_y_0_xxyy_yyyzz[k] * cd_x[k] + g_y_0_xxyy_xyyyzz[k];

                g_y_0_xxxyy_yyzzz[k] = -g_y_0_xxyy_yyzzz[k] * cd_x[k] + g_y_0_xxyy_xyyzzz[k];

                g_y_0_xxxyy_yzzzz[k] = -g_y_0_xxyy_yzzzz[k] * cd_x[k] + g_y_0_xxyy_xyzzzz[k];

                g_y_0_xxxyy_zzzzz[k] = -g_y_0_xxyy_zzzzz[k] * cd_x[k] + g_y_0_xxyy_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 84);

            auto g_y_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 85);

            auto g_y_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 86);

            auto g_y_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 87);

            auto g_y_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 88);

            auto g_y_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 89);

            auto g_y_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 90);

            auto g_y_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 91);

            auto g_y_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 92);

            auto g_y_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 93);

            auto g_y_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 94);

            auto g_y_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 95);

            auto g_y_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 96);

            auto g_y_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 97);

            auto g_y_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 98);

            auto g_y_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 99);

            auto g_y_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 100);

            auto g_y_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 101);

            auto g_y_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 102);

            auto g_y_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 103);

            auto g_y_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyz_xxxxx, g_y_0_xxxyz_xxxxy, g_y_0_xxxyz_xxxxz, g_y_0_xxxyz_xxxyy, g_y_0_xxxyz_xxxyz, g_y_0_xxxyz_xxxzz, g_y_0_xxxyz_xxyyy, g_y_0_xxxyz_xxyyz, g_y_0_xxxyz_xxyzz, g_y_0_xxxyz_xxzzz, g_y_0_xxxyz_xyyyy, g_y_0_xxxyz_xyyyz, g_y_0_xxxyz_xyyzz, g_y_0_xxxyz_xyzzz, g_y_0_xxxyz_xzzzz, g_y_0_xxxyz_yyyyy, g_y_0_xxxyz_yyyyz, g_y_0_xxxyz_yyyzz, g_y_0_xxxyz_yyzzz, g_y_0_xxxyz_yzzzz, g_y_0_xxxyz_zzzzz, g_y_0_xxyz_xxxxx, g_y_0_xxyz_xxxxxx, g_y_0_xxyz_xxxxxy, g_y_0_xxyz_xxxxxz, g_y_0_xxyz_xxxxy, g_y_0_xxyz_xxxxyy, g_y_0_xxyz_xxxxyz, g_y_0_xxyz_xxxxz, g_y_0_xxyz_xxxxzz, g_y_0_xxyz_xxxyy, g_y_0_xxyz_xxxyyy, g_y_0_xxyz_xxxyyz, g_y_0_xxyz_xxxyz, g_y_0_xxyz_xxxyzz, g_y_0_xxyz_xxxzz, g_y_0_xxyz_xxxzzz, g_y_0_xxyz_xxyyy, g_y_0_xxyz_xxyyyy, g_y_0_xxyz_xxyyyz, g_y_0_xxyz_xxyyz, g_y_0_xxyz_xxyyzz, g_y_0_xxyz_xxyzz, g_y_0_xxyz_xxyzzz, g_y_0_xxyz_xxzzz, g_y_0_xxyz_xxzzzz, g_y_0_xxyz_xyyyy, g_y_0_xxyz_xyyyyy, g_y_0_xxyz_xyyyyz, g_y_0_xxyz_xyyyz, g_y_0_xxyz_xyyyzz, g_y_0_xxyz_xyyzz, g_y_0_xxyz_xyyzzz, g_y_0_xxyz_xyzzz, g_y_0_xxyz_xyzzzz, g_y_0_xxyz_xzzzz, g_y_0_xxyz_xzzzzz, g_y_0_xxyz_yyyyy, g_y_0_xxyz_yyyyz, g_y_0_xxyz_yyyzz, g_y_0_xxyz_yyzzz, g_y_0_xxyz_yzzzz, g_y_0_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xxxxx[k] = -g_y_0_xxyz_xxxxx[k] * cd_x[k] + g_y_0_xxyz_xxxxxx[k];

                g_y_0_xxxyz_xxxxy[k] = -g_y_0_xxyz_xxxxy[k] * cd_x[k] + g_y_0_xxyz_xxxxxy[k];

                g_y_0_xxxyz_xxxxz[k] = -g_y_0_xxyz_xxxxz[k] * cd_x[k] + g_y_0_xxyz_xxxxxz[k];

                g_y_0_xxxyz_xxxyy[k] = -g_y_0_xxyz_xxxyy[k] * cd_x[k] + g_y_0_xxyz_xxxxyy[k];

                g_y_0_xxxyz_xxxyz[k] = -g_y_0_xxyz_xxxyz[k] * cd_x[k] + g_y_0_xxyz_xxxxyz[k];

                g_y_0_xxxyz_xxxzz[k] = -g_y_0_xxyz_xxxzz[k] * cd_x[k] + g_y_0_xxyz_xxxxzz[k];

                g_y_0_xxxyz_xxyyy[k] = -g_y_0_xxyz_xxyyy[k] * cd_x[k] + g_y_0_xxyz_xxxyyy[k];

                g_y_0_xxxyz_xxyyz[k] = -g_y_0_xxyz_xxyyz[k] * cd_x[k] + g_y_0_xxyz_xxxyyz[k];

                g_y_0_xxxyz_xxyzz[k] = -g_y_0_xxyz_xxyzz[k] * cd_x[k] + g_y_0_xxyz_xxxyzz[k];

                g_y_0_xxxyz_xxzzz[k] = -g_y_0_xxyz_xxzzz[k] * cd_x[k] + g_y_0_xxyz_xxxzzz[k];

                g_y_0_xxxyz_xyyyy[k] = -g_y_0_xxyz_xyyyy[k] * cd_x[k] + g_y_0_xxyz_xxyyyy[k];

                g_y_0_xxxyz_xyyyz[k] = -g_y_0_xxyz_xyyyz[k] * cd_x[k] + g_y_0_xxyz_xxyyyz[k];

                g_y_0_xxxyz_xyyzz[k] = -g_y_0_xxyz_xyyzz[k] * cd_x[k] + g_y_0_xxyz_xxyyzz[k];

                g_y_0_xxxyz_xyzzz[k] = -g_y_0_xxyz_xyzzz[k] * cd_x[k] + g_y_0_xxyz_xxyzzz[k];

                g_y_0_xxxyz_xzzzz[k] = -g_y_0_xxyz_xzzzz[k] * cd_x[k] + g_y_0_xxyz_xxzzzz[k];

                g_y_0_xxxyz_yyyyy[k] = -g_y_0_xxyz_yyyyy[k] * cd_x[k] + g_y_0_xxyz_xyyyyy[k];

                g_y_0_xxxyz_yyyyz[k] = -g_y_0_xxyz_yyyyz[k] * cd_x[k] + g_y_0_xxyz_xyyyyz[k];

                g_y_0_xxxyz_yyyzz[k] = -g_y_0_xxyz_yyyzz[k] * cd_x[k] + g_y_0_xxyz_xyyyzz[k];

                g_y_0_xxxyz_yyzzz[k] = -g_y_0_xxyz_yyzzz[k] * cd_x[k] + g_y_0_xxyz_xyyzzz[k];

                g_y_0_xxxyz_yzzzz[k] = -g_y_0_xxyz_yzzzz[k] * cd_x[k] + g_y_0_xxyz_xyzzzz[k];

                g_y_0_xxxyz_zzzzz[k] = -g_y_0_xxyz_zzzzz[k] * cd_x[k] + g_y_0_xxyz_xzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 105);

            auto g_y_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 106);

            auto g_y_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 107);

            auto g_y_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 108);

            auto g_y_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 109);

            auto g_y_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 110);

            auto g_y_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 111);

            auto g_y_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 112);

            auto g_y_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 113);

            auto g_y_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 114);

            auto g_y_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 115);

            auto g_y_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 116);

            auto g_y_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 117);

            auto g_y_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 118);

            auto g_y_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 119);

            auto g_y_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 120);

            auto g_y_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 121);

            auto g_y_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 122);

            auto g_y_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 123);

            auto g_y_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 124);

            auto g_y_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzz_xxxxx, g_y_0_xxxzz_xxxxy, g_y_0_xxxzz_xxxxz, g_y_0_xxxzz_xxxyy, g_y_0_xxxzz_xxxyz, g_y_0_xxxzz_xxxzz, g_y_0_xxxzz_xxyyy, g_y_0_xxxzz_xxyyz, g_y_0_xxxzz_xxyzz, g_y_0_xxxzz_xxzzz, g_y_0_xxxzz_xyyyy, g_y_0_xxxzz_xyyyz, g_y_0_xxxzz_xyyzz, g_y_0_xxxzz_xyzzz, g_y_0_xxxzz_xzzzz, g_y_0_xxxzz_yyyyy, g_y_0_xxxzz_yyyyz, g_y_0_xxxzz_yyyzz, g_y_0_xxxzz_yyzzz, g_y_0_xxxzz_yzzzz, g_y_0_xxxzz_zzzzz, g_y_0_xxzz_xxxxx, g_y_0_xxzz_xxxxxx, g_y_0_xxzz_xxxxxy, g_y_0_xxzz_xxxxxz, g_y_0_xxzz_xxxxy, g_y_0_xxzz_xxxxyy, g_y_0_xxzz_xxxxyz, g_y_0_xxzz_xxxxz, g_y_0_xxzz_xxxxzz, g_y_0_xxzz_xxxyy, g_y_0_xxzz_xxxyyy, g_y_0_xxzz_xxxyyz, g_y_0_xxzz_xxxyz, g_y_0_xxzz_xxxyzz, g_y_0_xxzz_xxxzz, g_y_0_xxzz_xxxzzz, g_y_0_xxzz_xxyyy, g_y_0_xxzz_xxyyyy, g_y_0_xxzz_xxyyyz, g_y_0_xxzz_xxyyz, g_y_0_xxzz_xxyyzz, g_y_0_xxzz_xxyzz, g_y_0_xxzz_xxyzzz, g_y_0_xxzz_xxzzz, g_y_0_xxzz_xxzzzz, g_y_0_xxzz_xyyyy, g_y_0_xxzz_xyyyyy, g_y_0_xxzz_xyyyyz, g_y_0_xxzz_xyyyz, g_y_0_xxzz_xyyyzz, g_y_0_xxzz_xyyzz, g_y_0_xxzz_xyyzzz, g_y_0_xxzz_xyzzz, g_y_0_xxzz_xyzzzz, g_y_0_xxzz_xzzzz, g_y_0_xxzz_xzzzzz, g_y_0_xxzz_yyyyy, g_y_0_xxzz_yyyyz, g_y_0_xxzz_yyyzz, g_y_0_xxzz_yyzzz, g_y_0_xxzz_yzzzz, g_y_0_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xxxxx[k] = -g_y_0_xxzz_xxxxx[k] * cd_x[k] + g_y_0_xxzz_xxxxxx[k];

                g_y_0_xxxzz_xxxxy[k] = -g_y_0_xxzz_xxxxy[k] * cd_x[k] + g_y_0_xxzz_xxxxxy[k];

                g_y_0_xxxzz_xxxxz[k] = -g_y_0_xxzz_xxxxz[k] * cd_x[k] + g_y_0_xxzz_xxxxxz[k];

                g_y_0_xxxzz_xxxyy[k] = -g_y_0_xxzz_xxxyy[k] * cd_x[k] + g_y_0_xxzz_xxxxyy[k];

                g_y_0_xxxzz_xxxyz[k] = -g_y_0_xxzz_xxxyz[k] * cd_x[k] + g_y_0_xxzz_xxxxyz[k];

                g_y_0_xxxzz_xxxzz[k] = -g_y_0_xxzz_xxxzz[k] * cd_x[k] + g_y_0_xxzz_xxxxzz[k];

                g_y_0_xxxzz_xxyyy[k] = -g_y_0_xxzz_xxyyy[k] * cd_x[k] + g_y_0_xxzz_xxxyyy[k];

                g_y_0_xxxzz_xxyyz[k] = -g_y_0_xxzz_xxyyz[k] * cd_x[k] + g_y_0_xxzz_xxxyyz[k];

                g_y_0_xxxzz_xxyzz[k] = -g_y_0_xxzz_xxyzz[k] * cd_x[k] + g_y_0_xxzz_xxxyzz[k];

                g_y_0_xxxzz_xxzzz[k] = -g_y_0_xxzz_xxzzz[k] * cd_x[k] + g_y_0_xxzz_xxxzzz[k];

                g_y_0_xxxzz_xyyyy[k] = -g_y_0_xxzz_xyyyy[k] * cd_x[k] + g_y_0_xxzz_xxyyyy[k];

                g_y_0_xxxzz_xyyyz[k] = -g_y_0_xxzz_xyyyz[k] * cd_x[k] + g_y_0_xxzz_xxyyyz[k];

                g_y_0_xxxzz_xyyzz[k] = -g_y_0_xxzz_xyyzz[k] * cd_x[k] + g_y_0_xxzz_xxyyzz[k];

                g_y_0_xxxzz_xyzzz[k] = -g_y_0_xxzz_xyzzz[k] * cd_x[k] + g_y_0_xxzz_xxyzzz[k];

                g_y_0_xxxzz_xzzzz[k] = -g_y_0_xxzz_xzzzz[k] * cd_x[k] + g_y_0_xxzz_xxzzzz[k];

                g_y_0_xxxzz_yyyyy[k] = -g_y_0_xxzz_yyyyy[k] * cd_x[k] + g_y_0_xxzz_xyyyyy[k];

                g_y_0_xxxzz_yyyyz[k] = -g_y_0_xxzz_yyyyz[k] * cd_x[k] + g_y_0_xxzz_xyyyyz[k];

                g_y_0_xxxzz_yyyzz[k] = -g_y_0_xxzz_yyyzz[k] * cd_x[k] + g_y_0_xxzz_xyyyzz[k];

                g_y_0_xxxzz_yyzzz[k] = -g_y_0_xxzz_yyzzz[k] * cd_x[k] + g_y_0_xxzz_xyyzzz[k];

                g_y_0_xxxzz_yzzzz[k] = -g_y_0_xxzz_yzzzz[k] * cd_x[k] + g_y_0_xxzz_xyzzzz[k];

                g_y_0_xxxzz_zzzzz[k] = -g_y_0_xxzz_zzzzz[k] * cd_x[k] + g_y_0_xxzz_xzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 126);

            auto g_y_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 127);

            auto g_y_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 128);

            auto g_y_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 129);

            auto g_y_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 130);

            auto g_y_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 131);

            auto g_y_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 132);

            auto g_y_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 133);

            auto g_y_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 134);

            auto g_y_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 135);

            auto g_y_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 136);

            auto g_y_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 137);

            auto g_y_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 138);

            auto g_y_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 139);

            auto g_y_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 140);

            auto g_y_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 141);

            auto g_y_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 142);

            auto g_y_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 143);

            auto g_y_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 144);

            auto g_y_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 145);

            auto g_y_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyy_xxxxx, g_y_0_xxyyy_xxxxy, g_y_0_xxyyy_xxxxz, g_y_0_xxyyy_xxxyy, g_y_0_xxyyy_xxxyz, g_y_0_xxyyy_xxxzz, g_y_0_xxyyy_xxyyy, g_y_0_xxyyy_xxyyz, g_y_0_xxyyy_xxyzz, g_y_0_xxyyy_xxzzz, g_y_0_xxyyy_xyyyy, g_y_0_xxyyy_xyyyz, g_y_0_xxyyy_xyyzz, g_y_0_xxyyy_xyzzz, g_y_0_xxyyy_xzzzz, g_y_0_xxyyy_yyyyy, g_y_0_xxyyy_yyyyz, g_y_0_xxyyy_yyyzz, g_y_0_xxyyy_yyzzz, g_y_0_xxyyy_yzzzz, g_y_0_xxyyy_zzzzz, g_y_0_xyyy_xxxxx, g_y_0_xyyy_xxxxxx, g_y_0_xyyy_xxxxxy, g_y_0_xyyy_xxxxxz, g_y_0_xyyy_xxxxy, g_y_0_xyyy_xxxxyy, g_y_0_xyyy_xxxxyz, g_y_0_xyyy_xxxxz, g_y_0_xyyy_xxxxzz, g_y_0_xyyy_xxxyy, g_y_0_xyyy_xxxyyy, g_y_0_xyyy_xxxyyz, g_y_0_xyyy_xxxyz, g_y_0_xyyy_xxxyzz, g_y_0_xyyy_xxxzz, g_y_0_xyyy_xxxzzz, g_y_0_xyyy_xxyyy, g_y_0_xyyy_xxyyyy, g_y_0_xyyy_xxyyyz, g_y_0_xyyy_xxyyz, g_y_0_xyyy_xxyyzz, g_y_0_xyyy_xxyzz, g_y_0_xyyy_xxyzzz, g_y_0_xyyy_xxzzz, g_y_0_xyyy_xxzzzz, g_y_0_xyyy_xyyyy, g_y_0_xyyy_xyyyyy, g_y_0_xyyy_xyyyyz, g_y_0_xyyy_xyyyz, g_y_0_xyyy_xyyyzz, g_y_0_xyyy_xyyzz, g_y_0_xyyy_xyyzzz, g_y_0_xyyy_xyzzz, g_y_0_xyyy_xyzzzz, g_y_0_xyyy_xzzzz, g_y_0_xyyy_xzzzzz, g_y_0_xyyy_yyyyy, g_y_0_xyyy_yyyyz, g_y_0_xyyy_yyyzz, g_y_0_xyyy_yyzzz, g_y_0_xyyy_yzzzz, g_y_0_xyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xxxxx[k] = -g_y_0_xyyy_xxxxx[k] * cd_x[k] + g_y_0_xyyy_xxxxxx[k];

                g_y_0_xxyyy_xxxxy[k] = -g_y_0_xyyy_xxxxy[k] * cd_x[k] + g_y_0_xyyy_xxxxxy[k];

                g_y_0_xxyyy_xxxxz[k] = -g_y_0_xyyy_xxxxz[k] * cd_x[k] + g_y_0_xyyy_xxxxxz[k];

                g_y_0_xxyyy_xxxyy[k] = -g_y_0_xyyy_xxxyy[k] * cd_x[k] + g_y_0_xyyy_xxxxyy[k];

                g_y_0_xxyyy_xxxyz[k] = -g_y_0_xyyy_xxxyz[k] * cd_x[k] + g_y_0_xyyy_xxxxyz[k];

                g_y_0_xxyyy_xxxzz[k] = -g_y_0_xyyy_xxxzz[k] * cd_x[k] + g_y_0_xyyy_xxxxzz[k];

                g_y_0_xxyyy_xxyyy[k] = -g_y_0_xyyy_xxyyy[k] * cd_x[k] + g_y_0_xyyy_xxxyyy[k];

                g_y_0_xxyyy_xxyyz[k] = -g_y_0_xyyy_xxyyz[k] * cd_x[k] + g_y_0_xyyy_xxxyyz[k];

                g_y_0_xxyyy_xxyzz[k] = -g_y_0_xyyy_xxyzz[k] * cd_x[k] + g_y_0_xyyy_xxxyzz[k];

                g_y_0_xxyyy_xxzzz[k] = -g_y_0_xyyy_xxzzz[k] * cd_x[k] + g_y_0_xyyy_xxxzzz[k];

                g_y_0_xxyyy_xyyyy[k] = -g_y_0_xyyy_xyyyy[k] * cd_x[k] + g_y_0_xyyy_xxyyyy[k];

                g_y_0_xxyyy_xyyyz[k] = -g_y_0_xyyy_xyyyz[k] * cd_x[k] + g_y_0_xyyy_xxyyyz[k];

                g_y_0_xxyyy_xyyzz[k] = -g_y_0_xyyy_xyyzz[k] * cd_x[k] + g_y_0_xyyy_xxyyzz[k];

                g_y_0_xxyyy_xyzzz[k] = -g_y_0_xyyy_xyzzz[k] * cd_x[k] + g_y_0_xyyy_xxyzzz[k];

                g_y_0_xxyyy_xzzzz[k] = -g_y_0_xyyy_xzzzz[k] * cd_x[k] + g_y_0_xyyy_xxzzzz[k];

                g_y_0_xxyyy_yyyyy[k] = -g_y_0_xyyy_yyyyy[k] * cd_x[k] + g_y_0_xyyy_xyyyyy[k];

                g_y_0_xxyyy_yyyyz[k] = -g_y_0_xyyy_yyyyz[k] * cd_x[k] + g_y_0_xyyy_xyyyyz[k];

                g_y_0_xxyyy_yyyzz[k] = -g_y_0_xyyy_yyyzz[k] * cd_x[k] + g_y_0_xyyy_xyyyzz[k];

                g_y_0_xxyyy_yyzzz[k] = -g_y_0_xyyy_yyzzz[k] * cd_x[k] + g_y_0_xyyy_xyyzzz[k];

                g_y_0_xxyyy_yzzzz[k] = -g_y_0_xyyy_yzzzz[k] * cd_x[k] + g_y_0_xyyy_xyzzzz[k];

                g_y_0_xxyyy_zzzzz[k] = -g_y_0_xyyy_zzzzz[k] * cd_x[k] + g_y_0_xyyy_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 147);

            auto g_y_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 148);

            auto g_y_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 149);

            auto g_y_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 150);

            auto g_y_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 151);

            auto g_y_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 152);

            auto g_y_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 153);

            auto g_y_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 154);

            auto g_y_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 155);

            auto g_y_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 156);

            auto g_y_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 157);

            auto g_y_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 158);

            auto g_y_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 159);

            auto g_y_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 160);

            auto g_y_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 161);

            auto g_y_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 162);

            auto g_y_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 163);

            auto g_y_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 164);

            auto g_y_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 165);

            auto g_y_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 166);

            auto g_y_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyz_xxxxx, g_y_0_xxyyz_xxxxy, g_y_0_xxyyz_xxxxz, g_y_0_xxyyz_xxxyy, g_y_0_xxyyz_xxxyz, g_y_0_xxyyz_xxxzz, g_y_0_xxyyz_xxyyy, g_y_0_xxyyz_xxyyz, g_y_0_xxyyz_xxyzz, g_y_0_xxyyz_xxzzz, g_y_0_xxyyz_xyyyy, g_y_0_xxyyz_xyyyz, g_y_0_xxyyz_xyyzz, g_y_0_xxyyz_xyzzz, g_y_0_xxyyz_xzzzz, g_y_0_xxyyz_yyyyy, g_y_0_xxyyz_yyyyz, g_y_0_xxyyz_yyyzz, g_y_0_xxyyz_yyzzz, g_y_0_xxyyz_yzzzz, g_y_0_xxyyz_zzzzz, g_y_0_xyyz_xxxxx, g_y_0_xyyz_xxxxxx, g_y_0_xyyz_xxxxxy, g_y_0_xyyz_xxxxxz, g_y_0_xyyz_xxxxy, g_y_0_xyyz_xxxxyy, g_y_0_xyyz_xxxxyz, g_y_0_xyyz_xxxxz, g_y_0_xyyz_xxxxzz, g_y_0_xyyz_xxxyy, g_y_0_xyyz_xxxyyy, g_y_0_xyyz_xxxyyz, g_y_0_xyyz_xxxyz, g_y_0_xyyz_xxxyzz, g_y_0_xyyz_xxxzz, g_y_0_xyyz_xxxzzz, g_y_0_xyyz_xxyyy, g_y_0_xyyz_xxyyyy, g_y_0_xyyz_xxyyyz, g_y_0_xyyz_xxyyz, g_y_0_xyyz_xxyyzz, g_y_0_xyyz_xxyzz, g_y_0_xyyz_xxyzzz, g_y_0_xyyz_xxzzz, g_y_0_xyyz_xxzzzz, g_y_0_xyyz_xyyyy, g_y_0_xyyz_xyyyyy, g_y_0_xyyz_xyyyyz, g_y_0_xyyz_xyyyz, g_y_0_xyyz_xyyyzz, g_y_0_xyyz_xyyzz, g_y_0_xyyz_xyyzzz, g_y_0_xyyz_xyzzz, g_y_0_xyyz_xyzzzz, g_y_0_xyyz_xzzzz, g_y_0_xyyz_xzzzzz, g_y_0_xyyz_yyyyy, g_y_0_xyyz_yyyyz, g_y_0_xyyz_yyyzz, g_y_0_xyyz_yyzzz, g_y_0_xyyz_yzzzz, g_y_0_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xxxxx[k] = -g_y_0_xyyz_xxxxx[k] * cd_x[k] + g_y_0_xyyz_xxxxxx[k];

                g_y_0_xxyyz_xxxxy[k] = -g_y_0_xyyz_xxxxy[k] * cd_x[k] + g_y_0_xyyz_xxxxxy[k];

                g_y_0_xxyyz_xxxxz[k] = -g_y_0_xyyz_xxxxz[k] * cd_x[k] + g_y_0_xyyz_xxxxxz[k];

                g_y_0_xxyyz_xxxyy[k] = -g_y_0_xyyz_xxxyy[k] * cd_x[k] + g_y_0_xyyz_xxxxyy[k];

                g_y_0_xxyyz_xxxyz[k] = -g_y_0_xyyz_xxxyz[k] * cd_x[k] + g_y_0_xyyz_xxxxyz[k];

                g_y_0_xxyyz_xxxzz[k] = -g_y_0_xyyz_xxxzz[k] * cd_x[k] + g_y_0_xyyz_xxxxzz[k];

                g_y_0_xxyyz_xxyyy[k] = -g_y_0_xyyz_xxyyy[k] * cd_x[k] + g_y_0_xyyz_xxxyyy[k];

                g_y_0_xxyyz_xxyyz[k] = -g_y_0_xyyz_xxyyz[k] * cd_x[k] + g_y_0_xyyz_xxxyyz[k];

                g_y_0_xxyyz_xxyzz[k] = -g_y_0_xyyz_xxyzz[k] * cd_x[k] + g_y_0_xyyz_xxxyzz[k];

                g_y_0_xxyyz_xxzzz[k] = -g_y_0_xyyz_xxzzz[k] * cd_x[k] + g_y_0_xyyz_xxxzzz[k];

                g_y_0_xxyyz_xyyyy[k] = -g_y_0_xyyz_xyyyy[k] * cd_x[k] + g_y_0_xyyz_xxyyyy[k];

                g_y_0_xxyyz_xyyyz[k] = -g_y_0_xyyz_xyyyz[k] * cd_x[k] + g_y_0_xyyz_xxyyyz[k];

                g_y_0_xxyyz_xyyzz[k] = -g_y_0_xyyz_xyyzz[k] * cd_x[k] + g_y_0_xyyz_xxyyzz[k];

                g_y_0_xxyyz_xyzzz[k] = -g_y_0_xyyz_xyzzz[k] * cd_x[k] + g_y_0_xyyz_xxyzzz[k];

                g_y_0_xxyyz_xzzzz[k] = -g_y_0_xyyz_xzzzz[k] * cd_x[k] + g_y_0_xyyz_xxzzzz[k];

                g_y_0_xxyyz_yyyyy[k] = -g_y_0_xyyz_yyyyy[k] * cd_x[k] + g_y_0_xyyz_xyyyyy[k];

                g_y_0_xxyyz_yyyyz[k] = -g_y_0_xyyz_yyyyz[k] * cd_x[k] + g_y_0_xyyz_xyyyyz[k];

                g_y_0_xxyyz_yyyzz[k] = -g_y_0_xyyz_yyyzz[k] * cd_x[k] + g_y_0_xyyz_xyyyzz[k];

                g_y_0_xxyyz_yyzzz[k] = -g_y_0_xyyz_yyzzz[k] * cd_x[k] + g_y_0_xyyz_xyyzzz[k];

                g_y_0_xxyyz_yzzzz[k] = -g_y_0_xyyz_yzzzz[k] * cd_x[k] + g_y_0_xyyz_xyzzzz[k];

                g_y_0_xxyyz_zzzzz[k] = -g_y_0_xyyz_zzzzz[k] * cd_x[k] + g_y_0_xyyz_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 168);

            auto g_y_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 169);

            auto g_y_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 170);

            auto g_y_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 171);

            auto g_y_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 172);

            auto g_y_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 173);

            auto g_y_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 174);

            auto g_y_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 175);

            auto g_y_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 176);

            auto g_y_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 177);

            auto g_y_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 178);

            auto g_y_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 179);

            auto g_y_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 180);

            auto g_y_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 181);

            auto g_y_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 182);

            auto g_y_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 183);

            auto g_y_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 184);

            auto g_y_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 185);

            auto g_y_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 186);

            auto g_y_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 187);

            auto g_y_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzz_xxxxx, g_y_0_xxyzz_xxxxy, g_y_0_xxyzz_xxxxz, g_y_0_xxyzz_xxxyy, g_y_0_xxyzz_xxxyz, g_y_0_xxyzz_xxxzz, g_y_0_xxyzz_xxyyy, g_y_0_xxyzz_xxyyz, g_y_0_xxyzz_xxyzz, g_y_0_xxyzz_xxzzz, g_y_0_xxyzz_xyyyy, g_y_0_xxyzz_xyyyz, g_y_0_xxyzz_xyyzz, g_y_0_xxyzz_xyzzz, g_y_0_xxyzz_xzzzz, g_y_0_xxyzz_yyyyy, g_y_0_xxyzz_yyyyz, g_y_0_xxyzz_yyyzz, g_y_0_xxyzz_yyzzz, g_y_0_xxyzz_yzzzz, g_y_0_xxyzz_zzzzz, g_y_0_xyzz_xxxxx, g_y_0_xyzz_xxxxxx, g_y_0_xyzz_xxxxxy, g_y_0_xyzz_xxxxxz, g_y_0_xyzz_xxxxy, g_y_0_xyzz_xxxxyy, g_y_0_xyzz_xxxxyz, g_y_0_xyzz_xxxxz, g_y_0_xyzz_xxxxzz, g_y_0_xyzz_xxxyy, g_y_0_xyzz_xxxyyy, g_y_0_xyzz_xxxyyz, g_y_0_xyzz_xxxyz, g_y_0_xyzz_xxxyzz, g_y_0_xyzz_xxxzz, g_y_0_xyzz_xxxzzz, g_y_0_xyzz_xxyyy, g_y_0_xyzz_xxyyyy, g_y_0_xyzz_xxyyyz, g_y_0_xyzz_xxyyz, g_y_0_xyzz_xxyyzz, g_y_0_xyzz_xxyzz, g_y_0_xyzz_xxyzzz, g_y_0_xyzz_xxzzz, g_y_0_xyzz_xxzzzz, g_y_0_xyzz_xyyyy, g_y_0_xyzz_xyyyyy, g_y_0_xyzz_xyyyyz, g_y_0_xyzz_xyyyz, g_y_0_xyzz_xyyyzz, g_y_0_xyzz_xyyzz, g_y_0_xyzz_xyyzzz, g_y_0_xyzz_xyzzz, g_y_0_xyzz_xyzzzz, g_y_0_xyzz_xzzzz, g_y_0_xyzz_xzzzzz, g_y_0_xyzz_yyyyy, g_y_0_xyzz_yyyyz, g_y_0_xyzz_yyyzz, g_y_0_xyzz_yyzzz, g_y_0_xyzz_yzzzz, g_y_0_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xxxxx[k] = -g_y_0_xyzz_xxxxx[k] * cd_x[k] + g_y_0_xyzz_xxxxxx[k];

                g_y_0_xxyzz_xxxxy[k] = -g_y_0_xyzz_xxxxy[k] * cd_x[k] + g_y_0_xyzz_xxxxxy[k];

                g_y_0_xxyzz_xxxxz[k] = -g_y_0_xyzz_xxxxz[k] * cd_x[k] + g_y_0_xyzz_xxxxxz[k];

                g_y_0_xxyzz_xxxyy[k] = -g_y_0_xyzz_xxxyy[k] * cd_x[k] + g_y_0_xyzz_xxxxyy[k];

                g_y_0_xxyzz_xxxyz[k] = -g_y_0_xyzz_xxxyz[k] * cd_x[k] + g_y_0_xyzz_xxxxyz[k];

                g_y_0_xxyzz_xxxzz[k] = -g_y_0_xyzz_xxxzz[k] * cd_x[k] + g_y_0_xyzz_xxxxzz[k];

                g_y_0_xxyzz_xxyyy[k] = -g_y_0_xyzz_xxyyy[k] * cd_x[k] + g_y_0_xyzz_xxxyyy[k];

                g_y_0_xxyzz_xxyyz[k] = -g_y_0_xyzz_xxyyz[k] * cd_x[k] + g_y_0_xyzz_xxxyyz[k];

                g_y_0_xxyzz_xxyzz[k] = -g_y_0_xyzz_xxyzz[k] * cd_x[k] + g_y_0_xyzz_xxxyzz[k];

                g_y_0_xxyzz_xxzzz[k] = -g_y_0_xyzz_xxzzz[k] * cd_x[k] + g_y_0_xyzz_xxxzzz[k];

                g_y_0_xxyzz_xyyyy[k] = -g_y_0_xyzz_xyyyy[k] * cd_x[k] + g_y_0_xyzz_xxyyyy[k];

                g_y_0_xxyzz_xyyyz[k] = -g_y_0_xyzz_xyyyz[k] * cd_x[k] + g_y_0_xyzz_xxyyyz[k];

                g_y_0_xxyzz_xyyzz[k] = -g_y_0_xyzz_xyyzz[k] * cd_x[k] + g_y_0_xyzz_xxyyzz[k];

                g_y_0_xxyzz_xyzzz[k] = -g_y_0_xyzz_xyzzz[k] * cd_x[k] + g_y_0_xyzz_xxyzzz[k];

                g_y_0_xxyzz_xzzzz[k] = -g_y_0_xyzz_xzzzz[k] * cd_x[k] + g_y_0_xyzz_xxzzzz[k];

                g_y_0_xxyzz_yyyyy[k] = -g_y_0_xyzz_yyyyy[k] * cd_x[k] + g_y_0_xyzz_xyyyyy[k];

                g_y_0_xxyzz_yyyyz[k] = -g_y_0_xyzz_yyyyz[k] * cd_x[k] + g_y_0_xyzz_xyyyyz[k];

                g_y_0_xxyzz_yyyzz[k] = -g_y_0_xyzz_yyyzz[k] * cd_x[k] + g_y_0_xyzz_xyyyzz[k];

                g_y_0_xxyzz_yyzzz[k] = -g_y_0_xyzz_yyzzz[k] * cd_x[k] + g_y_0_xyzz_xyyzzz[k];

                g_y_0_xxyzz_yzzzz[k] = -g_y_0_xyzz_yzzzz[k] * cd_x[k] + g_y_0_xyzz_xyzzzz[k];

                g_y_0_xxyzz_zzzzz[k] = -g_y_0_xyzz_zzzzz[k] * cd_x[k] + g_y_0_xyzz_xzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 189);

            auto g_y_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 190);

            auto g_y_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 191);

            auto g_y_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 192);

            auto g_y_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 193);

            auto g_y_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 194);

            auto g_y_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 195);

            auto g_y_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 196);

            auto g_y_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 197);

            auto g_y_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 198);

            auto g_y_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 199);

            auto g_y_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 200);

            auto g_y_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 201);

            auto g_y_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 202);

            auto g_y_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 203);

            auto g_y_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 204);

            auto g_y_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 205);

            auto g_y_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 206);

            auto g_y_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 207);

            auto g_y_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 208);

            auto g_y_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzz_xxxxx, g_y_0_xxzzz_xxxxy, g_y_0_xxzzz_xxxxz, g_y_0_xxzzz_xxxyy, g_y_0_xxzzz_xxxyz, g_y_0_xxzzz_xxxzz, g_y_0_xxzzz_xxyyy, g_y_0_xxzzz_xxyyz, g_y_0_xxzzz_xxyzz, g_y_0_xxzzz_xxzzz, g_y_0_xxzzz_xyyyy, g_y_0_xxzzz_xyyyz, g_y_0_xxzzz_xyyzz, g_y_0_xxzzz_xyzzz, g_y_0_xxzzz_xzzzz, g_y_0_xxzzz_yyyyy, g_y_0_xxzzz_yyyyz, g_y_0_xxzzz_yyyzz, g_y_0_xxzzz_yyzzz, g_y_0_xxzzz_yzzzz, g_y_0_xxzzz_zzzzz, g_y_0_xzzz_xxxxx, g_y_0_xzzz_xxxxxx, g_y_0_xzzz_xxxxxy, g_y_0_xzzz_xxxxxz, g_y_0_xzzz_xxxxy, g_y_0_xzzz_xxxxyy, g_y_0_xzzz_xxxxyz, g_y_0_xzzz_xxxxz, g_y_0_xzzz_xxxxzz, g_y_0_xzzz_xxxyy, g_y_0_xzzz_xxxyyy, g_y_0_xzzz_xxxyyz, g_y_0_xzzz_xxxyz, g_y_0_xzzz_xxxyzz, g_y_0_xzzz_xxxzz, g_y_0_xzzz_xxxzzz, g_y_0_xzzz_xxyyy, g_y_0_xzzz_xxyyyy, g_y_0_xzzz_xxyyyz, g_y_0_xzzz_xxyyz, g_y_0_xzzz_xxyyzz, g_y_0_xzzz_xxyzz, g_y_0_xzzz_xxyzzz, g_y_0_xzzz_xxzzz, g_y_0_xzzz_xxzzzz, g_y_0_xzzz_xyyyy, g_y_0_xzzz_xyyyyy, g_y_0_xzzz_xyyyyz, g_y_0_xzzz_xyyyz, g_y_0_xzzz_xyyyzz, g_y_0_xzzz_xyyzz, g_y_0_xzzz_xyyzzz, g_y_0_xzzz_xyzzz, g_y_0_xzzz_xyzzzz, g_y_0_xzzz_xzzzz, g_y_0_xzzz_xzzzzz, g_y_0_xzzz_yyyyy, g_y_0_xzzz_yyyyz, g_y_0_xzzz_yyyzz, g_y_0_xzzz_yyzzz, g_y_0_xzzz_yzzzz, g_y_0_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xxxxx[k] = -g_y_0_xzzz_xxxxx[k] * cd_x[k] + g_y_0_xzzz_xxxxxx[k];

                g_y_0_xxzzz_xxxxy[k] = -g_y_0_xzzz_xxxxy[k] * cd_x[k] + g_y_0_xzzz_xxxxxy[k];

                g_y_0_xxzzz_xxxxz[k] = -g_y_0_xzzz_xxxxz[k] * cd_x[k] + g_y_0_xzzz_xxxxxz[k];

                g_y_0_xxzzz_xxxyy[k] = -g_y_0_xzzz_xxxyy[k] * cd_x[k] + g_y_0_xzzz_xxxxyy[k];

                g_y_0_xxzzz_xxxyz[k] = -g_y_0_xzzz_xxxyz[k] * cd_x[k] + g_y_0_xzzz_xxxxyz[k];

                g_y_0_xxzzz_xxxzz[k] = -g_y_0_xzzz_xxxzz[k] * cd_x[k] + g_y_0_xzzz_xxxxzz[k];

                g_y_0_xxzzz_xxyyy[k] = -g_y_0_xzzz_xxyyy[k] * cd_x[k] + g_y_0_xzzz_xxxyyy[k];

                g_y_0_xxzzz_xxyyz[k] = -g_y_0_xzzz_xxyyz[k] * cd_x[k] + g_y_0_xzzz_xxxyyz[k];

                g_y_0_xxzzz_xxyzz[k] = -g_y_0_xzzz_xxyzz[k] * cd_x[k] + g_y_0_xzzz_xxxyzz[k];

                g_y_0_xxzzz_xxzzz[k] = -g_y_0_xzzz_xxzzz[k] * cd_x[k] + g_y_0_xzzz_xxxzzz[k];

                g_y_0_xxzzz_xyyyy[k] = -g_y_0_xzzz_xyyyy[k] * cd_x[k] + g_y_0_xzzz_xxyyyy[k];

                g_y_0_xxzzz_xyyyz[k] = -g_y_0_xzzz_xyyyz[k] * cd_x[k] + g_y_0_xzzz_xxyyyz[k];

                g_y_0_xxzzz_xyyzz[k] = -g_y_0_xzzz_xyyzz[k] * cd_x[k] + g_y_0_xzzz_xxyyzz[k];

                g_y_0_xxzzz_xyzzz[k] = -g_y_0_xzzz_xyzzz[k] * cd_x[k] + g_y_0_xzzz_xxyzzz[k];

                g_y_0_xxzzz_xzzzz[k] = -g_y_0_xzzz_xzzzz[k] * cd_x[k] + g_y_0_xzzz_xxzzzz[k];

                g_y_0_xxzzz_yyyyy[k] = -g_y_0_xzzz_yyyyy[k] * cd_x[k] + g_y_0_xzzz_xyyyyy[k];

                g_y_0_xxzzz_yyyyz[k] = -g_y_0_xzzz_yyyyz[k] * cd_x[k] + g_y_0_xzzz_xyyyyz[k];

                g_y_0_xxzzz_yyyzz[k] = -g_y_0_xzzz_yyyzz[k] * cd_x[k] + g_y_0_xzzz_xyyyzz[k];

                g_y_0_xxzzz_yyzzz[k] = -g_y_0_xzzz_yyzzz[k] * cd_x[k] + g_y_0_xzzz_xyyzzz[k];

                g_y_0_xxzzz_yzzzz[k] = -g_y_0_xzzz_yzzzz[k] * cd_x[k] + g_y_0_xzzz_xyzzzz[k];

                g_y_0_xxzzz_zzzzz[k] = -g_y_0_xzzz_zzzzz[k] * cd_x[k] + g_y_0_xzzz_xzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 210);

            auto g_y_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 211);

            auto g_y_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 212);

            auto g_y_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 213);

            auto g_y_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 214);

            auto g_y_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 215);

            auto g_y_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 216);

            auto g_y_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 217);

            auto g_y_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 218);

            auto g_y_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 219);

            auto g_y_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 220);

            auto g_y_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 221);

            auto g_y_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 222);

            auto g_y_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 223);

            auto g_y_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 224);

            auto g_y_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 225);

            auto g_y_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 226);

            auto g_y_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 227);

            auto g_y_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 228);

            auto g_y_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 229);

            auto g_y_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyy_xxxxx, g_y_0_xyyyy_xxxxy, g_y_0_xyyyy_xxxxz, g_y_0_xyyyy_xxxyy, g_y_0_xyyyy_xxxyz, g_y_0_xyyyy_xxxzz, g_y_0_xyyyy_xxyyy, g_y_0_xyyyy_xxyyz, g_y_0_xyyyy_xxyzz, g_y_0_xyyyy_xxzzz, g_y_0_xyyyy_xyyyy, g_y_0_xyyyy_xyyyz, g_y_0_xyyyy_xyyzz, g_y_0_xyyyy_xyzzz, g_y_0_xyyyy_xzzzz, g_y_0_xyyyy_yyyyy, g_y_0_xyyyy_yyyyz, g_y_0_xyyyy_yyyzz, g_y_0_xyyyy_yyzzz, g_y_0_xyyyy_yzzzz, g_y_0_xyyyy_zzzzz, g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xxxxx[k] = -g_y_0_yyyy_xxxxx[k] * cd_x[k] + g_y_0_yyyy_xxxxxx[k];

                g_y_0_xyyyy_xxxxy[k] = -g_y_0_yyyy_xxxxy[k] * cd_x[k] + g_y_0_yyyy_xxxxxy[k];

                g_y_0_xyyyy_xxxxz[k] = -g_y_0_yyyy_xxxxz[k] * cd_x[k] + g_y_0_yyyy_xxxxxz[k];

                g_y_0_xyyyy_xxxyy[k] = -g_y_0_yyyy_xxxyy[k] * cd_x[k] + g_y_0_yyyy_xxxxyy[k];

                g_y_0_xyyyy_xxxyz[k] = -g_y_0_yyyy_xxxyz[k] * cd_x[k] + g_y_0_yyyy_xxxxyz[k];

                g_y_0_xyyyy_xxxzz[k] = -g_y_0_yyyy_xxxzz[k] * cd_x[k] + g_y_0_yyyy_xxxxzz[k];

                g_y_0_xyyyy_xxyyy[k] = -g_y_0_yyyy_xxyyy[k] * cd_x[k] + g_y_0_yyyy_xxxyyy[k];

                g_y_0_xyyyy_xxyyz[k] = -g_y_0_yyyy_xxyyz[k] * cd_x[k] + g_y_0_yyyy_xxxyyz[k];

                g_y_0_xyyyy_xxyzz[k] = -g_y_0_yyyy_xxyzz[k] * cd_x[k] + g_y_0_yyyy_xxxyzz[k];

                g_y_0_xyyyy_xxzzz[k] = -g_y_0_yyyy_xxzzz[k] * cd_x[k] + g_y_0_yyyy_xxxzzz[k];

                g_y_0_xyyyy_xyyyy[k] = -g_y_0_yyyy_xyyyy[k] * cd_x[k] + g_y_0_yyyy_xxyyyy[k];

                g_y_0_xyyyy_xyyyz[k] = -g_y_0_yyyy_xyyyz[k] * cd_x[k] + g_y_0_yyyy_xxyyyz[k];

                g_y_0_xyyyy_xyyzz[k] = -g_y_0_yyyy_xyyzz[k] * cd_x[k] + g_y_0_yyyy_xxyyzz[k];

                g_y_0_xyyyy_xyzzz[k] = -g_y_0_yyyy_xyzzz[k] * cd_x[k] + g_y_0_yyyy_xxyzzz[k];

                g_y_0_xyyyy_xzzzz[k] = -g_y_0_yyyy_xzzzz[k] * cd_x[k] + g_y_0_yyyy_xxzzzz[k];

                g_y_0_xyyyy_yyyyy[k] = -g_y_0_yyyy_yyyyy[k] * cd_x[k] + g_y_0_yyyy_xyyyyy[k];

                g_y_0_xyyyy_yyyyz[k] = -g_y_0_yyyy_yyyyz[k] * cd_x[k] + g_y_0_yyyy_xyyyyz[k];

                g_y_0_xyyyy_yyyzz[k] = -g_y_0_yyyy_yyyzz[k] * cd_x[k] + g_y_0_yyyy_xyyyzz[k];

                g_y_0_xyyyy_yyzzz[k] = -g_y_0_yyyy_yyzzz[k] * cd_x[k] + g_y_0_yyyy_xyyzzz[k];

                g_y_0_xyyyy_yzzzz[k] = -g_y_0_yyyy_yzzzz[k] * cd_x[k] + g_y_0_yyyy_xyzzzz[k];

                g_y_0_xyyyy_zzzzz[k] = -g_y_0_yyyy_zzzzz[k] * cd_x[k] + g_y_0_yyyy_xzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 231);

            auto g_y_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 232);

            auto g_y_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 233);

            auto g_y_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 234);

            auto g_y_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 235);

            auto g_y_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 236);

            auto g_y_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 237);

            auto g_y_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 238);

            auto g_y_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 239);

            auto g_y_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 240);

            auto g_y_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 241);

            auto g_y_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 242);

            auto g_y_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 243);

            auto g_y_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 244);

            auto g_y_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 245);

            auto g_y_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 246);

            auto g_y_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 247);

            auto g_y_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 248);

            auto g_y_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 249);

            auto g_y_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 250);

            auto g_y_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyz_xxxxx, g_y_0_xyyyz_xxxxy, g_y_0_xyyyz_xxxxz, g_y_0_xyyyz_xxxyy, g_y_0_xyyyz_xxxyz, g_y_0_xyyyz_xxxzz, g_y_0_xyyyz_xxyyy, g_y_0_xyyyz_xxyyz, g_y_0_xyyyz_xxyzz, g_y_0_xyyyz_xxzzz, g_y_0_xyyyz_xyyyy, g_y_0_xyyyz_xyyyz, g_y_0_xyyyz_xyyzz, g_y_0_xyyyz_xyzzz, g_y_0_xyyyz_xzzzz, g_y_0_xyyyz_yyyyy, g_y_0_xyyyz_yyyyz, g_y_0_xyyyz_yyyzz, g_y_0_xyyyz_yyzzz, g_y_0_xyyyz_yzzzz, g_y_0_xyyyz_zzzzz, g_y_0_yyyz_xxxxx, g_y_0_yyyz_xxxxxx, g_y_0_yyyz_xxxxxy, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxy, g_y_0_yyyz_xxxxyy, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxyy, g_y_0_yyyz_xxxyyy, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxyyy, g_y_0_yyyz_xxyyyy, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xyyyy, g_y_0_yyyz_xyyyyy, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_yyyyy, g_y_0_yyyz_yyyyz, g_y_0_yyyz_yyyzz, g_y_0_yyyz_yyzzz, g_y_0_yyyz_yzzzz, g_y_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xxxxx[k] = -g_y_0_yyyz_xxxxx[k] * cd_x[k] + g_y_0_yyyz_xxxxxx[k];

                g_y_0_xyyyz_xxxxy[k] = -g_y_0_yyyz_xxxxy[k] * cd_x[k] + g_y_0_yyyz_xxxxxy[k];

                g_y_0_xyyyz_xxxxz[k] = -g_y_0_yyyz_xxxxz[k] * cd_x[k] + g_y_0_yyyz_xxxxxz[k];

                g_y_0_xyyyz_xxxyy[k] = -g_y_0_yyyz_xxxyy[k] * cd_x[k] + g_y_0_yyyz_xxxxyy[k];

                g_y_0_xyyyz_xxxyz[k] = -g_y_0_yyyz_xxxyz[k] * cd_x[k] + g_y_0_yyyz_xxxxyz[k];

                g_y_0_xyyyz_xxxzz[k] = -g_y_0_yyyz_xxxzz[k] * cd_x[k] + g_y_0_yyyz_xxxxzz[k];

                g_y_0_xyyyz_xxyyy[k] = -g_y_0_yyyz_xxyyy[k] * cd_x[k] + g_y_0_yyyz_xxxyyy[k];

                g_y_0_xyyyz_xxyyz[k] = -g_y_0_yyyz_xxyyz[k] * cd_x[k] + g_y_0_yyyz_xxxyyz[k];

                g_y_0_xyyyz_xxyzz[k] = -g_y_0_yyyz_xxyzz[k] * cd_x[k] + g_y_0_yyyz_xxxyzz[k];

                g_y_0_xyyyz_xxzzz[k] = -g_y_0_yyyz_xxzzz[k] * cd_x[k] + g_y_0_yyyz_xxxzzz[k];

                g_y_0_xyyyz_xyyyy[k] = -g_y_0_yyyz_xyyyy[k] * cd_x[k] + g_y_0_yyyz_xxyyyy[k];

                g_y_0_xyyyz_xyyyz[k] = -g_y_0_yyyz_xyyyz[k] * cd_x[k] + g_y_0_yyyz_xxyyyz[k];

                g_y_0_xyyyz_xyyzz[k] = -g_y_0_yyyz_xyyzz[k] * cd_x[k] + g_y_0_yyyz_xxyyzz[k];

                g_y_0_xyyyz_xyzzz[k] = -g_y_0_yyyz_xyzzz[k] * cd_x[k] + g_y_0_yyyz_xxyzzz[k];

                g_y_0_xyyyz_xzzzz[k] = -g_y_0_yyyz_xzzzz[k] * cd_x[k] + g_y_0_yyyz_xxzzzz[k];

                g_y_0_xyyyz_yyyyy[k] = -g_y_0_yyyz_yyyyy[k] * cd_x[k] + g_y_0_yyyz_xyyyyy[k];

                g_y_0_xyyyz_yyyyz[k] = -g_y_0_yyyz_yyyyz[k] * cd_x[k] + g_y_0_yyyz_xyyyyz[k];

                g_y_0_xyyyz_yyyzz[k] = -g_y_0_yyyz_yyyzz[k] * cd_x[k] + g_y_0_yyyz_xyyyzz[k];

                g_y_0_xyyyz_yyzzz[k] = -g_y_0_yyyz_yyzzz[k] * cd_x[k] + g_y_0_yyyz_xyyzzz[k];

                g_y_0_xyyyz_yzzzz[k] = -g_y_0_yyyz_yzzzz[k] * cd_x[k] + g_y_0_yyyz_xyzzzz[k];

                g_y_0_xyyyz_zzzzz[k] = -g_y_0_yyyz_zzzzz[k] * cd_x[k] + g_y_0_yyyz_xzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 252);

            auto g_y_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 253);

            auto g_y_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 254);

            auto g_y_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 255);

            auto g_y_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 256);

            auto g_y_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 257);

            auto g_y_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 258);

            auto g_y_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 259);

            auto g_y_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 260);

            auto g_y_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 261);

            auto g_y_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 262);

            auto g_y_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 263);

            auto g_y_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 264);

            auto g_y_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 265);

            auto g_y_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 266);

            auto g_y_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 267);

            auto g_y_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 268);

            auto g_y_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 269);

            auto g_y_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 270);

            auto g_y_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 271);

            auto g_y_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzz_xxxxx, g_y_0_xyyzz_xxxxy, g_y_0_xyyzz_xxxxz, g_y_0_xyyzz_xxxyy, g_y_0_xyyzz_xxxyz, g_y_0_xyyzz_xxxzz, g_y_0_xyyzz_xxyyy, g_y_0_xyyzz_xxyyz, g_y_0_xyyzz_xxyzz, g_y_0_xyyzz_xxzzz, g_y_0_xyyzz_xyyyy, g_y_0_xyyzz_xyyyz, g_y_0_xyyzz_xyyzz, g_y_0_xyyzz_xyzzz, g_y_0_xyyzz_xzzzz, g_y_0_xyyzz_yyyyy, g_y_0_xyyzz_yyyyz, g_y_0_xyyzz_yyyzz, g_y_0_xyyzz_yyzzz, g_y_0_xyyzz_yzzzz, g_y_0_xyyzz_zzzzz, g_y_0_yyzz_xxxxx, g_y_0_yyzz_xxxxxx, g_y_0_yyzz_xxxxxy, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxy, g_y_0_yyzz_xxxxyy, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxyy, g_y_0_yyzz_xxxyyy, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxyyy, g_y_0_yyzz_xxyyyy, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xyyyy, g_y_0_yyzz_xyyyyy, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_yyyyy, g_y_0_yyzz_yyyyz, g_y_0_yyzz_yyyzz, g_y_0_yyzz_yyzzz, g_y_0_yyzz_yzzzz, g_y_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xxxxx[k] = -g_y_0_yyzz_xxxxx[k] * cd_x[k] + g_y_0_yyzz_xxxxxx[k];

                g_y_0_xyyzz_xxxxy[k] = -g_y_0_yyzz_xxxxy[k] * cd_x[k] + g_y_0_yyzz_xxxxxy[k];

                g_y_0_xyyzz_xxxxz[k] = -g_y_0_yyzz_xxxxz[k] * cd_x[k] + g_y_0_yyzz_xxxxxz[k];

                g_y_0_xyyzz_xxxyy[k] = -g_y_0_yyzz_xxxyy[k] * cd_x[k] + g_y_0_yyzz_xxxxyy[k];

                g_y_0_xyyzz_xxxyz[k] = -g_y_0_yyzz_xxxyz[k] * cd_x[k] + g_y_0_yyzz_xxxxyz[k];

                g_y_0_xyyzz_xxxzz[k] = -g_y_0_yyzz_xxxzz[k] * cd_x[k] + g_y_0_yyzz_xxxxzz[k];

                g_y_0_xyyzz_xxyyy[k] = -g_y_0_yyzz_xxyyy[k] * cd_x[k] + g_y_0_yyzz_xxxyyy[k];

                g_y_0_xyyzz_xxyyz[k] = -g_y_0_yyzz_xxyyz[k] * cd_x[k] + g_y_0_yyzz_xxxyyz[k];

                g_y_0_xyyzz_xxyzz[k] = -g_y_0_yyzz_xxyzz[k] * cd_x[k] + g_y_0_yyzz_xxxyzz[k];

                g_y_0_xyyzz_xxzzz[k] = -g_y_0_yyzz_xxzzz[k] * cd_x[k] + g_y_0_yyzz_xxxzzz[k];

                g_y_0_xyyzz_xyyyy[k] = -g_y_0_yyzz_xyyyy[k] * cd_x[k] + g_y_0_yyzz_xxyyyy[k];

                g_y_0_xyyzz_xyyyz[k] = -g_y_0_yyzz_xyyyz[k] * cd_x[k] + g_y_0_yyzz_xxyyyz[k];

                g_y_0_xyyzz_xyyzz[k] = -g_y_0_yyzz_xyyzz[k] * cd_x[k] + g_y_0_yyzz_xxyyzz[k];

                g_y_0_xyyzz_xyzzz[k] = -g_y_0_yyzz_xyzzz[k] * cd_x[k] + g_y_0_yyzz_xxyzzz[k];

                g_y_0_xyyzz_xzzzz[k] = -g_y_0_yyzz_xzzzz[k] * cd_x[k] + g_y_0_yyzz_xxzzzz[k];

                g_y_0_xyyzz_yyyyy[k] = -g_y_0_yyzz_yyyyy[k] * cd_x[k] + g_y_0_yyzz_xyyyyy[k];

                g_y_0_xyyzz_yyyyz[k] = -g_y_0_yyzz_yyyyz[k] * cd_x[k] + g_y_0_yyzz_xyyyyz[k];

                g_y_0_xyyzz_yyyzz[k] = -g_y_0_yyzz_yyyzz[k] * cd_x[k] + g_y_0_yyzz_xyyyzz[k];

                g_y_0_xyyzz_yyzzz[k] = -g_y_0_yyzz_yyzzz[k] * cd_x[k] + g_y_0_yyzz_xyyzzz[k];

                g_y_0_xyyzz_yzzzz[k] = -g_y_0_yyzz_yzzzz[k] * cd_x[k] + g_y_0_yyzz_xyzzzz[k];

                g_y_0_xyyzz_zzzzz[k] = -g_y_0_yyzz_zzzzz[k] * cd_x[k] + g_y_0_yyzz_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 273);

            auto g_y_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 274);

            auto g_y_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 275);

            auto g_y_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 276);

            auto g_y_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 277);

            auto g_y_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 278);

            auto g_y_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 279);

            auto g_y_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 280);

            auto g_y_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 281);

            auto g_y_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 282);

            auto g_y_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 283);

            auto g_y_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 284);

            auto g_y_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 285);

            auto g_y_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 286);

            auto g_y_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 287);

            auto g_y_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 288);

            auto g_y_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 289);

            auto g_y_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 290);

            auto g_y_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 291);

            auto g_y_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 292);

            auto g_y_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzz_xxxxx, g_y_0_xyzzz_xxxxy, g_y_0_xyzzz_xxxxz, g_y_0_xyzzz_xxxyy, g_y_0_xyzzz_xxxyz, g_y_0_xyzzz_xxxzz, g_y_0_xyzzz_xxyyy, g_y_0_xyzzz_xxyyz, g_y_0_xyzzz_xxyzz, g_y_0_xyzzz_xxzzz, g_y_0_xyzzz_xyyyy, g_y_0_xyzzz_xyyyz, g_y_0_xyzzz_xyyzz, g_y_0_xyzzz_xyzzz, g_y_0_xyzzz_xzzzz, g_y_0_xyzzz_yyyyy, g_y_0_xyzzz_yyyyz, g_y_0_xyzzz_yyyzz, g_y_0_xyzzz_yyzzz, g_y_0_xyzzz_yzzzz, g_y_0_xyzzz_zzzzz, g_y_0_yzzz_xxxxx, g_y_0_yzzz_xxxxxx, g_y_0_yzzz_xxxxxy, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxy, g_y_0_yzzz_xxxxyy, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxyy, g_y_0_yzzz_xxxyyy, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxyyy, g_y_0_yzzz_xxyyyy, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xyyyy, g_y_0_yzzz_xyyyyy, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_yyyyy, g_y_0_yzzz_yyyyz, g_y_0_yzzz_yyyzz, g_y_0_yzzz_yyzzz, g_y_0_yzzz_yzzzz, g_y_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xxxxx[k] = -g_y_0_yzzz_xxxxx[k] * cd_x[k] + g_y_0_yzzz_xxxxxx[k];

                g_y_0_xyzzz_xxxxy[k] = -g_y_0_yzzz_xxxxy[k] * cd_x[k] + g_y_0_yzzz_xxxxxy[k];

                g_y_0_xyzzz_xxxxz[k] = -g_y_0_yzzz_xxxxz[k] * cd_x[k] + g_y_0_yzzz_xxxxxz[k];

                g_y_0_xyzzz_xxxyy[k] = -g_y_0_yzzz_xxxyy[k] * cd_x[k] + g_y_0_yzzz_xxxxyy[k];

                g_y_0_xyzzz_xxxyz[k] = -g_y_0_yzzz_xxxyz[k] * cd_x[k] + g_y_0_yzzz_xxxxyz[k];

                g_y_0_xyzzz_xxxzz[k] = -g_y_0_yzzz_xxxzz[k] * cd_x[k] + g_y_0_yzzz_xxxxzz[k];

                g_y_0_xyzzz_xxyyy[k] = -g_y_0_yzzz_xxyyy[k] * cd_x[k] + g_y_0_yzzz_xxxyyy[k];

                g_y_0_xyzzz_xxyyz[k] = -g_y_0_yzzz_xxyyz[k] * cd_x[k] + g_y_0_yzzz_xxxyyz[k];

                g_y_0_xyzzz_xxyzz[k] = -g_y_0_yzzz_xxyzz[k] * cd_x[k] + g_y_0_yzzz_xxxyzz[k];

                g_y_0_xyzzz_xxzzz[k] = -g_y_0_yzzz_xxzzz[k] * cd_x[k] + g_y_0_yzzz_xxxzzz[k];

                g_y_0_xyzzz_xyyyy[k] = -g_y_0_yzzz_xyyyy[k] * cd_x[k] + g_y_0_yzzz_xxyyyy[k];

                g_y_0_xyzzz_xyyyz[k] = -g_y_0_yzzz_xyyyz[k] * cd_x[k] + g_y_0_yzzz_xxyyyz[k];

                g_y_0_xyzzz_xyyzz[k] = -g_y_0_yzzz_xyyzz[k] * cd_x[k] + g_y_0_yzzz_xxyyzz[k];

                g_y_0_xyzzz_xyzzz[k] = -g_y_0_yzzz_xyzzz[k] * cd_x[k] + g_y_0_yzzz_xxyzzz[k];

                g_y_0_xyzzz_xzzzz[k] = -g_y_0_yzzz_xzzzz[k] * cd_x[k] + g_y_0_yzzz_xxzzzz[k];

                g_y_0_xyzzz_yyyyy[k] = -g_y_0_yzzz_yyyyy[k] * cd_x[k] + g_y_0_yzzz_xyyyyy[k];

                g_y_0_xyzzz_yyyyz[k] = -g_y_0_yzzz_yyyyz[k] * cd_x[k] + g_y_0_yzzz_xyyyyz[k];

                g_y_0_xyzzz_yyyzz[k] = -g_y_0_yzzz_yyyzz[k] * cd_x[k] + g_y_0_yzzz_xyyyzz[k];

                g_y_0_xyzzz_yyzzz[k] = -g_y_0_yzzz_yyzzz[k] * cd_x[k] + g_y_0_yzzz_xyyzzz[k];

                g_y_0_xyzzz_yzzzz[k] = -g_y_0_yzzz_yzzzz[k] * cd_x[k] + g_y_0_yzzz_xyzzzz[k];

                g_y_0_xyzzz_zzzzz[k] = -g_y_0_yzzz_zzzzz[k] * cd_x[k] + g_y_0_yzzz_xzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 294);

            auto g_y_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 295);

            auto g_y_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 296);

            auto g_y_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 297);

            auto g_y_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 298);

            auto g_y_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 299);

            auto g_y_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 300);

            auto g_y_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 301);

            auto g_y_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 302);

            auto g_y_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 303);

            auto g_y_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 304);

            auto g_y_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 305);

            auto g_y_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 306);

            auto g_y_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 307);

            auto g_y_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 308);

            auto g_y_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 309);

            auto g_y_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 310);

            auto g_y_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 311);

            auto g_y_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 312);

            auto g_y_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 313);

            auto g_y_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzz_xxxxx, g_y_0_xzzzz_xxxxy, g_y_0_xzzzz_xxxxz, g_y_0_xzzzz_xxxyy, g_y_0_xzzzz_xxxyz, g_y_0_xzzzz_xxxzz, g_y_0_xzzzz_xxyyy, g_y_0_xzzzz_xxyyz, g_y_0_xzzzz_xxyzz, g_y_0_xzzzz_xxzzz, g_y_0_xzzzz_xyyyy, g_y_0_xzzzz_xyyyz, g_y_0_xzzzz_xyyzz, g_y_0_xzzzz_xyzzz, g_y_0_xzzzz_xzzzz, g_y_0_xzzzz_yyyyy, g_y_0_xzzzz_yyyyz, g_y_0_xzzzz_yyyzz, g_y_0_xzzzz_yyzzz, g_y_0_xzzzz_yzzzz, g_y_0_xzzzz_zzzzz, g_y_0_zzzz_xxxxx, g_y_0_zzzz_xxxxxx, g_y_0_zzzz_xxxxxy, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxy, g_y_0_zzzz_xxxxyy, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxyy, g_y_0_zzzz_xxxyyy, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxyyy, g_y_0_zzzz_xxyyyy, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xyyyy, g_y_0_zzzz_xyyyyy, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_yyyyy, g_y_0_zzzz_yyyyz, g_y_0_zzzz_yyyzz, g_y_0_zzzz_yyzzz, g_y_0_zzzz_yzzzz, g_y_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xxxxx[k] = -g_y_0_zzzz_xxxxx[k] * cd_x[k] + g_y_0_zzzz_xxxxxx[k];

                g_y_0_xzzzz_xxxxy[k] = -g_y_0_zzzz_xxxxy[k] * cd_x[k] + g_y_0_zzzz_xxxxxy[k];

                g_y_0_xzzzz_xxxxz[k] = -g_y_0_zzzz_xxxxz[k] * cd_x[k] + g_y_0_zzzz_xxxxxz[k];

                g_y_0_xzzzz_xxxyy[k] = -g_y_0_zzzz_xxxyy[k] * cd_x[k] + g_y_0_zzzz_xxxxyy[k];

                g_y_0_xzzzz_xxxyz[k] = -g_y_0_zzzz_xxxyz[k] * cd_x[k] + g_y_0_zzzz_xxxxyz[k];

                g_y_0_xzzzz_xxxzz[k] = -g_y_0_zzzz_xxxzz[k] * cd_x[k] + g_y_0_zzzz_xxxxzz[k];

                g_y_0_xzzzz_xxyyy[k] = -g_y_0_zzzz_xxyyy[k] * cd_x[k] + g_y_0_zzzz_xxxyyy[k];

                g_y_0_xzzzz_xxyyz[k] = -g_y_0_zzzz_xxyyz[k] * cd_x[k] + g_y_0_zzzz_xxxyyz[k];

                g_y_0_xzzzz_xxyzz[k] = -g_y_0_zzzz_xxyzz[k] * cd_x[k] + g_y_0_zzzz_xxxyzz[k];

                g_y_0_xzzzz_xxzzz[k] = -g_y_0_zzzz_xxzzz[k] * cd_x[k] + g_y_0_zzzz_xxxzzz[k];

                g_y_0_xzzzz_xyyyy[k] = -g_y_0_zzzz_xyyyy[k] * cd_x[k] + g_y_0_zzzz_xxyyyy[k];

                g_y_0_xzzzz_xyyyz[k] = -g_y_0_zzzz_xyyyz[k] * cd_x[k] + g_y_0_zzzz_xxyyyz[k];

                g_y_0_xzzzz_xyyzz[k] = -g_y_0_zzzz_xyyzz[k] * cd_x[k] + g_y_0_zzzz_xxyyzz[k];

                g_y_0_xzzzz_xyzzz[k] = -g_y_0_zzzz_xyzzz[k] * cd_x[k] + g_y_0_zzzz_xxyzzz[k];

                g_y_0_xzzzz_xzzzz[k] = -g_y_0_zzzz_xzzzz[k] * cd_x[k] + g_y_0_zzzz_xxzzzz[k];

                g_y_0_xzzzz_yyyyy[k] = -g_y_0_zzzz_yyyyy[k] * cd_x[k] + g_y_0_zzzz_xyyyyy[k];

                g_y_0_xzzzz_yyyyz[k] = -g_y_0_zzzz_yyyyz[k] * cd_x[k] + g_y_0_zzzz_xyyyyz[k];

                g_y_0_xzzzz_yyyzz[k] = -g_y_0_zzzz_yyyzz[k] * cd_x[k] + g_y_0_zzzz_xyyyzz[k];

                g_y_0_xzzzz_yyzzz[k] = -g_y_0_zzzz_yyzzz[k] * cd_x[k] + g_y_0_zzzz_xyyzzz[k];

                g_y_0_xzzzz_yzzzz[k] = -g_y_0_zzzz_yzzzz[k] * cd_x[k] + g_y_0_zzzz_xyzzzz[k];

                g_y_0_xzzzz_zzzzz[k] = -g_y_0_zzzz_zzzzz[k] * cd_x[k] + g_y_0_zzzz_xzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 315);

            auto g_y_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 316);

            auto g_y_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 317);

            auto g_y_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 318);

            auto g_y_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 319);

            auto g_y_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 320);

            auto g_y_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 321);

            auto g_y_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 322);

            auto g_y_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 323);

            auto g_y_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 324);

            auto g_y_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 325);

            auto g_y_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 326);

            auto g_y_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 327);

            auto g_y_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 328);

            auto g_y_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 329);

            auto g_y_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 330);

            auto g_y_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 331);

            auto g_y_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 332);

            auto g_y_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 333);

            auto g_y_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 334);

            auto g_y_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_y, g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_zzzzz, g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_zzzzz, g_yyyy_xxxxx, g_yyyy_xxxxy, g_yyyy_xxxxz, g_yyyy_xxxyy, g_yyyy_xxxyz, g_yyyy_xxxzz, g_yyyy_xxyyy, g_yyyy_xxyyz, g_yyyy_xxyzz, g_yyyy_xxzzz, g_yyyy_xyyyy, g_yyyy_xyyyz, g_yyyy_xyyzz, g_yyyy_xyzzz, g_yyyy_xzzzz, g_yyyy_yyyyy, g_yyyy_yyyyz, g_yyyy_yyyzz, g_yyyy_yyzzz, g_yyyy_yzzzz, g_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xxxxx[k] = -g_yyyy_xxxxx[k] - g_y_0_yyyy_xxxxx[k] * cd_y[k] + g_y_0_yyyy_xxxxxy[k];

                g_y_0_yyyyy_xxxxy[k] = -g_yyyy_xxxxy[k] - g_y_0_yyyy_xxxxy[k] * cd_y[k] + g_y_0_yyyy_xxxxyy[k];

                g_y_0_yyyyy_xxxxz[k] = -g_yyyy_xxxxz[k] - g_y_0_yyyy_xxxxz[k] * cd_y[k] + g_y_0_yyyy_xxxxyz[k];

                g_y_0_yyyyy_xxxyy[k] = -g_yyyy_xxxyy[k] - g_y_0_yyyy_xxxyy[k] * cd_y[k] + g_y_0_yyyy_xxxyyy[k];

                g_y_0_yyyyy_xxxyz[k] = -g_yyyy_xxxyz[k] - g_y_0_yyyy_xxxyz[k] * cd_y[k] + g_y_0_yyyy_xxxyyz[k];

                g_y_0_yyyyy_xxxzz[k] = -g_yyyy_xxxzz[k] - g_y_0_yyyy_xxxzz[k] * cd_y[k] + g_y_0_yyyy_xxxyzz[k];

                g_y_0_yyyyy_xxyyy[k] = -g_yyyy_xxyyy[k] - g_y_0_yyyy_xxyyy[k] * cd_y[k] + g_y_0_yyyy_xxyyyy[k];

                g_y_0_yyyyy_xxyyz[k] = -g_yyyy_xxyyz[k] - g_y_0_yyyy_xxyyz[k] * cd_y[k] + g_y_0_yyyy_xxyyyz[k];

                g_y_0_yyyyy_xxyzz[k] = -g_yyyy_xxyzz[k] - g_y_0_yyyy_xxyzz[k] * cd_y[k] + g_y_0_yyyy_xxyyzz[k];

                g_y_0_yyyyy_xxzzz[k] = -g_yyyy_xxzzz[k] - g_y_0_yyyy_xxzzz[k] * cd_y[k] + g_y_0_yyyy_xxyzzz[k];

                g_y_0_yyyyy_xyyyy[k] = -g_yyyy_xyyyy[k] - g_y_0_yyyy_xyyyy[k] * cd_y[k] + g_y_0_yyyy_xyyyyy[k];

                g_y_0_yyyyy_xyyyz[k] = -g_yyyy_xyyyz[k] - g_y_0_yyyy_xyyyz[k] * cd_y[k] + g_y_0_yyyy_xyyyyz[k];

                g_y_0_yyyyy_xyyzz[k] = -g_yyyy_xyyzz[k] - g_y_0_yyyy_xyyzz[k] * cd_y[k] + g_y_0_yyyy_xyyyzz[k];

                g_y_0_yyyyy_xyzzz[k] = -g_yyyy_xyzzz[k] - g_y_0_yyyy_xyzzz[k] * cd_y[k] + g_y_0_yyyy_xyyzzz[k];

                g_y_0_yyyyy_xzzzz[k] = -g_yyyy_xzzzz[k] - g_y_0_yyyy_xzzzz[k] * cd_y[k] + g_y_0_yyyy_xyzzzz[k];

                g_y_0_yyyyy_yyyyy[k] = -g_yyyy_yyyyy[k] - g_y_0_yyyy_yyyyy[k] * cd_y[k] + g_y_0_yyyy_yyyyyy[k];

                g_y_0_yyyyy_yyyyz[k] = -g_yyyy_yyyyz[k] - g_y_0_yyyy_yyyyz[k] * cd_y[k] + g_y_0_yyyy_yyyyyz[k];

                g_y_0_yyyyy_yyyzz[k] = -g_yyyy_yyyzz[k] - g_y_0_yyyy_yyyzz[k] * cd_y[k] + g_y_0_yyyy_yyyyzz[k];

                g_y_0_yyyyy_yyzzz[k] = -g_yyyy_yyzzz[k] - g_y_0_yyyy_yyzzz[k] * cd_y[k] + g_y_0_yyyy_yyyzzz[k];

                g_y_0_yyyyy_yzzzz[k] = -g_yyyy_yzzzz[k] - g_y_0_yyyy_yzzzz[k] * cd_y[k] + g_y_0_yyyy_yyzzzz[k];

                g_y_0_yyyyy_zzzzz[k] = -g_yyyy_zzzzz[k] - g_y_0_yyyy_zzzzz[k] * cd_y[k] + g_y_0_yyyy_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 336);

            auto g_y_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 337);

            auto g_y_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 338);

            auto g_y_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 339);

            auto g_y_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 340);

            auto g_y_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 341);

            auto g_y_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 342);

            auto g_y_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 343);

            auto g_y_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 344);

            auto g_y_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 345);

            auto g_y_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 346);

            auto g_y_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 347);

            auto g_y_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 348);

            auto g_y_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 349);

            auto g_y_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 350);

            auto g_y_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 351);

            auto g_y_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 352);

            auto g_y_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 353);

            auto g_y_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 354);

            auto g_y_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 355);

            auto g_y_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 356);

            #pragma omp simd aligned(cd_z, g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_zzzzz, g_y_0_yyyy_zzzzzz, g_y_0_yyyyz_xxxxx, g_y_0_yyyyz_xxxxy, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxyy, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxyyy, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xyyyy, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_yyyyy, g_y_0_yyyyz_yyyyz, g_y_0_yyyyz_yyyzz, g_y_0_yyyyz_yyzzz, g_y_0_yyyyz_yzzzz, g_y_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xxxxx[k] = -g_y_0_yyyy_xxxxx[k] * cd_z[k] + g_y_0_yyyy_xxxxxz[k];

                g_y_0_yyyyz_xxxxy[k] = -g_y_0_yyyy_xxxxy[k] * cd_z[k] + g_y_0_yyyy_xxxxyz[k];

                g_y_0_yyyyz_xxxxz[k] = -g_y_0_yyyy_xxxxz[k] * cd_z[k] + g_y_0_yyyy_xxxxzz[k];

                g_y_0_yyyyz_xxxyy[k] = -g_y_0_yyyy_xxxyy[k] * cd_z[k] + g_y_0_yyyy_xxxyyz[k];

                g_y_0_yyyyz_xxxyz[k] = -g_y_0_yyyy_xxxyz[k] * cd_z[k] + g_y_0_yyyy_xxxyzz[k];

                g_y_0_yyyyz_xxxzz[k] = -g_y_0_yyyy_xxxzz[k] * cd_z[k] + g_y_0_yyyy_xxxzzz[k];

                g_y_0_yyyyz_xxyyy[k] = -g_y_0_yyyy_xxyyy[k] * cd_z[k] + g_y_0_yyyy_xxyyyz[k];

                g_y_0_yyyyz_xxyyz[k] = -g_y_0_yyyy_xxyyz[k] * cd_z[k] + g_y_0_yyyy_xxyyzz[k];

                g_y_0_yyyyz_xxyzz[k] = -g_y_0_yyyy_xxyzz[k] * cd_z[k] + g_y_0_yyyy_xxyzzz[k];

                g_y_0_yyyyz_xxzzz[k] = -g_y_0_yyyy_xxzzz[k] * cd_z[k] + g_y_0_yyyy_xxzzzz[k];

                g_y_0_yyyyz_xyyyy[k] = -g_y_0_yyyy_xyyyy[k] * cd_z[k] + g_y_0_yyyy_xyyyyz[k];

                g_y_0_yyyyz_xyyyz[k] = -g_y_0_yyyy_xyyyz[k] * cd_z[k] + g_y_0_yyyy_xyyyzz[k];

                g_y_0_yyyyz_xyyzz[k] = -g_y_0_yyyy_xyyzz[k] * cd_z[k] + g_y_0_yyyy_xyyzzz[k];

                g_y_0_yyyyz_xyzzz[k] = -g_y_0_yyyy_xyzzz[k] * cd_z[k] + g_y_0_yyyy_xyzzzz[k];

                g_y_0_yyyyz_xzzzz[k] = -g_y_0_yyyy_xzzzz[k] * cd_z[k] + g_y_0_yyyy_xzzzzz[k];

                g_y_0_yyyyz_yyyyy[k] = -g_y_0_yyyy_yyyyy[k] * cd_z[k] + g_y_0_yyyy_yyyyyz[k];

                g_y_0_yyyyz_yyyyz[k] = -g_y_0_yyyy_yyyyz[k] * cd_z[k] + g_y_0_yyyy_yyyyzz[k];

                g_y_0_yyyyz_yyyzz[k] = -g_y_0_yyyy_yyyzz[k] * cd_z[k] + g_y_0_yyyy_yyyzzz[k];

                g_y_0_yyyyz_yyzzz[k] = -g_y_0_yyyy_yyzzz[k] * cd_z[k] + g_y_0_yyyy_yyzzzz[k];

                g_y_0_yyyyz_yzzzz[k] = -g_y_0_yyyy_yzzzz[k] * cd_z[k] + g_y_0_yyyy_yzzzzz[k];

                g_y_0_yyyyz_zzzzz[k] = -g_y_0_yyyy_zzzzz[k] * cd_z[k] + g_y_0_yyyy_zzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 357);

            auto g_y_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 358);

            auto g_y_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 359);

            auto g_y_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 360);

            auto g_y_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 361);

            auto g_y_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 362);

            auto g_y_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 363);

            auto g_y_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 364);

            auto g_y_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 365);

            auto g_y_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 366);

            auto g_y_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 367);

            auto g_y_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 368);

            auto g_y_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 369);

            auto g_y_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 370);

            auto g_y_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 371);

            auto g_y_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 372);

            auto g_y_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 373);

            auto g_y_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 374);

            auto g_y_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 375);

            auto g_y_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 376);

            auto g_y_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 377);

            #pragma omp simd aligned(cd_z, g_y_0_yyyz_xxxxx, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxy, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxyy, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxyyy, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xyyyy, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_yyyyy, g_y_0_yyyz_yyyyyz, g_y_0_yyyz_yyyyz, g_y_0_yyyz_yyyyzz, g_y_0_yyyz_yyyzz, g_y_0_yyyz_yyyzzz, g_y_0_yyyz_yyzzz, g_y_0_yyyz_yyzzzz, g_y_0_yyyz_yzzzz, g_y_0_yyyz_yzzzzz, g_y_0_yyyz_zzzzz, g_y_0_yyyz_zzzzzz, g_y_0_yyyzz_xxxxx, g_y_0_yyyzz_xxxxy, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxyy, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxyyy, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xyyyy, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_yyyyy, g_y_0_yyyzz_yyyyz, g_y_0_yyyzz_yyyzz, g_y_0_yyyzz_yyzzz, g_y_0_yyyzz_yzzzz, g_y_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xxxxx[k] = -g_y_0_yyyz_xxxxx[k] * cd_z[k] + g_y_0_yyyz_xxxxxz[k];

                g_y_0_yyyzz_xxxxy[k] = -g_y_0_yyyz_xxxxy[k] * cd_z[k] + g_y_0_yyyz_xxxxyz[k];

                g_y_0_yyyzz_xxxxz[k] = -g_y_0_yyyz_xxxxz[k] * cd_z[k] + g_y_0_yyyz_xxxxzz[k];

                g_y_0_yyyzz_xxxyy[k] = -g_y_0_yyyz_xxxyy[k] * cd_z[k] + g_y_0_yyyz_xxxyyz[k];

                g_y_0_yyyzz_xxxyz[k] = -g_y_0_yyyz_xxxyz[k] * cd_z[k] + g_y_0_yyyz_xxxyzz[k];

                g_y_0_yyyzz_xxxzz[k] = -g_y_0_yyyz_xxxzz[k] * cd_z[k] + g_y_0_yyyz_xxxzzz[k];

                g_y_0_yyyzz_xxyyy[k] = -g_y_0_yyyz_xxyyy[k] * cd_z[k] + g_y_0_yyyz_xxyyyz[k];

                g_y_0_yyyzz_xxyyz[k] = -g_y_0_yyyz_xxyyz[k] * cd_z[k] + g_y_0_yyyz_xxyyzz[k];

                g_y_0_yyyzz_xxyzz[k] = -g_y_0_yyyz_xxyzz[k] * cd_z[k] + g_y_0_yyyz_xxyzzz[k];

                g_y_0_yyyzz_xxzzz[k] = -g_y_0_yyyz_xxzzz[k] * cd_z[k] + g_y_0_yyyz_xxzzzz[k];

                g_y_0_yyyzz_xyyyy[k] = -g_y_0_yyyz_xyyyy[k] * cd_z[k] + g_y_0_yyyz_xyyyyz[k];

                g_y_0_yyyzz_xyyyz[k] = -g_y_0_yyyz_xyyyz[k] * cd_z[k] + g_y_0_yyyz_xyyyzz[k];

                g_y_0_yyyzz_xyyzz[k] = -g_y_0_yyyz_xyyzz[k] * cd_z[k] + g_y_0_yyyz_xyyzzz[k];

                g_y_0_yyyzz_xyzzz[k] = -g_y_0_yyyz_xyzzz[k] * cd_z[k] + g_y_0_yyyz_xyzzzz[k];

                g_y_0_yyyzz_xzzzz[k] = -g_y_0_yyyz_xzzzz[k] * cd_z[k] + g_y_0_yyyz_xzzzzz[k];

                g_y_0_yyyzz_yyyyy[k] = -g_y_0_yyyz_yyyyy[k] * cd_z[k] + g_y_0_yyyz_yyyyyz[k];

                g_y_0_yyyzz_yyyyz[k] = -g_y_0_yyyz_yyyyz[k] * cd_z[k] + g_y_0_yyyz_yyyyzz[k];

                g_y_0_yyyzz_yyyzz[k] = -g_y_0_yyyz_yyyzz[k] * cd_z[k] + g_y_0_yyyz_yyyzzz[k];

                g_y_0_yyyzz_yyzzz[k] = -g_y_0_yyyz_yyzzz[k] * cd_z[k] + g_y_0_yyyz_yyzzzz[k];

                g_y_0_yyyzz_yzzzz[k] = -g_y_0_yyyz_yzzzz[k] * cd_z[k] + g_y_0_yyyz_yzzzzz[k];

                g_y_0_yyyzz_zzzzz[k] = -g_y_0_yyyz_zzzzz[k] * cd_z[k] + g_y_0_yyyz_zzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 378);

            auto g_y_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 379);

            auto g_y_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 380);

            auto g_y_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 381);

            auto g_y_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 382);

            auto g_y_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 383);

            auto g_y_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 384);

            auto g_y_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 385);

            auto g_y_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 386);

            auto g_y_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 387);

            auto g_y_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 388);

            auto g_y_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 389);

            auto g_y_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 390);

            auto g_y_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 391);

            auto g_y_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 392);

            auto g_y_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 393);

            auto g_y_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 394);

            auto g_y_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 395);

            auto g_y_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 396);

            auto g_y_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 397);

            auto g_y_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 398);

            #pragma omp simd aligned(cd_z, g_y_0_yyzz_xxxxx, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxy, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxyy, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxyyy, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xyyyy, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_yyyyy, g_y_0_yyzz_yyyyyz, g_y_0_yyzz_yyyyz, g_y_0_yyzz_yyyyzz, g_y_0_yyzz_yyyzz, g_y_0_yyzz_yyyzzz, g_y_0_yyzz_yyzzz, g_y_0_yyzz_yyzzzz, g_y_0_yyzz_yzzzz, g_y_0_yyzz_yzzzzz, g_y_0_yyzz_zzzzz, g_y_0_yyzz_zzzzzz, g_y_0_yyzzz_xxxxx, g_y_0_yyzzz_xxxxy, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxyy, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxyyy, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xyyyy, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_yyyyy, g_y_0_yyzzz_yyyyz, g_y_0_yyzzz_yyyzz, g_y_0_yyzzz_yyzzz, g_y_0_yyzzz_yzzzz, g_y_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xxxxx[k] = -g_y_0_yyzz_xxxxx[k] * cd_z[k] + g_y_0_yyzz_xxxxxz[k];

                g_y_0_yyzzz_xxxxy[k] = -g_y_0_yyzz_xxxxy[k] * cd_z[k] + g_y_0_yyzz_xxxxyz[k];

                g_y_0_yyzzz_xxxxz[k] = -g_y_0_yyzz_xxxxz[k] * cd_z[k] + g_y_0_yyzz_xxxxzz[k];

                g_y_0_yyzzz_xxxyy[k] = -g_y_0_yyzz_xxxyy[k] * cd_z[k] + g_y_0_yyzz_xxxyyz[k];

                g_y_0_yyzzz_xxxyz[k] = -g_y_0_yyzz_xxxyz[k] * cd_z[k] + g_y_0_yyzz_xxxyzz[k];

                g_y_0_yyzzz_xxxzz[k] = -g_y_0_yyzz_xxxzz[k] * cd_z[k] + g_y_0_yyzz_xxxzzz[k];

                g_y_0_yyzzz_xxyyy[k] = -g_y_0_yyzz_xxyyy[k] * cd_z[k] + g_y_0_yyzz_xxyyyz[k];

                g_y_0_yyzzz_xxyyz[k] = -g_y_0_yyzz_xxyyz[k] * cd_z[k] + g_y_0_yyzz_xxyyzz[k];

                g_y_0_yyzzz_xxyzz[k] = -g_y_0_yyzz_xxyzz[k] * cd_z[k] + g_y_0_yyzz_xxyzzz[k];

                g_y_0_yyzzz_xxzzz[k] = -g_y_0_yyzz_xxzzz[k] * cd_z[k] + g_y_0_yyzz_xxzzzz[k];

                g_y_0_yyzzz_xyyyy[k] = -g_y_0_yyzz_xyyyy[k] * cd_z[k] + g_y_0_yyzz_xyyyyz[k];

                g_y_0_yyzzz_xyyyz[k] = -g_y_0_yyzz_xyyyz[k] * cd_z[k] + g_y_0_yyzz_xyyyzz[k];

                g_y_0_yyzzz_xyyzz[k] = -g_y_0_yyzz_xyyzz[k] * cd_z[k] + g_y_0_yyzz_xyyzzz[k];

                g_y_0_yyzzz_xyzzz[k] = -g_y_0_yyzz_xyzzz[k] * cd_z[k] + g_y_0_yyzz_xyzzzz[k];

                g_y_0_yyzzz_xzzzz[k] = -g_y_0_yyzz_xzzzz[k] * cd_z[k] + g_y_0_yyzz_xzzzzz[k];

                g_y_0_yyzzz_yyyyy[k] = -g_y_0_yyzz_yyyyy[k] * cd_z[k] + g_y_0_yyzz_yyyyyz[k];

                g_y_0_yyzzz_yyyyz[k] = -g_y_0_yyzz_yyyyz[k] * cd_z[k] + g_y_0_yyzz_yyyyzz[k];

                g_y_0_yyzzz_yyyzz[k] = -g_y_0_yyzz_yyyzz[k] * cd_z[k] + g_y_0_yyzz_yyyzzz[k];

                g_y_0_yyzzz_yyzzz[k] = -g_y_0_yyzz_yyzzz[k] * cd_z[k] + g_y_0_yyzz_yyzzzz[k];

                g_y_0_yyzzz_yzzzz[k] = -g_y_0_yyzz_yzzzz[k] * cd_z[k] + g_y_0_yyzz_yzzzzz[k];

                g_y_0_yyzzz_zzzzz[k] = -g_y_0_yyzz_zzzzz[k] * cd_z[k] + g_y_0_yyzz_zzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 399);

            auto g_y_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 400);

            auto g_y_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 401);

            auto g_y_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 402);

            auto g_y_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 403);

            auto g_y_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 404);

            auto g_y_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 405);

            auto g_y_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 406);

            auto g_y_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 407);

            auto g_y_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 408);

            auto g_y_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 409);

            auto g_y_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 410);

            auto g_y_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 411);

            auto g_y_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 412);

            auto g_y_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 413);

            auto g_y_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 414);

            auto g_y_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 415);

            auto g_y_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 416);

            auto g_y_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 417);

            auto g_y_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 418);

            auto g_y_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_z, g_y_0_yzzz_xxxxx, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxy, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxyy, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxyyy, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xyyyy, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_yyyyy, g_y_0_yzzz_yyyyyz, g_y_0_yzzz_yyyyz, g_y_0_yzzz_yyyyzz, g_y_0_yzzz_yyyzz, g_y_0_yzzz_yyyzzz, g_y_0_yzzz_yyzzz, g_y_0_yzzz_yyzzzz, g_y_0_yzzz_yzzzz, g_y_0_yzzz_yzzzzz, g_y_0_yzzz_zzzzz, g_y_0_yzzz_zzzzzz, g_y_0_yzzzz_xxxxx, g_y_0_yzzzz_xxxxy, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxyy, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxyyy, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xyyyy, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_yyyyy, g_y_0_yzzzz_yyyyz, g_y_0_yzzzz_yyyzz, g_y_0_yzzzz_yyzzz, g_y_0_yzzzz_yzzzz, g_y_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xxxxx[k] = -g_y_0_yzzz_xxxxx[k] * cd_z[k] + g_y_0_yzzz_xxxxxz[k];

                g_y_0_yzzzz_xxxxy[k] = -g_y_0_yzzz_xxxxy[k] * cd_z[k] + g_y_0_yzzz_xxxxyz[k];

                g_y_0_yzzzz_xxxxz[k] = -g_y_0_yzzz_xxxxz[k] * cd_z[k] + g_y_0_yzzz_xxxxzz[k];

                g_y_0_yzzzz_xxxyy[k] = -g_y_0_yzzz_xxxyy[k] * cd_z[k] + g_y_0_yzzz_xxxyyz[k];

                g_y_0_yzzzz_xxxyz[k] = -g_y_0_yzzz_xxxyz[k] * cd_z[k] + g_y_0_yzzz_xxxyzz[k];

                g_y_0_yzzzz_xxxzz[k] = -g_y_0_yzzz_xxxzz[k] * cd_z[k] + g_y_0_yzzz_xxxzzz[k];

                g_y_0_yzzzz_xxyyy[k] = -g_y_0_yzzz_xxyyy[k] * cd_z[k] + g_y_0_yzzz_xxyyyz[k];

                g_y_0_yzzzz_xxyyz[k] = -g_y_0_yzzz_xxyyz[k] * cd_z[k] + g_y_0_yzzz_xxyyzz[k];

                g_y_0_yzzzz_xxyzz[k] = -g_y_0_yzzz_xxyzz[k] * cd_z[k] + g_y_0_yzzz_xxyzzz[k];

                g_y_0_yzzzz_xxzzz[k] = -g_y_0_yzzz_xxzzz[k] * cd_z[k] + g_y_0_yzzz_xxzzzz[k];

                g_y_0_yzzzz_xyyyy[k] = -g_y_0_yzzz_xyyyy[k] * cd_z[k] + g_y_0_yzzz_xyyyyz[k];

                g_y_0_yzzzz_xyyyz[k] = -g_y_0_yzzz_xyyyz[k] * cd_z[k] + g_y_0_yzzz_xyyyzz[k];

                g_y_0_yzzzz_xyyzz[k] = -g_y_0_yzzz_xyyzz[k] * cd_z[k] + g_y_0_yzzz_xyyzzz[k];

                g_y_0_yzzzz_xyzzz[k] = -g_y_0_yzzz_xyzzz[k] * cd_z[k] + g_y_0_yzzz_xyzzzz[k];

                g_y_0_yzzzz_xzzzz[k] = -g_y_0_yzzz_xzzzz[k] * cd_z[k] + g_y_0_yzzz_xzzzzz[k];

                g_y_0_yzzzz_yyyyy[k] = -g_y_0_yzzz_yyyyy[k] * cd_z[k] + g_y_0_yzzz_yyyyyz[k];

                g_y_0_yzzzz_yyyyz[k] = -g_y_0_yzzz_yyyyz[k] * cd_z[k] + g_y_0_yzzz_yyyyzz[k];

                g_y_0_yzzzz_yyyzz[k] = -g_y_0_yzzz_yyyzz[k] * cd_z[k] + g_y_0_yzzz_yyyzzz[k];

                g_y_0_yzzzz_yyzzz[k] = -g_y_0_yzzz_yyzzz[k] * cd_z[k] + g_y_0_yzzz_yyzzzz[k];

                g_y_0_yzzzz_yzzzz[k] = -g_y_0_yzzz_yzzzz[k] * cd_z[k] + g_y_0_yzzz_yzzzzz[k];

                g_y_0_yzzzz_zzzzz[k] = -g_y_0_yzzz_zzzzz[k] * cd_z[k] + g_y_0_yzzz_zzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 420);

            auto g_y_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 421);

            auto g_y_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 422);

            auto g_y_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 423);

            auto g_y_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 424);

            auto g_y_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 425);

            auto g_y_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 426);

            auto g_y_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 427);

            auto g_y_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 428);

            auto g_y_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 429);

            auto g_y_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 430);

            auto g_y_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 431);

            auto g_y_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 432);

            auto g_y_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 433);

            auto g_y_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 434);

            auto g_y_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 435);

            auto g_y_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 436);

            auto g_y_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 437);

            auto g_y_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 438);

            auto g_y_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 439);

            auto g_y_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 441 * acomps * bcomps + 440);

            #pragma omp simd aligned(cd_z, g_y_0_zzzz_xxxxx, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxy, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxyy, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxyyy, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xyyyy, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_yyyyy, g_y_0_zzzz_yyyyyz, g_y_0_zzzz_yyyyz, g_y_0_zzzz_yyyyzz, g_y_0_zzzz_yyyzz, g_y_0_zzzz_yyyzzz, g_y_0_zzzz_yyzzz, g_y_0_zzzz_yyzzzz, g_y_0_zzzz_yzzzz, g_y_0_zzzz_yzzzzz, g_y_0_zzzz_zzzzz, g_y_0_zzzz_zzzzzz, g_y_0_zzzzz_xxxxx, g_y_0_zzzzz_xxxxy, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxyy, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxyyy, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xyyyy, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_yyyyy, g_y_0_zzzzz_yyyyz, g_y_0_zzzzz_yyyzz, g_y_0_zzzzz_yyzzz, g_y_0_zzzzz_yzzzz, g_y_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xxxxx[k] = -g_y_0_zzzz_xxxxx[k] * cd_z[k] + g_y_0_zzzz_xxxxxz[k];

                g_y_0_zzzzz_xxxxy[k] = -g_y_0_zzzz_xxxxy[k] * cd_z[k] + g_y_0_zzzz_xxxxyz[k];

                g_y_0_zzzzz_xxxxz[k] = -g_y_0_zzzz_xxxxz[k] * cd_z[k] + g_y_0_zzzz_xxxxzz[k];

                g_y_0_zzzzz_xxxyy[k] = -g_y_0_zzzz_xxxyy[k] * cd_z[k] + g_y_0_zzzz_xxxyyz[k];

                g_y_0_zzzzz_xxxyz[k] = -g_y_0_zzzz_xxxyz[k] * cd_z[k] + g_y_0_zzzz_xxxyzz[k];

                g_y_0_zzzzz_xxxzz[k] = -g_y_0_zzzz_xxxzz[k] * cd_z[k] + g_y_0_zzzz_xxxzzz[k];

                g_y_0_zzzzz_xxyyy[k] = -g_y_0_zzzz_xxyyy[k] * cd_z[k] + g_y_0_zzzz_xxyyyz[k];

                g_y_0_zzzzz_xxyyz[k] = -g_y_0_zzzz_xxyyz[k] * cd_z[k] + g_y_0_zzzz_xxyyzz[k];

                g_y_0_zzzzz_xxyzz[k] = -g_y_0_zzzz_xxyzz[k] * cd_z[k] + g_y_0_zzzz_xxyzzz[k];

                g_y_0_zzzzz_xxzzz[k] = -g_y_0_zzzz_xxzzz[k] * cd_z[k] + g_y_0_zzzz_xxzzzz[k];

                g_y_0_zzzzz_xyyyy[k] = -g_y_0_zzzz_xyyyy[k] * cd_z[k] + g_y_0_zzzz_xyyyyz[k];

                g_y_0_zzzzz_xyyyz[k] = -g_y_0_zzzz_xyyyz[k] * cd_z[k] + g_y_0_zzzz_xyyyzz[k];

                g_y_0_zzzzz_xyyzz[k] = -g_y_0_zzzz_xyyzz[k] * cd_z[k] + g_y_0_zzzz_xyyzzz[k];

                g_y_0_zzzzz_xyzzz[k] = -g_y_0_zzzz_xyzzz[k] * cd_z[k] + g_y_0_zzzz_xyzzzz[k];

                g_y_0_zzzzz_xzzzz[k] = -g_y_0_zzzz_xzzzz[k] * cd_z[k] + g_y_0_zzzz_xzzzzz[k];

                g_y_0_zzzzz_yyyyy[k] = -g_y_0_zzzz_yyyyy[k] * cd_z[k] + g_y_0_zzzz_yyyyyz[k];

                g_y_0_zzzzz_yyyyz[k] = -g_y_0_zzzz_yyyyz[k] * cd_z[k] + g_y_0_zzzz_yyyyzz[k];

                g_y_0_zzzzz_yyyzz[k] = -g_y_0_zzzz_yyyzz[k] * cd_z[k] + g_y_0_zzzz_yyyzzz[k];

                g_y_0_zzzzz_yyzzz[k] = -g_y_0_zzzz_yyzzz[k] * cd_z[k] + g_y_0_zzzz_yyzzzz[k];

                g_y_0_zzzzz_yzzzz[k] = -g_y_0_zzzz_yzzzz[k] * cd_z[k] + g_y_0_zzzz_yzzzzz[k];

                g_y_0_zzzzz_zzzzz[k] = -g_y_0_zzzz_zzzzz[k] * cd_z[k] + g_y_0_zzzz_zzzzzz[k];
            }
            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 5);

            auto g_z_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 6);

            auto g_z_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 7);

            auto g_z_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 8);

            auto g_z_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 9);

            auto g_z_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 10);

            auto g_z_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 11);

            auto g_z_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 12);

            auto g_z_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 13);

            auto g_z_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 14);

            auto g_z_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 15);

            auto g_z_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 16);

            auto g_z_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 17);

            auto g_z_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 18);

            auto g_z_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 19);

            auto g_z_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_z_0_xxxx_xxxxx, g_z_0_xxxx_xxxxxx, g_z_0_xxxx_xxxxxy, g_z_0_xxxx_xxxxxz, g_z_0_xxxx_xxxxy, g_z_0_xxxx_xxxxyy, g_z_0_xxxx_xxxxyz, g_z_0_xxxx_xxxxz, g_z_0_xxxx_xxxxzz, g_z_0_xxxx_xxxyy, g_z_0_xxxx_xxxyyy, g_z_0_xxxx_xxxyyz, g_z_0_xxxx_xxxyz, g_z_0_xxxx_xxxyzz, g_z_0_xxxx_xxxzz, g_z_0_xxxx_xxxzzz, g_z_0_xxxx_xxyyy, g_z_0_xxxx_xxyyyy, g_z_0_xxxx_xxyyyz, g_z_0_xxxx_xxyyz, g_z_0_xxxx_xxyyzz, g_z_0_xxxx_xxyzz, g_z_0_xxxx_xxyzzz, g_z_0_xxxx_xxzzz, g_z_0_xxxx_xxzzzz, g_z_0_xxxx_xyyyy, g_z_0_xxxx_xyyyyy, g_z_0_xxxx_xyyyyz, g_z_0_xxxx_xyyyz, g_z_0_xxxx_xyyyzz, g_z_0_xxxx_xyyzz, g_z_0_xxxx_xyyzzz, g_z_0_xxxx_xyzzz, g_z_0_xxxx_xyzzzz, g_z_0_xxxx_xzzzz, g_z_0_xxxx_xzzzzz, g_z_0_xxxx_yyyyy, g_z_0_xxxx_yyyyz, g_z_0_xxxx_yyyzz, g_z_0_xxxx_yyzzz, g_z_0_xxxx_yzzzz, g_z_0_xxxx_zzzzz, g_z_0_xxxxx_xxxxx, g_z_0_xxxxx_xxxxy, g_z_0_xxxxx_xxxxz, g_z_0_xxxxx_xxxyy, g_z_0_xxxxx_xxxyz, g_z_0_xxxxx_xxxzz, g_z_0_xxxxx_xxyyy, g_z_0_xxxxx_xxyyz, g_z_0_xxxxx_xxyzz, g_z_0_xxxxx_xxzzz, g_z_0_xxxxx_xyyyy, g_z_0_xxxxx_xyyyz, g_z_0_xxxxx_xyyzz, g_z_0_xxxxx_xyzzz, g_z_0_xxxxx_xzzzz, g_z_0_xxxxx_yyyyy, g_z_0_xxxxx_yyyyz, g_z_0_xxxxx_yyyzz, g_z_0_xxxxx_yyzzz, g_z_0_xxxxx_yzzzz, g_z_0_xxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xxxxx[k] = -g_z_0_xxxx_xxxxx[k] * cd_x[k] + g_z_0_xxxx_xxxxxx[k];

                g_z_0_xxxxx_xxxxy[k] = -g_z_0_xxxx_xxxxy[k] * cd_x[k] + g_z_0_xxxx_xxxxxy[k];

                g_z_0_xxxxx_xxxxz[k] = -g_z_0_xxxx_xxxxz[k] * cd_x[k] + g_z_0_xxxx_xxxxxz[k];

                g_z_0_xxxxx_xxxyy[k] = -g_z_0_xxxx_xxxyy[k] * cd_x[k] + g_z_0_xxxx_xxxxyy[k];

                g_z_0_xxxxx_xxxyz[k] = -g_z_0_xxxx_xxxyz[k] * cd_x[k] + g_z_0_xxxx_xxxxyz[k];

                g_z_0_xxxxx_xxxzz[k] = -g_z_0_xxxx_xxxzz[k] * cd_x[k] + g_z_0_xxxx_xxxxzz[k];

                g_z_0_xxxxx_xxyyy[k] = -g_z_0_xxxx_xxyyy[k] * cd_x[k] + g_z_0_xxxx_xxxyyy[k];

                g_z_0_xxxxx_xxyyz[k] = -g_z_0_xxxx_xxyyz[k] * cd_x[k] + g_z_0_xxxx_xxxyyz[k];

                g_z_0_xxxxx_xxyzz[k] = -g_z_0_xxxx_xxyzz[k] * cd_x[k] + g_z_0_xxxx_xxxyzz[k];

                g_z_0_xxxxx_xxzzz[k] = -g_z_0_xxxx_xxzzz[k] * cd_x[k] + g_z_0_xxxx_xxxzzz[k];

                g_z_0_xxxxx_xyyyy[k] = -g_z_0_xxxx_xyyyy[k] * cd_x[k] + g_z_0_xxxx_xxyyyy[k];

                g_z_0_xxxxx_xyyyz[k] = -g_z_0_xxxx_xyyyz[k] * cd_x[k] + g_z_0_xxxx_xxyyyz[k];

                g_z_0_xxxxx_xyyzz[k] = -g_z_0_xxxx_xyyzz[k] * cd_x[k] + g_z_0_xxxx_xxyyzz[k];

                g_z_0_xxxxx_xyzzz[k] = -g_z_0_xxxx_xyzzz[k] * cd_x[k] + g_z_0_xxxx_xxyzzz[k];

                g_z_0_xxxxx_xzzzz[k] = -g_z_0_xxxx_xzzzz[k] * cd_x[k] + g_z_0_xxxx_xxzzzz[k];

                g_z_0_xxxxx_yyyyy[k] = -g_z_0_xxxx_yyyyy[k] * cd_x[k] + g_z_0_xxxx_xyyyyy[k];

                g_z_0_xxxxx_yyyyz[k] = -g_z_0_xxxx_yyyyz[k] * cd_x[k] + g_z_0_xxxx_xyyyyz[k];

                g_z_0_xxxxx_yyyzz[k] = -g_z_0_xxxx_yyyzz[k] * cd_x[k] + g_z_0_xxxx_xyyyzz[k];

                g_z_0_xxxxx_yyzzz[k] = -g_z_0_xxxx_yyzzz[k] * cd_x[k] + g_z_0_xxxx_xyyzzz[k];

                g_z_0_xxxxx_yzzzz[k] = -g_z_0_xxxx_yzzzz[k] * cd_x[k] + g_z_0_xxxx_xyzzzz[k];

                g_z_0_xxxxx_zzzzz[k] = -g_z_0_xxxx_zzzzz[k] * cd_x[k] + g_z_0_xxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 21);

            auto g_z_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 22);

            auto g_z_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 23);

            auto g_z_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 24);

            auto g_z_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 25);

            auto g_z_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 26);

            auto g_z_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 27);

            auto g_z_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 28);

            auto g_z_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 29);

            auto g_z_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 30);

            auto g_z_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 31);

            auto g_z_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 32);

            auto g_z_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 33);

            auto g_z_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 34);

            auto g_z_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 35);

            auto g_z_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 36);

            auto g_z_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 37);

            auto g_z_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 38);

            auto g_z_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 39);

            auto g_z_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 40);

            auto g_z_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxy_xxxxx, g_z_0_xxxxy_xxxxy, g_z_0_xxxxy_xxxxz, g_z_0_xxxxy_xxxyy, g_z_0_xxxxy_xxxyz, g_z_0_xxxxy_xxxzz, g_z_0_xxxxy_xxyyy, g_z_0_xxxxy_xxyyz, g_z_0_xxxxy_xxyzz, g_z_0_xxxxy_xxzzz, g_z_0_xxxxy_xyyyy, g_z_0_xxxxy_xyyyz, g_z_0_xxxxy_xyyzz, g_z_0_xxxxy_xyzzz, g_z_0_xxxxy_xzzzz, g_z_0_xxxxy_yyyyy, g_z_0_xxxxy_yyyyz, g_z_0_xxxxy_yyyzz, g_z_0_xxxxy_yyzzz, g_z_0_xxxxy_yzzzz, g_z_0_xxxxy_zzzzz, g_z_0_xxxy_xxxxx, g_z_0_xxxy_xxxxxx, g_z_0_xxxy_xxxxxy, g_z_0_xxxy_xxxxxz, g_z_0_xxxy_xxxxy, g_z_0_xxxy_xxxxyy, g_z_0_xxxy_xxxxyz, g_z_0_xxxy_xxxxz, g_z_0_xxxy_xxxxzz, g_z_0_xxxy_xxxyy, g_z_0_xxxy_xxxyyy, g_z_0_xxxy_xxxyyz, g_z_0_xxxy_xxxyz, g_z_0_xxxy_xxxyzz, g_z_0_xxxy_xxxzz, g_z_0_xxxy_xxxzzz, g_z_0_xxxy_xxyyy, g_z_0_xxxy_xxyyyy, g_z_0_xxxy_xxyyyz, g_z_0_xxxy_xxyyz, g_z_0_xxxy_xxyyzz, g_z_0_xxxy_xxyzz, g_z_0_xxxy_xxyzzz, g_z_0_xxxy_xxzzz, g_z_0_xxxy_xxzzzz, g_z_0_xxxy_xyyyy, g_z_0_xxxy_xyyyyy, g_z_0_xxxy_xyyyyz, g_z_0_xxxy_xyyyz, g_z_0_xxxy_xyyyzz, g_z_0_xxxy_xyyzz, g_z_0_xxxy_xyyzzz, g_z_0_xxxy_xyzzz, g_z_0_xxxy_xyzzzz, g_z_0_xxxy_xzzzz, g_z_0_xxxy_xzzzzz, g_z_0_xxxy_yyyyy, g_z_0_xxxy_yyyyz, g_z_0_xxxy_yyyzz, g_z_0_xxxy_yyzzz, g_z_0_xxxy_yzzzz, g_z_0_xxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xxxxx[k] = -g_z_0_xxxy_xxxxx[k] * cd_x[k] + g_z_0_xxxy_xxxxxx[k];

                g_z_0_xxxxy_xxxxy[k] = -g_z_0_xxxy_xxxxy[k] * cd_x[k] + g_z_0_xxxy_xxxxxy[k];

                g_z_0_xxxxy_xxxxz[k] = -g_z_0_xxxy_xxxxz[k] * cd_x[k] + g_z_0_xxxy_xxxxxz[k];

                g_z_0_xxxxy_xxxyy[k] = -g_z_0_xxxy_xxxyy[k] * cd_x[k] + g_z_0_xxxy_xxxxyy[k];

                g_z_0_xxxxy_xxxyz[k] = -g_z_0_xxxy_xxxyz[k] * cd_x[k] + g_z_0_xxxy_xxxxyz[k];

                g_z_0_xxxxy_xxxzz[k] = -g_z_0_xxxy_xxxzz[k] * cd_x[k] + g_z_0_xxxy_xxxxzz[k];

                g_z_0_xxxxy_xxyyy[k] = -g_z_0_xxxy_xxyyy[k] * cd_x[k] + g_z_0_xxxy_xxxyyy[k];

                g_z_0_xxxxy_xxyyz[k] = -g_z_0_xxxy_xxyyz[k] * cd_x[k] + g_z_0_xxxy_xxxyyz[k];

                g_z_0_xxxxy_xxyzz[k] = -g_z_0_xxxy_xxyzz[k] * cd_x[k] + g_z_0_xxxy_xxxyzz[k];

                g_z_0_xxxxy_xxzzz[k] = -g_z_0_xxxy_xxzzz[k] * cd_x[k] + g_z_0_xxxy_xxxzzz[k];

                g_z_0_xxxxy_xyyyy[k] = -g_z_0_xxxy_xyyyy[k] * cd_x[k] + g_z_0_xxxy_xxyyyy[k];

                g_z_0_xxxxy_xyyyz[k] = -g_z_0_xxxy_xyyyz[k] * cd_x[k] + g_z_0_xxxy_xxyyyz[k];

                g_z_0_xxxxy_xyyzz[k] = -g_z_0_xxxy_xyyzz[k] * cd_x[k] + g_z_0_xxxy_xxyyzz[k];

                g_z_0_xxxxy_xyzzz[k] = -g_z_0_xxxy_xyzzz[k] * cd_x[k] + g_z_0_xxxy_xxyzzz[k];

                g_z_0_xxxxy_xzzzz[k] = -g_z_0_xxxy_xzzzz[k] * cd_x[k] + g_z_0_xxxy_xxzzzz[k];

                g_z_0_xxxxy_yyyyy[k] = -g_z_0_xxxy_yyyyy[k] * cd_x[k] + g_z_0_xxxy_xyyyyy[k];

                g_z_0_xxxxy_yyyyz[k] = -g_z_0_xxxy_yyyyz[k] * cd_x[k] + g_z_0_xxxy_xyyyyz[k];

                g_z_0_xxxxy_yyyzz[k] = -g_z_0_xxxy_yyyzz[k] * cd_x[k] + g_z_0_xxxy_xyyyzz[k];

                g_z_0_xxxxy_yyzzz[k] = -g_z_0_xxxy_yyzzz[k] * cd_x[k] + g_z_0_xxxy_xyyzzz[k];

                g_z_0_xxxxy_yzzzz[k] = -g_z_0_xxxy_yzzzz[k] * cd_x[k] + g_z_0_xxxy_xyzzzz[k];

                g_z_0_xxxxy_zzzzz[k] = -g_z_0_xxxy_zzzzz[k] * cd_x[k] + g_z_0_xxxy_xzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 42);

            auto g_z_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 43);

            auto g_z_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 44);

            auto g_z_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 45);

            auto g_z_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 46);

            auto g_z_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 47);

            auto g_z_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 48);

            auto g_z_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 49);

            auto g_z_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 50);

            auto g_z_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 51);

            auto g_z_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 52);

            auto g_z_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 53);

            auto g_z_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 54);

            auto g_z_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 55);

            auto g_z_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 56);

            auto g_z_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 57);

            auto g_z_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 58);

            auto g_z_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 59);

            auto g_z_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 60);

            auto g_z_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 61);

            auto g_z_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxz_xxxxx, g_z_0_xxxxz_xxxxy, g_z_0_xxxxz_xxxxz, g_z_0_xxxxz_xxxyy, g_z_0_xxxxz_xxxyz, g_z_0_xxxxz_xxxzz, g_z_0_xxxxz_xxyyy, g_z_0_xxxxz_xxyyz, g_z_0_xxxxz_xxyzz, g_z_0_xxxxz_xxzzz, g_z_0_xxxxz_xyyyy, g_z_0_xxxxz_xyyyz, g_z_0_xxxxz_xyyzz, g_z_0_xxxxz_xyzzz, g_z_0_xxxxz_xzzzz, g_z_0_xxxxz_yyyyy, g_z_0_xxxxz_yyyyz, g_z_0_xxxxz_yyyzz, g_z_0_xxxxz_yyzzz, g_z_0_xxxxz_yzzzz, g_z_0_xxxxz_zzzzz, g_z_0_xxxz_xxxxx, g_z_0_xxxz_xxxxxx, g_z_0_xxxz_xxxxxy, g_z_0_xxxz_xxxxxz, g_z_0_xxxz_xxxxy, g_z_0_xxxz_xxxxyy, g_z_0_xxxz_xxxxyz, g_z_0_xxxz_xxxxz, g_z_0_xxxz_xxxxzz, g_z_0_xxxz_xxxyy, g_z_0_xxxz_xxxyyy, g_z_0_xxxz_xxxyyz, g_z_0_xxxz_xxxyz, g_z_0_xxxz_xxxyzz, g_z_0_xxxz_xxxzz, g_z_0_xxxz_xxxzzz, g_z_0_xxxz_xxyyy, g_z_0_xxxz_xxyyyy, g_z_0_xxxz_xxyyyz, g_z_0_xxxz_xxyyz, g_z_0_xxxz_xxyyzz, g_z_0_xxxz_xxyzz, g_z_0_xxxz_xxyzzz, g_z_0_xxxz_xxzzz, g_z_0_xxxz_xxzzzz, g_z_0_xxxz_xyyyy, g_z_0_xxxz_xyyyyy, g_z_0_xxxz_xyyyyz, g_z_0_xxxz_xyyyz, g_z_0_xxxz_xyyyzz, g_z_0_xxxz_xyyzz, g_z_0_xxxz_xyyzzz, g_z_0_xxxz_xyzzz, g_z_0_xxxz_xyzzzz, g_z_0_xxxz_xzzzz, g_z_0_xxxz_xzzzzz, g_z_0_xxxz_yyyyy, g_z_0_xxxz_yyyyz, g_z_0_xxxz_yyyzz, g_z_0_xxxz_yyzzz, g_z_0_xxxz_yzzzz, g_z_0_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xxxxx[k] = -g_z_0_xxxz_xxxxx[k] * cd_x[k] + g_z_0_xxxz_xxxxxx[k];

                g_z_0_xxxxz_xxxxy[k] = -g_z_0_xxxz_xxxxy[k] * cd_x[k] + g_z_0_xxxz_xxxxxy[k];

                g_z_0_xxxxz_xxxxz[k] = -g_z_0_xxxz_xxxxz[k] * cd_x[k] + g_z_0_xxxz_xxxxxz[k];

                g_z_0_xxxxz_xxxyy[k] = -g_z_0_xxxz_xxxyy[k] * cd_x[k] + g_z_0_xxxz_xxxxyy[k];

                g_z_0_xxxxz_xxxyz[k] = -g_z_0_xxxz_xxxyz[k] * cd_x[k] + g_z_0_xxxz_xxxxyz[k];

                g_z_0_xxxxz_xxxzz[k] = -g_z_0_xxxz_xxxzz[k] * cd_x[k] + g_z_0_xxxz_xxxxzz[k];

                g_z_0_xxxxz_xxyyy[k] = -g_z_0_xxxz_xxyyy[k] * cd_x[k] + g_z_0_xxxz_xxxyyy[k];

                g_z_0_xxxxz_xxyyz[k] = -g_z_0_xxxz_xxyyz[k] * cd_x[k] + g_z_0_xxxz_xxxyyz[k];

                g_z_0_xxxxz_xxyzz[k] = -g_z_0_xxxz_xxyzz[k] * cd_x[k] + g_z_0_xxxz_xxxyzz[k];

                g_z_0_xxxxz_xxzzz[k] = -g_z_0_xxxz_xxzzz[k] * cd_x[k] + g_z_0_xxxz_xxxzzz[k];

                g_z_0_xxxxz_xyyyy[k] = -g_z_0_xxxz_xyyyy[k] * cd_x[k] + g_z_0_xxxz_xxyyyy[k];

                g_z_0_xxxxz_xyyyz[k] = -g_z_0_xxxz_xyyyz[k] * cd_x[k] + g_z_0_xxxz_xxyyyz[k];

                g_z_0_xxxxz_xyyzz[k] = -g_z_0_xxxz_xyyzz[k] * cd_x[k] + g_z_0_xxxz_xxyyzz[k];

                g_z_0_xxxxz_xyzzz[k] = -g_z_0_xxxz_xyzzz[k] * cd_x[k] + g_z_0_xxxz_xxyzzz[k];

                g_z_0_xxxxz_xzzzz[k] = -g_z_0_xxxz_xzzzz[k] * cd_x[k] + g_z_0_xxxz_xxzzzz[k];

                g_z_0_xxxxz_yyyyy[k] = -g_z_0_xxxz_yyyyy[k] * cd_x[k] + g_z_0_xxxz_xyyyyy[k];

                g_z_0_xxxxz_yyyyz[k] = -g_z_0_xxxz_yyyyz[k] * cd_x[k] + g_z_0_xxxz_xyyyyz[k];

                g_z_0_xxxxz_yyyzz[k] = -g_z_0_xxxz_yyyzz[k] * cd_x[k] + g_z_0_xxxz_xyyyzz[k];

                g_z_0_xxxxz_yyzzz[k] = -g_z_0_xxxz_yyzzz[k] * cd_x[k] + g_z_0_xxxz_xyyzzz[k];

                g_z_0_xxxxz_yzzzz[k] = -g_z_0_xxxz_yzzzz[k] * cd_x[k] + g_z_0_xxxz_xyzzzz[k];

                g_z_0_xxxxz_zzzzz[k] = -g_z_0_xxxz_zzzzz[k] * cd_x[k] + g_z_0_xxxz_xzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 63);

            auto g_z_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 64);

            auto g_z_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 65);

            auto g_z_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 66);

            auto g_z_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 67);

            auto g_z_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 68);

            auto g_z_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 69);

            auto g_z_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 70);

            auto g_z_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 71);

            auto g_z_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 72);

            auto g_z_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 73);

            auto g_z_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 74);

            auto g_z_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 75);

            auto g_z_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 76);

            auto g_z_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 77);

            auto g_z_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 78);

            auto g_z_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 79);

            auto g_z_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 80);

            auto g_z_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 81);

            auto g_z_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 82);

            auto g_z_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyy_xxxxx, g_z_0_xxxyy_xxxxy, g_z_0_xxxyy_xxxxz, g_z_0_xxxyy_xxxyy, g_z_0_xxxyy_xxxyz, g_z_0_xxxyy_xxxzz, g_z_0_xxxyy_xxyyy, g_z_0_xxxyy_xxyyz, g_z_0_xxxyy_xxyzz, g_z_0_xxxyy_xxzzz, g_z_0_xxxyy_xyyyy, g_z_0_xxxyy_xyyyz, g_z_0_xxxyy_xyyzz, g_z_0_xxxyy_xyzzz, g_z_0_xxxyy_xzzzz, g_z_0_xxxyy_yyyyy, g_z_0_xxxyy_yyyyz, g_z_0_xxxyy_yyyzz, g_z_0_xxxyy_yyzzz, g_z_0_xxxyy_yzzzz, g_z_0_xxxyy_zzzzz, g_z_0_xxyy_xxxxx, g_z_0_xxyy_xxxxxx, g_z_0_xxyy_xxxxxy, g_z_0_xxyy_xxxxxz, g_z_0_xxyy_xxxxy, g_z_0_xxyy_xxxxyy, g_z_0_xxyy_xxxxyz, g_z_0_xxyy_xxxxz, g_z_0_xxyy_xxxxzz, g_z_0_xxyy_xxxyy, g_z_0_xxyy_xxxyyy, g_z_0_xxyy_xxxyyz, g_z_0_xxyy_xxxyz, g_z_0_xxyy_xxxyzz, g_z_0_xxyy_xxxzz, g_z_0_xxyy_xxxzzz, g_z_0_xxyy_xxyyy, g_z_0_xxyy_xxyyyy, g_z_0_xxyy_xxyyyz, g_z_0_xxyy_xxyyz, g_z_0_xxyy_xxyyzz, g_z_0_xxyy_xxyzz, g_z_0_xxyy_xxyzzz, g_z_0_xxyy_xxzzz, g_z_0_xxyy_xxzzzz, g_z_0_xxyy_xyyyy, g_z_0_xxyy_xyyyyy, g_z_0_xxyy_xyyyyz, g_z_0_xxyy_xyyyz, g_z_0_xxyy_xyyyzz, g_z_0_xxyy_xyyzz, g_z_0_xxyy_xyyzzz, g_z_0_xxyy_xyzzz, g_z_0_xxyy_xyzzzz, g_z_0_xxyy_xzzzz, g_z_0_xxyy_xzzzzz, g_z_0_xxyy_yyyyy, g_z_0_xxyy_yyyyz, g_z_0_xxyy_yyyzz, g_z_0_xxyy_yyzzz, g_z_0_xxyy_yzzzz, g_z_0_xxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xxxxx[k] = -g_z_0_xxyy_xxxxx[k] * cd_x[k] + g_z_0_xxyy_xxxxxx[k];

                g_z_0_xxxyy_xxxxy[k] = -g_z_0_xxyy_xxxxy[k] * cd_x[k] + g_z_0_xxyy_xxxxxy[k];

                g_z_0_xxxyy_xxxxz[k] = -g_z_0_xxyy_xxxxz[k] * cd_x[k] + g_z_0_xxyy_xxxxxz[k];

                g_z_0_xxxyy_xxxyy[k] = -g_z_0_xxyy_xxxyy[k] * cd_x[k] + g_z_0_xxyy_xxxxyy[k];

                g_z_0_xxxyy_xxxyz[k] = -g_z_0_xxyy_xxxyz[k] * cd_x[k] + g_z_0_xxyy_xxxxyz[k];

                g_z_0_xxxyy_xxxzz[k] = -g_z_0_xxyy_xxxzz[k] * cd_x[k] + g_z_0_xxyy_xxxxzz[k];

                g_z_0_xxxyy_xxyyy[k] = -g_z_0_xxyy_xxyyy[k] * cd_x[k] + g_z_0_xxyy_xxxyyy[k];

                g_z_0_xxxyy_xxyyz[k] = -g_z_0_xxyy_xxyyz[k] * cd_x[k] + g_z_0_xxyy_xxxyyz[k];

                g_z_0_xxxyy_xxyzz[k] = -g_z_0_xxyy_xxyzz[k] * cd_x[k] + g_z_0_xxyy_xxxyzz[k];

                g_z_0_xxxyy_xxzzz[k] = -g_z_0_xxyy_xxzzz[k] * cd_x[k] + g_z_0_xxyy_xxxzzz[k];

                g_z_0_xxxyy_xyyyy[k] = -g_z_0_xxyy_xyyyy[k] * cd_x[k] + g_z_0_xxyy_xxyyyy[k];

                g_z_0_xxxyy_xyyyz[k] = -g_z_0_xxyy_xyyyz[k] * cd_x[k] + g_z_0_xxyy_xxyyyz[k];

                g_z_0_xxxyy_xyyzz[k] = -g_z_0_xxyy_xyyzz[k] * cd_x[k] + g_z_0_xxyy_xxyyzz[k];

                g_z_0_xxxyy_xyzzz[k] = -g_z_0_xxyy_xyzzz[k] * cd_x[k] + g_z_0_xxyy_xxyzzz[k];

                g_z_0_xxxyy_xzzzz[k] = -g_z_0_xxyy_xzzzz[k] * cd_x[k] + g_z_0_xxyy_xxzzzz[k];

                g_z_0_xxxyy_yyyyy[k] = -g_z_0_xxyy_yyyyy[k] * cd_x[k] + g_z_0_xxyy_xyyyyy[k];

                g_z_0_xxxyy_yyyyz[k] = -g_z_0_xxyy_yyyyz[k] * cd_x[k] + g_z_0_xxyy_xyyyyz[k];

                g_z_0_xxxyy_yyyzz[k] = -g_z_0_xxyy_yyyzz[k] * cd_x[k] + g_z_0_xxyy_xyyyzz[k];

                g_z_0_xxxyy_yyzzz[k] = -g_z_0_xxyy_yyzzz[k] * cd_x[k] + g_z_0_xxyy_xyyzzz[k];

                g_z_0_xxxyy_yzzzz[k] = -g_z_0_xxyy_yzzzz[k] * cd_x[k] + g_z_0_xxyy_xyzzzz[k];

                g_z_0_xxxyy_zzzzz[k] = -g_z_0_xxyy_zzzzz[k] * cd_x[k] + g_z_0_xxyy_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 84);

            auto g_z_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 85);

            auto g_z_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 86);

            auto g_z_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 87);

            auto g_z_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 88);

            auto g_z_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 89);

            auto g_z_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 90);

            auto g_z_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 91);

            auto g_z_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 92);

            auto g_z_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 93);

            auto g_z_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 94);

            auto g_z_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 95);

            auto g_z_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 96);

            auto g_z_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 97);

            auto g_z_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 98);

            auto g_z_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 99);

            auto g_z_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 100);

            auto g_z_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 101);

            auto g_z_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 102);

            auto g_z_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 103);

            auto g_z_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyz_xxxxx, g_z_0_xxxyz_xxxxy, g_z_0_xxxyz_xxxxz, g_z_0_xxxyz_xxxyy, g_z_0_xxxyz_xxxyz, g_z_0_xxxyz_xxxzz, g_z_0_xxxyz_xxyyy, g_z_0_xxxyz_xxyyz, g_z_0_xxxyz_xxyzz, g_z_0_xxxyz_xxzzz, g_z_0_xxxyz_xyyyy, g_z_0_xxxyz_xyyyz, g_z_0_xxxyz_xyyzz, g_z_0_xxxyz_xyzzz, g_z_0_xxxyz_xzzzz, g_z_0_xxxyz_yyyyy, g_z_0_xxxyz_yyyyz, g_z_0_xxxyz_yyyzz, g_z_0_xxxyz_yyzzz, g_z_0_xxxyz_yzzzz, g_z_0_xxxyz_zzzzz, g_z_0_xxyz_xxxxx, g_z_0_xxyz_xxxxxx, g_z_0_xxyz_xxxxxy, g_z_0_xxyz_xxxxxz, g_z_0_xxyz_xxxxy, g_z_0_xxyz_xxxxyy, g_z_0_xxyz_xxxxyz, g_z_0_xxyz_xxxxz, g_z_0_xxyz_xxxxzz, g_z_0_xxyz_xxxyy, g_z_0_xxyz_xxxyyy, g_z_0_xxyz_xxxyyz, g_z_0_xxyz_xxxyz, g_z_0_xxyz_xxxyzz, g_z_0_xxyz_xxxzz, g_z_0_xxyz_xxxzzz, g_z_0_xxyz_xxyyy, g_z_0_xxyz_xxyyyy, g_z_0_xxyz_xxyyyz, g_z_0_xxyz_xxyyz, g_z_0_xxyz_xxyyzz, g_z_0_xxyz_xxyzz, g_z_0_xxyz_xxyzzz, g_z_0_xxyz_xxzzz, g_z_0_xxyz_xxzzzz, g_z_0_xxyz_xyyyy, g_z_0_xxyz_xyyyyy, g_z_0_xxyz_xyyyyz, g_z_0_xxyz_xyyyz, g_z_0_xxyz_xyyyzz, g_z_0_xxyz_xyyzz, g_z_0_xxyz_xyyzzz, g_z_0_xxyz_xyzzz, g_z_0_xxyz_xyzzzz, g_z_0_xxyz_xzzzz, g_z_0_xxyz_xzzzzz, g_z_0_xxyz_yyyyy, g_z_0_xxyz_yyyyz, g_z_0_xxyz_yyyzz, g_z_0_xxyz_yyzzz, g_z_0_xxyz_yzzzz, g_z_0_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xxxxx[k] = -g_z_0_xxyz_xxxxx[k] * cd_x[k] + g_z_0_xxyz_xxxxxx[k];

                g_z_0_xxxyz_xxxxy[k] = -g_z_0_xxyz_xxxxy[k] * cd_x[k] + g_z_0_xxyz_xxxxxy[k];

                g_z_0_xxxyz_xxxxz[k] = -g_z_0_xxyz_xxxxz[k] * cd_x[k] + g_z_0_xxyz_xxxxxz[k];

                g_z_0_xxxyz_xxxyy[k] = -g_z_0_xxyz_xxxyy[k] * cd_x[k] + g_z_0_xxyz_xxxxyy[k];

                g_z_0_xxxyz_xxxyz[k] = -g_z_0_xxyz_xxxyz[k] * cd_x[k] + g_z_0_xxyz_xxxxyz[k];

                g_z_0_xxxyz_xxxzz[k] = -g_z_0_xxyz_xxxzz[k] * cd_x[k] + g_z_0_xxyz_xxxxzz[k];

                g_z_0_xxxyz_xxyyy[k] = -g_z_0_xxyz_xxyyy[k] * cd_x[k] + g_z_0_xxyz_xxxyyy[k];

                g_z_0_xxxyz_xxyyz[k] = -g_z_0_xxyz_xxyyz[k] * cd_x[k] + g_z_0_xxyz_xxxyyz[k];

                g_z_0_xxxyz_xxyzz[k] = -g_z_0_xxyz_xxyzz[k] * cd_x[k] + g_z_0_xxyz_xxxyzz[k];

                g_z_0_xxxyz_xxzzz[k] = -g_z_0_xxyz_xxzzz[k] * cd_x[k] + g_z_0_xxyz_xxxzzz[k];

                g_z_0_xxxyz_xyyyy[k] = -g_z_0_xxyz_xyyyy[k] * cd_x[k] + g_z_0_xxyz_xxyyyy[k];

                g_z_0_xxxyz_xyyyz[k] = -g_z_0_xxyz_xyyyz[k] * cd_x[k] + g_z_0_xxyz_xxyyyz[k];

                g_z_0_xxxyz_xyyzz[k] = -g_z_0_xxyz_xyyzz[k] * cd_x[k] + g_z_0_xxyz_xxyyzz[k];

                g_z_0_xxxyz_xyzzz[k] = -g_z_0_xxyz_xyzzz[k] * cd_x[k] + g_z_0_xxyz_xxyzzz[k];

                g_z_0_xxxyz_xzzzz[k] = -g_z_0_xxyz_xzzzz[k] * cd_x[k] + g_z_0_xxyz_xxzzzz[k];

                g_z_0_xxxyz_yyyyy[k] = -g_z_0_xxyz_yyyyy[k] * cd_x[k] + g_z_0_xxyz_xyyyyy[k];

                g_z_0_xxxyz_yyyyz[k] = -g_z_0_xxyz_yyyyz[k] * cd_x[k] + g_z_0_xxyz_xyyyyz[k];

                g_z_0_xxxyz_yyyzz[k] = -g_z_0_xxyz_yyyzz[k] * cd_x[k] + g_z_0_xxyz_xyyyzz[k];

                g_z_0_xxxyz_yyzzz[k] = -g_z_0_xxyz_yyzzz[k] * cd_x[k] + g_z_0_xxyz_xyyzzz[k];

                g_z_0_xxxyz_yzzzz[k] = -g_z_0_xxyz_yzzzz[k] * cd_x[k] + g_z_0_xxyz_xyzzzz[k];

                g_z_0_xxxyz_zzzzz[k] = -g_z_0_xxyz_zzzzz[k] * cd_x[k] + g_z_0_xxyz_xzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 105);

            auto g_z_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 106);

            auto g_z_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 107);

            auto g_z_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 108);

            auto g_z_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 109);

            auto g_z_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 110);

            auto g_z_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 111);

            auto g_z_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 112);

            auto g_z_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 113);

            auto g_z_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 114);

            auto g_z_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 115);

            auto g_z_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 116);

            auto g_z_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 117);

            auto g_z_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 118);

            auto g_z_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 119);

            auto g_z_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 120);

            auto g_z_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 121);

            auto g_z_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 122);

            auto g_z_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 123);

            auto g_z_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 124);

            auto g_z_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzz_xxxxx, g_z_0_xxxzz_xxxxy, g_z_0_xxxzz_xxxxz, g_z_0_xxxzz_xxxyy, g_z_0_xxxzz_xxxyz, g_z_0_xxxzz_xxxzz, g_z_0_xxxzz_xxyyy, g_z_0_xxxzz_xxyyz, g_z_0_xxxzz_xxyzz, g_z_0_xxxzz_xxzzz, g_z_0_xxxzz_xyyyy, g_z_0_xxxzz_xyyyz, g_z_0_xxxzz_xyyzz, g_z_0_xxxzz_xyzzz, g_z_0_xxxzz_xzzzz, g_z_0_xxxzz_yyyyy, g_z_0_xxxzz_yyyyz, g_z_0_xxxzz_yyyzz, g_z_0_xxxzz_yyzzz, g_z_0_xxxzz_yzzzz, g_z_0_xxxzz_zzzzz, g_z_0_xxzz_xxxxx, g_z_0_xxzz_xxxxxx, g_z_0_xxzz_xxxxxy, g_z_0_xxzz_xxxxxz, g_z_0_xxzz_xxxxy, g_z_0_xxzz_xxxxyy, g_z_0_xxzz_xxxxyz, g_z_0_xxzz_xxxxz, g_z_0_xxzz_xxxxzz, g_z_0_xxzz_xxxyy, g_z_0_xxzz_xxxyyy, g_z_0_xxzz_xxxyyz, g_z_0_xxzz_xxxyz, g_z_0_xxzz_xxxyzz, g_z_0_xxzz_xxxzz, g_z_0_xxzz_xxxzzz, g_z_0_xxzz_xxyyy, g_z_0_xxzz_xxyyyy, g_z_0_xxzz_xxyyyz, g_z_0_xxzz_xxyyz, g_z_0_xxzz_xxyyzz, g_z_0_xxzz_xxyzz, g_z_0_xxzz_xxyzzz, g_z_0_xxzz_xxzzz, g_z_0_xxzz_xxzzzz, g_z_0_xxzz_xyyyy, g_z_0_xxzz_xyyyyy, g_z_0_xxzz_xyyyyz, g_z_0_xxzz_xyyyz, g_z_0_xxzz_xyyyzz, g_z_0_xxzz_xyyzz, g_z_0_xxzz_xyyzzz, g_z_0_xxzz_xyzzz, g_z_0_xxzz_xyzzzz, g_z_0_xxzz_xzzzz, g_z_0_xxzz_xzzzzz, g_z_0_xxzz_yyyyy, g_z_0_xxzz_yyyyz, g_z_0_xxzz_yyyzz, g_z_0_xxzz_yyzzz, g_z_0_xxzz_yzzzz, g_z_0_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xxxxx[k] = -g_z_0_xxzz_xxxxx[k] * cd_x[k] + g_z_0_xxzz_xxxxxx[k];

                g_z_0_xxxzz_xxxxy[k] = -g_z_0_xxzz_xxxxy[k] * cd_x[k] + g_z_0_xxzz_xxxxxy[k];

                g_z_0_xxxzz_xxxxz[k] = -g_z_0_xxzz_xxxxz[k] * cd_x[k] + g_z_0_xxzz_xxxxxz[k];

                g_z_0_xxxzz_xxxyy[k] = -g_z_0_xxzz_xxxyy[k] * cd_x[k] + g_z_0_xxzz_xxxxyy[k];

                g_z_0_xxxzz_xxxyz[k] = -g_z_0_xxzz_xxxyz[k] * cd_x[k] + g_z_0_xxzz_xxxxyz[k];

                g_z_0_xxxzz_xxxzz[k] = -g_z_0_xxzz_xxxzz[k] * cd_x[k] + g_z_0_xxzz_xxxxzz[k];

                g_z_0_xxxzz_xxyyy[k] = -g_z_0_xxzz_xxyyy[k] * cd_x[k] + g_z_0_xxzz_xxxyyy[k];

                g_z_0_xxxzz_xxyyz[k] = -g_z_0_xxzz_xxyyz[k] * cd_x[k] + g_z_0_xxzz_xxxyyz[k];

                g_z_0_xxxzz_xxyzz[k] = -g_z_0_xxzz_xxyzz[k] * cd_x[k] + g_z_0_xxzz_xxxyzz[k];

                g_z_0_xxxzz_xxzzz[k] = -g_z_0_xxzz_xxzzz[k] * cd_x[k] + g_z_0_xxzz_xxxzzz[k];

                g_z_0_xxxzz_xyyyy[k] = -g_z_0_xxzz_xyyyy[k] * cd_x[k] + g_z_0_xxzz_xxyyyy[k];

                g_z_0_xxxzz_xyyyz[k] = -g_z_0_xxzz_xyyyz[k] * cd_x[k] + g_z_0_xxzz_xxyyyz[k];

                g_z_0_xxxzz_xyyzz[k] = -g_z_0_xxzz_xyyzz[k] * cd_x[k] + g_z_0_xxzz_xxyyzz[k];

                g_z_0_xxxzz_xyzzz[k] = -g_z_0_xxzz_xyzzz[k] * cd_x[k] + g_z_0_xxzz_xxyzzz[k];

                g_z_0_xxxzz_xzzzz[k] = -g_z_0_xxzz_xzzzz[k] * cd_x[k] + g_z_0_xxzz_xxzzzz[k];

                g_z_0_xxxzz_yyyyy[k] = -g_z_0_xxzz_yyyyy[k] * cd_x[k] + g_z_0_xxzz_xyyyyy[k];

                g_z_0_xxxzz_yyyyz[k] = -g_z_0_xxzz_yyyyz[k] * cd_x[k] + g_z_0_xxzz_xyyyyz[k];

                g_z_0_xxxzz_yyyzz[k] = -g_z_0_xxzz_yyyzz[k] * cd_x[k] + g_z_0_xxzz_xyyyzz[k];

                g_z_0_xxxzz_yyzzz[k] = -g_z_0_xxzz_yyzzz[k] * cd_x[k] + g_z_0_xxzz_xyyzzz[k];

                g_z_0_xxxzz_yzzzz[k] = -g_z_0_xxzz_yzzzz[k] * cd_x[k] + g_z_0_xxzz_xyzzzz[k];

                g_z_0_xxxzz_zzzzz[k] = -g_z_0_xxzz_zzzzz[k] * cd_x[k] + g_z_0_xxzz_xzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 126);

            auto g_z_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 127);

            auto g_z_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 128);

            auto g_z_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 129);

            auto g_z_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 130);

            auto g_z_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 131);

            auto g_z_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 132);

            auto g_z_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 133);

            auto g_z_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 134);

            auto g_z_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 135);

            auto g_z_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 136);

            auto g_z_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 137);

            auto g_z_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 138);

            auto g_z_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 139);

            auto g_z_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 140);

            auto g_z_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 141);

            auto g_z_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 142);

            auto g_z_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 143);

            auto g_z_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 144);

            auto g_z_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 145);

            auto g_z_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyy_xxxxx, g_z_0_xxyyy_xxxxy, g_z_0_xxyyy_xxxxz, g_z_0_xxyyy_xxxyy, g_z_0_xxyyy_xxxyz, g_z_0_xxyyy_xxxzz, g_z_0_xxyyy_xxyyy, g_z_0_xxyyy_xxyyz, g_z_0_xxyyy_xxyzz, g_z_0_xxyyy_xxzzz, g_z_0_xxyyy_xyyyy, g_z_0_xxyyy_xyyyz, g_z_0_xxyyy_xyyzz, g_z_0_xxyyy_xyzzz, g_z_0_xxyyy_xzzzz, g_z_0_xxyyy_yyyyy, g_z_0_xxyyy_yyyyz, g_z_0_xxyyy_yyyzz, g_z_0_xxyyy_yyzzz, g_z_0_xxyyy_yzzzz, g_z_0_xxyyy_zzzzz, g_z_0_xyyy_xxxxx, g_z_0_xyyy_xxxxxx, g_z_0_xyyy_xxxxxy, g_z_0_xyyy_xxxxxz, g_z_0_xyyy_xxxxy, g_z_0_xyyy_xxxxyy, g_z_0_xyyy_xxxxyz, g_z_0_xyyy_xxxxz, g_z_0_xyyy_xxxxzz, g_z_0_xyyy_xxxyy, g_z_0_xyyy_xxxyyy, g_z_0_xyyy_xxxyyz, g_z_0_xyyy_xxxyz, g_z_0_xyyy_xxxyzz, g_z_0_xyyy_xxxzz, g_z_0_xyyy_xxxzzz, g_z_0_xyyy_xxyyy, g_z_0_xyyy_xxyyyy, g_z_0_xyyy_xxyyyz, g_z_0_xyyy_xxyyz, g_z_0_xyyy_xxyyzz, g_z_0_xyyy_xxyzz, g_z_0_xyyy_xxyzzz, g_z_0_xyyy_xxzzz, g_z_0_xyyy_xxzzzz, g_z_0_xyyy_xyyyy, g_z_0_xyyy_xyyyyy, g_z_0_xyyy_xyyyyz, g_z_0_xyyy_xyyyz, g_z_0_xyyy_xyyyzz, g_z_0_xyyy_xyyzz, g_z_0_xyyy_xyyzzz, g_z_0_xyyy_xyzzz, g_z_0_xyyy_xyzzzz, g_z_0_xyyy_xzzzz, g_z_0_xyyy_xzzzzz, g_z_0_xyyy_yyyyy, g_z_0_xyyy_yyyyz, g_z_0_xyyy_yyyzz, g_z_0_xyyy_yyzzz, g_z_0_xyyy_yzzzz, g_z_0_xyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xxxxx[k] = -g_z_0_xyyy_xxxxx[k] * cd_x[k] + g_z_0_xyyy_xxxxxx[k];

                g_z_0_xxyyy_xxxxy[k] = -g_z_0_xyyy_xxxxy[k] * cd_x[k] + g_z_0_xyyy_xxxxxy[k];

                g_z_0_xxyyy_xxxxz[k] = -g_z_0_xyyy_xxxxz[k] * cd_x[k] + g_z_0_xyyy_xxxxxz[k];

                g_z_0_xxyyy_xxxyy[k] = -g_z_0_xyyy_xxxyy[k] * cd_x[k] + g_z_0_xyyy_xxxxyy[k];

                g_z_0_xxyyy_xxxyz[k] = -g_z_0_xyyy_xxxyz[k] * cd_x[k] + g_z_0_xyyy_xxxxyz[k];

                g_z_0_xxyyy_xxxzz[k] = -g_z_0_xyyy_xxxzz[k] * cd_x[k] + g_z_0_xyyy_xxxxzz[k];

                g_z_0_xxyyy_xxyyy[k] = -g_z_0_xyyy_xxyyy[k] * cd_x[k] + g_z_0_xyyy_xxxyyy[k];

                g_z_0_xxyyy_xxyyz[k] = -g_z_0_xyyy_xxyyz[k] * cd_x[k] + g_z_0_xyyy_xxxyyz[k];

                g_z_0_xxyyy_xxyzz[k] = -g_z_0_xyyy_xxyzz[k] * cd_x[k] + g_z_0_xyyy_xxxyzz[k];

                g_z_0_xxyyy_xxzzz[k] = -g_z_0_xyyy_xxzzz[k] * cd_x[k] + g_z_0_xyyy_xxxzzz[k];

                g_z_0_xxyyy_xyyyy[k] = -g_z_0_xyyy_xyyyy[k] * cd_x[k] + g_z_0_xyyy_xxyyyy[k];

                g_z_0_xxyyy_xyyyz[k] = -g_z_0_xyyy_xyyyz[k] * cd_x[k] + g_z_0_xyyy_xxyyyz[k];

                g_z_0_xxyyy_xyyzz[k] = -g_z_0_xyyy_xyyzz[k] * cd_x[k] + g_z_0_xyyy_xxyyzz[k];

                g_z_0_xxyyy_xyzzz[k] = -g_z_0_xyyy_xyzzz[k] * cd_x[k] + g_z_0_xyyy_xxyzzz[k];

                g_z_0_xxyyy_xzzzz[k] = -g_z_0_xyyy_xzzzz[k] * cd_x[k] + g_z_0_xyyy_xxzzzz[k];

                g_z_0_xxyyy_yyyyy[k] = -g_z_0_xyyy_yyyyy[k] * cd_x[k] + g_z_0_xyyy_xyyyyy[k];

                g_z_0_xxyyy_yyyyz[k] = -g_z_0_xyyy_yyyyz[k] * cd_x[k] + g_z_0_xyyy_xyyyyz[k];

                g_z_0_xxyyy_yyyzz[k] = -g_z_0_xyyy_yyyzz[k] * cd_x[k] + g_z_0_xyyy_xyyyzz[k];

                g_z_0_xxyyy_yyzzz[k] = -g_z_0_xyyy_yyzzz[k] * cd_x[k] + g_z_0_xyyy_xyyzzz[k];

                g_z_0_xxyyy_yzzzz[k] = -g_z_0_xyyy_yzzzz[k] * cd_x[k] + g_z_0_xyyy_xyzzzz[k];

                g_z_0_xxyyy_zzzzz[k] = -g_z_0_xyyy_zzzzz[k] * cd_x[k] + g_z_0_xyyy_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 147);

            auto g_z_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 148);

            auto g_z_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 149);

            auto g_z_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 150);

            auto g_z_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 151);

            auto g_z_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 152);

            auto g_z_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 153);

            auto g_z_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 154);

            auto g_z_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 155);

            auto g_z_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 156);

            auto g_z_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 157);

            auto g_z_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 158);

            auto g_z_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 159);

            auto g_z_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 160);

            auto g_z_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 161);

            auto g_z_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 162);

            auto g_z_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 163);

            auto g_z_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 164);

            auto g_z_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 165);

            auto g_z_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 166);

            auto g_z_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyz_xxxxx, g_z_0_xxyyz_xxxxy, g_z_0_xxyyz_xxxxz, g_z_0_xxyyz_xxxyy, g_z_0_xxyyz_xxxyz, g_z_0_xxyyz_xxxzz, g_z_0_xxyyz_xxyyy, g_z_0_xxyyz_xxyyz, g_z_0_xxyyz_xxyzz, g_z_0_xxyyz_xxzzz, g_z_0_xxyyz_xyyyy, g_z_0_xxyyz_xyyyz, g_z_0_xxyyz_xyyzz, g_z_0_xxyyz_xyzzz, g_z_0_xxyyz_xzzzz, g_z_0_xxyyz_yyyyy, g_z_0_xxyyz_yyyyz, g_z_0_xxyyz_yyyzz, g_z_0_xxyyz_yyzzz, g_z_0_xxyyz_yzzzz, g_z_0_xxyyz_zzzzz, g_z_0_xyyz_xxxxx, g_z_0_xyyz_xxxxxx, g_z_0_xyyz_xxxxxy, g_z_0_xyyz_xxxxxz, g_z_0_xyyz_xxxxy, g_z_0_xyyz_xxxxyy, g_z_0_xyyz_xxxxyz, g_z_0_xyyz_xxxxz, g_z_0_xyyz_xxxxzz, g_z_0_xyyz_xxxyy, g_z_0_xyyz_xxxyyy, g_z_0_xyyz_xxxyyz, g_z_0_xyyz_xxxyz, g_z_0_xyyz_xxxyzz, g_z_0_xyyz_xxxzz, g_z_0_xyyz_xxxzzz, g_z_0_xyyz_xxyyy, g_z_0_xyyz_xxyyyy, g_z_0_xyyz_xxyyyz, g_z_0_xyyz_xxyyz, g_z_0_xyyz_xxyyzz, g_z_0_xyyz_xxyzz, g_z_0_xyyz_xxyzzz, g_z_0_xyyz_xxzzz, g_z_0_xyyz_xxzzzz, g_z_0_xyyz_xyyyy, g_z_0_xyyz_xyyyyy, g_z_0_xyyz_xyyyyz, g_z_0_xyyz_xyyyz, g_z_0_xyyz_xyyyzz, g_z_0_xyyz_xyyzz, g_z_0_xyyz_xyyzzz, g_z_0_xyyz_xyzzz, g_z_0_xyyz_xyzzzz, g_z_0_xyyz_xzzzz, g_z_0_xyyz_xzzzzz, g_z_0_xyyz_yyyyy, g_z_0_xyyz_yyyyz, g_z_0_xyyz_yyyzz, g_z_0_xyyz_yyzzz, g_z_0_xyyz_yzzzz, g_z_0_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xxxxx[k] = -g_z_0_xyyz_xxxxx[k] * cd_x[k] + g_z_0_xyyz_xxxxxx[k];

                g_z_0_xxyyz_xxxxy[k] = -g_z_0_xyyz_xxxxy[k] * cd_x[k] + g_z_0_xyyz_xxxxxy[k];

                g_z_0_xxyyz_xxxxz[k] = -g_z_0_xyyz_xxxxz[k] * cd_x[k] + g_z_0_xyyz_xxxxxz[k];

                g_z_0_xxyyz_xxxyy[k] = -g_z_0_xyyz_xxxyy[k] * cd_x[k] + g_z_0_xyyz_xxxxyy[k];

                g_z_0_xxyyz_xxxyz[k] = -g_z_0_xyyz_xxxyz[k] * cd_x[k] + g_z_0_xyyz_xxxxyz[k];

                g_z_0_xxyyz_xxxzz[k] = -g_z_0_xyyz_xxxzz[k] * cd_x[k] + g_z_0_xyyz_xxxxzz[k];

                g_z_0_xxyyz_xxyyy[k] = -g_z_0_xyyz_xxyyy[k] * cd_x[k] + g_z_0_xyyz_xxxyyy[k];

                g_z_0_xxyyz_xxyyz[k] = -g_z_0_xyyz_xxyyz[k] * cd_x[k] + g_z_0_xyyz_xxxyyz[k];

                g_z_0_xxyyz_xxyzz[k] = -g_z_0_xyyz_xxyzz[k] * cd_x[k] + g_z_0_xyyz_xxxyzz[k];

                g_z_0_xxyyz_xxzzz[k] = -g_z_0_xyyz_xxzzz[k] * cd_x[k] + g_z_0_xyyz_xxxzzz[k];

                g_z_0_xxyyz_xyyyy[k] = -g_z_0_xyyz_xyyyy[k] * cd_x[k] + g_z_0_xyyz_xxyyyy[k];

                g_z_0_xxyyz_xyyyz[k] = -g_z_0_xyyz_xyyyz[k] * cd_x[k] + g_z_0_xyyz_xxyyyz[k];

                g_z_0_xxyyz_xyyzz[k] = -g_z_0_xyyz_xyyzz[k] * cd_x[k] + g_z_0_xyyz_xxyyzz[k];

                g_z_0_xxyyz_xyzzz[k] = -g_z_0_xyyz_xyzzz[k] * cd_x[k] + g_z_0_xyyz_xxyzzz[k];

                g_z_0_xxyyz_xzzzz[k] = -g_z_0_xyyz_xzzzz[k] * cd_x[k] + g_z_0_xyyz_xxzzzz[k];

                g_z_0_xxyyz_yyyyy[k] = -g_z_0_xyyz_yyyyy[k] * cd_x[k] + g_z_0_xyyz_xyyyyy[k];

                g_z_0_xxyyz_yyyyz[k] = -g_z_0_xyyz_yyyyz[k] * cd_x[k] + g_z_0_xyyz_xyyyyz[k];

                g_z_0_xxyyz_yyyzz[k] = -g_z_0_xyyz_yyyzz[k] * cd_x[k] + g_z_0_xyyz_xyyyzz[k];

                g_z_0_xxyyz_yyzzz[k] = -g_z_0_xyyz_yyzzz[k] * cd_x[k] + g_z_0_xyyz_xyyzzz[k];

                g_z_0_xxyyz_yzzzz[k] = -g_z_0_xyyz_yzzzz[k] * cd_x[k] + g_z_0_xyyz_xyzzzz[k];

                g_z_0_xxyyz_zzzzz[k] = -g_z_0_xyyz_zzzzz[k] * cd_x[k] + g_z_0_xyyz_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 168);

            auto g_z_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 169);

            auto g_z_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 170);

            auto g_z_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 171);

            auto g_z_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 172);

            auto g_z_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 173);

            auto g_z_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 174);

            auto g_z_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 175);

            auto g_z_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 176);

            auto g_z_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 177);

            auto g_z_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 178);

            auto g_z_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 179);

            auto g_z_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 180);

            auto g_z_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 181);

            auto g_z_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 182);

            auto g_z_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 183);

            auto g_z_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 184);

            auto g_z_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 185);

            auto g_z_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 186);

            auto g_z_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 187);

            auto g_z_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzz_xxxxx, g_z_0_xxyzz_xxxxy, g_z_0_xxyzz_xxxxz, g_z_0_xxyzz_xxxyy, g_z_0_xxyzz_xxxyz, g_z_0_xxyzz_xxxzz, g_z_0_xxyzz_xxyyy, g_z_0_xxyzz_xxyyz, g_z_0_xxyzz_xxyzz, g_z_0_xxyzz_xxzzz, g_z_0_xxyzz_xyyyy, g_z_0_xxyzz_xyyyz, g_z_0_xxyzz_xyyzz, g_z_0_xxyzz_xyzzz, g_z_0_xxyzz_xzzzz, g_z_0_xxyzz_yyyyy, g_z_0_xxyzz_yyyyz, g_z_0_xxyzz_yyyzz, g_z_0_xxyzz_yyzzz, g_z_0_xxyzz_yzzzz, g_z_0_xxyzz_zzzzz, g_z_0_xyzz_xxxxx, g_z_0_xyzz_xxxxxx, g_z_0_xyzz_xxxxxy, g_z_0_xyzz_xxxxxz, g_z_0_xyzz_xxxxy, g_z_0_xyzz_xxxxyy, g_z_0_xyzz_xxxxyz, g_z_0_xyzz_xxxxz, g_z_0_xyzz_xxxxzz, g_z_0_xyzz_xxxyy, g_z_0_xyzz_xxxyyy, g_z_0_xyzz_xxxyyz, g_z_0_xyzz_xxxyz, g_z_0_xyzz_xxxyzz, g_z_0_xyzz_xxxzz, g_z_0_xyzz_xxxzzz, g_z_0_xyzz_xxyyy, g_z_0_xyzz_xxyyyy, g_z_0_xyzz_xxyyyz, g_z_0_xyzz_xxyyz, g_z_0_xyzz_xxyyzz, g_z_0_xyzz_xxyzz, g_z_0_xyzz_xxyzzz, g_z_0_xyzz_xxzzz, g_z_0_xyzz_xxzzzz, g_z_0_xyzz_xyyyy, g_z_0_xyzz_xyyyyy, g_z_0_xyzz_xyyyyz, g_z_0_xyzz_xyyyz, g_z_0_xyzz_xyyyzz, g_z_0_xyzz_xyyzz, g_z_0_xyzz_xyyzzz, g_z_0_xyzz_xyzzz, g_z_0_xyzz_xyzzzz, g_z_0_xyzz_xzzzz, g_z_0_xyzz_xzzzzz, g_z_0_xyzz_yyyyy, g_z_0_xyzz_yyyyz, g_z_0_xyzz_yyyzz, g_z_0_xyzz_yyzzz, g_z_0_xyzz_yzzzz, g_z_0_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xxxxx[k] = -g_z_0_xyzz_xxxxx[k] * cd_x[k] + g_z_0_xyzz_xxxxxx[k];

                g_z_0_xxyzz_xxxxy[k] = -g_z_0_xyzz_xxxxy[k] * cd_x[k] + g_z_0_xyzz_xxxxxy[k];

                g_z_0_xxyzz_xxxxz[k] = -g_z_0_xyzz_xxxxz[k] * cd_x[k] + g_z_0_xyzz_xxxxxz[k];

                g_z_0_xxyzz_xxxyy[k] = -g_z_0_xyzz_xxxyy[k] * cd_x[k] + g_z_0_xyzz_xxxxyy[k];

                g_z_0_xxyzz_xxxyz[k] = -g_z_0_xyzz_xxxyz[k] * cd_x[k] + g_z_0_xyzz_xxxxyz[k];

                g_z_0_xxyzz_xxxzz[k] = -g_z_0_xyzz_xxxzz[k] * cd_x[k] + g_z_0_xyzz_xxxxzz[k];

                g_z_0_xxyzz_xxyyy[k] = -g_z_0_xyzz_xxyyy[k] * cd_x[k] + g_z_0_xyzz_xxxyyy[k];

                g_z_0_xxyzz_xxyyz[k] = -g_z_0_xyzz_xxyyz[k] * cd_x[k] + g_z_0_xyzz_xxxyyz[k];

                g_z_0_xxyzz_xxyzz[k] = -g_z_0_xyzz_xxyzz[k] * cd_x[k] + g_z_0_xyzz_xxxyzz[k];

                g_z_0_xxyzz_xxzzz[k] = -g_z_0_xyzz_xxzzz[k] * cd_x[k] + g_z_0_xyzz_xxxzzz[k];

                g_z_0_xxyzz_xyyyy[k] = -g_z_0_xyzz_xyyyy[k] * cd_x[k] + g_z_0_xyzz_xxyyyy[k];

                g_z_0_xxyzz_xyyyz[k] = -g_z_0_xyzz_xyyyz[k] * cd_x[k] + g_z_0_xyzz_xxyyyz[k];

                g_z_0_xxyzz_xyyzz[k] = -g_z_0_xyzz_xyyzz[k] * cd_x[k] + g_z_0_xyzz_xxyyzz[k];

                g_z_0_xxyzz_xyzzz[k] = -g_z_0_xyzz_xyzzz[k] * cd_x[k] + g_z_0_xyzz_xxyzzz[k];

                g_z_0_xxyzz_xzzzz[k] = -g_z_0_xyzz_xzzzz[k] * cd_x[k] + g_z_0_xyzz_xxzzzz[k];

                g_z_0_xxyzz_yyyyy[k] = -g_z_0_xyzz_yyyyy[k] * cd_x[k] + g_z_0_xyzz_xyyyyy[k];

                g_z_0_xxyzz_yyyyz[k] = -g_z_0_xyzz_yyyyz[k] * cd_x[k] + g_z_0_xyzz_xyyyyz[k];

                g_z_0_xxyzz_yyyzz[k] = -g_z_0_xyzz_yyyzz[k] * cd_x[k] + g_z_0_xyzz_xyyyzz[k];

                g_z_0_xxyzz_yyzzz[k] = -g_z_0_xyzz_yyzzz[k] * cd_x[k] + g_z_0_xyzz_xyyzzz[k];

                g_z_0_xxyzz_yzzzz[k] = -g_z_0_xyzz_yzzzz[k] * cd_x[k] + g_z_0_xyzz_xyzzzz[k];

                g_z_0_xxyzz_zzzzz[k] = -g_z_0_xyzz_zzzzz[k] * cd_x[k] + g_z_0_xyzz_xzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 189);

            auto g_z_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 190);

            auto g_z_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 191);

            auto g_z_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 192);

            auto g_z_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 193);

            auto g_z_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 194);

            auto g_z_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 195);

            auto g_z_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 196);

            auto g_z_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 197);

            auto g_z_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 198);

            auto g_z_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 199);

            auto g_z_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 200);

            auto g_z_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 201);

            auto g_z_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 202);

            auto g_z_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 203);

            auto g_z_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 204);

            auto g_z_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 205);

            auto g_z_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 206);

            auto g_z_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 207);

            auto g_z_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 208);

            auto g_z_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzz_xxxxx, g_z_0_xxzzz_xxxxy, g_z_0_xxzzz_xxxxz, g_z_0_xxzzz_xxxyy, g_z_0_xxzzz_xxxyz, g_z_0_xxzzz_xxxzz, g_z_0_xxzzz_xxyyy, g_z_0_xxzzz_xxyyz, g_z_0_xxzzz_xxyzz, g_z_0_xxzzz_xxzzz, g_z_0_xxzzz_xyyyy, g_z_0_xxzzz_xyyyz, g_z_0_xxzzz_xyyzz, g_z_0_xxzzz_xyzzz, g_z_0_xxzzz_xzzzz, g_z_0_xxzzz_yyyyy, g_z_0_xxzzz_yyyyz, g_z_0_xxzzz_yyyzz, g_z_0_xxzzz_yyzzz, g_z_0_xxzzz_yzzzz, g_z_0_xxzzz_zzzzz, g_z_0_xzzz_xxxxx, g_z_0_xzzz_xxxxxx, g_z_0_xzzz_xxxxxy, g_z_0_xzzz_xxxxxz, g_z_0_xzzz_xxxxy, g_z_0_xzzz_xxxxyy, g_z_0_xzzz_xxxxyz, g_z_0_xzzz_xxxxz, g_z_0_xzzz_xxxxzz, g_z_0_xzzz_xxxyy, g_z_0_xzzz_xxxyyy, g_z_0_xzzz_xxxyyz, g_z_0_xzzz_xxxyz, g_z_0_xzzz_xxxyzz, g_z_0_xzzz_xxxzz, g_z_0_xzzz_xxxzzz, g_z_0_xzzz_xxyyy, g_z_0_xzzz_xxyyyy, g_z_0_xzzz_xxyyyz, g_z_0_xzzz_xxyyz, g_z_0_xzzz_xxyyzz, g_z_0_xzzz_xxyzz, g_z_0_xzzz_xxyzzz, g_z_0_xzzz_xxzzz, g_z_0_xzzz_xxzzzz, g_z_0_xzzz_xyyyy, g_z_0_xzzz_xyyyyy, g_z_0_xzzz_xyyyyz, g_z_0_xzzz_xyyyz, g_z_0_xzzz_xyyyzz, g_z_0_xzzz_xyyzz, g_z_0_xzzz_xyyzzz, g_z_0_xzzz_xyzzz, g_z_0_xzzz_xyzzzz, g_z_0_xzzz_xzzzz, g_z_0_xzzz_xzzzzz, g_z_0_xzzz_yyyyy, g_z_0_xzzz_yyyyz, g_z_0_xzzz_yyyzz, g_z_0_xzzz_yyzzz, g_z_0_xzzz_yzzzz, g_z_0_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xxxxx[k] = -g_z_0_xzzz_xxxxx[k] * cd_x[k] + g_z_0_xzzz_xxxxxx[k];

                g_z_0_xxzzz_xxxxy[k] = -g_z_0_xzzz_xxxxy[k] * cd_x[k] + g_z_0_xzzz_xxxxxy[k];

                g_z_0_xxzzz_xxxxz[k] = -g_z_0_xzzz_xxxxz[k] * cd_x[k] + g_z_0_xzzz_xxxxxz[k];

                g_z_0_xxzzz_xxxyy[k] = -g_z_0_xzzz_xxxyy[k] * cd_x[k] + g_z_0_xzzz_xxxxyy[k];

                g_z_0_xxzzz_xxxyz[k] = -g_z_0_xzzz_xxxyz[k] * cd_x[k] + g_z_0_xzzz_xxxxyz[k];

                g_z_0_xxzzz_xxxzz[k] = -g_z_0_xzzz_xxxzz[k] * cd_x[k] + g_z_0_xzzz_xxxxzz[k];

                g_z_0_xxzzz_xxyyy[k] = -g_z_0_xzzz_xxyyy[k] * cd_x[k] + g_z_0_xzzz_xxxyyy[k];

                g_z_0_xxzzz_xxyyz[k] = -g_z_0_xzzz_xxyyz[k] * cd_x[k] + g_z_0_xzzz_xxxyyz[k];

                g_z_0_xxzzz_xxyzz[k] = -g_z_0_xzzz_xxyzz[k] * cd_x[k] + g_z_0_xzzz_xxxyzz[k];

                g_z_0_xxzzz_xxzzz[k] = -g_z_0_xzzz_xxzzz[k] * cd_x[k] + g_z_0_xzzz_xxxzzz[k];

                g_z_0_xxzzz_xyyyy[k] = -g_z_0_xzzz_xyyyy[k] * cd_x[k] + g_z_0_xzzz_xxyyyy[k];

                g_z_0_xxzzz_xyyyz[k] = -g_z_0_xzzz_xyyyz[k] * cd_x[k] + g_z_0_xzzz_xxyyyz[k];

                g_z_0_xxzzz_xyyzz[k] = -g_z_0_xzzz_xyyzz[k] * cd_x[k] + g_z_0_xzzz_xxyyzz[k];

                g_z_0_xxzzz_xyzzz[k] = -g_z_0_xzzz_xyzzz[k] * cd_x[k] + g_z_0_xzzz_xxyzzz[k];

                g_z_0_xxzzz_xzzzz[k] = -g_z_0_xzzz_xzzzz[k] * cd_x[k] + g_z_0_xzzz_xxzzzz[k];

                g_z_0_xxzzz_yyyyy[k] = -g_z_0_xzzz_yyyyy[k] * cd_x[k] + g_z_0_xzzz_xyyyyy[k];

                g_z_0_xxzzz_yyyyz[k] = -g_z_0_xzzz_yyyyz[k] * cd_x[k] + g_z_0_xzzz_xyyyyz[k];

                g_z_0_xxzzz_yyyzz[k] = -g_z_0_xzzz_yyyzz[k] * cd_x[k] + g_z_0_xzzz_xyyyzz[k];

                g_z_0_xxzzz_yyzzz[k] = -g_z_0_xzzz_yyzzz[k] * cd_x[k] + g_z_0_xzzz_xyyzzz[k];

                g_z_0_xxzzz_yzzzz[k] = -g_z_0_xzzz_yzzzz[k] * cd_x[k] + g_z_0_xzzz_xyzzzz[k];

                g_z_0_xxzzz_zzzzz[k] = -g_z_0_xzzz_zzzzz[k] * cd_x[k] + g_z_0_xzzz_xzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 210);

            auto g_z_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 211);

            auto g_z_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 212);

            auto g_z_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 213);

            auto g_z_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 214);

            auto g_z_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 215);

            auto g_z_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 216);

            auto g_z_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 217);

            auto g_z_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 218);

            auto g_z_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 219);

            auto g_z_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 220);

            auto g_z_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 221);

            auto g_z_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 222);

            auto g_z_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 223);

            auto g_z_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 224);

            auto g_z_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 225);

            auto g_z_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 226);

            auto g_z_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 227);

            auto g_z_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 228);

            auto g_z_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 229);

            auto g_z_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyy_xxxxx, g_z_0_xyyyy_xxxxy, g_z_0_xyyyy_xxxxz, g_z_0_xyyyy_xxxyy, g_z_0_xyyyy_xxxyz, g_z_0_xyyyy_xxxzz, g_z_0_xyyyy_xxyyy, g_z_0_xyyyy_xxyyz, g_z_0_xyyyy_xxyzz, g_z_0_xyyyy_xxzzz, g_z_0_xyyyy_xyyyy, g_z_0_xyyyy_xyyyz, g_z_0_xyyyy_xyyzz, g_z_0_xyyyy_xyzzz, g_z_0_xyyyy_xzzzz, g_z_0_xyyyy_yyyyy, g_z_0_xyyyy_yyyyz, g_z_0_xyyyy_yyyzz, g_z_0_xyyyy_yyzzz, g_z_0_xyyyy_yzzzz, g_z_0_xyyyy_zzzzz, g_z_0_yyyy_xxxxx, g_z_0_yyyy_xxxxxx, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxxz, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxz, g_z_0_yyyy_xxxxzz, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxzz, g_z_0_yyyy_xxxzzz, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxzzz, g_z_0_yyyy_xxzzzz, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xzzzz, g_z_0_yyyy_xzzzzz, g_z_0_yyyy_yyyyy, g_z_0_yyyy_yyyyz, g_z_0_yyyy_yyyzz, g_z_0_yyyy_yyzzz, g_z_0_yyyy_yzzzz, g_z_0_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xxxxx[k] = -g_z_0_yyyy_xxxxx[k] * cd_x[k] + g_z_0_yyyy_xxxxxx[k];

                g_z_0_xyyyy_xxxxy[k] = -g_z_0_yyyy_xxxxy[k] * cd_x[k] + g_z_0_yyyy_xxxxxy[k];

                g_z_0_xyyyy_xxxxz[k] = -g_z_0_yyyy_xxxxz[k] * cd_x[k] + g_z_0_yyyy_xxxxxz[k];

                g_z_0_xyyyy_xxxyy[k] = -g_z_0_yyyy_xxxyy[k] * cd_x[k] + g_z_0_yyyy_xxxxyy[k];

                g_z_0_xyyyy_xxxyz[k] = -g_z_0_yyyy_xxxyz[k] * cd_x[k] + g_z_0_yyyy_xxxxyz[k];

                g_z_0_xyyyy_xxxzz[k] = -g_z_0_yyyy_xxxzz[k] * cd_x[k] + g_z_0_yyyy_xxxxzz[k];

                g_z_0_xyyyy_xxyyy[k] = -g_z_0_yyyy_xxyyy[k] * cd_x[k] + g_z_0_yyyy_xxxyyy[k];

                g_z_0_xyyyy_xxyyz[k] = -g_z_0_yyyy_xxyyz[k] * cd_x[k] + g_z_0_yyyy_xxxyyz[k];

                g_z_0_xyyyy_xxyzz[k] = -g_z_0_yyyy_xxyzz[k] * cd_x[k] + g_z_0_yyyy_xxxyzz[k];

                g_z_0_xyyyy_xxzzz[k] = -g_z_0_yyyy_xxzzz[k] * cd_x[k] + g_z_0_yyyy_xxxzzz[k];

                g_z_0_xyyyy_xyyyy[k] = -g_z_0_yyyy_xyyyy[k] * cd_x[k] + g_z_0_yyyy_xxyyyy[k];

                g_z_0_xyyyy_xyyyz[k] = -g_z_0_yyyy_xyyyz[k] * cd_x[k] + g_z_0_yyyy_xxyyyz[k];

                g_z_0_xyyyy_xyyzz[k] = -g_z_0_yyyy_xyyzz[k] * cd_x[k] + g_z_0_yyyy_xxyyzz[k];

                g_z_0_xyyyy_xyzzz[k] = -g_z_0_yyyy_xyzzz[k] * cd_x[k] + g_z_0_yyyy_xxyzzz[k];

                g_z_0_xyyyy_xzzzz[k] = -g_z_0_yyyy_xzzzz[k] * cd_x[k] + g_z_0_yyyy_xxzzzz[k];

                g_z_0_xyyyy_yyyyy[k] = -g_z_0_yyyy_yyyyy[k] * cd_x[k] + g_z_0_yyyy_xyyyyy[k];

                g_z_0_xyyyy_yyyyz[k] = -g_z_0_yyyy_yyyyz[k] * cd_x[k] + g_z_0_yyyy_xyyyyz[k];

                g_z_0_xyyyy_yyyzz[k] = -g_z_0_yyyy_yyyzz[k] * cd_x[k] + g_z_0_yyyy_xyyyzz[k];

                g_z_0_xyyyy_yyzzz[k] = -g_z_0_yyyy_yyzzz[k] * cd_x[k] + g_z_0_yyyy_xyyzzz[k];

                g_z_0_xyyyy_yzzzz[k] = -g_z_0_yyyy_yzzzz[k] * cd_x[k] + g_z_0_yyyy_xyzzzz[k];

                g_z_0_xyyyy_zzzzz[k] = -g_z_0_yyyy_zzzzz[k] * cd_x[k] + g_z_0_yyyy_xzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 231);

            auto g_z_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 232);

            auto g_z_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 233);

            auto g_z_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 234);

            auto g_z_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 235);

            auto g_z_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 236);

            auto g_z_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 237);

            auto g_z_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 238);

            auto g_z_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 239);

            auto g_z_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 240);

            auto g_z_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 241);

            auto g_z_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 242);

            auto g_z_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 243);

            auto g_z_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 244);

            auto g_z_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 245);

            auto g_z_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 246);

            auto g_z_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 247);

            auto g_z_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 248);

            auto g_z_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 249);

            auto g_z_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 250);

            auto g_z_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyz_xxxxx, g_z_0_xyyyz_xxxxy, g_z_0_xyyyz_xxxxz, g_z_0_xyyyz_xxxyy, g_z_0_xyyyz_xxxyz, g_z_0_xyyyz_xxxzz, g_z_0_xyyyz_xxyyy, g_z_0_xyyyz_xxyyz, g_z_0_xyyyz_xxyzz, g_z_0_xyyyz_xxzzz, g_z_0_xyyyz_xyyyy, g_z_0_xyyyz_xyyyz, g_z_0_xyyyz_xyyzz, g_z_0_xyyyz_xyzzz, g_z_0_xyyyz_xzzzz, g_z_0_xyyyz_yyyyy, g_z_0_xyyyz_yyyyz, g_z_0_xyyyz_yyyzz, g_z_0_xyyyz_yyzzz, g_z_0_xyyyz_yzzzz, g_z_0_xyyyz_zzzzz, g_z_0_yyyz_xxxxx, g_z_0_yyyz_xxxxxx, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxxz, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxz, g_z_0_yyyz_xxxxzz, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxzz, g_z_0_yyyz_xxxzzz, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxzzz, g_z_0_yyyz_xxzzzz, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xzzzz, g_z_0_yyyz_xzzzzz, g_z_0_yyyz_yyyyy, g_z_0_yyyz_yyyyz, g_z_0_yyyz_yyyzz, g_z_0_yyyz_yyzzz, g_z_0_yyyz_yzzzz, g_z_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xxxxx[k] = -g_z_0_yyyz_xxxxx[k] * cd_x[k] + g_z_0_yyyz_xxxxxx[k];

                g_z_0_xyyyz_xxxxy[k] = -g_z_0_yyyz_xxxxy[k] * cd_x[k] + g_z_0_yyyz_xxxxxy[k];

                g_z_0_xyyyz_xxxxz[k] = -g_z_0_yyyz_xxxxz[k] * cd_x[k] + g_z_0_yyyz_xxxxxz[k];

                g_z_0_xyyyz_xxxyy[k] = -g_z_0_yyyz_xxxyy[k] * cd_x[k] + g_z_0_yyyz_xxxxyy[k];

                g_z_0_xyyyz_xxxyz[k] = -g_z_0_yyyz_xxxyz[k] * cd_x[k] + g_z_0_yyyz_xxxxyz[k];

                g_z_0_xyyyz_xxxzz[k] = -g_z_0_yyyz_xxxzz[k] * cd_x[k] + g_z_0_yyyz_xxxxzz[k];

                g_z_0_xyyyz_xxyyy[k] = -g_z_0_yyyz_xxyyy[k] * cd_x[k] + g_z_0_yyyz_xxxyyy[k];

                g_z_0_xyyyz_xxyyz[k] = -g_z_0_yyyz_xxyyz[k] * cd_x[k] + g_z_0_yyyz_xxxyyz[k];

                g_z_0_xyyyz_xxyzz[k] = -g_z_0_yyyz_xxyzz[k] * cd_x[k] + g_z_0_yyyz_xxxyzz[k];

                g_z_0_xyyyz_xxzzz[k] = -g_z_0_yyyz_xxzzz[k] * cd_x[k] + g_z_0_yyyz_xxxzzz[k];

                g_z_0_xyyyz_xyyyy[k] = -g_z_0_yyyz_xyyyy[k] * cd_x[k] + g_z_0_yyyz_xxyyyy[k];

                g_z_0_xyyyz_xyyyz[k] = -g_z_0_yyyz_xyyyz[k] * cd_x[k] + g_z_0_yyyz_xxyyyz[k];

                g_z_0_xyyyz_xyyzz[k] = -g_z_0_yyyz_xyyzz[k] * cd_x[k] + g_z_0_yyyz_xxyyzz[k];

                g_z_0_xyyyz_xyzzz[k] = -g_z_0_yyyz_xyzzz[k] * cd_x[k] + g_z_0_yyyz_xxyzzz[k];

                g_z_0_xyyyz_xzzzz[k] = -g_z_0_yyyz_xzzzz[k] * cd_x[k] + g_z_0_yyyz_xxzzzz[k];

                g_z_0_xyyyz_yyyyy[k] = -g_z_0_yyyz_yyyyy[k] * cd_x[k] + g_z_0_yyyz_xyyyyy[k];

                g_z_0_xyyyz_yyyyz[k] = -g_z_0_yyyz_yyyyz[k] * cd_x[k] + g_z_0_yyyz_xyyyyz[k];

                g_z_0_xyyyz_yyyzz[k] = -g_z_0_yyyz_yyyzz[k] * cd_x[k] + g_z_0_yyyz_xyyyzz[k];

                g_z_0_xyyyz_yyzzz[k] = -g_z_0_yyyz_yyzzz[k] * cd_x[k] + g_z_0_yyyz_xyyzzz[k];

                g_z_0_xyyyz_yzzzz[k] = -g_z_0_yyyz_yzzzz[k] * cd_x[k] + g_z_0_yyyz_xyzzzz[k];

                g_z_0_xyyyz_zzzzz[k] = -g_z_0_yyyz_zzzzz[k] * cd_x[k] + g_z_0_yyyz_xzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 252);

            auto g_z_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 253);

            auto g_z_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 254);

            auto g_z_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 255);

            auto g_z_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 256);

            auto g_z_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 257);

            auto g_z_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 258);

            auto g_z_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 259);

            auto g_z_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 260);

            auto g_z_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 261);

            auto g_z_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 262);

            auto g_z_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 263);

            auto g_z_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 264);

            auto g_z_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 265);

            auto g_z_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 266);

            auto g_z_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 267);

            auto g_z_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 268);

            auto g_z_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 269);

            auto g_z_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 270);

            auto g_z_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 271);

            auto g_z_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzz_xxxxx, g_z_0_xyyzz_xxxxy, g_z_0_xyyzz_xxxxz, g_z_0_xyyzz_xxxyy, g_z_0_xyyzz_xxxyz, g_z_0_xyyzz_xxxzz, g_z_0_xyyzz_xxyyy, g_z_0_xyyzz_xxyyz, g_z_0_xyyzz_xxyzz, g_z_0_xyyzz_xxzzz, g_z_0_xyyzz_xyyyy, g_z_0_xyyzz_xyyyz, g_z_0_xyyzz_xyyzz, g_z_0_xyyzz_xyzzz, g_z_0_xyyzz_xzzzz, g_z_0_xyyzz_yyyyy, g_z_0_xyyzz_yyyyz, g_z_0_xyyzz_yyyzz, g_z_0_xyyzz_yyzzz, g_z_0_xyyzz_yzzzz, g_z_0_xyyzz_zzzzz, g_z_0_yyzz_xxxxx, g_z_0_yyzz_xxxxxx, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxxz, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxz, g_z_0_yyzz_xxxxzz, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxzz, g_z_0_yyzz_xxxzzz, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxzzz, g_z_0_yyzz_xxzzzz, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xzzzz, g_z_0_yyzz_xzzzzz, g_z_0_yyzz_yyyyy, g_z_0_yyzz_yyyyz, g_z_0_yyzz_yyyzz, g_z_0_yyzz_yyzzz, g_z_0_yyzz_yzzzz, g_z_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xxxxx[k] = -g_z_0_yyzz_xxxxx[k] * cd_x[k] + g_z_0_yyzz_xxxxxx[k];

                g_z_0_xyyzz_xxxxy[k] = -g_z_0_yyzz_xxxxy[k] * cd_x[k] + g_z_0_yyzz_xxxxxy[k];

                g_z_0_xyyzz_xxxxz[k] = -g_z_0_yyzz_xxxxz[k] * cd_x[k] + g_z_0_yyzz_xxxxxz[k];

                g_z_0_xyyzz_xxxyy[k] = -g_z_0_yyzz_xxxyy[k] * cd_x[k] + g_z_0_yyzz_xxxxyy[k];

                g_z_0_xyyzz_xxxyz[k] = -g_z_0_yyzz_xxxyz[k] * cd_x[k] + g_z_0_yyzz_xxxxyz[k];

                g_z_0_xyyzz_xxxzz[k] = -g_z_0_yyzz_xxxzz[k] * cd_x[k] + g_z_0_yyzz_xxxxzz[k];

                g_z_0_xyyzz_xxyyy[k] = -g_z_0_yyzz_xxyyy[k] * cd_x[k] + g_z_0_yyzz_xxxyyy[k];

                g_z_0_xyyzz_xxyyz[k] = -g_z_0_yyzz_xxyyz[k] * cd_x[k] + g_z_0_yyzz_xxxyyz[k];

                g_z_0_xyyzz_xxyzz[k] = -g_z_0_yyzz_xxyzz[k] * cd_x[k] + g_z_0_yyzz_xxxyzz[k];

                g_z_0_xyyzz_xxzzz[k] = -g_z_0_yyzz_xxzzz[k] * cd_x[k] + g_z_0_yyzz_xxxzzz[k];

                g_z_0_xyyzz_xyyyy[k] = -g_z_0_yyzz_xyyyy[k] * cd_x[k] + g_z_0_yyzz_xxyyyy[k];

                g_z_0_xyyzz_xyyyz[k] = -g_z_0_yyzz_xyyyz[k] * cd_x[k] + g_z_0_yyzz_xxyyyz[k];

                g_z_0_xyyzz_xyyzz[k] = -g_z_0_yyzz_xyyzz[k] * cd_x[k] + g_z_0_yyzz_xxyyzz[k];

                g_z_0_xyyzz_xyzzz[k] = -g_z_0_yyzz_xyzzz[k] * cd_x[k] + g_z_0_yyzz_xxyzzz[k];

                g_z_0_xyyzz_xzzzz[k] = -g_z_0_yyzz_xzzzz[k] * cd_x[k] + g_z_0_yyzz_xxzzzz[k];

                g_z_0_xyyzz_yyyyy[k] = -g_z_0_yyzz_yyyyy[k] * cd_x[k] + g_z_0_yyzz_xyyyyy[k];

                g_z_0_xyyzz_yyyyz[k] = -g_z_0_yyzz_yyyyz[k] * cd_x[k] + g_z_0_yyzz_xyyyyz[k];

                g_z_0_xyyzz_yyyzz[k] = -g_z_0_yyzz_yyyzz[k] * cd_x[k] + g_z_0_yyzz_xyyyzz[k];

                g_z_0_xyyzz_yyzzz[k] = -g_z_0_yyzz_yyzzz[k] * cd_x[k] + g_z_0_yyzz_xyyzzz[k];

                g_z_0_xyyzz_yzzzz[k] = -g_z_0_yyzz_yzzzz[k] * cd_x[k] + g_z_0_yyzz_xyzzzz[k];

                g_z_0_xyyzz_zzzzz[k] = -g_z_0_yyzz_zzzzz[k] * cd_x[k] + g_z_0_yyzz_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 273);

            auto g_z_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 274);

            auto g_z_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 275);

            auto g_z_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 276);

            auto g_z_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 277);

            auto g_z_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 278);

            auto g_z_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 279);

            auto g_z_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 280);

            auto g_z_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 281);

            auto g_z_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 282);

            auto g_z_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 283);

            auto g_z_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 284);

            auto g_z_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 285);

            auto g_z_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 286);

            auto g_z_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 287);

            auto g_z_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 288);

            auto g_z_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 289);

            auto g_z_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 290);

            auto g_z_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 291);

            auto g_z_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 292);

            auto g_z_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzz_xxxxx, g_z_0_xyzzz_xxxxy, g_z_0_xyzzz_xxxxz, g_z_0_xyzzz_xxxyy, g_z_0_xyzzz_xxxyz, g_z_0_xyzzz_xxxzz, g_z_0_xyzzz_xxyyy, g_z_0_xyzzz_xxyyz, g_z_0_xyzzz_xxyzz, g_z_0_xyzzz_xxzzz, g_z_0_xyzzz_xyyyy, g_z_0_xyzzz_xyyyz, g_z_0_xyzzz_xyyzz, g_z_0_xyzzz_xyzzz, g_z_0_xyzzz_xzzzz, g_z_0_xyzzz_yyyyy, g_z_0_xyzzz_yyyyz, g_z_0_xyzzz_yyyzz, g_z_0_xyzzz_yyzzz, g_z_0_xyzzz_yzzzz, g_z_0_xyzzz_zzzzz, g_z_0_yzzz_xxxxx, g_z_0_yzzz_xxxxxx, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxxz, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxz, g_z_0_yzzz_xxxxzz, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxzz, g_z_0_yzzz_xxxzzz, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxzzz, g_z_0_yzzz_xxzzzz, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xzzzz, g_z_0_yzzz_xzzzzz, g_z_0_yzzz_yyyyy, g_z_0_yzzz_yyyyz, g_z_0_yzzz_yyyzz, g_z_0_yzzz_yyzzz, g_z_0_yzzz_yzzzz, g_z_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xxxxx[k] = -g_z_0_yzzz_xxxxx[k] * cd_x[k] + g_z_0_yzzz_xxxxxx[k];

                g_z_0_xyzzz_xxxxy[k] = -g_z_0_yzzz_xxxxy[k] * cd_x[k] + g_z_0_yzzz_xxxxxy[k];

                g_z_0_xyzzz_xxxxz[k] = -g_z_0_yzzz_xxxxz[k] * cd_x[k] + g_z_0_yzzz_xxxxxz[k];

                g_z_0_xyzzz_xxxyy[k] = -g_z_0_yzzz_xxxyy[k] * cd_x[k] + g_z_0_yzzz_xxxxyy[k];

                g_z_0_xyzzz_xxxyz[k] = -g_z_0_yzzz_xxxyz[k] * cd_x[k] + g_z_0_yzzz_xxxxyz[k];

                g_z_0_xyzzz_xxxzz[k] = -g_z_0_yzzz_xxxzz[k] * cd_x[k] + g_z_0_yzzz_xxxxzz[k];

                g_z_0_xyzzz_xxyyy[k] = -g_z_0_yzzz_xxyyy[k] * cd_x[k] + g_z_0_yzzz_xxxyyy[k];

                g_z_0_xyzzz_xxyyz[k] = -g_z_0_yzzz_xxyyz[k] * cd_x[k] + g_z_0_yzzz_xxxyyz[k];

                g_z_0_xyzzz_xxyzz[k] = -g_z_0_yzzz_xxyzz[k] * cd_x[k] + g_z_0_yzzz_xxxyzz[k];

                g_z_0_xyzzz_xxzzz[k] = -g_z_0_yzzz_xxzzz[k] * cd_x[k] + g_z_0_yzzz_xxxzzz[k];

                g_z_0_xyzzz_xyyyy[k] = -g_z_0_yzzz_xyyyy[k] * cd_x[k] + g_z_0_yzzz_xxyyyy[k];

                g_z_0_xyzzz_xyyyz[k] = -g_z_0_yzzz_xyyyz[k] * cd_x[k] + g_z_0_yzzz_xxyyyz[k];

                g_z_0_xyzzz_xyyzz[k] = -g_z_0_yzzz_xyyzz[k] * cd_x[k] + g_z_0_yzzz_xxyyzz[k];

                g_z_0_xyzzz_xyzzz[k] = -g_z_0_yzzz_xyzzz[k] * cd_x[k] + g_z_0_yzzz_xxyzzz[k];

                g_z_0_xyzzz_xzzzz[k] = -g_z_0_yzzz_xzzzz[k] * cd_x[k] + g_z_0_yzzz_xxzzzz[k];

                g_z_0_xyzzz_yyyyy[k] = -g_z_0_yzzz_yyyyy[k] * cd_x[k] + g_z_0_yzzz_xyyyyy[k];

                g_z_0_xyzzz_yyyyz[k] = -g_z_0_yzzz_yyyyz[k] * cd_x[k] + g_z_0_yzzz_xyyyyz[k];

                g_z_0_xyzzz_yyyzz[k] = -g_z_0_yzzz_yyyzz[k] * cd_x[k] + g_z_0_yzzz_xyyyzz[k];

                g_z_0_xyzzz_yyzzz[k] = -g_z_0_yzzz_yyzzz[k] * cd_x[k] + g_z_0_yzzz_xyyzzz[k];

                g_z_0_xyzzz_yzzzz[k] = -g_z_0_yzzz_yzzzz[k] * cd_x[k] + g_z_0_yzzz_xyzzzz[k];

                g_z_0_xyzzz_zzzzz[k] = -g_z_0_yzzz_zzzzz[k] * cd_x[k] + g_z_0_yzzz_xzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 294);

            auto g_z_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 295);

            auto g_z_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 296);

            auto g_z_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 297);

            auto g_z_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 298);

            auto g_z_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 299);

            auto g_z_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 300);

            auto g_z_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 301);

            auto g_z_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 302);

            auto g_z_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 303);

            auto g_z_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 304);

            auto g_z_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 305);

            auto g_z_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 306);

            auto g_z_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 307);

            auto g_z_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 308);

            auto g_z_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 309);

            auto g_z_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 310);

            auto g_z_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 311);

            auto g_z_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 312);

            auto g_z_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 313);

            auto g_z_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzz_xxxxx, g_z_0_xzzzz_xxxxy, g_z_0_xzzzz_xxxxz, g_z_0_xzzzz_xxxyy, g_z_0_xzzzz_xxxyz, g_z_0_xzzzz_xxxzz, g_z_0_xzzzz_xxyyy, g_z_0_xzzzz_xxyyz, g_z_0_xzzzz_xxyzz, g_z_0_xzzzz_xxzzz, g_z_0_xzzzz_xyyyy, g_z_0_xzzzz_xyyyz, g_z_0_xzzzz_xyyzz, g_z_0_xzzzz_xyzzz, g_z_0_xzzzz_xzzzz, g_z_0_xzzzz_yyyyy, g_z_0_xzzzz_yyyyz, g_z_0_xzzzz_yyyzz, g_z_0_xzzzz_yyzzz, g_z_0_xzzzz_yzzzz, g_z_0_xzzzz_zzzzz, g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xxxxx[k] = -g_z_0_zzzz_xxxxx[k] * cd_x[k] + g_z_0_zzzz_xxxxxx[k];

                g_z_0_xzzzz_xxxxy[k] = -g_z_0_zzzz_xxxxy[k] * cd_x[k] + g_z_0_zzzz_xxxxxy[k];

                g_z_0_xzzzz_xxxxz[k] = -g_z_0_zzzz_xxxxz[k] * cd_x[k] + g_z_0_zzzz_xxxxxz[k];

                g_z_0_xzzzz_xxxyy[k] = -g_z_0_zzzz_xxxyy[k] * cd_x[k] + g_z_0_zzzz_xxxxyy[k];

                g_z_0_xzzzz_xxxyz[k] = -g_z_0_zzzz_xxxyz[k] * cd_x[k] + g_z_0_zzzz_xxxxyz[k];

                g_z_0_xzzzz_xxxzz[k] = -g_z_0_zzzz_xxxzz[k] * cd_x[k] + g_z_0_zzzz_xxxxzz[k];

                g_z_0_xzzzz_xxyyy[k] = -g_z_0_zzzz_xxyyy[k] * cd_x[k] + g_z_0_zzzz_xxxyyy[k];

                g_z_0_xzzzz_xxyyz[k] = -g_z_0_zzzz_xxyyz[k] * cd_x[k] + g_z_0_zzzz_xxxyyz[k];

                g_z_0_xzzzz_xxyzz[k] = -g_z_0_zzzz_xxyzz[k] * cd_x[k] + g_z_0_zzzz_xxxyzz[k];

                g_z_0_xzzzz_xxzzz[k] = -g_z_0_zzzz_xxzzz[k] * cd_x[k] + g_z_0_zzzz_xxxzzz[k];

                g_z_0_xzzzz_xyyyy[k] = -g_z_0_zzzz_xyyyy[k] * cd_x[k] + g_z_0_zzzz_xxyyyy[k];

                g_z_0_xzzzz_xyyyz[k] = -g_z_0_zzzz_xyyyz[k] * cd_x[k] + g_z_0_zzzz_xxyyyz[k];

                g_z_0_xzzzz_xyyzz[k] = -g_z_0_zzzz_xyyzz[k] * cd_x[k] + g_z_0_zzzz_xxyyzz[k];

                g_z_0_xzzzz_xyzzz[k] = -g_z_0_zzzz_xyzzz[k] * cd_x[k] + g_z_0_zzzz_xxyzzz[k];

                g_z_0_xzzzz_xzzzz[k] = -g_z_0_zzzz_xzzzz[k] * cd_x[k] + g_z_0_zzzz_xxzzzz[k];

                g_z_0_xzzzz_yyyyy[k] = -g_z_0_zzzz_yyyyy[k] * cd_x[k] + g_z_0_zzzz_xyyyyy[k];

                g_z_0_xzzzz_yyyyz[k] = -g_z_0_zzzz_yyyyz[k] * cd_x[k] + g_z_0_zzzz_xyyyyz[k];

                g_z_0_xzzzz_yyyzz[k] = -g_z_0_zzzz_yyyzz[k] * cd_x[k] + g_z_0_zzzz_xyyyzz[k];

                g_z_0_xzzzz_yyzzz[k] = -g_z_0_zzzz_yyzzz[k] * cd_x[k] + g_z_0_zzzz_xyyzzz[k];

                g_z_0_xzzzz_yzzzz[k] = -g_z_0_zzzz_yzzzz[k] * cd_x[k] + g_z_0_zzzz_xyzzzz[k];

                g_z_0_xzzzz_zzzzz[k] = -g_z_0_zzzz_zzzzz[k] * cd_x[k] + g_z_0_zzzz_xzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 315);

            auto g_z_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 316);

            auto g_z_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 317);

            auto g_z_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 318);

            auto g_z_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 319);

            auto g_z_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 320);

            auto g_z_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 321);

            auto g_z_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 322);

            auto g_z_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 323);

            auto g_z_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 324);

            auto g_z_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 325);

            auto g_z_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 326);

            auto g_z_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 327);

            auto g_z_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 328);

            auto g_z_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 329);

            auto g_z_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 330);

            auto g_z_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 331);

            auto g_z_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 332);

            auto g_z_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 333);

            auto g_z_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 334);

            auto g_z_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_y, g_z_0_yyyy_xxxxx, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxz, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxzz, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxzzz, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xzzzz, g_z_0_yyyy_yyyyy, g_z_0_yyyy_yyyyyy, g_z_0_yyyy_yyyyyz, g_z_0_yyyy_yyyyz, g_z_0_yyyy_yyyyzz, g_z_0_yyyy_yyyzz, g_z_0_yyyy_yyyzzz, g_z_0_yyyy_yyzzz, g_z_0_yyyy_yyzzzz, g_z_0_yyyy_yzzzz, g_z_0_yyyy_yzzzzz, g_z_0_yyyy_zzzzz, g_z_0_yyyyy_xxxxx, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxxz, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxzz, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxzzz, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xzzzz, g_z_0_yyyyy_yyyyy, g_z_0_yyyyy_yyyyz, g_z_0_yyyyy_yyyzz, g_z_0_yyyyy_yyzzz, g_z_0_yyyyy_yzzzz, g_z_0_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xxxxx[k] = -g_z_0_yyyy_xxxxx[k] * cd_y[k] + g_z_0_yyyy_xxxxxy[k];

                g_z_0_yyyyy_xxxxy[k] = -g_z_0_yyyy_xxxxy[k] * cd_y[k] + g_z_0_yyyy_xxxxyy[k];

                g_z_0_yyyyy_xxxxz[k] = -g_z_0_yyyy_xxxxz[k] * cd_y[k] + g_z_0_yyyy_xxxxyz[k];

                g_z_0_yyyyy_xxxyy[k] = -g_z_0_yyyy_xxxyy[k] * cd_y[k] + g_z_0_yyyy_xxxyyy[k];

                g_z_0_yyyyy_xxxyz[k] = -g_z_0_yyyy_xxxyz[k] * cd_y[k] + g_z_0_yyyy_xxxyyz[k];

                g_z_0_yyyyy_xxxzz[k] = -g_z_0_yyyy_xxxzz[k] * cd_y[k] + g_z_0_yyyy_xxxyzz[k];

                g_z_0_yyyyy_xxyyy[k] = -g_z_0_yyyy_xxyyy[k] * cd_y[k] + g_z_0_yyyy_xxyyyy[k];

                g_z_0_yyyyy_xxyyz[k] = -g_z_0_yyyy_xxyyz[k] * cd_y[k] + g_z_0_yyyy_xxyyyz[k];

                g_z_0_yyyyy_xxyzz[k] = -g_z_0_yyyy_xxyzz[k] * cd_y[k] + g_z_0_yyyy_xxyyzz[k];

                g_z_0_yyyyy_xxzzz[k] = -g_z_0_yyyy_xxzzz[k] * cd_y[k] + g_z_0_yyyy_xxyzzz[k];

                g_z_0_yyyyy_xyyyy[k] = -g_z_0_yyyy_xyyyy[k] * cd_y[k] + g_z_0_yyyy_xyyyyy[k];

                g_z_0_yyyyy_xyyyz[k] = -g_z_0_yyyy_xyyyz[k] * cd_y[k] + g_z_0_yyyy_xyyyyz[k];

                g_z_0_yyyyy_xyyzz[k] = -g_z_0_yyyy_xyyzz[k] * cd_y[k] + g_z_0_yyyy_xyyyzz[k];

                g_z_0_yyyyy_xyzzz[k] = -g_z_0_yyyy_xyzzz[k] * cd_y[k] + g_z_0_yyyy_xyyzzz[k];

                g_z_0_yyyyy_xzzzz[k] = -g_z_0_yyyy_xzzzz[k] * cd_y[k] + g_z_0_yyyy_xyzzzz[k];

                g_z_0_yyyyy_yyyyy[k] = -g_z_0_yyyy_yyyyy[k] * cd_y[k] + g_z_0_yyyy_yyyyyy[k];

                g_z_0_yyyyy_yyyyz[k] = -g_z_0_yyyy_yyyyz[k] * cd_y[k] + g_z_0_yyyy_yyyyyz[k];

                g_z_0_yyyyy_yyyzz[k] = -g_z_0_yyyy_yyyzz[k] * cd_y[k] + g_z_0_yyyy_yyyyzz[k];

                g_z_0_yyyyy_yyzzz[k] = -g_z_0_yyyy_yyzzz[k] * cd_y[k] + g_z_0_yyyy_yyyzzz[k];

                g_z_0_yyyyy_yzzzz[k] = -g_z_0_yyyy_yzzzz[k] * cd_y[k] + g_z_0_yyyy_yyzzzz[k];

                g_z_0_yyyyy_zzzzz[k] = -g_z_0_yyyy_zzzzz[k] * cd_y[k] + g_z_0_yyyy_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 336);

            auto g_z_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 337);

            auto g_z_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 338);

            auto g_z_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 339);

            auto g_z_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 340);

            auto g_z_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 341);

            auto g_z_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 342);

            auto g_z_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 343);

            auto g_z_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 344);

            auto g_z_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 345);

            auto g_z_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 346);

            auto g_z_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 347);

            auto g_z_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 348);

            auto g_z_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 349);

            auto g_z_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 350);

            auto g_z_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 351);

            auto g_z_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 352);

            auto g_z_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 353);

            auto g_z_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 354);

            auto g_z_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 355);

            auto g_z_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 356);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyz_xxxxx, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxxz, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxzz, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxzzz, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xzzzz, g_z_0_yyyyz_yyyyy, g_z_0_yyyyz_yyyyz, g_z_0_yyyyz_yyyzz, g_z_0_yyyyz_yyzzz, g_z_0_yyyyz_yzzzz, g_z_0_yyyyz_zzzzz, g_z_0_yyyz_xxxxx, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxz, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxzz, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxzzz, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xzzzz, g_z_0_yyyz_yyyyy, g_z_0_yyyz_yyyyyy, g_z_0_yyyz_yyyyyz, g_z_0_yyyz_yyyyz, g_z_0_yyyz_yyyyzz, g_z_0_yyyz_yyyzz, g_z_0_yyyz_yyyzzz, g_z_0_yyyz_yyzzz, g_z_0_yyyz_yyzzzz, g_z_0_yyyz_yzzzz, g_z_0_yyyz_yzzzzz, g_z_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xxxxx[k] = -g_z_0_yyyz_xxxxx[k] * cd_y[k] + g_z_0_yyyz_xxxxxy[k];

                g_z_0_yyyyz_xxxxy[k] = -g_z_0_yyyz_xxxxy[k] * cd_y[k] + g_z_0_yyyz_xxxxyy[k];

                g_z_0_yyyyz_xxxxz[k] = -g_z_0_yyyz_xxxxz[k] * cd_y[k] + g_z_0_yyyz_xxxxyz[k];

                g_z_0_yyyyz_xxxyy[k] = -g_z_0_yyyz_xxxyy[k] * cd_y[k] + g_z_0_yyyz_xxxyyy[k];

                g_z_0_yyyyz_xxxyz[k] = -g_z_0_yyyz_xxxyz[k] * cd_y[k] + g_z_0_yyyz_xxxyyz[k];

                g_z_0_yyyyz_xxxzz[k] = -g_z_0_yyyz_xxxzz[k] * cd_y[k] + g_z_0_yyyz_xxxyzz[k];

                g_z_0_yyyyz_xxyyy[k] = -g_z_0_yyyz_xxyyy[k] * cd_y[k] + g_z_0_yyyz_xxyyyy[k];

                g_z_0_yyyyz_xxyyz[k] = -g_z_0_yyyz_xxyyz[k] * cd_y[k] + g_z_0_yyyz_xxyyyz[k];

                g_z_0_yyyyz_xxyzz[k] = -g_z_0_yyyz_xxyzz[k] * cd_y[k] + g_z_0_yyyz_xxyyzz[k];

                g_z_0_yyyyz_xxzzz[k] = -g_z_0_yyyz_xxzzz[k] * cd_y[k] + g_z_0_yyyz_xxyzzz[k];

                g_z_0_yyyyz_xyyyy[k] = -g_z_0_yyyz_xyyyy[k] * cd_y[k] + g_z_0_yyyz_xyyyyy[k];

                g_z_0_yyyyz_xyyyz[k] = -g_z_0_yyyz_xyyyz[k] * cd_y[k] + g_z_0_yyyz_xyyyyz[k];

                g_z_0_yyyyz_xyyzz[k] = -g_z_0_yyyz_xyyzz[k] * cd_y[k] + g_z_0_yyyz_xyyyzz[k];

                g_z_0_yyyyz_xyzzz[k] = -g_z_0_yyyz_xyzzz[k] * cd_y[k] + g_z_0_yyyz_xyyzzz[k];

                g_z_0_yyyyz_xzzzz[k] = -g_z_0_yyyz_xzzzz[k] * cd_y[k] + g_z_0_yyyz_xyzzzz[k];

                g_z_0_yyyyz_yyyyy[k] = -g_z_0_yyyz_yyyyy[k] * cd_y[k] + g_z_0_yyyz_yyyyyy[k];

                g_z_0_yyyyz_yyyyz[k] = -g_z_0_yyyz_yyyyz[k] * cd_y[k] + g_z_0_yyyz_yyyyyz[k];

                g_z_0_yyyyz_yyyzz[k] = -g_z_0_yyyz_yyyzz[k] * cd_y[k] + g_z_0_yyyz_yyyyzz[k];

                g_z_0_yyyyz_yyzzz[k] = -g_z_0_yyyz_yyzzz[k] * cd_y[k] + g_z_0_yyyz_yyyzzz[k];

                g_z_0_yyyyz_yzzzz[k] = -g_z_0_yyyz_yzzzz[k] * cd_y[k] + g_z_0_yyyz_yyzzzz[k];

                g_z_0_yyyyz_zzzzz[k] = -g_z_0_yyyz_zzzzz[k] * cd_y[k] + g_z_0_yyyz_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 357);

            auto g_z_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 358);

            auto g_z_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 359);

            auto g_z_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 360);

            auto g_z_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 361);

            auto g_z_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 362);

            auto g_z_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 363);

            auto g_z_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 364);

            auto g_z_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 365);

            auto g_z_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 366);

            auto g_z_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 367);

            auto g_z_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 368);

            auto g_z_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 369);

            auto g_z_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 370);

            auto g_z_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 371);

            auto g_z_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 372);

            auto g_z_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 373);

            auto g_z_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 374);

            auto g_z_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 375);

            auto g_z_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 376);

            auto g_z_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 377);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzz_xxxxx, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxxz, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxzz, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxzzz, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xzzzz, g_z_0_yyyzz_yyyyy, g_z_0_yyyzz_yyyyz, g_z_0_yyyzz_yyyzz, g_z_0_yyyzz_yyzzz, g_z_0_yyyzz_yzzzz, g_z_0_yyyzz_zzzzz, g_z_0_yyzz_xxxxx, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxz, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxzz, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxzzz, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xzzzz, g_z_0_yyzz_yyyyy, g_z_0_yyzz_yyyyyy, g_z_0_yyzz_yyyyyz, g_z_0_yyzz_yyyyz, g_z_0_yyzz_yyyyzz, g_z_0_yyzz_yyyzz, g_z_0_yyzz_yyyzzz, g_z_0_yyzz_yyzzz, g_z_0_yyzz_yyzzzz, g_z_0_yyzz_yzzzz, g_z_0_yyzz_yzzzzz, g_z_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xxxxx[k] = -g_z_0_yyzz_xxxxx[k] * cd_y[k] + g_z_0_yyzz_xxxxxy[k];

                g_z_0_yyyzz_xxxxy[k] = -g_z_0_yyzz_xxxxy[k] * cd_y[k] + g_z_0_yyzz_xxxxyy[k];

                g_z_0_yyyzz_xxxxz[k] = -g_z_0_yyzz_xxxxz[k] * cd_y[k] + g_z_0_yyzz_xxxxyz[k];

                g_z_0_yyyzz_xxxyy[k] = -g_z_0_yyzz_xxxyy[k] * cd_y[k] + g_z_0_yyzz_xxxyyy[k];

                g_z_0_yyyzz_xxxyz[k] = -g_z_0_yyzz_xxxyz[k] * cd_y[k] + g_z_0_yyzz_xxxyyz[k];

                g_z_0_yyyzz_xxxzz[k] = -g_z_0_yyzz_xxxzz[k] * cd_y[k] + g_z_0_yyzz_xxxyzz[k];

                g_z_0_yyyzz_xxyyy[k] = -g_z_0_yyzz_xxyyy[k] * cd_y[k] + g_z_0_yyzz_xxyyyy[k];

                g_z_0_yyyzz_xxyyz[k] = -g_z_0_yyzz_xxyyz[k] * cd_y[k] + g_z_0_yyzz_xxyyyz[k];

                g_z_0_yyyzz_xxyzz[k] = -g_z_0_yyzz_xxyzz[k] * cd_y[k] + g_z_0_yyzz_xxyyzz[k];

                g_z_0_yyyzz_xxzzz[k] = -g_z_0_yyzz_xxzzz[k] * cd_y[k] + g_z_0_yyzz_xxyzzz[k];

                g_z_0_yyyzz_xyyyy[k] = -g_z_0_yyzz_xyyyy[k] * cd_y[k] + g_z_0_yyzz_xyyyyy[k];

                g_z_0_yyyzz_xyyyz[k] = -g_z_0_yyzz_xyyyz[k] * cd_y[k] + g_z_0_yyzz_xyyyyz[k];

                g_z_0_yyyzz_xyyzz[k] = -g_z_0_yyzz_xyyzz[k] * cd_y[k] + g_z_0_yyzz_xyyyzz[k];

                g_z_0_yyyzz_xyzzz[k] = -g_z_0_yyzz_xyzzz[k] * cd_y[k] + g_z_0_yyzz_xyyzzz[k];

                g_z_0_yyyzz_xzzzz[k] = -g_z_0_yyzz_xzzzz[k] * cd_y[k] + g_z_0_yyzz_xyzzzz[k];

                g_z_0_yyyzz_yyyyy[k] = -g_z_0_yyzz_yyyyy[k] * cd_y[k] + g_z_0_yyzz_yyyyyy[k];

                g_z_0_yyyzz_yyyyz[k] = -g_z_0_yyzz_yyyyz[k] * cd_y[k] + g_z_0_yyzz_yyyyyz[k];

                g_z_0_yyyzz_yyyzz[k] = -g_z_0_yyzz_yyyzz[k] * cd_y[k] + g_z_0_yyzz_yyyyzz[k];

                g_z_0_yyyzz_yyzzz[k] = -g_z_0_yyzz_yyzzz[k] * cd_y[k] + g_z_0_yyzz_yyyzzz[k];

                g_z_0_yyyzz_yzzzz[k] = -g_z_0_yyzz_yzzzz[k] * cd_y[k] + g_z_0_yyzz_yyzzzz[k];

                g_z_0_yyyzz_zzzzz[k] = -g_z_0_yyzz_zzzzz[k] * cd_y[k] + g_z_0_yyzz_yzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 378);

            auto g_z_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 379);

            auto g_z_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 380);

            auto g_z_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 381);

            auto g_z_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 382);

            auto g_z_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 383);

            auto g_z_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 384);

            auto g_z_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 385);

            auto g_z_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 386);

            auto g_z_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 387);

            auto g_z_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 388);

            auto g_z_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 389);

            auto g_z_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 390);

            auto g_z_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 391);

            auto g_z_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 392);

            auto g_z_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 393);

            auto g_z_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 394);

            auto g_z_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 395);

            auto g_z_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 396);

            auto g_z_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 397);

            auto g_z_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 398);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzz_xxxxx, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxxz, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxzz, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxzzz, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xzzzz, g_z_0_yyzzz_yyyyy, g_z_0_yyzzz_yyyyz, g_z_0_yyzzz_yyyzz, g_z_0_yyzzz_yyzzz, g_z_0_yyzzz_yzzzz, g_z_0_yyzzz_zzzzz, g_z_0_yzzz_xxxxx, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxz, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxzz, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxzzz, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xzzzz, g_z_0_yzzz_yyyyy, g_z_0_yzzz_yyyyyy, g_z_0_yzzz_yyyyyz, g_z_0_yzzz_yyyyz, g_z_0_yzzz_yyyyzz, g_z_0_yzzz_yyyzz, g_z_0_yzzz_yyyzzz, g_z_0_yzzz_yyzzz, g_z_0_yzzz_yyzzzz, g_z_0_yzzz_yzzzz, g_z_0_yzzz_yzzzzz, g_z_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xxxxx[k] = -g_z_0_yzzz_xxxxx[k] * cd_y[k] + g_z_0_yzzz_xxxxxy[k];

                g_z_0_yyzzz_xxxxy[k] = -g_z_0_yzzz_xxxxy[k] * cd_y[k] + g_z_0_yzzz_xxxxyy[k];

                g_z_0_yyzzz_xxxxz[k] = -g_z_0_yzzz_xxxxz[k] * cd_y[k] + g_z_0_yzzz_xxxxyz[k];

                g_z_0_yyzzz_xxxyy[k] = -g_z_0_yzzz_xxxyy[k] * cd_y[k] + g_z_0_yzzz_xxxyyy[k];

                g_z_0_yyzzz_xxxyz[k] = -g_z_0_yzzz_xxxyz[k] * cd_y[k] + g_z_0_yzzz_xxxyyz[k];

                g_z_0_yyzzz_xxxzz[k] = -g_z_0_yzzz_xxxzz[k] * cd_y[k] + g_z_0_yzzz_xxxyzz[k];

                g_z_0_yyzzz_xxyyy[k] = -g_z_0_yzzz_xxyyy[k] * cd_y[k] + g_z_0_yzzz_xxyyyy[k];

                g_z_0_yyzzz_xxyyz[k] = -g_z_0_yzzz_xxyyz[k] * cd_y[k] + g_z_0_yzzz_xxyyyz[k];

                g_z_0_yyzzz_xxyzz[k] = -g_z_0_yzzz_xxyzz[k] * cd_y[k] + g_z_0_yzzz_xxyyzz[k];

                g_z_0_yyzzz_xxzzz[k] = -g_z_0_yzzz_xxzzz[k] * cd_y[k] + g_z_0_yzzz_xxyzzz[k];

                g_z_0_yyzzz_xyyyy[k] = -g_z_0_yzzz_xyyyy[k] * cd_y[k] + g_z_0_yzzz_xyyyyy[k];

                g_z_0_yyzzz_xyyyz[k] = -g_z_0_yzzz_xyyyz[k] * cd_y[k] + g_z_0_yzzz_xyyyyz[k];

                g_z_0_yyzzz_xyyzz[k] = -g_z_0_yzzz_xyyzz[k] * cd_y[k] + g_z_0_yzzz_xyyyzz[k];

                g_z_0_yyzzz_xyzzz[k] = -g_z_0_yzzz_xyzzz[k] * cd_y[k] + g_z_0_yzzz_xyyzzz[k];

                g_z_0_yyzzz_xzzzz[k] = -g_z_0_yzzz_xzzzz[k] * cd_y[k] + g_z_0_yzzz_xyzzzz[k];

                g_z_0_yyzzz_yyyyy[k] = -g_z_0_yzzz_yyyyy[k] * cd_y[k] + g_z_0_yzzz_yyyyyy[k];

                g_z_0_yyzzz_yyyyz[k] = -g_z_0_yzzz_yyyyz[k] * cd_y[k] + g_z_0_yzzz_yyyyyz[k];

                g_z_0_yyzzz_yyyzz[k] = -g_z_0_yzzz_yyyzz[k] * cd_y[k] + g_z_0_yzzz_yyyyzz[k];

                g_z_0_yyzzz_yyzzz[k] = -g_z_0_yzzz_yyzzz[k] * cd_y[k] + g_z_0_yzzz_yyyzzz[k];

                g_z_0_yyzzz_yzzzz[k] = -g_z_0_yzzz_yzzzz[k] * cd_y[k] + g_z_0_yzzz_yyzzzz[k];

                g_z_0_yyzzz_zzzzz[k] = -g_z_0_yzzz_zzzzz[k] * cd_y[k] + g_z_0_yzzz_yzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 399);

            auto g_z_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 400);

            auto g_z_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 401);

            auto g_z_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 402);

            auto g_z_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 403);

            auto g_z_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 404);

            auto g_z_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 405);

            auto g_z_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 406);

            auto g_z_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 407);

            auto g_z_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 408);

            auto g_z_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 409);

            auto g_z_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 410);

            auto g_z_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 411);

            auto g_z_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 412);

            auto g_z_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 413);

            auto g_z_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 414);

            auto g_z_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 415);

            auto g_z_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 416);

            auto g_z_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 417);

            auto g_z_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 418);

            auto g_z_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzz_xxxxx, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxxz, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxzz, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxzzz, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xzzzz, g_z_0_yzzzz_yyyyy, g_z_0_yzzzz_yyyyz, g_z_0_yzzzz_yyyzz, g_z_0_yzzzz_yyzzz, g_z_0_yzzzz_yzzzz, g_z_0_yzzzz_zzzzz, g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xxxxx[k] = -g_z_0_zzzz_xxxxx[k] * cd_y[k] + g_z_0_zzzz_xxxxxy[k];

                g_z_0_yzzzz_xxxxy[k] = -g_z_0_zzzz_xxxxy[k] * cd_y[k] + g_z_0_zzzz_xxxxyy[k];

                g_z_0_yzzzz_xxxxz[k] = -g_z_0_zzzz_xxxxz[k] * cd_y[k] + g_z_0_zzzz_xxxxyz[k];

                g_z_0_yzzzz_xxxyy[k] = -g_z_0_zzzz_xxxyy[k] * cd_y[k] + g_z_0_zzzz_xxxyyy[k];

                g_z_0_yzzzz_xxxyz[k] = -g_z_0_zzzz_xxxyz[k] * cd_y[k] + g_z_0_zzzz_xxxyyz[k];

                g_z_0_yzzzz_xxxzz[k] = -g_z_0_zzzz_xxxzz[k] * cd_y[k] + g_z_0_zzzz_xxxyzz[k];

                g_z_0_yzzzz_xxyyy[k] = -g_z_0_zzzz_xxyyy[k] * cd_y[k] + g_z_0_zzzz_xxyyyy[k];

                g_z_0_yzzzz_xxyyz[k] = -g_z_0_zzzz_xxyyz[k] * cd_y[k] + g_z_0_zzzz_xxyyyz[k];

                g_z_0_yzzzz_xxyzz[k] = -g_z_0_zzzz_xxyzz[k] * cd_y[k] + g_z_0_zzzz_xxyyzz[k];

                g_z_0_yzzzz_xxzzz[k] = -g_z_0_zzzz_xxzzz[k] * cd_y[k] + g_z_0_zzzz_xxyzzz[k];

                g_z_0_yzzzz_xyyyy[k] = -g_z_0_zzzz_xyyyy[k] * cd_y[k] + g_z_0_zzzz_xyyyyy[k];

                g_z_0_yzzzz_xyyyz[k] = -g_z_0_zzzz_xyyyz[k] * cd_y[k] + g_z_0_zzzz_xyyyyz[k];

                g_z_0_yzzzz_xyyzz[k] = -g_z_0_zzzz_xyyzz[k] * cd_y[k] + g_z_0_zzzz_xyyyzz[k];

                g_z_0_yzzzz_xyzzz[k] = -g_z_0_zzzz_xyzzz[k] * cd_y[k] + g_z_0_zzzz_xyyzzz[k];

                g_z_0_yzzzz_xzzzz[k] = -g_z_0_zzzz_xzzzz[k] * cd_y[k] + g_z_0_zzzz_xyzzzz[k];

                g_z_0_yzzzz_yyyyy[k] = -g_z_0_zzzz_yyyyy[k] * cd_y[k] + g_z_0_zzzz_yyyyyy[k];

                g_z_0_yzzzz_yyyyz[k] = -g_z_0_zzzz_yyyyz[k] * cd_y[k] + g_z_0_zzzz_yyyyyz[k];

                g_z_0_yzzzz_yyyzz[k] = -g_z_0_zzzz_yyyzz[k] * cd_y[k] + g_z_0_zzzz_yyyyzz[k];

                g_z_0_yzzzz_yyzzz[k] = -g_z_0_zzzz_yyzzz[k] * cd_y[k] + g_z_0_zzzz_yyyzzz[k];

                g_z_0_yzzzz_yzzzz[k] = -g_z_0_zzzz_yzzzz[k] * cd_y[k] + g_z_0_zzzz_yyzzzz[k];

                g_z_0_yzzzz_zzzzz[k] = -g_z_0_zzzz_zzzzz[k] * cd_y[k] + g_z_0_zzzz_yzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 420);

            auto g_z_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 421);

            auto g_z_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 422);

            auto g_z_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 423);

            auto g_z_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 424);

            auto g_z_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 425);

            auto g_z_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 426);

            auto g_z_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 427);

            auto g_z_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 428);

            auto g_z_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 429);

            auto g_z_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 430);

            auto g_z_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 431);

            auto g_z_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 432);

            auto g_z_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 433);

            auto g_z_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 434);

            auto g_z_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 435);

            auto g_z_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 436);

            auto g_z_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 437);

            auto g_z_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 438);

            auto g_z_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 439);

            auto g_z_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 882 * acomps * bcomps + 440);

            #pragma omp simd aligned(cd_z, g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_zzzzz, g_z_0_zzzz_zzzzzz, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_zzzzz, g_zzzz_xxxxx, g_zzzz_xxxxy, g_zzzz_xxxxz, g_zzzz_xxxyy, g_zzzz_xxxyz, g_zzzz_xxxzz, g_zzzz_xxyyy, g_zzzz_xxyyz, g_zzzz_xxyzz, g_zzzz_xxzzz, g_zzzz_xyyyy, g_zzzz_xyyyz, g_zzzz_xyyzz, g_zzzz_xyzzz, g_zzzz_xzzzz, g_zzzz_yyyyy, g_zzzz_yyyyz, g_zzzz_yyyzz, g_zzzz_yyzzz, g_zzzz_yzzzz, g_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xxxxx[k] = -g_zzzz_xxxxx[k] - g_z_0_zzzz_xxxxx[k] * cd_z[k] + g_z_0_zzzz_xxxxxz[k];

                g_z_0_zzzzz_xxxxy[k] = -g_zzzz_xxxxy[k] - g_z_0_zzzz_xxxxy[k] * cd_z[k] + g_z_0_zzzz_xxxxyz[k];

                g_z_0_zzzzz_xxxxz[k] = -g_zzzz_xxxxz[k] - g_z_0_zzzz_xxxxz[k] * cd_z[k] + g_z_0_zzzz_xxxxzz[k];

                g_z_0_zzzzz_xxxyy[k] = -g_zzzz_xxxyy[k] - g_z_0_zzzz_xxxyy[k] * cd_z[k] + g_z_0_zzzz_xxxyyz[k];

                g_z_0_zzzzz_xxxyz[k] = -g_zzzz_xxxyz[k] - g_z_0_zzzz_xxxyz[k] * cd_z[k] + g_z_0_zzzz_xxxyzz[k];

                g_z_0_zzzzz_xxxzz[k] = -g_zzzz_xxxzz[k] - g_z_0_zzzz_xxxzz[k] * cd_z[k] + g_z_0_zzzz_xxxzzz[k];

                g_z_0_zzzzz_xxyyy[k] = -g_zzzz_xxyyy[k] - g_z_0_zzzz_xxyyy[k] * cd_z[k] + g_z_0_zzzz_xxyyyz[k];

                g_z_0_zzzzz_xxyyz[k] = -g_zzzz_xxyyz[k] - g_z_0_zzzz_xxyyz[k] * cd_z[k] + g_z_0_zzzz_xxyyzz[k];

                g_z_0_zzzzz_xxyzz[k] = -g_zzzz_xxyzz[k] - g_z_0_zzzz_xxyzz[k] * cd_z[k] + g_z_0_zzzz_xxyzzz[k];

                g_z_0_zzzzz_xxzzz[k] = -g_zzzz_xxzzz[k] - g_z_0_zzzz_xxzzz[k] * cd_z[k] + g_z_0_zzzz_xxzzzz[k];

                g_z_0_zzzzz_xyyyy[k] = -g_zzzz_xyyyy[k] - g_z_0_zzzz_xyyyy[k] * cd_z[k] + g_z_0_zzzz_xyyyyz[k];

                g_z_0_zzzzz_xyyyz[k] = -g_zzzz_xyyyz[k] - g_z_0_zzzz_xyyyz[k] * cd_z[k] + g_z_0_zzzz_xyyyzz[k];

                g_z_0_zzzzz_xyyzz[k] = -g_zzzz_xyyzz[k] - g_z_0_zzzz_xyyzz[k] * cd_z[k] + g_z_0_zzzz_xyyzzz[k];

                g_z_0_zzzzz_xyzzz[k] = -g_zzzz_xyzzz[k] - g_z_0_zzzz_xyzzz[k] * cd_z[k] + g_z_0_zzzz_xyzzzz[k];

                g_z_0_zzzzz_xzzzz[k] = -g_zzzz_xzzzz[k] - g_z_0_zzzz_xzzzz[k] * cd_z[k] + g_z_0_zzzz_xzzzzz[k];

                g_z_0_zzzzz_yyyyy[k] = -g_zzzz_yyyyy[k] - g_z_0_zzzz_yyyyy[k] * cd_z[k] + g_z_0_zzzz_yyyyyz[k];

                g_z_0_zzzzz_yyyyz[k] = -g_zzzz_yyyyz[k] - g_z_0_zzzz_yyyyz[k] * cd_z[k] + g_z_0_zzzz_yyyyzz[k];

                g_z_0_zzzzz_yyyzz[k] = -g_zzzz_yyyzz[k] - g_z_0_zzzz_yyyzz[k] * cd_z[k] + g_z_0_zzzz_yyyzzz[k];

                g_z_0_zzzzz_yyzzz[k] = -g_zzzz_yyzzz[k] - g_z_0_zzzz_yyzzz[k] * cd_z[k] + g_z_0_zzzz_yyzzzz[k];

                g_z_0_zzzzz_yzzzz[k] = -g_zzzz_yzzzz[k] - g_z_0_zzzz_yzzzz[k] * cd_z[k] + g_z_0_zzzz_yzzzzz[k];

                g_z_0_zzzzz_zzzzz[k] = -g_zzzz_zzzzz[k] - g_z_0_zzzz_zzzzz[k] * cd_z[k] + g_z_0_zzzz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

