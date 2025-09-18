#include "ElectronRepulsionGeom0010ContrRecXXIH.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxih(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxih,
                                            const size_t idx_xxhh,
                                            const size_t idx_geom_10_xxhh,
                                            const size_t idx_geom_10_xxhi,
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
            /// Set up components of auxilary buffer : SSHH

            const auto hh_off = idx_xxhh + (i * bcomps + j) * 441;

            auto g_xxxxx_xxxxx = cbuffer.data(hh_off + 0);

            auto g_xxxxx_xxxxy = cbuffer.data(hh_off + 1);

            auto g_xxxxx_xxxxz = cbuffer.data(hh_off + 2);

            auto g_xxxxx_xxxyy = cbuffer.data(hh_off + 3);

            auto g_xxxxx_xxxyz = cbuffer.data(hh_off + 4);

            auto g_xxxxx_xxxzz = cbuffer.data(hh_off + 5);

            auto g_xxxxx_xxyyy = cbuffer.data(hh_off + 6);

            auto g_xxxxx_xxyyz = cbuffer.data(hh_off + 7);

            auto g_xxxxx_xxyzz = cbuffer.data(hh_off + 8);

            auto g_xxxxx_xxzzz = cbuffer.data(hh_off + 9);

            auto g_xxxxx_xyyyy = cbuffer.data(hh_off + 10);

            auto g_xxxxx_xyyyz = cbuffer.data(hh_off + 11);

            auto g_xxxxx_xyyzz = cbuffer.data(hh_off + 12);

            auto g_xxxxx_xyzzz = cbuffer.data(hh_off + 13);

            auto g_xxxxx_xzzzz = cbuffer.data(hh_off + 14);

            auto g_xxxxx_yyyyy = cbuffer.data(hh_off + 15);

            auto g_xxxxx_yyyyz = cbuffer.data(hh_off + 16);

            auto g_xxxxx_yyyzz = cbuffer.data(hh_off + 17);

            auto g_xxxxx_yyzzz = cbuffer.data(hh_off + 18);

            auto g_xxxxx_yzzzz = cbuffer.data(hh_off + 19);

            auto g_xxxxx_zzzzz = cbuffer.data(hh_off + 20);

            auto g_xxxxy_xxxxx = cbuffer.data(hh_off + 21);

            auto g_xxxxy_xxxxy = cbuffer.data(hh_off + 22);

            auto g_xxxxy_xxxxz = cbuffer.data(hh_off + 23);

            auto g_xxxxy_xxxyy = cbuffer.data(hh_off + 24);

            auto g_xxxxy_xxxyz = cbuffer.data(hh_off + 25);

            auto g_xxxxy_xxxzz = cbuffer.data(hh_off + 26);

            auto g_xxxxy_xxyyy = cbuffer.data(hh_off + 27);

            auto g_xxxxy_xxyyz = cbuffer.data(hh_off + 28);

            auto g_xxxxy_xxyzz = cbuffer.data(hh_off + 29);

            auto g_xxxxy_xxzzz = cbuffer.data(hh_off + 30);

            auto g_xxxxy_xyyyy = cbuffer.data(hh_off + 31);

            auto g_xxxxy_xyyyz = cbuffer.data(hh_off + 32);

            auto g_xxxxy_xyyzz = cbuffer.data(hh_off + 33);

            auto g_xxxxy_xyzzz = cbuffer.data(hh_off + 34);

            auto g_xxxxy_xzzzz = cbuffer.data(hh_off + 35);

            auto g_xxxxy_yyyyy = cbuffer.data(hh_off + 36);

            auto g_xxxxy_yyyyz = cbuffer.data(hh_off + 37);

            auto g_xxxxy_yyyzz = cbuffer.data(hh_off + 38);

            auto g_xxxxy_yyzzz = cbuffer.data(hh_off + 39);

            auto g_xxxxy_yzzzz = cbuffer.data(hh_off + 40);

            auto g_xxxxy_zzzzz = cbuffer.data(hh_off + 41);

            auto g_xxxxz_xxxxx = cbuffer.data(hh_off + 42);

            auto g_xxxxz_xxxxy = cbuffer.data(hh_off + 43);

            auto g_xxxxz_xxxxz = cbuffer.data(hh_off + 44);

            auto g_xxxxz_xxxyy = cbuffer.data(hh_off + 45);

            auto g_xxxxz_xxxyz = cbuffer.data(hh_off + 46);

            auto g_xxxxz_xxxzz = cbuffer.data(hh_off + 47);

            auto g_xxxxz_xxyyy = cbuffer.data(hh_off + 48);

            auto g_xxxxz_xxyyz = cbuffer.data(hh_off + 49);

            auto g_xxxxz_xxyzz = cbuffer.data(hh_off + 50);

            auto g_xxxxz_xxzzz = cbuffer.data(hh_off + 51);

            auto g_xxxxz_xyyyy = cbuffer.data(hh_off + 52);

            auto g_xxxxz_xyyyz = cbuffer.data(hh_off + 53);

            auto g_xxxxz_xyyzz = cbuffer.data(hh_off + 54);

            auto g_xxxxz_xyzzz = cbuffer.data(hh_off + 55);

            auto g_xxxxz_xzzzz = cbuffer.data(hh_off + 56);

            auto g_xxxxz_yyyyy = cbuffer.data(hh_off + 57);

            auto g_xxxxz_yyyyz = cbuffer.data(hh_off + 58);

            auto g_xxxxz_yyyzz = cbuffer.data(hh_off + 59);

            auto g_xxxxz_yyzzz = cbuffer.data(hh_off + 60);

            auto g_xxxxz_yzzzz = cbuffer.data(hh_off + 61);

            auto g_xxxxz_zzzzz = cbuffer.data(hh_off + 62);

            auto g_xxxyy_xxxxx = cbuffer.data(hh_off + 63);

            auto g_xxxyy_xxxxy = cbuffer.data(hh_off + 64);

            auto g_xxxyy_xxxxz = cbuffer.data(hh_off + 65);

            auto g_xxxyy_xxxyy = cbuffer.data(hh_off + 66);

            auto g_xxxyy_xxxyz = cbuffer.data(hh_off + 67);

            auto g_xxxyy_xxxzz = cbuffer.data(hh_off + 68);

            auto g_xxxyy_xxyyy = cbuffer.data(hh_off + 69);

            auto g_xxxyy_xxyyz = cbuffer.data(hh_off + 70);

            auto g_xxxyy_xxyzz = cbuffer.data(hh_off + 71);

            auto g_xxxyy_xxzzz = cbuffer.data(hh_off + 72);

            auto g_xxxyy_xyyyy = cbuffer.data(hh_off + 73);

            auto g_xxxyy_xyyyz = cbuffer.data(hh_off + 74);

            auto g_xxxyy_xyyzz = cbuffer.data(hh_off + 75);

            auto g_xxxyy_xyzzz = cbuffer.data(hh_off + 76);

            auto g_xxxyy_xzzzz = cbuffer.data(hh_off + 77);

            auto g_xxxyy_yyyyy = cbuffer.data(hh_off + 78);

            auto g_xxxyy_yyyyz = cbuffer.data(hh_off + 79);

            auto g_xxxyy_yyyzz = cbuffer.data(hh_off + 80);

            auto g_xxxyy_yyzzz = cbuffer.data(hh_off + 81);

            auto g_xxxyy_yzzzz = cbuffer.data(hh_off + 82);

            auto g_xxxyy_zzzzz = cbuffer.data(hh_off + 83);

            auto g_xxxyz_xxxxx = cbuffer.data(hh_off + 84);

            auto g_xxxyz_xxxxy = cbuffer.data(hh_off + 85);

            auto g_xxxyz_xxxxz = cbuffer.data(hh_off + 86);

            auto g_xxxyz_xxxyy = cbuffer.data(hh_off + 87);

            auto g_xxxyz_xxxyz = cbuffer.data(hh_off + 88);

            auto g_xxxyz_xxxzz = cbuffer.data(hh_off + 89);

            auto g_xxxyz_xxyyy = cbuffer.data(hh_off + 90);

            auto g_xxxyz_xxyyz = cbuffer.data(hh_off + 91);

            auto g_xxxyz_xxyzz = cbuffer.data(hh_off + 92);

            auto g_xxxyz_xxzzz = cbuffer.data(hh_off + 93);

            auto g_xxxyz_xyyyy = cbuffer.data(hh_off + 94);

            auto g_xxxyz_xyyyz = cbuffer.data(hh_off + 95);

            auto g_xxxyz_xyyzz = cbuffer.data(hh_off + 96);

            auto g_xxxyz_xyzzz = cbuffer.data(hh_off + 97);

            auto g_xxxyz_xzzzz = cbuffer.data(hh_off + 98);

            auto g_xxxyz_yyyyy = cbuffer.data(hh_off + 99);

            auto g_xxxyz_yyyyz = cbuffer.data(hh_off + 100);

            auto g_xxxyz_yyyzz = cbuffer.data(hh_off + 101);

            auto g_xxxyz_yyzzz = cbuffer.data(hh_off + 102);

            auto g_xxxyz_yzzzz = cbuffer.data(hh_off + 103);

            auto g_xxxyz_zzzzz = cbuffer.data(hh_off + 104);

            auto g_xxxzz_xxxxx = cbuffer.data(hh_off + 105);

            auto g_xxxzz_xxxxy = cbuffer.data(hh_off + 106);

            auto g_xxxzz_xxxxz = cbuffer.data(hh_off + 107);

            auto g_xxxzz_xxxyy = cbuffer.data(hh_off + 108);

            auto g_xxxzz_xxxyz = cbuffer.data(hh_off + 109);

            auto g_xxxzz_xxxzz = cbuffer.data(hh_off + 110);

            auto g_xxxzz_xxyyy = cbuffer.data(hh_off + 111);

            auto g_xxxzz_xxyyz = cbuffer.data(hh_off + 112);

            auto g_xxxzz_xxyzz = cbuffer.data(hh_off + 113);

            auto g_xxxzz_xxzzz = cbuffer.data(hh_off + 114);

            auto g_xxxzz_xyyyy = cbuffer.data(hh_off + 115);

            auto g_xxxzz_xyyyz = cbuffer.data(hh_off + 116);

            auto g_xxxzz_xyyzz = cbuffer.data(hh_off + 117);

            auto g_xxxzz_xyzzz = cbuffer.data(hh_off + 118);

            auto g_xxxzz_xzzzz = cbuffer.data(hh_off + 119);

            auto g_xxxzz_yyyyy = cbuffer.data(hh_off + 120);

            auto g_xxxzz_yyyyz = cbuffer.data(hh_off + 121);

            auto g_xxxzz_yyyzz = cbuffer.data(hh_off + 122);

            auto g_xxxzz_yyzzz = cbuffer.data(hh_off + 123);

            auto g_xxxzz_yzzzz = cbuffer.data(hh_off + 124);

            auto g_xxxzz_zzzzz = cbuffer.data(hh_off + 125);

            auto g_xxyyy_xxxxx = cbuffer.data(hh_off + 126);

            auto g_xxyyy_xxxxy = cbuffer.data(hh_off + 127);

            auto g_xxyyy_xxxxz = cbuffer.data(hh_off + 128);

            auto g_xxyyy_xxxyy = cbuffer.data(hh_off + 129);

            auto g_xxyyy_xxxyz = cbuffer.data(hh_off + 130);

            auto g_xxyyy_xxxzz = cbuffer.data(hh_off + 131);

            auto g_xxyyy_xxyyy = cbuffer.data(hh_off + 132);

            auto g_xxyyy_xxyyz = cbuffer.data(hh_off + 133);

            auto g_xxyyy_xxyzz = cbuffer.data(hh_off + 134);

            auto g_xxyyy_xxzzz = cbuffer.data(hh_off + 135);

            auto g_xxyyy_xyyyy = cbuffer.data(hh_off + 136);

            auto g_xxyyy_xyyyz = cbuffer.data(hh_off + 137);

            auto g_xxyyy_xyyzz = cbuffer.data(hh_off + 138);

            auto g_xxyyy_xyzzz = cbuffer.data(hh_off + 139);

            auto g_xxyyy_xzzzz = cbuffer.data(hh_off + 140);

            auto g_xxyyy_yyyyy = cbuffer.data(hh_off + 141);

            auto g_xxyyy_yyyyz = cbuffer.data(hh_off + 142);

            auto g_xxyyy_yyyzz = cbuffer.data(hh_off + 143);

            auto g_xxyyy_yyzzz = cbuffer.data(hh_off + 144);

            auto g_xxyyy_yzzzz = cbuffer.data(hh_off + 145);

            auto g_xxyyy_zzzzz = cbuffer.data(hh_off + 146);

            auto g_xxyyz_xxxxx = cbuffer.data(hh_off + 147);

            auto g_xxyyz_xxxxy = cbuffer.data(hh_off + 148);

            auto g_xxyyz_xxxxz = cbuffer.data(hh_off + 149);

            auto g_xxyyz_xxxyy = cbuffer.data(hh_off + 150);

            auto g_xxyyz_xxxyz = cbuffer.data(hh_off + 151);

            auto g_xxyyz_xxxzz = cbuffer.data(hh_off + 152);

            auto g_xxyyz_xxyyy = cbuffer.data(hh_off + 153);

            auto g_xxyyz_xxyyz = cbuffer.data(hh_off + 154);

            auto g_xxyyz_xxyzz = cbuffer.data(hh_off + 155);

            auto g_xxyyz_xxzzz = cbuffer.data(hh_off + 156);

            auto g_xxyyz_xyyyy = cbuffer.data(hh_off + 157);

            auto g_xxyyz_xyyyz = cbuffer.data(hh_off + 158);

            auto g_xxyyz_xyyzz = cbuffer.data(hh_off + 159);

            auto g_xxyyz_xyzzz = cbuffer.data(hh_off + 160);

            auto g_xxyyz_xzzzz = cbuffer.data(hh_off + 161);

            auto g_xxyyz_yyyyy = cbuffer.data(hh_off + 162);

            auto g_xxyyz_yyyyz = cbuffer.data(hh_off + 163);

            auto g_xxyyz_yyyzz = cbuffer.data(hh_off + 164);

            auto g_xxyyz_yyzzz = cbuffer.data(hh_off + 165);

            auto g_xxyyz_yzzzz = cbuffer.data(hh_off + 166);

            auto g_xxyyz_zzzzz = cbuffer.data(hh_off + 167);

            auto g_xxyzz_xxxxx = cbuffer.data(hh_off + 168);

            auto g_xxyzz_xxxxy = cbuffer.data(hh_off + 169);

            auto g_xxyzz_xxxxz = cbuffer.data(hh_off + 170);

            auto g_xxyzz_xxxyy = cbuffer.data(hh_off + 171);

            auto g_xxyzz_xxxyz = cbuffer.data(hh_off + 172);

            auto g_xxyzz_xxxzz = cbuffer.data(hh_off + 173);

            auto g_xxyzz_xxyyy = cbuffer.data(hh_off + 174);

            auto g_xxyzz_xxyyz = cbuffer.data(hh_off + 175);

            auto g_xxyzz_xxyzz = cbuffer.data(hh_off + 176);

            auto g_xxyzz_xxzzz = cbuffer.data(hh_off + 177);

            auto g_xxyzz_xyyyy = cbuffer.data(hh_off + 178);

            auto g_xxyzz_xyyyz = cbuffer.data(hh_off + 179);

            auto g_xxyzz_xyyzz = cbuffer.data(hh_off + 180);

            auto g_xxyzz_xyzzz = cbuffer.data(hh_off + 181);

            auto g_xxyzz_xzzzz = cbuffer.data(hh_off + 182);

            auto g_xxyzz_yyyyy = cbuffer.data(hh_off + 183);

            auto g_xxyzz_yyyyz = cbuffer.data(hh_off + 184);

            auto g_xxyzz_yyyzz = cbuffer.data(hh_off + 185);

            auto g_xxyzz_yyzzz = cbuffer.data(hh_off + 186);

            auto g_xxyzz_yzzzz = cbuffer.data(hh_off + 187);

            auto g_xxyzz_zzzzz = cbuffer.data(hh_off + 188);

            auto g_xxzzz_xxxxx = cbuffer.data(hh_off + 189);

            auto g_xxzzz_xxxxy = cbuffer.data(hh_off + 190);

            auto g_xxzzz_xxxxz = cbuffer.data(hh_off + 191);

            auto g_xxzzz_xxxyy = cbuffer.data(hh_off + 192);

            auto g_xxzzz_xxxyz = cbuffer.data(hh_off + 193);

            auto g_xxzzz_xxxzz = cbuffer.data(hh_off + 194);

            auto g_xxzzz_xxyyy = cbuffer.data(hh_off + 195);

            auto g_xxzzz_xxyyz = cbuffer.data(hh_off + 196);

            auto g_xxzzz_xxyzz = cbuffer.data(hh_off + 197);

            auto g_xxzzz_xxzzz = cbuffer.data(hh_off + 198);

            auto g_xxzzz_xyyyy = cbuffer.data(hh_off + 199);

            auto g_xxzzz_xyyyz = cbuffer.data(hh_off + 200);

            auto g_xxzzz_xyyzz = cbuffer.data(hh_off + 201);

            auto g_xxzzz_xyzzz = cbuffer.data(hh_off + 202);

            auto g_xxzzz_xzzzz = cbuffer.data(hh_off + 203);

            auto g_xxzzz_yyyyy = cbuffer.data(hh_off + 204);

            auto g_xxzzz_yyyyz = cbuffer.data(hh_off + 205);

            auto g_xxzzz_yyyzz = cbuffer.data(hh_off + 206);

            auto g_xxzzz_yyzzz = cbuffer.data(hh_off + 207);

            auto g_xxzzz_yzzzz = cbuffer.data(hh_off + 208);

            auto g_xxzzz_zzzzz = cbuffer.data(hh_off + 209);

            auto g_xyyyy_xxxxx = cbuffer.data(hh_off + 210);

            auto g_xyyyy_xxxxy = cbuffer.data(hh_off + 211);

            auto g_xyyyy_xxxxz = cbuffer.data(hh_off + 212);

            auto g_xyyyy_xxxyy = cbuffer.data(hh_off + 213);

            auto g_xyyyy_xxxyz = cbuffer.data(hh_off + 214);

            auto g_xyyyy_xxxzz = cbuffer.data(hh_off + 215);

            auto g_xyyyy_xxyyy = cbuffer.data(hh_off + 216);

            auto g_xyyyy_xxyyz = cbuffer.data(hh_off + 217);

            auto g_xyyyy_xxyzz = cbuffer.data(hh_off + 218);

            auto g_xyyyy_xxzzz = cbuffer.data(hh_off + 219);

            auto g_xyyyy_xyyyy = cbuffer.data(hh_off + 220);

            auto g_xyyyy_xyyyz = cbuffer.data(hh_off + 221);

            auto g_xyyyy_xyyzz = cbuffer.data(hh_off + 222);

            auto g_xyyyy_xyzzz = cbuffer.data(hh_off + 223);

            auto g_xyyyy_xzzzz = cbuffer.data(hh_off + 224);

            auto g_xyyyy_yyyyy = cbuffer.data(hh_off + 225);

            auto g_xyyyy_yyyyz = cbuffer.data(hh_off + 226);

            auto g_xyyyy_yyyzz = cbuffer.data(hh_off + 227);

            auto g_xyyyy_yyzzz = cbuffer.data(hh_off + 228);

            auto g_xyyyy_yzzzz = cbuffer.data(hh_off + 229);

            auto g_xyyyy_zzzzz = cbuffer.data(hh_off + 230);

            auto g_xyyyz_xxxxx = cbuffer.data(hh_off + 231);

            auto g_xyyyz_xxxxy = cbuffer.data(hh_off + 232);

            auto g_xyyyz_xxxxz = cbuffer.data(hh_off + 233);

            auto g_xyyyz_xxxyy = cbuffer.data(hh_off + 234);

            auto g_xyyyz_xxxyz = cbuffer.data(hh_off + 235);

            auto g_xyyyz_xxxzz = cbuffer.data(hh_off + 236);

            auto g_xyyyz_xxyyy = cbuffer.data(hh_off + 237);

            auto g_xyyyz_xxyyz = cbuffer.data(hh_off + 238);

            auto g_xyyyz_xxyzz = cbuffer.data(hh_off + 239);

            auto g_xyyyz_xxzzz = cbuffer.data(hh_off + 240);

            auto g_xyyyz_xyyyy = cbuffer.data(hh_off + 241);

            auto g_xyyyz_xyyyz = cbuffer.data(hh_off + 242);

            auto g_xyyyz_xyyzz = cbuffer.data(hh_off + 243);

            auto g_xyyyz_xyzzz = cbuffer.data(hh_off + 244);

            auto g_xyyyz_xzzzz = cbuffer.data(hh_off + 245);

            auto g_xyyyz_yyyyy = cbuffer.data(hh_off + 246);

            auto g_xyyyz_yyyyz = cbuffer.data(hh_off + 247);

            auto g_xyyyz_yyyzz = cbuffer.data(hh_off + 248);

            auto g_xyyyz_yyzzz = cbuffer.data(hh_off + 249);

            auto g_xyyyz_yzzzz = cbuffer.data(hh_off + 250);

            auto g_xyyyz_zzzzz = cbuffer.data(hh_off + 251);

            auto g_xyyzz_xxxxx = cbuffer.data(hh_off + 252);

            auto g_xyyzz_xxxxy = cbuffer.data(hh_off + 253);

            auto g_xyyzz_xxxxz = cbuffer.data(hh_off + 254);

            auto g_xyyzz_xxxyy = cbuffer.data(hh_off + 255);

            auto g_xyyzz_xxxyz = cbuffer.data(hh_off + 256);

            auto g_xyyzz_xxxzz = cbuffer.data(hh_off + 257);

            auto g_xyyzz_xxyyy = cbuffer.data(hh_off + 258);

            auto g_xyyzz_xxyyz = cbuffer.data(hh_off + 259);

            auto g_xyyzz_xxyzz = cbuffer.data(hh_off + 260);

            auto g_xyyzz_xxzzz = cbuffer.data(hh_off + 261);

            auto g_xyyzz_xyyyy = cbuffer.data(hh_off + 262);

            auto g_xyyzz_xyyyz = cbuffer.data(hh_off + 263);

            auto g_xyyzz_xyyzz = cbuffer.data(hh_off + 264);

            auto g_xyyzz_xyzzz = cbuffer.data(hh_off + 265);

            auto g_xyyzz_xzzzz = cbuffer.data(hh_off + 266);

            auto g_xyyzz_yyyyy = cbuffer.data(hh_off + 267);

            auto g_xyyzz_yyyyz = cbuffer.data(hh_off + 268);

            auto g_xyyzz_yyyzz = cbuffer.data(hh_off + 269);

            auto g_xyyzz_yyzzz = cbuffer.data(hh_off + 270);

            auto g_xyyzz_yzzzz = cbuffer.data(hh_off + 271);

            auto g_xyyzz_zzzzz = cbuffer.data(hh_off + 272);

            auto g_xyzzz_xxxxx = cbuffer.data(hh_off + 273);

            auto g_xyzzz_xxxxy = cbuffer.data(hh_off + 274);

            auto g_xyzzz_xxxxz = cbuffer.data(hh_off + 275);

            auto g_xyzzz_xxxyy = cbuffer.data(hh_off + 276);

            auto g_xyzzz_xxxyz = cbuffer.data(hh_off + 277);

            auto g_xyzzz_xxxzz = cbuffer.data(hh_off + 278);

            auto g_xyzzz_xxyyy = cbuffer.data(hh_off + 279);

            auto g_xyzzz_xxyyz = cbuffer.data(hh_off + 280);

            auto g_xyzzz_xxyzz = cbuffer.data(hh_off + 281);

            auto g_xyzzz_xxzzz = cbuffer.data(hh_off + 282);

            auto g_xyzzz_xyyyy = cbuffer.data(hh_off + 283);

            auto g_xyzzz_xyyyz = cbuffer.data(hh_off + 284);

            auto g_xyzzz_xyyzz = cbuffer.data(hh_off + 285);

            auto g_xyzzz_xyzzz = cbuffer.data(hh_off + 286);

            auto g_xyzzz_xzzzz = cbuffer.data(hh_off + 287);

            auto g_xyzzz_yyyyy = cbuffer.data(hh_off + 288);

            auto g_xyzzz_yyyyz = cbuffer.data(hh_off + 289);

            auto g_xyzzz_yyyzz = cbuffer.data(hh_off + 290);

            auto g_xyzzz_yyzzz = cbuffer.data(hh_off + 291);

            auto g_xyzzz_yzzzz = cbuffer.data(hh_off + 292);

            auto g_xyzzz_zzzzz = cbuffer.data(hh_off + 293);

            auto g_xzzzz_xxxxx = cbuffer.data(hh_off + 294);

            auto g_xzzzz_xxxxy = cbuffer.data(hh_off + 295);

            auto g_xzzzz_xxxxz = cbuffer.data(hh_off + 296);

            auto g_xzzzz_xxxyy = cbuffer.data(hh_off + 297);

            auto g_xzzzz_xxxyz = cbuffer.data(hh_off + 298);

            auto g_xzzzz_xxxzz = cbuffer.data(hh_off + 299);

            auto g_xzzzz_xxyyy = cbuffer.data(hh_off + 300);

            auto g_xzzzz_xxyyz = cbuffer.data(hh_off + 301);

            auto g_xzzzz_xxyzz = cbuffer.data(hh_off + 302);

            auto g_xzzzz_xxzzz = cbuffer.data(hh_off + 303);

            auto g_xzzzz_xyyyy = cbuffer.data(hh_off + 304);

            auto g_xzzzz_xyyyz = cbuffer.data(hh_off + 305);

            auto g_xzzzz_xyyzz = cbuffer.data(hh_off + 306);

            auto g_xzzzz_xyzzz = cbuffer.data(hh_off + 307);

            auto g_xzzzz_xzzzz = cbuffer.data(hh_off + 308);

            auto g_xzzzz_yyyyy = cbuffer.data(hh_off + 309);

            auto g_xzzzz_yyyyz = cbuffer.data(hh_off + 310);

            auto g_xzzzz_yyyzz = cbuffer.data(hh_off + 311);

            auto g_xzzzz_yyzzz = cbuffer.data(hh_off + 312);

            auto g_xzzzz_yzzzz = cbuffer.data(hh_off + 313);

            auto g_xzzzz_zzzzz = cbuffer.data(hh_off + 314);

            auto g_yyyyy_xxxxx = cbuffer.data(hh_off + 315);

            auto g_yyyyy_xxxxy = cbuffer.data(hh_off + 316);

            auto g_yyyyy_xxxxz = cbuffer.data(hh_off + 317);

            auto g_yyyyy_xxxyy = cbuffer.data(hh_off + 318);

            auto g_yyyyy_xxxyz = cbuffer.data(hh_off + 319);

            auto g_yyyyy_xxxzz = cbuffer.data(hh_off + 320);

            auto g_yyyyy_xxyyy = cbuffer.data(hh_off + 321);

            auto g_yyyyy_xxyyz = cbuffer.data(hh_off + 322);

            auto g_yyyyy_xxyzz = cbuffer.data(hh_off + 323);

            auto g_yyyyy_xxzzz = cbuffer.data(hh_off + 324);

            auto g_yyyyy_xyyyy = cbuffer.data(hh_off + 325);

            auto g_yyyyy_xyyyz = cbuffer.data(hh_off + 326);

            auto g_yyyyy_xyyzz = cbuffer.data(hh_off + 327);

            auto g_yyyyy_xyzzz = cbuffer.data(hh_off + 328);

            auto g_yyyyy_xzzzz = cbuffer.data(hh_off + 329);

            auto g_yyyyy_yyyyy = cbuffer.data(hh_off + 330);

            auto g_yyyyy_yyyyz = cbuffer.data(hh_off + 331);

            auto g_yyyyy_yyyzz = cbuffer.data(hh_off + 332);

            auto g_yyyyy_yyzzz = cbuffer.data(hh_off + 333);

            auto g_yyyyy_yzzzz = cbuffer.data(hh_off + 334);

            auto g_yyyyy_zzzzz = cbuffer.data(hh_off + 335);

            auto g_yyyyz_xxxxx = cbuffer.data(hh_off + 336);

            auto g_yyyyz_xxxxy = cbuffer.data(hh_off + 337);

            auto g_yyyyz_xxxxz = cbuffer.data(hh_off + 338);

            auto g_yyyyz_xxxyy = cbuffer.data(hh_off + 339);

            auto g_yyyyz_xxxyz = cbuffer.data(hh_off + 340);

            auto g_yyyyz_xxxzz = cbuffer.data(hh_off + 341);

            auto g_yyyyz_xxyyy = cbuffer.data(hh_off + 342);

            auto g_yyyyz_xxyyz = cbuffer.data(hh_off + 343);

            auto g_yyyyz_xxyzz = cbuffer.data(hh_off + 344);

            auto g_yyyyz_xxzzz = cbuffer.data(hh_off + 345);

            auto g_yyyyz_xyyyy = cbuffer.data(hh_off + 346);

            auto g_yyyyz_xyyyz = cbuffer.data(hh_off + 347);

            auto g_yyyyz_xyyzz = cbuffer.data(hh_off + 348);

            auto g_yyyyz_xyzzz = cbuffer.data(hh_off + 349);

            auto g_yyyyz_xzzzz = cbuffer.data(hh_off + 350);

            auto g_yyyyz_yyyyy = cbuffer.data(hh_off + 351);

            auto g_yyyyz_yyyyz = cbuffer.data(hh_off + 352);

            auto g_yyyyz_yyyzz = cbuffer.data(hh_off + 353);

            auto g_yyyyz_yyzzz = cbuffer.data(hh_off + 354);

            auto g_yyyyz_yzzzz = cbuffer.data(hh_off + 355);

            auto g_yyyyz_zzzzz = cbuffer.data(hh_off + 356);

            auto g_yyyzz_xxxxx = cbuffer.data(hh_off + 357);

            auto g_yyyzz_xxxxy = cbuffer.data(hh_off + 358);

            auto g_yyyzz_xxxxz = cbuffer.data(hh_off + 359);

            auto g_yyyzz_xxxyy = cbuffer.data(hh_off + 360);

            auto g_yyyzz_xxxyz = cbuffer.data(hh_off + 361);

            auto g_yyyzz_xxxzz = cbuffer.data(hh_off + 362);

            auto g_yyyzz_xxyyy = cbuffer.data(hh_off + 363);

            auto g_yyyzz_xxyyz = cbuffer.data(hh_off + 364);

            auto g_yyyzz_xxyzz = cbuffer.data(hh_off + 365);

            auto g_yyyzz_xxzzz = cbuffer.data(hh_off + 366);

            auto g_yyyzz_xyyyy = cbuffer.data(hh_off + 367);

            auto g_yyyzz_xyyyz = cbuffer.data(hh_off + 368);

            auto g_yyyzz_xyyzz = cbuffer.data(hh_off + 369);

            auto g_yyyzz_xyzzz = cbuffer.data(hh_off + 370);

            auto g_yyyzz_xzzzz = cbuffer.data(hh_off + 371);

            auto g_yyyzz_yyyyy = cbuffer.data(hh_off + 372);

            auto g_yyyzz_yyyyz = cbuffer.data(hh_off + 373);

            auto g_yyyzz_yyyzz = cbuffer.data(hh_off + 374);

            auto g_yyyzz_yyzzz = cbuffer.data(hh_off + 375);

            auto g_yyyzz_yzzzz = cbuffer.data(hh_off + 376);

            auto g_yyyzz_zzzzz = cbuffer.data(hh_off + 377);

            auto g_yyzzz_xxxxx = cbuffer.data(hh_off + 378);

            auto g_yyzzz_xxxxy = cbuffer.data(hh_off + 379);

            auto g_yyzzz_xxxxz = cbuffer.data(hh_off + 380);

            auto g_yyzzz_xxxyy = cbuffer.data(hh_off + 381);

            auto g_yyzzz_xxxyz = cbuffer.data(hh_off + 382);

            auto g_yyzzz_xxxzz = cbuffer.data(hh_off + 383);

            auto g_yyzzz_xxyyy = cbuffer.data(hh_off + 384);

            auto g_yyzzz_xxyyz = cbuffer.data(hh_off + 385);

            auto g_yyzzz_xxyzz = cbuffer.data(hh_off + 386);

            auto g_yyzzz_xxzzz = cbuffer.data(hh_off + 387);

            auto g_yyzzz_xyyyy = cbuffer.data(hh_off + 388);

            auto g_yyzzz_xyyyz = cbuffer.data(hh_off + 389);

            auto g_yyzzz_xyyzz = cbuffer.data(hh_off + 390);

            auto g_yyzzz_xyzzz = cbuffer.data(hh_off + 391);

            auto g_yyzzz_xzzzz = cbuffer.data(hh_off + 392);

            auto g_yyzzz_yyyyy = cbuffer.data(hh_off + 393);

            auto g_yyzzz_yyyyz = cbuffer.data(hh_off + 394);

            auto g_yyzzz_yyyzz = cbuffer.data(hh_off + 395);

            auto g_yyzzz_yyzzz = cbuffer.data(hh_off + 396);

            auto g_yyzzz_yzzzz = cbuffer.data(hh_off + 397);

            auto g_yyzzz_zzzzz = cbuffer.data(hh_off + 398);

            auto g_yzzzz_xxxxx = cbuffer.data(hh_off + 399);

            auto g_yzzzz_xxxxy = cbuffer.data(hh_off + 400);

            auto g_yzzzz_xxxxz = cbuffer.data(hh_off + 401);

            auto g_yzzzz_xxxyy = cbuffer.data(hh_off + 402);

            auto g_yzzzz_xxxyz = cbuffer.data(hh_off + 403);

            auto g_yzzzz_xxxzz = cbuffer.data(hh_off + 404);

            auto g_yzzzz_xxyyy = cbuffer.data(hh_off + 405);

            auto g_yzzzz_xxyyz = cbuffer.data(hh_off + 406);

            auto g_yzzzz_xxyzz = cbuffer.data(hh_off + 407);

            auto g_yzzzz_xxzzz = cbuffer.data(hh_off + 408);

            auto g_yzzzz_xyyyy = cbuffer.data(hh_off + 409);

            auto g_yzzzz_xyyyz = cbuffer.data(hh_off + 410);

            auto g_yzzzz_xyyzz = cbuffer.data(hh_off + 411);

            auto g_yzzzz_xyzzz = cbuffer.data(hh_off + 412);

            auto g_yzzzz_xzzzz = cbuffer.data(hh_off + 413);

            auto g_yzzzz_yyyyy = cbuffer.data(hh_off + 414);

            auto g_yzzzz_yyyyz = cbuffer.data(hh_off + 415);

            auto g_yzzzz_yyyzz = cbuffer.data(hh_off + 416);

            auto g_yzzzz_yyzzz = cbuffer.data(hh_off + 417);

            auto g_yzzzz_yzzzz = cbuffer.data(hh_off + 418);

            auto g_yzzzz_zzzzz = cbuffer.data(hh_off + 419);

            auto g_zzzzz_xxxxx = cbuffer.data(hh_off + 420);

            auto g_zzzzz_xxxxy = cbuffer.data(hh_off + 421);

            auto g_zzzzz_xxxxz = cbuffer.data(hh_off + 422);

            auto g_zzzzz_xxxyy = cbuffer.data(hh_off + 423);

            auto g_zzzzz_xxxyz = cbuffer.data(hh_off + 424);

            auto g_zzzzz_xxxzz = cbuffer.data(hh_off + 425);

            auto g_zzzzz_xxyyy = cbuffer.data(hh_off + 426);

            auto g_zzzzz_xxyyz = cbuffer.data(hh_off + 427);

            auto g_zzzzz_xxyzz = cbuffer.data(hh_off + 428);

            auto g_zzzzz_xxzzz = cbuffer.data(hh_off + 429);

            auto g_zzzzz_xyyyy = cbuffer.data(hh_off + 430);

            auto g_zzzzz_xyyyz = cbuffer.data(hh_off + 431);

            auto g_zzzzz_xyyzz = cbuffer.data(hh_off + 432);

            auto g_zzzzz_xyzzz = cbuffer.data(hh_off + 433);

            auto g_zzzzz_xzzzz = cbuffer.data(hh_off + 434);

            auto g_zzzzz_yyyyy = cbuffer.data(hh_off + 435);

            auto g_zzzzz_yyyyz = cbuffer.data(hh_off + 436);

            auto g_zzzzz_yyyzz = cbuffer.data(hh_off + 437);

            auto g_zzzzz_yyzzz = cbuffer.data(hh_off + 438);

            auto g_zzzzz_yzzzz = cbuffer.data(hh_off + 439);

            auto g_zzzzz_zzzzz = cbuffer.data(hh_off + 440);

            /// Set up components of auxilary buffer : SSHH

            const auto hh_geom_10_off = idx_geom_10_xxhh + (i * bcomps + j) * 441;

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

            /// Set up components of auxilary buffer : SSHI

            const auto hi_geom_10_off = idx_geom_10_xxhi + (i * bcomps + j) * 588;

            auto g_x_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_x_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 335);

            auto g_x_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 363);

            auto g_x_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 377);

            auto g_x_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 391);

            auto g_x_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 398);

            auto g_x_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 419);

            auto g_x_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 420);

            auto g_x_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 421);

            auto g_x_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 422);

            auto g_x_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 423);

            auto g_x_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 424);

            auto g_x_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 425);

            auto g_x_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 426);

            auto g_x_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 427);

            auto g_x_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 428);

            auto g_x_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 429);

            auto g_x_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 430);

            auto g_x_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 431);

            auto g_x_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 432);

            auto g_x_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 433);

            auto g_x_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 434);

            auto g_x_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 435);

            auto g_x_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 436);

            auto g_x_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 437);

            auto g_x_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 438);

            auto g_x_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 439);

            auto g_x_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 440);

            auto g_x_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 441);

            auto g_x_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 442);

            auto g_x_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 443);

            auto g_x_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 444);

            auto g_x_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 445);

            auto g_x_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 446);

            auto g_x_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 447);

            auto g_x_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 448);

            auto g_x_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 449);

            auto g_x_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 450);

            auto g_x_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 451);

            auto g_x_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 452);

            auto g_x_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 453);

            auto g_x_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 454);

            auto g_x_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 455);

            auto g_x_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 456);

            auto g_x_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 457);

            auto g_x_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 458);

            auto g_x_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 459);

            auto g_x_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 460);

            auto g_x_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 461);

            auto g_x_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 462);

            auto g_x_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 463);

            auto g_x_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 464);

            auto g_x_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 465);

            auto g_x_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 466);

            auto g_x_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 467);

            auto g_x_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 468);

            auto g_x_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 469);

            auto g_x_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 470);

            auto g_x_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 471);

            auto g_x_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 472);

            auto g_x_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 473);

            auto g_x_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 474);

            auto g_x_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 475);

            auto g_x_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 476);

            auto g_x_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 477);

            auto g_x_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 478);

            auto g_x_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 479);

            auto g_x_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 480);

            auto g_x_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 481);

            auto g_x_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 482);

            auto g_x_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 483);

            auto g_x_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 484);

            auto g_x_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 485);

            auto g_x_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 486);

            auto g_x_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 487);

            auto g_x_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 488);

            auto g_x_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 489);

            auto g_x_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 490);

            auto g_x_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 491);

            auto g_x_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 492);

            auto g_x_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 493);

            auto g_x_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 494);

            auto g_x_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 495);

            auto g_x_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 496);

            auto g_x_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 497);

            auto g_x_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 498);

            auto g_x_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 499);

            auto g_x_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 500);

            auto g_x_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 501);

            auto g_x_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 502);

            auto g_x_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 503);

            auto g_x_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 504);

            auto g_x_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 505);

            auto g_x_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 506);

            auto g_x_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 507);

            auto g_x_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 508);

            auto g_x_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 509);

            auto g_x_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 510);

            auto g_x_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 511);

            auto g_x_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 512);

            auto g_x_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 513);

            auto g_x_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 514);

            auto g_x_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 515);

            auto g_x_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 516);

            auto g_x_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 517);

            auto g_x_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 518);

            auto g_x_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 519);

            auto g_x_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 520);

            auto g_x_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 521);

            auto g_x_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 522);

            auto g_x_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 523);

            auto g_x_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 524);

            auto g_x_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 525);

            auto g_x_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 526);

            auto g_x_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 527);

            auto g_x_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 528);

            auto g_x_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 529);

            auto g_x_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 530);

            auto g_x_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 531);

            auto g_x_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 532);

            auto g_x_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 533);

            auto g_x_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 534);

            auto g_x_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 535);

            auto g_x_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 536);

            auto g_x_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 537);

            auto g_x_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 538);

            auto g_x_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 539);

            auto g_x_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 540);

            auto g_x_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 541);

            auto g_x_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 542);

            auto g_x_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 543);

            auto g_x_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 544);

            auto g_x_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 545);

            auto g_x_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 546);

            auto g_x_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 547);

            auto g_x_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 548);

            auto g_x_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 549);

            auto g_x_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 550);

            auto g_x_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 551);

            auto g_x_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 552);

            auto g_x_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 553);

            auto g_x_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 554);

            auto g_x_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 555);

            auto g_x_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 556);

            auto g_x_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 557);

            auto g_x_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 558);

            auto g_x_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 559);

            auto g_x_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 560);

            auto g_x_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 561);

            auto g_x_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 562);

            auto g_x_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 563);

            auto g_x_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 564);

            auto g_x_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 565);

            auto g_x_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 566);

            auto g_x_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 567);

            auto g_x_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 568);

            auto g_x_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 569);

            auto g_x_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 570);

            auto g_x_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 571);

            auto g_x_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 572);

            auto g_x_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 573);

            auto g_x_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 574);

            auto g_x_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 575);

            auto g_x_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 576);

            auto g_x_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 577);

            auto g_x_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 578);

            auto g_x_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 579);

            auto g_x_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 580);

            auto g_x_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 581);

            auto g_x_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 582);

            auto g_x_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 583);

            auto g_x_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 584);

            auto g_x_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 585);

            auto g_x_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 586);

            auto g_x_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 0 * acomps * bcomps + 587);

            auto g_y_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 5);

            auto g_y_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 6);

            auto g_y_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 7);

            auto g_y_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 8);

            auto g_y_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 9);

            auto g_y_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 10);

            auto g_y_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 11);

            auto g_y_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 12);

            auto g_y_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 13);

            auto g_y_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 14);

            auto g_y_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 15);

            auto g_y_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 16);

            auto g_y_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 17);

            auto g_y_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 18);

            auto g_y_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 19);

            auto g_y_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 20);

            auto g_y_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 21);

            auto g_y_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 22);

            auto g_y_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 23);

            auto g_y_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 24);

            auto g_y_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 25);

            auto g_y_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 26);

            auto g_y_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 27);

            auto g_y_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 28);

            auto g_y_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 29);

            auto g_y_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 30);

            auto g_y_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 31);

            auto g_y_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 32);

            auto g_y_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 33);

            auto g_y_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 34);

            auto g_y_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 35);

            auto g_y_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 36);

            auto g_y_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 37);

            auto g_y_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 38);

            auto g_y_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 39);

            auto g_y_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 40);

            auto g_y_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 41);

            auto g_y_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 42);

            auto g_y_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 43);

            auto g_y_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 44);

            auto g_y_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 45);

            auto g_y_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 46);

            auto g_y_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 47);

            auto g_y_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 48);

            auto g_y_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 49);

            auto g_y_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 50);

            auto g_y_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 51);

            auto g_y_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 52);

            auto g_y_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 53);

            auto g_y_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 54);

            auto g_y_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 55);

            auto g_y_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 56);

            auto g_y_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 57);

            auto g_y_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 58);

            auto g_y_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 59);

            auto g_y_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 60);

            auto g_y_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 61);

            auto g_y_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 62);

            auto g_y_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 63);

            auto g_y_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 64);

            auto g_y_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 65);

            auto g_y_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 66);

            auto g_y_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 67);

            auto g_y_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 68);

            auto g_y_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 69);

            auto g_y_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 70);

            auto g_y_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 71);

            auto g_y_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 72);

            auto g_y_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 73);

            auto g_y_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 74);

            auto g_y_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 75);

            auto g_y_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 76);

            auto g_y_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 77);

            auto g_y_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 78);

            auto g_y_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 79);

            auto g_y_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 80);

            auto g_y_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 81);

            auto g_y_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 82);

            auto g_y_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 83);

            auto g_y_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 84);

            auto g_y_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 85);

            auto g_y_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 86);

            auto g_y_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 87);

            auto g_y_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 88);

            auto g_y_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 89);

            auto g_y_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 90);

            auto g_y_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 91);

            auto g_y_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 92);

            auto g_y_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 93);

            auto g_y_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 94);

            auto g_y_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 95);

            auto g_y_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 96);

            auto g_y_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 97);

            auto g_y_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 98);

            auto g_y_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 99);

            auto g_y_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 100);

            auto g_y_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 101);

            auto g_y_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 102);

            auto g_y_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 103);

            auto g_y_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 104);

            auto g_y_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 105);

            auto g_y_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 106);

            auto g_y_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 107);

            auto g_y_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 108);

            auto g_y_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 109);

            auto g_y_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 110);

            auto g_y_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 111);

            auto g_y_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 112);

            auto g_y_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 113);

            auto g_y_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 114);

            auto g_y_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 115);

            auto g_y_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 116);

            auto g_y_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 117);

            auto g_y_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 118);

            auto g_y_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 119);

            auto g_y_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 120);

            auto g_y_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 121);

            auto g_y_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 122);

            auto g_y_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 123);

            auto g_y_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 124);

            auto g_y_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 125);

            auto g_y_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 126);

            auto g_y_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 127);

            auto g_y_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 128);

            auto g_y_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 129);

            auto g_y_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 130);

            auto g_y_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 131);

            auto g_y_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 132);

            auto g_y_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 133);

            auto g_y_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 134);

            auto g_y_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 135);

            auto g_y_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 136);

            auto g_y_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 137);

            auto g_y_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 138);

            auto g_y_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 139);

            auto g_y_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 140);

            auto g_y_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 141);

            auto g_y_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 142);

            auto g_y_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 143);

            auto g_y_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 144);

            auto g_y_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 145);

            auto g_y_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 146);

            auto g_y_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 147);

            auto g_y_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 148);

            auto g_y_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 149);

            auto g_y_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 150);

            auto g_y_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 151);

            auto g_y_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 152);

            auto g_y_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 153);

            auto g_y_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 154);

            auto g_y_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 155);

            auto g_y_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 156);

            auto g_y_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 157);

            auto g_y_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 158);

            auto g_y_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 159);

            auto g_y_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 160);

            auto g_y_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 161);

            auto g_y_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 162);

            auto g_y_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 163);

            auto g_y_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 164);

            auto g_y_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 165);

            auto g_y_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 166);

            auto g_y_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 167);

            auto g_y_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 168);

            auto g_y_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 169);

            auto g_y_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 170);

            auto g_y_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 171);

            auto g_y_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 172);

            auto g_y_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 173);

            auto g_y_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 174);

            auto g_y_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 175);

            auto g_y_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 176);

            auto g_y_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 177);

            auto g_y_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 178);

            auto g_y_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 179);

            auto g_y_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 180);

            auto g_y_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 181);

            auto g_y_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 182);

            auto g_y_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 183);

            auto g_y_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 184);

            auto g_y_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 185);

            auto g_y_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 186);

            auto g_y_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 187);

            auto g_y_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 188);

            auto g_y_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 189);

            auto g_y_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 190);

            auto g_y_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 191);

            auto g_y_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 192);

            auto g_y_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 193);

            auto g_y_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 194);

            auto g_y_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 195);

            auto g_y_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 196);

            auto g_y_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 197);

            auto g_y_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 198);

            auto g_y_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 199);

            auto g_y_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 200);

            auto g_y_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 201);

            auto g_y_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 202);

            auto g_y_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 203);

            auto g_y_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 204);

            auto g_y_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 205);

            auto g_y_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 206);

            auto g_y_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 207);

            auto g_y_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 208);

            auto g_y_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 209);

            auto g_y_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 210);

            auto g_y_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 211);

            auto g_y_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 212);

            auto g_y_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 213);

            auto g_y_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 214);

            auto g_y_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 215);

            auto g_y_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 216);

            auto g_y_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 217);

            auto g_y_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 218);

            auto g_y_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 219);

            auto g_y_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 220);

            auto g_y_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 221);

            auto g_y_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 222);

            auto g_y_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 223);

            auto g_y_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 224);

            auto g_y_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 225);

            auto g_y_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 226);

            auto g_y_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 227);

            auto g_y_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 228);

            auto g_y_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 229);

            auto g_y_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 230);

            auto g_y_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 231);

            auto g_y_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 232);

            auto g_y_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 233);

            auto g_y_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 234);

            auto g_y_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 235);

            auto g_y_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 236);

            auto g_y_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 237);

            auto g_y_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 238);

            auto g_y_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 239);

            auto g_y_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 240);

            auto g_y_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 241);

            auto g_y_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 242);

            auto g_y_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 243);

            auto g_y_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 244);

            auto g_y_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 245);

            auto g_y_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 246);

            auto g_y_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 247);

            auto g_y_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 248);

            auto g_y_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 249);

            auto g_y_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 250);

            auto g_y_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 251);

            auto g_y_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 252);

            auto g_y_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 253);

            auto g_y_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 254);

            auto g_y_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 255);

            auto g_y_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 256);

            auto g_y_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 257);

            auto g_y_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 258);

            auto g_y_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 259);

            auto g_y_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 260);

            auto g_y_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 261);

            auto g_y_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 262);

            auto g_y_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 263);

            auto g_y_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 264);

            auto g_y_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 265);

            auto g_y_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 266);

            auto g_y_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 267);

            auto g_y_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 268);

            auto g_y_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 269);

            auto g_y_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 270);

            auto g_y_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 271);

            auto g_y_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 272);

            auto g_y_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 273);

            auto g_y_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 274);

            auto g_y_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 275);

            auto g_y_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 276);

            auto g_y_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 277);

            auto g_y_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 278);

            auto g_y_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 279);

            auto g_y_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 280);

            auto g_y_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 281);

            auto g_y_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 282);

            auto g_y_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 283);

            auto g_y_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 284);

            auto g_y_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 285);

            auto g_y_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 286);

            auto g_y_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 287);

            auto g_y_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 288);

            auto g_y_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 289);

            auto g_y_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 290);

            auto g_y_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 291);

            auto g_y_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 292);

            auto g_y_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 293);

            auto g_y_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 294);

            auto g_y_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 295);

            auto g_y_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 296);

            auto g_y_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 297);

            auto g_y_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 298);

            auto g_y_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 299);

            auto g_y_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 300);

            auto g_y_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 301);

            auto g_y_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 302);

            auto g_y_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 303);

            auto g_y_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 304);

            auto g_y_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 305);

            auto g_y_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 306);

            auto g_y_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 307);

            auto g_y_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 308);

            auto g_y_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 309);

            auto g_y_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 310);

            auto g_y_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 311);

            auto g_y_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 312);

            auto g_y_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 313);

            auto g_y_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 314);

            auto g_y_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 315);

            auto g_y_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 316);

            auto g_y_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 317);

            auto g_y_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 318);

            auto g_y_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 319);

            auto g_y_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 320);

            auto g_y_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 321);

            auto g_y_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 322);

            auto g_y_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 323);

            auto g_y_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 324);

            auto g_y_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 325);

            auto g_y_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 326);

            auto g_y_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 327);

            auto g_y_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 328);

            auto g_y_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 329);

            auto g_y_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 330);

            auto g_y_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 331);

            auto g_y_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 332);

            auto g_y_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 333);

            auto g_y_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 334);

            auto g_y_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 335);

            auto g_y_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 336);

            auto g_y_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 337);

            auto g_y_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 338);

            auto g_y_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 339);

            auto g_y_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 340);

            auto g_y_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 341);

            auto g_y_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 342);

            auto g_y_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 343);

            auto g_y_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 344);

            auto g_y_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 345);

            auto g_y_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 346);

            auto g_y_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 347);

            auto g_y_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 348);

            auto g_y_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 349);

            auto g_y_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 350);

            auto g_y_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 351);

            auto g_y_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 352);

            auto g_y_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 353);

            auto g_y_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 354);

            auto g_y_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 355);

            auto g_y_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 356);

            auto g_y_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 357);

            auto g_y_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 358);

            auto g_y_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 359);

            auto g_y_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 360);

            auto g_y_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 361);

            auto g_y_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 362);

            auto g_y_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 363);

            auto g_y_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 364);

            auto g_y_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 365);

            auto g_y_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 366);

            auto g_y_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 367);

            auto g_y_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 368);

            auto g_y_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 369);

            auto g_y_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 370);

            auto g_y_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 371);

            auto g_y_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 372);

            auto g_y_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 373);

            auto g_y_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 374);

            auto g_y_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 375);

            auto g_y_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 376);

            auto g_y_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 377);

            auto g_y_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 378);

            auto g_y_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 379);

            auto g_y_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 380);

            auto g_y_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 381);

            auto g_y_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 382);

            auto g_y_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 383);

            auto g_y_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 384);

            auto g_y_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 385);

            auto g_y_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 386);

            auto g_y_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 387);

            auto g_y_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 388);

            auto g_y_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 389);

            auto g_y_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 390);

            auto g_y_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 391);

            auto g_y_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 392);

            auto g_y_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 393);

            auto g_y_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 394);

            auto g_y_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 395);

            auto g_y_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 396);

            auto g_y_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 397);

            auto g_y_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 398);

            auto g_y_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 399);

            auto g_y_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 400);

            auto g_y_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 401);

            auto g_y_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 402);

            auto g_y_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 403);

            auto g_y_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 404);

            auto g_y_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 405);

            auto g_y_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 406);

            auto g_y_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 407);

            auto g_y_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 408);

            auto g_y_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 409);

            auto g_y_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 410);

            auto g_y_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 411);

            auto g_y_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 412);

            auto g_y_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 413);

            auto g_y_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 414);

            auto g_y_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 415);

            auto g_y_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 416);

            auto g_y_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 417);

            auto g_y_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 418);

            auto g_y_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 419);

            auto g_y_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 420);

            auto g_y_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 421);

            auto g_y_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 422);

            auto g_y_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 423);

            auto g_y_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 424);

            auto g_y_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 425);

            auto g_y_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 426);

            auto g_y_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 427);

            auto g_y_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 428);

            auto g_y_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 429);

            auto g_y_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 430);

            auto g_y_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 431);

            auto g_y_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 432);

            auto g_y_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 433);

            auto g_y_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 434);

            auto g_y_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 435);

            auto g_y_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 436);

            auto g_y_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 437);

            auto g_y_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 438);

            auto g_y_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 439);

            auto g_y_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 440);

            auto g_y_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 441);

            auto g_y_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 442);

            auto g_y_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 443);

            auto g_y_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 444);

            auto g_y_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 445);

            auto g_y_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 446);

            auto g_y_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 447);

            auto g_y_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 448);

            auto g_y_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 449);

            auto g_y_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 450);

            auto g_y_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 451);

            auto g_y_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 452);

            auto g_y_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 453);

            auto g_y_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 454);

            auto g_y_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 455);

            auto g_y_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 456);

            auto g_y_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 457);

            auto g_y_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 458);

            auto g_y_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 459);

            auto g_y_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 460);

            auto g_y_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 461);

            auto g_y_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 462);

            auto g_y_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 463);

            auto g_y_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 464);

            auto g_y_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 465);

            auto g_y_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 466);

            auto g_y_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 467);

            auto g_y_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 468);

            auto g_y_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 469);

            auto g_y_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 470);

            auto g_y_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 471);

            auto g_y_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 472);

            auto g_y_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 473);

            auto g_y_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 474);

            auto g_y_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 475);

            auto g_y_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 476);

            auto g_y_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 477);

            auto g_y_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 478);

            auto g_y_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 479);

            auto g_y_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 480);

            auto g_y_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 481);

            auto g_y_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 482);

            auto g_y_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 483);

            auto g_y_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 484);

            auto g_y_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 485);

            auto g_y_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 486);

            auto g_y_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 487);

            auto g_y_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 488);

            auto g_y_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 489);

            auto g_y_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 490);

            auto g_y_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 491);

            auto g_y_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 492);

            auto g_y_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 493);

            auto g_y_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 494);

            auto g_y_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 495);

            auto g_y_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 496);

            auto g_y_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 497);

            auto g_y_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 498);

            auto g_y_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 499);

            auto g_y_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 500);

            auto g_y_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 501);

            auto g_y_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 502);

            auto g_y_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 503);

            auto g_y_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 504);

            auto g_y_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 505);

            auto g_y_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 506);

            auto g_y_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 507);

            auto g_y_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 508);

            auto g_y_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 509);

            auto g_y_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 510);

            auto g_y_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 511);

            auto g_y_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 512);

            auto g_y_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 513);

            auto g_y_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 514);

            auto g_y_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 515);

            auto g_y_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 516);

            auto g_y_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 517);

            auto g_y_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 518);

            auto g_y_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 519);

            auto g_y_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 520);

            auto g_y_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 521);

            auto g_y_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 522);

            auto g_y_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 523);

            auto g_y_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 524);

            auto g_y_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 525);

            auto g_y_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 526);

            auto g_y_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 527);

            auto g_y_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 528);

            auto g_y_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 529);

            auto g_y_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 530);

            auto g_y_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 531);

            auto g_y_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 532);

            auto g_y_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 533);

            auto g_y_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 534);

            auto g_y_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 535);

            auto g_y_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 536);

            auto g_y_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 537);

            auto g_y_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 538);

            auto g_y_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 539);

            auto g_y_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 540);

            auto g_y_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 541);

            auto g_y_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 542);

            auto g_y_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 543);

            auto g_y_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 544);

            auto g_y_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 545);

            auto g_y_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 546);

            auto g_y_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 547);

            auto g_y_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 548);

            auto g_y_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 549);

            auto g_y_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 550);

            auto g_y_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 551);

            auto g_y_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 552);

            auto g_y_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 553);

            auto g_y_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 554);

            auto g_y_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 555);

            auto g_y_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 556);

            auto g_y_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 557);

            auto g_y_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 558);

            auto g_y_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 559);

            auto g_y_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 560);

            auto g_y_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 561);

            auto g_y_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 562);

            auto g_y_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 563);

            auto g_y_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 564);

            auto g_y_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 565);

            auto g_y_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 566);

            auto g_y_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 567);

            auto g_y_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 568);

            auto g_y_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 569);

            auto g_y_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 570);

            auto g_y_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 571);

            auto g_y_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 572);

            auto g_y_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 573);

            auto g_y_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 574);

            auto g_y_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 575);

            auto g_y_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 576);

            auto g_y_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 577);

            auto g_y_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 578);

            auto g_y_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 579);

            auto g_y_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 580);

            auto g_y_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 581);

            auto g_y_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 582);

            auto g_y_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 583);

            auto g_y_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 584);

            auto g_y_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 585);

            auto g_y_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 586);

            auto g_y_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 588 * acomps * bcomps + 587);

            auto g_z_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 5);

            auto g_z_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 6);

            auto g_z_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 7);

            auto g_z_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 8);

            auto g_z_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 9);

            auto g_z_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 10);

            auto g_z_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 11);

            auto g_z_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 12);

            auto g_z_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 13);

            auto g_z_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 14);

            auto g_z_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 15);

            auto g_z_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 16);

            auto g_z_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 17);

            auto g_z_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 18);

            auto g_z_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 19);

            auto g_z_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 20);

            auto g_z_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 21);

            auto g_z_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 22);

            auto g_z_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 23);

            auto g_z_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 24);

            auto g_z_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 25);

            auto g_z_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 26);

            auto g_z_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 27);

            auto g_z_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 28);

            auto g_z_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 29);

            auto g_z_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 30);

            auto g_z_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 31);

            auto g_z_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 32);

            auto g_z_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 33);

            auto g_z_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 34);

            auto g_z_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 35);

            auto g_z_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 36);

            auto g_z_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 37);

            auto g_z_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 38);

            auto g_z_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 39);

            auto g_z_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 40);

            auto g_z_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 41);

            auto g_z_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 42);

            auto g_z_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 43);

            auto g_z_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 44);

            auto g_z_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 45);

            auto g_z_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 46);

            auto g_z_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 47);

            auto g_z_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 48);

            auto g_z_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 49);

            auto g_z_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 50);

            auto g_z_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 51);

            auto g_z_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 52);

            auto g_z_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 53);

            auto g_z_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 54);

            auto g_z_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 55);

            auto g_z_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 56);

            auto g_z_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 57);

            auto g_z_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 58);

            auto g_z_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 59);

            auto g_z_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 60);

            auto g_z_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 61);

            auto g_z_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 62);

            auto g_z_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 63);

            auto g_z_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 64);

            auto g_z_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 65);

            auto g_z_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 66);

            auto g_z_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 67);

            auto g_z_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 68);

            auto g_z_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 69);

            auto g_z_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 70);

            auto g_z_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 71);

            auto g_z_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 72);

            auto g_z_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 73);

            auto g_z_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 74);

            auto g_z_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 75);

            auto g_z_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 76);

            auto g_z_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 77);

            auto g_z_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 78);

            auto g_z_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 79);

            auto g_z_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 80);

            auto g_z_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 81);

            auto g_z_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 82);

            auto g_z_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 83);

            auto g_z_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 84);

            auto g_z_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 85);

            auto g_z_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 86);

            auto g_z_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 87);

            auto g_z_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 88);

            auto g_z_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 89);

            auto g_z_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 90);

            auto g_z_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 91);

            auto g_z_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 92);

            auto g_z_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 93);

            auto g_z_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 94);

            auto g_z_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 95);

            auto g_z_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 96);

            auto g_z_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 97);

            auto g_z_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 98);

            auto g_z_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 99);

            auto g_z_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 100);

            auto g_z_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 101);

            auto g_z_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 102);

            auto g_z_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 103);

            auto g_z_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 104);

            auto g_z_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 105);

            auto g_z_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 106);

            auto g_z_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 107);

            auto g_z_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 108);

            auto g_z_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 109);

            auto g_z_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 110);

            auto g_z_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 111);

            auto g_z_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 112);

            auto g_z_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 113);

            auto g_z_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 114);

            auto g_z_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 115);

            auto g_z_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 116);

            auto g_z_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 117);

            auto g_z_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 118);

            auto g_z_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 119);

            auto g_z_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 120);

            auto g_z_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 121);

            auto g_z_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 122);

            auto g_z_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 123);

            auto g_z_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 124);

            auto g_z_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 125);

            auto g_z_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 126);

            auto g_z_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 127);

            auto g_z_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 128);

            auto g_z_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 129);

            auto g_z_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 130);

            auto g_z_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 131);

            auto g_z_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 132);

            auto g_z_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 133);

            auto g_z_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 134);

            auto g_z_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 135);

            auto g_z_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 136);

            auto g_z_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 137);

            auto g_z_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 138);

            auto g_z_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 139);

            auto g_z_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 140);

            auto g_z_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 141);

            auto g_z_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 142);

            auto g_z_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 143);

            auto g_z_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 144);

            auto g_z_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 145);

            auto g_z_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 146);

            auto g_z_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 147);

            auto g_z_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 148);

            auto g_z_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 149);

            auto g_z_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 150);

            auto g_z_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 151);

            auto g_z_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 152);

            auto g_z_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 153);

            auto g_z_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 154);

            auto g_z_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 155);

            auto g_z_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 156);

            auto g_z_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 157);

            auto g_z_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 158);

            auto g_z_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 159);

            auto g_z_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 160);

            auto g_z_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 161);

            auto g_z_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 162);

            auto g_z_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 163);

            auto g_z_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 164);

            auto g_z_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 165);

            auto g_z_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 166);

            auto g_z_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 167);

            auto g_z_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 168);

            auto g_z_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 169);

            auto g_z_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 170);

            auto g_z_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 171);

            auto g_z_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 172);

            auto g_z_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 173);

            auto g_z_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 174);

            auto g_z_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 175);

            auto g_z_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 176);

            auto g_z_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 177);

            auto g_z_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 178);

            auto g_z_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 179);

            auto g_z_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 180);

            auto g_z_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 181);

            auto g_z_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 182);

            auto g_z_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 183);

            auto g_z_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 184);

            auto g_z_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 185);

            auto g_z_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 186);

            auto g_z_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 187);

            auto g_z_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 188);

            auto g_z_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 189);

            auto g_z_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 190);

            auto g_z_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 191);

            auto g_z_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 192);

            auto g_z_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 193);

            auto g_z_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 194);

            auto g_z_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 195);

            auto g_z_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 196);

            auto g_z_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 197);

            auto g_z_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 198);

            auto g_z_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 199);

            auto g_z_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 200);

            auto g_z_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 201);

            auto g_z_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 202);

            auto g_z_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 203);

            auto g_z_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 204);

            auto g_z_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 205);

            auto g_z_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 206);

            auto g_z_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 207);

            auto g_z_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 208);

            auto g_z_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 209);

            auto g_z_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 210);

            auto g_z_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 211);

            auto g_z_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 212);

            auto g_z_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 213);

            auto g_z_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 214);

            auto g_z_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 215);

            auto g_z_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 216);

            auto g_z_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 217);

            auto g_z_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 218);

            auto g_z_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 219);

            auto g_z_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 220);

            auto g_z_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 221);

            auto g_z_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 222);

            auto g_z_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 223);

            auto g_z_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 224);

            auto g_z_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 225);

            auto g_z_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 226);

            auto g_z_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 227);

            auto g_z_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 228);

            auto g_z_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 229);

            auto g_z_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 230);

            auto g_z_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 231);

            auto g_z_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 232);

            auto g_z_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 233);

            auto g_z_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 234);

            auto g_z_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 235);

            auto g_z_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 236);

            auto g_z_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 237);

            auto g_z_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 238);

            auto g_z_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 239);

            auto g_z_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 240);

            auto g_z_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 241);

            auto g_z_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 242);

            auto g_z_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 243);

            auto g_z_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 244);

            auto g_z_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 245);

            auto g_z_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 246);

            auto g_z_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 247);

            auto g_z_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 248);

            auto g_z_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 249);

            auto g_z_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 250);

            auto g_z_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 251);

            auto g_z_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 252);

            auto g_z_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 253);

            auto g_z_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 254);

            auto g_z_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 255);

            auto g_z_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 256);

            auto g_z_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 257);

            auto g_z_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 258);

            auto g_z_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 259);

            auto g_z_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 260);

            auto g_z_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 261);

            auto g_z_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 262);

            auto g_z_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 263);

            auto g_z_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 264);

            auto g_z_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 265);

            auto g_z_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 266);

            auto g_z_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 267);

            auto g_z_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 268);

            auto g_z_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 269);

            auto g_z_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 270);

            auto g_z_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 271);

            auto g_z_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 272);

            auto g_z_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 273);

            auto g_z_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 274);

            auto g_z_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 275);

            auto g_z_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 276);

            auto g_z_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 277);

            auto g_z_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 278);

            auto g_z_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 279);

            auto g_z_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 280);

            auto g_z_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 281);

            auto g_z_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 282);

            auto g_z_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 283);

            auto g_z_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 284);

            auto g_z_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 285);

            auto g_z_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 286);

            auto g_z_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 287);

            auto g_z_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 288);

            auto g_z_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 289);

            auto g_z_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 290);

            auto g_z_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 291);

            auto g_z_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 292);

            auto g_z_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 293);

            auto g_z_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 294);

            auto g_z_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 295);

            auto g_z_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 296);

            auto g_z_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 297);

            auto g_z_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 298);

            auto g_z_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 299);

            auto g_z_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 300);

            auto g_z_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 301);

            auto g_z_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 302);

            auto g_z_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 303);

            auto g_z_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 304);

            auto g_z_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 305);

            auto g_z_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 306);

            auto g_z_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 307);

            auto g_z_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 308);

            auto g_z_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 309);

            auto g_z_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 310);

            auto g_z_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 311);

            auto g_z_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 312);

            auto g_z_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 313);

            auto g_z_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 314);

            auto g_z_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 315);

            auto g_z_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 316);

            auto g_z_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 317);

            auto g_z_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 318);

            auto g_z_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 319);

            auto g_z_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 320);

            auto g_z_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 321);

            auto g_z_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 322);

            auto g_z_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 323);

            auto g_z_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 324);

            auto g_z_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 325);

            auto g_z_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 326);

            auto g_z_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 327);

            auto g_z_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 328);

            auto g_z_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 329);

            auto g_z_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 330);

            auto g_z_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 331);

            auto g_z_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 332);

            auto g_z_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 333);

            auto g_z_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 334);

            auto g_z_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 335);

            auto g_z_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 336);

            auto g_z_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 337);

            auto g_z_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 338);

            auto g_z_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 339);

            auto g_z_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 340);

            auto g_z_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 341);

            auto g_z_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 342);

            auto g_z_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 343);

            auto g_z_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 344);

            auto g_z_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 345);

            auto g_z_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 346);

            auto g_z_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 347);

            auto g_z_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 348);

            auto g_z_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 349);

            auto g_z_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 350);

            auto g_z_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 351);

            auto g_z_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 352);

            auto g_z_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 353);

            auto g_z_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 354);

            auto g_z_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 355);

            auto g_z_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 356);

            auto g_z_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 357);

            auto g_z_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 358);

            auto g_z_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 359);

            auto g_z_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 360);

            auto g_z_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 361);

            auto g_z_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 362);

            auto g_z_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 363);

            auto g_z_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 364);

            auto g_z_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 365);

            auto g_z_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 366);

            auto g_z_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 367);

            auto g_z_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 368);

            auto g_z_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 369);

            auto g_z_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 370);

            auto g_z_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 371);

            auto g_z_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 372);

            auto g_z_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 373);

            auto g_z_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 374);

            auto g_z_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 375);

            auto g_z_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 376);

            auto g_z_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 377);

            auto g_z_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 378);

            auto g_z_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 379);

            auto g_z_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 380);

            auto g_z_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 381);

            auto g_z_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 382);

            auto g_z_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 383);

            auto g_z_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 384);

            auto g_z_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 385);

            auto g_z_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 386);

            auto g_z_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 387);

            auto g_z_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 388);

            auto g_z_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 389);

            auto g_z_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 390);

            auto g_z_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 391);

            auto g_z_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 392);

            auto g_z_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 393);

            auto g_z_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 394);

            auto g_z_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 395);

            auto g_z_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 396);

            auto g_z_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 397);

            auto g_z_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 398);

            auto g_z_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 399);

            auto g_z_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 400);

            auto g_z_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 401);

            auto g_z_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 402);

            auto g_z_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 403);

            auto g_z_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 404);

            auto g_z_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 405);

            auto g_z_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 406);

            auto g_z_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 407);

            auto g_z_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 408);

            auto g_z_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 409);

            auto g_z_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 410);

            auto g_z_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 411);

            auto g_z_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 412);

            auto g_z_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 413);

            auto g_z_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 414);

            auto g_z_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 415);

            auto g_z_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 416);

            auto g_z_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 417);

            auto g_z_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 418);

            auto g_z_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 419);

            auto g_z_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 420);

            auto g_z_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 421);

            auto g_z_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 422);

            auto g_z_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 423);

            auto g_z_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 424);

            auto g_z_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 425);

            auto g_z_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 426);

            auto g_z_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 427);

            auto g_z_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 428);

            auto g_z_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 429);

            auto g_z_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 430);

            auto g_z_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 431);

            auto g_z_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 432);

            auto g_z_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 433);

            auto g_z_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 434);

            auto g_z_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 435);

            auto g_z_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 436);

            auto g_z_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 437);

            auto g_z_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 438);

            auto g_z_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 439);

            auto g_z_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 440);

            auto g_z_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 441);

            auto g_z_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 442);

            auto g_z_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 443);

            auto g_z_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 444);

            auto g_z_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 445);

            auto g_z_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 446);

            auto g_z_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 447);

            auto g_z_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 448);

            auto g_z_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 449);

            auto g_z_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 450);

            auto g_z_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 451);

            auto g_z_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 452);

            auto g_z_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 453);

            auto g_z_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 454);

            auto g_z_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 455);

            auto g_z_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 456);

            auto g_z_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 457);

            auto g_z_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 458);

            auto g_z_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 459);

            auto g_z_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 460);

            auto g_z_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 461);

            auto g_z_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 462);

            auto g_z_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 463);

            auto g_z_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 464);

            auto g_z_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 465);

            auto g_z_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 466);

            auto g_z_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 467);

            auto g_z_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 468);

            auto g_z_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 469);

            auto g_z_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 470);

            auto g_z_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 471);

            auto g_z_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 472);

            auto g_z_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 473);

            auto g_z_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 474);

            auto g_z_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 475);

            auto g_z_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 476);

            auto g_z_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 477);

            auto g_z_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 478);

            auto g_z_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 479);

            auto g_z_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 480);

            auto g_z_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 481);

            auto g_z_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 482);

            auto g_z_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 483);

            auto g_z_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 484);

            auto g_z_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 485);

            auto g_z_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 486);

            auto g_z_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 487);

            auto g_z_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 488);

            auto g_z_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 489);

            auto g_z_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 490);

            auto g_z_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 491);

            auto g_z_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 492);

            auto g_z_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 493);

            auto g_z_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 494);

            auto g_z_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 495);

            auto g_z_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 496);

            auto g_z_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 497);

            auto g_z_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 498);

            auto g_z_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 499);

            auto g_z_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 500);

            auto g_z_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 501);

            auto g_z_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 502);

            auto g_z_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 503);

            auto g_z_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 504);

            auto g_z_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 505);

            auto g_z_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 506);

            auto g_z_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 507);

            auto g_z_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 508);

            auto g_z_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 509);

            auto g_z_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 510);

            auto g_z_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 511);

            auto g_z_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 512);

            auto g_z_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 513);

            auto g_z_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 514);

            auto g_z_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 515);

            auto g_z_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 516);

            auto g_z_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 517);

            auto g_z_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 518);

            auto g_z_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 519);

            auto g_z_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 520);

            auto g_z_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 521);

            auto g_z_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 522);

            auto g_z_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 523);

            auto g_z_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 524);

            auto g_z_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 525);

            auto g_z_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 526);

            auto g_z_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 527);

            auto g_z_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 528);

            auto g_z_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 529);

            auto g_z_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 530);

            auto g_z_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 531);

            auto g_z_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 532);

            auto g_z_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 533);

            auto g_z_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 534);

            auto g_z_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 535);

            auto g_z_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 536);

            auto g_z_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 537);

            auto g_z_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 538);

            auto g_z_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 539);

            auto g_z_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 540);

            auto g_z_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 541);

            auto g_z_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 542);

            auto g_z_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 543);

            auto g_z_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 544);

            auto g_z_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 545);

            auto g_z_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 546);

            auto g_z_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 547);

            auto g_z_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 548);

            auto g_z_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 549);

            auto g_z_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 550);

            auto g_z_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 551);

            auto g_z_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 552);

            auto g_z_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 553);

            auto g_z_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 554);

            auto g_z_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 555);

            auto g_z_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 556);

            auto g_z_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 557);

            auto g_z_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 558);

            auto g_z_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 559);

            auto g_z_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 560);

            auto g_z_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 561);

            auto g_z_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 562);

            auto g_z_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 563);

            auto g_z_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 564);

            auto g_z_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 565);

            auto g_z_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 566);

            auto g_z_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 567);

            auto g_z_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 568);

            auto g_z_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 569);

            auto g_z_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 570);

            auto g_z_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 571);

            auto g_z_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 572);

            auto g_z_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 573);

            auto g_z_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 574);

            auto g_z_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 575);

            auto g_z_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 576);

            auto g_z_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 577);

            auto g_z_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 578);

            auto g_z_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 579);

            auto g_z_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 580);

            auto g_z_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 581);

            auto g_z_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 582);

            auto g_z_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 583);

            auto g_z_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 584);

            auto g_z_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 585);

            auto g_z_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 586);

            auto g_z_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1176 * acomps * bcomps + 587);

            /// set up bra offset for contr_buffer_xxih

            const auto ih_geom_10_off = idx_geom_10_xxih + (i * bcomps + j) * 588;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxxx_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxxx_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxxx_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxxx_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxxx_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxxx_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxxx_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxxx_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxxx_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxxx_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxxx_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxxx_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxxx_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxxx_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxxx_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxxx_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxxx_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxxx_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxxx_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxxx_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxxx, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_zzzzz, g_x_0_xxxxxx_xxxxx, g_x_0_xxxxxx_xxxxy, g_x_0_xxxxxx_xxxxz, g_x_0_xxxxxx_xxxyy, g_x_0_xxxxxx_xxxyz, g_x_0_xxxxxx_xxxzz, g_x_0_xxxxxx_xxyyy, g_x_0_xxxxxx_xxyyz, g_x_0_xxxxxx_xxyzz, g_x_0_xxxxxx_xxzzz, g_x_0_xxxxxx_xyyyy, g_x_0_xxxxxx_xyyyz, g_x_0_xxxxxx_xyyzz, g_x_0_xxxxxx_xyzzz, g_x_0_xxxxxx_xzzzz, g_x_0_xxxxxx_yyyyy, g_x_0_xxxxxx_yyyyz, g_x_0_xxxxxx_yyyzz, g_x_0_xxxxxx_yyzzz, g_x_0_xxxxxx_yzzzz, g_x_0_xxxxxx_zzzzz, g_xxxxx_xxxxx, g_xxxxx_xxxxy, g_xxxxx_xxxxz, g_xxxxx_xxxyy, g_xxxxx_xxxyz, g_xxxxx_xxxzz, g_xxxxx_xxyyy, g_xxxxx_xxyyz, g_xxxxx_xxyzz, g_xxxxx_xxzzz, g_xxxxx_xyyyy, g_xxxxx_xyyyz, g_xxxxx_xyyzz, g_xxxxx_xyzzz, g_xxxxx_xzzzz, g_xxxxx_yyyyy, g_xxxxx_yyyyz, g_xxxxx_yyyzz, g_xxxxx_yyzzz, g_xxxxx_yzzzz, g_xxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_xxxxx[k] = -g_xxxxx_xxxxx[k] - g_x_0_xxxxx_xxxxx[k] * cd_x[k] + g_x_0_xxxxx_xxxxxx[k];

                g_x_0_xxxxxx_xxxxy[k] = -g_xxxxx_xxxxy[k] - g_x_0_xxxxx_xxxxy[k] * cd_x[k] + g_x_0_xxxxx_xxxxxy[k];

                g_x_0_xxxxxx_xxxxz[k] = -g_xxxxx_xxxxz[k] - g_x_0_xxxxx_xxxxz[k] * cd_x[k] + g_x_0_xxxxx_xxxxxz[k];

                g_x_0_xxxxxx_xxxyy[k] = -g_xxxxx_xxxyy[k] - g_x_0_xxxxx_xxxyy[k] * cd_x[k] + g_x_0_xxxxx_xxxxyy[k];

                g_x_0_xxxxxx_xxxyz[k] = -g_xxxxx_xxxyz[k] - g_x_0_xxxxx_xxxyz[k] * cd_x[k] + g_x_0_xxxxx_xxxxyz[k];

                g_x_0_xxxxxx_xxxzz[k] = -g_xxxxx_xxxzz[k] - g_x_0_xxxxx_xxxzz[k] * cd_x[k] + g_x_0_xxxxx_xxxxzz[k];

                g_x_0_xxxxxx_xxyyy[k] = -g_xxxxx_xxyyy[k] - g_x_0_xxxxx_xxyyy[k] * cd_x[k] + g_x_0_xxxxx_xxxyyy[k];

                g_x_0_xxxxxx_xxyyz[k] = -g_xxxxx_xxyyz[k] - g_x_0_xxxxx_xxyyz[k] * cd_x[k] + g_x_0_xxxxx_xxxyyz[k];

                g_x_0_xxxxxx_xxyzz[k] = -g_xxxxx_xxyzz[k] - g_x_0_xxxxx_xxyzz[k] * cd_x[k] + g_x_0_xxxxx_xxxyzz[k];

                g_x_0_xxxxxx_xxzzz[k] = -g_xxxxx_xxzzz[k] - g_x_0_xxxxx_xxzzz[k] * cd_x[k] + g_x_0_xxxxx_xxxzzz[k];

                g_x_0_xxxxxx_xyyyy[k] = -g_xxxxx_xyyyy[k] - g_x_0_xxxxx_xyyyy[k] * cd_x[k] + g_x_0_xxxxx_xxyyyy[k];

                g_x_0_xxxxxx_xyyyz[k] = -g_xxxxx_xyyyz[k] - g_x_0_xxxxx_xyyyz[k] * cd_x[k] + g_x_0_xxxxx_xxyyyz[k];

                g_x_0_xxxxxx_xyyzz[k] = -g_xxxxx_xyyzz[k] - g_x_0_xxxxx_xyyzz[k] * cd_x[k] + g_x_0_xxxxx_xxyyzz[k];

                g_x_0_xxxxxx_xyzzz[k] = -g_xxxxx_xyzzz[k] - g_x_0_xxxxx_xyzzz[k] * cd_x[k] + g_x_0_xxxxx_xxyzzz[k];

                g_x_0_xxxxxx_xzzzz[k] = -g_xxxxx_xzzzz[k] - g_x_0_xxxxx_xzzzz[k] * cd_x[k] + g_x_0_xxxxx_xxzzzz[k];

                g_x_0_xxxxxx_yyyyy[k] = -g_xxxxx_yyyyy[k] - g_x_0_xxxxx_yyyyy[k] * cd_x[k] + g_x_0_xxxxx_xyyyyy[k];

                g_x_0_xxxxxx_yyyyz[k] = -g_xxxxx_yyyyz[k] - g_x_0_xxxxx_yyyyz[k] * cd_x[k] + g_x_0_xxxxx_xyyyyz[k];

                g_x_0_xxxxxx_yyyzz[k] = -g_xxxxx_yyyzz[k] - g_x_0_xxxxx_yyyzz[k] * cd_x[k] + g_x_0_xxxxx_xyyyzz[k];

                g_x_0_xxxxxx_yyzzz[k] = -g_xxxxx_yyzzz[k] - g_x_0_xxxxx_yyzzz[k] * cd_x[k] + g_x_0_xxxxx_xyyzzz[k];

                g_x_0_xxxxxx_yzzzz[k] = -g_xxxxx_yzzzz[k] - g_x_0_xxxxx_yzzzz[k] * cd_x[k] + g_x_0_xxxxx_xyzzzz[k];

                g_x_0_xxxxxx_zzzzz[k] = -g_xxxxx_zzzzz[k] - g_x_0_xxxxx_zzzzz[k] * cd_x[k] + g_x_0_xxxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxxy_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxxy_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxxy_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxxy_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxxy_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxxy_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxxy_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxxy_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxxxy_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxxy_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxxy_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxxy_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxxy_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxxy_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxxy_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxxy_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxxy_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxxy_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxxxy_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxxy_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyyy, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_zzzzz, g_x_0_xxxxxy_xxxxx, g_x_0_xxxxxy_xxxxy, g_x_0_xxxxxy_xxxxz, g_x_0_xxxxxy_xxxyy, g_x_0_xxxxxy_xxxyz, g_x_0_xxxxxy_xxxzz, g_x_0_xxxxxy_xxyyy, g_x_0_xxxxxy_xxyyz, g_x_0_xxxxxy_xxyzz, g_x_0_xxxxxy_xxzzz, g_x_0_xxxxxy_xyyyy, g_x_0_xxxxxy_xyyyz, g_x_0_xxxxxy_xyyzz, g_x_0_xxxxxy_xyzzz, g_x_0_xxxxxy_xzzzz, g_x_0_xxxxxy_yyyyy, g_x_0_xxxxxy_yyyyz, g_x_0_xxxxxy_yyyzz, g_x_0_xxxxxy_yyzzz, g_x_0_xxxxxy_yzzzz, g_x_0_xxxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_xxxxx[k] = -g_x_0_xxxxx_xxxxx[k] * cd_y[k] + g_x_0_xxxxx_xxxxxy[k];

                g_x_0_xxxxxy_xxxxy[k] = -g_x_0_xxxxx_xxxxy[k] * cd_y[k] + g_x_0_xxxxx_xxxxyy[k];

                g_x_0_xxxxxy_xxxxz[k] = -g_x_0_xxxxx_xxxxz[k] * cd_y[k] + g_x_0_xxxxx_xxxxyz[k];

                g_x_0_xxxxxy_xxxyy[k] = -g_x_0_xxxxx_xxxyy[k] * cd_y[k] + g_x_0_xxxxx_xxxyyy[k];

                g_x_0_xxxxxy_xxxyz[k] = -g_x_0_xxxxx_xxxyz[k] * cd_y[k] + g_x_0_xxxxx_xxxyyz[k];

                g_x_0_xxxxxy_xxxzz[k] = -g_x_0_xxxxx_xxxzz[k] * cd_y[k] + g_x_0_xxxxx_xxxyzz[k];

                g_x_0_xxxxxy_xxyyy[k] = -g_x_0_xxxxx_xxyyy[k] * cd_y[k] + g_x_0_xxxxx_xxyyyy[k];

                g_x_0_xxxxxy_xxyyz[k] = -g_x_0_xxxxx_xxyyz[k] * cd_y[k] + g_x_0_xxxxx_xxyyyz[k];

                g_x_0_xxxxxy_xxyzz[k] = -g_x_0_xxxxx_xxyzz[k] * cd_y[k] + g_x_0_xxxxx_xxyyzz[k];

                g_x_0_xxxxxy_xxzzz[k] = -g_x_0_xxxxx_xxzzz[k] * cd_y[k] + g_x_0_xxxxx_xxyzzz[k];

                g_x_0_xxxxxy_xyyyy[k] = -g_x_0_xxxxx_xyyyy[k] * cd_y[k] + g_x_0_xxxxx_xyyyyy[k];

                g_x_0_xxxxxy_xyyyz[k] = -g_x_0_xxxxx_xyyyz[k] * cd_y[k] + g_x_0_xxxxx_xyyyyz[k];

                g_x_0_xxxxxy_xyyzz[k] = -g_x_0_xxxxx_xyyzz[k] * cd_y[k] + g_x_0_xxxxx_xyyyzz[k];

                g_x_0_xxxxxy_xyzzz[k] = -g_x_0_xxxxx_xyzzz[k] * cd_y[k] + g_x_0_xxxxx_xyyzzz[k];

                g_x_0_xxxxxy_xzzzz[k] = -g_x_0_xxxxx_xzzzz[k] * cd_y[k] + g_x_0_xxxxx_xyzzzz[k];

                g_x_0_xxxxxy_yyyyy[k] = -g_x_0_xxxxx_yyyyy[k] * cd_y[k] + g_x_0_xxxxx_yyyyyy[k];

                g_x_0_xxxxxy_yyyyz[k] = -g_x_0_xxxxx_yyyyz[k] * cd_y[k] + g_x_0_xxxxx_yyyyyz[k];

                g_x_0_xxxxxy_yyyzz[k] = -g_x_0_xxxxx_yyyzz[k] * cd_y[k] + g_x_0_xxxxx_yyyyzz[k];

                g_x_0_xxxxxy_yyzzz[k] = -g_x_0_xxxxx_yyzzz[k] * cd_y[k] + g_x_0_xxxxx_yyyzzz[k];

                g_x_0_xxxxxy_yzzzz[k] = -g_x_0_xxxxx_yzzzz[k] * cd_y[k] + g_x_0_xxxxx_yyzzzz[k];

                g_x_0_xxxxxy_zzzzz[k] = -g_x_0_xxxxx_zzzzz[k] * cd_y[k] + g_x_0_xxxxx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxxz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxxz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxxxz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxxxz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxxxz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxxxz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxxxz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxxxz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxxxz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxxxz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxxxz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxxxz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxxxz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxxxz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxxxz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxxxz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxxxz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxxxz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxxxz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxxxz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_zzzzz, g_x_0_xxxxx_zzzzzz, g_x_0_xxxxxz_xxxxx, g_x_0_xxxxxz_xxxxy, g_x_0_xxxxxz_xxxxz, g_x_0_xxxxxz_xxxyy, g_x_0_xxxxxz_xxxyz, g_x_0_xxxxxz_xxxzz, g_x_0_xxxxxz_xxyyy, g_x_0_xxxxxz_xxyyz, g_x_0_xxxxxz_xxyzz, g_x_0_xxxxxz_xxzzz, g_x_0_xxxxxz_xyyyy, g_x_0_xxxxxz_xyyyz, g_x_0_xxxxxz_xyyzz, g_x_0_xxxxxz_xyzzz, g_x_0_xxxxxz_xzzzz, g_x_0_xxxxxz_yyyyy, g_x_0_xxxxxz_yyyyz, g_x_0_xxxxxz_yyyzz, g_x_0_xxxxxz_yyzzz, g_x_0_xxxxxz_yzzzz, g_x_0_xxxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_xxxxx[k] = -g_x_0_xxxxx_xxxxx[k] * cd_z[k] + g_x_0_xxxxx_xxxxxz[k];

                g_x_0_xxxxxz_xxxxy[k] = -g_x_0_xxxxx_xxxxy[k] * cd_z[k] + g_x_0_xxxxx_xxxxyz[k];

                g_x_0_xxxxxz_xxxxz[k] = -g_x_0_xxxxx_xxxxz[k] * cd_z[k] + g_x_0_xxxxx_xxxxzz[k];

                g_x_0_xxxxxz_xxxyy[k] = -g_x_0_xxxxx_xxxyy[k] * cd_z[k] + g_x_0_xxxxx_xxxyyz[k];

                g_x_0_xxxxxz_xxxyz[k] = -g_x_0_xxxxx_xxxyz[k] * cd_z[k] + g_x_0_xxxxx_xxxyzz[k];

                g_x_0_xxxxxz_xxxzz[k] = -g_x_0_xxxxx_xxxzz[k] * cd_z[k] + g_x_0_xxxxx_xxxzzz[k];

                g_x_0_xxxxxz_xxyyy[k] = -g_x_0_xxxxx_xxyyy[k] * cd_z[k] + g_x_0_xxxxx_xxyyyz[k];

                g_x_0_xxxxxz_xxyyz[k] = -g_x_0_xxxxx_xxyyz[k] * cd_z[k] + g_x_0_xxxxx_xxyyzz[k];

                g_x_0_xxxxxz_xxyzz[k] = -g_x_0_xxxxx_xxyzz[k] * cd_z[k] + g_x_0_xxxxx_xxyzzz[k];

                g_x_0_xxxxxz_xxzzz[k] = -g_x_0_xxxxx_xxzzz[k] * cd_z[k] + g_x_0_xxxxx_xxzzzz[k];

                g_x_0_xxxxxz_xyyyy[k] = -g_x_0_xxxxx_xyyyy[k] * cd_z[k] + g_x_0_xxxxx_xyyyyz[k];

                g_x_0_xxxxxz_xyyyz[k] = -g_x_0_xxxxx_xyyyz[k] * cd_z[k] + g_x_0_xxxxx_xyyyzz[k];

                g_x_0_xxxxxz_xyyzz[k] = -g_x_0_xxxxx_xyyzz[k] * cd_z[k] + g_x_0_xxxxx_xyyzzz[k];

                g_x_0_xxxxxz_xyzzz[k] = -g_x_0_xxxxx_xyzzz[k] * cd_z[k] + g_x_0_xxxxx_xyzzzz[k];

                g_x_0_xxxxxz_xzzzz[k] = -g_x_0_xxxxx_xzzzz[k] * cd_z[k] + g_x_0_xxxxx_xzzzzz[k];

                g_x_0_xxxxxz_yyyyy[k] = -g_x_0_xxxxx_yyyyy[k] * cd_z[k] + g_x_0_xxxxx_yyyyyz[k];

                g_x_0_xxxxxz_yyyyz[k] = -g_x_0_xxxxx_yyyyz[k] * cd_z[k] + g_x_0_xxxxx_yyyyzz[k];

                g_x_0_xxxxxz_yyyzz[k] = -g_x_0_xxxxx_yyyzz[k] * cd_z[k] + g_x_0_xxxxx_yyyzzz[k];

                g_x_0_xxxxxz_yyzzz[k] = -g_x_0_xxxxx_yyzzz[k] * cd_z[k] + g_x_0_xxxxx_yyzzzz[k];

                g_x_0_xxxxxz_yzzzz[k] = -g_x_0_xxxxx_yzzzz[k] * cd_z[k] + g_x_0_xxxxx_yzzzzz[k];

                g_x_0_xxxxxz_zzzzz[k] = -g_x_0_xxxxx_zzzzz[k] * cd_z[k] + g_x_0_xxxxx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxxyy_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxxyy_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxxyy_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxxyy_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxxyy_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxxyy_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxxyy_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxxyy_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxxyy_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxxyy_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxxyy_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxxyy_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxxyy_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxxyy_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxxyy_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxxyy_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxxyy_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxxyy_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxxyy_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxxyy_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxy_xxxxx, g_x_0_xxxxy_xxxxxy, g_x_0_xxxxy_xxxxy, g_x_0_xxxxy_xxxxyy, g_x_0_xxxxy_xxxxyz, g_x_0_xxxxy_xxxxz, g_x_0_xxxxy_xxxyy, g_x_0_xxxxy_xxxyyy, g_x_0_xxxxy_xxxyyz, g_x_0_xxxxy_xxxyz, g_x_0_xxxxy_xxxyzz, g_x_0_xxxxy_xxxzz, g_x_0_xxxxy_xxyyy, g_x_0_xxxxy_xxyyyy, g_x_0_xxxxy_xxyyyz, g_x_0_xxxxy_xxyyz, g_x_0_xxxxy_xxyyzz, g_x_0_xxxxy_xxyzz, g_x_0_xxxxy_xxyzzz, g_x_0_xxxxy_xxzzz, g_x_0_xxxxy_xyyyy, g_x_0_xxxxy_xyyyyy, g_x_0_xxxxy_xyyyyz, g_x_0_xxxxy_xyyyz, g_x_0_xxxxy_xyyyzz, g_x_0_xxxxy_xyyzz, g_x_0_xxxxy_xyyzzz, g_x_0_xxxxy_xyzzz, g_x_0_xxxxy_xyzzzz, g_x_0_xxxxy_xzzzz, g_x_0_xxxxy_yyyyy, g_x_0_xxxxy_yyyyyy, g_x_0_xxxxy_yyyyyz, g_x_0_xxxxy_yyyyz, g_x_0_xxxxy_yyyyzz, g_x_0_xxxxy_yyyzz, g_x_0_xxxxy_yyyzzz, g_x_0_xxxxy_yyzzz, g_x_0_xxxxy_yyzzzz, g_x_0_xxxxy_yzzzz, g_x_0_xxxxy_yzzzzz, g_x_0_xxxxy_zzzzz, g_x_0_xxxxyy_xxxxx, g_x_0_xxxxyy_xxxxy, g_x_0_xxxxyy_xxxxz, g_x_0_xxxxyy_xxxyy, g_x_0_xxxxyy_xxxyz, g_x_0_xxxxyy_xxxzz, g_x_0_xxxxyy_xxyyy, g_x_0_xxxxyy_xxyyz, g_x_0_xxxxyy_xxyzz, g_x_0_xxxxyy_xxzzz, g_x_0_xxxxyy_xyyyy, g_x_0_xxxxyy_xyyyz, g_x_0_xxxxyy_xyyzz, g_x_0_xxxxyy_xyzzz, g_x_0_xxxxyy_xzzzz, g_x_0_xxxxyy_yyyyy, g_x_0_xxxxyy_yyyyz, g_x_0_xxxxyy_yyyzz, g_x_0_xxxxyy_yyzzz, g_x_0_xxxxyy_yzzzz, g_x_0_xxxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_xxxxx[k] = -g_x_0_xxxxy_xxxxx[k] * cd_y[k] + g_x_0_xxxxy_xxxxxy[k];

                g_x_0_xxxxyy_xxxxy[k] = -g_x_0_xxxxy_xxxxy[k] * cd_y[k] + g_x_0_xxxxy_xxxxyy[k];

                g_x_0_xxxxyy_xxxxz[k] = -g_x_0_xxxxy_xxxxz[k] * cd_y[k] + g_x_0_xxxxy_xxxxyz[k];

                g_x_0_xxxxyy_xxxyy[k] = -g_x_0_xxxxy_xxxyy[k] * cd_y[k] + g_x_0_xxxxy_xxxyyy[k];

                g_x_0_xxxxyy_xxxyz[k] = -g_x_0_xxxxy_xxxyz[k] * cd_y[k] + g_x_0_xxxxy_xxxyyz[k];

                g_x_0_xxxxyy_xxxzz[k] = -g_x_0_xxxxy_xxxzz[k] * cd_y[k] + g_x_0_xxxxy_xxxyzz[k];

                g_x_0_xxxxyy_xxyyy[k] = -g_x_0_xxxxy_xxyyy[k] * cd_y[k] + g_x_0_xxxxy_xxyyyy[k];

                g_x_0_xxxxyy_xxyyz[k] = -g_x_0_xxxxy_xxyyz[k] * cd_y[k] + g_x_0_xxxxy_xxyyyz[k];

                g_x_0_xxxxyy_xxyzz[k] = -g_x_0_xxxxy_xxyzz[k] * cd_y[k] + g_x_0_xxxxy_xxyyzz[k];

                g_x_0_xxxxyy_xxzzz[k] = -g_x_0_xxxxy_xxzzz[k] * cd_y[k] + g_x_0_xxxxy_xxyzzz[k];

                g_x_0_xxxxyy_xyyyy[k] = -g_x_0_xxxxy_xyyyy[k] * cd_y[k] + g_x_0_xxxxy_xyyyyy[k];

                g_x_0_xxxxyy_xyyyz[k] = -g_x_0_xxxxy_xyyyz[k] * cd_y[k] + g_x_0_xxxxy_xyyyyz[k];

                g_x_0_xxxxyy_xyyzz[k] = -g_x_0_xxxxy_xyyzz[k] * cd_y[k] + g_x_0_xxxxy_xyyyzz[k];

                g_x_0_xxxxyy_xyzzz[k] = -g_x_0_xxxxy_xyzzz[k] * cd_y[k] + g_x_0_xxxxy_xyyzzz[k];

                g_x_0_xxxxyy_xzzzz[k] = -g_x_0_xxxxy_xzzzz[k] * cd_y[k] + g_x_0_xxxxy_xyzzzz[k];

                g_x_0_xxxxyy_yyyyy[k] = -g_x_0_xxxxy_yyyyy[k] * cd_y[k] + g_x_0_xxxxy_yyyyyy[k];

                g_x_0_xxxxyy_yyyyz[k] = -g_x_0_xxxxy_yyyyz[k] * cd_y[k] + g_x_0_xxxxy_yyyyyz[k];

                g_x_0_xxxxyy_yyyzz[k] = -g_x_0_xxxxy_yyyzz[k] * cd_y[k] + g_x_0_xxxxy_yyyyzz[k];

                g_x_0_xxxxyy_yyzzz[k] = -g_x_0_xxxxy_yyzzz[k] * cd_y[k] + g_x_0_xxxxy_yyyzzz[k];

                g_x_0_xxxxyy_yzzzz[k] = -g_x_0_xxxxy_yzzzz[k] * cd_y[k] + g_x_0_xxxxy_yyzzzz[k];

                g_x_0_xxxxyy_zzzzz[k] = -g_x_0_xxxxy_zzzzz[k] * cd_y[k] + g_x_0_xxxxy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxxyz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxxyz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxxyz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxxyz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxxyz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxxxyz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxxxyz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxxxyz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxxxyz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxxxyz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxxxyz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxxxyz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxxxyz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxxxyz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxxxyz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxxxyz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxxxyz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxxxyz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxxxyz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxxxyz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxyz_xxxxx, g_x_0_xxxxyz_xxxxy, g_x_0_xxxxyz_xxxxz, g_x_0_xxxxyz_xxxyy, g_x_0_xxxxyz_xxxyz, g_x_0_xxxxyz_xxxzz, g_x_0_xxxxyz_xxyyy, g_x_0_xxxxyz_xxyyz, g_x_0_xxxxyz_xxyzz, g_x_0_xxxxyz_xxzzz, g_x_0_xxxxyz_xyyyy, g_x_0_xxxxyz_xyyyz, g_x_0_xxxxyz_xyyzz, g_x_0_xxxxyz_xyzzz, g_x_0_xxxxyz_xzzzz, g_x_0_xxxxyz_yyyyy, g_x_0_xxxxyz_yyyyz, g_x_0_xxxxyz_yyyzz, g_x_0_xxxxyz_yyzzz, g_x_0_xxxxyz_yzzzz, g_x_0_xxxxyz_zzzzz, g_x_0_xxxxz_xxxxx, g_x_0_xxxxz_xxxxxy, g_x_0_xxxxz_xxxxy, g_x_0_xxxxz_xxxxyy, g_x_0_xxxxz_xxxxyz, g_x_0_xxxxz_xxxxz, g_x_0_xxxxz_xxxyy, g_x_0_xxxxz_xxxyyy, g_x_0_xxxxz_xxxyyz, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxyzz, g_x_0_xxxxz_xxxzz, g_x_0_xxxxz_xxyyy, g_x_0_xxxxz_xxyyyy, g_x_0_xxxxz_xxyyyz, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyyzz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxyzzz, g_x_0_xxxxz_xxzzz, g_x_0_xxxxz_xyyyy, g_x_0_xxxxz_xyyyyy, g_x_0_xxxxz_xyyyyz, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyyzz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyyzzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xyzzzz, g_x_0_xxxxz_xzzzz, g_x_0_xxxxz_yyyyy, g_x_0_xxxxz_yyyyyy, g_x_0_xxxxz_yyyyyz, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyyzz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyyzzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yyzzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_yzzzzz, g_x_0_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_xxxxx[k] = -g_x_0_xxxxz_xxxxx[k] * cd_y[k] + g_x_0_xxxxz_xxxxxy[k];

                g_x_0_xxxxyz_xxxxy[k] = -g_x_0_xxxxz_xxxxy[k] * cd_y[k] + g_x_0_xxxxz_xxxxyy[k];

                g_x_0_xxxxyz_xxxxz[k] = -g_x_0_xxxxz_xxxxz[k] * cd_y[k] + g_x_0_xxxxz_xxxxyz[k];

                g_x_0_xxxxyz_xxxyy[k] = -g_x_0_xxxxz_xxxyy[k] * cd_y[k] + g_x_0_xxxxz_xxxyyy[k];

                g_x_0_xxxxyz_xxxyz[k] = -g_x_0_xxxxz_xxxyz[k] * cd_y[k] + g_x_0_xxxxz_xxxyyz[k];

                g_x_0_xxxxyz_xxxzz[k] = -g_x_0_xxxxz_xxxzz[k] * cd_y[k] + g_x_0_xxxxz_xxxyzz[k];

                g_x_0_xxxxyz_xxyyy[k] = -g_x_0_xxxxz_xxyyy[k] * cd_y[k] + g_x_0_xxxxz_xxyyyy[k];

                g_x_0_xxxxyz_xxyyz[k] = -g_x_0_xxxxz_xxyyz[k] * cd_y[k] + g_x_0_xxxxz_xxyyyz[k];

                g_x_0_xxxxyz_xxyzz[k] = -g_x_0_xxxxz_xxyzz[k] * cd_y[k] + g_x_0_xxxxz_xxyyzz[k];

                g_x_0_xxxxyz_xxzzz[k] = -g_x_0_xxxxz_xxzzz[k] * cd_y[k] + g_x_0_xxxxz_xxyzzz[k];

                g_x_0_xxxxyz_xyyyy[k] = -g_x_0_xxxxz_xyyyy[k] * cd_y[k] + g_x_0_xxxxz_xyyyyy[k];

                g_x_0_xxxxyz_xyyyz[k] = -g_x_0_xxxxz_xyyyz[k] * cd_y[k] + g_x_0_xxxxz_xyyyyz[k];

                g_x_0_xxxxyz_xyyzz[k] = -g_x_0_xxxxz_xyyzz[k] * cd_y[k] + g_x_0_xxxxz_xyyyzz[k];

                g_x_0_xxxxyz_xyzzz[k] = -g_x_0_xxxxz_xyzzz[k] * cd_y[k] + g_x_0_xxxxz_xyyzzz[k];

                g_x_0_xxxxyz_xzzzz[k] = -g_x_0_xxxxz_xzzzz[k] * cd_y[k] + g_x_0_xxxxz_xyzzzz[k];

                g_x_0_xxxxyz_yyyyy[k] = -g_x_0_xxxxz_yyyyy[k] * cd_y[k] + g_x_0_xxxxz_yyyyyy[k];

                g_x_0_xxxxyz_yyyyz[k] = -g_x_0_xxxxz_yyyyz[k] * cd_y[k] + g_x_0_xxxxz_yyyyyz[k];

                g_x_0_xxxxyz_yyyzz[k] = -g_x_0_xxxxz_yyyzz[k] * cd_y[k] + g_x_0_xxxxz_yyyyzz[k];

                g_x_0_xxxxyz_yyzzz[k] = -g_x_0_xxxxz_yyzzz[k] * cd_y[k] + g_x_0_xxxxz_yyyzzz[k];

                g_x_0_xxxxyz_yzzzz[k] = -g_x_0_xxxxz_yzzzz[k] * cd_y[k] + g_x_0_xxxxz_yyzzzz[k];

                g_x_0_xxxxyz_zzzzz[k] = -g_x_0_xxxxz_zzzzz[k] * cd_y[k] + g_x_0_xxxxz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxxxzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxxxzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxxxzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxxxzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxxxzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxxxzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxxxzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxxxzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxxxzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxxxzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxxxzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxxxzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxxxzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxxxzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxxxzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxxxzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxxxzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxxxzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxxxzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxxxzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxz_xxxxx, g_x_0_xxxxz_xxxxxz, g_x_0_xxxxz_xxxxy, g_x_0_xxxxz_xxxxyz, g_x_0_xxxxz_xxxxz, g_x_0_xxxxz_xxxxzz, g_x_0_xxxxz_xxxyy, g_x_0_xxxxz_xxxyyz, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxyzz, g_x_0_xxxxz_xxxzz, g_x_0_xxxxz_xxxzzz, g_x_0_xxxxz_xxyyy, g_x_0_xxxxz_xxyyyz, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyyzz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxyzzz, g_x_0_xxxxz_xxzzz, g_x_0_xxxxz_xxzzzz, g_x_0_xxxxz_xyyyy, g_x_0_xxxxz_xyyyyz, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyyzz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyyzzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xyzzzz, g_x_0_xxxxz_xzzzz, g_x_0_xxxxz_xzzzzz, g_x_0_xxxxz_yyyyy, g_x_0_xxxxz_yyyyyz, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyyzz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyyzzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yyzzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_yzzzzz, g_x_0_xxxxz_zzzzz, g_x_0_xxxxz_zzzzzz, g_x_0_xxxxzz_xxxxx, g_x_0_xxxxzz_xxxxy, g_x_0_xxxxzz_xxxxz, g_x_0_xxxxzz_xxxyy, g_x_0_xxxxzz_xxxyz, g_x_0_xxxxzz_xxxzz, g_x_0_xxxxzz_xxyyy, g_x_0_xxxxzz_xxyyz, g_x_0_xxxxzz_xxyzz, g_x_0_xxxxzz_xxzzz, g_x_0_xxxxzz_xyyyy, g_x_0_xxxxzz_xyyyz, g_x_0_xxxxzz_xyyzz, g_x_0_xxxxzz_xyzzz, g_x_0_xxxxzz_xzzzz, g_x_0_xxxxzz_yyyyy, g_x_0_xxxxzz_yyyyz, g_x_0_xxxxzz_yyyzz, g_x_0_xxxxzz_yyzzz, g_x_0_xxxxzz_yzzzz, g_x_0_xxxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_xxxxx[k] = -g_x_0_xxxxz_xxxxx[k] * cd_z[k] + g_x_0_xxxxz_xxxxxz[k];

                g_x_0_xxxxzz_xxxxy[k] = -g_x_0_xxxxz_xxxxy[k] * cd_z[k] + g_x_0_xxxxz_xxxxyz[k];

                g_x_0_xxxxzz_xxxxz[k] = -g_x_0_xxxxz_xxxxz[k] * cd_z[k] + g_x_0_xxxxz_xxxxzz[k];

                g_x_0_xxxxzz_xxxyy[k] = -g_x_0_xxxxz_xxxyy[k] * cd_z[k] + g_x_0_xxxxz_xxxyyz[k];

                g_x_0_xxxxzz_xxxyz[k] = -g_x_0_xxxxz_xxxyz[k] * cd_z[k] + g_x_0_xxxxz_xxxyzz[k];

                g_x_0_xxxxzz_xxxzz[k] = -g_x_0_xxxxz_xxxzz[k] * cd_z[k] + g_x_0_xxxxz_xxxzzz[k];

                g_x_0_xxxxzz_xxyyy[k] = -g_x_0_xxxxz_xxyyy[k] * cd_z[k] + g_x_0_xxxxz_xxyyyz[k];

                g_x_0_xxxxzz_xxyyz[k] = -g_x_0_xxxxz_xxyyz[k] * cd_z[k] + g_x_0_xxxxz_xxyyzz[k];

                g_x_0_xxxxzz_xxyzz[k] = -g_x_0_xxxxz_xxyzz[k] * cd_z[k] + g_x_0_xxxxz_xxyzzz[k];

                g_x_0_xxxxzz_xxzzz[k] = -g_x_0_xxxxz_xxzzz[k] * cd_z[k] + g_x_0_xxxxz_xxzzzz[k];

                g_x_0_xxxxzz_xyyyy[k] = -g_x_0_xxxxz_xyyyy[k] * cd_z[k] + g_x_0_xxxxz_xyyyyz[k];

                g_x_0_xxxxzz_xyyyz[k] = -g_x_0_xxxxz_xyyyz[k] * cd_z[k] + g_x_0_xxxxz_xyyyzz[k];

                g_x_0_xxxxzz_xyyzz[k] = -g_x_0_xxxxz_xyyzz[k] * cd_z[k] + g_x_0_xxxxz_xyyzzz[k];

                g_x_0_xxxxzz_xyzzz[k] = -g_x_0_xxxxz_xyzzz[k] * cd_z[k] + g_x_0_xxxxz_xyzzzz[k];

                g_x_0_xxxxzz_xzzzz[k] = -g_x_0_xxxxz_xzzzz[k] * cd_z[k] + g_x_0_xxxxz_xzzzzz[k];

                g_x_0_xxxxzz_yyyyy[k] = -g_x_0_xxxxz_yyyyy[k] * cd_z[k] + g_x_0_xxxxz_yyyyyz[k];

                g_x_0_xxxxzz_yyyyz[k] = -g_x_0_xxxxz_yyyyz[k] * cd_z[k] + g_x_0_xxxxz_yyyyzz[k];

                g_x_0_xxxxzz_yyyzz[k] = -g_x_0_xxxxz_yyyzz[k] * cd_z[k] + g_x_0_xxxxz_yyyzzz[k];

                g_x_0_xxxxzz_yyzzz[k] = -g_x_0_xxxxz_yyzzz[k] * cd_z[k] + g_x_0_xxxxz_yyzzzz[k];

                g_x_0_xxxxzz_yzzzz[k] = -g_x_0_xxxxz_yzzzz[k] * cd_z[k] + g_x_0_xxxxz_yzzzzz[k];

                g_x_0_xxxxzz_zzzzz[k] = -g_x_0_xxxxz_zzzzz[k] * cd_z[k] + g_x_0_xxxxz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxxyyy_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxxyyy_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxxyyy_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxxyyy_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxxyyy_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxxyyy_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxxyyy_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxxyyy_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxxyyy_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxxyyy_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxxyyy_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxxyyy_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxxyyy_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxxyyy_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxxyyy_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxxyyy_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxxyyy_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxxyyy_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxxyyy_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxxyyy_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyy_xxxxx, g_x_0_xxxyy_xxxxxy, g_x_0_xxxyy_xxxxy, g_x_0_xxxyy_xxxxyy, g_x_0_xxxyy_xxxxyz, g_x_0_xxxyy_xxxxz, g_x_0_xxxyy_xxxyy, g_x_0_xxxyy_xxxyyy, g_x_0_xxxyy_xxxyyz, g_x_0_xxxyy_xxxyz, g_x_0_xxxyy_xxxyzz, g_x_0_xxxyy_xxxzz, g_x_0_xxxyy_xxyyy, g_x_0_xxxyy_xxyyyy, g_x_0_xxxyy_xxyyyz, g_x_0_xxxyy_xxyyz, g_x_0_xxxyy_xxyyzz, g_x_0_xxxyy_xxyzz, g_x_0_xxxyy_xxyzzz, g_x_0_xxxyy_xxzzz, g_x_0_xxxyy_xyyyy, g_x_0_xxxyy_xyyyyy, g_x_0_xxxyy_xyyyyz, g_x_0_xxxyy_xyyyz, g_x_0_xxxyy_xyyyzz, g_x_0_xxxyy_xyyzz, g_x_0_xxxyy_xyyzzz, g_x_0_xxxyy_xyzzz, g_x_0_xxxyy_xyzzzz, g_x_0_xxxyy_xzzzz, g_x_0_xxxyy_yyyyy, g_x_0_xxxyy_yyyyyy, g_x_0_xxxyy_yyyyyz, g_x_0_xxxyy_yyyyz, g_x_0_xxxyy_yyyyzz, g_x_0_xxxyy_yyyzz, g_x_0_xxxyy_yyyzzz, g_x_0_xxxyy_yyzzz, g_x_0_xxxyy_yyzzzz, g_x_0_xxxyy_yzzzz, g_x_0_xxxyy_yzzzzz, g_x_0_xxxyy_zzzzz, g_x_0_xxxyyy_xxxxx, g_x_0_xxxyyy_xxxxy, g_x_0_xxxyyy_xxxxz, g_x_0_xxxyyy_xxxyy, g_x_0_xxxyyy_xxxyz, g_x_0_xxxyyy_xxxzz, g_x_0_xxxyyy_xxyyy, g_x_0_xxxyyy_xxyyz, g_x_0_xxxyyy_xxyzz, g_x_0_xxxyyy_xxzzz, g_x_0_xxxyyy_xyyyy, g_x_0_xxxyyy_xyyyz, g_x_0_xxxyyy_xyyzz, g_x_0_xxxyyy_xyzzz, g_x_0_xxxyyy_xzzzz, g_x_0_xxxyyy_yyyyy, g_x_0_xxxyyy_yyyyz, g_x_0_xxxyyy_yyyzz, g_x_0_xxxyyy_yyzzz, g_x_0_xxxyyy_yzzzz, g_x_0_xxxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_xxxxx[k] = -g_x_0_xxxyy_xxxxx[k] * cd_y[k] + g_x_0_xxxyy_xxxxxy[k];

                g_x_0_xxxyyy_xxxxy[k] = -g_x_0_xxxyy_xxxxy[k] * cd_y[k] + g_x_0_xxxyy_xxxxyy[k];

                g_x_0_xxxyyy_xxxxz[k] = -g_x_0_xxxyy_xxxxz[k] * cd_y[k] + g_x_0_xxxyy_xxxxyz[k];

                g_x_0_xxxyyy_xxxyy[k] = -g_x_0_xxxyy_xxxyy[k] * cd_y[k] + g_x_0_xxxyy_xxxyyy[k];

                g_x_0_xxxyyy_xxxyz[k] = -g_x_0_xxxyy_xxxyz[k] * cd_y[k] + g_x_0_xxxyy_xxxyyz[k];

                g_x_0_xxxyyy_xxxzz[k] = -g_x_0_xxxyy_xxxzz[k] * cd_y[k] + g_x_0_xxxyy_xxxyzz[k];

                g_x_0_xxxyyy_xxyyy[k] = -g_x_0_xxxyy_xxyyy[k] * cd_y[k] + g_x_0_xxxyy_xxyyyy[k];

                g_x_0_xxxyyy_xxyyz[k] = -g_x_0_xxxyy_xxyyz[k] * cd_y[k] + g_x_0_xxxyy_xxyyyz[k];

                g_x_0_xxxyyy_xxyzz[k] = -g_x_0_xxxyy_xxyzz[k] * cd_y[k] + g_x_0_xxxyy_xxyyzz[k];

                g_x_0_xxxyyy_xxzzz[k] = -g_x_0_xxxyy_xxzzz[k] * cd_y[k] + g_x_0_xxxyy_xxyzzz[k];

                g_x_0_xxxyyy_xyyyy[k] = -g_x_0_xxxyy_xyyyy[k] * cd_y[k] + g_x_0_xxxyy_xyyyyy[k];

                g_x_0_xxxyyy_xyyyz[k] = -g_x_0_xxxyy_xyyyz[k] * cd_y[k] + g_x_0_xxxyy_xyyyyz[k];

                g_x_0_xxxyyy_xyyzz[k] = -g_x_0_xxxyy_xyyzz[k] * cd_y[k] + g_x_0_xxxyy_xyyyzz[k];

                g_x_0_xxxyyy_xyzzz[k] = -g_x_0_xxxyy_xyzzz[k] * cd_y[k] + g_x_0_xxxyy_xyyzzz[k];

                g_x_0_xxxyyy_xzzzz[k] = -g_x_0_xxxyy_xzzzz[k] * cd_y[k] + g_x_0_xxxyy_xyzzzz[k];

                g_x_0_xxxyyy_yyyyy[k] = -g_x_0_xxxyy_yyyyy[k] * cd_y[k] + g_x_0_xxxyy_yyyyyy[k];

                g_x_0_xxxyyy_yyyyz[k] = -g_x_0_xxxyy_yyyyz[k] * cd_y[k] + g_x_0_xxxyy_yyyyyz[k];

                g_x_0_xxxyyy_yyyzz[k] = -g_x_0_xxxyy_yyyzz[k] * cd_y[k] + g_x_0_xxxyy_yyyyzz[k];

                g_x_0_xxxyyy_yyzzz[k] = -g_x_0_xxxyy_yyzzz[k] * cd_y[k] + g_x_0_xxxyy_yyyzzz[k];

                g_x_0_xxxyyy_yzzzz[k] = -g_x_0_xxxyy_yzzzz[k] * cd_y[k] + g_x_0_xxxyy_yyzzzz[k];

                g_x_0_xxxyyy_zzzzz[k] = -g_x_0_xxxyy_zzzzz[k] * cd_y[k] + g_x_0_xxxyy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxxyyz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxxyyz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxxyyz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxxyyz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxxyyz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxxyyz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxxyyz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxxyyz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxxyyz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxxyyz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxxyyz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxxyyz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxxyyz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxxyyz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxxyyz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxxyyz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxxyyz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxxyyz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxxyyz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxxyyz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyyz_xxxxx, g_x_0_xxxyyz_xxxxy, g_x_0_xxxyyz_xxxxz, g_x_0_xxxyyz_xxxyy, g_x_0_xxxyyz_xxxyz, g_x_0_xxxyyz_xxxzz, g_x_0_xxxyyz_xxyyy, g_x_0_xxxyyz_xxyyz, g_x_0_xxxyyz_xxyzz, g_x_0_xxxyyz_xxzzz, g_x_0_xxxyyz_xyyyy, g_x_0_xxxyyz_xyyyz, g_x_0_xxxyyz_xyyzz, g_x_0_xxxyyz_xyzzz, g_x_0_xxxyyz_xzzzz, g_x_0_xxxyyz_yyyyy, g_x_0_xxxyyz_yyyyz, g_x_0_xxxyyz_yyyzz, g_x_0_xxxyyz_yyzzz, g_x_0_xxxyyz_yzzzz, g_x_0_xxxyyz_zzzzz, g_x_0_xxxyz_xxxxx, g_x_0_xxxyz_xxxxxy, g_x_0_xxxyz_xxxxy, g_x_0_xxxyz_xxxxyy, g_x_0_xxxyz_xxxxyz, g_x_0_xxxyz_xxxxz, g_x_0_xxxyz_xxxyy, g_x_0_xxxyz_xxxyyy, g_x_0_xxxyz_xxxyyz, g_x_0_xxxyz_xxxyz, g_x_0_xxxyz_xxxyzz, g_x_0_xxxyz_xxxzz, g_x_0_xxxyz_xxyyy, g_x_0_xxxyz_xxyyyy, g_x_0_xxxyz_xxyyyz, g_x_0_xxxyz_xxyyz, g_x_0_xxxyz_xxyyzz, g_x_0_xxxyz_xxyzz, g_x_0_xxxyz_xxyzzz, g_x_0_xxxyz_xxzzz, g_x_0_xxxyz_xyyyy, g_x_0_xxxyz_xyyyyy, g_x_0_xxxyz_xyyyyz, g_x_0_xxxyz_xyyyz, g_x_0_xxxyz_xyyyzz, g_x_0_xxxyz_xyyzz, g_x_0_xxxyz_xyyzzz, g_x_0_xxxyz_xyzzz, g_x_0_xxxyz_xyzzzz, g_x_0_xxxyz_xzzzz, g_x_0_xxxyz_yyyyy, g_x_0_xxxyz_yyyyyy, g_x_0_xxxyz_yyyyyz, g_x_0_xxxyz_yyyyz, g_x_0_xxxyz_yyyyzz, g_x_0_xxxyz_yyyzz, g_x_0_xxxyz_yyyzzz, g_x_0_xxxyz_yyzzz, g_x_0_xxxyz_yyzzzz, g_x_0_xxxyz_yzzzz, g_x_0_xxxyz_yzzzzz, g_x_0_xxxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_xxxxx[k] = -g_x_0_xxxyz_xxxxx[k] * cd_y[k] + g_x_0_xxxyz_xxxxxy[k];

                g_x_0_xxxyyz_xxxxy[k] = -g_x_0_xxxyz_xxxxy[k] * cd_y[k] + g_x_0_xxxyz_xxxxyy[k];

                g_x_0_xxxyyz_xxxxz[k] = -g_x_0_xxxyz_xxxxz[k] * cd_y[k] + g_x_0_xxxyz_xxxxyz[k];

                g_x_0_xxxyyz_xxxyy[k] = -g_x_0_xxxyz_xxxyy[k] * cd_y[k] + g_x_0_xxxyz_xxxyyy[k];

                g_x_0_xxxyyz_xxxyz[k] = -g_x_0_xxxyz_xxxyz[k] * cd_y[k] + g_x_0_xxxyz_xxxyyz[k];

                g_x_0_xxxyyz_xxxzz[k] = -g_x_0_xxxyz_xxxzz[k] * cd_y[k] + g_x_0_xxxyz_xxxyzz[k];

                g_x_0_xxxyyz_xxyyy[k] = -g_x_0_xxxyz_xxyyy[k] * cd_y[k] + g_x_0_xxxyz_xxyyyy[k];

                g_x_0_xxxyyz_xxyyz[k] = -g_x_0_xxxyz_xxyyz[k] * cd_y[k] + g_x_0_xxxyz_xxyyyz[k];

                g_x_0_xxxyyz_xxyzz[k] = -g_x_0_xxxyz_xxyzz[k] * cd_y[k] + g_x_0_xxxyz_xxyyzz[k];

                g_x_0_xxxyyz_xxzzz[k] = -g_x_0_xxxyz_xxzzz[k] * cd_y[k] + g_x_0_xxxyz_xxyzzz[k];

                g_x_0_xxxyyz_xyyyy[k] = -g_x_0_xxxyz_xyyyy[k] * cd_y[k] + g_x_0_xxxyz_xyyyyy[k];

                g_x_0_xxxyyz_xyyyz[k] = -g_x_0_xxxyz_xyyyz[k] * cd_y[k] + g_x_0_xxxyz_xyyyyz[k];

                g_x_0_xxxyyz_xyyzz[k] = -g_x_0_xxxyz_xyyzz[k] * cd_y[k] + g_x_0_xxxyz_xyyyzz[k];

                g_x_0_xxxyyz_xyzzz[k] = -g_x_0_xxxyz_xyzzz[k] * cd_y[k] + g_x_0_xxxyz_xyyzzz[k];

                g_x_0_xxxyyz_xzzzz[k] = -g_x_0_xxxyz_xzzzz[k] * cd_y[k] + g_x_0_xxxyz_xyzzzz[k];

                g_x_0_xxxyyz_yyyyy[k] = -g_x_0_xxxyz_yyyyy[k] * cd_y[k] + g_x_0_xxxyz_yyyyyy[k];

                g_x_0_xxxyyz_yyyyz[k] = -g_x_0_xxxyz_yyyyz[k] * cd_y[k] + g_x_0_xxxyz_yyyyyz[k];

                g_x_0_xxxyyz_yyyzz[k] = -g_x_0_xxxyz_yyyzz[k] * cd_y[k] + g_x_0_xxxyz_yyyyzz[k];

                g_x_0_xxxyyz_yyzzz[k] = -g_x_0_xxxyz_yyzzz[k] * cd_y[k] + g_x_0_xxxyz_yyyzzz[k];

                g_x_0_xxxyyz_yzzzz[k] = -g_x_0_xxxyz_yzzzz[k] * cd_y[k] + g_x_0_xxxyz_yyzzzz[k];

                g_x_0_xxxyyz_zzzzz[k] = -g_x_0_xxxyz_zzzzz[k] * cd_y[k] + g_x_0_xxxyz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xxxyzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xxxyzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xxxyzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xxxyzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xxxyzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xxxyzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xxxyzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xxxyzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xxxyzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xxxyzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xxxyzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xxxyzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xxxyzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xxxyzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xxxyzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xxxyzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xxxyzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xxxyzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xxxyzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xxxyzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyzz_xxxxx, g_x_0_xxxyzz_xxxxy, g_x_0_xxxyzz_xxxxz, g_x_0_xxxyzz_xxxyy, g_x_0_xxxyzz_xxxyz, g_x_0_xxxyzz_xxxzz, g_x_0_xxxyzz_xxyyy, g_x_0_xxxyzz_xxyyz, g_x_0_xxxyzz_xxyzz, g_x_0_xxxyzz_xxzzz, g_x_0_xxxyzz_xyyyy, g_x_0_xxxyzz_xyyyz, g_x_0_xxxyzz_xyyzz, g_x_0_xxxyzz_xyzzz, g_x_0_xxxyzz_xzzzz, g_x_0_xxxyzz_yyyyy, g_x_0_xxxyzz_yyyyz, g_x_0_xxxyzz_yyyzz, g_x_0_xxxyzz_yyzzz, g_x_0_xxxyzz_yzzzz, g_x_0_xxxyzz_zzzzz, g_x_0_xxxzz_xxxxx, g_x_0_xxxzz_xxxxxy, g_x_0_xxxzz_xxxxy, g_x_0_xxxzz_xxxxyy, g_x_0_xxxzz_xxxxyz, g_x_0_xxxzz_xxxxz, g_x_0_xxxzz_xxxyy, g_x_0_xxxzz_xxxyyy, g_x_0_xxxzz_xxxyyz, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxyzz, g_x_0_xxxzz_xxxzz, g_x_0_xxxzz_xxyyy, g_x_0_xxxzz_xxyyyy, g_x_0_xxxzz_xxyyyz, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyyzz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxyzzz, g_x_0_xxxzz_xxzzz, g_x_0_xxxzz_xyyyy, g_x_0_xxxzz_xyyyyy, g_x_0_xxxzz_xyyyyz, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyyzz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyyzzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xyzzzz, g_x_0_xxxzz_xzzzz, g_x_0_xxxzz_yyyyy, g_x_0_xxxzz_yyyyyy, g_x_0_xxxzz_yyyyyz, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyyzz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyyzzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yyzzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_yzzzzz, g_x_0_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_xxxxx[k] = -g_x_0_xxxzz_xxxxx[k] * cd_y[k] + g_x_0_xxxzz_xxxxxy[k];

                g_x_0_xxxyzz_xxxxy[k] = -g_x_0_xxxzz_xxxxy[k] * cd_y[k] + g_x_0_xxxzz_xxxxyy[k];

                g_x_0_xxxyzz_xxxxz[k] = -g_x_0_xxxzz_xxxxz[k] * cd_y[k] + g_x_0_xxxzz_xxxxyz[k];

                g_x_0_xxxyzz_xxxyy[k] = -g_x_0_xxxzz_xxxyy[k] * cd_y[k] + g_x_0_xxxzz_xxxyyy[k];

                g_x_0_xxxyzz_xxxyz[k] = -g_x_0_xxxzz_xxxyz[k] * cd_y[k] + g_x_0_xxxzz_xxxyyz[k];

                g_x_0_xxxyzz_xxxzz[k] = -g_x_0_xxxzz_xxxzz[k] * cd_y[k] + g_x_0_xxxzz_xxxyzz[k];

                g_x_0_xxxyzz_xxyyy[k] = -g_x_0_xxxzz_xxyyy[k] * cd_y[k] + g_x_0_xxxzz_xxyyyy[k];

                g_x_0_xxxyzz_xxyyz[k] = -g_x_0_xxxzz_xxyyz[k] * cd_y[k] + g_x_0_xxxzz_xxyyyz[k];

                g_x_0_xxxyzz_xxyzz[k] = -g_x_0_xxxzz_xxyzz[k] * cd_y[k] + g_x_0_xxxzz_xxyyzz[k];

                g_x_0_xxxyzz_xxzzz[k] = -g_x_0_xxxzz_xxzzz[k] * cd_y[k] + g_x_0_xxxzz_xxyzzz[k];

                g_x_0_xxxyzz_xyyyy[k] = -g_x_0_xxxzz_xyyyy[k] * cd_y[k] + g_x_0_xxxzz_xyyyyy[k];

                g_x_0_xxxyzz_xyyyz[k] = -g_x_0_xxxzz_xyyyz[k] * cd_y[k] + g_x_0_xxxzz_xyyyyz[k];

                g_x_0_xxxyzz_xyyzz[k] = -g_x_0_xxxzz_xyyzz[k] * cd_y[k] + g_x_0_xxxzz_xyyyzz[k];

                g_x_0_xxxyzz_xyzzz[k] = -g_x_0_xxxzz_xyzzz[k] * cd_y[k] + g_x_0_xxxzz_xyyzzz[k];

                g_x_0_xxxyzz_xzzzz[k] = -g_x_0_xxxzz_xzzzz[k] * cd_y[k] + g_x_0_xxxzz_xyzzzz[k];

                g_x_0_xxxyzz_yyyyy[k] = -g_x_0_xxxzz_yyyyy[k] * cd_y[k] + g_x_0_xxxzz_yyyyyy[k];

                g_x_0_xxxyzz_yyyyz[k] = -g_x_0_xxxzz_yyyyz[k] * cd_y[k] + g_x_0_xxxzz_yyyyyz[k];

                g_x_0_xxxyzz_yyyzz[k] = -g_x_0_xxxzz_yyyzz[k] * cd_y[k] + g_x_0_xxxzz_yyyyzz[k];

                g_x_0_xxxyzz_yyzzz[k] = -g_x_0_xxxzz_yyzzz[k] * cd_y[k] + g_x_0_xxxzz_yyyzzz[k];

                g_x_0_xxxyzz_yzzzz[k] = -g_x_0_xxxzz_yzzzz[k] * cd_y[k] + g_x_0_xxxzz_yyzzzz[k];

                g_x_0_xxxyzz_zzzzz[k] = -g_x_0_xxxzz_zzzzz[k] * cd_y[k] + g_x_0_xxxzz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xxxzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xxxzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xxxzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xxxzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xxxzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xxxzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xxxzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xxxzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xxxzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xxxzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xxxzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xxxzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xxxzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xxxzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xxxzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xxxzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xxxzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xxxzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xxxzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xxxzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_z, g_x_0_xxxzz_xxxxx, g_x_0_xxxzz_xxxxxz, g_x_0_xxxzz_xxxxy, g_x_0_xxxzz_xxxxyz, g_x_0_xxxzz_xxxxz, g_x_0_xxxzz_xxxxzz, g_x_0_xxxzz_xxxyy, g_x_0_xxxzz_xxxyyz, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxyzz, g_x_0_xxxzz_xxxzz, g_x_0_xxxzz_xxxzzz, g_x_0_xxxzz_xxyyy, g_x_0_xxxzz_xxyyyz, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyyzz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxyzzz, g_x_0_xxxzz_xxzzz, g_x_0_xxxzz_xxzzzz, g_x_0_xxxzz_xyyyy, g_x_0_xxxzz_xyyyyz, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyyzz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyyzzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xyzzzz, g_x_0_xxxzz_xzzzz, g_x_0_xxxzz_xzzzzz, g_x_0_xxxzz_yyyyy, g_x_0_xxxzz_yyyyyz, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyyzz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyyzzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yyzzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_yzzzzz, g_x_0_xxxzz_zzzzz, g_x_0_xxxzz_zzzzzz, g_x_0_xxxzzz_xxxxx, g_x_0_xxxzzz_xxxxy, g_x_0_xxxzzz_xxxxz, g_x_0_xxxzzz_xxxyy, g_x_0_xxxzzz_xxxyz, g_x_0_xxxzzz_xxxzz, g_x_0_xxxzzz_xxyyy, g_x_0_xxxzzz_xxyyz, g_x_0_xxxzzz_xxyzz, g_x_0_xxxzzz_xxzzz, g_x_0_xxxzzz_xyyyy, g_x_0_xxxzzz_xyyyz, g_x_0_xxxzzz_xyyzz, g_x_0_xxxzzz_xyzzz, g_x_0_xxxzzz_xzzzz, g_x_0_xxxzzz_yyyyy, g_x_0_xxxzzz_yyyyz, g_x_0_xxxzzz_yyyzz, g_x_0_xxxzzz_yyzzz, g_x_0_xxxzzz_yzzzz, g_x_0_xxxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_xxxxx[k] = -g_x_0_xxxzz_xxxxx[k] * cd_z[k] + g_x_0_xxxzz_xxxxxz[k];

                g_x_0_xxxzzz_xxxxy[k] = -g_x_0_xxxzz_xxxxy[k] * cd_z[k] + g_x_0_xxxzz_xxxxyz[k];

                g_x_0_xxxzzz_xxxxz[k] = -g_x_0_xxxzz_xxxxz[k] * cd_z[k] + g_x_0_xxxzz_xxxxzz[k];

                g_x_0_xxxzzz_xxxyy[k] = -g_x_0_xxxzz_xxxyy[k] * cd_z[k] + g_x_0_xxxzz_xxxyyz[k];

                g_x_0_xxxzzz_xxxyz[k] = -g_x_0_xxxzz_xxxyz[k] * cd_z[k] + g_x_0_xxxzz_xxxyzz[k];

                g_x_0_xxxzzz_xxxzz[k] = -g_x_0_xxxzz_xxxzz[k] * cd_z[k] + g_x_0_xxxzz_xxxzzz[k];

                g_x_0_xxxzzz_xxyyy[k] = -g_x_0_xxxzz_xxyyy[k] * cd_z[k] + g_x_0_xxxzz_xxyyyz[k];

                g_x_0_xxxzzz_xxyyz[k] = -g_x_0_xxxzz_xxyyz[k] * cd_z[k] + g_x_0_xxxzz_xxyyzz[k];

                g_x_0_xxxzzz_xxyzz[k] = -g_x_0_xxxzz_xxyzz[k] * cd_z[k] + g_x_0_xxxzz_xxyzzz[k];

                g_x_0_xxxzzz_xxzzz[k] = -g_x_0_xxxzz_xxzzz[k] * cd_z[k] + g_x_0_xxxzz_xxzzzz[k];

                g_x_0_xxxzzz_xyyyy[k] = -g_x_0_xxxzz_xyyyy[k] * cd_z[k] + g_x_0_xxxzz_xyyyyz[k];

                g_x_0_xxxzzz_xyyyz[k] = -g_x_0_xxxzz_xyyyz[k] * cd_z[k] + g_x_0_xxxzz_xyyyzz[k];

                g_x_0_xxxzzz_xyyzz[k] = -g_x_0_xxxzz_xyyzz[k] * cd_z[k] + g_x_0_xxxzz_xyyzzz[k];

                g_x_0_xxxzzz_xyzzz[k] = -g_x_0_xxxzz_xyzzz[k] * cd_z[k] + g_x_0_xxxzz_xyzzzz[k];

                g_x_0_xxxzzz_xzzzz[k] = -g_x_0_xxxzz_xzzzz[k] * cd_z[k] + g_x_0_xxxzz_xzzzzz[k];

                g_x_0_xxxzzz_yyyyy[k] = -g_x_0_xxxzz_yyyyy[k] * cd_z[k] + g_x_0_xxxzz_yyyyyz[k];

                g_x_0_xxxzzz_yyyyz[k] = -g_x_0_xxxzz_yyyyz[k] * cd_z[k] + g_x_0_xxxzz_yyyyzz[k];

                g_x_0_xxxzzz_yyyzz[k] = -g_x_0_xxxzz_yyyzz[k] * cd_z[k] + g_x_0_xxxzz_yyyzzz[k];

                g_x_0_xxxzzz_yyzzz[k] = -g_x_0_xxxzz_yyzzz[k] * cd_z[k] + g_x_0_xxxzz_yyzzzz[k];

                g_x_0_xxxzzz_yzzzz[k] = -g_x_0_xxxzz_yzzzz[k] * cd_z[k] + g_x_0_xxxzz_yzzzzz[k];

                g_x_0_xxxzzz_zzzzz[k] = -g_x_0_xxxzz_zzzzz[k] * cd_z[k] + g_x_0_xxxzz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xxyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xxyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xxyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xxyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xxyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xxyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xxyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xxyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xxyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xxyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xxyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xxyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xxyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xxyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xxyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xxyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xxyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xxyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xxyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xxyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyy_xxxxx, g_x_0_xxyyy_xxxxxy, g_x_0_xxyyy_xxxxy, g_x_0_xxyyy_xxxxyy, g_x_0_xxyyy_xxxxyz, g_x_0_xxyyy_xxxxz, g_x_0_xxyyy_xxxyy, g_x_0_xxyyy_xxxyyy, g_x_0_xxyyy_xxxyyz, g_x_0_xxyyy_xxxyz, g_x_0_xxyyy_xxxyzz, g_x_0_xxyyy_xxxzz, g_x_0_xxyyy_xxyyy, g_x_0_xxyyy_xxyyyy, g_x_0_xxyyy_xxyyyz, g_x_0_xxyyy_xxyyz, g_x_0_xxyyy_xxyyzz, g_x_0_xxyyy_xxyzz, g_x_0_xxyyy_xxyzzz, g_x_0_xxyyy_xxzzz, g_x_0_xxyyy_xyyyy, g_x_0_xxyyy_xyyyyy, g_x_0_xxyyy_xyyyyz, g_x_0_xxyyy_xyyyz, g_x_0_xxyyy_xyyyzz, g_x_0_xxyyy_xyyzz, g_x_0_xxyyy_xyyzzz, g_x_0_xxyyy_xyzzz, g_x_0_xxyyy_xyzzzz, g_x_0_xxyyy_xzzzz, g_x_0_xxyyy_yyyyy, g_x_0_xxyyy_yyyyyy, g_x_0_xxyyy_yyyyyz, g_x_0_xxyyy_yyyyz, g_x_0_xxyyy_yyyyzz, g_x_0_xxyyy_yyyzz, g_x_0_xxyyy_yyyzzz, g_x_0_xxyyy_yyzzz, g_x_0_xxyyy_yyzzzz, g_x_0_xxyyy_yzzzz, g_x_0_xxyyy_yzzzzz, g_x_0_xxyyy_zzzzz, g_x_0_xxyyyy_xxxxx, g_x_0_xxyyyy_xxxxy, g_x_0_xxyyyy_xxxxz, g_x_0_xxyyyy_xxxyy, g_x_0_xxyyyy_xxxyz, g_x_0_xxyyyy_xxxzz, g_x_0_xxyyyy_xxyyy, g_x_0_xxyyyy_xxyyz, g_x_0_xxyyyy_xxyzz, g_x_0_xxyyyy_xxzzz, g_x_0_xxyyyy_xyyyy, g_x_0_xxyyyy_xyyyz, g_x_0_xxyyyy_xyyzz, g_x_0_xxyyyy_xyzzz, g_x_0_xxyyyy_xzzzz, g_x_0_xxyyyy_yyyyy, g_x_0_xxyyyy_yyyyz, g_x_0_xxyyyy_yyyzz, g_x_0_xxyyyy_yyzzz, g_x_0_xxyyyy_yzzzz, g_x_0_xxyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_xxxxx[k] = -g_x_0_xxyyy_xxxxx[k] * cd_y[k] + g_x_0_xxyyy_xxxxxy[k];

                g_x_0_xxyyyy_xxxxy[k] = -g_x_0_xxyyy_xxxxy[k] * cd_y[k] + g_x_0_xxyyy_xxxxyy[k];

                g_x_0_xxyyyy_xxxxz[k] = -g_x_0_xxyyy_xxxxz[k] * cd_y[k] + g_x_0_xxyyy_xxxxyz[k];

                g_x_0_xxyyyy_xxxyy[k] = -g_x_0_xxyyy_xxxyy[k] * cd_y[k] + g_x_0_xxyyy_xxxyyy[k];

                g_x_0_xxyyyy_xxxyz[k] = -g_x_0_xxyyy_xxxyz[k] * cd_y[k] + g_x_0_xxyyy_xxxyyz[k];

                g_x_0_xxyyyy_xxxzz[k] = -g_x_0_xxyyy_xxxzz[k] * cd_y[k] + g_x_0_xxyyy_xxxyzz[k];

                g_x_0_xxyyyy_xxyyy[k] = -g_x_0_xxyyy_xxyyy[k] * cd_y[k] + g_x_0_xxyyy_xxyyyy[k];

                g_x_0_xxyyyy_xxyyz[k] = -g_x_0_xxyyy_xxyyz[k] * cd_y[k] + g_x_0_xxyyy_xxyyyz[k];

                g_x_0_xxyyyy_xxyzz[k] = -g_x_0_xxyyy_xxyzz[k] * cd_y[k] + g_x_0_xxyyy_xxyyzz[k];

                g_x_0_xxyyyy_xxzzz[k] = -g_x_0_xxyyy_xxzzz[k] * cd_y[k] + g_x_0_xxyyy_xxyzzz[k];

                g_x_0_xxyyyy_xyyyy[k] = -g_x_0_xxyyy_xyyyy[k] * cd_y[k] + g_x_0_xxyyy_xyyyyy[k];

                g_x_0_xxyyyy_xyyyz[k] = -g_x_0_xxyyy_xyyyz[k] * cd_y[k] + g_x_0_xxyyy_xyyyyz[k];

                g_x_0_xxyyyy_xyyzz[k] = -g_x_0_xxyyy_xyyzz[k] * cd_y[k] + g_x_0_xxyyy_xyyyzz[k];

                g_x_0_xxyyyy_xyzzz[k] = -g_x_0_xxyyy_xyzzz[k] * cd_y[k] + g_x_0_xxyyy_xyyzzz[k];

                g_x_0_xxyyyy_xzzzz[k] = -g_x_0_xxyyy_xzzzz[k] * cd_y[k] + g_x_0_xxyyy_xyzzzz[k];

                g_x_0_xxyyyy_yyyyy[k] = -g_x_0_xxyyy_yyyyy[k] * cd_y[k] + g_x_0_xxyyy_yyyyyy[k];

                g_x_0_xxyyyy_yyyyz[k] = -g_x_0_xxyyy_yyyyz[k] * cd_y[k] + g_x_0_xxyyy_yyyyyz[k];

                g_x_0_xxyyyy_yyyzz[k] = -g_x_0_xxyyy_yyyzz[k] * cd_y[k] + g_x_0_xxyyy_yyyyzz[k];

                g_x_0_xxyyyy_yyzzz[k] = -g_x_0_xxyyy_yyzzz[k] * cd_y[k] + g_x_0_xxyyy_yyyzzz[k];

                g_x_0_xxyyyy_yzzzz[k] = -g_x_0_xxyyy_yzzzz[k] * cd_y[k] + g_x_0_xxyyy_yyzzzz[k];

                g_x_0_xxyyyy_zzzzz[k] = -g_x_0_xxyyy_zzzzz[k] * cd_y[k] + g_x_0_xxyyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xxyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xxyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xxyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xxyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xxyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xxyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xxyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xxyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xxyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xxyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xxyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xxyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xxyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xxyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xxyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xxyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xxyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xxyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xxyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xxyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyyz_xxxxx, g_x_0_xxyyyz_xxxxy, g_x_0_xxyyyz_xxxxz, g_x_0_xxyyyz_xxxyy, g_x_0_xxyyyz_xxxyz, g_x_0_xxyyyz_xxxzz, g_x_0_xxyyyz_xxyyy, g_x_0_xxyyyz_xxyyz, g_x_0_xxyyyz_xxyzz, g_x_0_xxyyyz_xxzzz, g_x_0_xxyyyz_xyyyy, g_x_0_xxyyyz_xyyyz, g_x_0_xxyyyz_xyyzz, g_x_0_xxyyyz_xyzzz, g_x_0_xxyyyz_xzzzz, g_x_0_xxyyyz_yyyyy, g_x_0_xxyyyz_yyyyz, g_x_0_xxyyyz_yyyzz, g_x_0_xxyyyz_yyzzz, g_x_0_xxyyyz_yzzzz, g_x_0_xxyyyz_zzzzz, g_x_0_xxyyz_xxxxx, g_x_0_xxyyz_xxxxxy, g_x_0_xxyyz_xxxxy, g_x_0_xxyyz_xxxxyy, g_x_0_xxyyz_xxxxyz, g_x_0_xxyyz_xxxxz, g_x_0_xxyyz_xxxyy, g_x_0_xxyyz_xxxyyy, g_x_0_xxyyz_xxxyyz, g_x_0_xxyyz_xxxyz, g_x_0_xxyyz_xxxyzz, g_x_0_xxyyz_xxxzz, g_x_0_xxyyz_xxyyy, g_x_0_xxyyz_xxyyyy, g_x_0_xxyyz_xxyyyz, g_x_0_xxyyz_xxyyz, g_x_0_xxyyz_xxyyzz, g_x_0_xxyyz_xxyzz, g_x_0_xxyyz_xxyzzz, g_x_0_xxyyz_xxzzz, g_x_0_xxyyz_xyyyy, g_x_0_xxyyz_xyyyyy, g_x_0_xxyyz_xyyyyz, g_x_0_xxyyz_xyyyz, g_x_0_xxyyz_xyyyzz, g_x_0_xxyyz_xyyzz, g_x_0_xxyyz_xyyzzz, g_x_0_xxyyz_xyzzz, g_x_0_xxyyz_xyzzzz, g_x_0_xxyyz_xzzzz, g_x_0_xxyyz_yyyyy, g_x_0_xxyyz_yyyyyy, g_x_0_xxyyz_yyyyyz, g_x_0_xxyyz_yyyyz, g_x_0_xxyyz_yyyyzz, g_x_0_xxyyz_yyyzz, g_x_0_xxyyz_yyyzzz, g_x_0_xxyyz_yyzzz, g_x_0_xxyyz_yyzzzz, g_x_0_xxyyz_yzzzz, g_x_0_xxyyz_yzzzzz, g_x_0_xxyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_xxxxx[k] = -g_x_0_xxyyz_xxxxx[k] * cd_y[k] + g_x_0_xxyyz_xxxxxy[k];

                g_x_0_xxyyyz_xxxxy[k] = -g_x_0_xxyyz_xxxxy[k] * cd_y[k] + g_x_0_xxyyz_xxxxyy[k];

                g_x_0_xxyyyz_xxxxz[k] = -g_x_0_xxyyz_xxxxz[k] * cd_y[k] + g_x_0_xxyyz_xxxxyz[k];

                g_x_0_xxyyyz_xxxyy[k] = -g_x_0_xxyyz_xxxyy[k] * cd_y[k] + g_x_0_xxyyz_xxxyyy[k];

                g_x_0_xxyyyz_xxxyz[k] = -g_x_0_xxyyz_xxxyz[k] * cd_y[k] + g_x_0_xxyyz_xxxyyz[k];

                g_x_0_xxyyyz_xxxzz[k] = -g_x_0_xxyyz_xxxzz[k] * cd_y[k] + g_x_0_xxyyz_xxxyzz[k];

                g_x_0_xxyyyz_xxyyy[k] = -g_x_0_xxyyz_xxyyy[k] * cd_y[k] + g_x_0_xxyyz_xxyyyy[k];

                g_x_0_xxyyyz_xxyyz[k] = -g_x_0_xxyyz_xxyyz[k] * cd_y[k] + g_x_0_xxyyz_xxyyyz[k];

                g_x_0_xxyyyz_xxyzz[k] = -g_x_0_xxyyz_xxyzz[k] * cd_y[k] + g_x_0_xxyyz_xxyyzz[k];

                g_x_0_xxyyyz_xxzzz[k] = -g_x_0_xxyyz_xxzzz[k] * cd_y[k] + g_x_0_xxyyz_xxyzzz[k];

                g_x_0_xxyyyz_xyyyy[k] = -g_x_0_xxyyz_xyyyy[k] * cd_y[k] + g_x_0_xxyyz_xyyyyy[k];

                g_x_0_xxyyyz_xyyyz[k] = -g_x_0_xxyyz_xyyyz[k] * cd_y[k] + g_x_0_xxyyz_xyyyyz[k];

                g_x_0_xxyyyz_xyyzz[k] = -g_x_0_xxyyz_xyyzz[k] * cd_y[k] + g_x_0_xxyyz_xyyyzz[k];

                g_x_0_xxyyyz_xyzzz[k] = -g_x_0_xxyyz_xyzzz[k] * cd_y[k] + g_x_0_xxyyz_xyyzzz[k];

                g_x_0_xxyyyz_xzzzz[k] = -g_x_0_xxyyz_xzzzz[k] * cd_y[k] + g_x_0_xxyyz_xyzzzz[k];

                g_x_0_xxyyyz_yyyyy[k] = -g_x_0_xxyyz_yyyyy[k] * cd_y[k] + g_x_0_xxyyz_yyyyyy[k];

                g_x_0_xxyyyz_yyyyz[k] = -g_x_0_xxyyz_yyyyz[k] * cd_y[k] + g_x_0_xxyyz_yyyyyz[k];

                g_x_0_xxyyyz_yyyzz[k] = -g_x_0_xxyyz_yyyzz[k] * cd_y[k] + g_x_0_xxyyz_yyyyzz[k];

                g_x_0_xxyyyz_yyzzz[k] = -g_x_0_xxyyz_yyzzz[k] * cd_y[k] + g_x_0_xxyyz_yyyzzz[k];

                g_x_0_xxyyyz_yzzzz[k] = -g_x_0_xxyyz_yzzzz[k] * cd_y[k] + g_x_0_xxyyz_yyzzzz[k];

                g_x_0_xxyyyz_zzzzz[k] = -g_x_0_xxyyz_zzzzz[k] * cd_y[k] + g_x_0_xxyyz_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xxyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xxyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xxyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xxyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xxyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xxyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xxyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xxyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xxyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xxyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xxyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xxyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xxyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xxyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xxyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xxyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xxyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xxyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xxyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xxyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyzz_xxxxx, g_x_0_xxyyzz_xxxxy, g_x_0_xxyyzz_xxxxz, g_x_0_xxyyzz_xxxyy, g_x_0_xxyyzz_xxxyz, g_x_0_xxyyzz_xxxzz, g_x_0_xxyyzz_xxyyy, g_x_0_xxyyzz_xxyyz, g_x_0_xxyyzz_xxyzz, g_x_0_xxyyzz_xxzzz, g_x_0_xxyyzz_xyyyy, g_x_0_xxyyzz_xyyyz, g_x_0_xxyyzz_xyyzz, g_x_0_xxyyzz_xyzzz, g_x_0_xxyyzz_xzzzz, g_x_0_xxyyzz_yyyyy, g_x_0_xxyyzz_yyyyz, g_x_0_xxyyzz_yyyzz, g_x_0_xxyyzz_yyzzz, g_x_0_xxyyzz_yzzzz, g_x_0_xxyyzz_zzzzz, g_x_0_xxyzz_xxxxx, g_x_0_xxyzz_xxxxxy, g_x_0_xxyzz_xxxxy, g_x_0_xxyzz_xxxxyy, g_x_0_xxyzz_xxxxyz, g_x_0_xxyzz_xxxxz, g_x_0_xxyzz_xxxyy, g_x_0_xxyzz_xxxyyy, g_x_0_xxyzz_xxxyyz, g_x_0_xxyzz_xxxyz, g_x_0_xxyzz_xxxyzz, g_x_0_xxyzz_xxxzz, g_x_0_xxyzz_xxyyy, g_x_0_xxyzz_xxyyyy, g_x_0_xxyzz_xxyyyz, g_x_0_xxyzz_xxyyz, g_x_0_xxyzz_xxyyzz, g_x_0_xxyzz_xxyzz, g_x_0_xxyzz_xxyzzz, g_x_0_xxyzz_xxzzz, g_x_0_xxyzz_xyyyy, g_x_0_xxyzz_xyyyyy, g_x_0_xxyzz_xyyyyz, g_x_0_xxyzz_xyyyz, g_x_0_xxyzz_xyyyzz, g_x_0_xxyzz_xyyzz, g_x_0_xxyzz_xyyzzz, g_x_0_xxyzz_xyzzz, g_x_0_xxyzz_xyzzzz, g_x_0_xxyzz_xzzzz, g_x_0_xxyzz_yyyyy, g_x_0_xxyzz_yyyyyy, g_x_0_xxyzz_yyyyyz, g_x_0_xxyzz_yyyyz, g_x_0_xxyzz_yyyyzz, g_x_0_xxyzz_yyyzz, g_x_0_xxyzz_yyyzzz, g_x_0_xxyzz_yyzzz, g_x_0_xxyzz_yyzzzz, g_x_0_xxyzz_yzzzz, g_x_0_xxyzz_yzzzzz, g_x_0_xxyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_xxxxx[k] = -g_x_0_xxyzz_xxxxx[k] * cd_y[k] + g_x_0_xxyzz_xxxxxy[k];

                g_x_0_xxyyzz_xxxxy[k] = -g_x_0_xxyzz_xxxxy[k] * cd_y[k] + g_x_0_xxyzz_xxxxyy[k];

                g_x_0_xxyyzz_xxxxz[k] = -g_x_0_xxyzz_xxxxz[k] * cd_y[k] + g_x_0_xxyzz_xxxxyz[k];

                g_x_0_xxyyzz_xxxyy[k] = -g_x_0_xxyzz_xxxyy[k] * cd_y[k] + g_x_0_xxyzz_xxxyyy[k];

                g_x_0_xxyyzz_xxxyz[k] = -g_x_0_xxyzz_xxxyz[k] * cd_y[k] + g_x_0_xxyzz_xxxyyz[k];

                g_x_0_xxyyzz_xxxzz[k] = -g_x_0_xxyzz_xxxzz[k] * cd_y[k] + g_x_0_xxyzz_xxxyzz[k];

                g_x_0_xxyyzz_xxyyy[k] = -g_x_0_xxyzz_xxyyy[k] * cd_y[k] + g_x_0_xxyzz_xxyyyy[k];

                g_x_0_xxyyzz_xxyyz[k] = -g_x_0_xxyzz_xxyyz[k] * cd_y[k] + g_x_0_xxyzz_xxyyyz[k];

                g_x_0_xxyyzz_xxyzz[k] = -g_x_0_xxyzz_xxyzz[k] * cd_y[k] + g_x_0_xxyzz_xxyyzz[k];

                g_x_0_xxyyzz_xxzzz[k] = -g_x_0_xxyzz_xxzzz[k] * cd_y[k] + g_x_0_xxyzz_xxyzzz[k];

                g_x_0_xxyyzz_xyyyy[k] = -g_x_0_xxyzz_xyyyy[k] * cd_y[k] + g_x_0_xxyzz_xyyyyy[k];

                g_x_0_xxyyzz_xyyyz[k] = -g_x_0_xxyzz_xyyyz[k] * cd_y[k] + g_x_0_xxyzz_xyyyyz[k];

                g_x_0_xxyyzz_xyyzz[k] = -g_x_0_xxyzz_xyyzz[k] * cd_y[k] + g_x_0_xxyzz_xyyyzz[k];

                g_x_0_xxyyzz_xyzzz[k] = -g_x_0_xxyzz_xyzzz[k] * cd_y[k] + g_x_0_xxyzz_xyyzzz[k];

                g_x_0_xxyyzz_xzzzz[k] = -g_x_0_xxyzz_xzzzz[k] * cd_y[k] + g_x_0_xxyzz_xyzzzz[k];

                g_x_0_xxyyzz_yyyyy[k] = -g_x_0_xxyzz_yyyyy[k] * cd_y[k] + g_x_0_xxyzz_yyyyyy[k];

                g_x_0_xxyyzz_yyyyz[k] = -g_x_0_xxyzz_yyyyz[k] * cd_y[k] + g_x_0_xxyzz_yyyyyz[k];

                g_x_0_xxyyzz_yyyzz[k] = -g_x_0_xxyzz_yyyzz[k] * cd_y[k] + g_x_0_xxyzz_yyyyzz[k];

                g_x_0_xxyyzz_yyzzz[k] = -g_x_0_xxyzz_yyzzz[k] * cd_y[k] + g_x_0_xxyzz_yyyzzz[k];

                g_x_0_xxyyzz_yzzzz[k] = -g_x_0_xxyzz_yzzzz[k] * cd_y[k] + g_x_0_xxyzz_yyzzzz[k];

                g_x_0_xxyyzz_zzzzz[k] = -g_x_0_xxyzz_zzzzz[k] * cd_y[k] + g_x_0_xxyzz_yzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xxyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xxyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xxyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xxyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xxyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xxyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_xxyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_xxyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_xxyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_xxyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_xxyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_xxyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_xxyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_xxyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_xxyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_xxyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_xxyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_xxyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_xxyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_xxyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzzz_xxxxx, g_x_0_xxyzzz_xxxxy, g_x_0_xxyzzz_xxxxz, g_x_0_xxyzzz_xxxyy, g_x_0_xxyzzz_xxxyz, g_x_0_xxyzzz_xxxzz, g_x_0_xxyzzz_xxyyy, g_x_0_xxyzzz_xxyyz, g_x_0_xxyzzz_xxyzz, g_x_0_xxyzzz_xxzzz, g_x_0_xxyzzz_xyyyy, g_x_0_xxyzzz_xyyyz, g_x_0_xxyzzz_xyyzz, g_x_0_xxyzzz_xyzzz, g_x_0_xxyzzz_xzzzz, g_x_0_xxyzzz_yyyyy, g_x_0_xxyzzz_yyyyz, g_x_0_xxyzzz_yyyzz, g_x_0_xxyzzz_yyzzz, g_x_0_xxyzzz_yzzzz, g_x_0_xxyzzz_zzzzz, g_x_0_xxzzz_xxxxx, g_x_0_xxzzz_xxxxxy, g_x_0_xxzzz_xxxxy, g_x_0_xxzzz_xxxxyy, g_x_0_xxzzz_xxxxyz, g_x_0_xxzzz_xxxxz, g_x_0_xxzzz_xxxyy, g_x_0_xxzzz_xxxyyy, g_x_0_xxzzz_xxxyyz, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxyzz, g_x_0_xxzzz_xxxzz, g_x_0_xxzzz_xxyyy, g_x_0_xxzzz_xxyyyy, g_x_0_xxzzz_xxyyyz, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyyzz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxyzzz, g_x_0_xxzzz_xxzzz, g_x_0_xxzzz_xyyyy, g_x_0_xxzzz_xyyyyy, g_x_0_xxzzz_xyyyyz, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyyzz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyyzzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xyzzzz, g_x_0_xxzzz_xzzzz, g_x_0_xxzzz_yyyyy, g_x_0_xxzzz_yyyyyy, g_x_0_xxzzz_yyyyyz, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyyzz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyyzzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yyzzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_yzzzzz, g_x_0_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_xxxxx[k] = -g_x_0_xxzzz_xxxxx[k] * cd_y[k] + g_x_0_xxzzz_xxxxxy[k];

                g_x_0_xxyzzz_xxxxy[k] = -g_x_0_xxzzz_xxxxy[k] * cd_y[k] + g_x_0_xxzzz_xxxxyy[k];

                g_x_0_xxyzzz_xxxxz[k] = -g_x_0_xxzzz_xxxxz[k] * cd_y[k] + g_x_0_xxzzz_xxxxyz[k];

                g_x_0_xxyzzz_xxxyy[k] = -g_x_0_xxzzz_xxxyy[k] * cd_y[k] + g_x_0_xxzzz_xxxyyy[k];

                g_x_0_xxyzzz_xxxyz[k] = -g_x_0_xxzzz_xxxyz[k] * cd_y[k] + g_x_0_xxzzz_xxxyyz[k];

                g_x_0_xxyzzz_xxxzz[k] = -g_x_0_xxzzz_xxxzz[k] * cd_y[k] + g_x_0_xxzzz_xxxyzz[k];

                g_x_0_xxyzzz_xxyyy[k] = -g_x_0_xxzzz_xxyyy[k] * cd_y[k] + g_x_0_xxzzz_xxyyyy[k];

                g_x_0_xxyzzz_xxyyz[k] = -g_x_0_xxzzz_xxyyz[k] * cd_y[k] + g_x_0_xxzzz_xxyyyz[k];

                g_x_0_xxyzzz_xxyzz[k] = -g_x_0_xxzzz_xxyzz[k] * cd_y[k] + g_x_0_xxzzz_xxyyzz[k];

                g_x_0_xxyzzz_xxzzz[k] = -g_x_0_xxzzz_xxzzz[k] * cd_y[k] + g_x_0_xxzzz_xxyzzz[k];

                g_x_0_xxyzzz_xyyyy[k] = -g_x_0_xxzzz_xyyyy[k] * cd_y[k] + g_x_0_xxzzz_xyyyyy[k];

                g_x_0_xxyzzz_xyyyz[k] = -g_x_0_xxzzz_xyyyz[k] * cd_y[k] + g_x_0_xxzzz_xyyyyz[k];

                g_x_0_xxyzzz_xyyzz[k] = -g_x_0_xxzzz_xyyzz[k] * cd_y[k] + g_x_0_xxzzz_xyyyzz[k];

                g_x_0_xxyzzz_xyzzz[k] = -g_x_0_xxzzz_xyzzz[k] * cd_y[k] + g_x_0_xxzzz_xyyzzz[k];

                g_x_0_xxyzzz_xzzzz[k] = -g_x_0_xxzzz_xzzzz[k] * cd_y[k] + g_x_0_xxzzz_xyzzzz[k];

                g_x_0_xxyzzz_yyyyy[k] = -g_x_0_xxzzz_yyyyy[k] * cd_y[k] + g_x_0_xxzzz_yyyyyy[k];

                g_x_0_xxyzzz_yyyyz[k] = -g_x_0_xxzzz_yyyyz[k] * cd_y[k] + g_x_0_xxzzz_yyyyyz[k];

                g_x_0_xxyzzz_yyyzz[k] = -g_x_0_xxzzz_yyyzz[k] * cd_y[k] + g_x_0_xxzzz_yyyyzz[k];

                g_x_0_xxyzzz_yyzzz[k] = -g_x_0_xxzzz_yyzzz[k] * cd_y[k] + g_x_0_xxzzz_yyyzzz[k];

                g_x_0_xxyzzz_yzzzz[k] = -g_x_0_xxzzz_yzzzz[k] * cd_y[k] + g_x_0_xxzzz_yyzzzz[k];

                g_x_0_xxyzzz_zzzzz[k] = -g_x_0_xxzzz_zzzzz[k] * cd_y[k] + g_x_0_xxzzz_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_xxzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_xxzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_xxzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_xxzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_xxzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_xxzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_xxzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_xxzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_xxzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_xxzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_xxzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_xxzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_xxzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_xxzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_xxzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_xxzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_xxzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_xxzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_xxzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_xxzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_x_0_xxzzz_xxxxx, g_x_0_xxzzz_xxxxxz, g_x_0_xxzzz_xxxxy, g_x_0_xxzzz_xxxxyz, g_x_0_xxzzz_xxxxz, g_x_0_xxzzz_xxxxzz, g_x_0_xxzzz_xxxyy, g_x_0_xxzzz_xxxyyz, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxyzz, g_x_0_xxzzz_xxxzz, g_x_0_xxzzz_xxxzzz, g_x_0_xxzzz_xxyyy, g_x_0_xxzzz_xxyyyz, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyyzz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxyzzz, g_x_0_xxzzz_xxzzz, g_x_0_xxzzz_xxzzzz, g_x_0_xxzzz_xyyyy, g_x_0_xxzzz_xyyyyz, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyyzz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyyzzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xyzzzz, g_x_0_xxzzz_xzzzz, g_x_0_xxzzz_xzzzzz, g_x_0_xxzzz_yyyyy, g_x_0_xxzzz_yyyyyz, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyyzz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyyzzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yyzzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_yzzzzz, g_x_0_xxzzz_zzzzz, g_x_0_xxzzz_zzzzzz, g_x_0_xxzzzz_xxxxx, g_x_0_xxzzzz_xxxxy, g_x_0_xxzzzz_xxxxz, g_x_0_xxzzzz_xxxyy, g_x_0_xxzzzz_xxxyz, g_x_0_xxzzzz_xxxzz, g_x_0_xxzzzz_xxyyy, g_x_0_xxzzzz_xxyyz, g_x_0_xxzzzz_xxyzz, g_x_0_xxzzzz_xxzzz, g_x_0_xxzzzz_xyyyy, g_x_0_xxzzzz_xyyyz, g_x_0_xxzzzz_xyyzz, g_x_0_xxzzzz_xyzzz, g_x_0_xxzzzz_xzzzz, g_x_0_xxzzzz_yyyyy, g_x_0_xxzzzz_yyyyz, g_x_0_xxzzzz_yyyzz, g_x_0_xxzzzz_yyzzz, g_x_0_xxzzzz_yzzzz, g_x_0_xxzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_xxxxx[k] = -g_x_0_xxzzz_xxxxx[k] * cd_z[k] + g_x_0_xxzzz_xxxxxz[k];

                g_x_0_xxzzzz_xxxxy[k] = -g_x_0_xxzzz_xxxxy[k] * cd_z[k] + g_x_0_xxzzz_xxxxyz[k];

                g_x_0_xxzzzz_xxxxz[k] = -g_x_0_xxzzz_xxxxz[k] * cd_z[k] + g_x_0_xxzzz_xxxxzz[k];

                g_x_0_xxzzzz_xxxyy[k] = -g_x_0_xxzzz_xxxyy[k] * cd_z[k] + g_x_0_xxzzz_xxxyyz[k];

                g_x_0_xxzzzz_xxxyz[k] = -g_x_0_xxzzz_xxxyz[k] * cd_z[k] + g_x_0_xxzzz_xxxyzz[k];

                g_x_0_xxzzzz_xxxzz[k] = -g_x_0_xxzzz_xxxzz[k] * cd_z[k] + g_x_0_xxzzz_xxxzzz[k];

                g_x_0_xxzzzz_xxyyy[k] = -g_x_0_xxzzz_xxyyy[k] * cd_z[k] + g_x_0_xxzzz_xxyyyz[k];

                g_x_0_xxzzzz_xxyyz[k] = -g_x_0_xxzzz_xxyyz[k] * cd_z[k] + g_x_0_xxzzz_xxyyzz[k];

                g_x_0_xxzzzz_xxyzz[k] = -g_x_0_xxzzz_xxyzz[k] * cd_z[k] + g_x_0_xxzzz_xxyzzz[k];

                g_x_0_xxzzzz_xxzzz[k] = -g_x_0_xxzzz_xxzzz[k] * cd_z[k] + g_x_0_xxzzz_xxzzzz[k];

                g_x_0_xxzzzz_xyyyy[k] = -g_x_0_xxzzz_xyyyy[k] * cd_z[k] + g_x_0_xxzzz_xyyyyz[k];

                g_x_0_xxzzzz_xyyyz[k] = -g_x_0_xxzzz_xyyyz[k] * cd_z[k] + g_x_0_xxzzz_xyyyzz[k];

                g_x_0_xxzzzz_xyyzz[k] = -g_x_0_xxzzz_xyyzz[k] * cd_z[k] + g_x_0_xxzzz_xyyzzz[k];

                g_x_0_xxzzzz_xyzzz[k] = -g_x_0_xxzzz_xyzzz[k] * cd_z[k] + g_x_0_xxzzz_xyzzzz[k];

                g_x_0_xxzzzz_xzzzz[k] = -g_x_0_xxzzz_xzzzz[k] * cd_z[k] + g_x_0_xxzzz_xzzzzz[k];

                g_x_0_xxzzzz_yyyyy[k] = -g_x_0_xxzzz_yyyyy[k] * cd_z[k] + g_x_0_xxzzz_yyyyyz[k];

                g_x_0_xxzzzz_yyyyz[k] = -g_x_0_xxzzz_yyyyz[k] * cd_z[k] + g_x_0_xxzzz_yyyyzz[k];

                g_x_0_xxzzzz_yyyzz[k] = -g_x_0_xxzzz_yyyzz[k] * cd_z[k] + g_x_0_xxzzz_yyyzzz[k];

                g_x_0_xxzzzz_yyzzz[k] = -g_x_0_xxzzz_yyzzz[k] * cd_z[k] + g_x_0_xxzzz_yyzzzz[k];

                g_x_0_xxzzzz_yzzzz[k] = -g_x_0_xxzzz_yzzzz[k] * cd_z[k] + g_x_0_xxzzz_yzzzzz[k];

                g_x_0_xxzzzz_zzzzz[k] = -g_x_0_xxzzz_zzzzz[k] * cd_z[k] + g_x_0_xxzzz_zzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_xyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_xyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_xyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_xyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_xyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_xyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_xyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_xyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_xyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_xyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_xyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_xyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_xyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_xyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_xyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_xyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_xyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_xyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_xyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_xyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyy_xxxxx, g_x_0_xyyyy_xxxxxy, g_x_0_xyyyy_xxxxy, g_x_0_xyyyy_xxxxyy, g_x_0_xyyyy_xxxxyz, g_x_0_xyyyy_xxxxz, g_x_0_xyyyy_xxxyy, g_x_0_xyyyy_xxxyyy, g_x_0_xyyyy_xxxyyz, g_x_0_xyyyy_xxxyz, g_x_0_xyyyy_xxxyzz, g_x_0_xyyyy_xxxzz, g_x_0_xyyyy_xxyyy, g_x_0_xyyyy_xxyyyy, g_x_0_xyyyy_xxyyyz, g_x_0_xyyyy_xxyyz, g_x_0_xyyyy_xxyyzz, g_x_0_xyyyy_xxyzz, g_x_0_xyyyy_xxyzzz, g_x_0_xyyyy_xxzzz, g_x_0_xyyyy_xyyyy, g_x_0_xyyyy_xyyyyy, g_x_0_xyyyy_xyyyyz, g_x_0_xyyyy_xyyyz, g_x_0_xyyyy_xyyyzz, g_x_0_xyyyy_xyyzz, g_x_0_xyyyy_xyyzzz, g_x_0_xyyyy_xyzzz, g_x_0_xyyyy_xyzzzz, g_x_0_xyyyy_xzzzz, g_x_0_xyyyy_yyyyy, g_x_0_xyyyy_yyyyyy, g_x_0_xyyyy_yyyyyz, g_x_0_xyyyy_yyyyz, g_x_0_xyyyy_yyyyzz, g_x_0_xyyyy_yyyzz, g_x_0_xyyyy_yyyzzz, g_x_0_xyyyy_yyzzz, g_x_0_xyyyy_yyzzzz, g_x_0_xyyyy_yzzzz, g_x_0_xyyyy_yzzzzz, g_x_0_xyyyy_zzzzz, g_x_0_xyyyyy_xxxxx, g_x_0_xyyyyy_xxxxy, g_x_0_xyyyyy_xxxxz, g_x_0_xyyyyy_xxxyy, g_x_0_xyyyyy_xxxyz, g_x_0_xyyyyy_xxxzz, g_x_0_xyyyyy_xxyyy, g_x_0_xyyyyy_xxyyz, g_x_0_xyyyyy_xxyzz, g_x_0_xyyyyy_xxzzz, g_x_0_xyyyyy_xyyyy, g_x_0_xyyyyy_xyyyz, g_x_0_xyyyyy_xyyzz, g_x_0_xyyyyy_xyzzz, g_x_0_xyyyyy_xzzzz, g_x_0_xyyyyy_yyyyy, g_x_0_xyyyyy_yyyyz, g_x_0_xyyyyy_yyyzz, g_x_0_xyyyyy_yyzzz, g_x_0_xyyyyy_yzzzz, g_x_0_xyyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_xxxxx[k] = -g_x_0_xyyyy_xxxxx[k] * cd_y[k] + g_x_0_xyyyy_xxxxxy[k];

                g_x_0_xyyyyy_xxxxy[k] = -g_x_0_xyyyy_xxxxy[k] * cd_y[k] + g_x_0_xyyyy_xxxxyy[k];

                g_x_0_xyyyyy_xxxxz[k] = -g_x_0_xyyyy_xxxxz[k] * cd_y[k] + g_x_0_xyyyy_xxxxyz[k];

                g_x_0_xyyyyy_xxxyy[k] = -g_x_0_xyyyy_xxxyy[k] * cd_y[k] + g_x_0_xyyyy_xxxyyy[k];

                g_x_0_xyyyyy_xxxyz[k] = -g_x_0_xyyyy_xxxyz[k] * cd_y[k] + g_x_0_xyyyy_xxxyyz[k];

                g_x_0_xyyyyy_xxxzz[k] = -g_x_0_xyyyy_xxxzz[k] * cd_y[k] + g_x_0_xyyyy_xxxyzz[k];

                g_x_0_xyyyyy_xxyyy[k] = -g_x_0_xyyyy_xxyyy[k] * cd_y[k] + g_x_0_xyyyy_xxyyyy[k];

                g_x_0_xyyyyy_xxyyz[k] = -g_x_0_xyyyy_xxyyz[k] * cd_y[k] + g_x_0_xyyyy_xxyyyz[k];

                g_x_0_xyyyyy_xxyzz[k] = -g_x_0_xyyyy_xxyzz[k] * cd_y[k] + g_x_0_xyyyy_xxyyzz[k];

                g_x_0_xyyyyy_xxzzz[k] = -g_x_0_xyyyy_xxzzz[k] * cd_y[k] + g_x_0_xyyyy_xxyzzz[k];

                g_x_0_xyyyyy_xyyyy[k] = -g_x_0_xyyyy_xyyyy[k] * cd_y[k] + g_x_0_xyyyy_xyyyyy[k];

                g_x_0_xyyyyy_xyyyz[k] = -g_x_0_xyyyy_xyyyz[k] * cd_y[k] + g_x_0_xyyyy_xyyyyz[k];

                g_x_0_xyyyyy_xyyzz[k] = -g_x_0_xyyyy_xyyzz[k] * cd_y[k] + g_x_0_xyyyy_xyyyzz[k];

                g_x_0_xyyyyy_xyzzz[k] = -g_x_0_xyyyy_xyzzz[k] * cd_y[k] + g_x_0_xyyyy_xyyzzz[k];

                g_x_0_xyyyyy_xzzzz[k] = -g_x_0_xyyyy_xzzzz[k] * cd_y[k] + g_x_0_xyyyy_xyzzzz[k];

                g_x_0_xyyyyy_yyyyy[k] = -g_x_0_xyyyy_yyyyy[k] * cd_y[k] + g_x_0_xyyyy_yyyyyy[k];

                g_x_0_xyyyyy_yyyyz[k] = -g_x_0_xyyyy_yyyyz[k] * cd_y[k] + g_x_0_xyyyy_yyyyyz[k];

                g_x_0_xyyyyy_yyyzz[k] = -g_x_0_xyyyy_yyyzz[k] * cd_y[k] + g_x_0_xyyyy_yyyyzz[k];

                g_x_0_xyyyyy_yyzzz[k] = -g_x_0_xyyyy_yyzzz[k] * cd_y[k] + g_x_0_xyyyy_yyyzzz[k];

                g_x_0_xyyyyy_yzzzz[k] = -g_x_0_xyyyy_yzzzz[k] * cd_y[k] + g_x_0_xyyyy_yyzzzz[k];

                g_x_0_xyyyyy_zzzzz[k] = -g_x_0_xyyyy_zzzzz[k] * cd_y[k] + g_x_0_xyyyy_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_xyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_xyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_xyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_xyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_xyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_xyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_xyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_xyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_xyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_xyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_xyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_xyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_xyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_xyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_xyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_xyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_xyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_xyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_xyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_xyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 356);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyyz_xxxxx, g_x_0_xyyyyz_xxxxy, g_x_0_xyyyyz_xxxxz, g_x_0_xyyyyz_xxxyy, g_x_0_xyyyyz_xxxyz, g_x_0_xyyyyz_xxxzz, g_x_0_xyyyyz_xxyyy, g_x_0_xyyyyz_xxyyz, g_x_0_xyyyyz_xxyzz, g_x_0_xyyyyz_xxzzz, g_x_0_xyyyyz_xyyyy, g_x_0_xyyyyz_xyyyz, g_x_0_xyyyyz_xyyzz, g_x_0_xyyyyz_xyzzz, g_x_0_xyyyyz_xzzzz, g_x_0_xyyyyz_yyyyy, g_x_0_xyyyyz_yyyyz, g_x_0_xyyyyz_yyyzz, g_x_0_xyyyyz_yyzzz, g_x_0_xyyyyz_yzzzz, g_x_0_xyyyyz_zzzzz, g_x_0_xyyyz_xxxxx, g_x_0_xyyyz_xxxxxy, g_x_0_xyyyz_xxxxy, g_x_0_xyyyz_xxxxyy, g_x_0_xyyyz_xxxxyz, g_x_0_xyyyz_xxxxz, g_x_0_xyyyz_xxxyy, g_x_0_xyyyz_xxxyyy, g_x_0_xyyyz_xxxyyz, g_x_0_xyyyz_xxxyz, g_x_0_xyyyz_xxxyzz, g_x_0_xyyyz_xxxzz, g_x_0_xyyyz_xxyyy, g_x_0_xyyyz_xxyyyy, g_x_0_xyyyz_xxyyyz, g_x_0_xyyyz_xxyyz, g_x_0_xyyyz_xxyyzz, g_x_0_xyyyz_xxyzz, g_x_0_xyyyz_xxyzzz, g_x_0_xyyyz_xxzzz, g_x_0_xyyyz_xyyyy, g_x_0_xyyyz_xyyyyy, g_x_0_xyyyz_xyyyyz, g_x_0_xyyyz_xyyyz, g_x_0_xyyyz_xyyyzz, g_x_0_xyyyz_xyyzz, g_x_0_xyyyz_xyyzzz, g_x_0_xyyyz_xyzzz, g_x_0_xyyyz_xyzzzz, g_x_0_xyyyz_xzzzz, g_x_0_xyyyz_yyyyy, g_x_0_xyyyz_yyyyyy, g_x_0_xyyyz_yyyyyz, g_x_0_xyyyz_yyyyz, g_x_0_xyyyz_yyyyzz, g_x_0_xyyyz_yyyzz, g_x_0_xyyyz_yyyzzz, g_x_0_xyyyz_yyzzz, g_x_0_xyyyz_yyzzzz, g_x_0_xyyyz_yzzzz, g_x_0_xyyyz_yzzzzz, g_x_0_xyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_xxxxx[k] = -g_x_0_xyyyz_xxxxx[k] * cd_y[k] + g_x_0_xyyyz_xxxxxy[k];

                g_x_0_xyyyyz_xxxxy[k] = -g_x_0_xyyyz_xxxxy[k] * cd_y[k] + g_x_0_xyyyz_xxxxyy[k];

                g_x_0_xyyyyz_xxxxz[k] = -g_x_0_xyyyz_xxxxz[k] * cd_y[k] + g_x_0_xyyyz_xxxxyz[k];

                g_x_0_xyyyyz_xxxyy[k] = -g_x_0_xyyyz_xxxyy[k] * cd_y[k] + g_x_0_xyyyz_xxxyyy[k];

                g_x_0_xyyyyz_xxxyz[k] = -g_x_0_xyyyz_xxxyz[k] * cd_y[k] + g_x_0_xyyyz_xxxyyz[k];

                g_x_0_xyyyyz_xxxzz[k] = -g_x_0_xyyyz_xxxzz[k] * cd_y[k] + g_x_0_xyyyz_xxxyzz[k];

                g_x_0_xyyyyz_xxyyy[k] = -g_x_0_xyyyz_xxyyy[k] * cd_y[k] + g_x_0_xyyyz_xxyyyy[k];

                g_x_0_xyyyyz_xxyyz[k] = -g_x_0_xyyyz_xxyyz[k] * cd_y[k] + g_x_0_xyyyz_xxyyyz[k];

                g_x_0_xyyyyz_xxyzz[k] = -g_x_0_xyyyz_xxyzz[k] * cd_y[k] + g_x_0_xyyyz_xxyyzz[k];

                g_x_0_xyyyyz_xxzzz[k] = -g_x_0_xyyyz_xxzzz[k] * cd_y[k] + g_x_0_xyyyz_xxyzzz[k];

                g_x_0_xyyyyz_xyyyy[k] = -g_x_0_xyyyz_xyyyy[k] * cd_y[k] + g_x_0_xyyyz_xyyyyy[k];

                g_x_0_xyyyyz_xyyyz[k] = -g_x_0_xyyyz_xyyyz[k] * cd_y[k] + g_x_0_xyyyz_xyyyyz[k];

                g_x_0_xyyyyz_xyyzz[k] = -g_x_0_xyyyz_xyyzz[k] * cd_y[k] + g_x_0_xyyyz_xyyyzz[k];

                g_x_0_xyyyyz_xyzzz[k] = -g_x_0_xyyyz_xyzzz[k] * cd_y[k] + g_x_0_xyyyz_xyyzzz[k];

                g_x_0_xyyyyz_xzzzz[k] = -g_x_0_xyyyz_xzzzz[k] * cd_y[k] + g_x_0_xyyyz_xyzzzz[k];

                g_x_0_xyyyyz_yyyyy[k] = -g_x_0_xyyyz_yyyyy[k] * cd_y[k] + g_x_0_xyyyz_yyyyyy[k];

                g_x_0_xyyyyz_yyyyz[k] = -g_x_0_xyyyz_yyyyz[k] * cd_y[k] + g_x_0_xyyyz_yyyyyz[k];

                g_x_0_xyyyyz_yyyzz[k] = -g_x_0_xyyyz_yyyzz[k] * cd_y[k] + g_x_0_xyyyz_yyyyzz[k];

                g_x_0_xyyyyz_yyzzz[k] = -g_x_0_xyyyz_yyzzz[k] * cd_y[k] + g_x_0_xyyyz_yyyzzz[k];

                g_x_0_xyyyyz_yzzzz[k] = -g_x_0_xyyyz_yzzzz[k] * cd_y[k] + g_x_0_xyyyz_yyzzzz[k];

                g_x_0_xyyyyz_zzzzz[k] = -g_x_0_xyyyz_zzzzz[k] * cd_y[k] + g_x_0_xyyyz_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_xyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_xyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_xyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_xyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_xyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_xyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 363);

            auto g_x_0_xyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_xyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_xyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_xyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_xyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_xyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_xyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_xyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_xyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_xyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_xyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_xyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_xyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_xyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 377);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyzz_xxxxx, g_x_0_xyyyzz_xxxxy, g_x_0_xyyyzz_xxxxz, g_x_0_xyyyzz_xxxyy, g_x_0_xyyyzz_xxxyz, g_x_0_xyyyzz_xxxzz, g_x_0_xyyyzz_xxyyy, g_x_0_xyyyzz_xxyyz, g_x_0_xyyyzz_xxyzz, g_x_0_xyyyzz_xxzzz, g_x_0_xyyyzz_xyyyy, g_x_0_xyyyzz_xyyyz, g_x_0_xyyyzz_xyyzz, g_x_0_xyyyzz_xyzzz, g_x_0_xyyyzz_xzzzz, g_x_0_xyyyzz_yyyyy, g_x_0_xyyyzz_yyyyz, g_x_0_xyyyzz_yyyzz, g_x_0_xyyyzz_yyzzz, g_x_0_xyyyzz_yzzzz, g_x_0_xyyyzz_zzzzz, g_x_0_xyyzz_xxxxx, g_x_0_xyyzz_xxxxxy, g_x_0_xyyzz_xxxxy, g_x_0_xyyzz_xxxxyy, g_x_0_xyyzz_xxxxyz, g_x_0_xyyzz_xxxxz, g_x_0_xyyzz_xxxyy, g_x_0_xyyzz_xxxyyy, g_x_0_xyyzz_xxxyyz, g_x_0_xyyzz_xxxyz, g_x_0_xyyzz_xxxyzz, g_x_0_xyyzz_xxxzz, g_x_0_xyyzz_xxyyy, g_x_0_xyyzz_xxyyyy, g_x_0_xyyzz_xxyyyz, g_x_0_xyyzz_xxyyz, g_x_0_xyyzz_xxyyzz, g_x_0_xyyzz_xxyzz, g_x_0_xyyzz_xxyzzz, g_x_0_xyyzz_xxzzz, g_x_0_xyyzz_xyyyy, g_x_0_xyyzz_xyyyyy, g_x_0_xyyzz_xyyyyz, g_x_0_xyyzz_xyyyz, g_x_0_xyyzz_xyyyzz, g_x_0_xyyzz_xyyzz, g_x_0_xyyzz_xyyzzz, g_x_0_xyyzz_xyzzz, g_x_0_xyyzz_xyzzzz, g_x_0_xyyzz_xzzzz, g_x_0_xyyzz_yyyyy, g_x_0_xyyzz_yyyyyy, g_x_0_xyyzz_yyyyyz, g_x_0_xyyzz_yyyyz, g_x_0_xyyzz_yyyyzz, g_x_0_xyyzz_yyyzz, g_x_0_xyyzz_yyyzzz, g_x_0_xyyzz_yyzzz, g_x_0_xyyzz_yyzzzz, g_x_0_xyyzz_yzzzz, g_x_0_xyyzz_yzzzzz, g_x_0_xyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_xxxxx[k] = -g_x_0_xyyzz_xxxxx[k] * cd_y[k] + g_x_0_xyyzz_xxxxxy[k];

                g_x_0_xyyyzz_xxxxy[k] = -g_x_0_xyyzz_xxxxy[k] * cd_y[k] + g_x_0_xyyzz_xxxxyy[k];

                g_x_0_xyyyzz_xxxxz[k] = -g_x_0_xyyzz_xxxxz[k] * cd_y[k] + g_x_0_xyyzz_xxxxyz[k];

                g_x_0_xyyyzz_xxxyy[k] = -g_x_0_xyyzz_xxxyy[k] * cd_y[k] + g_x_0_xyyzz_xxxyyy[k];

                g_x_0_xyyyzz_xxxyz[k] = -g_x_0_xyyzz_xxxyz[k] * cd_y[k] + g_x_0_xyyzz_xxxyyz[k];

                g_x_0_xyyyzz_xxxzz[k] = -g_x_0_xyyzz_xxxzz[k] * cd_y[k] + g_x_0_xyyzz_xxxyzz[k];

                g_x_0_xyyyzz_xxyyy[k] = -g_x_0_xyyzz_xxyyy[k] * cd_y[k] + g_x_0_xyyzz_xxyyyy[k];

                g_x_0_xyyyzz_xxyyz[k] = -g_x_0_xyyzz_xxyyz[k] * cd_y[k] + g_x_0_xyyzz_xxyyyz[k];

                g_x_0_xyyyzz_xxyzz[k] = -g_x_0_xyyzz_xxyzz[k] * cd_y[k] + g_x_0_xyyzz_xxyyzz[k];

                g_x_0_xyyyzz_xxzzz[k] = -g_x_0_xyyzz_xxzzz[k] * cd_y[k] + g_x_0_xyyzz_xxyzzz[k];

                g_x_0_xyyyzz_xyyyy[k] = -g_x_0_xyyzz_xyyyy[k] * cd_y[k] + g_x_0_xyyzz_xyyyyy[k];

                g_x_0_xyyyzz_xyyyz[k] = -g_x_0_xyyzz_xyyyz[k] * cd_y[k] + g_x_0_xyyzz_xyyyyz[k];

                g_x_0_xyyyzz_xyyzz[k] = -g_x_0_xyyzz_xyyzz[k] * cd_y[k] + g_x_0_xyyzz_xyyyzz[k];

                g_x_0_xyyyzz_xyzzz[k] = -g_x_0_xyyzz_xyzzz[k] * cd_y[k] + g_x_0_xyyzz_xyyzzz[k];

                g_x_0_xyyyzz_xzzzz[k] = -g_x_0_xyyzz_xzzzz[k] * cd_y[k] + g_x_0_xyyzz_xyzzzz[k];

                g_x_0_xyyyzz_yyyyy[k] = -g_x_0_xyyzz_yyyyy[k] * cd_y[k] + g_x_0_xyyzz_yyyyyy[k];

                g_x_0_xyyyzz_yyyyz[k] = -g_x_0_xyyzz_yyyyz[k] * cd_y[k] + g_x_0_xyyzz_yyyyyz[k];

                g_x_0_xyyyzz_yyyzz[k] = -g_x_0_xyyzz_yyyzz[k] * cd_y[k] + g_x_0_xyyzz_yyyyzz[k];

                g_x_0_xyyyzz_yyzzz[k] = -g_x_0_xyyzz_yyzzz[k] * cd_y[k] + g_x_0_xyyzz_yyyzzz[k];

                g_x_0_xyyyzz_yzzzz[k] = -g_x_0_xyyzz_yzzzz[k] * cd_y[k] + g_x_0_xyyzz_yyzzzz[k];

                g_x_0_xyyyzz_zzzzz[k] = -g_x_0_xyyzz_zzzzz[k] * cd_y[k] + g_x_0_xyyzz_yzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_xyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_xyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_xyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_xyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_xyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_xyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_xyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_xyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_xyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_xyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_xyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_xyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_xyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 391);

            auto g_x_0_xyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_xyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_xyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_xyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_xyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_xyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_xyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 398);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzzz_xxxxx, g_x_0_xyyzzz_xxxxy, g_x_0_xyyzzz_xxxxz, g_x_0_xyyzzz_xxxyy, g_x_0_xyyzzz_xxxyz, g_x_0_xyyzzz_xxxzz, g_x_0_xyyzzz_xxyyy, g_x_0_xyyzzz_xxyyz, g_x_0_xyyzzz_xxyzz, g_x_0_xyyzzz_xxzzz, g_x_0_xyyzzz_xyyyy, g_x_0_xyyzzz_xyyyz, g_x_0_xyyzzz_xyyzz, g_x_0_xyyzzz_xyzzz, g_x_0_xyyzzz_xzzzz, g_x_0_xyyzzz_yyyyy, g_x_0_xyyzzz_yyyyz, g_x_0_xyyzzz_yyyzz, g_x_0_xyyzzz_yyzzz, g_x_0_xyyzzz_yzzzz, g_x_0_xyyzzz_zzzzz, g_x_0_xyzzz_xxxxx, g_x_0_xyzzz_xxxxxy, g_x_0_xyzzz_xxxxy, g_x_0_xyzzz_xxxxyy, g_x_0_xyzzz_xxxxyz, g_x_0_xyzzz_xxxxz, g_x_0_xyzzz_xxxyy, g_x_0_xyzzz_xxxyyy, g_x_0_xyzzz_xxxyyz, g_x_0_xyzzz_xxxyz, g_x_0_xyzzz_xxxyzz, g_x_0_xyzzz_xxxzz, g_x_0_xyzzz_xxyyy, g_x_0_xyzzz_xxyyyy, g_x_0_xyzzz_xxyyyz, g_x_0_xyzzz_xxyyz, g_x_0_xyzzz_xxyyzz, g_x_0_xyzzz_xxyzz, g_x_0_xyzzz_xxyzzz, g_x_0_xyzzz_xxzzz, g_x_0_xyzzz_xyyyy, g_x_0_xyzzz_xyyyyy, g_x_0_xyzzz_xyyyyz, g_x_0_xyzzz_xyyyz, g_x_0_xyzzz_xyyyzz, g_x_0_xyzzz_xyyzz, g_x_0_xyzzz_xyyzzz, g_x_0_xyzzz_xyzzz, g_x_0_xyzzz_xyzzzz, g_x_0_xyzzz_xzzzz, g_x_0_xyzzz_yyyyy, g_x_0_xyzzz_yyyyyy, g_x_0_xyzzz_yyyyyz, g_x_0_xyzzz_yyyyz, g_x_0_xyzzz_yyyyzz, g_x_0_xyzzz_yyyzz, g_x_0_xyzzz_yyyzzz, g_x_0_xyzzz_yyzzz, g_x_0_xyzzz_yyzzzz, g_x_0_xyzzz_yzzzz, g_x_0_xyzzz_yzzzzz, g_x_0_xyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_xxxxx[k] = -g_x_0_xyzzz_xxxxx[k] * cd_y[k] + g_x_0_xyzzz_xxxxxy[k];

                g_x_0_xyyzzz_xxxxy[k] = -g_x_0_xyzzz_xxxxy[k] * cd_y[k] + g_x_0_xyzzz_xxxxyy[k];

                g_x_0_xyyzzz_xxxxz[k] = -g_x_0_xyzzz_xxxxz[k] * cd_y[k] + g_x_0_xyzzz_xxxxyz[k];

                g_x_0_xyyzzz_xxxyy[k] = -g_x_0_xyzzz_xxxyy[k] * cd_y[k] + g_x_0_xyzzz_xxxyyy[k];

                g_x_0_xyyzzz_xxxyz[k] = -g_x_0_xyzzz_xxxyz[k] * cd_y[k] + g_x_0_xyzzz_xxxyyz[k];

                g_x_0_xyyzzz_xxxzz[k] = -g_x_0_xyzzz_xxxzz[k] * cd_y[k] + g_x_0_xyzzz_xxxyzz[k];

                g_x_0_xyyzzz_xxyyy[k] = -g_x_0_xyzzz_xxyyy[k] * cd_y[k] + g_x_0_xyzzz_xxyyyy[k];

                g_x_0_xyyzzz_xxyyz[k] = -g_x_0_xyzzz_xxyyz[k] * cd_y[k] + g_x_0_xyzzz_xxyyyz[k];

                g_x_0_xyyzzz_xxyzz[k] = -g_x_0_xyzzz_xxyzz[k] * cd_y[k] + g_x_0_xyzzz_xxyyzz[k];

                g_x_0_xyyzzz_xxzzz[k] = -g_x_0_xyzzz_xxzzz[k] * cd_y[k] + g_x_0_xyzzz_xxyzzz[k];

                g_x_0_xyyzzz_xyyyy[k] = -g_x_0_xyzzz_xyyyy[k] * cd_y[k] + g_x_0_xyzzz_xyyyyy[k];

                g_x_0_xyyzzz_xyyyz[k] = -g_x_0_xyzzz_xyyyz[k] * cd_y[k] + g_x_0_xyzzz_xyyyyz[k];

                g_x_0_xyyzzz_xyyzz[k] = -g_x_0_xyzzz_xyyzz[k] * cd_y[k] + g_x_0_xyzzz_xyyyzz[k];

                g_x_0_xyyzzz_xyzzz[k] = -g_x_0_xyzzz_xyzzz[k] * cd_y[k] + g_x_0_xyzzz_xyyzzz[k];

                g_x_0_xyyzzz_xzzzz[k] = -g_x_0_xyzzz_xzzzz[k] * cd_y[k] + g_x_0_xyzzz_xyzzzz[k];

                g_x_0_xyyzzz_yyyyy[k] = -g_x_0_xyzzz_yyyyy[k] * cd_y[k] + g_x_0_xyzzz_yyyyyy[k];

                g_x_0_xyyzzz_yyyyz[k] = -g_x_0_xyzzz_yyyyz[k] * cd_y[k] + g_x_0_xyzzz_yyyyyz[k];

                g_x_0_xyyzzz_yyyzz[k] = -g_x_0_xyzzz_yyyzz[k] * cd_y[k] + g_x_0_xyzzz_yyyyzz[k];

                g_x_0_xyyzzz_yyzzz[k] = -g_x_0_xyzzz_yyzzz[k] * cd_y[k] + g_x_0_xyzzz_yyyzzz[k];

                g_x_0_xyyzzz_yzzzz[k] = -g_x_0_xyzzz_yzzzz[k] * cd_y[k] + g_x_0_xyzzz_yyzzzz[k];

                g_x_0_xyyzzz_zzzzz[k] = -g_x_0_xyzzz_zzzzz[k] * cd_y[k] + g_x_0_xyzzz_yzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_xyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_xyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_xyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_xyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_xyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_xyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_xyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_xyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_xyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_xyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_xyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_xyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_xyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_xyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_xyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_xyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_xyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_xyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_xyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_xyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzzz_xxxxx, g_x_0_xyzzzz_xxxxy, g_x_0_xyzzzz_xxxxz, g_x_0_xyzzzz_xxxyy, g_x_0_xyzzzz_xxxyz, g_x_0_xyzzzz_xxxzz, g_x_0_xyzzzz_xxyyy, g_x_0_xyzzzz_xxyyz, g_x_0_xyzzzz_xxyzz, g_x_0_xyzzzz_xxzzz, g_x_0_xyzzzz_xyyyy, g_x_0_xyzzzz_xyyyz, g_x_0_xyzzzz_xyyzz, g_x_0_xyzzzz_xyzzz, g_x_0_xyzzzz_xzzzz, g_x_0_xyzzzz_yyyyy, g_x_0_xyzzzz_yyyyz, g_x_0_xyzzzz_yyyzz, g_x_0_xyzzzz_yyzzz, g_x_0_xyzzzz_yzzzz, g_x_0_xyzzzz_zzzzz, g_x_0_xzzzz_xxxxx, g_x_0_xzzzz_xxxxxy, g_x_0_xzzzz_xxxxy, g_x_0_xzzzz_xxxxyy, g_x_0_xzzzz_xxxxyz, g_x_0_xzzzz_xxxxz, g_x_0_xzzzz_xxxyy, g_x_0_xzzzz_xxxyyy, g_x_0_xzzzz_xxxyyz, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxyzz, g_x_0_xzzzz_xxxzz, g_x_0_xzzzz_xxyyy, g_x_0_xzzzz_xxyyyy, g_x_0_xzzzz_xxyyyz, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyyzz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxyzzz, g_x_0_xzzzz_xxzzz, g_x_0_xzzzz_xyyyy, g_x_0_xzzzz_xyyyyy, g_x_0_xzzzz_xyyyyz, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyyzz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyyzzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xyzzzz, g_x_0_xzzzz_xzzzz, g_x_0_xzzzz_yyyyy, g_x_0_xzzzz_yyyyyy, g_x_0_xzzzz_yyyyyz, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyyzz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyyzzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yyzzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_yzzzzz, g_x_0_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_xxxxx[k] = -g_x_0_xzzzz_xxxxx[k] * cd_y[k] + g_x_0_xzzzz_xxxxxy[k];

                g_x_0_xyzzzz_xxxxy[k] = -g_x_0_xzzzz_xxxxy[k] * cd_y[k] + g_x_0_xzzzz_xxxxyy[k];

                g_x_0_xyzzzz_xxxxz[k] = -g_x_0_xzzzz_xxxxz[k] * cd_y[k] + g_x_0_xzzzz_xxxxyz[k];

                g_x_0_xyzzzz_xxxyy[k] = -g_x_0_xzzzz_xxxyy[k] * cd_y[k] + g_x_0_xzzzz_xxxyyy[k];

                g_x_0_xyzzzz_xxxyz[k] = -g_x_0_xzzzz_xxxyz[k] * cd_y[k] + g_x_0_xzzzz_xxxyyz[k];

                g_x_0_xyzzzz_xxxzz[k] = -g_x_0_xzzzz_xxxzz[k] * cd_y[k] + g_x_0_xzzzz_xxxyzz[k];

                g_x_0_xyzzzz_xxyyy[k] = -g_x_0_xzzzz_xxyyy[k] * cd_y[k] + g_x_0_xzzzz_xxyyyy[k];

                g_x_0_xyzzzz_xxyyz[k] = -g_x_0_xzzzz_xxyyz[k] * cd_y[k] + g_x_0_xzzzz_xxyyyz[k];

                g_x_0_xyzzzz_xxyzz[k] = -g_x_0_xzzzz_xxyzz[k] * cd_y[k] + g_x_0_xzzzz_xxyyzz[k];

                g_x_0_xyzzzz_xxzzz[k] = -g_x_0_xzzzz_xxzzz[k] * cd_y[k] + g_x_0_xzzzz_xxyzzz[k];

                g_x_0_xyzzzz_xyyyy[k] = -g_x_0_xzzzz_xyyyy[k] * cd_y[k] + g_x_0_xzzzz_xyyyyy[k];

                g_x_0_xyzzzz_xyyyz[k] = -g_x_0_xzzzz_xyyyz[k] * cd_y[k] + g_x_0_xzzzz_xyyyyz[k];

                g_x_0_xyzzzz_xyyzz[k] = -g_x_0_xzzzz_xyyzz[k] * cd_y[k] + g_x_0_xzzzz_xyyyzz[k];

                g_x_0_xyzzzz_xyzzz[k] = -g_x_0_xzzzz_xyzzz[k] * cd_y[k] + g_x_0_xzzzz_xyyzzz[k];

                g_x_0_xyzzzz_xzzzz[k] = -g_x_0_xzzzz_xzzzz[k] * cd_y[k] + g_x_0_xzzzz_xyzzzz[k];

                g_x_0_xyzzzz_yyyyy[k] = -g_x_0_xzzzz_yyyyy[k] * cd_y[k] + g_x_0_xzzzz_yyyyyy[k];

                g_x_0_xyzzzz_yyyyz[k] = -g_x_0_xzzzz_yyyyz[k] * cd_y[k] + g_x_0_xzzzz_yyyyyz[k];

                g_x_0_xyzzzz_yyyzz[k] = -g_x_0_xzzzz_yyyzz[k] * cd_y[k] + g_x_0_xzzzz_yyyyzz[k];

                g_x_0_xyzzzz_yyzzz[k] = -g_x_0_xzzzz_yyzzz[k] * cd_y[k] + g_x_0_xzzzz_yyyzzz[k];

                g_x_0_xyzzzz_yzzzz[k] = -g_x_0_xzzzz_yzzzz[k] * cd_y[k] + g_x_0_xzzzz_yyzzzz[k];

                g_x_0_xyzzzz_zzzzz[k] = -g_x_0_xzzzz_zzzzz[k] * cd_y[k] + g_x_0_xzzzz_yzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 420);

            auto g_x_0_xzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 421);

            auto g_x_0_xzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 422);

            auto g_x_0_xzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 423);

            auto g_x_0_xzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 424);

            auto g_x_0_xzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 425);

            auto g_x_0_xzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 426);

            auto g_x_0_xzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 427);

            auto g_x_0_xzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 428);

            auto g_x_0_xzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 429);

            auto g_x_0_xzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 430);

            auto g_x_0_xzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 431);

            auto g_x_0_xzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 432);

            auto g_x_0_xzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 433);

            auto g_x_0_xzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 434);

            auto g_x_0_xzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 435);

            auto g_x_0_xzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 436);

            auto g_x_0_xzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 437);

            auto g_x_0_xzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 438);

            auto g_x_0_xzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 439);

            auto g_x_0_xzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 440);

            #pragma omp simd aligned(cd_z, g_x_0_xzzzz_xxxxx, g_x_0_xzzzz_xxxxxz, g_x_0_xzzzz_xxxxy, g_x_0_xzzzz_xxxxyz, g_x_0_xzzzz_xxxxz, g_x_0_xzzzz_xxxxzz, g_x_0_xzzzz_xxxyy, g_x_0_xzzzz_xxxyyz, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxyzz, g_x_0_xzzzz_xxxzz, g_x_0_xzzzz_xxxzzz, g_x_0_xzzzz_xxyyy, g_x_0_xzzzz_xxyyyz, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyyzz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxyzzz, g_x_0_xzzzz_xxzzz, g_x_0_xzzzz_xxzzzz, g_x_0_xzzzz_xyyyy, g_x_0_xzzzz_xyyyyz, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyyzz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyyzzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xyzzzz, g_x_0_xzzzz_xzzzz, g_x_0_xzzzz_xzzzzz, g_x_0_xzzzz_yyyyy, g_x_0_xzzzz_yyyyyz, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyyzz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyyzzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yyzzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_yzzzzz, g_x_0_xzzzz_zzzzz, g_x_0_xzzzz_zzzzzz, g_x_0_xzzzzz_xxxxx, g_x_0_xzzzzz_xxxxy, g_x_0_xzzzzz_xxxxz, g_x_0_xzzzzz_xxxyy, g_x_0_xzzzzz_xxxyz, g_x_0_xzzzzz_xxxzz, g_x_0_xzzzzz_xxyyy, g_x_0_xzzzzz_xxyyz, g_x_0_xzzzzz_xxyzz, g_x_0_xzzzzz_xxzzz, g_x_0_xzzzzz_xyyyy, g_x_0_xzzzzz_xyyyz, g_x_0_xzzzzz_xyyzz, g_x_0_xzzzzz_xyzzz, g_x_0_xzzzzz_xzzzz, g_x_0_xzzzzz_yyyyy, g_x_0_xzzzzz_yyyyz, g_x_0_xzzzzz_yyyzz, g_x_0_xzzzzz_yyzzz, g_x_0_xzzzzz_yzzzz, g_x_0_xzzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_xxxxx[k] = -g_x_0_xzzzz_xxxxx[k] * cd_z[k] + g_x_0_xzzzz_xxxxxz[k];

                g_x_0_xzzzzz_xxxxy[k] = -g_x_0_xzzzz_xxxxy[k] * cd_z[k] + g_x_0_xzzzz_xxxxyz[k];

                g_x_0_xzzzzz_xxxxz[k] = -g_x_0_xzzzz_xxxxz[k] * cd_z[k] + g_x_0_xzzzz_xxxxzz[k];

                g_x_0_xzzzzz_xxxyy[k] = -g_x_0_xzzzz_xxxyy[k] * cd_z[k] + g_x_0_xzzzz_xxxyyz[k];

                g_x_0_xzzzzz_xxxyz[k] = -g_x_0_xzzzz_xxxyz[k] * cd_z[k] + g_x_0_xzzzz_xxxyzz[k];

                g_x_0_xzzzzz_xxxzz[k] = -g_x_0_xzzzz_xxxzz[k] * cd_z[k] + g_x_0_xzzzz_xxxzzz[k];

                g_x_0_xzzzzz_xxyyy[k] = -g_x_0_xzzzz_xxyyy[k] * cd_z[k] + g_x_0_xzzzz_xxyyyz[k];

                g_x_0_xzzzzz_xxyyz[k] = -g_x_0_xzzzz_xxyyz[k] * cd_z[k] + g_x_0_xzzzz_xxyyzz[k];

                g_x_0_xzzzzz_xxyzz[k] = -g_x_0_xzzzz_xxyzz[k] * cd_z[k] + g_x_0_xzzzz_xxyzzz[k];

                g_x_0_xzzzzz_xxzzz[k] = -g_x_0_xzzzz_xxzzz[k] * cd_z[k] + g_x_0_xzzzz_xxzzzz[k];

                g_x_0_xzzzzz_xyyyy[k] = -g_x_0_xzzzz_xyyyy[k] * cd_z[k] + g_x_0_xzzzz_xyyyyz[k];

                g_x_0_xzzzzz_xyyyz[k] = -g_x_0_xzzzz_xyyyz[k] * cd_z[k] + g_x_0_xzzzz_xyyyzz[k];

                g_x_0_xzzzzz_xyyzz[k] = -g_x_0_xzzzz_xyyzz[k] * cd_z[k] + g_x_0_xzzzz_xyyzzz[k];

                g_x_0_xzzzzz_xyzzz[k] = -g_x_0_xzzzz_xyzzz[k] * cd_z[k] + g_x_0_xzzzz_xyzzzz[k];

                g_x_0_xzzzzz_xzzzz[k] = -g_x_0_xzzzz_xzzzz[k] * cd_z[k] + g_x_0_xzzzz_xzzzzz[k];

                g_x_0_xzzzzz_yyyyy[k] = -g_x_0_xzzzz_yyyyy[k] * cd_z[k] + g_x_0_xzzzz_yyyyyz[k];

                g_x_0_xzzzzz_yyyyz[k] = -g_x_0_xzzzz_yyyyz[k] * cd_z[k] + g_x_0_xzzzz_yyyyzz[k];

                g_x_0_xzzzzz_yyyzz[k] = -g_x_0_xzzzz_yyyzz[k] * cd_z[k] + g_x_0_xzzzz_yyyzzz[k];

                g_x_0_xzzzzz_yyzzz[k] = -g_x_0_xzzzz_yyzzz[k] * cd_z[k] + g_x_0_xzzzz_yyzzzz[k];

                g_x_0_xzzzzz_yzzzz[k] = -g_x_0_xzzzz_yzzzz[k] * cd_z[k] + g_x_0_xzzzz_yzzzzz[k];

                g_x_0_xzzzzz_zzzzz[k] = -g_x_0_xzzzz_zzzzz[k] * cd_z[k] + g_x_0_xzzzz_zzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 441);

            auto g_x_0_yyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 442);

            auto g_x_0_yyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 443);

            auto g_x_0_yyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 444);

            auto g_x_0_yyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 445);

            auto g_x_0_yyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 446);

            auto g_x_0_yyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 447);

            auto g_x_0_yyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 448);

            auto g_x_0_yyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 449);

            auto g_x_0_yyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 450);

            auto g_x_0_yyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 451);

            auto g_x_0_yyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 452);

            auto g_x_0_yyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 453);

            auto g_x_0_yyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 454);

            auto g_x_0_yyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 455);

            auto g_x_0_yyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 456);

            auto g_x_0_yyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 457);

            auto g_x_0_yyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 458);

            auto g_x_0_yyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 459);

            auto g_x_0_yyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 460);

            auto g_x_0_yyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 461);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyy_xxxxx, g_x_0_yyyyy_xxxxxy, g_x_0_yyyyy_xxxxy, g_x_0_yyyyy_xxxxyy, g_x_0_yyyyy_xxxxyz, g_x_0_yyyyy_xxxxz, g_x_0_yyyyy_xxxyy, g_x_0_yyyyy_xxxyyy, g_x_0_yyyyy_xxxyyz, g_x_0_yyyyy_xxxyz, g_x_0_yyyyy_xxxyzz, g_x_0_yyyyy_xxxzz, g_x_0_yyyyy_xxyyy, g_x_0_yyyyy_xxyyyy, g_x_0_yyyyy_xxyyyz, g_x_0_yyyyy_xxyyz, g_x_0_yyyyy_xxyyzz, g_x_0_yyyyy_xxyzz, g_x_0_yyyyy_xxyzzz, g_x_0_yyyyy_xxzzz, g_x_0_yyyyy_xyyyy, g_x_0_yyyyy_xyyyyy, g_x_0_yyyyy_xyyyyz, g_x_0_yyyyy_xyyyz, g_x_0_yyyyy_xyyyzz, g_x_0_yyyyy_xyyzz, g_x_0_yyyyy_xyyzzz, g_x_0_yyyyy_xyzzz, g_x_0_yyyyy_xyzzzz, g_x_0_yyyyy_xzzzz, g_x_0_yyyyy_yyyyy, g_x_0_yyyyy_yyyyyy, g_x_0_yyyyy_yyyyyz, g_x_0_yyyyy_yyyyz, g_x_0_yyyyy_yyyyzz, g_x_0_yyyyy_yyyzz, g_x_0_yyyyy_yyyzzz, g_x_0_yyyyy_yyzzz, g_x_0_yyyyy_yyzzzz, g_x_0_yyyyy_yzzzz, g_x_0_yyyyy_yzzzzz, g_x_0_yyyyy_zzzzz, g_x_0_yyyyyy_xxxxx, g_x_0_yyyyyy_xxxxy, g_x_0_yyyyyy_xxxxz, g_x_0_yyyyyy_xxxyy, g_x_0_yyyyyy_xxxyz, g_x_0_yyyyyy_xxxzz, g_x_0_yyyyyy_xxyyy, g_x_0_yyyyyy_xxyyz, g_x_0_yyyyyy_xxyzz, g_x_0_yyyyyy_xxzzz, g_x_0_yyyyyy_xyyyy, g_x_0_yyyyyy_xyyyz, g_x_0_yyyyyy_xyyzz, g_x_0_yyyyyy_xyzzz, g_x_0_yyyyyy_xzzzz, g_x_0_yyyyyy_yyyyy, g_x_0_yyyyyy_yyyyz, g_x_0_yyyyyy_yyyzz, g_x_0_yyyyyy_yyzzz, g_x_0_yyyyyy_yzzzz, g_x_0_yyyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_xxxxx[k] = -g_x_0_yyyyy_xxxxx[k] * cd_y[k] + g_x_0_yyyyy_xxxxxy[k];

                g_x_0_yyyyyy_xxxxy[k] = -g_x_0_yyyyy_xxxxy[k] * cd_y[k] + g_x_0_yyyyy_xxxxyy[k];

                g_x_0_yyyyyy_xxxxz[k] = -g_x_0_yyyyy_xxxxz[k] * cd_y[k] + g_x_0_yyyyy_xxxxyz[k];

                g_x_0_yyyyyy_xxxyy[k] = -g_x_0_yyyyy_xxxyy[k] * cd_y[k] + g_x_0_yyyyy_xxxyyy[k];

                g_x_0_yyyyyy_xxxyz[k] = -g_x_0_yyyyy_xxxyz[k] * cd_y[k] + g_x_0_yyyyy_xxxyyz[k];

                g_x_0_yyyyyy_xxxzz[k] = -g_x_0_yyyyy_xxxzz[k] * cd_y[k] + g_x_0_yyyyy_xxxyzz[k];

                g_x_0_yyyyyy_xxyyy[k] = -g_x_0_yyyyy_xxyyy[k] * cd_y[k] + g_x_0_yyyyy_xxyyyy[k];

                g_x_0_yyyyyy_xxyyz[k] = -g_x_0_yyyyy_xxyyz[k] * cd_y[k] + g_x_0_yyyyy_xxyyyz[k];

                g_x_0_yyyyyy_xxyzz[k] = -g_x_0_yyyyy_xxyzz[k] * cd_y[k] + g_x_0_yyyyy_xxyyzz[k];

                g_x_0_yyyyyy_xxzzz[k] = -g_x_0_yyyyy_xxzzz[k] * cd_y[k] + g_x_0_yyyyy_xxyzzz[k];

                g_x_0_yyyyyy_xyyyy[k] = -g_x_0_yyyyy_xyyyy[k] * cd_y[k] + g_x_0_yyyyy_xyyyyy[k];

                g_x_0_yyyyyy_xyyyz[k] = -g_x_0_yyyyy_xyyyz[k] * cd_y[k] + g_x_0_yyyyy_xyyyyz[k];

                g_x_0_yyyyyy_xyyzz[k] = -g_x_0_yyyyy_xyyzz[k] * cd_y[k] + g_x_0_yyyyy_xyyyzz[k];

                g_x_0_yyyyyy_xyzzz[k] = -g_x_0_yyyyy_xyzzz[k] * cd_y[k] + g_x_0_yyyyy_xyyzzz[k];

                g_x_0_yyyyyy_xzzzz[k] = -g_x_0_yyyyy_xzzzz[k] * cd_y[k] + g_x_0_yyyyy_xyzzzz[k];

                g_x_0_yyyyyy_yyyyy[k] = -g_x_0_yyyyy_yyyyy[k] * cd_y[k] + g_x_0_yyyyy_yyyyyy[k];

                g_x_0_yyyyyy_yyyyz[k] = -g_x_0_yyyyy_yyyyz[k] * cd_y[k] + g_x_0_yyyyy_yyyyyz[k];

                g_x_0_yyyyyy_yyyzz[k] = -g_x_0_yyyyy_yyyzz[k] * cd_y[k] + g_x_0_yyyyy_yyyyzz[k];

                g_x_0_yyyyyy_yyzzz[k] = -g_x_0_yyyyy_yyzzz[k] * cd_y[k] + g_x_0_yyyyy_yyyzzz[k];

                g_x_0_yyyyyy_yzzzz[k] = -g_x_0_yyyyy_yzzzz[k] * cd_y[k] + g_x_0_yyyyy_yyzzzz[k];

                g_x_0_yyyyyy_zzzzz[k] = -g_x_0_yyyyy_zzzzz[k] * cd_y[k] + g_x_0_yyyyy_yzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 462);

            auto g_x_0_yyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 463);

            auto g_x_0_yyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 464);

            auto g_x_0_yyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 465);

            auto g_x_0_yyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 466);

            auto g_x_0_yyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 467);

            auto g_x_0_yyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 468);

            auto g_x_0_yyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 469);

            auto g_x_0_yyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 470);

            auto g_x_0_yyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 471);

            auto g_x_0_yyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 472);

            auto g_x_0_yyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 473);

            auto g_x_0_yyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 474);

            auto g_x_0_yyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 475);

            auto g_x_0_yyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 476);

            auto g_x_0_yyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 477);

            auto g_x_0_yyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 478);

            auto g_x_0_yyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 479);

            auto g_x_0_yyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 480);

            auto g_x_0_yyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 481);

            auto g_x_0_yyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 482);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyyz_xxxxx, g_x_0_yyyyyz_xxxxy, g_x_0_yyyyyz_xxxxz, g_x_0_yyyyyz_xxxyy, g_x_0_yyyyyz_xxxyz, g_x_0_yyyyyz_xxxzz, g_x_0_yyyyyz_xxyyy, g_x_0_yyyyyz_xxyyz, g_x_0_yyyyyz_xxyzz, g_x_0_yyyyyz_xxzzz, g_x_0_yyyyyz_xyyyy, g_x_0_yyyyyz_xyyyz, g_x_0_yyyyyz_xyyzz, g_x_0_yyyyyz_xyzzz, g_x_0_yyyyyz_xzzzz, g_x_0_yyyyyz_yyyyy, g_x_0_yyyyyz_yyyyz, g_x_0_yyyyyz_yyyzz, g_x_0_yyyyyz_yyzzz, g_x_0_yyyyyz_yzzzz, g_x_0_yyyyyz_zzzzz, g_x_0_yyyyz_xxxxx, g_x_0_yyyyz_xxxxxy, g_x_0_yyyyz_xxxxy, g_x_0_yyyyz_xxxxyy, g_x_0_yyyyz_xxxxyz, g_x_0_yyyyz_xxxxz, g_x_0_yyyyz_xxxyy, g_x_0_yyyyz_xxxyyy, g_x_0_yyyyz_xxxyyz, g_x_0_yyyyz_xxxyz, g_x_0_yyyyz_xxxyzz, g_x_0_yyyyz_xxxzz, g_x_0_yyyyz_xxyyy, g_x_0_yyyyz_xxyyyy, g_x_0_yyyyz_xxyyyz, g_x_0_yyyyz_xxyyz, g_x_0_yyyyz_xxyyzz, g_x_0_yyyyz_xxyzz, g_x_0_yyyyz_xxyzzz, g_x_0_yyyyz_xxzzz, g_x_0_yyyyz_xyyyy, g_x_0_yyyyz_xyyyyy, g_x_0_yyyyz_xyyyyz, g_x_0_yyyyz_xyyyz, g_x_0_yyyyz_xyyyzz, g_x_0_yyyyz_xyyzz, g_x_0_yyyyz_xyyzzz, g_x_0_yyyyz_xyzzz, g_x_0_yyyyz_xyzzzz, g_x_0_yyyyz_xzzzz, g_x_0_yyyyz_yyyyy, g_x_0_yyyyz_yyyyyy, g_x_0_yyyyz_yyyyyz, g_x_0_yyyyz_yyyyz, g_x_0_yyyyz_yyyyzz, g_x_0_yyyyz_yyyzz, g_x_0_yyyyz_yyyzzz, g_x_0_yyyyz_yyzzz, g_x_0_yyyyz_yyzzzz, g_x_0_yyyyz_yzzzz, g_x_0_yyyyz_yzzzzz, g_x_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_xxxxx[k] = -g_x_0_yyyyz_xxxxx[k] * cd_y[k] + g_x_0_yyyyz_xxxxxy[k];

                g_x_0_yyyyyz_xxxxy[k] = -g_x_0_yyyyz_xxxxy[k] * cd_y[k] + g_x_0_yyyyz_xxxxyy[k];

                g_x_0_yyyyyz_xxxxz[k] = -g_x_0_yyyyz_xxxxz[k] * cd_y[k] + g_x_0_yyyyz_xxxxyz[k];

                g_x_0_yyyyyz_xxxyy[k] = -g_x_0_yyyyz_xxxyy[k] * cd_y[k] + g_x_0_yyyyz_xxxyyy[k];

                g_x_0_yyyyyz_xxxyz[k] = -g_x_0_yyyyz_xxxyz[k] * cd_y[k] + g_x_0_yyyyz_xxxyyz[k];

                g_x_0_yyyyyz_xxxzz[k] = -g_x_0_yyyyz_xxxzz[k] * cd_y[k] + g_x_0_yyyyz_xxxyzz[k];

                g_x_0_yyyyyz_xxyyy[k] = -g_x_0_yyyyz_xxyyy[k] * cd_y[k] + g_x_0_yyyyz_xxyyyy[k];

                g_x_0_yyyyyz_xxyyz[k] = -g_x_0_yyyyz_xxyyz[k] * cd_y[k] + g_x_0_yyyyz_xxyyyz[k];

                g_x_0_yyyyyz_xxyzz[k] = -g_x_0_yyyyz_xxyzz[k] * cd_y[k] + g_x_0_yyyyz_xxyyzz[k];

                g_x_0_yyyyyz_xxzzz[k] = -g_x_0_yyyyz_xxzzz[k] * cd_y[k] + g_x_0_yyyyz_xxyzzz[k];

                g_x_0_yyyyyz_xyyyy[k] = -g_x_0_yyyyz_xyyyy[k] * cd_y[k] + g_x_0_yyyyz_xyyyyy[k];

                g_x_0_yyyyyz_xyyyz[k] = -g_x_0_yyyyz_xyyyz[k] * cd_y[k] + g_x_0_yyyyz_xyyyyz[k];

                g_x_0_yyyyyz_xyyzz[k] = -g_x_0_yyyyz_xyyzz[k] * cd_y[k] + g_x_0_yyyyz_xyyyzz[k];

                g_x_0_yyyyyz_xyzzz[k] = -g_x_0_yyyyz_xyzzz[k] * cd_y[k] + g_x_0_yyyyz_xyyzzz[k];

                g_x_0_yyyyyz_xzzzz[k] = -g_x_0_yyyyz_xzzzz[k] * cd_y[k] + g_x_0_yyyyz_xyzzzz[k];

                g_x_0_yyyyyz_yyyyy[k] = -g_x_0_yyyyz_yyyyy[k] * cd_y[k] + g_x_0_yyyyz_yyyyyy[k];

                g_x_0_yyyyyz_yyyyz[k] = -g_x_0_yyyyz_yyyyz[k] * cd_y[k] + g_x_0_yyyyz_yyyyyz[k];

                g_x_0_yyyyyz_yyyzz[k] = -g_x_0_yyyyz_yyyzz[k] * cd_y[k] + g_x_0_yyyyz_yyyyzz[k];

                g_x_0_yyyyyz_yyzzz[k] = -g_x_0_yyyyz_yyzzz[k] * cd_y[k] + g_x_0_yyyyz_yyyzzz[k];

                g_x_0_yyyyyz_yzzzz[k] = -g_x_0_yyyyz_yzzzz[k] * cd_y[k] + g_x_0_yyyyz_yyzzzz[k];

                g_x_0_yyyyyz_zzzzz[k] = -g_x_0_yyyyz_zzzzz[k] * cd_y[k] + g_x_0_yyyyz_yzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 483);

            auto g_x_0_yyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 484);

            auto g_x_0_yyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 485);

            auto g_x_0_yyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 486);

            auto g_x_0_yyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 487);

            auto g_x_0_yyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 488);

            auto g_x_0_yyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 489);

            auto g_x_0_yyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 490);

            auto g_x_0_yyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 491);

            auto g_x_0_yyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 492);

            auto g_x_0_yyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 493);

            auto g_x_0_yyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 494);

            auto g_x_0_yyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 495);

            auto g_x_0_yyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 496);

            auto g_x_0_yyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 497);

            auto g_x_0_yyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 498);

            auto g_x_0_yyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 499);

            auto g_x_0_yyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 500);

            auto g_x_0_yyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 501);

            auto g_x_0_yyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 502);

            auto g_x_0_yyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 503);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyzz_xxxxx, g_x_0_yyyyzz_xxxxy, g_x_0_yyyyzz_xxxxz, g_x_0_yyyyzz_xxxyy, g_x_0_yyyyzz_xxxyz, g_x_0_yyyyzz_xxxzz, g_x_0_yyyyzz_xxyyy, g_x_0_yyyyzz_xxyyz, g_x_0_yyyyzz_xxyzz, g_x_0_yyyyzz_xxzzz, g_x_0_yyyyzz_xyyyy, g_x_0_yyyyzz_xyyyz, g_x_0_yyyyzz_xyyzz, g_x_0_yyyyzz_xyzzz, g_x_0_yyyyzz_xzzzz, g_x_0_yyyyzz_yyyyy, g_x_0_yyyyzz_yyyyz, g_x_0_yyyyzz_yyyzz, g_x_0_yyyyzz_yyzzz, g_x_0_yyyyzz_yzzzz, g_x_0_yyyyzz_zzzzz, g_x_0_yyyzz_xxxxx, g_x_0_yyyzz_xxxxxy, g_x_0_yyyzz_xxxxy, g_x_0_yyyzz_xxxxyy, g_x_0_yyyzz_xxxxyz, g_x_0_yyyzz_xxxxz, g_x_0_yyyzz_xxxyy, g_x_0_yyyzz_xxxyyy, g_x_0_yyyzz_xxxyyz, g_x_0_yyyzz_xxxyz, g_x_0_yyyzz_xxxyzz, g_x_0_yyyzz_xxxzz, g_x_0_yyyzz_xxyyy, g_x_0_yyyzz_xxyyyy, g_x_0_yyyzz_xxyyyz, g_x_0_yyyzz_xxyyz, g_x_0_yyyzz_xxyyzz, g_x_0_yyyzz_xxyzz, g_x_0_yyyzz_xxyzzz, g_x_0_yyyzz_xxzzz, g_x_0_yyyzz_xyyyy, g_x_0_yyyzz_xyyyyy, g_x_0_yyyzz_xyyyyz, g_x_0_yyyzz_xyyyz, g_x_0_yyyzz_xyyyzz, g_x_0_yyyzz_xyyzz, g_x_0_yyyzz_xyyzzz, g_x_0_yyyzz_xyzzz, g_x_0_yyyzz_xyzzzz, g_x_0_yyyzz_xzzzz, g_x_0_yyyzz_yyyyy, g_x_0_yyyzz_yyyyyy, g_x_0_yyyzz_yyyyyz, g_x_0_yyyzz_yyyyz, g_x_0_yyyzz_yyyyzz, g_x_0_yyyzz_yyyzz, g_x_0_yyyzz_yyyzzz, g_x_0_yyyzz_yyzzz, g_x_0_yyyzz_yyzzzz, g_x_0_yyyzz_yzzzz, g_x_0_yyyzz_yzzzzz, g_x_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_xxxxx[k] = -g_x_0_yyyzz_xxxxx[k] * cd_y[k] + g_x_0_yyyzz_xxxxxy[k];

                g_x_0_yyyyzz_xxxxy[k] = -g_x_0_yyyzz_xxxxy[k] * cd_y[k] + g_x_0_yyyzz_xxxxyy[k];

                g_x_0_yyyyzz_xxxxz[k] = -g_x_0_yyyzz_xxxxz[k] * cd_y[k] + g_x_0_yyyzz_xxxxyz[k];

                g_x_0_yyyyzz_xxxyy[k] = -g_x_0_yyyzz_xxxyy[k] * cd_y[k] + g_x_0_yyyzz_xxxyyy[k];

                g_x_0_yyyyzz_xxxyz[k] = -g_x_0_yyyzz_xxxyz[k] * cd_y[k] + g_x_0_yyyzz_xxxyyz[k];

                g_x_0_yyyyzz_xxxzz[k] = -g_x_0_yyyzz_xxxzz[k] * cd_y[k] + g_x_0_yyyzz_xxxyzz[k];

                g_x_0_yyyyzz_xxyyy[k] = -g_x_0_yyyzz_xxyyy[k] * cd_y[k] + g_x_0_yyyzz_xxyyyy[k];

                g_x_0_yyyyzz_xxyyz[k] = -g_x_0_yyyzz_xxyyz[k] * cd_y[k] + g_x_0_yyyzz_xxyyyz[k];

                g_x_0_yyyyzz_xxyzz[k] = -g_x_0_yyyzz_xxyzz[k] * cd_y[k] + g_x_0_yyyzz_xxyyzz[k];

                g_x_0_yyyyzz_xxzzz[k] = -g_x_0_yyyzz_xxzzz[k] * cd_y[k] + g_x_0_yyyzz_xxyzzz[k];

                g_x_0_yyyyzz_xyyyy[k] = -g_x_0_yyyzz_xyyyy[k] * cd_y[k] + g_x_0_yyyzz_xyyyyy[k];

                g_x_0_yyyyzz_xyyyz[k] = -g_x_0_yyyzz_xyyyz[k] * cd_y[k] + g_x_0_yyyzz_xyyyyz[k];

                g_x_0_yyyyzz_xyyzz[k] = -g_x_0_yyyzz_xyyzz[k] * cd_y[k] + g_x_0_yyyzz_xyyyzz[k];

                g_x_0_yyyyzz_xyzzz[k] = -g_x_0_yyyzz_xyzzz[k] * cd_y[k] + g_x_0_yyyzz_xyyzzz[k];

                g_x_0_yyyyzz_xzzzz[k] = -g_x_0_yyyzz_xzzzz[k] * cd_y[k] + g_x_0_yyyzz_xyzzzz[k];

                g_x_0_yyyyzz_yyyyy[k] = -g_x_0_yyyzz_yyyyy[k] * cd_y[k] + g_x_0_yyyzz_yyyyyy[k];

                g_x_0_yyyyzz_yyyyz[k] = -g_x_0_yyyzz_yyyyz[k] * cd_y[k] + g_x_0_yyyzz_yyyyyz[k];

                g_x_0_yyyyzz_yyyzz[k] = -g_x_0_yyyzz_yyyzz[k] * cd_y[k] + g_x_0_yyyzz_yyyyzz[k];

                g_x_0_yyyyzz_yyzzz[k] = -g_x_0_yyyzz_yyzzz[k] * cd_y[k] + g_x_0_yyyzz_yyyzzz[k];

                g_x_0_yyyyzz_yzzzz[k] = -g_x_0_yyyzz_yzzzz[k] * cd_y[k] + g_x_0_yyyzz_yyzzzz[k];

                g_x_0_yyyyzz_zzzzz[k] = -g_x_0_yyyzz_zzzzz[k] * cd_y[k] + g_x_0_yyyzz_yzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 504);

            auto g_x_0_yyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 505);

            auto g_x_0_yyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 506);

            auto g_x_0_yyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 507);

            auto g_x_0_yyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 508);

            auto g_x_0_yyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 509);

            auto g_x_0_yyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 510);

            auto g_x_0_yyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 511);

            auto g_x_0_yyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 512);

            auto g_x_0_yyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 513);

            auto g_x_0_yyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 514);

            auto g_x_0_yyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 515);

            auto g_x_0_yyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 516);

            auto g_x_0_yyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 517);

            auto g_x_0_yyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 518);

            auto g_x_0_yyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 519);

            auto g_x_0_yyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 520);

            auto g_x_0_yyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 521);

            auto g_x_0_yyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 522);

            auto g_x_0_yyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 523);

            auto g_x_0_yyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 524);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzzz_xxxxx, g_x_0_yyyzzz_xxxxy, g_x_0_yyyzzz_xxxxz, g_x_0_yyyzzz_xxxyy, g_x_0_yyyzzz_xxxyz, g_x_0_yyyzzz_xxxzz, g_x_0_yyyzzz_xxyyy, g_x_0_yyyzzz_xxyyz, g_x_0_yyyzzz_xxyzz, g_x_0_yyyzzz_xxzzz, g_x_0_yyyzzz_xyyyy, g_x_0_yyyzzz_xyyyz, g_x_0_yyyzzz_xyyzz, g_x_0_yyyzzz_xyzzz, g_x_0_yyyzzz_xzzzz, g_x_0_yyyzzz_yyyyy, g_x_0_yyyzzz_yyyyz, g_x_0_yyyzzz_yyyzz, g_x_0_yyyzzz_yyzzz, g_x_0_yyyzzz_yzzzz, g_x_0_yyyzzz_zzzzz, g_x_0_yyzzz_xxxxx, g_x_0_yyzzz_xxxxxy, g_x_0_yyzzz_xxxxy, g_x_0_yyzzz_xxxxyy, g_x_0_yyzzz_xxxxyz, g_x_0_yyzzz_xxxxz, g_x_0_yyzzz_xxxyy, g_x_0_yyzzz_xxxyyy, g_x_0_yyzzz_xxxyyz, g_x_0_yyzzz_xxxyz, g_x_0_yyzzz_xxxyzz, g_x_0_yyzzz_xxxzz, g_x_0_yyzzz_xxyyy, g_x_0_yyzzz_xxyyyy, g_x_0_yyzzz_xxyyyz, g_x_0_yyzzz_xxyyz, g_x_0_yyzzz_xxyyzz, g_x_0_yyzzz_xxyzz, g_x_0_yyzzz_xxyzzz, g_x_0_yyzzz_xxzzz, g_x_0_yyzzz_xyyyy, g_x_0_yyzzz_xyyyyy, g_x_0_yyzzz_xyyyyz, g_x_0_yyzzz_xyyyz, g_x_0_yyzzz_xyyyzz, g_x_0_yyzzz_xyyzz, g_x_0_yyzzz_xyyzzz, g_x_0_yyzzz_xyzzz, g_x_0_yyzzz_xyzzzz, g_x_0_yyzzz_xzzzz, g_x_0_yyzzz_yyyyy, g_x_0_yyzzz_yyyyyy, g_x_0_yyzzz_yyyyyz, g_x_0_yyzzz_yyyyz, g_x_0_yyzzz_yyyyzz, g_x_0_yyzzz_yyyzz, g_x_0_yyzzz_yyyzzz, g_x_0_yyzzz_yyzzz, g_x_0_yyzzz_yyzzzz, g_x_0_yyzzz_yzzzz, g_x_0_yyzzz_yzzzzz, g_x_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_xxxxx[k] = -g_x_0_yyzzz_xxxxx[k] * cd_y[k] + g_x_0_yyzzz_xxxxxy[k];

                g_x_0_yyyzzz_xxxxy[k] = -g_x_0_yyzzz_xxxxy[k] * cd_y[k] + g_x_0_yyzzz_xxxxyy[k];

                g_x_0_yyyzzz_xxxxz[k] = -g_x_0_yyzzz_xxxxz[k] * cd_y[k] + g_x_0_yyzzz_xxxxyz[k];

                g_x_0_yyyzzz_xxxyy[k] = -g_x_0_yyzzz_xxxyy[k] * cd_y[k] + g_x_0_yyzzz_xxxyyy[k];

                g_x_0_yyyzzz_xxxyz[k] = -g_x_0_yyzzz_xxxyz[k] * cd_y[k] + g_x_0_yyzzz_xxxyyz[k];

                g_x_0_yyyzzz_xxxzz[k] = -g_x_0_yyzzz_xxxzz[k] * cd_y[k] + g_x_0_yyzzz_xxxyzz[k];

                g_x_0_yyyzzz_xxyyy[k] = -g_x_0_yyzzz_xxyyy[k] * cd_y[k] + g_x_0_yyzzz_xxyyyy[k];

                g_x_0_yyyzzz_xxyyz[k] = -g_x_0_yyzzz_xxyyz[k] * cd_y[k] + g_x_0_yyzzz_xxyyyz[k];

                g_x_0_yyyzzz_xxyzz[k] = -g_x_0_yyzzz_xxyzz[k] * cd_y[k] + g_x_0_yyzzz_xxyyzz[k];

                g_x_0_yyyzzz_xxzzz[k] = -g_x_0_yyzzz_xxzzz[k] * cd_y[k] + g_x_0_yyzzz_xxyzzz[k];

                g_x_0_yyyzzz_xyyyy[k] = -g_x_0_yyzzz_xyyyy[k] * cd_y[k] + g_x_0_yyzzz_xyyyyy[k];

                g_x_0_yyyzzz_xyyyz[k] = -g_x_0_yyzzz_xyyyz[k] * cd_y[k] + g_x_0_yyzzz_xyyyyz[k];

                g_x_0_yyyzzz_xyyzz[k] = -g_x_0_yyzzz_xyyzz[k] * cd_y[k] + g_x_0_yyzzz_xyyyzz[k];

                g_x_0_yyyzzz_xyzzz[k] = -g_x_0_yyzzz_xyzzz[k] * cd_y[k] + g_x_0_yyzzz_xyyzzz[k];

                g_x_0_yyyzzz_xzzzz[k] = -g_x_0_yyzzz_xzzzz[k] * cd_y[k] + g_x_0_yyzzz_xyzzzz[k];

                g_x_0_yyyzzz_yyyyy[k] = -g_x_0_yyzzz_yyyyy[k] * cd_y[k] + g_x_0_yyzzz_yyyyyy[k];

                g_x_0_yyyzzz_yyyyz[k] = -g_x_0_yyzzz_yyyyz[k] * cd_y[k] + g_x_0_yyzzz_yyyyyz[k];

                g_x_0_yyyzzz_yyyzz[k] = -g_x_0_yyzzz_yyyzz[k] * cd_y[k] + g_x_0_yyzzz_yyyyzz[k];

                g_x_0_yyyzzz_yyzzz[k] = -g_x_0_yyzzz_yyzzz[k] * cd_y[k] + g_x_0_yyzzz_yyyzzz[k];

                g_x_0_yyyzzz_yzzzz[k] = -g_x_0_yyzzz_yzzzz[k] * cd_y[k] + g_x_0_yyzzz_yyzzzz[k];

                g_x_0_yyyzzz_zzzzz[k] = -g_x_0_yyzzz_zzzzz[k] * cd_y[k] + g_x_0_yyzzz_yzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 525);

            auto g_x_0_yyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 526);

            auto g_x_0_yyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 527);

            auto g_x_0_yyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 528);

            auto g_x_0_yyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 529);

            auto g_x_0_yyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 530);

            auto g_x_0_yyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 531);

            auto g_x_0_yyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 532);

            auto g_x_0_yyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 533);

            auto g_x_0_yyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 534);

            auto g_x_0_yyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 535);

            auto g_x_0_yyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 536);

            auto g_x_0_yyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 537);

            auto g_x_0_yyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 538);

            auto g_x_0_yyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 539);

            auto g_x_0_yyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 540);

            auto g_x_0_yyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 541);

            auto g_x_0_yyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 542);

            auto g_x_0_yyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 543);

            auto g_x_0_yyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 544);

            auto g_x_0_yyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 545);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzzz_xxxxx, g_x_0_yyzzzz_xxxxy, g_x_0_yyzzzz_xxxxz, g_x_0_yyzzzz_xxxyy, g_x_0_yyzzzz_xxxyz, g_x_0_yyzzzz_xxxzz, g_x_0_yyzzzz_xxyyy, g_x_0_yyzzzz_xxyyz, g_x_0_yyzzzz_xxyzz, g_x_0_yyzzzz_xxzzz, g_x_0_yyzzzz_xyyyy, g_x_0_yyzzzz_xyyyz, g_x_0_yyzzzz_xyyzz, g_x_0_yyzzzz_xyzzz, g_x_0_yyzzzz_xzzzz, g_x_0_yyzzzz_yyyyy, g_x_0_yyzzzz_yyyyz, g_x_0_yyzzzz_yyyzz, g_x_0_yyzzzz_yyzzz, g_x_0_yyzzzz_yzzzz, g_x_0_yyzzzz_zzzzz, g_x_0_yzzzz_xxxxx, g_x_0_yzzzz_xxxxxy, g_x_0_yzzzz_xxxxy, g_x_0_yzzzz_xxxxyy, g_x_0_yzzzz_xxxxyz, g_x_0_yzzzz_xxxxz, g_x_0_yzzzz_xxxyy, g_x_0_yzzzz_xxxyyy, g_x_0_yzzzz_xxxyyz, g_x_0_yzzzz_xxxyz, g_x_0_yzzzz_xxxyzz, g_x_0_yzzzz_xxxzz, g_x_0_yzzzz_xxyyy, g_x_0_yzzzz_xxyyyy, g_x_0_yzzzz_xxyyyz, g_x_0_yzzzz_xxyyz, g_x_0_yzzzz_xxyyzz, g_x_0_yzzzz_xxyzz, g_x_0_yzzzz_xxyzzz, g_x_0_yzzzz_xxzzz, g_x_0_yzzzz_xyyyy, g_x_0_yzzzz_xyyyyy, g_x_0_yzzzz_xyyyyz, g_x_0_yzzzz_xyyyz, g_x_0_yzzzz_xyyyzz, g_x_0_yzzzz_xyyzz, g_x_0_yzzzz_xyyzzz, g_x_0_yzzzz_xyzzz, g_x_0_yzzzz_xyzzzz, g_x_0_yzzzz_xzzzz, g_x_0_yzzzz_yyyyy, g_x_0_yzzzz_yyyyyy, g_x_0_yzzzz_yyyyyz, g_x_0_yzzzz_yyyyz, g_x_0_yzzzz_yyyyzz, g_x_0_yzzzz_yyyzz, g_x_0_yzzzz_yyyzzz, g_x_0_yzzzz_yyzzz, g_x_0_yzzzz_yyzzzz, g_x_0_yzzzz_yzzzz, g_x_0_yzzzz_yzzzzz, g_x_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_xxxxx[k] = -g_x_0_yzzzz_xxxxx[k] * cd_y[k] + g_x_0_yzzzz_xxxxxy[k];

                g_x_0_yyzzzz_xxxxy[k] = -g_x_0_yzzzz_xxxxy[k] * cd_y[k] + g_x_0_yzzzz_xxxxyy[k];

                g_x_0_yyzzzz_xxxxz[k] = -g_x_0_yzzzz_xxxxz[k] * cd_y[k] + g_x_0_yzzzz_xxxxyz[k];

                g_x_0_yyzzzz_xxxyy[k] = -g_x_0_yzzzz_xxxyy[k] * cd_y[k] + g_x_0_yzzzz_xxxyyy[k];

                g_x_0_yyzzzz_xxxyz[k] = -g_x_0_yzzzz_xxxyz[k] * cd_y[k] + g_x_0_yzzzz_xxxyyz[k];

                g_x_0_yyzzzz_xxxzz[k] = -g_x_0_yzzzz_xxxzz[k] * cd_y[k] + g_x_0_yzzzz_xxxyzz[k];

                g_x_0_yyzzzz_xxyyy[k] = -g_x_0_yzzzz_xxyyy[k] * cd_y[k] + g_x_0_yzzzz_xxyyyy[k];

                g_x_0_yyzzzz_xxyyz[k] = -g_x_0_yzzzz_xxyyz[k] * cd_y[k] + g_x_0_yzzzz_xxyyyz[k];

                g_x_0_yyzzzz_xxyzz[k] = -g_x_0_yzzzz_xxyzz[k] * cd_y[k] + g_x_0_yzzzz_xxyyzz[k];

                g_x_0_yyzzzz_xxzzz[k] = -g_x_0_yzzzz_xxzzz[k] * cd_y[k] + g_x_0_yzzzz_xxyzzz[k];

                g_x_0_yyzzzz_xyyyy[k] = -g_x_0_yzzzz_xyyyy[k] * cd_y[k] + g_x_0_yzzzz_xyyyyy[k];

                g_x_0_yyzzzz_xyyyz[k] = -g_x_0_yzzzz_xyyyz[k] * cd_y[k] + g_x_0_yzzzz_xyyyyz[k];

                g_x_0_yyzzzz_xyyzz[k] = -g_x_0_yzzzz_xyyzz[k] * cd_y[k] + g_x_0_yzzzz_xyyyzz[k];

                g_x_0_yyzzzz_xyzzz[k] = -g_x_0_yzzzz_xyzzz[k] * cd_y[k] + g_x_0_yzzzz_xyyzzz[k];

                g_x_0_yyzzzz_xzzzz[k] = -g_x_0_yzzzz_xzzzz[k] * cd_y[k] + g_x_0_yzzzz_xyzzzz[k];

                g_x_0_yyzzzz_yyyyy[k] = -g_x_0_yzzzz_yyyyy[k] * cd_y[k] + g_x_0_yzzzz_yyyyyy[k];

                g_x_0_yyzzzz_yyyyz[k] = -g_x_0_yzzzz_yyyyz[k] * cd_y[k] + g_x_0_yzzzz_yyyyyz[k];

                g_x_0_yyzzzz_yyyzz[k] = -g_x_0_yzzzz_yyyzz[k] * cd_y[k] + g_x_0_yzzzz_yyyyzz[k];

                g_x_0_yyzzzz_yyzzz[k] = -g_x_0_yzzzz_yyzzz[k] * cd_y[k] + g_x_0_yzzzz_yyyzzz[k];

                g_x_0_yyzzzz_yzzzz[k] = -g_x_0_yzzzz_yzzzz[k] * cd_y[k] + g_x_0_yzzzz_yyzzzz[k];

                g_x_0_yyzzzz_zzzzz[k] = -g_x_0_yzzzz_zzzzz[k] * cd_y[k] + g_x_0_yzzzz_yzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 546);

            auto g_x_0_yzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 547);

            auto g_x_0_yzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 548);

            auto g_x_0_yzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 549);

            auto g_x_0_yzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 550);

            auto g_x_0_yzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 551);

            auto g_x_0_yzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 552);

            auto g_x_0_yzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 553);

            auto g_x_0_yzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 554);

            auto g_x_0_yzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 555);

            auto g_x_0_yzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 556);

            auto g_x_0_yzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 557);

            auto g_x_0_yzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 558);

            auto g_x_0_yzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 559);

            auto g_x_0_yzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 560);

            auto g_x_0_yzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 561);

            auto g_x_0_yzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 562);

            auto g_x_0_yzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 563);

            auto g_x_0_yzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 564);

            auto g_x_0_yzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 565);

            auto g_x_0_yzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 566);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzzz_xxxxx, g_x_0_yzzzzz_xxxxy, g_x_0_yzzzzz_xxxxz, g_x_0_yzzzzz_xxxyy, g_x_0_yzzzzz_xxxyz, g_x_0_yzzzzz_xxxzz, g_x_0_yzzzzz_xxyyy, g_x_0_yzzzzz_xxyyz, g_x_0_yzzzzz_xxyzz, g_x_0_yzzzzz_xxzzz, g_x_0_yzzzzz_xyyyy, g_x_0_yzzzzz_xyyyz, g_x_0_yzzzzz_xyyzz, g_x_0_yzzzzz_xyzzz, g_x_0_yzzzzz_xzzzz, g_x_0_yzzzzz_yyyyy, g_x_0_yzzzzz_yyyyz, g_x_0_yzzzzz_yyyzz, g_x_0_yzzzzz_yyzzz, g_x_0_yzzzzz_yzzzz, g_x_0_yzzzzz_zzzzz, g_x_0_zzzzz_xxxxx, g_x_0_zzzzz_xxxxxy, g_x_0_zzzzz_xxxxy, g_x_0_zzzzz_xxxxyy, g_x_0_zzzzz_xxxxyz, g_x_0_zzzzz_xxxxz, g_x_0_zzzzz_xxxyy, g_x_0_zzzzz_xxxyyy, g_x_0_zzzzz_xxxyyz, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxyzz, g_x_0_zzzzz_xxxzz, g_x_0_zzzzz_xxyyy, g_x_0_zzzzz_xxyyyy, g_x_0_zzzzz_xxyyyz, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyyzz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxyzzz, g_x_0_zzzzz_xxzzz, g_x_0_zzzzz_xyyyy, g_x_0_zzzzz_xyyyyy, g_x_0_zzzzz_xyyyyz, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyyzz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyyzzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xyzzzz, g_x_0_zzzzz_xzzzz, g_x_0_zzzzz_yyyyy, g_x_0_zzzzz_yyyyyy, g_x_0_zzzzz_yyyyyz, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyyzz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyyzzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yyzzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_yzzzzz, g_x_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_xxxxx[k] = -g_x_0_zzzzz_xxxxx[k] * cd_y[k] + g_x_0_zzzzz_xxxxxy[k];

                g_x_0_yzzzzz_xxxxy[k] = -g_x_0_zzzzz_xxxxy[k] * cd_y[k] + g_x_0_zzzzz_xxxxyy[k];

                g_x_0_yzzzzz_xxxxz[k] = -g_x_0_zzzzz_xxxxz[k] * cd_y[k] + g_x_0_zzzzz_xxxxyz[k];

                g_x_0_yzzzzz_xxxyy[k] = -g_x_0_zzzzz_xxxyy[k] * cd_y[k] + g_x_0_zzzzz_xxxyyy[k];

                g_x_0_yzzzzz_xxxyz[k] = -g_x_0_zzzzz_xxxyz[k] * cd_y[k] + g_x_0_zzzzz_xxxyyz[k];

                g_x_0_yzzzzz_xxxzz[k] = -g_x_0_zzzzz_xxxzz[k] * cd_y[k] + g_x_0_zzzzz_xxxyzz[k];

                g_x_0_yzzzzz_xxyyy[k] = -g_x_0_zzzzz_xxyyy[k] * cd_y[k] + g_x_0_zzzzz_xxyyyy[k];

                g_x_0_yzzzzz_xxyyz[k] = -g_x_0_zzzzz_xxyyz[k] * cd_y[k] + g_x_0_zzzzz_xxyyyz[k];

                g_x_0_yzzzzz_xxyzz[k] = -g_x_0_zzzzz_xxyzz[k] * cd_y[k] + g_x_0_zzzzz_xxyyzz[k];

                g_x_0_yzzzzz_xxzzz[k] = -g_x_0_zzzzz_xxzzz[k] * cd_y[k] + g_x_0_zzzzz_xxyzzz[k];

                g_x_0_yzzzzz_xyyyy[k] = -g_x_0_zzzzz_xyyyy[k] * cd_y[k] + g_x_0_zzzzz_xyyyyy[k];

                g_x_0_yzzzzz_xyyyz[k] = -g_x_0_zzzzz_xyyyz[k] * cd_y[k] + g_x_0_zzzzz_xyyyyz[k];

                g_x_0_yzzzzz_xyyzz[k] = -g_x_0_zzzzz_xyyzz[k] * cd_y[k] + g_x_0_zzzzz_xyyyzz[k];

                g_x_0_yzzzzz_xyzzz[k] = -g_x_0_zzzzz_xyzzz[k] * cd_y[k] + g_x_0_zzzzz_xyyzzz[k];

                g_x_0_yzzzzz_xzzzz[k] = -g_x_0_zzzzz_xzzzz[k] * cd_y[k] + g_x_0_zzzzz_xyzzzz[k];

                g_x_0_yzzzzz_yyyyy[k] = -g_x_0_zzzzz_yyyyy[k] * cd_y[k] + g_x_0_zzzzz_yyyyyy[k];

                g_x_0_yzzzzz_yyyyz[k] = -g_x_0_zzzzz_yyyyz[k] * cd_y[k] + g_x_0_zzzzz_yyyyyz[k];

                g_x_0_yzzzzz_yyyzz[k] = -g_x_0_zzzzz_yyyzz[k] * cd_y[k] + g_x_0_zzzzz_yyyyzz[k];

                g_x_0_yzzzzz_yyzzz[k] = -g_x_0_zzzzz_yyzzz[k] * cd_y[k] + g_x_0_zzzzz_yyyzzz[k];

                g_x_0_yzzzzz_yzzzz[k] = -g_x_0_zzzzz_yzzzz[k] * cd_y[k] + g_x_0_zzzzz_yyzzzz[k];

                g_x_0_yzzzzz_zzzzz[k] = -g_x_0_zzzzz_zzzzz[k] * cd_y[k] + g_x_0_zzzzz_yzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 567);

            auto g_x_0_zzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 568);

            auto g_x_0_zzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 569);

            auto g_x_0_zzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 570);

            auto g_x_0_zzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 571);

            auto g_x_0_zzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 572);

            auto g_x_0_zzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 573);

            auto g_x_0_zzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 574);

            auto g_x_0_zzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 575);

            auto g_x_0_zzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 576);

            auto g_x_0_zzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 577);

            auto g_x_0_zzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 578);

            auto g_x_0_zzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 579);

            auto g_x_0_zzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 580);

            auto g_x_0_zzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 581);

            auto g_x_0_zzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 582);

            auto g_x_0_zzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 583);

            auto g_x_0_zzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 584);

            auto g_x_0_zzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 585);

            auto g_x_0_zzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 586);

            auto g_x_0_zzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 0 * acomps * bcomps + 587);

            #pragma omp simd aligned(cd_z, g_x_0_zzzzz_xxxxx, g_x_0_zzzzz_xxxxxz, g_x_0_zzzzz_xxxxy, g_x_0_zzzzz_xxxxyz, g_x_0_zzzzz_xxxxz, g_x_0_zzzzz_xxxxzz, g_x_0_zzzzz_xxxyy, g_x_0_zzzzz_xxxyyz, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxyzz, g_x_0_zzzzz_xxxzz, g_x_0_zzzzz_xxxzzz, g_x_0_zzzzz_xxyyy, g_x_0_zzzzz_xxyyyz, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyyzz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxyzzz, g_x_0_zzzzz_xxzzz, g_x_0_zzzzz_xxzzzz, g_x_0_zzzzz_xyyyy, g_x_0_zzzzz_xyyyyz, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyyzz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyyzzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xyzzzz, g_x_0_zzzzz_xzzzz, g_x_0_zzzzz_xzzzzz, g_x_0_zzzzz_yyyyy, g_x_0_zzzzz_yyyyyz, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyyzz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyyzzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yyzzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_yzzzzz, g_x_0_zzzzz_zzzzz, g_x_0_zzzzz_zzzzzz, g_x_0_zzzzzz_xxxxx, g_x_0_zzzzzz_xxxxy, g_x_0_zzzzzz_xxxxz, g_x_0_zzzzzz_xxxyy, g_x_0_zzzzzz_xxxyz, g_x_0_zzzzzz_xxxzz, g_x_0_zzzzzz_xxyyy, g_x_0_zzzzzz_xxyyz, g_x_0_zzzzzz_xxyzz, g_x_0_zzzzzz_xxzzz, g_x_0_zzzzzz_xyyyy, g_x_0_zzzzzz_xyyyz, g_x_0_zzzzzz_xyyzz, g_x_0_zzzzzz_xyzzz, g_x_0_zzzzzz_xzzzz, g_x_0_zzzzzz_yyyyy, g_x_0_zzzzzz_yyyyz, g_x_0_zzzzzz_yyyzz, g_x_0_zzzzzz_yyzzz, g_x_0_zzzzzz_yzzzz, g_x_0_zzzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_xxxxx[k] = -g_x_0_zzzzz_xxxxx[k] * cd_z[k] + g_x_0_zzzzz_xxxxxz[k];

                g_x_0_zzzzzz_xxxxy[k] = -g_x_0_zzzzz_xxxxy[k] * cd_z[k] + g_x_0_zzzzz_xxxxyz[k];

                g_x_0_zzzzzz_xxxxz[k] = -g_x_0_zzzzz_xxxxz[k] * cd_z[k] + g_x_0_zzzzz_xxxxzz[k];

                g_x_0_zzzzzz_xxxyy[k] = -g_x_0_zzzzz_xxxyy[k] * cd_z[k] + g_x_0_zzzzz_xxxyyz[k];

                g_x_0_zzzzzz_xxxyz[k] = -g_x_0_zzzzz_xxxyz[k] * cd_z[k] + g_x_0_zzzzz_xxxyzz[k];

                g_x_0_zzzzzz_xxxzz[k] = -g_x_0_zzzzz_xxxzz[k] * cd_z[k] + g_x_0_zzzzz_xxxzzz[k];

                g_x_0_zzzzzz_xxyyy[k] = -g_x_0_zzzzz_xxyyy[k] * cd_z[k] + g_x_0_zzzzz_xxyyyz[k];

                g_x_0_zzzzzz_xxyyz[k] = -g_x_0_zzzzz_xxyyz[k] * cd_z[k] + g_x_0_zzzzz_xxyyzz[k];

                g_x_0_zzzzzz_xxyzz[k] = -g_x_0_zzzzz_xxyzz[k] * cd_z[k] + g_x_0_zzzzz_xxyzzz[k];

                g_x_0_zzzzzz_xxzzz[k] = -g_x_0_zzzzz_xxzzz[k] * cd_z[k] + g_x_0_zzzzz_xxzzzz[k];

                g_x_0_zzzzzz_xyyyy[k] = -g_x_0_zzzzz_xyyyy[k] * cd_z[k] + g_x_0_zzzzz_xyyyyz[k];

                g_x_0_zzzzzz_xyyyz[k] = -g_x_0_zzzzz_xyyyz[k] * cd_z[k] + g_x_0_zzzzz_xyyyzz[k];

                g_x_0_zzzzzz_xyyzz[k] = -g_x_0_zzzzz_xyyzz[k] * cd_z[k] + g_x_0_zzzzz_xyyzzz[k];

                g_x_0_zzzzzz_xyzzz[k] = -g_x_0_zzzzz_xyzzz[k] * cd_z[k] + g_x_0_zzzzz_xyzzzz[k];

                g_x_0_zzzzzz_xzzzz[k] = -g_x_0_zzzzz_xzzzz[k] * cd_z[k] + g_x_0_zzzzz_xzzzzz[k];

                g_x_0_zzzzzz_yyyyy[k] = -g_x_0_zzzzz_yyyyy[k] * cd_z[k] + g_x_0_zzzzz_yyyyyz[k];

                g_x_0_zzzzzz_yyyyz[k] = -g_x_0_zzzzz_yyyyz[k] * cd_z[k] + g_x_0_zzzzz_yyyyzz[k];

                g_x_0_zzzzzz_yyyzz[k] = -g_x_0_zzzzz_yyyzz[k] * cd_z[k] + g_x_0_zzzzz_yyyzzz[k];

                g_x_0_zzzzzz_yyzzz[k] = -g_x_0_zzzzz_yyzzz[k] * cd_z[k] + g_x_0_zzzzz_yyzzzz[k];

                g_x_0_zzzzzz_yzzzz[k] = -g_x_0_zzzzz_yzzzz[k] * cd_z[k] + g_x_0_zzzzz_yzzzzz[k];

                g_x_0_zzzzzz_zzzzz[k] = -g_x_0_zzzzz_zzzzz[k] * cd_z[k] + g_x_0_zzzzz_zzzzzz[k];
            }
            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 0);

            auto g_y_0_xxxxxx_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 1);

            auto g_y_0_xxxxxx_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 2);

            auto g_y_0_xxxxxx_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 3);

            auto g_y_0_xxxxxx_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 4);

            auto g_y_0_xxxxxx_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 5);

            auto g_y_0_xxxxxx_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 6);

            auto g_y_0_xxxxxx_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 7);

            auto g_y_0_xxxxxx_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 8);

            auto g_y_0_xxxxxx_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 9);

            auto g_y_0_xxxxxx_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 10);

            auto g_y_0_xxxxxx_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 11);

            auto g_y_0_xxxxxx_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 12);

            auto g_y_0_xxxxxx_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 13);

            auto g_y_0_xxxxxx_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 14);

            auto g_y_0_xxxxxx_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 15);

            auto g_y_0_xxxxxx_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 16);

            auto g_y_0_xxxxxx_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 17);

            auto g_y_0_xxxxxx_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 18);

            auto g_y_0_xxxxxx_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 19);

            auto g_y_0_xxxxxx_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxx_xxxxx, g_y_0_xxxxx_xxxxxx, g_y_0_xxxxx_xxxxxy, g_y_0_xxxxx_xxxxxz, g_y_0_xxxxx_xxxxy, g_y_0_xxxxx_xxxxyy, g_y_0_xxxxx_xxxxyz, g_y_0_xxxxx_xxxxz, g_y_0_xxxxx_xxxxzz, g_y_0_xxxxx_xxxyy, g_y_0_xxxxx_xxxyyy, g_y_0_xxxxx_xxxyyz, g_y_0_xxxxx_xxxyz, g_y_0_xxxxx_xxxyzz, g_y_0_xxxxx_xxxzz, g_y_0_xxxxx_xxxzzz, g_y_0_xxxxx_xxyyy, g_y_0_xxxxx_xxyyyy, g_y_0_xxxxx_xxyyyz, g_y_0_xxxxx_xxyyz, g_y_0_xxxxx_xxyyzz, g_y_0_xxxxx_xxyzz, g_y_0_xxxxx_xxyzzz, g_y_0_xxxxx_xxzzz, g_y_0_xxxxx_xxzzzz, g_y_0_xxxxx_xyyyy, g_y_0_xxxxx_xyyyyy, g_y_0_xxxxx_xyyyyz, g_y_0_xxxxx_xyyyz, g_y_0_xxxxx_xyyyzz, g_y_0_xxxxx_xyyzz, g_y_0_xxxxx_xyyzzz, g_y_0_xxxxx_xyzzz, g_y_0_xxxxx_xyzzzz, g_y_0_xxxxx_xzzzz, g_y_0_xxxxx_xzzzzz, g_y_0_xxxxx_yyyyy, g_y_0_xxxxx_yyyyz, g_y_0_xxxxx_yyyzz, g_y_0_xxxxx_yyzzz, g_y_0_xxxxx_yzzzz, g_y_0_xxxxx_zzzzz, g_y_0_xxxxxx_xxxxx, g_y_0_xxxxxx_xxxxy, g_y_0_xxxxxx_xxxxz, g_y_0_xxxxxx_xxxyy, g_y_0_xxxxxx_xxxyz, g_y_0_xxxxxx_xxxzz, g_y_0_xxxxxx_xxyyy, g_y_0_xxxxxx_xxyyz, g_y_0_xxxxxx_xxyzz, g_y_0_xxxxxx_xxzzz, g_y_0_xxxxxx_xyyyy, g_y_0_xxxxxx_xyyyz, g_y_0_xxxxxx_xyyzz, g_y_0_xxxxxx_xyzzz, g_y_0_xxxxxx_xzzzz, g_y_0_xxxxxx_yyyyy, g_y_0_xxxxxx_yyyyz, g_y_0_xxxxxx_yyyzz, g_y_0_xxxxxx_yyzzz, g_y_0_xxxxxx_yzzzz, g_y_0_xxxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_xxxxx[k] = -g_y_0_xxxxx_xxxxx[k] * cd_x[k] + g_y_0_xxxxx_xxxxxx[k];

                g_y_0_xxxxxx_xxxxy[k] = -g_y_0_xxxxx_xxxxy[k] * cd_x[k] + g_y_0_xxxxx_xxxxxy[k];

                g_y_0_xxxxxx_xxxxz[k] = -g_y_0_xxxxx_xxxxz[k] * cd_x[k] + g_y_0_xxxxx_xxxxxz[k];

                g_y_0_xxxxxx_xxxyy[k] = -g_y_0_xxxxx_xxxyy[k] * cd_x[k] + g_y_0_xxxxx_xxxxyy[k];

                g_y_0_xxxxxx_xxxyz[k] = -g_y_0_xxxxx_xxxyz[k] * cd_x[k] + g_y_0_xxxxx_xxxxyz[k];

                g_y_0_xxxxxx_xxxzz[k] = -g_y_0_xxxxx_xxxzz[k] * cd_x[k] + g_y_0_xxxxx_xxxxzz[k];

                g_y_0_xxxxxx_xxyyy[k] = -g_y_0_xxxxx_xxyyy[k] * cd_x[k] + g_y_0_xxxxx_xxxyyy[k];

                g_y_0_xxxxxx_xxyyz[k] = -g_y_0_xxxxx_xxyyz[k] * cd_x[k] + g_y_0_xxxxx_xxxyyz[k];

                g_y_0_xxxxxx_xxyzz[k] = -g_y_0_xxxxx_xxyzz[k] * cd_x[k] + g_y_0_xxxxx_xxxyzz[k];

                g_y_0_xxxxxx_xxzzz[k] = -g_y_0_xxxxx_xxzzz[k] * cd_x[k] + g_y_0_xxxxx_xxxzzz[k];

                g_y_0_xxxxxx_xyyyy[k] = -g_y_0_xxxxx_xyyyy[k] * cd_x[k] + g_y_0_xxxxx_xxyyyy[k];

                g_y_0_xxxxxx_xyyyz[k] = -g_y_0_xxxxx_xyyyz[k] * cd_x[k] + g_y_0_xxxxx_xxyyyz[k];

                g_y_0_xxxxxx_xyyzz[k] = -g_y_0_xxxxx_xyyzz[k] * cd_x[k] + g_y_0_xxxxx_xxyyzz[k];

                g_y_0_xxxxxx_xyzzz[k] = -g_y_0_xxxxx_xyzzz[k] * cd_x[k] + g_y_0_xxxxx_xxyzzz[k];

                g_y_0_xxxxxx_xzzzz[k] = -g_y_0_xxxxx_xzzzz[k] * cd_x[k] + g_y_0_xxxxx_xxzzzz[k];

                g_y_0_xxxxxx_yyyyy[k] = -g_y_0_xxxxx_yyyyy[k] * cd_x[k] + g_y_0_xxxxx_xyyyyy[k];

                g_y_0_xxxxxx_yyyyz[k] = -g_y_0_xxxxx_yyyyz[k] * cd_x[k] + g_y_0_xxxxx_xyyyyz[k];

                g_y_0_xxxxxx_yyyzz[k] = -g_y_0_xxxxx_yyyzz[k] * cd_x[k] + g_y_0_xxxxx_xyyyzz[k];

                g_y_0_xxxxxx_yyzzz[k] = -g_y_0_xxxxx_yyzzz[k] * cd_x[k] + g_y_0_xxxxx_xyyzzz[k];

                g_y_0_xxxxxx_yzzzz[k] = -g_y_0_xxxxx_yzzzz[k] * cd_x[k] + g_y_0_xxxxx_xyzzzz[k];

                g_y_0_xxxxxx_zzzzz[k] = -g_y_0_xxxxx_zzzzz[k] * cd_x[k] + g_y_0_xxxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 21);

            auto g_y_0_xxxxxy_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 22);

            auto g_y_0_xxxxxy_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 23);

            auto g_y_0_xxxxxy_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 24);

            auto g_y_0_xxxxxy_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 25);

            auto g_y_0_xxxxxy_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 26);

            auto g_y_0_xxxxxy_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 27);

            auto g_y_0_xxxxxy_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 28);

            auto g_y_0_xxxxxy_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 29);

            auto g_y_0_xxxxxy_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 30);

            auto g_y_0_xxxxxy_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 31);

            auto g_y_0_xxxxxy_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 32);

            auto g_y_0_xxxxxy_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 33);

            auto g_y_0_xxxxxy_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 34);

            auto g_y_0_xxxxxy_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 35);

            auto g_y_0_xxxxxy_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 36);

            auto g_y_0_xxxxxy_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 37);

            auto g_y_0_xxxxxy_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 38);

            auto g_y_0_xxxxxy_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 39);

            auto g_y_0_xxxxxy_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 40);

            auto g_y_0_xxxxxy_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxy_xxxxx, g_y_0_xxxxxy_xxxxy, g_y_0_xxxxxy_xxxxz, g_y_0_xxxxxy_xxxyy, g_y_0_xxxxxy_xxxyz, g_y_0_xxxxxy_xxxzz, g_y_0_xxxxxy_xxyyy, g_y_0_xxxxxy_xxyyz, g_y_0_xxxxxy_xxyzz, g_y_0_xxxxxy_xxzzz, g_y_0_xxxxxy_xyyyy, g_y_0_xxxxxy_xyyyz, g_y_0_xxxxxy_xyyzz, g_y_0_xxxxxy_xyzzz, g_y_0_xxxxxy_xzzzz, g_y_0_xxxxxy_yyyyy, g_y_0_xxxxxy_yyyyz, g_y_0_xxxxxy_yyyzz, g_y_0_xxxxxy_yyzzz, g_y_0_xxxxxy_yzzzz, g_y_0_xxxxxy_zzzzz, g_y_0_xxxxy_xxxxx, g_y_0_xxxxy_xxxxxx, g_y_0_xxxxy_xxxxxy, g_y_0_xxxxy_xxxxxz, g_y_0_xxxxy_xxxxy, g_y_0_xxxxy_xxxxyy, g_y_0_xxxxy_xxxxyz, g_y_0_xxxxy_xxxxz, g_y_0_xxxxy_xxxxzz, g_y_0_xxxxy_xxxyy, g_y_0_xxxxy_xxxyyy, g_y_0_xxxxy_xxxyyz, g_y_0_xxxxy_xxxyz, g_y_0_xxxxy_xxxyzz, g_y_0_xxxxy_xxxzz, g_y_0_xxxxy_xxxzzz, g_y_0_xxxxy_xxyyy, g_y_0_xxxxy_xxyyyy, g_y_0_xxxxy_xxyyyz, g_y_0_xxxxy_xxyyz, g_y_0_xxxxy_xxyyzz, g_y_0_xxxxy_xxyzz, g_y_0_xxxxy_xxyzzz, g_y_0_xxxxy_xxzzz, g_y_0_xxxxy_xxzzzz, g_y_0_xxxxy_xyyyy, g_y_0_xxxxy_xyyyyy, g_y_0_xxxxy_xyyyyz, g_y_0_xxxxy_xyyyz, g_y_0_xxxxy_xyyyzz, g_y_0_xxxxy_xyyzz, g_y_0_xxxxy_xyyzzz, g_y_0_xxxxy_xyzzz, g_y_0_xxxxy_xyzzzz, g_y_0_xxxxy_xzzzz, g_y_0_xxxxy_xzzzzz, g_y_0_xxxxy_yyyyy, g_y_0_xxxxy_yyyyz, g_y_0_xxxxy_yyyzz, g_y_0_xxxxy_yyzzz, g_y_0_xxxxy_yzzzz, g_y_0_xxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_xxxxx[k] = -g_y_0_xxxxy_xxxxx[k] * cd_x[k] + g_y_0_xxxxy_xxxxxx[k];

                g_y_0_xxxxxy_xxxxy[k] = -g_y_0_xxxxy_xxxxy[k] * cd_x[k] + g_y_0_xxxxy_xxxxxy[k];

                g_y_0_xxxxxy_xxxxz[k] = -g_y_0_xxxxy_xxxxz[k] * cd_x[k] + g_y_0_xxxxy_xxxxxz[k];

                g_y_0_xxxxxy_xxxyy[k] = -g_y_0_xxxxy_xxxyy[k] * cd_x[k] + g_y_0_xxxxy_xxxxyy[k];

                g_y_0_xxxxxy_xxxyz[k] = -g_y_0_xxxxy_xxxyz[k] * cd_x[k] + g_y_0_xxxxy_xxxxyz[k];

                g_y_0_xxxxxy_xxxzz[k] = -g_y_0_xxxxy_xxxzz[k] * cd_x[k] + g_y_0_xxxxy_xxxxzz[k];

                g_y_0_xxxxxy_xxyyy[k] = -g_y_0_xxxxy_xxyyy[k] * cd_x[k] + g_y_0_xxxxy_xxxyyy[k];

                g_y_0_xxxxxy_xxyyz[k] = -g_y_0_xxxxy_xxyyz[k] * cd_x[k] + g_y_0_xxxxy_xxxyyz[k];

                g_y_0_xxxxxy_xxyzz[k] = -g_y_0_xxxxy_xxyzz[k] * cd_x[k] + g_y_0_xxxxy_xxxyzz[k];

                g_y_0_xxxxxy_xxzzz[k] = -g_y_0_xxxxy_xxzzz[k] * cd_x[k] + g_y_0_xxxxy_xxxzzz[k];

                g_y_0_xxxxxy_xyyyy[k] = -g_y_0_xxxxy_xyyyy[k] * cd_x[k] + g_y_0_xxxxy_xxyyyy[k];

                g_y_0_xxxxxy_xyyyz[k] = -g_y_0_xxxxy_xyyyz[k] * cd_x[k] + g_y_0_xxxxy_xxyyyz[k];

                g_y_0_xxxxxy_xyyzz[k] = -g_y_0_xxxxy_xyyzz[k] * cd_x[k] + g_y_0_xxxxy_xxyyzz[k];

                g_y_0_xxxxxy_xyzzz[k] = -g_y_0_xxxxy_xyzzz[k] * cd_x[k] + g_y_0_xxxxy_xxyzzz[k];

                g_y_0_xxxxxy_xzzzz[k] = -g_y_0_xxxxy_xzzzz[k] * cd_x[k] + g_y_0_xxxxy_xxzzzz[k];

                g_y_0_xxxxxy_yyyyy[k] = -g_y_0_xxxxy_yyyyy[k] * cd_x[k] + g_y_0_xxxxy_xyyyyy[k];

                g_y_0_xxxxxy_yyyyz[k] = -g_y_0_xxxxy_yyyyz[k] * cd_x[k] + g_y_0_xxxxy_xyyyyz[k];

                g_y_0_xxxxxy_yyyzz[k] = -g_y_0_xxxxy_yyyzz[k] * cd_x[k] + g_y_0_xxxxy_xyyyzz[k];

                g_y_0_xxxxxy_yyzzz[k] = -g_y_0_xxxxy_yyzzz[k] * cd_x[k] + g_y_0_xxxxy_xyyzzz[k];

                g_y_0_xxxxxy_yzzzz[k] = -g_y_0_xxxxy_yzzzz[k] * cd_x[k] + g_y_0_xxxxy_xyzzzz[k];

                g_y_0_xxxxxy_zzzzz[k] = -g_y_0_xxxxy_zzzzz[k] * cd_x[k] + g_y_0_xxxxy_xzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 42);

            auto g_y_0_xxxxxz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 43);

            auto g_y_0_xxxxxz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 44);

            auto g_y_0_xxxxxz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 45);

            auto g_y_0_xxxxxz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 46);

            auto g_y_0_xxxxxz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 47);

            auto g_y_0_xxxxxz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 48);

            auto g_y_0_xxxxxz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 49);

            auto g_y_0_xxxxxz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 50);

            auto g_y_0_xxxxxz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 51);

            auto g_y_0_xxxxxz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 52);

            auto g_y_0_xxxxxz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 53);

            auto g_y_0_xxxxxz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 54);

            auto g_y_0_xxxxxz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 55);

            auto g_y_0_xxxxxz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 56);

            auto g_y_0_xxxxxz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 57);

            auto g_y_0_xxxxxz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 58);

            auto g_y_0_xxxxxz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 59);

            auto g_y_0_xxxxxz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 60);

            auto g_y_0_xxxxxz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 61);

            auto g_y_0_xxxxxz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxz_xxxxx, g_y_0_xxxxxz_xxxxy, g_y_0_xxxxxz_xxxxz, g_y_0_xxxxxz_xxxyy, g_y_0_xxxxxz_xxxyz, g_y_0_xxxxxz_xxxzz, g_y_0_xxxxxz_xxyyy, g_y_0_xxxxxz_xxyyz, g_y_0_xxxxxz_xxyzz, g_y_0_xxxxxz_xxzzz, g_y_0_xxxxxz_xyyyy, g_y_0_xxxxxz_xyyyz, g_y_0_xxxxxz_xyyzz, g_y_0_xxxxxz_xyzzz, g_y_0_xxxxxz_xzzzz, g_y_0_xxxxxz_yyyyy, g_y_0_xxxxxz_yyyyz, g_y_0_xxxxxz_yyyzz, g_y_0_xxxxxz_yyzzz, g_y_0_xxxxxz_yzzzz, g_y_0_xxxxxz_zzzzz, g_y_0_xxxxz_xxxxx, g_y_0_xxxxz_xxxxxx, g_y_0_xxxxz_xxxxxy, g_y_0_xxxxz_xxxxxz, g_y_0_xxxxz_xxxxy, g_y_0_xxxxz_xxxxyy, g_y_0_xxxxz_xxxxyz, g_y_0_xxxxz_xxxxz, g_y_0_xxxxz_xxxxzz, g_y_0_xxxxz_xxxyy, g_y_0_xxxxz_xxxyyy, g_y_0_xxxxz_xxxyyz, g_y_0_xxxxz_xxxyz, g_y_0_xxxxz_xxxyzz, g_y_0_xxxxz_xxxzz, g_y_0_xxxxz_xxxzzz, g_y_0_xxxxz_xxyyy, g_y_0_xxxxz_xxyyyy, g_y_0_xxxxz_xxyyyz, g_y_0_xxxxz_xxyyz, g_y_0_xxxxz_xxyyzz, g_y_0_xxxxz_xxyzz, g_y_0_xxxxz_xxyzzz, g_y_0_xxxxz_xxzzz, g_y_0_xxxxz_xxzzzz, g_y_0_xxxxz_xyyyy, g_y_0_xxxxz_xyyyyy, g_y_0_xxxxz_xyyyyz, g_y_0_xxxxz_xyyyz, g_y_0_xxxxz_xyyyzz, g_y_0_xxxxz_xyyzz, g_y_0_xxxxz_xyyzzz, g_y_0_xxxxz_xyzzz, g_y_0_xxxxz_xyzzzz, g_y_0_xxxxz_xzzzz, g_y_0_xxxxz_xzzzzz, g_y_0_xxxxz_yyyyy, g_y_0_xxxxz_yyyyz, g_y_0_xxxxz_yyyzz, g_y_0_xxxxz_yyzzz, g_y_0_xxxxz_yzzzz, g_y_0_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_xxxxx[k] = -g_y_0_xxxxz_xxxxx[k] * cd_x[k] + g_y_0_xxxxz_xxxxxx[k];

                g_y_0_xxxxxz_xxxxy[k] = -g_y_0_xxxxz_xxxxy[k] * cd_x[k] + g_y_0_xxxxz_xxxxxy[k];

                g_y_0_xxxxxz_xxxxz[k] = -g_y_0_xxxxz_xxxxz[k] * cd_x[k] + g_y_0_xxxxz_xxxxxz[k];

                g_y_0_xxxxxz_xxxyy[k] = -g_y_0_xxxxz_xxxyy[k] * cd_x[k] + g_y_0_xxxxz_xxxxyy[k];

                g_y_0_xxxxxz_xxxyz[k] = -g_y_0_xxxxz_xxxyz[k] * cd_x[k] + g_y_0_xxxxz_xxxxyz[k];

                g_y_0_xxxxxz_xxxzz[k] = -g_y_0_xxxxz_xxxzz[k] * cd_x[k] + g_y_0_xxxxz_xxxxzz[k];

                g_y_0_xxxxxz_xxyyy[k] = -g_y_0_xxxxz_xxyyy[k] * cd_x[k] + g_y_0_xxxxz_xxxyyy[k];

                g_y_0_xxxxxz_xxyyz[k] = -g_y_0_xxxxz_xxyyz[k] * cd_x[k] + g_y_0_xxxxz_xxxyyz[k];

                g_y_0_xxxxxz_xxyzz[k] = -g_y_0_xxxxz_xxyzz[k] * cd_x[k] + g_y_0_xxxxz_xxxyzz[k];

                g_y_0_xxxxxz_xxzzz[k] = -g_y_0_xxxxz_xxzzz[k] * cd_x[k] + g_y_0_xxxxz_xxxzzz[k];

                g_y_0_xxxxxz_xyyyy[k] = -g_y_0_xxxxz_xyyyy[k] * cd_x[k] + g_y_0_xxxxz_xxyyyy[k];

                g_y_0_xxxxxz_xyyyz[k] = -g_y_0_xxxxz_xyyyz[k] * cd_x[k] + g_y_0_xxxxz_xxyyyz[k];

                g_y_0_xxxxxz_xyyzz[k] = -g_y_0_xxxxz_xyyzz[k] * cd_x[k] + g_y_0_xxxxz_xxyyzz[k];

                g_y_0_xxxxxz_xyzzz[k] = -g_y_0_xxxxz_xyzzz[k] * cd_x[k] + g_y_0_xxxxz_xxyzzz[k];

                g_y_0_xxxxxz_xzzzz[k] = -g_y_0_xxxxz_xzzzz[k] * cd_x[k] + g_y_0_xxxxz_xxzzzz[k];

                g_y_0_xxxxxz_yyyyy[k] = -g_y_0_xxxxz_yyyyy[k] * cd_x[k] + g_y_0_xxxxz_xyyyyy[k];

                g_y_0_xxxxxz_yyyyz[k] = -g_y_0_xxxxz_yyyyz[k] * cd_x[k] + g_y_0_xxxxz_xyyyyz[k];

                g_y_0_xxxxxz_yyyzz[k] = -g_y_0_xxxxz_yyyzz[k] * cd_x[k] + g_y_0_xxxxz_xyyyzz[k];

                g_y_0_xxxxxz_yyzzz[k] = -g_y_0_xxxxz_yyzzz[k] * cd_x[k] + g_y_0_xxxxz_xyyzzz[k];

                g_y_0_xxxxxz_yzzzz[k] = -g_y_0_xxxxz_yzzzz[k] * cd_x[k] + g_y_0_xxxxz_xyzzzz[k];

                g_y_0_xxxxxz_zzzzz[k] = -g_y_0_xxxxz_zzzzz[k] * cd_x[k] + g_y_0_xxxxz_xzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 63);

            auto g_y_0_xxxxyy_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 64);

            auto g_y_0_xxxxyy_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 65);

            auto g_y_0_xxxxyy_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 66);

            auto g_y_0_xxxxyy_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 67);

            auto g_y_0_xxxxyy_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 68);

            auto g_y_0_xxxxyy_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 69);

            auto g_y_0_xxxxyy_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 70);

            auto g_y_0_xxxxyy_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 71);

            auto g_y_0_xxxxyy_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 72);

            auto g_y_0_xxxxyy_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 73);

            auto g_y_0_xxxxyy_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 74);

            auto g_y_0_xxxxyy_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 75);

            auto g_y_0_xxxxyy_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 76);

            auto g_y_0_xxxxyy_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 77);

            auto g_y_0_xxxxyy_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 78);

            auto g_y_0_xxxxyy_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 79);

            auto g_y_0_xxxxyy_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 80);

            auto g_y_0_xxxxyy_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 81);

            auto g_y_0_xxxxyy_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 82);

            auto g_y_0_xxxxyy_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyy_xxxxx, g_y_0_xxxxyy_xxxxy, g_y_0_xxxxyy_xxxxz, g_y_0_xxxxyy_xxxyy, g_y_0_xxxxyy_xxxyz, g_y_0_xxxxyy_xxxzz, g_y_0_xxxxyy_xxyyy, g_y_0_xxxxyy_xxyyz, g_y_0_xxxxyy_xxyzz, g_y_0_xxxxyy_xxzzz, g_y_0_xxxxyy_xyyyy, g_y_0_xxxxyy_xyyyz, g_y_0_xxxxyy_xyyzz, g_y_0_xxxxyy_xyzzz, g_y_0_xxxxyy_xzzzz, g_y_0_xxxxyy_yyyyy, g_y_0_xxxxyy_yyyyz, g_y_0_xxxxyy_yyyzz, g_y_0_xxxxyy_yyzzz, g_y_0_xxxxyy_yzzzz, g_y_0_xxxxyy_zzzzz, g_y_0_xxxyy_xxxxx, g_y_0_xxxyy_xxxxxx, g_y_0_xxxyy_xxxxxy, g_y_0_xxxyy_xxxxxz, g_y_0_xxxyy_xxxxy, g_y_0_xxxyy_xxxxyy, g_y_0_xxxyy_xxxxyz, g_y_0_xxxyy_xxxxz, g_y_0_xxxyy_xxxxzz, g_y_0_xxxyy_xxxyy, g_y_0_xxxyy_xxxyyy, g_y_0_xxxyy_xxxyyz, g_y_0_xxxyy_xxxyz, g_y_0_xxxyy_xxxyzz, g_y_0_xxxyy_xxxzz, g_y_0_xxxyy_xxxzzz, g_y_0_xxxyy_xxyyy, g_y_0_xxxyy_xxyyyy, g_y_0_xxxyy_xxyyyz, g_y_0_xxxyy_xxyyz, g_y_0_xxxyy_xxyyzz, g_y_0_xxxyy_xxyzz, g_y_0_xxxyy_xxyzzz, g_y_0_xxxyy_xxzzz, g_y_0_xxxyy_xxzzzz, g_y_0_xxxyy_xyyyy, g_y_0_xxxyy_xyyyyy, g_y_0_xxxyy_xyyyyz, g_y_0_xxxyy_xyyyz, g_y_0_xxxyy_xyyyzz, g_y_0_xxxyy_xyyzz, g_y_0_xxxyy_xyyzzz, g_y_0_xxxyy_xyzzz, g_y_0_xxxyy_xyzzzz, g_y_0_xxxyy_xzzzz, g_y_0_xxxyy_xzzzzz, g_y_0_xxxyy_yyyyy, g_y_0_xxxyy_yyyyz, g_y_0_xxxyy_yyyzz, g_y_0_xxxyy_yyzzz, g_y_0_xxxyy_yzzzz, g_y_0_xxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_xxxxx[k] = -g_y_0_xxxyy_xxxxx[k] * cd_x[k] + g_y_0_xxxyy_xxxxxx[k];

                g_y_0_xxxxyy_xxxxy[k] = -g_y_0_xxxyy_xxxxy[k] * cd_x[k] + g_y_0_xxxyy_xxxxxy[k];

                g_y_0_xxxxyy_xxxxz[k] = -g_y_0_xxxyy_xxxxz[k] * cd_x[k] + g_y_0_xxxyy_xxxxxz[k];

                g_y_0_xxxxyy_xxxyy[k] = -g_y_0_xxxyy_xxxyy[k] * cd_x[k] + g_y_0_xxxyy_xxxxyy[k];

                g_y_0_xxxxyy_xxxyz[k] = -g_y_0_xxxyy_xxxyz[k] * cd_x[k] + g_y_0_xxxyy_xxxxyz[k];

                g_y_0_xxxxyy_xxxzz[k] = -g_y_0_xxxyy_xxxzz[k] * cd_x[k] + g_y_0_xxxyy_xxxxzz[k];

                g_y_0_xxxxyy_xxyyy[k] = -g_y_0_xxxyy_xxyyy[k] * cd_x[k] + g_y_0_xxxyy_xxxyyy[k];

                g_y_0_xxxxyy_xxyyz[k] = -g_y_0_xxxyy_xxyyz[k] * cd_x[k] + g_y_0_xxxyy_xxxyyz[k];

                g_y_0_xxxxyy_xxyzz[k] = -g_y_0_xxxyy_xxyzz[k] * cd_x[k] + g_y_0_xxxyy_xxxyzz[k];

                g_y_0_xxxxyy_xxzzz[k] = -g_y_0_xxxyy_xxzzz[k] * cd_x[k] + g_y_0_xxxyy_xxxzzz[k];

                g_y_0_xxxxyy_xyyyy[k] = -g_y_0_xxxyy_xyyyy[k] * cd_x[k] + g_y_0_xxxyy_xxyyyy[k];

                g_y_0_xxxxyy_xyyyz[k] = -g_y_0_xxxyy_xyyyz[k] * cd_x[k] + g_y_0_xxxyy_xxyyyz[k];

                g_y_0_xxxxyy_xyyzz[k] = -g_y_0_xxxyy_xyyzz[k] * cd_x[k] + g_y_0_xxxyy_xxyyzz[k];

                g_y_0_xxxxyy_xyzzz[k] = -g_y_0_xxxyy_xyzzz[k] * cd_x[k] + g_y_0_xxxyy_xxyzzz[k];

                g_y_0_xxxxyy_xzzzz[k] = -g_y_0_xxxyy_xzzzz[k] * cd_x[k] + g_y_0_xxxyy_xxzzzz[k];

                g_y_0_xxxxyy_yyyyy[k] = -g_y_0_xxxyy_yyyyy[k] * cd_x[k] + g_y_0_xxxyy_xyyyyy[k];

                g_y_0_xxxxyy_yyyyz[k] = -g_y_0_xxxyy_yyyyz[k] * cd_x[k] + g_y_0_xxxyy_xyyyyz[k];

                g_y_0_xxxxyy_yyyzz[k] = -g_y_0_xxxyy_yyyzz[k] * cd_x[k] + g_y_0_xxxyy_xyyyzz[k];

                g_y_0_xxxxyy_yyzzz[k] = -g_y_0_xxxyy_yyzzz[k] * cd_x[k] + g_y_0_xxxyy_xyyzzz[k];

                g_y_0_xxxxyy_yzzzz[k] = -g_y_0_xxxyy_yzzzz[k] * cd_x[k] + g_y_0_xxxyy_xyzzzz[k];

                g_y_0_xxxxyy_zzzzz[k] = -g_y_0_xxxyy_zzzzz[k] * cd_x[k] + g_y_0_xxxyy_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 84);

            auto g_y_0_xxxxyz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 85);

            auto g_y_0_xxxxyz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 86);

            auto g_y_0_xxxxyz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 87);

            auto g_y_0_xxxxyz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 88);

            auto g_y_0_xxxxyz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 89);

            auto g_y_0_xxxxyz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 90);

            auto g_y_0_xxxxyz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 91);

            auto g_y_0_xxxxyz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 92);

            auto g_y_0_xxxxyz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 93);

            auto g_y_0_xxxxyz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 94);

            auto g_y_0_xxxxyz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 95);

            auto g_y_0_xxxxyz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 96);

            auto g_y_0_xxxxyz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 97);

            auto g_y_0_xxxxyz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 98);

            auto g_y_0_xxxxyz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 99);

            auto g_y_0_xxxxyz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 100);

            auto g_y_0_xxxxyz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 101);

            auto g_y_0_xxxxyz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 102);

            auto g_y_0_xxxxyz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 103);

            auto g_y_0_xxxxyz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyz_xxxxx, g_y_0_xxxxyz_xxxxy, g_y_0_xxxxyz_xxxxz, g_y_0_xxxxyz_xxxyy, g_y_0_xxxxyz_xxxyz, g_y_0_xxxxyz_xxxzz, g_y_0_xxxxyz_xxyyy, g_y_0_xxxxyz_xxyyz, g_y_0_xxxxyz_xxyzz, g_y_0_xxxxyz_xxzzz, g_y_0_xxxxyz_xyyyy, g_y_0_xxxxyz_xyyyz, g_y_0_xxxxyz_xyyzz, g_y_0_xxxxyz_xyzzz, g_y_0_xxxxyz_xzzzz, g_y_0_xxxxyz_yyyyy, g_y_0_xxxxyz_yyyyz, g_y_0_xxxxyz_yyyzz, g_y_0_xxxxyz_yyzzz, g_y_0_xxxxyz_yzzzz, g_y_0_xxxxyz_zzzzz, g_y_0_xxxyz_xxxxx, g_y_0_xxxyz_xxxxxx, g_y_0_xxxyz_xxxxxy, g_y_0_xxxyz_xxxxxz, g_y_0_xxxyz_xxxxy, g_y_0_xxxyz_xxxxyy, g_y_0_xxxyz_xxxxyz, g_y_0_xxxyz_xxxxz, g_y_0_xxxyz_xxxxzz, g_y_0_xxxyz_xxxyy, g_y_0_xxxyz_xxxyyy, g_y_0_xxxyz_xxxyyz, g_y_0_xxxyz_xxxyz, g_y_0_xxxyz_xxxyzz, g_y_0_xxxyz_xxxzz, g_y_0_xxxyz_xxxzzz, g_y_0_xxxyz_xxyyy, g_y_0_xxxyz_xxyyyy, g_y_0_xxxyz_xxyyyz, g_y_0_xxxyz_xxyyz, g_y_0_xxxyz_xxyyzz, g_y_0_xxxyz_xxyzz, g_y_0_xxxyz_xxyzzz, g_y_0_xxxyz_xxzzz, g_y_0_xxxyz_xxzzzz, g_y_0_xxxyz_xyyyy, g_y_0_xxxyz_xyyyyy, g_y_0_xxxyz_xyyyyz, g_y_0_xxxyz_xyyyz, g_y_0_xxxyz_xyyyzz, g_y_0_xxxyz_xyyzz, g_y_0_xxxyz_xyyzzz, g_y_0_xxxyz_xyzzz, g_y_0_xxxyz_xyzzzz, g_y_0_xxxyz_xzzzz, g_y_0_xxxyz_xzzzzz, g_y_0_xxxyz_yyyyy, g_y_0_xxxyz_yyyyz, g_y_0_xxxyz_yyyzz, g_y_0_xxxyz_yyzzz, g_y_0_xxxyz_yzzzz, g_y_0_xxxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_xxxxx[k] = -g_y_0_xxxyz_xxxxx[k] * cd_x[k] + g_y_0_xxxyz_xxxxxx[k];

                g_y_0_xxxxyz_xxxxy[k] = -g_y_0_xxxyz_xxxxy[k] * cd_x[k] + g_y_0_xxxyz_xxxxxy[k];

                g_y_0_xxxxyz_xxxxz[k] = -g_y_0_xxxyz_xxxxz[k] * cd_x[k] + g_y_0_xxxyz_xxxxxz[k];

                g_y_0_xxxxyz_xxxyy[k] = -g_y_0_xxxyz_xxxyy[k] * cd_x[k] + g_y_0_xxxyz_xxxxyy[k];

                g_y_0_xxxxyz_xxxyz[k] = -g_y_0_xxxyz_xxxyz[k] * cd_x[k] + g_y_0_xxxyz_xxxxyz[k];

                g_y_0_xxxxyz_xxxzz[k] = -g_y_0_xxxyz_xxxzz[k] * cd_x[k] + g_y_0_xxxyz_xxxxzz[k];

                g_y_0_xxxxyz_xxyyy[k] = -g_y_0_xxxyz_xxyyy[k] * cd_x[k] + g_y_0_xxxyz_xxxyyy[k];

                g_y_0_xxxxyz_xxyyz[k] = -g_y_0_xxxyz_xxyyz[k] * cd_x[k] + g_y_0_xxxyz_xxxyyz[k];

                g_y_0_xxxxyz_xxyzz[k] = -g_y_0_xxxyz_xxyzz[k] * cd_x[k] + g_y_0_xxxyz_xxxyzz[k];

                g_y_0_xxxxyz_xxzzz[k] = -g_y_0_xxxyz_xxzzz[k] * cd_x[k] + g_y_0_xxxyz_xxxzzz[k];

                g_y_0_xxxxyz_xyyyy[k] = -g_y_0_xxxyz_xyyyy[k] * cd_x[k] + g_y_0_xxxyz_xxyyyy[k];

                g_y_0_xxxxyz_xyyyz[k] = -g_y_0_xxxyz_xyyyz[k] * cd_x[k] + g_y_0_xxxyz_xxyyyz[k];

                g_y_0_xxxxyz_xyyzz[k] = -g_y_0_xxxyz_xyyzz[k] * cd_x[k] + g_y_0_xxxyz_xxyyzz[k];

                g_y_0_xxxxyz_xyzzz[k] = -g_y_0_xxxyz_xyzzz[k] * cd_x[k] + g_y_0_xxxyz_xxyzzz[k];

                g_y_0_xxxxyz_xzzzz[k] = -g_y_0_xxxyz_xzzzz[k] * cd_x[k] + g_y_0_xxxyz_xxzzzz[k];

                g_y_0_xxxxyz_yyyyy[k] = -g_y_0_xxxyz_yyyyy[k] * cd_x[k] + g_y_0_xxxyz_xyyyyy[k];

                g_y_0_xxxxyz_yyyyz[k] = -g_y_0_xxxyz_yyyyz[k] * cd_x[k] + g_y_0_xxxyz_xyyyyz[k];

                g_y_0_xxxxyz_yyyzz[k] = -g_y_0_xxxyz_yyyzz[k] * cd_x[k] + g_y_0_xxxyz_xyyyzz[k];

                g_y_0_xxxxyz_yyzzz[k] = -g_y_0_xxxyz_yyzzz[k] * cd_x[k] + g_y_0_xxxyz_xyyzzz[k];

                g_y_0_xxxxyz_yzzzz[k] = -g_y_0_xxxyz_yzzzz[k] * cd_x[k] + g_y_0_xxxyz_xyzzzz[k];

                g_y_0_xxxxyz_zzzzz[k] = -g_y_0_xxxyz_zzzzz[k] * cd_x[k] + g_y_0_xxxyz_xzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 105);

            auto g_y_0_xxxxzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 106);

            auto g_y_0_xxxxzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 107);

            auto g_y_0_xxxxzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 108);

            auto g_y_0_xxxxzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 109);

            auto g_y_0_xxxxzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 110);

            auto g_y_0_xxxxzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 111);

            auto g_y_0_xxxxzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 112);

            auto g_y_0_xxxxzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 113);

            auto g_y_0_xxxxzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 114);

            auto g_y_0_xxxxzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 115);

            auto g_y_0_xxxxzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 116);

            auto g_y_0_xxxxzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 117);

            auto g_y_0_xxxxzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 118);

            auto g_y_0_xxxxzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 119);

            auto g_y_0_xxxxzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 120);

            auto g_y_0_xxxxzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 121);

            auto g_y_0_xxxxzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 122);

            auto g_y_0_xxxxzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 123);

            auto g_y_0_xxxxzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 124);

            auto g_y_0_xxxxzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxzz_xxxxx, g_y_0_xxxxzz_xxxxy, g_y_0_xxxxzz_xxxxz, g_y_0_xxxxzz_xxxyy, g_y_0_xxxxzz_xxxyz, g_y_0_xxxxzz_xxxzz, g_y_0_xxxxzz_xxyyy, g_y_0_xxxxzz_xxyyz, g_y_0_xxxxzz_xxyzz, g_y_0_xxxxzz_xxzzz, g_y_0_xxxxzz_xyyyy, g_y_0_xxxxzz_xyyyz, g_y_0_xxxxzz_xyyzz, g_y_0_xxxxzz_xyzzz, g_y_0_xxxxzz_xzzzz, g_y_0_xxxxzz_yyyyy, g_y_0_xxxxzz_yyyyz, g_y_0_xxxxzz_yyyzz, g_y_0_xxxxzz_yyzzz, g_y_0_xxxxzz_yzzzz, g_y_0_xxxxzz_zzzzz, g_y_0_xxxzz_xxxxx, g_y_0_xxxzz_xxxxxx, g_y_0_xxxzz_xxxxxy, g_y_0_xxxzz_xxxxxz, g_y_0_xxxzz_xxxxy, g_y_0_xxxzz_xxxxyy, g_y_0_xxxzz_xxxxyz, g_y_0_xxxzz_xxxxz, g_y_0_xxxzz_xxxxzz, g_y_0_xxxzz_xxxyy, g_y_0_xxxzz_xxxyyy, g_y_0_xxxzz_xxxyyz, g_y_0_xxxzz_xxxyz, g_y_0_xxxzz_xxxyzz, g_y_0_xxxzz_xxxzz, g_y_0_xxxzz_xxxzzz, g_y_0_xxxzz_xxyyy, g_y_0_xxxzz_xxyyyy, g_y_0_xxxzz_xxyyyz, g_y_0_xxxzz_xxyyz, g_y_0_xxxzz_xxyyzz, g_y_0_xxxzz_xxyzz, g_y_0_xxxzz_xxyzzz, g_y_0_xxxzz_xxzzz, g_y_0_xxxzz_xxzzzz, g_y_0_xxxzz_xyyyy, g_y_0_xxxzz_xyyyyy, g_y_0_xxxzz_xyyyyz, g_y_0_xxxzz_xyyyz, g_y_0_xxxzz_xyyyzz, g_y_0_xxxzz_xyyzz, g_y_0_xxxzz_xyyzzz, g_y_0_xxxzz_xyzzz, g_y_0_xxxzz_xyzzzz, g_y_0_xxxzz_xzzzz, g_y_0_xxxzz_xzzzzz, g_y_0_xxxzz_yyyyy, g_y_0_xxxzz_yyyyz, g_y_0_xxxzz_yyyzz, g_y_0_xxxzz_yyzzz, g_y_0_xxxzz_yzzzz, g_y_0_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_xxxxx[k] = -g_y_0_xxxzz_xxxxx[k] * cd_x[k] + g_y_0_xxxzz_xxxxxx[k];

                g_y_0_xxxxzz_xxxxy[k] = -g_y_0_xxxzz_xxxxy[k] * cd_x[k] + g_y_0_xxxzz_xxxxxy[k];

                g_y_0_xxxxzz_xxxxz[k] = -g_y_0_xxxzz_xxxxz[k] * cd_x[k] + g_y_0_xxxzz_xxxxxz[k];

                g_y_0_xxxxzz_xxxyy[k] = -g_y_0_xxxzz_xxxyy[k] * cd_x[k] + g_y_0_xxxzz_xxxxyy[k];

                g_y_0_xxxxzz_xxxyz[k] = -g_y_0_xxxzz_xxxyz[k] * cd_x[k] + g_y_0_xxxzz_xxxxyz[k];

                g_y_0_xxxxzz_xxxzz[k] = -g_y_0_xxxzz_xxxzz[k] * cd_x[k] + g_y_0_xxxzz_xxxxzz[k];

                g_y_0_xxxxzz_xxyyy[k] = -g_y_0_xxxzz_xxyyy[k] * cd_x[k] + g_y_0_xxxzz_xxxyyy[k];

                g_y_0_xxxxzz_xxyyz[k] = -g_y_0_xxxzz_xxyyz[k] * cd_x[k] + g_y_0_xxxzz_xxxyyz[k];

                g_y_0_xxxxzz_xxyzz[k] = -g_y_0_xxxzz_xxyzz[k] * cd_x[k] + g_y_0_xxxzz_xxxyzz[k];

                g_y_0_xxxxzz_xxzzz[k] = -g_y_0_xxxzz_xxzzz[k] * cd_x[k] + g_y_0_xxxzz_xxxzzz[k];

                g_y_0_xxxxzz_xyyyy[k] = -g_y_0_xxxzz_xyyyy[k] * cd_x[k] + g_y_0_xxxzz_xxyyyy[k];

                g_y_0_xxxxzz_xyyyz[k] = -g_y_0_xxxzz_xyyyz[k] * cd_x[k] + g_y_0_xxxzz_xxyyyz[k];

                g_y_0_xxxxzz_xyyzz[k] = -g_y_0_xxxzz_xyyzz[k] * cd_x[k] + g_y_0_xxxzz_xxyyzz[k];

                g_y_0_xxxxzz_xyzzz[k] = -g_y_0_xxxzz_xyzzz[k] * cd_x[k] + g_y_0_xxxzz_xxyzzz[k];

                g_y_0_xxxxzz_xzzzz[k] = -g_y_0_xxxzz_xzzzz[k] * cd_x[k] + g_y_0_xxxzz_xxzzzz[k];

                g_y_0_xxxxzz_yyyyy[k] = -g_y_0_xxxzz_yyyyy[k] * cd_x[k] + g_y_0_xxxzz_xyyyyy[k];

                g_y_0_xxxxzz_yyyyz[k] = -g_y_0_xxxzz_yyyyz[k] * cd_x[k] + g_y_0_xxxzz_xyyyyz[k];

                g_y_0_xxxxzz_yyyzz[k] = -g_y_0_xxxzz_yyyzz[k] * cd_x[k] + g_y_0_xxxzz_xyyyzz[k];

                g_y_0_xxxxzz_yyzzz[k] = -g_y_0_xxxzz_yyzzz[k] * cd_x[k] + g_y_0_xxxzz_xyyzzz[k];

                g_y_0_xxxxzz_yzzzz[k] = -g_y_0_xxxzz_yzzzz[k] * cd_x[k] + g_y_0_xxxzz_xyzzzz[k];

                g_y_0_xxxxzz_zzzzz[k] = -g_y_0_xxxzz_zzzzz[k] * cd_x[k] + g_y_0_xxxzz_xzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 126);

            auto g_y_0_xxxyyy_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 127);

            auto g_y_0_xxxyyy_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 128);

            auto g_y_0_xxxyyy_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 129);

            auto g_y_0_xxxyyy_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 130);

            auto g_y_0_xxxyyy_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 131);

            auto g_y_0_xxxyyy_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 132);

            auto g_y_0_xxxyyy_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 133);

            auto g_y_0_xxxyyy_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 134);

            auto g_y_0_xxxyyy_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 135);

            auto g_y_0_xxxyyy_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 136);

            auto g_y_0_xxxyyy_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 137);

            auto g_y_0_xxxyyy_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 138);

            auto g_y_0_xxxyyy_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 139);

            auto g_y_0_xxxyyy_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 140);

            auto g_y_0_xxxyyy_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 141);

            auto g_y_0_xxxyyy_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 142);

            auto g_y_0_xxxyyy_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 143);

            auto g_y_0_xxxyyy_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 144);

            auto g_y_0_xxxyyy_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 145);

            auto g_y_0_xxxyyy_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyy_xxxxx, g_y_0_xxxyyy_xxxxy, g_y_0_xxxyyy_xxxxz, g_y_0_xxxyyy_xxxyy, g_y_0_xxxyyy_xxxyz, g_y_0_xxxyyy_xxxzz, g_y_0_xxxyyy_xxyyy, g_y_0_xxxyyy_xxyyz, g_y_0_xxxyyy_xxyzz, g_y_0_xxxyyy_xxzzz, g_y_0_xxxyyy_xyyyy, g_y_0_xxxyyy_xyyyz, g_y_0_xxxyyy_xyyzz, g_y_0_xxxyyy_xyzzz, g_y_0_xxxyyy_xzzzz, g_y_0_xxxyyy_yyyyy, g_y_0_xxxyyy_yyyyz, g_y_0_xxxyyy_yyyzz, g_y_0_xxxyyy_yyzzz, g_y_0_xxxyyy_yzzzz, g_y_0_xxxyyy_zzzzz, g_y_0_xxyyy_xxxxx, g_y_0_xxyyy_xxxxxx, g_y_0_xxyyy_xxxxxy, g_y_0_xxyyy_xxxxxz, g_y_0_xxyyy_xxxxy, g_y_0_xxyyy_xxxxyy, g_y_0_xxyyy_xxxxyz, g_y_0_xxyyy_xxxxz, g_y_0_xxyyy_xxxxzz, g_y_0_xxyyy_xxxyy, g_y_0_xxyyy_xxxyyy, g_y_0_xxyyy_xxxyyz, g_y_0_xxyyy_xxxyz, g_y_0_xxyyy_xxxyzz, g_y_0_xxyyy_xxxzz, g_y_0_xxyyy_xxxzzz, g_y_0_xxyyy_xxyyy, g_y_0_xxyyy_xxyyyy, g_y_0_xxyyy_xxyyyz, g_y_0_xxyyy_xxyyz, g_y_0_xxyyy_xxyyzz, g_y_0_xxyyy_xxyzz, g_y_0_xxyyy_xxyzzz, g_y_0_xxyyy_xxzzz, g_y_0_xxyyy_xxzzzz, g_y_0_xxyyy_xyyyy, g_y_0_xxyyy_xyyyyy, g_y_0_xxyyy_xyyyyz, g_y_0_xxyyy_xyyyz, g_y_0_xxyyy_xyyyzz, g_y_0_xxyyy_xyyzz, g_y_0_xxyyy_xyyzzz, g_y_0_xxyyy_xyzzz, g_y_0_xxyyy_xyzzzz, g_y_0_xxyyy_xzzzz, g_y_0_xxyyy_xzzzzz, g_y_0_xxyyy_yyyyy, g_y_0_xxyyy_yyyyz, g_y_0_xxyyy_yyyzz, g_y_0_xxyyy_yyzzz, g_y_0_xxyyy_yzzzz, g_y_0_xxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_xxxxx[k] = -g_y_0_xxyyy_xxxxx[k] * cd_x[k] + g_y_0_xxyyy_xxxxxx[k];

                g_y_0_xxxyyy_xxxxy[k] = -g_y_0_xxyyy_xxxxy[k] * cd_x[k] + g_y_0_xxyyy_xxxxxy[k];

                g_y_0_xxxyyy_xxxxz[k] = -g_y_0_xxyyy_xxxxz[k] * cd_x[k] + g_y_0_xxyyy_xxxxxz[k];

                g_y_0_xxxyyy_xxxyy[k] = -g_y_0_xxyyy_xxxyy[k] * cd_x[k] + g_y_0_xxyyy_xxxxyy[k];

                g_y_0_xxxyyy_xxxyz[k] = -g_y_0_xxyyy_xxxyz[k] * cd_x[k] + g_y_0_xxyyy_xxxxyz[k];

                g_y_0_xxxyyy_xxxzz[k] = -g_y_0_xxyyy_xxxzz[k] * cd_x[k] + g_y_0_xxyyy_xxxxzz[k];

                g_y_0_xxxyyy_xxyyy[k] = -g_y_0_xxyyy_xxyyy[k] * cd_x[k] + g_y_0_xxyyy_xxxyyy[k];

                g_y_0_xxxyyy_xxyyz[k] = -g_y_0_xxyyy_xxyyz[k] * cd_x[k] + g_y_0_xxyyy_xxxyyz[k];

                g_y_0_xxxyyy_xxyzz[k] = -g_y_0_xxyyy_xxyzz[k] * cd_x[k] + g_y_0_xxyyy_xxxyzz[k];

                g_y_0_xxxyyy_xxzzz[k] = -g_y_0_xxyyy_xxzzz[k] * cd_x[k] + g_y_0_xxyyy_xxxzzz[k];

                g_y_0_xxxyyy_xyyyy[k] = -g_y_0_xxyyy_xyyyy[k] * cd_x[k] + g_y_0_xxyyy_xxyyyy[k];

                g_y_0_xxxyyy_xyyyz[k] = -g_y_0_xxyyy_xyyyz[k] * cd_x[k] + g_y_0_xxyyy_xxyyyz[k];

                g_y_0_xxxyyy_xyyzz[k] = -g_y_0_xxyyy_xyyzz[k] * cd_x[k] + g_y_0_xxyyy_xxyyzz[k];

                g_y_0_xxxyyy_xyzzz[k] = -g_y_0_xxyyy_xyzzz[k] * cd_x[k] + g_y_0_xxyyy_xxyzzz[k];

                g_y_0_xxxyyy_xzzzz[k] = -g_y_0_xxyyy_xzzzz[k] * cd_x[k] + g_y_0_xxyyy_xxzzzz[k];

                g_y_0_xxxyyy_yyyyy[k] = -g_y_0_xxyyy_yyyyy[k] * cd_x[k] + g_y_0_xxyyy_xyyyyy[k];

                g_y_0_xxxyyy_yyyyz[k] = -g_y_0_xxyyy_yyyyz[k] * cd_x[k] + g_y_0_xxyyy_xyyyyz[k];

                g_y_0_xxxyyy_yyyzz[k] = -g_y_0_xxyyy_yyyzz[k] * cd_x[k] + g_y_0_xxyyy_xyyyzz[k];

                g_y_0_xxxyyy_yyzzz[k] = -g_y_0_xxyyy_yyzzz[k] * cd_x[k] + g_y_0_xxyyy_xyyzzz[k];

                g_y_0_xxxyyy_yzzzz[k] = -g_y_0_xxyyy_yzzzz[k] * cd_x[k] + g_y_0_xxyyy_xyzzzz[k];

                g_y_0_xxxyyy_zzzzz[k] = -g_y_0_xxyyy_zzzzz[k] * cd_x[k] + g_y_0_xxyyy_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 147);

            auto g_y_0_xxxyyz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 148);

            auto g_y_0_xxxyyz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 149);

            auto g_y_0_xxxyyz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 150);

            auto g_y_0_xxxyyz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 151);

            auto g_y_0_xxxyyz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 152);

            auto g_y_0_xxxyyz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 153);

            auto g_y_0_xxxyyz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 154);

            auto g_y_0_xxxyyz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 155);

            auto g_y_0_xxxyyz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 156);

            auto g_y_0_xxxyyz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 157);

            auto g_y_0_xxxyyz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 158);

            auto g_y_0_xxxyyz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 159);

            auto g_y_0_xxxyyz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 160);

            auto g_y_0_xxxyyz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 161);

            auto g_y_0_xxxyyz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 162);

            auto g_y_0_xxxyyz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 163);

            auto g_y_0_xxxyyz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 164);

            auto g_y_0_xxxyyz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 165);

            auto g_y_0_xxxyyz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 166);

            auto g_y_0_xxxyyz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyz_xxxxx, g_y_0_xxxyyz_xxxxy, g_y_0_xxxyyz_xxxxz, g_y_0_xxxyyz_xxxyy, g_y_0_xxxyyz_xxxyz, g_y_0_xxxyyz_xxxzz, g_y_0_xxxyyz_xxyyy, g_y_0_xxxyyz_xxyyz, g_y_0_xxxyyz_xxyzz, g_y_0_xxxyyz_xxzzz, g_y_0_xxxyyz_xyyyy, g_y_0_xxxyyz_xyyyz, g_y_0_xxxyyz_xyyzz, g_y_0_xxxyyz_xyzzz, g_y_0_xxxyyz_xzzzz, g_y_0_xxxyyz_yyyyy, g_y_0_xxxyyz_yyyyz, g_y_0_xxxyyz_yyyzz, g_y_0_xxxyyz_yyzzz, g_y_0_xxxyyz_yzzzz, g_y_0_xxxyyz_zzzzz, g_y_0_xxyyz_xxxxx, g_y_0_xxyyz_xxxxxx, g_y_0_xxyyz_xxxxxy, g_y_0_xxyyz_xxxxxz, g_y_0_xxyyz_xxxxy, g_y_0_xxyyz_xxxxyy, g_y_0_xxyyz_xxxxyz, g_y_0_xxyyz_xxxxz, g_y_0_xxyyz_xxxxzz, g_y_0_xxyyz_xxxyy, g_y_0_xxyyz_xxxyyy, g_y_0_xxyyz_xxxyyz, g_y_0_xxyyz_xxxyz, g_y_0_xxyyz_xxxyzz, g_y_0_xxyyz_xxxzz, g_y_0_xxyyz_xxxzzz, g_y_0_xxyyz_xxyyy, g_y_0_xxyyz_xxyyyy, g_y_0_xxyyz_xxyyyz, g_y_0_xxyyz_xxyyz, g_y_0_xxyyz_xxyyzz, g_y_0_xxyyz_xxyzz, g_y_0_xxyyz_xxyzzz, g_y_0_xxyyz_xxzzz, g_y_0_xxyyz_xxzzzz, g_y_0_xxyyz_xyyyy, g_y_0_xxyyz_xyyyyy, g_y_0_xxyyz_xyyyyz, g_y_0_xxyyz_xyyyz, g_y_0_xxyyz_xyyyzz, g_y_0_xxyyz_xyyzz, g_y_0_xxyyz_xyyzzz, g_y_0_xxyyz_xyzzz, g_y_0_xxyyz_xyzzzz, g_y_0_xxyyz_xzzzz, g_y_0_xxyyz_xzzzzz, g_y_0_xxyyz_yyyyy, g_y_0_xxyyz_yyyyz, g_y_0_xxyyz_yyyzz, g_y_0_xxyyz_yyzzz, g_y_0_xxyyz_yzzzz, g_y_0_xxyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_xxxxx[k] = -g_y_0_xxyyz_xxxxx[k] * cd_x[k] + g_y_0_xxyyz_xxxxxx[k];

                g_y_0_xxxyyz_xxxxy[k] = -g_y_0_xxyyz_xxxxy[k] * cd_x[k] + g_y_0_xxyyz_xxxxxy[k];

                g_y_0_xxxyyz_xxxxz[k] = -g_y_0_xxyyz_xxxxz[k] * cd_x[k] + g_y_0_xxyyz_xxxxxz[k];

                g_y_0_xxxyyz_xxxyy[k] = -g_y_0_xxyyz_xxxyy[k] * cd_x[k] + g_y_0_xxyyz_xxxxyy[k];

                g_y_0_xxxyyz_xxxyz[k] = -g_y_0_xxyyz_xxxyz[k] * cd_x[k] + g_y_0_xxyyz_xxxxyz[k];

                g_y_0_xxxyyz_xxxzz[k] = -g_y_0_xxyyz_xxxzz[k] * cd_x[k] + g_y_0_xxyyz_xxxxzz[k];

                g_y_0_xxxyyz_xxyyy[k] = -g_y_0_xxyyz_xxyyy[k] * cd_x[k] + g_y_0_xxyyz_xxxyyy[k];

                g_y_0_xxxyyz_xxyyz[k] = -g_y_0_xxyyz_xxyyz[k] * cd_x[k] + g_y_0_xxyyz_xxxyyz[k];

                g_y_0_xxxyyz_xxyzz[k] = -g_y_0_xxyyz_xxyzz[k] * cd_x[k] + g_y_0_xxyyz_xxxyzz[k];

                g_y_0_xxxyyz_xxzzz[k] = -g_y_0_xxyyz_xxzzz[k] * cd_x[k] + g_y_0_xxyyz_xxxzzz[k];

                g_y_0_xxxyyz_xyyyy[k] = -g_y_0_xxyyz_xyyyy[k] * cd_x[k] + g_y_0_xxyyz_xxyyyy[k];

                g_y_0_xxxyyz_xyyyz[k] = -g_y_0_xxyyz_xyyyz[k] * cd_x[k] + g_y_0_xxyyz_xxyyyz[k];

                g_y_0_xxxyyz_xyyzz[k] = -g_y_0_xxyyz_xyyzz[k] * cd_x[k] + g_y_0_xxyyz_xxyyzz[k];

                g_y_0_xxxyyz_xyzzz[k] = -g_y_0_xxyyz_xyzzz[k] * cd_x[k] + g_y_0_xxyyz_xxyzzz[k];

                g_y_0_xxxyyz_xzzzz[k] = -g_y_0_xxyyz_xzzzz[k] * cd_x[k] + g_y_0_xxyyz_xxzzzz[k];

                g_y_0_xxxyyz_yyyyy[k] = -g_y_0_xxyyz_yyyyy[k] * cd_x[k] + g_y_0_xxyyz_xyyyyy[k];

                g_y_0_xxxyyz_yyyyz[k] = -g_y_0_xxyyz_yyyyz[k] * cd_x[k] + g_y_0_xxyyz_xyyyyz[k];

                g_y_0_xxxyyz_yyyzz[k] = -g_y_0_xxyyz_yyyzz[k] * cd_x[k] + g_y_0_xxyyz_xyyyzz[k];

                g_y_0_xxxyyz_yyzzz[k] = -g_y_0_xxyyz_yyzzz[k] * cd_x[k] + g_y_0_xxyyz_xyyzzz[k];

                g_y_0_xxxyyz_yzzzz[k] = -g_y_0_xxyyz_yzzzz[k] * cd_x[k] + g_y_0_xxyyz_xyzzzz[k];

                g_y_0_xxxyyz_zzzzz[k] = -g_y_0_xxyyz_zzzzz[k] * cd_x[k] + g_y_0_xxyyz_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 168);

            auto g_y_0_xxxyzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 169);

            auto g_y_0_xxxyzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 170);

            auto g_y_0_xxxyzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 171);

            auto g_y_0_xxxyzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 172);

            auto g_y_0_xxxyzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 173);

            auto g_y_0_xxxyzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 174);

            auto g_y_0_xxxyzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 175);

            auto g_y_0_xxxyzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 176);

            auto g_y_0_xxxyzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 177);

            auto g_y_0_xxxyzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 178);

            auto g_y_0_xxxyzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 179);

            auto g_y_0_xxxyzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 180);

            auto g_y_0_xxxyzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 181);

            auto g_y_0_xxxyzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 182);

            auto g_y_0_xxxyzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 183);

            auto g_y_0_xxxyzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 184);

            auto g_y_0_xxxyzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 185);

            auto g_y_0_xxxyzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 186);

            auto g_y_0_xxxyzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 187);

            auto g_y_0_xxxyzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyzz_xxxxx, g_y_0_xxxyzz_xxxxy, g_y_0_xxxyzz_xxxxz, g_y_0_xxxyzz_xxxyy, g_y_0_xxxyzz_xxxyz, g_y_0_xxxyzz_xxxzz, g_y_0_xxxyzz_xxyyy, g_y_0_xxxyzz_xxyyz, g_y_0_xxxyzz_xxyzz, g_y_0_xxxyzz_xxzzz, g_y_0_xxxyzz_xyyyy, g_y_0_xxxyzz_xyyyz, g_y_0_xxxyzz_xyyzz, g_y_0_xxxyzz_xyzzz, g_y_0_xxxyzz_xzzzz, g_y_0_xxxyzz_yyyyy, g_y_0_xxxyzz_yyyyz, g_y_0_xxxyzz_yyyzz, g_y_0_xxxyzz_yyzzz, g_y_0_xxxyzz_yzzzz, g_y_0_xxxyzz_zzzzz, g_y_0_xxyzz_xxxxx, g_y_0_xxyzz_xxxxxx, g_y_0_xxyzz_xxxxxy, g_y_0_xxyzz_xxxxxz, g_y_0_xxyzz_xxxxy, g_y_0_xxyzz_xxxxyy, g_y_0_xxyzz_xxxxyz, g_y_0_xxyzz_xxxxz, g_y_0_xxyzz_xxxxzz, g_y_0_xxyzz_xxxyy, g_y_0_xxyzz_xxxyyy, g_y_0_xxyzz_xxxyyz, g_y_0_xxyzz_xxxyz, g_y_0_xxyzz_xxxyzz, g_y_0_xxyzz_xxxzz, g_y_0_xxyzz_xxxzzz, g_y_0_xxyzz_xxyyy, g_y_0_xxyzz_xxyyyy, g_y_0_xxyzz_xxyyyz, g_y_0_xxyzz_xxyyz, g_y_0_xxyzz_xxyyzz, g_y_0_xxyzz_xxyzz, g_y_0_xxyzz_xxyzzz, g_y_0_xxyzz_xxzzz, g_y_0_xxyzz_xxzzzz, g_y_0_xxyzz_xyyyy, g_y_0_xxyzz_xyyyyy, g_y_0_xxyzz_xyyyyz, g_y_0_xxyzz_xyyyz, g_y_0_xxyzz_xyyyzz, g_y_0_xxyzz_xyyzz, g_y_0_xxyzz_xyyzzz, g_y_0_xxyzz_xyzzz, g_y_0_xxyzz_xyzzzz, g_y_0_xxyzz_xzzzz, g_y_0_xxyzz_xzzzzz, g_y_0_xxyzz_yyyyy, g_y_0_xxyzz_yyyyz, g_y_0_xxyzz_yyyzz, g_y_0_xxyzz_yyzzz, g_y_0_xxyzz_yzzzz, g_y_0_xxyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_xxxxx[k] = -g_y_0_xxyzz_xxxxx[k] * cd_x[k] + g_y_0_xxyzz_xxxxxx[k];

                g_y_0_xxxyzz_xxxxy[k] = -g_y_0_xxyzz_xxxxy[k] * cd_x[k] + g_y_0_xxyzz_xxxxxy[k];

                g_y_0_xxxyzz_xxxxz[k] = -g_y_0_xxyzz_xxxxz[k] * cd_x[k] + g_y_0_xxyzz_xxxxxz[k];

                g_y_0_xxxyzz_xxxyy[k] = -g_y_0_xxyzz_xxxyy[k] * cd_x[k] + g_y_0_xxyzz_xxxxyy[k];

                g_y_0_xxxyzz_xxxyz[k] = -g_y_0_xxyzz_xxxyz[k] * cd_x[k] + g_y_0_xxyzz_xxxxyz[k];

                g_y_0_xxxyzz_xxxzz[k] = -g_y_0_xxyzz_xxxzz[k] * cd_x[k] + g_y_0_xxyzz_xxxxzz[k];

                g_y_0_xxxyzz_xxyyy[k] = -g_y_0_xxyzz_xxyyy[k] * cd_x[k] + g_y_0_xxyzz_xxxyyy[k];

                g_y_0_xxxyzz_xxyyz[k] = -g_y_0_xxyzz_xxyyz[k] * cd_x[k] + g_y_0_xxyzz_xxxyyz[k];

                g_y_0_xxxyzz_xxyzz[k] = -g_y_0_xxyzz_xxyzz[k] * cd_x[k] + g_y_0_xxyzz_xxxyzz[k];

                g_y_0_xxxyzz_xxzzz[k] = -g_y_0_xxyzz_xxzzz[k] * cd_x[k] + g_y_0_xxyzz_xxxzzz[k];

                g_y_0_xxxyzz_xyyyy[k] = -g_y_0_xxyzz_xyyyy[k] * cd_x[k] + g_y_0_xxyzz_xxyyyy[k];

                g_y_0_xxxyzz_xyyyz[k] = -g_y_0_xxyzz_xyyyz[k] * cd_x[k] + g_y_0_xxyzz_xxyyyz[k];

                g_y_0_xxxyzz_xyyzz[k] = -g_y_0_xxyzz_xyyzz[k] * cd_x[k] + g_y_0_xxyzz_xxyyzz[k];

                g_y_0_xxxyzz_xyzzz[k] = -g_y_0_xxyzz_xyzzz[k] * cd_x[k] + g_y_0_xxyzz_xxyzzz[k];

                g_y_0_xxxyzz_xzzzz[k] = -g_y_0_xxyzz_xzzzz[k] * cd_x[k] + g_y_0_xxyzz_xxzzzz[k];

                g_y_0_xxxyzz_yyyyy[k] = -g_y_0_xxyzz_yyyyy[k] * cd_x[k] + g_y_0_xxyzz_xyyyyy[k];

                g_y_0_xxxyzz_yyyyz[k] = -g_y_0_xxyzz_yyyyz[k] * cd_x[k] + g_y_0_xxyzz_xyyyyz[k];

                g_y_0_xxxyzz_yyyzz[k] = -g_y_0_xxyzz_yyyzz[k] * cd_x[k] + g_y_0_xxyzz_xyyyzz[k];

                g_y_0_xxxyzz_yyzzz[k] = -g_y_0_xxyzz_yyzzz[k] * cd_x[k] + g_y_0_xxyzz_xyyzzz[k];

                g_y_0_xxxyzz_yzzzz[k] = -g_y_0_xxyzz_yzzzz[k] * cd_x[k] + g_y_0_xxyzz_xyzzzz[k];

                g_y_0_xxxyzz_zzzzz[k] = -g_y_0_xxyzz_zzzzz[k] * cd_x[k] + g_y_0_xxyzz_xzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 189);

            auto g_y_0_xxxzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 190);

            auto g_y_0_xxxzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 191);

            auto g_y_0_xxxzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 192);

            auto g_y_0_xxxzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 193);

            auto g_y_0_xxxzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 194);

            auto g_y_0_xxxzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 195);

            auto g_y_0_xxxzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 196);

            auto g_y_0_xxxzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 197);

            auto g_y_0_xxxzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 198);

            auto g_y_0_xxxzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 199);

            auto g_y_0_xxxzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 200);

            auto g_y_0_xxxzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 201);

            auto g_y_0_xxxzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 202);

            auto g_y_0_xxxzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 203);

            auto g_y_0_xxxzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 204);

            auto g_y_0_xxxzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 205);

            auto g_y_0_xxxzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 206);

            auto g_y_0_xxxzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 207);

            auto g_y_0_xxxzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 208);

            auto g_y_0_xxxzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzzz_xxxxx, g_y_0_xxxzzz_xxxxy, g_y_0_xxxzzz_xxxxz, g_y_0_xxxzzz_xxxyy, g_y_0_xxxzzz_xxxyz, g_y_0_xxxzzz_xxxzz, g_y_0_xxxzzz_xxyyy, g_y_0_xxxzzz_xxyyz, g_y_0_xxxzzz_xxyzz, g_y_0_xxxzzz_xxzzz, g_y_0_xxxzzz_xyyyy, g_y_0_xxxzzz_xyyyz, g_y_0_xxxzzz_xyyzz, g_y_0_xxxzzz_xyzzz, g_y_0_xxxzzz_xzzzz, g_y_0_xxxzzz_yyyyy, g_y_0_xxxzzz_yyyyz, g_y_0_xxxzzz_yyyzz, g_y_0_xxxzzz_yyzzz, g_y_0_xxxzzz_yzzzz, g_y_0_xxxzzz_zzzzz, g_y_0_xxzzz_xxxxx, g_y_0_xxzzz_xxxxxx, g_y_0_xxzzz_xxxxxy, g_y_0_xxzzz_xxxxxz, g_y_0_xxzzz_xxxxy, g_y_0_xxzzz_xxxxyy, g_y_0_xxzzz_xxxxyz, g_y_0_xxzzz_xxxxz, g_y_0_xxzzz_xxxxzz, g_y_0_xxzzz_xxxyy, g_y_0_xxzzz_xxxyyy, g_y_0_xxzzz_xxxyyz, g_y_0_xxzzz_xxxyz, g_y_0_xxzzz_xxxyzz, g_y_0_xxzzz_xxxzz, g_y_0_xxzzz_xxxzzz, g_y_0_xxzzz_xxyyy, g_y_0_xxzzz_xxyyyy, g_y_0_xxzzz_xxyyyz, g_y_0_xxzzz_xxyyz, g_y_0_xxzzz_xxyyzz, g_y_0_xxzzz_xxyzz, g_y_0_xxzzz_xxyzzz, g_y_0_xxzzz_xxzzz, g_y_0_xxzzz_xxzzzz, g_y_0_xxzzz_xyyyy, g_y_0_xxzzz_xyyyyy, g_y_0_xxzzz_xyyyyz, g_y_0_xxzzz_xyyyz, g_y_0_xxzzz_xyyyzz, g_y_0_xxzzz_xyyzz, g_y_0_xxzzz_xyyzzz, g_y_0_xxzzz_xyzzz, g_y_0_xxzzz_xyzzzz, g_y_0_xxzzz_xzzzz, g_y_0_xxzzz_xzzzzz, g_y_0_xxzzz_yyyyy, g_y_0_xxzzz_yyyyz, g_y_0_xxzzz_yyyzz, g_y_0_xxzzz_yyzzz, g_y_0_xxzzz_yzzzz, g_y_0_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_xxxxx[k] = -g_y_0_xxzzz_xxxxx[k] * cd_x[k] + g_y_0_xxzzz_xxxxxx[k];

                g_y_0_xxxzzz_xxxxy[k] = -g_y_0_xxzzz_xxxxy[k] * cd_x[k] + g_y_0_xxzzz_xxxxxy[k];

                g_y_0_xxxzzz_xxxxz[k] = -g_y_0_xxzzz_xxxxz[k] * cd_x[k] + g_y_0_xxzzz_xxxxxz[k];

                g_y_0_xxxzzz_xxxyy[k] = -g_y_0_xxzzz_xxxyy[k] * cd_x[k] + g_y_0_xxzzz_xxxxyy[k];

                g_y_0_xxxzzz_xxxyz[k] = -g_y_0_xxzzz_xxxyz[k] * cd_x[k] + g_y_0_xxzzz_xxxxyz[k];

                g_y_0_xxxzzz_xxxzz[k] = -g_y_0_xxzzz_xxxzz[k] * cd_x[k] + g_y_0_xxzzz_xxxxzz[k];

                g_y_0_xxxzzz_xxyyy[k] = -g_y_0_xxzzz_xxyyy[k] * cd_x[k] + g_y_0_xxzzz_xxxyyy[k];

                g_y_0_xxxzzz_xxyyz[k] = -g_y_0_xxzzz_xxyyz[k] * cd_x[k] + g_y_0_xxzzz_xxxyyz[k];

                g_y_0_xxxzzz_xxyzz[k] = -g_y_0_xxzzz_xxyzz[k] * cd_x[k] + g_y_0_xxzzz_xxxyzz[k];

                g_y_0_xxxzzz_xxzzz[k] = -g_y_0_xxzzz_xxzzz[k] * cd_x[k] + g_y_0_xxzzz_xxxzzz[k];

                g_y_0_xxxzzz_xyyyy[k] = -g_y_0_xxzzz_xyyyy[k] * cd_x[k] + g_y_0_xxzzz_xxyyyy[k];

                g_y_0_xxxzzz_xyyyz[k] = -g_y_0_xxzzz_xyyyz[k] * cd_x[k] + g_y_0_xxzzz_xxyyyz[k];

                g_y_0_xxxzzz_xyyzz[k] = -g_y_0_xxzzz_xyyzz[k] * cd_x[k] + g_y_0_xxzzz_xxyyzz[k];

                g_y_0_xxxzzz_xyzzz[k] = -g_y_0_xxzzz_xyzzz[k] * cd_x[k] + g_y_0_xxzzz_xxyzzz[k];

                g_y_0_xxxzzz_xzzzz[k] = -g_y_0_xxzzz_xzzzz[k] * cd_x[k] + g_y_0_xxzzz_xxzzzz[k];

                g_y_0_xxxzzz_yyyyy[k] = -g_y_0_xxzzz_yyyyy[k] * cd_x[k] + g_y_0_xxzzz_xyyyyy[k];

                g_y_0_xxxzzz_yyyyz[k] = -g_y_0_xxzzz_yyyyz[k] * cd_x[k] + g_y_0_xxzzz_xyyyyz[k];

                g_y_0_xxxzzz_yyyzz[k] = -g_y_0_xxzzz_yyyzz[k] * cd_x[k] + g_y_0_xxzzz_xyyyzz[k];

                g_y_0_xxxzzz_yyzzz[k] = -g_y_0_xxzzz_yyzzz[k] * cd_x[k] + g_y_0_xxzzz_xyyzzz[k];

                g_y_0_xxxzzz_yzzzz[k] = -g_y_0_xxzzz_yzzzz[k] * cd_x[k] + g_y_0_xxzzz_xyzzzz[k];

                g_y_0_xxxzzz_zzzzz[k] = -g_y_0_xxzzz_zzzzz[k] * cd_x[k] + g_y_0_xxzzz_xzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 210);

            auto g_y_0_xxyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 211);

            auto g_y_0_xxyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 212);

            auto g_y_0_xxyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 213);

            auto g_y_0_xxyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 214);

            auto g_y_0_xxyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 215);

            auto g_y_0_xxyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 216);

            auto g_y_0_xxyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 217);

            auto g_y_0_xxyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 218);

            auto g_y_0_xxyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 219);

            auto g_y_0_xxyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 220);

            auto g_y_0_xxyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 221);

            auto g_y_0_xxyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 222);

            auto g_y_0_xxyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 223);

            auto g_y_0_xxyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 224);

            auto g_y_0_xxyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 225);

            auto g_y_0_xxyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 226);

            auto g_y_0_xxyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 227);

            auto g_y_0_xxyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 228);

            auto g_y_0_xxyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 229);

            auto g_y_0_xxyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyy_xxxxx, g_y_0_xxyyyy_xxxxy, g_y_0_xxyyyy_xxxxz, g_y_0_xxyyyy_xxxyy, g_y_0_xxyyyy_xxxyz, g_y_0_xxyyyy_xxxzz, g_y_0_xxyyyy_xxyyy, g_y_0_xxyyyy_xxyyz, g_y_0_xxyyyy_xxyzz, g_y_0_xxyyyy_xxzzz, g_y_0_xxyyyy_xyyyy, g_y_0_xxyyyy_xyyyz, g_y_0_xxyyyy_xyyzz, g_y_0_xxyyyy_xyzzz, g_y_0_xxyyyy_xzzzz, g_y_0_xxyyyy_yyyyy, g_y_0_xxyyyy_yyyyz, g_y_0_xxyyyy_yyyzz, g_y_0_xxyyyy_yyzzz, g_y_0_xxyyyy_yzzzz, g_y_0_xxyyyy_zzzzz, g_y_0_xyyyy_xxxxx, g_y_0_xyyyy_xxxxxx, g_y_0_xyyyy_xxxxxy, g_y_0_xyyyy_xxxxxz, g_y_0_xyyyy_xxxxy, g_y_0_xyyyy_xxxxyy, g_y_0_xyyyy_xxxxyz, g_y_0_xyyyy_xxxxz, g_y_0_xyyyy_xxxxzz, g_y_0_xyyyy_xxxyy, g_y_0_xyyyy_xxxyyy, g_y_0_xyyyy_xxxyyz, g_y_0_xyyyy_xxxyz, g_y_0_xyyyy_xxxyzz, g_y_0_xyyyy_xxxzz, g_y_0_xyyyy_xxxzzz, g_y_0_xyyyy_xxyyy, g_y_0_xyyyy_xxyyyy, g_y_0_xyyyy_xxyyyz, g_y_0_xyyyy_xxyyz, g_y_0_xyyyy_xxyyzz, g_y_0_xyyyy_xxyzz, g_y_0_xyyyy_xxyzzz, g_y_0_xyyyy_xxzzz, g_y_0_xyyyy_xxzzzz, g_y_0_xyyyy_xyyyy, g_y_0_xyyyy_xyyyyy, g_y_0_xyyyy_xyyyyz, g_y_0_xyyyy_xyyyz, g_y_0_xyyyy_xyyyzz, g_y_0_xyyyy_xyyzz, g_y_0_xyyyy_xyyzzz, g_y_0_xyyyy_xyzzz, g_y_0_xyyyy_xyzzzz, g_y_0_xyyyy_xzzzz, g_y_0_xyyyy_xzzzzz, g_y_0_xyyyy_yyyyy, g_y_0_xyyyy_yyyyz, g_y_0_xyyyy_yyyzz, g_y_0_xyyyy_yyzzz, g_y_0_xyyyy_yzzzz, g_y_0_xyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_xxxxx[k] = -g_y_0_xyyyy_xxxxx[k] * cd_x[k] + g_y_0_xyyyy_xxxxxx[k];

                g_y_0_xxyyyy_xxxxy[k] = -g_y_0_xyyyy_xxxxy[k] * cd_x[k] + g_y_0_xyyyy_xxxxxy[k];

                g_y_0_xxyyyy_xxxxz[k] = -g_y_0_xyyyy_xxxxz[k] * cd_x[k] + g_y_0_xyyyy_xxxxxz[k];

                g_y_0_xxyyyy_xxxyy[k] = -g_y_0_xyyyy_xxxyy[k] * cd_x[k] + g_y_0_xyyyy_xxxxyy[k];

                g_y_0_xxyyyy_xxxyz[k] = -g_y_0_xyyyy_xxxyz[k] * cd_x[k] + g_y_0_xyyyy_xxxxyz[k];

                g_y_0_xxyyyy_xxxzz[k] = -g_y_0_xyyyy_xxxzz[k] * cd_x[k] + g_y_0_xyyyy_xxxxzz[k];

                g_y_0_xxyyyy_xxyyy[k] = -g_y_0_xyyyy_xxyyy[k] * cd_x[k] + g_y_0_xyyyy_xxxyyy[k];

                g_y_0_xxyyyy_xxyyz[k] = -g_y_0_xyyyy_xxyyz[k] * cd_x[k] + g_y_0_xyyyy_xxxyyz[k];

                g_y_0_xxyyyy_xxyzz[k] = -g_y_0_xyyyy_xxyzz[k] * cd_x[k] + g_y_0_xyyyy_xxxyzz[k];

                g_y_0_xxyyyy_xxzzz[k] = -g_y_0_xyyyy_xxzzz[k] * cd_x[k] + g_y_0_xyyyy_xxxzzz[k];

                g_y_0_xxyyyy_xyyyy[k] = -g_y_0_xyyyy_xyyyy[k] * cd_x[k] + g_y_0_xyyyy_xxyyyy[k];

                g_y_0_xxyyyy_xyyyz[k] = -g_y_0_xyyyy_xyyyz[k] * cd_x[k] + g_y_0_xyyyy_xxyyyz[k];

                g_y_0_xxyyyy_xyyzz[k] = -g_y_0_xyyyy_xyyzz[k] * cd_x[k] + g_y_0_xyyyy_xxyyzz[k];

                g_y_0_xxyyyy_xyzzz[k] = -g_y_0_xyyyy_xyzzz[k] * cd_x[k] + g_y_0_xyyyy_xxyzzz[k];

                g_y_0_xxyyyy_xzzzz[k] = -g_y_0_xyyyy_xzzzz[k] * cd_x[k] + g_y_0_xyyyy_xxzzzz[k];

                g_y_0_xxyyyy_yyyyy[k] = -g_y_0_xyyyy_yyyyy[k] * cd_x[k] + g_y_0_xyyyy_xyyyyy[k];

                g_y_0_xxyyyy_yyyyz[k] = -g_y_0_xyyyy_yyyyz[k] * cd_x[k] + g_y_0_xyyyy_xyyyyz[k];

                g_y_0_xxyyyy_yyyzz[k] = -g_y_0_xyyyy_yyyzz[k] * cd_x[k] + g_y_0_xyyyy_xyyyzz[k];

                g_y_0_xxyyyy_yyzzz[k] = -g_y_0_xyyyy_yyzzz[k] * cd_x[k] + g_y_0_xyyyy_xyyzzz[k];

                g_y_0_xxyyyy_yzzzz[k] = -g_y_0_xyyyy_yzzzz[k] * cd_x[k] + g_y_0_xyyyy_xyzzzz[k];

                g_y_0_xxyyyy_zzzzz[k] = -g_y_0_xyyyy_zzzzz[k] * cd_x[k] + g_y_0_xyyyy_xzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 231);

            auto g_y_0_xxyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 232);

            auto g_y_0_xxyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 233);

            auto g_y_0_xxyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 234);

            auto g_y_0_xxyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 235);

            auto g_y_0_xxyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 236);

            auto g_y_0_xxyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 237);

            auto g_y_0_xxyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 238);

            auto g_y_0_xxyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 239);

            auto g_y_0_xxyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 240);

            auto g_y_0_xxyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 241);

            auto g_y_0_xxyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 242);

            auto g_y_0_xxyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 243);

            auto g_y_0_xxyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 244);

            auto g_y_0_xxyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 245);

            auto g_y_0_xxyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 246);

            auto g_y_0_xxyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 247);

            auto g_y_0_xxyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 248);

            auto g_y_0_xxyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 249);

            auto g_y_0_xxyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 250);

            auto g_y_0_xxyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyz_xxxxx, g_y_0_xxyyyz_xxxxy, g_y_0_xxyyyz_xxxxz, g_y_0_xxyyyz_xxxyy, g_y_0_xxyyyz_xxxyz, g_y_0_xxyyyz_xxxzz, g_y_0_xxyyyz_xxyyy, g_y_0_xxyyyz_xxyyz, g_y_0_xxyyyz_xxyzz, g_y_0_xxyyyz_xxzzz, g_y_0_xxyyyz_xyyyy, g_y_0_xxyyyz_xyyyz, g_y_0_xxyyyz_xyyzz, g_y_0_xxyyyz_xyzzz, g_y_0_xxyyyz_xzzzz, g_y_0_xxyyyz_yyyyy, g_y_0_xxyyyz_yyyyz, g_y_0_xxyyyz_yyyzz, g_y_0_xxyyyz_yyzzz, g_y_0_xxyyyz_yzzzz, g_y_0_xxyyyz_zzzzz, g_y_0_xyyyz_xxxxx, g_y_0_xyyyz_xxxxxx, g_y_0_xyyyz_xxxxxy, g_y_0_xyyyz_xxxxxz, g_y_0_xyyyz_xxxxy, g_y_0_xyyyz_xxxxyy, g_y_0_xyyyz_xxxxyz, g_y_0_xyyyz_xxxxz, g_y_0_xyyyz_xxxxzz, g_y_0_xyyyz_xxxyy, g_y_0_xyyyz_xxxyyy, g_y_0_xyyyz_xxxyyz, g_y_0_xyyyz_xxxyz, g_y_0_xyyyz_xxxyzz, g_y_0_xyyyz_xxxzz, g_y_0_xyyyz_xxxzzz, g_y_0_xyyyz_xxyyy, g_y_0_xyyyz_xxyyyy, g_y_0_xyyyz_xxyyyz, g_y_0_xyyyz_xxyyz, g_y_0_xyyyz_xxyyzz, g_y_0_xyyyz_xxyzz, g_y_0_xyyyz_xxyzzz, g_y_0_xyyyz_xxzzz, g_y_0_xyyyz_xxzzzz, g_y_0_xyyyz_xyyyy, g_y_0_xyyyz_xyyyyy, g_y_0_xyyyz_xyyyyz, g_y_0_xyyyz_xyyyz, g_y_0_xyyyz_xyyyzz, g_y_0_xyyyz_xyyzz, g_y_0_xyyyz_xyyzzz, g_y_0_xyyyz_xyzzz, g_y_0_xyyyz_xyzzzz, g_y_0_xyyyz_xzzzz, g_y_0_xyyyz_xzzzzz, g_y_0_xyyyz_yyyyy, g_y_0_xyyyz_yyyyz, g_y_0_xyyyz_yyyzz, g_y_0_xyyyz_yyzzz, g_y_0_xyyyz_yzzzz, g_y_0_xyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_xxxxx[k] = -g_y_0_xyyyz_xxxxx[k] * cd_x[k] + g_y_0_xyyyz_xxxxxx[k];

                g_y_0_xxyyyz_xxxxy[k] = -g_y_0_xyyyz_xxxxy[k] * cd_x[k] + g_y_0_xyyyz_xxxxxy[k];

                g_y_0_xxyyyz_xxxxz[k] = -g_y_0_xyyyz_xxxxz[k] * cd_x[k] + g_y_0_xyyyz_xxxxxz[k];

                g_y_0_xxyyyz_xxxyy[k] = -g_y_0_xyyyz_xxxyy[k] * cd_x[k] + g_y_0_xyyyz_xxxxyy[k];

                g_y_0_xxyyyz_xxxyz[k] = -g_y_0_xyyyz_xxxyz[k] * cd_x[k] + g_y_0_xyyyz_xxxxyz[k];

                g_y_0_xxyyyz_xxxzz[k] = -g_y_0_xyyyz_xxxzz[k] * cd_x[k] + g_y_0_xyyyz_xxxxzz[k];

                g_y_0_xxyyyz_xxyyy[k] = -g_y_0_xyyyz_xxyyy[k] * cd_x[k] + g_y_0_xyyyz_xxxyyy[k];

                g_y_0_xxyyyz_xxyyz[k] = -g_y_0_xyyyz_xxyyz[k] * cd_x[k] + g_y_0_xyyyz_xxxyyz[k];

                g_y_0_xxyyyz_xxyzz[k] = -g_y_0_xyyyz_xxyzz[k] * cd_x[k] + g_y_0_xyyyz_xxxyzz[k];

                g_y_0_xxyyyz_xxzzz[k] = -g_y_0_xyyyz_xxzzz[k] * cd_x[k] + g_y_0_xyyyz_xxxzzz[k];

                g_y_0_xxyyyz_xyyyy[k] = -g_y_0_xyyyz_xyyyy[k] * cd_x[k] + g_y_0_xyyyz_xxyyyy[k];

                g_y_0_xxyyyz_xyyyz[k] = -g_y_0_xyyyz_xyyyz[k] * cd_x[k] + g_y_0_xyyyz_xxyyyz[k];

                g_y_0_xxyyyz_xyyzz[k] = -g_y_0_xyyyz_xyyzz[k] * cd_x[k] + g_y_0_xyyyz_xxyyzz[k];

                g_y_0_xxyyyz_xyzzz[k] = -g_y_0_xyyyz_xyzzz[k] * cd_x[k] + g_y_0_xyyyz_xxyzzz[k];

                g_y_0_xxyyyz_xzzzz[k] = -g_y_0_xyyyz_xzzzz[k] * cd_x[k] + g_y_0_xyyyz_xxzzzz[k];

                g_y_0_xxyyyz_yyyyy[k] = -g_y_0_xyyyz_yyyyy[k] * cd_x[k] + g_y_0_xyyyz_xyyyyy[k];

                g_y_0_xxyyyz_yyyyz[k] = -g_y_0_xyyyz_yyyyz[k] * cd_x[k] + g_y_0_xyyyz_xyyyyz[k];

                g_y_0_xxyyyz_yyyzz[k] = -g_y_0_xyyyz_yyyzz[k] * cd_x[k] + g_y_0_xyyyz_xyyyzz[k];

                g_y_0_xxyyyz_yyzzz[k] = -g_y_0_xyyyz_yyzzz[k] * cd_x[k] + g_y_0_xyyyz_xyyzzz[k];

                g_y_0_xxyyyz_yzzzz[k] = -g_y_0_xyyyz_yzzzz[k] * cd_x[k] + g_y_0_xyyyz_xyzzzz[k];

                g_y_0_xxyyyz_zzzzz[k] = -g_y_0_xyyyz_zzzzz[k] * cd_x[k] + g_y_0_xyyyz_xzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 252);

            auto g_y_0_xxyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 253);

            auto g_y_0_xxyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 254);

            auto g_y_0_xxyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 255);

            auto g_y_0_xxyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 256);

            auto g_y_0_xxyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 257);

            auto g_y_0_xxyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 258);

            auto g_y_0_xxyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 259);

            auto g_y_0_xxyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 260);

            auto g_y_0_xxyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 261);

            auto g_y_0_xxyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 262);

            auto g_y_0_xxyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 263);

            auto g_y_0_xxyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 264);

            auto g_y_0_xxyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 265);

            auto g_y_0_xxyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 266);

            auto g_y_0_xxyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 267);

            auto g_y_0_xxyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 268);

            auto g_y_0_xxyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 269);

            auto g_y_0_xxyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 270);

            auto g_y_0_xxyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 271);

            auto g_y_0_xxyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyzz_xxxxx, g_y_0_xxyyzz_xxxxy, g_y_0_xxyyzz_xxxxz, g_y_0_xxyyzz_xxxyy, g_y_0_xxyyzz_xxxyz, g_y_0_xxyyzz_xxxzz, g_y_0_xxyyzz_xxyyy, g_y_0_xxyyzz_xxyyz, g_y_0_xxyyzz_xxyzz, g_y_0_xxyyzz_xxzzz, g_y_0_xxyyzz_xyyyy, g_y_0_xxyyzz_xyyyz, g_y_0_xxyyzz_xyyzz, g_y_0_xxyyzz_xyzzz, g_y_0_xxyyzz_xzzzz, g_y_0_xxyyzz_yyyyy, g_y_0_xxyyzz_yyyyz, g_y_0_xxyyzz_yyyzz, g_y_0_xxyyzz_yyzzz, g_y_0_xxyyzz_yzzzz, g_y_0_xxyyzz_zzzzz, g_y_0_xyyzz_xxxxx, g_y_0_xyyzz_xxxxxx, g_y_0_xyyzz_xxxxxy, g_y_0_xyyzz_xxxxxz, g_y_0_xyyzz_xxxxy, g_y_0_xyyzz_xxxxyy, g_y_0_xyyzz_xxxxyz, g_y_0_xyyzz_xxxxz, g_y_0_xyyzz_xxxxzz, g_y_0_xyyzz_xxxyy, g_y_0_xyyzz_xxxyyy, g_y_0_xyyzz_xxxyyz, g_y_0_xyyzz_xxxyz, g_y_0_xyyzz_xxxyzz, g_y_0_xyyzz_xxxzz, g_y_0_xyyzz_xxxzzz, g_y_0_xyyzz_xxyyy, g_y_0_xyyzz_xxyyyy, g_y_0_xyyzz_xxyyyz, g_y_0_xyyzz_xxyyz, g_y_0_xyyzz_xxyyzz, g_y_0_xyyzz_xxyzz, g_y_0_xyyzz_xxyzzz, g_y_0_xyyzz_xxzzz, g_y_0_xyyzz_xxzzzz, g_y_0_xyyzz_xyyyy, g_y_0_xyyzz_xyyyyy, g_y_0_xyyzz_xyyyyz, g_y_0_xyyzz_xyyyz, g_y_0_xyyzz_xyyyzz, g_y_0_xyyzz_xyyzz, g_y_0_xyyzz_xyyzzz, g_y_0_xyyzz_xyzzz, g_y_0_xyyzz_xyzzzz, g_y_0_xyyzz_xzzzz, g_y_0_xyyzz_xzzzzz, g_y_0_xyyzz_yyyyy, g_y_0_xyyzz_yyyyz, g_y_0_xyyzz_yyyzz, g_y_0_xyyzz_yyzzz, g_y_0_xyyzz_yzzzz, g_y_0_xyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_xxxxx[k] = -g_y_0_xyyzz_xxxxx[k] * cd_x[k] + g_y_0_xyyzz_xxxxxx[k];

                g_y_0_xxyyzz_xxxxy[k] = -g_y_0_xyyzz_xxxxy[k] * cd_x[k] + g_y_0_xyyzz_xxxxxy[k];

                g_y_0_xxyyzz_xxxxz[k] = -g_y_0_xyyzz_xxxxz[k] * cd_x[k] + g_y_0_xyyzz_xxxxxz[k];

                g_y_0_xxyyzz_xxxyy[k] = -g_y_0_xyyzz_xxxyy[k] * cd_x[k] + g_y_0_xyyzz_xxxxyy[k];

                g_y_0_xxyyzz_xxxyz[k] = -g_y_0_xyyzz_xxxyz[k] * cd_x[k] + g_y_0_xyyzz_xxxxyz[k];

                g_y_0_xxyyzz_xxxzz[k] = -g_y_0_xyyzz_xxxzz[k] * cd_x[k] + g_y_0_xyyzz_xxxxzz[k];

                g_y_0_xxyyzz_xxyyy[k] = -g_y_0_xyyzz_xxyyy[k] * cd_x[k] + g_y_0_xyyzz_xxxyyy[k];

                g_y_0_xxyyzz_xxyyz[k] = -g_y_0_xyyzz_xxyyz[k] * cd_x[k] + g_y_0_xyyzz_xxxyyz[k];

                g_y_0_xxyyzz_xxyzz[k] = -g_y_0_xyyzz_xxyzz[k] * cd_x[k] + g_y_0_xyyzz_xxxyzz[k];

                g_y_0_xxyyzz_xxzzz[k] = -g_y_0_xyyzz_xxzzz[k] * cd_x[k] + g_y_0_xyyzz_xxxzzz[k];

                g_y_0_xxyyzz_xyyyy[k] = -g_y_0_xyyzz_xyyyy[k] * cd_x[k] + g_y_0_xyyzz_xxyyyy[k];

                g_y_0_xxyyzz_xyyyz[k] = -g_y_0_xyyzz_xyyyz[k] * cd_x[k] + g_y_0_xyyzz_xxyyyz[k];

                g_y_0_xxyyzz_xyyzz[k] = -g_y_0_xyyzz_xyyzz[k] * cd_x[k] + g_y_0_xyyzz_xxyyzz[k];

                g_y_0_xxyyzz_xyzzz[k] = -g_y_0_xyyzz_xyzzz[k] * cd_x[k] + g_y_0_xyyzz_xxyzzz[k];

                g_y_0_xxyyzz_xzzzz[k] = -g_y_0_xyyzz_xzzzz[k] * cd_x[k] + g_y_0_xyyzz_xxzzzz[k];

                g_y_0_xxyyzz_yyyyy[k] = -g_y_0_xyyzz_yyyyy[k] * cd_x[k] + g_y_0_xyyzz_xyyyyy[k];

                g_y_0_xxyyzz_yyyyz[k] = -g_y_0_xyyzz_yyyyz[k] * cd_x[k] + g_y_0_xyyzz_xyyyyz[k];

                g_y_0_xxyyzz_yyyzz[k] = -g_y_0_xyyzz_yyyzz[k] * cd_x[k] + g_y_0_xyyzz_xyyyzz[k];

                g_y_0_xxyyzz_yyzzz[k] = -g_y_0_xyyzz_yyzzz[k] * cd_x[k] + g_y_0_xyyzz_xyyzzz[k];

                g_y_0_xxyyzz_yzzzz[k] = -g_y_0_xyyzz_yzzzz[k] * cd_x[k] + g_y_0_xyyzz_xyzzzz[k];

                g_y_0_xxyyzz_zzzzz[k] = -g_y_0_xyyzz_zzzzz[k] * cd_x[k] + g_y_0_xyyzz_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 273);

            auto g_y_0_xxyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 274);

            auto g_y_0_xxyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 275);

            auto g_y_0_xxyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 276);

            auto g_y_0_xxyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 277);

            auto g_y_0_xxyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 278);

            auto g_y_0_xxyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 279);

            auto g_y_0_xxyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 280);

            auto g_y_0_xxyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 281);

            auto g_y_0_xxyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 282);

            auto g_y_0_xxyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 283);

            auto g_y_0_xxyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 284);

            auto g_y_0_xxyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 285);

            auto g_y_0_xxyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 286);

            auto g_y_0_xxyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 287);

            auto g_y_0_xxyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 288);

            auto g_y_0_xxyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 289);

            auto g_y_0_xxyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 290);

            auto g_y_0_xxyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 291);

            auto g_y_0_xxyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 292);

            auto g_y_0_xxyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzzz_xxxxx, g_y_0_xxyzzz_xxxxy, g_y_0_xxyzzz_xxxxz, g_y_0_xxyzzz_xxxyy, g_y_0_xxyzzz_xxxyz, g_y_0_xxyzzz_xxxzz, g_y_0_xxyzzz_xxyyy, g_y_0_xxyzzz_xxyyz, g_y_0_xxyzzz_xxyzz, g_y_0_xxyzzz_xxzzz, g_y_0_xxyzzz_xyyyy, g_y_0_xxyzzz_xyyyz, g_y_0_xxyzzz_xyyzz, g_y_0_xxyzzz_xyzzz, g_y_0_xxyzzz_xzzzz, g_y_0_xxyzzz_yyyyy, g_y_0_xxyzzz_yyyyz, g_y_0_xxyzzz_yyyzz, g_y_0_xxyzzz_yyzzz, g_y_0_xxyzzz_yzzzz, g_y_0_xxyzzz_zzzzz, g_y_0_xyzzz_xxxxx, g_y_0_xyzzz_xxxxxx, g_y_0_xyzzz_xxxxxy, g_y_0_xyzzz_xxxxxz, g_y_0_xyzzz_xxxxy, g_y_0_xyzzz_xxxxyy, g_y_0_xyzzz_xxxxyz, g_y_0_xyzzz_xxxxz, g_y_0_xyzzz_xxxxzz, g_y_0_xyzzz_xxxyy, g_y_0_xyzzz_xxxyyy, g_y_0_xyzzz_xxxyyz, g_y_0_xyzzz_xxxyz, g_y_0_xyzzz_xxxyzz, g_y_0_xyzzz_xxxzz, g_y_0_xyzzz_xxxzzz, g_y_0_xyzzz_xxyyy, g_y_0_xyzzz_xxyyyy, g_y_0_xyzzz_xxyyyz, g_y_0_xyzzz_xxyyz, g_y_0_xyzzz_xxyyzz, g_y_0_xyzzz_xxyzz, g_y_0_xyzzz_xxyzzz, g_y_0_xyzzz_xxzzz, g_y_0_xyzzz_xxzzzz, g_y_0_xyzzz_xyyyy, g_y_0_xyzzz_xyyyyy, g_y_0_xyzzz_xyyyyz, g_y_0_xyzzz_xyyyz, g_y_0_xyzzz_xyyyzz, g_y_0_xyzzz_xyyzz, g_y_0_xyzzz_xyyzzz, g_y_0_xyzzz_xyzzz, g_y_0_xyzzz_xyzzzz, g_y_0_xyzzz_xzzzz, g_y_0_xyzzz_xzzzzz, g_y_0_xyzzz_yyyyy, g_y_0_xyzzz_yyyyz, g_y_0_xyzzz_yyyzz, g_y_0_xyzzz_yyzzz, g_y_0_xyzzz_yzzzz, g_y_0_xyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_xxxxx[k] = -g_y_0_xyzzz_xxxxx[k] * cd_x[k] + g_y_0_xyzzz_xxxxxx[k];

                g_y_0_xxyzzz_xxxxy[k] = -g_y_0_xyzzz_xxxxy[k] * cd_x[k] + g_y_0_xyzzz_xxxxxy[k];

                g_y_0_xxyzzz_xxxxz[k] = -g_y_0_xyzzz_xxxxz[k] * cd_x[k] + g_y_0_xyzzz_xxxxxz[k];

                g_y_0_xxyzzz_xxxyy[k] = -g_y_0_xyzzz_xxxyy[k] * cd_x[k] + g_y_0_xyzzz_xxxxyy[k];

                g_y_0_xxyzzz_xxxyz[k] = -g_y_0_xyzzz_xxxyz[k] * cd_x[k] + g_y_0_xyzzz_xxxxyz[k];

                g_y_0_xxyzzz_xxxzz[k] = -g_y_0_xyzzz_xxxzz[k] * cd_x[k] + g_y_0_xyzzz_xxxxzz[k];

                g_y_0_xxyzzz_xxyyy[k] = -g_y_0_xyzzz_xxyyy[k] * cd_x[k] + g_y_0_xyzzz_xxxyyy[k];

                g_y_0_xxyzzz_xxyyz[k] = -g_y_0_xyzzz_xxyyz[k] * cd_x[k] + g_y_0_xyzzz_xxxyyz[k];

                g_y_0_xxyzzz_xxyzz[k] = -g_y_0_xyzzz_xxyzz[k] * cd_x[k] + g_y_0_xyzzz_xxxyzz[k];

                g_y_0_xxyzzz_xxzzz[k] = -g_y_0_xyzzz_xxzzz[k] * cd_x[k] + g_y_0_xyzzz_xxxzzz[k];

                g_y_0_xxyzzz_xyyyy[k] = -g_y_0_xyzzz_xyyyy[k] * cd_x[k] + g_y_0_xyzzz_xxyyyy[k];

                g_y_0_xxyzzz_xyyyz[k] = -g_y_0_xyzzz_xyyyz[k] * cd_x[k] + g_y_0_xyzzz_xxyyyz[k];

                g_y_0_xxyzzz_xyyzz[k] = -g_y_0_xyzzz_xyyzz[k] * cd_x[k] + g_y_0_xyzzz_xxyyzz[k];

                g_y_0_xxyzzz_xyzzz[k] = -g_y_0_xyzzz_xyzzz[k] * cd_x[k] + g_y_0_xyzzz_xxyzzz[k];

                g_y_0_xxyzzz_xzzzz[k] = -g_y_0_xyzzz_xzzzz[k] * cd_x[k] + g_y_0_xyzzz_xxzzzz[k];

                g_y_0_xxyzzz_yyyyy[k] = -g_y_0_xyzzz_yyyyy[k] * cd_x[k] + g_y_0_xyzzz_xyyyyy[k];

                g_y_0_xxyzzz_yyyyz[k] = -g_y_0_xyzzz_yyyyz[k] * cd_x[k] + g_y_0_xyzzz_xyyyyz[k];

                g_y_0_xxyzzz_yyyzz[k] = -g_y_0_xyzzz_yyyzz[k] * cd_x[k] + g_y_0_xyzzz_xyyyzz[k];

                g_y_0_xxyzzz_yyzzz[k] = -g_y_0_xyzzz_yyzzz[k] * cd_x[k] + g_y_0_xyzzz_xyyzzz[k];

                g_y_0_xxyzzz_yzzzz[k] = -g_y_0_xyzzz_yzzzz[k] * cd_x[k] + g_y_0_xyzzz_xyzzzz[k];

                g_y_0_xxyzzz_zzzzz[k] = -g_y_0_xyzzz_zzzzz[k] * cd_x[k] + g_y_0_xyzzz_xzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 294);

            auto g_y_0_xxzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 295);

            auto g_y_0_xxzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 296);

            auto g_y_0_xxzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 297);

            auto g_y_0_xxzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 298);

            auto g_y_0_xxzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 299);

            auto g_y_0_xxzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 300);

            auto g_y_0_xxzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 301);

            auto g_y_0_xxzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 302);

            auto g_y_0_xxzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 303);

            auto g_y_0_xxzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 304);

            auto g_y_0_xxzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 305);

            auto g_y_0_xxzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 306);

            auto g_y_0_xxzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 307);

            auto g_y_0_xxzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 308);

            auto g_y_0_xxzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 309);

            auto g_y_0_xxzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 310);

            auto g_y_0_xxzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 311);

            auto g_y_0_xxzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 312);

            auto g_y_0_xxzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 313);

            auto g_y_0_xxzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzzz_xxxxx, g_y_0_xxzzzz_xxxxy, g_y_0_xxzzzz_xxxxz, g_y_0_xxzzzz_xxxyy, g_y_0_xxzzzz_xxxyz, g_y_0_xxzzzz_xxxzz, g_y_0_xxzzzz_xxyyy, g_y_0_xxzzzz_xxyyz, g_y_0_xxzzzz_xxyzz, g_y_0_xxzzzz_xxzzz, g_y_0_xxzzzz_xyyyy, g_y_0_xxzzzz_xyyyz, g_y_0_xxzzzz_xyyzz, g_y_0_xxzzzz_xyzzz, g_y_0_xxzzzz_xzzzz, g_y_0_xxzzzz_yyyyy, g_y_0_xxzzzz_yyyyz, g_y_0_xxzzzz_yyyzz, g_y_0_xxzzzz_yyzzz, g_y_0_xxzzzz_yzzzz, g_y_0_xxzzzz_zzzzz, g_y_0_xzzzz_xxxxx, g_y_0_xzzzz_xxxxxx, g_y_0_xzzzz_xxxxxy, g_y_0_xzzzz_xxxxxz, g_y_0_xzzzz_xxxxy, g_y_0_xzzzz_xxxxyy, g_y_0_xzzzz_xxxxyz, g_y_0_xzzzz_xxxxz, g_y_0_xzzzz_xxxxzz, g_y_0_xzzzz_xxxyy, g_y_0_xzzzz_xxxyyy, g_y_0_xzzzz_xxxyyz, g_y_0_xzzzz_xxxyz, g_y_0_xzzzz_xxxyzz, g_y_0_xzzzz_xxxzz, g_y_0_xzzzz_xxxzzz, g_y_0_xzzzz_xxyyy, g_y_0_xzzzz_xxyyyy, g_y_0_xzzzz_xxyyyz, g_y_0_xzzzz_xxyyz, g_y_0_xzzzz_xxyyzz, g_y_0_xzzzz_xxyzz, g_y_0_xzzzz_xxyzzz, g_y_0_xzzzz_xxzzz, g_y_0_xzzzz_xxzzzz, g_y_0_xzzzz_xyyyy, g_y_0_xzzzz_xyyyyy, g_y_0_xzzzz_xyyyyz, g_y_0_xzzzz_xyyyz, g_y_0_xzzzz_xyyyzz, g_y_0_xzzzz_xyyzz, g_y_0_xzzzz_xyyzzz, g_y_0_xzzzz_xyzzz, g_y_0_xzzzz_xyzzzz, g_y_0_xzzzz_xzzzz, g_y_0_xzzzz_xzzzzz, g_y_0_xzzzz_yyyyy, g_y_0_xzzzz_yyyyz, g_y_0_xzzzz_yyyzz, g_y_0_xzzzz_yyzzz, g_y_0_xzzzz_yzzzz, g_y_0_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_xxxxx[k] = -g_y_0_xzzzz_xxxxx[k] * cd_x[k] + g_y_0_xzzzz_xxxxxx[k];

                g_y_0_xxzzzz_xxxxy[k] = -g_y_0_xzzzz_xxxxy[k] * cd_x[k] + g_y_0_xzzzz_xxxxxy[k];

                g_y_0_xxzzzz_xxxxz[k] = -g_y_0_xzzzz_xxxxz[k] * cd_x[k] + g_y_0_xzzzz_xxxxxz[k];

                g_y_0_xxzzzz_xxxyy[k] = -g_y_0_xzzzz_xxxyy[k] * cd_x[k] + g_y_0_xzzzz_xxxxyy[k];

                g_y_0_xxzzzz_xxxyz[k] = -g_y_0_xzzzz_xxxyz[k] * cd_x[k] + g_y_0_xzzzz_xxxxyz[k];

                g_y_0_xxzzzz_xxxzz[k] = -g_y_0_xzzzz_xxxzz[k] * cd_x[k] + g_y_0_xzzzz_xxxxzz[k];

                g_y_0_xxzzzz_xxyyy[k] = -g_y_0_xzzzz_xxyyy[k] * cd_x[k] + g_y_0_xzzzz_xxxyyy[k];

                g_y_0_xxzzzz_xxyyz[k] = -g_y_0_xzzzz_xxyyz[k] * cd_x[k] + g_y_0_xzzzz_xxxyyz[k];

                g_y_0_xxzzzz_xxyzz[k] = -g_y_0_xzzzz_xxyzz[k] * cd_x[k] + g_y_0_xzzzz_xxxyzz[k];

                g_y_0_xxzzzz_xxzzz[k] = -g_y_0_xzzzz_xxzzz[k] * cd_x[k] + g_y_0_xzzzz_xxxzzz[k];

                g_y_0_xxzzzz_xyyyy[k] = -g_y_0_xzzzz_xyyyy[k] * cd_x[k] + g_y_0_xzzzz_xxyyyy[k];

                g_y_0_xxzzzz_xyyyz[k] = -g_y_0_xzzzz_xyyyz[k] * cd_x[k] + g_y_0_xzzzz_xxyyyz[k];

                g_y_0_xxzzzz_xyyzz[k] = -g_y_0_xzzzz_xyyzz[k] * cd_x[k] + g_y_0_xzzzz_xxyyzz[k];

                g_y_0_xxzzzz_xyzzz[k] = -g_y_0_xzzzz_xyzzz[k] * cd_x[k] + g_y_0_xzzzz_xxyzzz[k];

                g_y_0_xxzzzz_xzzzz[k] = -g_y_0_xzzzz_xzzzz[k] * cd_x[k] + g_y_0_xzzzz_xxzzzz[k];

                g_y_0_xxzzzz_yyyyy[k] = -g_y_0_xzzzz_yyyyy[k] * cd_x[k] + g_y_0_xzzzz_xyyyyy[k];

                g_y_0_xxzzzz_yyyyz[k] = -g_y_0_xzzzz_yyyyz[k] * cd_x[k] + g_y_0_xzzzz_xyyyyz[k];

                g_y_0_xxzzzz_yyyzz[k] = -g_y_0_xzzzz_yyyzz[k] * cd_x[k] + g_y_0_xzzzz_xyyyzz[k];

                g_y_0_xxzzzz_yyzzz[k] = -g_y_0_xzzzz_yyzzz[k] * cd_x[k] + g_y_0_xzzzz_xyyzzz[k];

                g_y_0_xxzzzz_yzzzz[k] = -g_y_0_xzzzz_yzzzz[k] * cd_x[k] + g_y_0_xzzzz_xyzzzz[k];

                g_y_0_xxzzzz_zzzzz[k] = -g_y_0_xzzzz_zzzzz[k] * cd_x[k] + g_y_0_xzzzz_xzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 315);

            auto g_y_0_xyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 316);

            auto g_y_0_xyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 317);

            auto g_y_0_xyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 318);

            auto g_y_0_xyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 319);

            auto g_y_0_xyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 320);

            auto g_y_0_xyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 321);

            auto g_y_0_xyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 322);

            auto g_y_0_xyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 323);

            auto g_y_0_xyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 324);

            auto g_y_0_xyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 325);

            auto g_y_0_xyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 326);

            auto g_y_0_xyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 327);

            auto g_y_0_xyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 328);

            auto g_y_0_xyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 329);

            auto g_y_0_xyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 330);

            auto g_y_0_xyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 331);

            auto g_y_0_xyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 332);

            auto g_y_0_xyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 333);

            auto g_y_0_xyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 334);

            auto g_y_0_xyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyy_xxxxx, g_y_0_xyyyyy_xxxxy, g_y_0_xyyyyy_xxxxz, g_y_0_xyyyyy_xxxyy, g_y_0_xyyyyy_xxxyz, g_y_0_xyyyyy_xxxzz, g_y_0_xyyyyy_xxyyy, g_y_0_xyyyyy_xxyyz, g_y_0_xyyyyy_xxyzz, g_y_0_xyyyyy_xxzzz, g_y_0_xyyyyy_xyyyy, g_y_0_xyyyyy_xyyyz, g_y_0_xyyyyy_xyyzz, g_y_0_xyyyyy_xyzzz, g_y_0_xyyyyy_xzzzz, g_y_0_xyyyyy_yyyyy, g_y_0_xyyyyy_yyyyz, g_y_0_xyyyyy_yyyzz, g_y_0_xyyyyy_yyzzz, g_y_0_xyyyyy_yzzzz, g_y_0_xyyyyy_zzzzz, g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxxx, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_xxxxx[k] = -g_y_0_yyyyy_xxxxx[k] * cd_x[k] + g_y_0_yyyyy_xxxxxx[k];

                g_y_0_xyyyyy_xxxxy[k] = -g_y_0_yyyyy_xxxxy[k] * cd_x[k] + g_y_0_yyyyy_xxxxxy[k];

                g_y_0_xyyyyy_xxxxz[k] = -g_y_0_yyyyy_xxxxz[k] * cd_x[k] + g_y_0_yyyyy_xxxxxz[k];

                g_y_0_xyyyyy_xxxyy[k] = -g_y_0_yyyyy_xxxyy[k] * cd_x[k] + g_y_0_yyyyy_xxxxyy[k];

                g_y_0_xyyyyy_xxxyz[k] = -g_y_0_yyyyy_xxxyz[k] * cd_x[k] + g_y_0_yyyyy_xxxxyz[k];

                g_y_0_xyyyyy_xxxzz[k] = -g_y_0_yyyyy_xxxzz[k] * cd_x[k] + g_y_0_yyyyy_xxxxzz[k];

                g_y_0_xyyyyy_xxyyy[k] = -g_y_0_yyyyy_xxyyy[k] * cd_x[k] + g_y_0_yyyyy_xxxyyy[k];

                g_y_0_xyyyyy_xxyyz[k] = -g_y_0_yyyyy_xxyyz[k] * cd_x[k] + g_y_0_yyyyy_xxxyyz[k];

                g_y_0_xyyyyy_xxyzz[k] = -g_y_0_yyyyy_xxyzz[k] * cd_x[k] + g_y_0_yyyyy_xxxyzz[k];

                g_y_0_xyyyyy_xxzzz[k] = -g_y_0_yyyyy_xxzzz[k] * cd_x[k] + g_y_0_yyyyy_xxxzzz[k];

                g_y_0_xyyyyy_xyyyy[k] = -g_y_0_yyyyy_xyyyy[k] * cd_x[k] + g_y_0_yyyyy_xxyyyy[k];

                g_y_0_xyyyyy_xyyyz[k] = -g_y_0_yyyyy_xyyyz[k] * cd_x[k] + g_y_0_yyyyy_xxyyyz[k];

                g_y_0_xyyyyy_xyyzz[k] = -g_y_0_yyyyy_xyyzz[k] * cd_x[k] + g_y_0_yyyyy_xxyyzz[k];

                g_y_0_xyyyyy_xyzzz[k] = -g_y_0_yyyyy_xyzzz[k] * cd_x[k] + g_y_0_yyyyy_xxyzzz[k];

                g_y_0_xyyyyy_xzzzz[k] = -g_y_0_yyyyy_xzzzz[k] * cd_x[k] + g_y_0_yyyyy_xxzzzz[k];

                g_y_0_xyyyyy_yyyyy[k] = -g_y_0_yyyyy_yyyyy[k] * cd_x[k] + g_y_0_yyyyy_xyyyyy[k];

                g_y_0_xyyyyy_yyyyz[k] = -g_y_0_yyyyy_yyyyz[k] * cd_x[k] + g_y_0_yyyyy_xyyyyz[k];

                g_y_0_xyyyyy_yyyzz[k] = -g_y_0_yyyyy_yyyzz[k] * cd_x[k] + g_y_0_yyyyy_xyyyzz[k];

                g_y_0_xyyyyy_yyzzz[k] = -g_y_0_yyyyy_yyzzz[k] * cd_x[k] + g_y_0_yyyyy_xyyzzz[k];

                g_y_0_xyyyyy_yzzzz[k] = -g_y_0_yyyyy_yzzzz[k] * cd_x[k] + g_y_0_yyyyy_xyzzzz[k];

                g_y_0_xyyyyy_zzzzz[k] = -g_y_0_yyyyy_zzzzz[k] * cd_x[k] + g_y_0_yyyyy_xzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 336);

            auto g_y_0_xyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 337);

            auto g_y_0_xyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 338);

            auto g_y_0_xyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 339);

            auto g_y_0_xyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 340);

            auto g_y_0_xyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 341);

            auto g_y_0_xyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 342);

            auto g_y_0_xyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 343);

            auto g_y_0_xyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 344);

            auto g_y_0_xyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 345);

            auto g_y_0_xyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 346);

            auto g_y_0_xyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 347);

            auto g_y_0_xyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 348);

            auto g_y_0_xyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 349);

            auto g_y_0_xyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 350);

            auto g_y_0_xyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 351);

            auto g_y_0_xyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 352);

            auto g_y_0_xyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 353);

            auto g_y_0_xyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 354);

            auto g_y_0_xyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 355);

            auto g_y_0_xyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 356);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyz_xxxxx, g_y_0_xyyyyz_xxxxy, g_y_0_xyyyyz_xxxxz, g_y_0_xyyyyz_xxxyy, g_y_0_xyyyyz_xxxyz, g_y_0_xyyyyz_xxxzz, g_y_0_xyyyyz_xxyyy, g_y_0_xyyyyz_xxyyz, g_y_0_xyyyyz_xxyzz, g_y_0_xyyyyz_xxzzz, g_y_0_xyyyyz_xyyyy, g_y_0_xyyyyz_xyyyz, g_y_0_xyyyyz_xyyzz, g_y_0_xyyyyz_xyzzz, g_y_0_xyyyyz_xzzzz, g_y_0_xyyyyz_yyyyy, g_y_0_xyyyyz_yyyyz, g_y_0_xyyyyz_yyyzz, g_y_0_xyyyyz_yyzzz, g_y_0_xyyyyz_yzzzz, g_y_0_xyyyyz_zzzzz, g_y_0_yyyyz_xxxxx, g_y_0_yyyyz_xxxxxx, g_y_0_yyyyz_xxxxxy, g_y_0_yyyyz_xxxxxz, g_y_0_yyyyz_xxxxy, g_y_0_yyyyz_xxxxyy, g_y_0_yyyyz_xxxxyz, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxxzz, g_y_0_yyyyz_xxxyy, g_y_0_yyyyz_xxxyyy, g_y_0_yyyyz_xxxyyz, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxyzz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxxzzz, g_y_0_yyyyz_xxyyy, g_y_0_yyyyz_xxyyyy, g_y_0_yyyyz_xxyyyz, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyyzz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxyzzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xxzzzz, g_y_0_yyyyz_xyyyy, g_y_0_yyyyz_xyyyyy, g_y_0_yyyyz_xyyyyz, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyyzz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyyzzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xyzzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_xzzzzz, g_y_0_yyyyz_yyyyy, g_y_0_yyyyz_yyyyz, g_y_0_yyyyz_yyyzz, g_y_0_yyyyz_yyzzz, g_y_0_yyyyz_yzzzz, g_y_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_xxxxx[k] = -g_y_0_yyyyz_xxxxx[k] * cd_x[k] + g_y_0_yyyyz_xxxxxx[k];

                g_y_0_xyyyyz_xxxxy[k] = -g_y_0_yyyyz_xxxxy[k] * cd_x[k] + g_y_0_yyyyz_xxxxxy[k];

                g_y_0_xyyyyz_xxxxz[k] = -g_y_0_yyyyz_xxxxz[k] * cd_x[k] + g_y_0_yyyyz_xxxxxz[k];

                g_y_0_xyyyyz_xxxyy[k] = -g_y_0_yyyyz_xxxyy[k] * cd_x[k] + g_y_0_yyyyz_xxxxyy[k];

                g_y_0_xyyyyz_xxxyz[k] = -g_y_0_yyyyz_xxxyz[k] * cd_x[k] + g_y_0_yyyyz_xxxxyz[k];

                g_y_0_xyyyyz_xxxzz[k] = -g_y_0_yyyyz_xxxzz[k] * cd_x[k] + g_y_0_yyyyz_xxxxzz[k];

                g_y_0_xyyyyz_xxyyy[k] = -g_y_0_yyyyz_xxyyy[k] * cd_x[k] + g_y_0_yyyyz_xxxyyy[k];

                g_y_0_xyyyyz_xxyyz[k] = -g_y_0_yyyyz_xxyyz[k] * cd_x[k] + g_y_0_yyyyz_xxxyyz[k];

                g_y_0_xyyyyz_xxyzz[k] = -g_y_0_yyyyz_xxyzz[k] * cd_x[k] + g_y_0_yyyyz_xxxyzz[k];

                g_y_0_xyyyyz_xxzzz[k] = -g_y_0_yyyyz_xxzzz[k] * cd_x[k] + g_y_0_yyyyz_xxxzzz[k];

                g_y_0_xyyyyz_xyyyy[k] = -g_y_0_yyyyz_xyyyy[k] * cd_x[k] + g_y_0_yyyyz_xxyyyy[k];

                g_y_0_xyyyyz_xyyyz[k] = -g_y_0_yyyyz_xyyyz[k] * cd_x[k] + g_y_0_yyyyz_xxyyyz[k];

                g_y_0_xyyyyz_xyyzz[k] = -g_y_0_yyyyz_xyyzz[k] * cd_x[k] + g_y_0_yyyyz_xxyyzz[k];

                g_y_0_xyyyyz_xyzzz[k] = -g_y_0_yyyyz_xyzzz[k] * cd_x[k] + g_y_0_yyyyz_xxyzzz[k];

                g_y_0_xyyyyz_xzzzz[k] = -g_y_0_yyyyz_xzzzz[k] * cd_x[k] + g_y_0_yyyyz_xxzzzz[k];

                g_y_0_xyyyyz_yyyyy[k] = -g_y_0_yyyyz_yyyyy[k] * cd_x[k] + g_y_0_yyyyz_xyyyyy[k];

                g_y_0_xyyyyz_yyyyz[k] = -g_y_0_yyyyz_yyyyz[k] * cd_x[k] + g_y_0_yyyyz_xyyyyz[k];

                g_y_0_xyyyyz_yyyzz[k] = -g_y_0_yyyyz_yyyzz[k] * cd_x[k] + g_y_0_yyyyz_xyyyzz[k];

                g_y_0_xyyyyz_yyzzz[k] = -g_y_0_yyyyz_yyzzz[k] * cd_x[k] + g_y_0_yyyyz_xyyzzz[k];

                g_y_0_xyyyyz_yzzzz[k] = -g_y_0_yyyyz_yzzzz[k] * cd_x[k] + g_y_0_yyyyz_xyzzzz[k];

                g_y_0_xyyyyz_zzzzz[k] = -g_y_0_yyyyz_zzzzz[k] * cd_x[k] + g_y_0_yyyyz_xzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 357);

            auto g_y_0_xyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 358);

            auto g_y_0_xyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 359);

            auto g_y_0_xyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 360);

            auto g_y_0_xyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 361);

            auto g_y_0_xyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 362);

            auto g_y_0_xyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 363);

            auto g_y_0_xyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 364);

            auto g_y_0_xyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 365);

            auto g_y_0_xyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 366);

            auto g_y_0_xyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 367);

            auto g_y_0_xyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 368);

            auto g_y_0_xyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 369);

            auto g_y_0_xyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 370);

            auto g_y_0_xyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 371);

            auto g_y_0_xyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 372);

            auto g_y_0_xyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 373);

            auto g_y_0_xyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 374);

            auto g_y_0_xyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 375);

            auto g_y_0_xyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 376);

            auto g_y_0_xyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 377);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyzz_xxxxx, g_y_0_xyyyzz_xxxxy, g_y_0_xyyyzz_xxxxz, g_y_0_xyyyzz_xxxyy, g_y_0_xyyyzz_xxxyz, g_y_0_xyyyzz_xxxzz, g_y_0_xyyyzz_xxyyy, g_y_0_xyyyzz_xxyyz, g_y_0_xyyyzz_xxyzz, g_y_0_xyyyzz_xxzzz, g_y_0_xyyyzz_xyyyy, g_y_0_xyyyzz_xyyyz, g_y_0_xyyyzz_xyyzz, g_y_0_xyyyzz_xyzzz, g_y_0_xyyyzz_xzzzz, g_y_0_xyyyzz_yyyyy, g_y_0_xyyyzz_yyyyz, g_y_0_xyyyzz_yyyzz, g_y_0_xyyyzz_yyzzz, g_y_0_xyyyzz_yzzzz, g_y_0_xyyyzz_zzzzz, g_y_0_yyyzz_xxxxx, g_y_0_yyyzz_xxxxxx, g_y_0_yyyzz_xxxxxy, g_y_0_yyyzz_xxxxxz, g_y_0_yyyzz_xxxxy, g_y_0_yyyzz_xxxxyy, g_y_0_yyyzz_xxxxyz, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxxzz, g_y_0_yyyzz_xxxyy, g_y_0_yyyzz_xxxyyy, g_y_0_yyyzz_xxxyyz, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxyzz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxxzzz, g_y_0_yyyzz_xxyyy, g_y_0_yyyzz_xxyyyy, g_y_0_yyyzz_xxyyyz, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyyzz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxyzzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xxzzzz, g_y_0_yyyzz_xyyyy, g_y_0_yyyzz_xyyyyy, g_y_0_yyyzz_xyyyyz, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyyzz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyyzzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xyzzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_xzzzzz, g_y_0_yyyzz_yyyyy, g_y_0_yyyzz_yyyyz, g_y_0_yyyzz_yyyzz, g_y_0_yyyzz_yyzzz, g_y_0_yyyzz_yzzzz, g_y_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_xxxxx[k] = -g_y_0_yyyzz_xxxxx[k] * cd_x[k] + g_y_0_yyyzz_xxxxxx[k];

                g_y_0_xyyyzz_xxxxy[k] = -g_y_0_yyyzz_xxxxy[k] * cd_x[k] + g_y_0_yyyzz_xxxxxy[k];

                g_y_0_xyyyzz_xxxxz[k] = -g_y_0_yyyzz_xxxxz[k] * cd_x[k] + g_y_0_yyyzz_xxxxxz[k];

                g_y_0_xyyyzz_xxxyy[k] = -g_y_0_yyyzz_xxxyy[k] * cd_x[k] + g_y_0_yyyzz_xxxxyy[k];

                g_y_0_xyyyzz_xxxyz[k] = -g_y_0_yyyzz_xxxyz[k] * cd_x[k] + g_y_0_yyyzz_xxxxyz[k];

                g_y_0_xyyyzz_xxxzz[k] = -g_y_0_yyyzz_xxxzz[k] * cd_x[k] + g_y_0_yyyzz_xxxxzz[k];

                g_y_0_xyyyzz_xxyyy[k] = -g_y_0_yyyzz_xxyyy[k] * cd_x[k] + g_y_0_yyyzz_xxxyyy[k];

                g_y_0_xyyyzz_xxyyz[k] = -g_y_0_yyyzz_xxyyz[k] * cd_x[k] + g_y_0_yyyzz_xxxyyz[k];

                g_y_0_xyyyzz_xxyzz[k] = -g_y_0_yyyzz_xxyzz[k] * cd_x[k] + g_y_0_yyyzz_xxxyzz[k];

                g_y_0_xyyyzz_xxzzz[k] = -g_y_0_yyyzz_xxzzz[k] * cd_x[k] + g_y_0_yyyzz_xxxzzz[k];

                g_y_0_xyyyzz_xyyyy[k] = -g_y_0_yyyzz_xyyyy[k] * cd_x[k] + g_y_0_yyyzz_xxyyyy[k];

                g_y_0_xyyyzz_xyyyz[k] = -g_y_0_yyyzz_xyyyz[k] * cd_x[k] + g_y_0_yyyzz_xxyyyz[k];

                g_y_0_xyyyzz_xyyzz[k] = -g_y_0_yyyzz_xyyzz[k] * cd_x[k] + g_y_0_yyyzz_xxyyzz[k];

                g_y_0_xyyyzz_xyzzz[k] = -g_y_0_yyyzz_xyzzz[k] * cd_x[k] + g_y_0_yyyzz_xxyzzz[k];

                g_y_0_xyyyzz_xzzzz[k] = -g_y_0_yyyzz_xzzzz[k] * cd_x[k] + g_y_0_yyyzz_xxzzzz[k];

                g_y_0_xyyyzz_yyyyy[k] = -g_y_0_yyyzz_yyyyy[k] * cd_x[k] + g_y_0_yyyzz_xyyyyy[k];

                g_y_0_xyyyzz_yyyyz[k] = -g_y_0_yyyzz_yyyyz[k] * cd_x[k] + g_y_0_yyyzz_xyyyyz[k];

                g_y_0_xyyyzz_yyyzz[k] = -g_y_0_yyyzz_yyyzz[k] * cd_x[k] + g_y_0_yyyzz_xyyyzz[k];

                g_y_0_xyyyzz_yyzzz[k] = -g_y_0_yyyzz_yyzzz[k] * cd_x[k] + g_y_0_yyyzz_xyyzzz[k];

                g_y_0_xyyyzz_yzzzz[k] = -g_y_0_yyyzz_yzzzz[k] * cd_x[k] + g_y_0_yyyzz_xyzzzz[k];

                g_y_0_xyyyzz_zzzzz[k] = -g_y_0_yyyzz_zzzzz[k] * cd_x[k] + g_y_0_yyyzz_xzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 378);

            auto g_y_0_xyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 379);

            auto g_y_0_xyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 380);

            auto g_y_0_xyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 381);

            auto g_y_0_xyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 382);

            auto g_y_0_xyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 383);

            auto g_y_0_xyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 384);

            auto g_y_0_xyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 385);

            auto g_y_0_xyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 386);

            auto g_y_0_xyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 387);

            auto g_y_0_xyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 388);

            auto g_y_0_xyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 389);

            auto g_y_0_xyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 390);

            auto g_y_0_xyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 391);

            auto g_y_0_xyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 392);

            auto g_y_0_xyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 393);

            auto g_y_0_xyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 394);

            auto g_y_0_xyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 395);

            auto g_y_0_xyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 396);

            auto g_y_0_xyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 397);

            auto g_y_0_xyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 398);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzzz_xxxxx, g_y_0_xyyzzz_xxxxy, g_y_0_xyyzzz_xxxxz, g_y_0_xyyzzz_xxxyy, g_y_0_xyyzzz_xxxyz, g_y_0_xyyzzz_xxxzz, g_y_0_xyyzzz_xxyyy, g_y_0_xyyzzz_xxyyz, g_y_0_xyyzzz_xxyzz, g_y_0_xyyzzz_xxzzz, g_y_0_xyyzzz_xyyyy, g_y_0_xyyzzz_xyyyz, g_y_0_xyyzzz_xyyzz, g_y_0_xyyzzz_xyzzz, g_y_0_xyyzzz_xzzzz, g_y_0_xyyzzz_yyyyy, g_y_0_xyyzzz_yyyyz, g_y_0_xyyzzz_yyyzz, g_y_0_xyyzzz_yyzzz, g_y_0_xyyzzz_yzzzz, g_y_0_xyyzzz_zzzzz, g_y_0_yyzzz_xxxxx, g_y_0_yyzzz_xxxxxx, g_y_0_yyzzz_xxxxxy, g_y_0_yyzzz_xxxxxz, g_y_0_yyzzz_xxxxy, g_y_0_yyzzz_xxxxyy, g_y_0_yyzzz_xxxxyz, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxxzz, g_y_0_yyzzz_xxxyy, g_y_0_yyzzz_xxxyyy, g_y_0_yyzzz_xxxyyz, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxyzz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxxzzz, g_y_0_yyzzz_xxyyy, g_y_0_yyzzz_xxyyyy, g_y_0_yyzzz_xxyyyz, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyyzz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxyzzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xxzzzz, g_y_0_yyzzz_xyyyy, g_y_0_yyzzz_xyyyyy, g_y_0_yyzzz_xyyyyz, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyyzz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyyzzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xyzzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_xzzzzz, g_y_0_yyzzz_yyyyy, g_y_0_yyzzz_yyyyz, g_y_0_yyzzz_yyyzz, g_y_0_yyzzz_yyzzz, g_y_0_yyzzz_yzzzz, g_y_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_xxxxx[k] = -g_y_0_yyzzz_xxxxx[k] * cd_x[k] + g_y_0_yyzzz_xxxxxx[k];

                g_y_0_xyyzzz_xxxxy[k] = -g_y_0_yyzzz_xxxxy[k] * cd_x[k] + g_y_0_yyzzz_xxxxxy[k];

                g_y_0_xyyzzz_xxxxz[k] = -g_y_0_yyzzz_xxxxz[k] * cd_x[k] + g_y_0_yyzzz_xxxxxz[k];

                g_y_0_xyyzzz_xxxyy[k] = -g_y_0_yyzzz_xxxyy[k] * cd_x[k] + g_y_0_yyzzz_xxxxyy[k];

                g_y_0_xyyzzz_xxxyz[k] = -g_y_0_yyzzz_xxxyz[k] * cd_x[k] + g_y_0_yyzzz_xxxxyz[k];

                g_y_0_xyyzzz_xxxzz[k] = -g_y_0_yyzzz_xxxzz[k] * cd_x[k] + g_y_0_yyzzz_xxxxzz[k];

                g_y_0_xyyzzz_xxyyy[k] = -g_y_0_yyzzz_xxyyy[k] * cd_x[k] + g_y_0_yyzzz_xxxyyy[k];

                g_y_0_xyyzzz_xxyyz[k] = -g_y_0_yyzzz_xxyyz[k] * cd_x[k] + g_y_0_yyzzz_xxxyyz[k];

                g_y_0_xyyzzz_xxyzz[k] = -g_y_0_yyzzz_xxyzz[k] * cd_x[k] + g_y_0_yyzzz_xxxyzz[k];

                g_y_0_xyyzzz_xxzzz[k] = -g_y_0_yyzzz_xxzzz[k] * cd_x[k] + g_y_0_yyzzz_xxxzzz[k];

                g_y_0_xyyzzz_xyyyy[k] = -g_y_0_yyzzz_xyyyy[k] * cd_x[k] + g_y_0_yyzzz_xxyyyy[k];

                g_y_0_xyyzzz_xyyyz[k] = -g_y_0_yyzzz_xyyyz[k] * cd_x[k] + g_y_0_yyzzz_xxyyyz[k];

                g_y_0_xyyzzz_xyyzz[k] = -g_y_0_yyzzz_xyyzz[k] * cd_x[k] + g_y_0_yyzzz_xxyyzz[k];

                g_y_0_xyyzzz_xyzzz[k] = -g_y_0_yyzzz_xyzzz[k] * cd_x[k] + g_y_0_yyzzz_xxyzzz[k];

                g_y_0_xyyzzz_xzzzz[k] = -g_y_0_yyzzz_xzzzz[k] * cd_x[k] + g_y_0_yyzzz_xxzzzz[k];

                g_y_0_xyyzzz_yyyyy[k] = -g_y_0_yyzzz_yyyyy[k] * cd_x[k] + g_y_0_yyzzz_xyyyyy[k];

                g_y_0_xyyzzz_yyyyz[k] = -g_y_0_yyzzz_yyyyz[k] * cd_x[k] + g_y_0_yyzzz_xyyyyz[k];

                g_y_0_xyyzzz_yyyzz[k] = -g_y_0_yyzzz_yyyzz[k] * cd_x[k] + g_y_0_yyzzz_xyyyzz[k];

                g_y_0_xyyzzz_yyzzz[k] = -g_y_0_yyzzz_yyzzz[k] * cd_x[k] + g_y_0_yyzzz_xyyzzz[k];

                g_y_0_xyyzzz_yzzzz[k] = -g_y_0_yyzzz_yzzzz[k] * cd_x[k] + g_y_0_yyzzz_xyzzzz[k];

                g_y_0_xyyzzz_zzzzz[k] = -g_y_0_yyzzz_zzzzz[k] * cd_x[k] + g_y_0_yyzzz_xzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 399);

            auto g_y_0_xyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 400);

            auto g_y_0_xyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 401);

            auto g_y_0_xyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 402);

            auto g_y_0_xyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 403);

            auto g_y_0_xyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 404);

            auto g_y_0_xyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 405);

            auto g_y_0_xyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 406);

            auto g_y_0_xyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 407);

            auto g_y_0_xyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 408);

            auto g_y_0_xyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 409);

            auto g_y_0_xyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 410);

            auto g_y_0_xyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 411);

            auto g_y_0_xyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 412);

            auto g_y_0_xyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 413);

            auto g_y_0_xyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 414);

            auto g_y_0_xyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 415);

            auto g_y_0_xyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 416);

            auto g_y_0_xyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 417);

            auto g_y_0_xyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 418);

            auto g_y_0_xyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzzz_xxxxx, g_y_0_xyzzzz_xxxxy, g_y_0_xyzzzz_xxxxz, g_y_0_xyzzzz_xxxyy, g_y_0_xyzzzz_xxxyz, g_y_0_xyzzzz_xxxzz, g_y_0_xyzzzz_xxyyy, g_y_0_xyzzzz_xxyyz, g_y_0_xyzzzz_xxyzz, g_y_0_xyzzzz_xxzzz, g_y_0_xyzzzz_xyyyy, g_y_0_xyzzzz_xyyyz, g_y_0_xyzzzz_xyyzz, g_y_0_xyzzzz_xyzzz, g_y_0_xyzzzz_xzzzz, g_y_0_xyzzzz_yyyyy, g_y_0_xyzzzz_yyyyz, g_y_0_xyzzzz_yyyzz, g_y_0_xyzzzz_yyzzz, g_y_0_xyzzzz_yzzzz, g_y_0_xyzzzz_zzzzz, g_y_0_yzzzz_xxxxx, g_y_0_yzzzz_xxxxxx, g_y_0_yzzzz_xxxxxy, g_y_0_yzzzz_xxxxxz, g_y_0_yzzzz_xxxxy, g_y_0_yzzzz_xxxxyy, g_y_0_yzzzz_xxxxyz, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxxzz, g_y_0_yzzzz_xxxyy, g_y_0_yzzzz_xxxyyy, g_y_0_yzzzz_xxxyyz, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxyzz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxxzzz, g_y_0_yzzzz_xxyyy, g_y_0_yzzzz_xxyyyy, g_y_0_yzzzz_xxyyyz, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyyzz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxyzzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xxzzzz, g_y_0_yzzzz_xyyyy, g_y_0_yzzzz_xyyyyy, g_y_0_yzzzz_xyyyyz, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyyzz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyyzzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xyzzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_xzzzzz, g_y_0_yzzzz_yyyyy, g_y_0_yzzzz_yyyyz, g_y_0_yzzzz_yyyzz, g_y_0_yzzzz_yyzzz, g_y_0_yzzzz_yzzzz, g_y_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_xxxxx[k] = -g_y_0_yzzzz_xxxxx[k] * cd_x[k] + g_y_0_yzzzz_xxxxxx[k];

                g_y_0_xyzzzz_xxxxy[k] = -g_y_0_yzzzz_xxxxy[k] * cd_x[k] + g_y_0_yzzzz_xxxxxy[k];

                g_y_0_xyzzzz_xxxxz[k] = -g_y_0_yzzzz_xxxxz[k] * cd_x[k] + g_y_0_yzzzz_xxxxxz[k];

                g_y_0_xyzzzz_xxxyy[k] = -g_y_0_yzzzz_xxxyy[k] * cd_x[k] + g_y_0_yzzzz_xxxxyy[k];

                g_y_0_xyzzzz_xxxyz[k] = -g_y_0_yzzzz_xxxyz[k] * cd_x[k] + g_y_0_yzzzz_xxxxyz[k];

                g_y_0_xyzzzz_xxxzz[k] = -g_y_0_yzzzz_xxxzz[k] * cd_x[k] + g_y_0_yzzzz_xxxxzz[k];

                g_y_0_xyzzzz_xxyyy[k] = -g_y_0_yzzzz_xxyyy[k] * cd_x[k] + g_y_0_yzzzz_xxxyyy[k];

                g_y_0_xyzzzz_xxyyz[k] = -g_y_0_yzzzz_xxyyz[k] * cd_x[k] + g_y_0_yzzzz_xxxyyz[k];

                g_y_0_xyzzzz_xxyzz[k] = -g_y_0_yzzzz_xxyzz[k] * cd_x[k] + g_y_0_yzzzz_xxxyzz[k];

                g_y_0_xyzzzz_xxzzz[k] = -g_y_0_yzzzz_xxzzz[k] * cd_x[k] + g_y_0_yzzzz_xxxzzz[k];

                g_y_0_xyzzzz_xyyyy[k] = -g_y_0_yzzzz_xyyyy[k] * cd_x[k] + g_y_0_yzzzz_xxyyyy[k];

                g_y_0_xyzzzz_xyyyz[k] = -g_y_0_yzzzz_xyyyz[k] * cd_x[k] + g_y_0_yzzzz_xxyyyz[k];

                g_y_0_xyzzzz_xyyzz[k] = -g_y_0_yzzzz_xyyzz[k] * cd_x[k] + g_y_0_yzzzz_xxyyzz[k];

                g_y_0_xyzzzz_xyzzz[k] = -g_y_0_yzzzz_xyzzz[k] * cd_x[k] + g_y_0_yzzzz_xxyzzz[k];

                g_y_0_xyzzzz_xzzzz[k] = -g_y_0_yzzzz_xzzzz[k] * cd_x[k] + g_y_0_yzzzz_xxzzzz[k];

                g_y_0_xyzzzz_yyyyy[k] = -g_y_0_yzzzz_yyyyy[k] * cd_x[k] + g_y_0_yzzzz_xyyyyy[k];

                g_y_0_xyzzzz_yyyyz[k] = -g_y_0_yzzzz_yyyyz[k] * cd_x[k] + g_y_0_yzzzz_xyyyyz[k];

                g_y_0_xyzzzz_yyyzz[k] = -g_y_0_yzzzz_yyyzz[k] * cd_x[k] + g_y_0_yzzzz_xyyyzz[k];

                g_y_0_xyzzzz_yyzzz[k] = -g_y_0_yzzzz_yyzzz[k] * cd_x[k] + g_y_0_yzzzz_xyyzzz[k];

                g_y_0_xyzzzz_yzzzz[k] = -g_y_0_yzzzz_yzzzz[k] * cd_x[k] + g_y_0_yzzzz_xyzzzz[k];

                g_y_0_xyzzzz_zzzzz[k] = -g_y_0_yzzzz_zzzzz[k] * cd_x[k] + g_y_0_yzzzz_xzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 420);

            auto g_y_0_xzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 421);

            auto g_y_0_xzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 422);

            auto g_y_0_xzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 423);

            auto g_y_0_xzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 424);

            auto g_y_0_xzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 425);

            auto g_y_0_xzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 426);

            auto g_y_0_xzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 427);

            auto g_y_0_xzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 428);

            auto g_y_0_xzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 429);

            auto g_y_0_xzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 430);

            auto g_y_0_xzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 431);

            auto g_y_0_xzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 432);

            auto g_y_0_xzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 433);

            auto g_y_0_xzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 434);

            auto g_y_0_xzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 435);

            auto g_y_0_xzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 436);

            auto g_y_0_xzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 437);

            auto g_y_0_xzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 438);

            auto g_y_0_xzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 439);

            auto g_y_0_xzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 440);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzzz_xxxxx, g_y_0_xzzzzz_xxxxy, g_y_0_xzzzzz_xxxxz, g_y_0_xzzzzz_xxxyy, g_y_0_xzzzzz_xxxyz, g_y_0_xzzzzz_xxxzz, g_y_0_xzzzzz_xxyyy, g_y_0_xzzzzz_xxyyz, g_y_0_xzzzzz_xxyzz, g_y_0_xzzzzz_xxzzz, g_y_0_xzzzzz_xyyyy, g_y_0_xzzzzz_xyyyz, g_y_0_xzzzzz_xyyzz, g_y_0_xzzzzz_xyzzz, g_y_0_xzzzzz_xzzzz, g_y_0_xzzzzz_yyyyy, g_y_0_xzzzzz_yyyyz, g_y_0_xzzzzz_yyyzz, g_y_0_xzzzzz_yyzzz, g_y_0_xzzzzz_yzzzz, g_y_0_xzzzzz_zzzzz, g_y_0_zzzzz_xxxxx, g_y_0_zzzzz_xxxxxx, g_y_0_zzzzz_xxxxxy, g_y_0_zzzzz_xxxxxz, g_y_0_zzzzz_xxxxy, g_y_0_zzzzz_xxxxyy, g_y_0_zzzzz_xxxxyz, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxxzz, g_y_0_zzzzz_xxxyy, g_y_0_zzzzz_xxxyyy, g_y_0_zzzzz_xxxyyz, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxyzz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxxzzz, g_y_0_zzzzz_xxyyy, g_y_0_zzzzz_xxyyyy, g_y_0_zzzzz_xxyyyz, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyyzz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxyzzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xxzzzz, g_y_0_zzzzz_xyyyy, g_y_0_zzzzz_xyyyyy, g_y_0_zzzzz_xyyyyz, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyyzz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyyzzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xyzzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_xzzzzz, g_y_0_zzzzz_yyyyy, g_y_0_zzzzz_yyyyz, g_y_0_zzzzz_yyyzz, g_y_0_zzzzz_yyzzz, g_y_0_zzzzz_yzzzz, g_y_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_xxxxx[k] = -g_y_0_zzzzz_xxxxx[k] * cd_x[k] + g_y_0_zzzzz_xxxxxx[k];

                g_y_0_xzzzzz_xxxxy[k] = -g_y_0_zzzzz_xxxxy[k] * cd_x[k] + g_y_0_zzzzz_xxxxxy[k];

                g_y_0_xzzzzz_xxxxz[k] = -g_y_0_zzzzz_xxxxz[k] * cd_x[k] + g_y_0_zzzzz_xxxxxz[k];

                g_y_0_xzzzzz_xxxyy[k] = -g_y_0_zzzzz_xxxyy[k] * cd_x[k] + g_y_0_zzzzz_xxxxyy[k];

                g_y_0_xzzzzz_xxxyz[k] = -g_y_0_zzzzz_xxxyz[k] * cd_x[k] + g_y_0_zzzzz_xxxxyz[k];

                g_y_0_xzzzzz_xxxzz[k] = -g_y_0_zzzzz_xxxzz[k] * cd_x[k] + g_y_0_zzzzz_xxxxzz[k];

                g_y_0_xzzzzz_xxyyy[k] = -g_y_0_zzzzz_xxyyy[k] * cd_x[k] + g_y_0_zzzzz_xxxyyy[k];

                g_y_0_xzzzzz_xxyyz[k] = -g_y_0_zzzzz_xxyyz[k] * cd_x[k] + g_y_0_zzzzz_xxxyyz[k];

                g_y_0_xzzzzz_xxyzz[k] = -g_y_0_zzzzz_xxyzz[k] * cd_x[k] + g_y_0_zzzzz_xxxyzz[k];

                g_y_0_xzzzzz_xxzzz[k] = -g_y_0_zzzzz_xxzzz[k] * cd_x[k] + g_y_0_zzzzz_xxxzzz[k];

                g_y_0_xzzzzz_xyyyy[k] = -g_y_0_zzzzz_xyyyy[k] * cd_x[k] + g_y_0_zzzzz_xxyyyy[k];

                g_y_0_xzzzzz_xyyyz[k] = -g_y_0_zzzzz_xyyyz[k] * cd_x[k] + g_y_0_zzzzz_xxyyyz[k];

                g_y_0_xzzzzz_xyyzz[k] = -g_y_0_zzzzz_xyyzz[k] * cd_x[k] + g_y_0_zzzzz_xxyyzz[k];

                g_y_0_xzzzzz_xyzzz[k] = -g_y_0_zzzzz_xyzzz[k] * cd_x[k] + g_y_0_zzzzz_xxyzzz[k];

                g_y_0_xzzzzz_xzzzz[k] = -g_y_0_zzzzz_xzzzz[k] * cd_x[k] + g_y_0_zzzzz_xxzzzz[k];

                g_y_0_xzzzzz_yyyyy[k] = -g_y_0_zzzzz_yyyyy[k] * cd_x[k] + g_y_0_zzzzz_xyyyyy[k];

                g_y_0_xzzzzz_yyyyz[k] = -g_y_0_zzzzz_yyyyz[k] * cd_x[k] + g_y_0_zzzzz_xyyyyz[k];

                g_y_0_xzzzzz_yyyzz[k] = -g_y_0_zzzzz_yyyzz[k] * cd_x[k] + g_y_0_zzzzz_xyyyzz[k];

                g_y_0_xzzzzz_yyzzz[k] = -g_y_0_zzzzz_yyzzz[k] * cd_x[k] + g_y_0_zzzzz_xyyzzz[k];

                g_y_0_xzzzzz_yzzzz[k] = -g_y_0_zzzzz_yzzzz[k] * cd_x[k] + g_y_0_zzzzz_xyzzzz[k];

                g_y_0_xzzzzz_zzzzz[k] = -g_y_0_zzzzz_zzzzz[k] * cd_x[k] + g_y_0_zzzzz_xzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 441);

            auto g_y_0_yyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 442);

            auto g_y_0_yyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 443);

            auto g_y_0_yyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 444);

            auto g_y_0_yyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 445);

            auto g_y_0_yyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 446);

            auto g_y_0_yyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 447);

            auto g_y_0_yyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 448);

            auto g_y_0_yyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 449);

            auto g_y_0_yyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 450);

            auto g_y_0_yyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 451);

            auto g_y_0_yyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 452);

            auto g_y_0_yyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 453);

            auto g_y_0_yyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 454);

            auto g_y_0_yyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 455);

            auto g_y_0_yyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 456);

            auto g_y_0_yyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 457);

            auto g_y_0_yyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 458);

            auto g_y_0_yyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 459);

            auto g_y_0_yyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 460);

            auto g_y_0_yyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 461);

            #pragma omp simd aligned(cd_y, g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyyy, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_zzzzz, g_y_0_yyyyyy_xxxxx, g_y_0_yyyyyy_xxxxy, g_y_0_yyyyyy_xxxxz, g_y_0_yyyyyy_xxxyy, g_y_0_yyyyyy_xxxyz, g_y_0_yyyyyy_xxxzz, g_y_0_yyyyyy_xxyyy, g_y_0_yyyyyy_xxyyz, g_y_0_yyyyyy_xxyzz, g_y_0_yyyyyy_xxzzz, g_y_0_yyyyyy_xyyyy, g_y_0_yyyyyy_xyyyz, g_y_0_yyyyyy_xyyzz, g_y_0_yyyyyy_xyzzz, g_y_0_yyyyyy_xzzzz, g_y_0_yyyyyy_yyyyy, g_y_0_yyyyyy_yyyyz, g_y_0_yyyyyy_yyyzz, g_y_0_yyyyyy_yyzzz, g_y_0_yyyyyy_yzzzz, g_y_0_yyyyyy_zzzzz, g_yyyyy_xxxxx, g_yyyyy_xxxxy, g_yyyyy_xxxxz, g_yyyyy_xxxyy, g_yyyyy_xxxyz, g_yyyyy_xxxzz, g_yyyyy_xxyyy, g_yyyyy_xxyyz, g_yyyyy_xxyzz, g_yyyyy_xxzzz, g_yyyyy_xyyyy, g_yyyyy_xyyyz, g_yyyyy_xyyzz, g_yyyyy_xyzzz, g_yyyyy_xzzzz, g_yyyyy_yyyyy, g_yyyyy_yyyyz, g_yyyyy_yyyzz, g_yyyyy_yyzzz, g_yyyyy_yzzzz, g_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_xxxxx[k] = -g_yyyyy_xxxxx[k] - g_y_0_yyyyy_xxxxx[k] * cd_y[k] + g_y_0_yyyyy_xxxxxy[k];

                g_y_0_yyyyyy_xxxxy[k] = -g_yyyyy_xxxxy[k] - g_y_0_yyyyy_xxxxy[k] * cd_y[k] + g_y_0_yyyyy_xxxxyy[k];

                g_y_0_yyyyyy_xxxxz[k] = -g_yyyyy_xxxxz[k] - g_y_0_yyyyy_xxxxz[k] * cd_y[k] + g_y_0_yyyyy_xxxxyz[k];

                g_y_0_yyyyyy_xxxyy[k] = -g_yyyyy_xxxyy[k] - g_y_0_yyyyy_xxxyy[k] * cd_y[k] + g_y_0_yyyyy_xxxyyy[k];

                g_y_0_yyyyyy_xxxyz[k] = -g_yyyyy_xxxyz[k] - g_y_0_yyyyy_xxxyz[k] * cd_y[k] + g_y_0_yyyyy_xxxyyz[k];

                g_y_0_yyyyyy_xxxzz[k] = -g_yyyyy_xxxzz[k] - g_y_0_yyyyy_xxxzz[k] * cd_y[k] + g_y_0_yyyyy_xxxyzz[k];

                g_y_0_yyyyyy_xxyyy[k] = -g_yyyyy_xxyyy[k] - g_y_0_yyyyy_xxyyy[k] * cd_y[k] + g_y_0_yyyyy_xxyyyy[k];

                g_y_0_yyyyyy_xxyyz[k] = -g_yyyyy_xxyyz[k] - g_y_0_yyyyy_xxyyz[k] * cd_y[k] + g_y_0_yyyyy_xxyyyz[k];

                g_y_0_yyyyyy_xxyzz[k] = -g_yyyyy_xxyzz[k] - g_y_0_yyyyy_xxyzz[k] * cd_y[k] + g_y_0_yyyyy_xxyyzz[k];

                g_y_0_yyyyyy_xxzzz[k] = -g_yyyyy_xxzzz[k] - g_y_0_yyyyy_xxzzz[k] * cd_y[k] + g_y_0_yyyyy_xxyzzz[k];

                g_y_0_yyyyyy_xyyyy[k] = -g_yyyyy_xyyyy[k] - g_y_0_yyyyy_xyyyy[k] * cd_y[k] + g_y_0_yyyyy_xyyyyy[k];

                g_y_0_yyyyyy_xyyyz[k] = -g_yyyyy_xyyyz[k] - g_y_0_yyyyy_xyyyz[k] * cd_y[k] + g_y_0_yyyyy_xyyyyz[k];

                g_y_0_yyyyyy_xyyzz[k] = -g_yyyyy_xyyzz[k] - g_y_0_yyyyy_xyyzz[k] * cd_y[k] + g_y_0_yyyyy_xyyyzz[k];

                g_y_0_yyyyyy_xyzzz[k] = -g_yyyyy_xyzzz[k] - g_y_0_yyyyy_xyzzz[k] * cd_y[k] + g_y_0_yyyyy_xyyzzz[k];

                g_y_0_yyyyyy_xzzzz[k] = -g_yyyyy_xzzzz[k] - g_y_0_yyyyy_xzzzz[k] * cd_y[k] + g_y_0_yyyyy_xyzzzz[k];

                g_y_0_yyyyyy_yyyyy[k] = -g_yyyyy_yyyyy[k] - g_y_0_yyyyy_yyyyy[k] * cd_y[k] + g_y_0_yyyyy_yyyyyy[k];

                g_y_0_yyyyyy_yyyyz[k] = -g_yyyyy_yyyyz[k] - g_y_0_yyyyy_yyyyz[k] * cd_y[k] + g_y_0_yyyyy_yyyyyz[k];

                g_y_0_yyyyyy_yyyzz[k] = -g_yyyyy_yyyzz[k] - g_y_0_yyyyy_yyyzz[k] * cd_y[k] + g_y_0_yyyyy_yyyyzz[k];

                g_y_0_yyyyyy_yyzzz[k] = -g_yyyyy_yyzzz[k] - g_y_0_yyyyy_yyzzz[k] * cd_y[k] + g_y_0_yyyyy_yyyzzz[k];

                g_y_0_yyyyyy_yzzzz[k] = -g_yyyyy_yzzzz[k] - g_y_0_yyyyy_yzzzz[k] * cd_y[k] + g_y_0_yyyyy_yyzzzz[k];

                g_y_0_yyyyyy_zzzzz[k] = -g_yyyyy_zzzzz[k] - g_y_0_yyyyy_zzzzz[k] * cd_y[k] + g_y_0_yyyyy_yzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 462);

            auto g_y_0_yyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 463);

            auto g_y_0_yyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 464);

            auto g_y_0_yyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 465);

            auto g_y_0_yyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 466);

            auto g_y_0_yyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 467);

            auto g_y_0_yyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 468);

            auto g_y_0_yyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 469);

            auto g_y_0_yyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 470);

            auto g_y_0_yyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 471);

            auto g_y_0_yyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 472);

            auto g_y_0_yyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 473);

            auto g_y_0_yyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 474);

            auto g_y_0_yyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 475);

            auto g_y_0_yyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 476);

            auto g_y_0_yyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 477);

            auto g_y_0_yyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 478);

            auto g_y_0_yyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 479);

            auto g_y_0_yyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 480);

            auto g_y_0_yyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 481);

            auto g_y_0_yyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 482);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_zzzzz, g_y_0_yyyyy_zzzzzz, g_y_0_yyyyyz_xxxxx, g_y_0_yyyyyz_xxxxy, g_y_0_yyyyyz_xxxxz, g_y_0_yyyyyz_xxxyy, g_y_0_yyyyyz_xxxyz, g_y_0_yyyyyz_xxxzz, g_y_0_yyyyyz_xxyyy, g_y_0_yyyyyz_xxyyz, g_y_0_yyyyyz_xxyzz, g_y_0_yyyyyz_xxzzz, g_y_0_yyyyyz_xyyyy, g_y_0_yyyyyz_xyyyz, g_y_0_yyyyyz_xyyzz, g_y_0_yyyyyz_xyzzz, g_y_0_yyyyyz_xzzzz, g_y_0_yyyyyz_yyyyy, g_y_0_yyyyyz_yyyyz, g_y_0_yyyyyz_yyyzz, g_y_0_yyyyyz_yyzzz, g_y_0_yyyyyz_yzzzz, g_y_0_yyyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_xxxxx[k] = -g_y_0_yyyyy_xxxxx[k] * cd_z[k] + g_y_0_yyyyy_xxxxxz[k];

                g_y_0_yyyyyz_xxxxy[k] = -g_y_0_yyyyy_xxxxy[k] * cd_z[k] + g_y_0_yyyyy_xxxxyz[k];

                g_y_0_yyyyyz_xxxxz[k] = -g_y_0_yyyyy_xxxxz[k] * cd_z[k] + g_y_0_yyyyy_xxxxzz[k];

                g_y_0_yyyyyz_xxxyy[k] = -g_y_0_yyyyy_xxxyy[k] * cd_z[k] + g_y_0_yyyyy_xxxyyz[k];

                g_y_0_yyyyyz_xxxyz[k] = -g_y_0_yyyyy_xxxyz[k] * cd_z[k] + g_y_0_yyyyy_xxxyzz[k];

                g_y_0_yyyyyz_xxxzz[k] = -g_y_0_yyyyy_xxxzz[k] * cd_z[k] + g_y_0_yyyyy_xxxzzz[k];

                g_y_0_yyyyyz_xxyyy[k] = -g_y_0_yyyyy_xxyyy[k] * cd_z[k] + g_y_0_yyyyy_xxyyyz[k];

                g_y_0_yyyyyz_xxyyz[k] = -g_y_0_yyyyy_xxyyz[k] * cd_z[k] + g_y_0_yyyyy_xxyyzz[k];

                g_y_0_yyyyyz_xxyzz[k] = -g_y_0_yyyyy_xxyzz[k] * cd_z[k] + g_y_0_yyyyy_xxyzzz[k];

                g_y_0_yyyyyz_xxzzz[k] = -g_y_0_yyyyy_xxzzz[k] * cd_z[k] + g_y_0_yyyyy_xxzzzz[k];

                g_y_0_yyyyyz_xyyyy[k] = -g_y_0_yyyyy_xyyyy[k] * cd_z[k] + g_y_0_yyyyy_xyyyyz[k];

                g_y_0_yyyyyz_xyyyz[k] = -g_y_0_yyyyy_xyyyz[k] * cd_z[k] + g_y_0_yyyyy_xyyyzz[k];

                g_y_0_yyyyyz_xyyzz[k] = -g_y_0_yyyyy_xyyzz[k] * cd_z[k] + g_y_0_yyyyy_xyyzzz[k];

                g_y_0_yyyyyz_xyzzz[k] = -g_y_0_yyyyy_xyzzz[k] * cd_z[k] + g_y_0_yyyyy_xyzzzz[k];

                g_y_0_yyyyyz_xzzzz[k] = -g_y_0_yyyyy_xzzzz[k] * cd_z[k] + g_y_0_yyyyy_xzzzzz[k];

                g_y_0_yyyyyz_yyyyy[k] = -g_y_0_yyyyy_yyyyy[k] * cd_z[k] + g_y_0_yyyyy_yyyyyz[k];

                g_y_0_yyyyyz_yyyyz[k] = -g_y_0_yyyyy_yyyyz[k] * cd_z[k] + g_y_0_yyyyy_yyyyzz[k];

                g_y_0_yyyyyz_yyyzz[k] = -g_y_0_yyyyy_yyyzz[k] * cd_z[k] + g_y_0_yyyyy_yyyzzz[k];

                g_y_0_yyyyyz_yyzzz[k] = -g_y_0_yyyyy_yyzzz[k] * cd_z[k] + g_y_0_yyyyy_yyzzzz[k];

                g_y_0_yyyyyz_yzzzz[k] = -g_y_0_yyyyy_yzzzz[k] * cd_z[k] + g_y_0_yyyyy_yzzzzz[k];

                g_y_0_yyyyyz_zzzzz[k] = -g_y_0_yyyyy_zzzzz[k] * cd_z[k] + g_y_0_yyyyy_zzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 483);

            auto g_y_0_yyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 484);

            auto g_y_0_yyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 485);

            auto g_y_0_yyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 486);

            auto g_y_0_yyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 487);

            auto g_y_0_yyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 488);

            auto g_y_0_yyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 489);

            auto g_y_0_yyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 490);

            auto g_y_0_yyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 491);

            auto g_y_0_yyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 492);

            auto g_y_0_yyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 493);

            auto g_y_0_yyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 494);

            auto g_y_0_yyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 495);

            auto g_y_0_yyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 496);

            auto g_y_0_yyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 497);

            auto g_y_0_yyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 498);

            auto g_y_0_yyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 499);

            auto g_y_0_yyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 500);

            auto g_y_0_yyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 501);

            auto g_y_0_yyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 502);

            auto g_y_0_yyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 503);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyz_xxxxx, g_y_0_yyyyz_xxxxxz, g_y_0_yyyyz_xxxxy, g_y_0_yyyyz_xxxxyz, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxxzz, g_y_0_yyyyz_xxxyy, g_y_0_yyyyz_xxxyyz, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxyzz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxxzzz, g_y_0_yyyyz_xxyyy, g_y_0_yyyyz_xxyyyz, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyyzz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxyzzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xxzzzz, g_y_0_yyyyz_xyyyy, g_y_0_yyyyz_xyyyyz, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyyzz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyyzzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xyzzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_xzzzzz, g_y_0_yyyyz_yyyyy, g_y_0_yyyyz_yyyyyz, g_y_0_yyyyz_yyyyz, g_y_0_yyyyz_yyyyzz, g_y_0_yyyyz_yyyzz, g_y_0_yyyyz_yyyzzz, g_y_0_yyyyz_yyzzz, g_y_0_yyyyz_yyzzzz, g_y_0_yyyyz_yzzzz, g_y_0_yyyyz_yzzzzz, g_y_0_yyyyz_zzzzz, g_y_0_yyyyz_zzzzzz, g_y_0_yyyyzz_xxxxx, g_y_0_yyyyzz_xxxxy, g_y_0_yyyyzz_xxxxz, g_y_0_yyyyzz_xxxyy, g_y_0_yyyyzz_xxxyz, g_y_0_yyyyzz_xxxzz, g_y_0_yyyyzz_xxyyy, g_y_0_yyyyzz_xxyyz, g_y_0_yyyyzz_xxyzz, g_y_0_yyyyzz_xxzzz, g_y_0_yyyyzz_xyyyy, g_y_0_yyyyzz_xyyyz, g_y_0_yyyyzz_xyyzz, g_y_0_yyyyzz_xyzzz, g_y_0_yyyyzz_xzzzz, g_y_0_yyyyzz_yyyyy, g_y_0_yyyyzz_yyyyz, g_y_0_yyyyzz_yyyzz, g_y_0_yyyyzz_yyzzz, g_y_0_yyyyzz_yzzzz, g_y_0_yyyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_xxxxx[k] = -g_y_0_yyyyz_xxxxx[k] * cd_z[k] + g_y_0_yyyyz_xxxxxz[k];

                g_y_0_yyyyzz_xxxxy[k] = -g_y_0_yyyyz_xxxxy[k] * cd_z[k] + g_y_0_yyyyz_xxxxyz[k];

                g_y_0_yyyyzz_xxxxz[k] = -g_y_0_yyyyz_xxxxz[k] * cd_z[k] + g_y_0_yyyyz_xxxxzz[k];

                g_y_0_yyyyzz_xxxyy[k] = -g_y_0_yyyyz_xxxyy[k] * cd_z[k] + g_y_0_yyyyz_xxxyyz[k];

                g_y_0_yyyyzz_xxxyz[k] = -g_y_0_yyyyz_xxxyz[k] * cd_z[k] + g_y_0_yyyyz_xxxyzz[k];

                g_y_0_yyyyzz_xxxzz[k] = -g_y_0_yyyyz_xxxzz[k] * cd_z[k] + g_y_0_yyyyz_xxxzzz[k];

                g_y_0_yyyyzz_xxyyy[k] = -g_y_0_yyyyz_xxyyy[k] * cd_z[k] + g_y_0_yyyyz_xxyyyz[k];

                g_y_0_yyyyzz_xxyyz[k] = -g_y_0_yyyyz_xxyyz[k] * cd_z[k] + g_y_0_yyyyz_xxyyzz[k];

                g_y_0_yyyyzz_xxyzz[k] = -g_y_0_yyyyz_xxyzz[k] * cd_z[k] + g_y_0_yyyyz_xxyzzz[k];

                g_y_0_yyyyzz_xxzzz[k] = -g_y_0_yyyyz_xxzzz[k] * cd_z[k] + g_y_0_yyyyz_xxzzzz[k];

                g_y_0_yyyyzz_xyyyy[k] = -g_y_0_yyyyz_xyyyy[k] * cd_z[k] + g_y_0_yyyyz_xyyyyz[k];

                g_y_0_yyyyzz_xyyyz[k] = -g_y_0_yyyyz_xyyyz[k] * cd_z[k] + g_y_0_yyyyz_xyyyzz[k];

                g_y_0_yyyyzz_xyyzz[k] = -g_y_0_yyyyz_xyyzz[k] * cd_z[k] + g_y_0_yyyyz_xyyzzz[k];

                g_y_0_yyyyzz_xyzzz[k] = -g_y_0_yyyyz_xyzzz[k] * cd_z[k] + g_y_0_yyyyz_xyzzzz[k];

                g_y_0_yyyyzz_xzzzz[k] = -g_y_0_yyyyz_xzzzz[k] * cd_z[k] + g_y_0_yyyyz_xzzzzz[k];

                g_y_0_yyyyzz_yyyyy[k] = -g_y_0_yyyyz_yyyyy[k] * cd_z[k] + g_y_0_yyyyz_yyyyyz[k];

                g_y_0_yyyyzz_yyyyz[k] = -g_y_0_yyyyz_yyyyz[k] * cd_z[k] + g_y_0_yyyyz_yyyyzz[k];

                g_y_0_yyyyzz_yyyzz[k] = -g_y_0_yyyyz_yyyzz[k] * cd_z[k] + g_y_0_yyyyz_yyyzzz[k];

                g_y_0_yyyyzz_yyzzz[k] = -g_y_0_yyyyz_yyzzz[k] * cd_z[k] + g_y_0_yyyyz_yyzzzz[k];

                g_y_0_yyyyzz_yzzzz[k] = -g_y_0_yyyyz_yzzzz[k] * cd_z[k] + g_y_0_yyyyz_yzzzzz[k];

                g_y_0_yyyyzz_zzzzz[k] = -g_y_0_yyyyz_zzzzz[k] * cd_z[k] + g_y_0_yyyyz_zzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 504);

            auto g_y_0_yyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 505);

            auto g_y_0_yyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 506);

            auto g_y_0_yyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 507);

            auto g_y_0_yyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 508);

            auto g_y_0_yyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 509);

            auto g_y_0_yyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 510);

            auto g_y_0_yyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 511);

            auto g_y_0_yyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 512);

            auto g_y_0_yyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 513);

            auto g_y_0_yyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 514);

            auto g_y_0_yyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 515);

            auto g_y_0_yyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 516);

            auto g_y_0_yyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 517);

            auto g_y_0_yyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 518);

            auto g_y_0_yyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 519);

            auto g_y_0_yyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 520);

            auto g_y_0_yyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 521);

            auto g_y_0_yyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 522);

            auto g_y_0_yyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 523);

            auto g_y_0_yyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 524);

            #pragma omp simd aligned(cd_z, g_y_0_yyyzz_xxxxx, g_y_0_yyyzz_xxxxxz, g_y_0_yyyzz_xxxxy, g_y_0_yyyzz_xxxxyz, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxxzz, g_y_0_yyyzz_xxxyy, g_y_0_yyyzz_xxxyyz, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxyzz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxxzzz, g_y_0_yyyzz_xxyyy, g_y_0_yyyzz_xxyyyz, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyyzz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxyzzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xxzzzz, g_y_0_yyyzz_xyyyy, g_y_0_yyyzz_xyyyyz, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyyzz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyyzzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xyzzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_xzzzzz, g_y_0_yyyzz_yyyyy, g_y_0_yyyzz_yyyyyz, g_y_0_yyyzz_yyyyz, g_y_0_yyyzz_yyyyzz, g_y_0_yyyzz_yyyzz, g_y_0_yyyzz_yyyzzz, g_y_0_yyyzz_yyzzz, g_y_0_yyyzz_yyzzzz, g_y_0_yyyzz_yzzzz, g_y_0_yyyzz_yzzzzz, g_y_0_yyyzz_zzzzz, g_y_0_yyyzz_zzzzzz, g_y_0_yyyzzz_xxxxx, g_y_0_yyyzzz_xxxxy, g_y_0_yyyzzz_xxxxz, g_y_0_yyyzzz_xxxyy, g_y_0_yyyzzz_xxxyz, g_y_0_yyyzzz_xxxzz, g_y_0_yyyzzz_xxyyy, g_y_0_yyyzzz_xxyyz, g_y_0_yyyzzz_xxyzz, g_y_0_yyyzzz_xxzzz, g_y_0_yyyzzz_xyyyy, g_y_0_yyyzzz_xyyyz, g_y_0_yyyzzz_xyyzz, g_y_0_yyyzzz_xyzzz, g_y_0_yyyzzz_xzzzz, g_y_0_yyyzzz_yyyyy, g_y_0_yyyzzz_yyyyz, g_y_0_yyyzzz_yyyzz, g_y_0_yyyzzz_yyzzz, g_y_0_yyyzzz_yzzzz, g_y_0_yyyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_xxxxx[k] = -g_y_0_yyyzz_xxxxx[k] * cd_z[k] + g_y_0_yyyzz_xxxxxz[k];

                g_y_0_yyyzzz_xxxxy[k] = -g_y_0_yyyzz_xxxxy[k] * cd_z[k] + g_y_0_yyyzz_xxxxyz[k];

                g_y_0_yyyzzz_xxxxz[k] = -g_y_0_yyyzz_xxxxz[k] * cd_z[k] + g_y_0_yyyzz_xxxxzz[k];

                g_y_0_yyyzzz_xxxyy[k] = -g_y_0_yyyzz_xxxyy[k] * cd_z[k] + g_y_0_yyyzz_xxxyyz[k];

                g_y_0_yyyzzz_xxxyz[k] = -g_y_0_yyyzz_xxxyz[k] * cd_z[k] + g_y_0_yyyzz_xxxyzz[k];

                g_y_0_yyyzzz_xxxzz[k] = -g_y_0_yyyzz_xxxzz[k] * cd_z[k] + g_y_0_yyyzz_xxxzzz[k];

                g_y_0_yyyzzz_xxyyy[k] = -g_y_0_yyyzz_xxyyy[k] * cd_z[k] + g_y_0_yyyzz_xxyyyz[k];

                g_y_0_yyyzzz_xxyyz[k] = -g_y_0_yyyzz_xxyyz[k] * cd_z[k] + g_y_0_yyyzz_xxyyzz[k];

                g_y_0_yyyzzz_xxyzz[k] = -g_y_0_yyyzz_xxyzz[k] * cd_z[k] + g_y_0_yyyzz_xxyzzz[k];

                g_y_0_yyyzzz_xxzzz[k] = -g_y_0_yyyzz_xxzzz[k] * cd_z[k] + g_y_0_yyyzz_xxzzzz[k];

                g_y_0_yyyzzz_xyyyy[k] = -g_y_0_yyyzz_xyyyy[k] * cd_z[k] + g_y_0_yyyzz_xyyyyz[k];

                g_y_0_yyyzzz_xyyyz[k] = -g_y_0_yyyzz_xyyyz[k] * cd_z[k] + g_y_0_yyyzz_xyyyzz[k];

                g_y_0_yyyzzz_xyyzz[k] = -g_y_0_yyyzz_xyyzz[k] * cd_z[k] + g_y_0_yyyzz_xyyzzz[k];

                g_y_0_yyyzzz_xyzzz[k] = -g_y_0_yyyzz_xyzzz[k] * cd_z[k] + g_y_0_yyyzz_xyzzzz[k];

                g_y_0_yyyzzz_xzzzz[k] = -g_y_0_yyyzz_xzzzz[k] * cd_z[k] + g_y_0_yyyzz_xzzzzz[k];

                g_y_0_yyyzzz_yyyyy[k] = -g_y_0_yyyzz_yyyyy[k] * cd_z[k] + g_y_0_yyyzz_yyyyyz[k];

                g_y_0_yyyzzz_yyyyz[k] = -g_y_0_yyyzz_yyyyz[k] * cd_z[k] + g_y_0_yyyzz_yyyyzz[k];

                g_y_0_yyyzzz_yyyzz[k] = -g_y_0_yyyzz_yyyzz[k] * cd_z[k] + g_y_0_yyyzz_yyyzzz[k];

                g_y_0_yyyzzz_yyzzz[k] = -g_y_0_yyyzz_yyzzz[k] * cd_z[k] + g_y_0_yyyzz_yyzzzz[k];

                g_y_0_yyyzzz_yzzzz[k] = -g_y_0_yyyzz_yzzzz[k] * cd_z[k] + g_y_0_yyyzz_yzzzzz[k];

                g_y_0_yyyzzz_zzzzz[k] = -g_y_0_yyyzz_zzzzz[k] * cd_z[k] + g_y_0_yyyzz_zzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 525);

            auto g_y_0_yyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 526);

            auto g_y_0_yyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 527);

            auto g_y_0_yyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 528);

            auto g_y_0_yyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 529);

            auto g_y_0_yyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 530);

            auto g_y_0_yyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 531);

            auto g_y_0_yyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 532);

            auto g_y_0_yyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 533);

            auto g_y_0_yyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 534);

            auto g_y_0_yyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 535);

            auto g_y_0_yyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 536);

            auto g_y_0_yyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 537);

            auto g_y_0_yyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 538);

            auto g_y_0_yyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 539);

            auto g_y_0_yyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 540);

            auto g_y_0_yyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 541);

            auto g_y_0_yyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 542);

            auto g_y_0_yyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 543);

            auto g_y_0_yyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 544);

            auto g_y_0_yyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 545);

            #pragma omp simd aligned(cd_z, g_y_0_yyzzz_xxxxx, g_y_0_yyzzz_xxxxxz, g_y_0_yyzzz_xxxxy, g_y_0_yyzzz_xxxxyz, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxxzz, g_y_0_yyzzz_xxxyy, g_y_0_yyzzz_xxxyyz, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxyzz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxxzzz, g_y_0_yyzzz_xxyyy, g_y_0_yyzzz_xxyyyz, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyyzz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxyzzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xxzzzz, g_y_0_yyzzz_xyyyy, g_y_0_yyzzz_xyyyyz, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyyzz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyyzzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xyzzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_xzzzzz, g_y_0_yyzzz_yyyyy, g_y_0_yyzzz_yyyyyz, g_y_0_yyzzz_yyyyz, g_y_0_yyzzz_yyyyzz, g_y_0_yyzzz_yyyzz, g_y_0_yyzzz_yyyzzz, g_y_0_yyzzz_yyzzz, g_y_0_yyzzz_yyzzzz, g_y_0_yyzzz_yzzzz, g_y_0_yyzzz_yzzzzz, g_y_0_yyzzz_zzzzz, g_y_0_yyzzz_zzzzzz, g_y_0_yyzzzz_xxxxx, g_y_0_yyzzzz_xxxxy, g_y_0_yyzzzz_xxxxz, g_y_0_yyzzzz_xxxyy, g_y_0_yyzzzz_xxxyz, g_y_0_yyzzzz_xxxzz, g_y_0_yyzzzz_xxyyy, g_y_0_yyzzzz_xxyyz, g_y_0_yyzzzz_xxyzz, g_y_0_yyzzzz_xxzzz, g_y_0_yyzzzz_xyyyy, g_y_0_yyzzzz_xyyyz, g_y_0_yyzzzz_xyyzz, g_y_0_yyzzzz_xyzzz, g_y_0_yyzzzz_xzzzz, g_y_0_yyzzzz_yyyyy, g_y_0_yyzzzz_yyyyz, g_y_0_yyzzzz_yyyzz, g_y_0_yyzzzz_yyzzz, g_y_0_yyzzzz_yzzzz, g_y_0_yyzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_xxxxx[k] = -g_y_0_yyzzz_xxxxx[k] * cd_z[k] + g_y_0_yyzzz_xxxxxz[k];

                g_y_0_yyzzzz_xxxxy[k] = -g_y_0_yyzzz_xxxxy[k] * cd_z[k] + g_y_0_yyzzz_xxxxyz[k];

                g_y_0_yyzzzz_xxxxz[k] = -g_y_0_yyzzz_xxxxz[k] * cd_z[k] + g_y_0_yyzzz_xxxxzz[k];

                g_y_0_yyzzzz_xxxyy[k] = -g_y_0_yyzzz_xxxyy[k] * cd_z[k] + g_y_0_yyzzz_xxxyyz[k];

                g_y_0_yyzzzz_xxxyz[k] = -g_y_0_yyzzz_xxxyz[k] * cd_z[k] + g_y_0_yyzzz_xxxyzz[k];

                g_y_0_yyzzzz_xxxzz[k] = -g_y_0_yyzzz_xxxzz[k] * cd_z[k] + g_y_0_yyzzz_xxxzzz[k];

                g_y_0_yyzzzz_xxyyy[k] = -g_y_0_yyzzz_xxyyy[k] * cd_z[k] + g_y_0_yyzzz_xxyyyz[k];

                g_y_0_yyzzzz_xxyyz[k] = -g_y_0_yyzzz_xxyyz[k] * cd_z[k] + g_y_0_yyzzz_xxyyzz[k];

                g_y_0_yyzzzz_xxyzz[k] = -g_y_0_yyzzz_xxyzz[k] * cd_z[k] + g_y_0_yyzzz_xxyzzz[k];

                g_y_0_yyzzzz_xxzzz[k] = -g_y_0_yyzzz_xxzzz[k] * cd_z[k] + g_y_0_yyzzz_xxzzzz[k];

                g_y_0_yyzzzz_xyyyy[k] = -g_y_0_yyzzz_xyyyy[k] * cd_z[k] + g_y_0_yyzzz_xyyyyz[k];

                g_y_0_yyzzzz_xyyyz[k] = -g_y_0_yyzzz_xyyyz[k] * cd_z[k] + g_y_0_yyzzz_xyyyzz[k];

                g_y_0_yyzzzz_xyyzz[k] = -g_y_0_yyzzz_xyyzz[k] * cd_z[k] + g_y_0_yyzzz_xyyzzz[k];

                g_y_0_yyzzzz_xyzzz[k] = -g_y_0_yyzzz_xyzzz[k] * cd_z[k] + g_y_0_yyzzz_xyzzzz[k];

                g_y_0_yyzzzz_xzzzz[k] = -g_y_0_yyzzz_xzzzz[k] * cd_z[k] + g_y_0_yyzzz_xzzzzz[k];

                g_y_0_yyzzzz_yyyyy[k] = -g_y_0_yyzzz_yyyyy[k] * cd_z[k] + g_y_0_yyzzz_yyyyyz[k];

                g_y_0_yyzzzz_yyyyz[k] = -g_y_0_yyzzz_yyyyz[k] * cd_z[k] + g_y_0_yyzzz_yyyyzz[k];

                g_y_0_yyzzzz_yyyzz[k] = -g_y_0_yyzzz_yyyzz[k] * cd_z[k] + g_y_0_yyzzz_yyyzzz[k];

                g_y_0_yyzzzz_yyzzz[k] = -g_y_0_yyzzz_yyzzz[k] * cd_z[k] + g_y_0_yyzzz_yyzzzz[k];

                g_y_0_yyzzzz_yzzzz[k] = -g_y_0_yyzzz_yzzzz[k] * cd_z[k] + g_y_0_yyzzz_yzzzzz[k];

                g_y_0_yyzzzz_zzzzz[k] = -g_y_0_yyzzz_zzzzz[k] * cd_z[k] + g_y_0_yyzzz_zzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 546);

            auto g_y_0_yzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 547);

            auto g_y_0_yzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 548);

            auto g_y_0_yzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 549);

            auto g_y_0_yzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 550);

            auto g_y_0_yzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 551);

            auto g_y_0_yzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 552);

            auto g_y_0_yzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 553);

            auto g_y_0_yzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 554);

            auto g_y_0_yzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 555);

            auto g_y_0_yzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 556);

            auto g_y_0_yzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 557);

            auto g_y_0_yzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 558);

            auto g_y_0_yzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 559);

            auto g_y_0_yzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 560);

            auto g_y_0_yzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 561);

            auto g_y_0_yzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 562);

            auto g_y_0_yzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 563);

            auto g_y_0_yzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 564);

            auto g_y_0_yzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 565);

            auto g_y_0_yzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 566);

            #pragma omp simd aligned(cd_z, g_y_0_yzzzz_xxxxx, g_y_0_yzzzz_xxxxxz, g_y_0_yzzzz_xxxxy, g_y_0_yzzzz_xxxxyz, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxxzz, g_y_0_yzzzz_xxxyy, g_y_0_yzzzz_xxxyyz, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxyzz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxxzzz, g_y_0_yzzzz_xxyyy, g_y_0_yzzzz_xxyyyz, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyyzz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxyzzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xxzzzz, g_y_0_yzzzz_xyyyy, g_y_0_yzzzz_xyyyyz, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyyzz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyyzzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xyzzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_xzzzzz, g_y_0_yzzzz_yyyyy, g_y_0_yzzzz_yyyyyz, g_y_0_yzzzz_yyyyz, g_y_0_yzzzz_yyyyzz, g_y_0_yzzzz_yyyzz, g_y_0_yzzzz_yyyzzz, g_y_0_yzzzz_yyzzz, g_y_0_yzzzz_yyzzzz, g_y_0_yzzzz_yzzzz, g_y_0_yzzzz_yzzzzz, g_y_0_yzzzz_zzzzz, g_y_0_yzzzz_zzzzzz, g_y_0_yzzzzz_xxxxx, g_y_0_yzzzzz_xxxxy, g_y_0_yzzzzz_xxxxz, g_y_0_yzzzzz_xxxyy, g_y_0_yzzzzz_xxxyz, g_y_0_yzzzzz_xxxzz, g_y_0_yzzzzz_xxyyy, g_y_0_yzzzzz_xxyyz, g_y_0_yzzzzz_xxyzz, g_y_0_yzzzzz_xxzzz, g_y_0_yzzzzz_xyyyy, g_y_0_yzzzzz_xyyyz, g_y_0_yzzzzz_xyyzz, g_y_0_yzzzzz_xyzzz, g_y_0_yzzzzz_xzzzz, g_y_0_yzzzzz_yyyyy, g_y_0_yzzzzz_yyyyz, g_y_0_yzzzzz_yyyzz, g_y_0_yzzzzz_yyzzz, g_y_0_yzzzzz_yzzzz, g_y_0_yzzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_xxxxx[k] = -g_y_0_yzzzz_xxxxx[k] * cd_z[k] + g_y_0_yzzzz_xxxxxz[k];

                g_y_0_yzzzzz_xxxxy[k] = -g_y_0_yzzzz_xxxxy[k] * cd_z[k] + g_y_0_yzzzz_xxxxyz[k];

                g_y_0_yzzzzz_xxxxz[k] = -g_y_0_yzzzz_xxxxz[k] * cd_z[k] + g_y_0_yzzzz_xxxxzz[k];

                g_y_0_yzzzzz_xxxyy[k] = -g_y_0_yzzzz_xxxyy[k] * cd_z[k] + g_y_0_yzzzz_xxxyyz[k];

                g_y_0_yzzzzz_xxxyz[k] = -g_y_0_yzzzz_xxxyz[k] * cd_z[k] + g_y_0_yzzzz_xxxyzz[k];

                g_y_0_yzzzzz_xxxzz[k] = -g_y_0_yzzzz_xxxzz[k] * cd_z[k] + g_y_0_yzzzz_xxxzzz[k];

                g_y_0_yzzzzz_xxyyy[k] = -g_y_0_yzzzz_xxyyy[k] * cd_z[k] + g_y_0_yzzzz_xxyyyz[k];

                g_y_0_yzzzzz_xxyyz[k] = -g_y_0_yzzzz_xxyyz[k] * cd_z[k] + g_y_0_yzzzz_xxyyzz[k];

                g_y_0_yzzzzz_xxyzz[k] = -g_y_0_yzzzz_xxyzz[k] * cd_z[k] + g_y_0_yzzzz_xxyzzz[k];

                g_y_0_yzzzzz_xxzzz[k] = -g_y_0_yzzzz_xxzzz[k] * cd_z[k] + g_y_0_yzzzz_xxzzzz[k];

                g_y_0_yzzzzz_xyyyy[k] = -g_y_0_yzzzz_xyyyy[k] * cd_z[k] + g_y_0_yzzzz_xyyyyz[k];

                g_y_0_yzzzzz_xyyyz[k] = -g_y_0_yzzzz_xyyyz[k] * cd_z[k] + g_y_0_yzzzz_xyyyzz[k];

                g_y_0_yzzzzz_xyyzz[k] = -g_y_0_yzzzz_xyyzz[k] * cd_z[k] + g_y_0_yzzzz_xyyzzz[k];

                g_y_0_yzzzzz_xyzzz[k] = -g_y_0_yzzzz_xyzzz[k] * cd_z[k] + g_y_0_yzzzz_xyzzzz[k];

                g_y_0_yzzzzz_xzzzz[k] = -g_y_0_yzzzz_xzzzz[k] * cd_z[k] + g_y_0_yzzzz_xzzzzz[k];

                g_y_0_yzzzzz_yyyyy[k] = -g_y_0_yzzzz_yyyyy[k] * cd_z[k] + g_y_0_yzzzz_yyyyyz[k];

                g_y_0_yzzzzz_yyyyz[k] = -g_y_0_yzzzz_yyyyz[k] * cd_z[k] + g_y_0_yzzzz_yyyyzz[k];

                g_y_0_yzzzzz_yyyzz[k] = -g_y_0_yzzzz_yyyzz[k] * cd_z[k] + g_y_0_yzzzz_yyyzzz[k];

                g_y_0_yzzzzz_yyzzz[k] = -g_y_0_yzzzz_yyzzz[k] * cd_z[k] + g_y_0_yzzzz_yyzzzz[k];

                g_y_0_yzzzzz_yzzzz[k] = -g_y_0_yzzzz_yzzzz[k] * cd_z[k] + g_y_0_yzzzz_yzzzzz[k];

                g_y_0_yzzzzz_zzzzz[k] = -g_y_0_yzzzz_zzzzz[k] * cd_z[k] + g_y_0_yzzzz_zzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 567);

            auto g_y_0_zzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 568);

            auto g_y_0_zzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 569);

            auto g_y_0_zzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 570);

            auto g_y_0_zzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 571);

            auto g_y_0_zzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 572);

            auto g_y_0_zzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 573);

            auto g_y_0_zzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 574);

            auto g_y_0_zzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 575);

            auto g_y_0_zzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 576);

            auto g_y_0_zzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 577);

            auto g_y_0_zzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 578);

            auto g_y_0_zzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 579);

            auto g_y_0_zzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 580);

            auto g_y_0_zzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 581);

            auto g_y_0_zzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 582);

            auto g_y_0_zzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 583);

            auto g_y_0_zzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 584);

            auto g_y_0_zzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 585);

            auto g_y_0_zzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 586);

            auto g_y_0_zzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 588 * acomps * bcomps + 587);

            #pragma omp simd aligned(cd_z, g_y_0_zzzzz_xxxxx, g_y_0_zzzzz_xxxxxz, g_y_0_zzzzz_xxxxy, g_y_0_zzzzz_xxxxyz, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxxzz, g_y_0_zzzzz_xxxyy, g_y_0_zzzzz_xxxyyz, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxyzz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxxzzz, g_y_0_zzzzz_xxyyy, g_y_0_zzzzz_xxyyyz, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyyzz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxyzzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xxzzzz, g_y_0_zzzzz_xyyyy, g_y_0_zzzzz_xyyyyz, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyyzz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyyzzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xyzzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_xzzzzz, g_y_0_zzzzz_yyyyy, g_y_0_zzzzz_yyyyyz, g_y_0_zzzzz_yyyyz, g_y_0_zzzzz_yyyyzz, g_y_0_zzzzz_yyyzz, g_y_0_zzzzz_yyyzzz, g_y_0_zzzzz_yyzzz, g_y_0_zzzzz_yyzzzz, g_y_0_zzzzz_yzzzz, g_y_0_zzzzz_yzzzzz, g_y_0_zzzzz_zzzzz, g_y_0_zzzzz_zzzzzz, g_y_0_zzzzzz_xxxxx, g_y_0_zzzzzz_xxxxy, g_y_0_zzzzzz_xxxxz, g_y_0_zzzzzz_xxxyy, g_y_0_zzzzzz_xxxyz, g_y_0_zzzzzz_xxxzz, g_y_0_zzzzzz_xxyyy, g_y_0_zzzzzz_xxyyz, g_y_0_zzzzzz_xxyzz, g_y_0_zzzzzz_xxzzz, g_y_0_zzzzzz_xyyyy, g_y_0_zzzzzz_xyyyz, g_y_0_zzzzzz_xyyzz, g_y_0_zzzzzz_xyzzz, g_y_0_zzzzzz_xzzzz, g_y_0_zzzzzz_yyyyy, g_y_0_zzzzzz_yyyyz, g_y_0_zzzzzz_yyyzz, g_y_0_zzzzzz_yyzzz, g_y_0_zzzzzz_yzzzz, g_y_0_zzzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_xxxxx[k] = -g_y_0_zzzzz_xxxxx[k] * cd_z[k] + g_y_0_zzzzz_xxxxxz[k];

                g_y_0_zzzzzz_xxxxy[k] = -g_y_0_zzzzz_xxxxy[k] * cd_z[k] + g_y_0_zzzzz_xxxxyz[k];

                g_y_0_zzzzzz_xxxxz[k] = -g_y_0_zzzzz_xxxxz[k] * cd_z[k] + g_y_0_zzzzz_xxxxzz[k];

                g_y_0_zzzzzz_xxxyy[k] = -g_y_0_zzzzz_xxxyy[k] * cd_z[k] + g_y_0_zzzzz_xxxyyz[k];

                g_y_0_zzzzzz_xxxyz[k] = -g_y_0_zzzzz_xxxyz[k] * cd_z[k] + g_y_0_zzzzz_xxxyzz[k];

                g_y_0_zzzzzz_xxxzz[k] = -g_y_0_zzzzz_xxxzz[k] * cd_z[k] + g_y_0_zzzzz_xxxzzz[k];

                g_y_0_zzzzzz_xxyyy[k] = -g_y_0_zzzzz_xxyyy[k] * cd_z[k] + g_y_0_zzzzz_xxyyyz[k];

                g_y_0_zzzzzz_xxyyz[k] = -g_y_0_zzzzz_xxyyz[k] * cd_z[k] + g_y_0_zzzzz_xxyyzz[k];

                g_y_0_zzzzzz_xxyzz[k] = -g_y_0_zzzzz_xxyzz[k] * cd_z[k] + g_y_0_zzzzz_xxyzzz[k];

                g_y_0_zzzzzz_xxzzz[k] = -g_y_0_zzzzz_xxzzz[k] * cd_z[k] + g_y_0_zzzzz_xxzzzz[k];

                g_y_0_zzzzzz_xyyyy[k] = -g_y_0_zzzzz_xyyyy[k] * cd_z[k] + g_y_0_zzzzz_xyyyyz[k];

                g_y_0_zzzzzz_xyyyz[k] = -g_y_0_zzzzz_xyyyz[k] * cd_z[k] + g_y_0_zzzzz_xyyyzz[k];

                g_y_0_zzzzzz_xyyzz[k] = -g_y_0_zzzzz_xyyzz[k] * cd_z[k] + g_y_0_zzzzz_xyyzzz[k];

                g_y_0_zzzzzz_xyzzz[k] = -g_y_0_zzzzz_xyzzz[k] * cd_z[k] + g_y_0_zzzzz_xyzzzz[k];

                g_y_0_zzzzzz_xzzzz[k] = -g_y_0_zzzzz_xzzzz[k] * cd_z[k] + g_y_0_zzzzz_xzzzzz[k];

                g_y_0_zzzzzz_yyyyy[k] = -g_y_0_zzzzz_yyyyy[k] * cd_z[k] + g_y_0_zzzzz_yyyyyz[k];

                g_y_0_zzzzzz_yyyyz[k] = -g_y_0_zzzzz_yyyyz[k] * cd_z[k] + g_y_0_zzzzz_yyyyzz[k];

                g_y_0_zzzzzz_yyyzz[k] = -g_y_0_zzzzz_yyyzz[k] * cd_z[k] + g_y_0_zzzzz_yyyzzz[k];

                g_y_0_zzzzzz_yyzzz[k] = -g_y_0_zzzzz_yyzzz[k] * cd_z[k] + g_y_0_zzzzz_yyzzzz[k];

                g_y_0_zzzzzz_yzzzz[k] = -g_y_0_zzzzz_yzzzz[k] * cd_z[k] + g_y_0_zzzzz_yzzzzz[k];

                g_y_0_zzzzzz_zzzzz[k] = -g_y_0_zzzzz_zzzzz[k] * cd_z[k] + g_y_0_zzzzz_zzzzzz[k];
            }
            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 0);

            auto g_z_0_xxxxxx_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 1);

            auto g_z_0_xxxxxx_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 2);

            auto g_z_0_xxxxxx_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 3);

            auto g_z_0_xxxxxx_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 4);

            auto g_z_0_xxxxxx_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 5);

            auto g_z_0_xxxxxx_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 6);

            auto g_z_0_xxxxxx_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 7);

            auto g_z_0_xxxxxx_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 8);

            auto g_z_0_xxxxxx_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 9);

            auto g_z_0_xxxxxx_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 10);

            auto g_z_0_xxxxxx_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 11);

            auto g_z_0_xxxxxx_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 12);

            auto g_z_0_xxxxxx_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 13);

            auto g_z_0_xxxxxx_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 14);

            auto g_z_0_xxxxxx_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 15);

            auto g_z_0_xxxxxx_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 16);

            auto g_z_0_xxxxxx_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 17);

            auto g_z_0_xxxxxx_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 18);

            auto g_z_0_xxxxxx_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 19);

            auto g_z_0_xxxxxx_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxx_xxxxx, g_z_0_xxxxx_xxxxxx, g_z_0_xxxxx_xxxxxy, g_z_0_xxxxx_xxxxxz, g_z_0_xxxxx_xxxxy, g_z_0_xxxxx_xxxxyy, g_z_0_xxxxx_xxxxyz, g_z_0_xxxxx_xxxxz, g_z_0_xxxxx_xxxxzz, g_z_0_xxxxx_xxxyy, g_z_0_xxxxx_xxxyyy, g_z_0_xxxxx_xxxyyz, g_z_0_xxxxx_xxxyz, g_z_0_xxxxx_xxxyzz, g_z_0_xxxxx_xxxzz, g_z_0_xxxxx_xxxzzz, g_z_0_xxxxx_xxyyy, g_z_0_xxxxx_xxyyyy, g_z_0_xxxxx_xxyyyz, g_z_0_xxxxx_xxyyz, g_z_0_xxxxx_xxyyzz, g_z_0_xxxxx_xxyzz, g_z_0_xxxxx_xxyzzz, g_z_0_xxxxx_xxzzz, g_z_0_xxxxx_xxzzzz, g_z_0_xxxxx_xyyyy, g_z_0_xxxxx_xyyyyy, g_z_0_xxxxx_xyyyyz, g_z_0_xxxxx_xyyyz, g_z_0_xxxxx_xyyyzz, g_z_0_xxxxx_xyyzz, g_z_0_xxxxx_xyyzzz, g_z_0_xxxxx_xyzzz, g_z_0_xxxxx_xyzzzz, g_z_0_xxxxx_xzzzz, g_z_0_xxxxx_xzzzzz, g_z_0_xxxxx_yyyyy, g_z_0_xxxxx_yyyyz, g_z_0_xxxxx_yyyzz, g_z_0_xxxxx_yyzzz, g_z_0_xxxxx_yzzzz, g_z_0_xxxxx_zzzzz, g_z_0_xxxxxx_xxxxx, g_z_0_xxxxxx_xxxxy, g_z_0_xxxxxx_xxxxz, g_z_0_xxxxxx_xxxyy, g_z_0_xxxxxx_xxxyz, g_z_0_xxxxxx_xxxzz, g_z_0_xxxxxx_xxyyy, g_z_0_xxxxxx_xxyyz, g_z_0_xxxxxx_xxyzz, g_z_0_xxxxxx_xxzzz, g_z_0_xxxxxx_xyyyy, g_z_0_xxxxxx_xyyyz, g_z_0_xxxxxx_xyyzz, g_z_0_xxxxxx_xyzzz, g_z_0_xxxxxx_xzzzz, g_z_0_xxxxxx_yyyyy, g_z_0_xxxxxx_yyyyz, g_z_0_xxxxxx_yyyzz, g_z_0_xxxxxx_yyzzz, g_z_0_xxxxxx_yzzzz, g_z_0_xxxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_xxxxx[k] = -g_z_0_xxxxx_xxxxx[k] * cd_x[k] + g_z_0_xxxxx_xxxxxx[k];

                g_z_0_xxxxxx_xxxxy[k] = -g_z_0_xxxxx_xxxxy[k] * cd_x[k] + g_z_0_xxxxx_xxxxxy[k];

                g_z_0_xxxxxx_xxxxz[k] = -g_z_0_xxxxx_xxxxz[k] * cd_x[k] + g_z_0_xxxxx_xxxxxz[k];

                g_z_0_xxxxxx_xxxyy[k] = -g_z_0_xxxxx_xxxyy[k] * cd_x[k] + g_z_0_xxxxx_xxxxyy[k];

                g_z_0_xxxxxx_xxxyz[k] = -g_z_0_xxxxx_xxxyz[k] * cd_x[k] + g_z_0_xxxxx_xxxxyz[k];

                g_z_0_xxxxxx_xxxzz[k] = -g_z_0_xxxxx_xxxzz[k] * cd_x[k] + g_z_0_xxxxx_xxxxzz[k];

                g_z_0_xxxxxx_xxyyy[k] = -g_z_0_xxxxx_xxyyy[k] * cd_x[k] + g_z_0_xxxxx_xxxyyy[k];

                g_z_0_xxxxxx_xxyyz[k] = -g_z_0_xxxxx_xxyyz[k] * cd_x[k] + g_z_0_xxxxx_xxxyyz[k];

                g_z_0_xxxxxx_xxyzz[k] = -g_z_0_xxxxx_xxyzz[k] * cd_x[k] + g_z_0_xxxxx_xxxyzz[k];

                g_z_0_xxxxxx_xxzzz[k] = -g_z_0_xxxxx_xxzzz[k] * cd_x[k] + g_z_0_xxxxx_xxxzzz[k];

                g_z_0_xxxxxx_xyyyy[k] = -g_z_0_xxxxx_xyyyy[k] * cd_x[k] + g_z_0_xxxxx_xxyyyy[k];

                g_z_0_xxxxxx_xyyyz[k] = -g_z_0_xxxxx_xyyyz[k] * cd_x[k] + g_z_0_xxxxx_xxyyyz[k];

                g_z_0_xxxxxx_xyyzz[k] = -g_z_0_xxxxx_xyyzz[k] * cd_x[k] + g_z_0_xxxxx_xxyyzz[k];

                g_z_0_xxxxxx_xyzzz[k] = -g_z_0_xxxxx_xyzzz[k] * cd_x[k] + g_z_0_xxxxx_xxyzzz[k];

                g_z_0_xxxxxx_xzzzz[k] = -g_z_0_xxxxx_xzzzz[k] * cd_x[k] + g_z_0_xxxxx_xxzzzz[k];

                g_z_0_xxxxxx_yyyyy[k] = -g_z_0_xxxxx_yyyyy[k] * cd_x[k] + g_z_0_xxxxx_xyyyyy[k];

                g_z_0_xxxxxx_yyyyz[k] = -g_z_0_xxxxx_yyyyz[k] * cd_x[k] + g_z_0_xxxxx_xyyyyz[k];

                g_z_0_xxxxxx_yyyzz[k] = -g_z_0_xxxxx_yyyzz[k] * cd_x[k] + g_z_0_xxxxx_xyyyzz[k];

                g_z_0_xxxxxx_yyzzz[k] = -g_z_0_xxxxx_yyzzz[k] * cd_x[k] + g_z_0_xxxxx_xyyzzz[k];

                g_z_0_xxxxxx_yzzzz[k] = -g_z_0_xxxxx_yzzzz[k] * cd_x[k] + g_z_0_xxxxx_xyzzzz[k];

                g_z_0_xxxxxx_zzzzz[k] = -g_z_0_xxxxx_zzzzz[k] * cd_x[k] + g_z_0_xxxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 21);

            auto g_z_0_xxxxxy_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 22);

            auto g_z_0_xxxxxy_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 23);

            auto g_z_0_xxxxxy_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 24);

            auto g_z_0_xxxxxy_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 25);

            auto g_z_0_xxxxxy_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 26);

            auto g_z_0_xxxxxy_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 27);

            auto g_z_0_xxxxxy_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 28);

            auto g_z_0_xxxxxy_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 29);

            auto g_z_0_xxxxxy_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 30);

            auto g_z_0_xxxxxy_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 31);

            auto g_z_0_xxxxxy_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 32);

            auto g_z_0_xxxxxy_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 33);

            auto g_z_0_xxxxxy_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 34);

            auto g_z_0_xxxxxy_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 35);

            auto g_z_0_xxxxxy_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 36);

            auto g_z_0_xxxxxy_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 37);

            auto g_z_0_xxxxxy_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 38);

            auto g_z_0_xxxxxy_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 39);

            auto g_z_0_xxxxxy_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 40);

            auto g_z_0_xxxxxy_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxy_xxxxx, g_z_0_xxxxxy_xxxxy, g_z_0_xxxxxy_xxxxz, g_z_0_xxxxxy_xxxyy, g_z_0_xxxxxy_xxxyz, g_z_0_xxxxxy_xxxzz, g_z_0_xxxxxy_xxyyy, g_z_0_xxxxxy_xxyyz, g_z_0_xxxxxy_xxyzz, g_z_0_xxxxxy_xxzzz, g_z_0_xxxxxy_xyyyy, g_z_0_xxxxxy_xyyyz, g_z_0_xxxxxy_xyyzz, g_z_0_xxxxxy_xyzzz, g_z_0_xxxxxy_xzzzz, g_z_0_xxxxxy_yyyyy, g_z_0_xxxxxy_yyyyz, g_z_0_xxxxxy_yyyzz, g_z_0_xxxxxy_yyzzz, g_z_0_xxxxxy_yzzzz, g_z_0_xxxxxy_zzzzz, g_z_0_xxxxy_xxxxx, g_z_0_xxxxy_xxxxxx, g_z_0_xxxxy_xxxxxy, g_z_0_xxxxy_xxxxxz, g_z_0_xxxxy_xxxxy, g_z_0_xxxxy_xxxxyy, g_z_0_xxxxy_xxxxyz, g_z_0_xxxxy_xxxxz, g_z_0_xxxxy_xxxxzz, g_z_0_xxxxy_xxxyy, g_z_0_xxxxy_xxxyyy, g_z_0_xxxxy_xxxyyz, g_z_0_xxxxy_xxxyz, g_z_0_xxxxy_xxxyzz, g_z_0_xxxxy_xxxzz, g_z_0_xxxxy_xxxzzz, g_z_0_xxxxy_xxyyy, g_z_0_xxxxy_xxyyyy, g_z_0_xxxxy_xxyyyz, g_z_0_xxxxy_xxyyz, g_z_0_xxxxy_xxyyzz, g_z_0_xxxxy_xxyzz, g_z_0_xxxxy_xxyzzz, g_z_0_xxxxy_xxzzz, g_z_0_xxxxy_xxzzzz, g_z_0_xxxxy_xyyyy, g_z_0_xxxxy_xyyyyy, g_z_0_xxxxy_xyyyyz, g_z_0_xxxxy_xyyyz, g_z_0_xxxxy_xyyyzz, g_z_0_xxxxy_xyyzz, g_z_0_xxxxy_xyyzzz, g_z_0_xxxxy_xyzzz, g_z_0_xxxxy_xyzzzz, g_z_0_xxxxy_xzzzz, g_z_0_xxxxy_xzzzzz, g_z_0_xxxxy_yyyyy, g_z_0_xxxxy_yyyyz, g_z_0_xxxxy_yyyzz, g_z_0_xxxxy_yyzzz, g_z_0_xxxxy_yzzzz, g_z_0_xxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_xxxxx[k] = -g_z_0_xxxxy_xxxxx[k] * cd_x[k] + g_z_0_xxxxy_xxxxxx[k];

                g_z_0_xxxxxy_xxxxy[k] = -g_z_0_xxxxy_xxxxy[k] * cd_x[k] + g_z_0_xxxxy_xxxxxy[k];

                g_z_0_xxxxxy_xxxxz[k] = -g_z_0_xxxxy_xxxxz[k] * cd_x[k] + g_z_0_xxxxy_xxxxxz[k];

                g_z_0_xxxxxy_xxxyy[k] = -g_z_0_xxxxy_xxxyy[k] * cd_x[k] + g_z_0_xxxxy_xxxxyy[k];

                g_z_0_xxxxxy_xxxyz[k] = -g_z_0_xxxxy_xxxyz[k] * cd_x[k] + g_z_0_xxxxy_xxxxyz[k];

                g_z_0_xxxxxy_xxxzz[k] = -g_z_0_xxxxy_xxxzz[k] * cd_x[k] + g_z_0_xxxxy_xxxxzz[k];

                g_z_0_xxxxxy_xxyyy[k] = -g_z_0_xxxxy_xxyyy[k] * cd_x[k] + g_z_0_xxxxy_xxxyyy[k];

                g_z_0_xxxxxy_xxyyz[k] = -g_z_0_xxxxy_xxyyz[k] * cd_x[k] + g_z_0_xxxxy_xxxyyz[k];

                g_z_0_xxxxxy_xxyzz[k] = -g_z_0_xxxxy_xxyzz[k] * cd_x[k] + g_z_0_xxxxy_xxxyzz[k];

                g_z_0_xxxxxy_xxzzz[k] = -g_z_0_xxxxy_xxzzz[k] * cd_x[k] + g_z_0_xxxxy_xxxzzz[k];

                g_z_0_xxxxxy_xyyyy[k] = -g_z_0_xxxxy_xyyyy[k] * cd_x[k] + g_z_0_xxxxy_xxyyyy[k];

                g_z_0_xxxxxy_xyyyz[k] = -g_z_0_xxxxy_xyyyz[k] * cd_x[k] + g_z_0_xxxxy_xxyyyz[k];

                g_z_0_xxxxxy_xyyzz[k] = -g_z_0_xxxxy_xyyzz[k] * cd_x[k] + g_z_0_xxxxy_xxyyzz[k];

                g_z_0_xxxxxy_xyzzz[k] = -g_z_0_xxxxy_xyzzz[k] * cd_x[k] + g_z_0_xxxxy_xxyzzz[k];

                g_z_0_xxxxxy_xzzzz[k] = -g_z_0_xxxxy_xzzzz[k] * cd_x[k] + g_z_0_xxxxy_xxzzzz[k];

                g_z_0_xxxxxy_yyyyy[k] = -g_z_0_xxxxy_yyyyy[k] * cd_x[k] + g_z_0_xxxxy_xyyyyy[k];

                g_z_0_xxxxxy_yyyyz[k] = -g_z_0_xxxxy_yyyyz[k] * cd_x[k] + g_z_0_xxxxy_xyyyyz[k];

                g_z_0_xxxxxy_yyyzz[k] = -g_z_0_xxxxy_yyyzz[k] * cd_x[k] + g_z_0_xxxxy_xyyyzz[k];

                g_z_0_xxxxxy_yyzzz[k] = -g_z_0_xxxxy_yyzzz[k] * cd_x[k] + g_z_0_xxxxy_xyyzzz[k];

                g_z_0_xxxxxy_yzzzz[k] = -g_z_0_xxxxy_yzzzz[k] * cd_x[k] + g_z_0_xxxxy_xyzzzz[k];

                g_z_0_xxxxxy_zzzzz[k] = -g_z_0_xxxxy_zzzzz[k] * cd_x[k] + g_z_0_xxxxy_xzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 42);

            auto g_z_0_xxxxxz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 43);

            auto g_z_0_xxxxxz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 44);

            auto g_z_0_xxxxxz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 45);

            auto g_z_0_xxxxxz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 46);

            auto g_z_0_xxxxxz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 47);

            auto g_z_0_xxxxxz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 48);

            auto g_z_0_xxxxxz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 49);

            auto g_z_0_xxxxxz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 50);

            auto g_z_0_xxxxxz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 51);

            auto g_z_0_xxxxxz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 52);

            auto g_z_0_xxxxxz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 53);

            auto g_z_0_xxxxxz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 54);

            auto g_z_0_xxxxxz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 55);

            auto g_z_0_xxxxxz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 56);

            auto g_z_0_xxxxxz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 57);

            auto g_z_0_xxxxxz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 58);

            auto g_z_0_xxxxxz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 59);

            auto g_z_0_xxxxxz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 60);

            auto g_z_0_xxxxxz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 61);

            auto g_z_0_xxxxxz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxz_xxxxx, g_z_0_xxxxxz_xxxxy, g_z_0_xxxxxz_xxxxz, g_z_0_xxxxxz_xxxyy, g_z_0_xxxxxz_xxxyz, g_z_0_xxxxxz_xxxzz, g_z_0_xxxxxz_xxyyy, g_z_0_xxxxxz_xxyyz, g_z_0_xxxxxz_xxyzz, g_z_0_xxxxxz_xxzzz, g_z_0_xxxxxz_xyyyy, g_z_0_xxxxxz_xyyyz, g_z_0_xxxxxz_xyyzz, g_z_0_xxxxxz_xyzzz, g_z_0_xxxxxz_xzzzz, g_z_0_xxxxxz_yyyyy, g_z_0_xxxxxz_yyyyz, g_z_0_xxxxxz_yyyzz, g_z_0_xxxxxz_yyzzz, g_z_0_xxxxxz_yzzzz, g_z_0_xxxxxz_zzzzz, g_z_0_xxxxz_xxxxx, g_z_0_xxxxz_xxxxxx, g_z_0_xxxxz_xxxxxy, g_z_0_xxxxz_xxxxxz, g_z_0_xxxxz_xxxxy, g_z_0_xxxxz_xxxxyy, g_z_0_xxxxz_xxxxyz, g_z_0_xxxxz_xxxxz, g_z_0_xxxxz_xxxxzz, g_z_0_xxxxz_xxxyy, g_z_0_xxxxz_xxxyyy, g_z_0_xxxxz_xxxyyz, g_z_0_xxxxz_xxxyz, g_z_0_xxxxz_xxxyzz, g_z_0_xxxxz_xxxzz, g_z_0_xxxxz_xxxzzz, g_z_0_xxxxz_xxyyy, g_z_0_xxxxz_xxyyyy, g_z_0_xxxxz_xxyyyz, g_z_0_xxxxz_xxyyz, g_z_0_xxxxz_xxyyzz, g_z_0_xxxxz_xxyzz, g_z_0_xxxxz_xxyzzz, g_z_0_xxxxz_xxzzz, g_z_0_xxxxz_xxzzzz, g_z_0_xxxxz_xyyyy, g_z_0_xxxxz_xyyyyy, g_z_0_xxxxz_xyyyyz, g_z_0_xxxxz_xyyyz, g_z_0_xxxxz_xyyyzz, g_z_0_xxxxz_xyyzz, g_z_0_xxxxz_xyyzzz, g_z_0_xxxxz_xyzzz, g_z_0_xxxxz_xyzzzz, g_z_0_xxxxz_xzzzz, g_z_0_xxxxz_xzzzzz, g_z_0_xxxxz_yyyyy, g_z_0_xxxxz_yyyyz, g_z_0_xxxxz_yyyzz, g_z_0_xxxxz_yyzzz, g_z_0_xxxxz_yzzzz, g_z_0_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_xxxxx[k] = -g_z_0_xxxxz_xxxxx[k] * cd_x[k] + g_z_0_xxxxz_xxxxxx[k];

                g_z_0_xxxxxz_xxxxy[k] = -g_z_0_xxxxz_xxxxy[k] * cd_x[k] + g_z_0_xxxxz_xxxxxy[k];

                g_z_0_xxxxxz_xxxxz[k] = -g_z_0_xxxxz_xxxxz[k] * cd_x[k] + g_z_0_xxxxz_xxxxxz[k];

                g_z_0_xxxxxz_xxxyy[k] = -g_z_0_xxxxz_xxxyy[k] * cd_x[k] + g_z_0_xxxxz_xxxxyy[k];

                g_z_0_xxxxxz_xxxyz[k] = -g_z_0_xxxxz_xxxyz[k] * cd_x[k] + g_z_0_xxxxz_xxxxyz[k];

                g_z_0_xxxxxz_xxxzz[k] = -g_z_0_xxxxz_xxxzz[k] * cd_x[k] + g_z_0_xxxxz_xxxxzz[k];

                g_z_0_xxxxxz_xxyyy[k] = -g_z_0_xxxxz_xxyyy[k] * cd_x[k] + g_z_0_xxxxz_xxxyyy[k];

                g_z_0_xxxxxz_xxyyz[k] = -g_z_0_xxxxz_xxyyz[k] * cd_x[k] + g_z_0_xxxxz_xxxyyz[k];

                g_z_0_xxxxxz_xxyzz[k] = -g_z_0_xxxxz_xxyzz[k] * cd_x[k] + g_z_0_xxxxz_xxxyzz[k];

                g_z_0_xxxxxz_xxzzz[k] = -g_z_0_xxxxz_xxzzz[k] * cd_x[k] + g_z_0_xxxxz_xxxzzz[k];

                g_z_0_xxxxxz_xyyyy[k] = -g_z_0_xxxxz_xyyyy[k] * cd_x[k] + g_z_0_xxxxz_xxyyyy[k];

                g_z_0_xxxxxz_xyyyz[k] = -g_z_0_xxxxz_xyyyz[k] * cd_x[k] + g_z_0_xxxxz_xxyyyz[k];

                g_z_0_xxxxxz_xyyzz[k] = -g_z_0_xxxxz_xyyzz[k] * cd_x[k] + g_z_0_xxxxz_xxyyzz[k];

                g_z_0_xxxxxz_xyzzz[k] = -g_z_0_xxxxz_xyzzz[k] * cd_x[k] + g_z_0_xxxxz_xxyzzz[k];

                g_z_0_xxxxxz_xzzzz[k] = -g_z_0_xxxxz_xzzzz[k] * cd_x[k] + g_z_0_xxxxz_xxzzzz[k];

                g_z_0_xxxxxz_yyyyy[k] = -g_z_0_xxxxz_yyyyy[k] * cd_x[k] + g_z_0_xxxxz_xyyyyy[k];

                g_z_0_xxxxxz_yyyyz[k] = -g_z_0_xxxxz_yyyyz[k] * cd_x[k] + g_z_0_xxxxz_xyyyyz[k];

                g_z_0_xxxxxz_yyyzz[k] = -g_z_0_xxxxz_yyyzz[k] * cd_x[k] + g_z_0_xxxxz_xyyyzz[k];

                g_z_0_xxxxxz_yyzzz[k] = -g_z_0_xxxxz_yyzzz[k] * cd_x[k] + g_z_0_xxxxz_xyyzzz[k];

                g_z_0_xxxxxz_yzzzz[k] = -g_z_0_xxxxz_yzzzz[k] * cd_x[k] + g_z_0_xxxxz_xyzzzz[k];

                g_z_0_xxxxxz_zzzzz[k] = -g_z_0_xxxxz_zzzzz[k] * cd_x[k] + g_z_0_xxxxz_xzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 63);

            auto g_z_0_xxxxyy_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 64);

            auto g_z_0_xxxxyy_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 65);

            auto g_z_0_xxxxyy_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 66);

            auto g_z_0_xxxxyy_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 67);

            auto g_z_0_xxxxyy_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 68);

            auto g_z_0_xxxxyy_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 69);

            auto g_z_0_xxxxyy_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 70);

            auto g_z_0_xxxxyy_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 71);

            auto g_z_0_xxxxyy_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 72);

            auto g_z_0_xxxxyy_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 73);

            auto g_z_0_xxxxyy_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 74);

            auto g_z_0_xxxxyy_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 75);

            auto g_z_0_xxxxyy_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 76);

            auto g_z_0_xxxxyy_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 77);

            auto g_z_0_xxxxyy_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 78);

            auto g_z_0_xxxxyy_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 79);

            auto g_z_0_xxxxyy_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 80);

            auto g_z_0_xxxxyy_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 81);

            auto g_z_0_xxxxyy_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 82);

            auto g_z_0_xxxxyy_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyy_xxxxx, g_z_0_xxxxyy_xxxxy, g_z_0_xxxxyy_xxxxz, g_z_0_xxxxyy_xxxyy, g_z_0_xxxxyy_xxxyz, g_z_0_xxxxyy_xxxzz, g_z_0_xxxxyy_xxyyy, g_z_0_xxxxyy_xxyyz, g_z_0_xxxxyy_xxyzz, g_z_0_xxxxyy_xxzzz, g_z_0_xxxxyy_xyyyy, g_z_0_xxxxyy_xyyyz, g_z_0_xxxxyy_xyyzz, g_z_0_xxxxyy_xyzzz, g_z_0_xxxxyy_xzzzz, g_z_0_xxxxyy_yyyyy, g_z_0_xxxxyy_yyyyz, g_z_0_xxxxyy_yyyzz, g_z_0_xxxxyy_yyzzz, g_z_0_xxxxyy_yzzzz, g_z_0_xxxxyy_zzzzz, g_z_0_xxxyy_xxxxx, g_z_0_xxxyy_xxxxxx, g_z_0_xxxyy_xxxxxy, g_z_0_xxxyy_xxxxxz, g_z_0_xxxyy_xxxxy, g_z_0_xxxyy_xxxxyy, g_z_0_xxxyy_xxxxyz, g_z_0_xxxyy_xxxxz, g_z_0_xxxyy_xxxxzz, g_z_0_xxxyy_xxxyy, g_z_0_xxxyy_xxxyyy, g_z_0_xxxyy_xxxyyz, g_z_0_xxxyy_xxxyz, g_z_0_xxxyy_xxxyzz, g_z_0_xxxyy_xxxzz, g_z_0_xxxyy_xxxzzz, g_z_0_xxxyy_xxyyy, g_z_0_xxxyy_xxyyyy, g_z_0_xxxyy_xxyyyz, g_z_0_xxxyy_xxyyz, g_z_0_xxxyy_xxyyzz, g_z_0_xxxyy_xxyzz, g_z_0_xxxyy_xxyzzz, g_z_0_xxxyy_xxzzz, g_z_0_xxxyy_xxzzzz, g_z_0_xxxyy_xyyyy, g_z_0_xxxyy_xyyyyy, g_z_0_xxxyy_xyyyyz, g_z_0_xxxyy_xyyyz, g_z_0_xxxyy_xyyyzz, g_z_0_xxxyy_xyyzz, g_z_0_xxxyy_xyyzzz, g_z_0_xxxyy_xyzzz, g_z_0_xxxyy_xyzzzz, g_z_0_xxxyy_xzzzz, g_z_0_xxxyy_xzzzzz, g_z_0_xxxyy_yyyyy, g_z_0_xxxyy_yyyyz, g_z_0_xxxyy_yyyzz, g_z_0_xxxyy_yyzzz, g_z_0_xxxyy_yzzzz, g_z_0_xxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_xxxxx[k] = -g_z_0_xxxyy_xxxxx[k] * cd_x[k] + g_z_0_xxxyy_xxxxxx[k];

                g_z_0_xxxxyy_xxxxy[k] = -g_z_0_xxxyy_xxxxy[k] * cd_x[k] + g_z_0_xxxyy_xxxxxy[k];

                g_z_0_xxxxyy_xxxxz[k] = -g_z_0_xxxyy_xxxxz[k] * cd_x[k] + g_z_0_xxxyy_xxxxxz[k];

                g_z_0_xxxxyy_xxxyy[k] = -g_z_0_xxxyy_xxxyy[k] * cd_x[k] + g_z_0_xxxyy_xxxxyy[k];

                g_z_0_xxxxyy_xxxyz[k] = -g_z_0_xxxyy_xxxyz[k] * cd_x[k] + g_z_0_xxxyy_xxxxyz[k];

                g_z_0_xxxxyy_xxxzz[k] = -g_z_0_xxxyy_xxxzz[k] * cd_x[k] + g_z_0_xxxyy_xxxxzz[k];

                g_z_0_xxxxyy_xxyyy[k] = -g_z_0_xxxyy_xxyyy[k] * cd_x[k] + g_z_0_xxxyy_xxxyyy[k];

                g_z_0_xxxxyy_xxyyz[k] = -g_z_0_xxxyy_xxyyz[k] * cd_x[k] + g_z_0_xxxyy_xxxyyz[k];

                g_z_0_xxxxyy_xxyzz[k] = -g_z_0_xxxyy_xxyzz[k] * cd_x[k] + g_z_0_xxxyy_xxxyzz[k];

                g_z_0_xxxxyy_xxzzz[k] = -g_z_0_xxxyy_xxzzz[k] * cd_x[k] + g_z_0_xxxyy_xxxzzz[k];

                g_z_0_xxxxyy_xyyyy[k] = -g_z_0_xxxyy_xyyyy[k] * cd_x[k] + g_z_0_xxxyy_xxyyyy[k];

                g_z_0_xxxxyy_xyyyz[k] = -g_z_0_xxxyy_xyyyz[k] * cd_x[k] + g_z_0_xxxyy_xxyyyz[k];

                g_z_0_xxxxyy_xyyzz[k] = -g_z_0_xxxyy_xyyzz[k] * cd_x[k] + g_z_0_xxxyy_xxyyzz[k];

                g_z_0_xxxxyy_xyzzz[k] = -g_z_0_xxxyy_xyzzz[k] * cd_x[k] + g_z_0_xxxyy_xxyzzz[k];

                g_z_0_xxxxyy_xzzzz[k] = -g_z_0_xxxyy_xzzzz[k] * cd_x[k] + g_z_0_xxxyy_xxzzzz[k];

                g_z_0_xxxxyy_yyyyy[k] = -g_z_0_xxxyy_yyyyy[k] * cd_x[k] + g_z_0_xxxyy_xyyyyy[k];

                g_z_0_xxxxyy_yyyyz[k] = -g_z_0_xxxyy_yyyyz[k] * cd_x[k] + g_z_0_xxxyy_xyyyyz[k];

                g_z_0_xxxxyy_yyyzz[k] = -g_z_0_xxxyy_yyyzz[k] * cd_x[k] + g_z_0_xxxyy_xyyyzz[k];

                g_z_0_xxxxyy_yyzzz[k] = -g_z_0_xxxyy_yyzzz[k] * cd_x[k] + g_z_0_xxxyy_xyyzzz[k];

                g_z_0_xxxxyy_yzzzz[k] = -g_z_0_xxxyy_yzzzz[k] * cd_x[k] + g_z_0_xxxyy_xyzzzz[k];

                g_z_0_xxxxyy_zzzzz[k] = -g_z_0_xxxyy_zzzzz[k] * cd_x[k] + g_z_0_xxxyy_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 84);

            auto g_z_0_xxxxyz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 85);

            auto g_z_0_xxxxyz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 86);

            auto g_z_0_xxxxyz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 87);

            auto g_z_0_xxxxyz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 88);

            auto g_z_0_xxxxyz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 89);

            auto g_z_0_xxxxyz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 90);

            auto g_z_0_xxxxyz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 91);

            auto g_z_0_xxxxyz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 92);

            auto g_z_0_xxxxyz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 93);

            auto g_z_0_xxxxyz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 94);

            auto g_z_0_xxxxyz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 95);

            auto g_z_0_xxxxyz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 96);

            auto g_z_0_xxxxyz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 97);

            auto g_z_0_xxxxyz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 98);

            auto g_z_0_xxxxyz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 99);

            auto g_z_0_xxxxyz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 100);

            auto g_z_0_xxxxyz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 101);

            auto g_z_0_xxxxyz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 102);

            auto g_z_0_xxxxyz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 103);

            auto g_z_0_xxxxyz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyz_xxxxx, g_z_0_xxxxyz_xxxxy, g_z_0_xxxxyz_xxxxz, g_z_0_xxxxyz_xxxyy, g_z_0_xxxxyz_xxxyz, g_z_0_xxxxyz_xxxzz, g_z_0_xxxxyz_xxyyy, g_z_0_xxxxyz_xxyyz, g_z_0_xxxxyz_xxyzz, g_z_0_xxxxyz_xxzzz, g_z_0_xxxxyz_xyyyy, g_z_0_xxxxyz_xyyyz, g_z_0_xxxxyz_xyyzz, g_z_0_xxxxyz_xyzzz, g_z_0_xxxxyz_xzzzz, g_z_0_xxxxyz_yyyyy, g_z_0_xxxxyz_yyyyz, g_z_0_xxxxyz_yyyzz, g_z_0_xxxxyz_yyzzz, g_z_0_xxxxyz_yzzzz, g_z_0_xxxxyz_zzzzz, g_z_0_xxxyz_xxxxx, g_z_0_xxxyz_xxxxxx, g_z_0_xxxyz_xxxxxy, g_z_0_xxxyz_xxxxxz, g_z_0_xxxyz_xxxxy, g_z_0_xxxyz_xxxxyy, g_z_0_xxxyz_xxxxyz, g_z_0_xxxyz_xxxxz, g_z_0_xxxyz_xxxxzz, g_z_0_xxxyz_xxxyy, g_z_0_xxxyz_xxxyyy, g_z_0_xxxyz_xxxyyz, g_z_0_xxxyz_xxxyz, g_z_0_xxxyz_xxxyzz, g_z_0_xxxyz_xxxzz, g_z_0_xxxyz_xxxzzz, g_z_0_xxxyz_xxyyy, g_z_0_xxxyz_xxyyyy, g_z_0_xxxyz_xxyyyz, g_z_0_xxxyz_xxyyz, g_z_0_xxxyz_xxyyzz, g_z_0_xxxyz_xxyzz, g_z_0_xxxyz_xxyzzz, g_z_0_xxxyz_xxzzz, g_z_0_xxxyz_xxzzzz, g_z_0_xxxyz_xyyyy, g_z_0_xxxyz_xyyyyy, g_z_0_xxxyz_xyyyyz, g_z_0_xxxyz_xyyyz, g_z_0_xxxyz_xyyyzz, g_z_0_xxxyz_xyyzz, g_z_0_xxxyz_xyyzzz, g_z_0_xxxyz_xyzzz, g_z_0_xxxyz_xyzzzz, g_z_0_xxxyz_xzzzz, g_z_0_xxxyz_xzzzzz, g_z_0_xxxyz_yyyyy, g_z_0_xxxyz_yyyyz, g_z_0_xxxyz_yyyzz, g_z_0_xxxyz_yyzzz, g_z_0_xxxyz_yzzzz, g_z_0_xxxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_xxxxx[k] = -g_z_0_xxxyz_xxxxx[k] * cd_x[k] + g_z_0_xxxyz_xxxxxx[k];

                g_z_0_xxxxyz_xxxxy[k] = -g_z_0_xxxyz_xxxxy[k] * cd_x[k] + g_z_0_xxxyz_xxxxxy[k];

                g_z_0_xxxxyz_xxxxz[k] = -g_z_0_xxxyz_xxxxz[k] * cd_x[k] + g_z_0_xxxyz_xxxxxz[k];

                g_z_0_xxxxyz_xxxyy[k] = -g_z_0_xxxyz_xxxyy[k] * cd_x[k] + g_z_0_xxxyz_xxxxyy[k];

                g_z_0_xxxxyz_xxxyz[k] = -g_z_0_xxxyz_xxxyz[k] * cd_x[k] + g_z_0_xxxyz_xxxxyz[k];

                g_z_0_xxxxyz_xxxzz[k] = -g_z_0_xxxyz_xxxzz[k] * cd_x[k] + g_z_0_xxxyz_xxxxzz[k];

                g_z_0_xxxxyz_xxyyy[k] = -g_z_0_xxxyz_xxyyy[k] * cd_x[k] + g_z_0_xxxyz_xxxyyy[k];

                g_z_0_xxxxyz_xxyyz[k] = -g_z_0_xxxyz_xxyyz[k] * cd_x[k] + g_z_0_xxxyz_xxxyyz[k];

                g_z_0_xxxxyz_xxyzz[k] = -g_z_0_xxxyz_xxyzz[k] * cd_x[k] + g_z_0_xxxyz_xxxyzz[k];

                g_z_0_xxxxyz_xxzzz[k] = -g_z_0_xxxyz_xxzzz[k] * cd_x[k] + g_z_0_xxxyz_xxxzzz[k];

                g_z_0_xxxxyz_xyyyy[k] = -g_z_0_xxxyz_xyyyy[k] * cd_x[k] + g_z_0_xxxyz_xxyyyy[k];

                g_z_0_xxxxyz_xyyyz[k] = -g_z_0_xxxyz_xyyyz[k] * cd_x[k] + g_z_0_xxxyz_xxyyyz[k];

                g_z_0_xxxxyz_xyyzz[k] = -g_z_0_xxxyz_xyyzz[k] * cd_x[k] + g_z_0_xxxyz_xxyyzz[k];

                g_z_0_xxxxyz_xyzzz[k] = -g_z_0_xxxyz_xyzzz[k] * cd_x[k] + g_z_0_xxxyz_xxyzzz[k];

                g_z_0_xxxxyz_xzzzz[k] = -g_z_0_xxxyz_xzzzz[k] * cd_x[k] + g_z_0_xxxyz_xxzzzz[k];

                g_z_0_xxxxyz_yyyyy[k] = -g_z_0_xxxyz_yyyyy[k] * cd_x[k] + g_z_0_xxxyz_xyyyyy[k];

                g_z_0_xxxxyz_yyyyz[k] = -g_z_0_xxxyz_yyyyz[k] * cd_x[k] + g_z_0_xxxyz_xyyyyz[k];

                g_z_0_xxxxyz_yyyzz[k] = -g_z_0_xxxyz_yyyzz[k] * cd_x[k] + g_z_0_xxxyz_xyyyzz[k];

                g_z_0_xxxxyz_yyzzz[k] = -g_z_0_xxxyz_yyzzz[k] * cd_x[k] + g_z_0_xxxyz_xyyzzz[k];

                g_z_0_xxxxyz_yzzzz[k] = -g_z_0_xxxyz_yzzzz[k] * cd_x[k] + g_z_0_xxxyz_xyzzzz[k];

                g_z_0_xxxxyz_zzzzz[k] = -g_z_0_xxxyz_zzzzz[k] * cd_x[k] + g_z_0_xxxyz_xzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 105);

            auto g_z_0_xxxxzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 106);

            auto g_z_0_xxxxzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 107);

            auto g_z_0_xxxxzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 108);

            auto g_z_0_xxxxzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 109);

            auto g_z_0_xxxxzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 110);

            auto g_z_0_xxxxzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 111);

            auto g_z_0_xxxxzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 112);

            auto g_z_0_xxxxzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 113);

            auto g_z_0_xxxxzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 114);

            auto g_z_0_xxxxzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 115);

            auto g_z_0_xxxxzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 116);

            auto g_z_0_xxxxzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 117);

            auto g_z_0_xxxxzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 118);

            auto g_z_0_xxxxzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 119);

            auto g_z_0_xxxxzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 120);

            auto g_z_0_xxxxzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 121);

            auto g_z_0_xxxxzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 122);

            auto g_z_0_xxxxzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 123);

            auto g_z_0_xxxxzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 124);

            auto g_z_0_xxxxzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxzz_xxxxx, g_z_0_xxxxzz_xxxxy, g_z_0_xxxxzz_xxxxz, g_z_0_xxxxzz_xxxyy, g_z_0_xxxxzz_xxxyz, g_z_0_xxxxzz_xxxzz, g_z_0_xxxxzz_xxyyy, g_z_0_xxxxzz_xxyyz, g_z_0_xxxxzz_xxyzz, g_z_0_xxxxzz_xxzzz, g_z_0_xxxxzz_xyyyy, g_z_0_xxxxzz_xyyyz, g_z_0_xxxxzz_xyyzz, g_z_0_xxxxzz_xyzzz, g_z_0_xxxxzz_xzzzz, g_z_0_xxxxzz_yyyyy, g_z_0_xxxxzz_yyyyz, g_z_0_xxxxzz_yyyzz, g_z_0_xxxxzz_yyzzz, g_z_0_xxxxzz_yzzzz, g_z_0_xxxxzz_zzzzz, g_z_0_xxxzz_xxxxx, g_z_0_xxxzz_xxxxxx, g_z_0_xxxzz_xxxxxy, g_z_0_xxxzz_xxxxxz, g_z_0_xxxzz_xxxxy, g_z_0_xxxzz_xxxxyy, g_z_0_xxxzz_xxxxyz, g_z_0_xxxzz_xxxxz, g_z_0_xxxzz_xxxxzz, g_z_0_xxxzz_xxxyy, g_z_0_xxxzz_xxxyyy, g_z_0_xxxzz_xxxyyz, g_z_0_xxxzz_xxxyz, g_z_0_xxxzz_xxxyzz, g_z_0_xxxzz_xxxzz, g_z_0_xxxzz_xxxzzz, g_z_0_xxxzz_xxyyy, g_z_0_xxxzz_xxyyyy, g_z_0_xxxzz_xxyyyz, g_z_0_xxxzz_xxyyz, g_z_0_xxxzz_xxyyzz, g_z_0_xxxzz_xxyzz, g_z_0_xxxzz_xxyzzz, g_z_0_xxxzz_xxzzz, g_z_0_xxxzz_xxzzzz, g_z_0_xxxzz_xyyyy, g_z_0_xxxzz_xyyyyy, g_z_0_xxxzz_xyyyyz, g_z_0_xxxzz_xyyyz, g_z_0_xxxzz_xyyyzz, g_z_0_xxxzz_xyyzz, g_z_0_xxxzz_xyyzzz, g_z_0_xxxzz_xyzzz, g_z_0_xxxzz_xyzzzz, g_z_0_xxxzz_xzzzz, g_z_0_xxxzz_xzzzzz, g_z_0_xxxzz_yyyyy, g_z_0_xxxzz_yyyyz, g_z_0_xxxzz_yyyzz, g_z_0_xxxzz_yyzzz, g_z_0_xxxzz_yzzzz, g_z_0_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_xxxxx[k] = -g_z_0_xxxzz_xxxxx[k] * cd_x[k] + g_z_0_xxxzz_xxxxxx[k];

                g_z_0_xxxxzz_xxxxy[k] = -g_z_0_xxxzz_xxxxy[k] * cd_x[k] + g_z_0_xxxzz_xxxxxy[k];

                g_z_0_xxxxzz_xxxxz[k] = -g_z_0_xxxzz_xxxxz[k] * cd_x[k] + g_z_0_xxxzz_xxxxxz[k];

                g_z_0_xxxxzz_xxxyy[k] = -g_z_0_xxxzz_xxxyy[k] * cd_x[k] + g_z_0_xxxzz_xxxxyy[k];

                g_z_0_xxxxzz_xxxyz[k] = -g_z_0_xxxzz_xxxyz[k] * cd_x[k] + g_z_0_xxxzz_xxxxyz[k];

                g_z_0_xxxxzz_xxxzz[k] = -g_z_0_xxxzz_xxxzz[k] * cd_x[k] + g_z_0_xxxzz_xxxxzz[k];

                g_z_0_xxxxzz_xxyyy[k] = -g_z_0_xxxzz_xxyyy[k] * cd_x[k] + g_z_0_xxxzz_xxxyyy[k];

                g_z_0_xxxxzz_xxyyz[k] = -g_z_0_xxxzz_xxyyz[k] * cd_x[k] + g_z_0_xxxzz_xxxyyz[k];

                g_z_0_xxxxzz_xxyzz[k] = -g_z_0_xxxzz_xxyzz[k] * cd_x[k] + g_z_0_xxxzz_xxxyzz[k];

                g_z_0_xxxxzz_xxzzz[k] = -g_z_0_xxxzz_xxzzz[k] * cd_x[k] + g_z_0_xxxzz_xxxzzz[k];

                g_z_0_xxxxzz_xyyyy[k] = -g_z_0_xxxzz_xyyyy[k] * cd_x[k] + g_z_0_xxxzz_xxyyyy[k];

                g_z_0_xxxxzz_xyyyz[k] = -g_z_0_xxxzz_xyyyz[k] * cd_x[k] + g_z_0_xxxzz_xxyyyz[k];

                g_z_0_xxxxzz_xyyzz[k] = -g_z_0_xxxzz_xyyzz[k] * cd_x[k] + g_z_0_xxxzz_xxyyzz[k];

                g_z_0_xxxxzz_xyzzz[k] = -g_z_0_xxxzz_xyzzz[k] * cd_x[k] + g_z_0_xxxzz_xxyzzz[k];

                g_z_0_xxxxzz_xzzzz[k] = -g_z_0_xxxzz_xzzzz[k] * cd_x[k] + g_z_0_xxxzz_xxzzzz[k];

                g_z_0_xxxxzz_yyyyy[k] = -g_z_0_xxxzz_yyyyy[k] * cd_x[k] + g_z_0_xxxzz_xyyyyy[k];

                g_z_0_xxxxzz_yyyyz[k] = -g_z_0_xxxzz_yyyyz[k] * cd_x[k] + g_z_0_xxxzz_xyyyyz[k];

                g_z_0_xxxxzz_yyyzz[k] = -g_z_0_xxxzz_yyyzz[k] * cd_x[k] + g_z_0_xxxzz_xyyyzz[k];

                g_z_0_xxxxzz_yyzzz[k] = -g_z_0_xxxzz_yyzzz[k] * cd_x[k] + g_z_0_xxxzz_xyyzzz[k];

                g_z_0_xxxxzz_yzzzz[k] = -g_z_0_xxxzz_yzzzz[k] * cd_x[k] + g_z_0_xxxzz_xyzzzz[k];

                g_z_0_xxxxzz_zzzzz[k] = -g_z_0_xxxzz_zzzzz[k] * cd_x[k] + g_z_0_xxxzz_xzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 126);

            auto g_z_0_xxxyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 127);

            auto g_z_0_xxxyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 128);

            auto g_z_0_xxxyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 129);

            auto g_z_0_xxxyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 130);

            auto g_z_0_xxxyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 131);

            auto g_z_0_xxxyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 132);

            auto g_z_0_xxxyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 133);

            auto g_z_0_xxxyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 134);

            auto g_z_0_xxxyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 135);

            auto g_z_0_xxxyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 136);

            auto g_z_0_xxxyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 137);

            auto g_z_0_xxxyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 138);

            auto g_z_0_xxxyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 139);

            auto g_z_0_xxxyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 140);

            auto g_z_0_xxxyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 141);

            auto g_z_0_xxxyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 142);

            auto g_z_0_xxxyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 143);

            auto g_z_0_xxxyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 144);

            auto g_z_0_xxxyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 145);

            auto g_z_0_xxxyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyy_xxxxx, g_z_0_xxxyyy_xxxxy, g_z_0_xxxyyy_xxxxz, g_z_0_xxxyyy_xxxyy, g_z_0_xxxyyy_xxxyz, g_z_0_xxxyyy_xxxzz, g_z_0_xxxyyy_xxyyy, g_z_0_xxxyyy_xxyyz, g_z_0_xxxyyy_xxyzz, g_z_0_xxxyyy_xxzzz, g_z_0_xxxyyy_xyyyy, g_z_0_xxxyyy_xyyyz, g_z_0_xxxyyy_xyyzz, g_z_0_xxxyyy_xyzzz, g_z_0_xxxyyy_xzzzz, g_z_0_xxxyyy_yyyyy, g_z_0_xxxyyy_yyyyz, g_z_0_xxxyyy_yyyzz, g_z_0_xxxyyy_yyzzz, g_z_0_xxxyyy_yzzzz, g_z_0_xxxyyy_zzzzz, g_z_0_xxyyy_xxxxx, g_z_0_xxyyy_xxxxxx, g_z_0_xxyyy_xxxxxy, g_z_0_xxyyy_xxxxxz, g_z_0_xxyyy_xxxxy, g_z_0_xxyyy_xxxxyy, g_z_0_xxyyy_xxxxyz, g_z_0_xxyyy_xxxxz, g_z_0_xxyyy_xxxxzz, g_z_0_xxyyy_xxxyy, g_z_0_xxyyy_xxxyyy, g_z_0_xxyyy_xxxyyz, g_z_0_xxyyy_xxxyz, g_z_0_xxyyy_xxxyzz, g_z_0_xxyyy_xxxzz, g_z_0_xxyyy_xxxzzz, g_z_0_xxyyy_xxyyy, g_z_0_xxyyy_xxyyyy, g_z_0_xxyyy_xxyyyz, g_z_0_xxyyy_xxyyz, g_z_0_xxyyy_xxyyzz, g_z_0_xxyyy_xxyzz, g_z_0_xxyyy_xxyzzz, g_z_0_xxyyy_xxzzz, g_z_0_xxyyy_xxzzzz, g_z_0_xxyyy_xyyyy, g_z_0_xxyyy_xyyyyy, g_z_0_xxyyy_xyyyyz, g_z_0_xxyyy_xyyyz, g_z_0_xxyyy_xyyyzz, g_z_0_xxyyy_xyyzz, g_z_0_xxyyy_xyyzzz, g_z_0_xxyyy_xyzzz, g_z_0_xxyyy_xyzzzz, g_z_0_xxyyy_xzzzz, g_z_0_xxyyy_xzzzzz, g_z_0_xxyyy_yyyyy, g_z_0_xxyyy_yyyyz, g_z_0_xxyyy_yyyzz, g_z_0_xxyyy_yyzzz, g_z_0_xxyyy_yzzzz, g_z_0_xxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_xxxxx[k] = -g_z_0_xxyyy_xxxxx[k] * cd_x[k] + g_z_0_xxyyy_xxxxxx[k];

                g_z_0_xxxyyy_xxxxy[k] = -g_z_0_xxyyy_xxxxy[k] * cd_x[k] + g_z_0_xxyyy_xxxxxy[k];

                g_z_0_xxxyyy_xxxxz[k] = -g_z_0_xxyyy_xxxxz[k] * cd_x[k] + g_z_0_xxyyy_xxxxxz[k];

                g_z_0_xxxyyy_xxxyy[k] = -g_z_0_xxyyy_xxxyy[k] * cd_x[k] + g_z_0_xxyyy_xxxxyy[k];

                g_z_0_xxxyyy_xxxyz[k] = -g_z_0_xxyyy_xxxyz[k] * cd_x[k] + g_z_0_xxyyy_xxxxyz[k];

                g_z_0_xxxyyy_xxxzz[k] = -g_z_0_xxyyy_xxxzz[k] * cd_x[k] + g_z_0_xxyyy_xxxxzz[k];

                g_z_0_xxxyyy_xxyyy[k] = -g_z_0_xxyyy_xxyyy[k] * cd_x[k] + g_z_0_xxyyy_xxxyyy[k];

                g_z_0_xxxyyy_xxyyz[k] = -g_z_0_xxyyy_xxyyz[k] * cd_x[k] + g_z_0_xxyyy_xxxyyz[k];

                g_z_0_xxxyyy_xxyzz[k] = -g_z_0_xxyyy_xxyzz[k] * cd_x[k] + g_z_0_xxyyy_xxxyzz[k];

                g_z_0_xxxyyy_xxzzz[k] = -g_z_0_xxyyy_xxzzz[k] * cd_x[k] + g_z_0_xxyyy_xxxzzz[k];

                g_z_0_xxxyyy_xyyyy[k] = -g_z_0_xxyyy_xyyyy[k] * cd_x[k] + g_z_0_xxyyy_xxyyyy[k];

                g_z_0_xxxyyy_xyyyz[k] = -g_z_0_xxyyy_xyyyz[k] * cd_x[k] + g_z_0_xxyyy_xxyyyz[k];

                g_z_0_xxxyyy_xyyzz[k] = -g_z_0_xxyyy_xyyzz[k] * cd_x[k] + g_z_0_xxyyy_xxyyzz[k];

                g_z_0_xxxyyy_xyzzz[k] = -g_z_0_xxyyy_xyzzz[k] * cd_x[k] + g_z_0_xxyyy_xxyzzz[k];

                g_z_0_xxxyyy_xzzzz[k] = -g_z_0_xxyyy_xzzzz[k] * cd_x[k] + g_z_0_xxyyy_xxzzzz[k];

                g_z_0_xxxyyy_yyyyy[k] = -g_z_0_xxyyy_yyyyy[k] * cd_x[k] + g_z_0_xxyyy_xyyyyy[k];

                g_z_0_xxxyyy_yyyyz[k] = -g_z_0_xxyyy_yyyyz[k] * cd_x[k] + g_z_0_xxyyy_xyyyyz[k];

                g_z_0_xxxyyy_yyyzz[k] = -g_z_0_xxyyy_yyyzz[k] * cd_x[k] + g_z_0_xxyyy_xyyyzz[k];

                g_z_0_xxxyyy_yyzzz[k] = -g_z_0_xxyyy_yyzzz[k] * cd_x[k] + g_z_0_xxyyy_xyyzzz[k];

                g_z_0_xxxyyy_yzzzz[k] = -g_z_0_xxyyy_yzzzz[k] * cd_x[k] + g_z_0_xxyyy_xyzzzz[k];

                g_z_0_xxxyyy_zzzzz[k] = -g_z_0_xxyyy_zzzzz[k] * cd_x[k] + g_z_0_xxyyy_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 147);

            auto g_z_0_xxxyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 148);

            auto g_z_0_xxxyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 149);

            auto g_z_0_xxxyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 150);

            auto g_z_0_xxxyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 151);

            auto g_z_0_xxxyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 152);

            auto g_z_0_xxxyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 153);

            auto g_z_0_xxxyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 154);

            auto g_z_0_xxxyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 155);

            auto g_z_0_xxxyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 156);

            auto g_z_0_xxxyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 157);

            auto g_z_0_xxxyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 158);

            auto g_z_0_xxxyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 159);

            auto g_z_0_xxxyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 160);

            auto g_z_0_xxxyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 161);

            auto g_z_0_xxxyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 162);

            auto g_z_0_xxxyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 163);

            auto g_z_0_xxxyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 164);

            auto g_z_0_xxxyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 165);

            auto g_z_0_xxxyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 166);

            auto g_z_0_xxxyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyz_xxxxx, g_z_0_xxxyyz_xxxxy, g_z_0_xxxyyz_xxxxz, g_z_0_xxxyyz_xxxyy, g_z_0_xxxyyz_xxxyz, g_z_0_xxxyyz_xxxzz, g_z_0_xxxyyz_xxyyy, g_z_0_xxxyyz_xxyyz, g_z_0_xxxyyz_xxyzz, g_z_0_xxxyyz_xxzzz, g_z_0_xxxyyz_xyyyy, g_z_0_xxxyyz_xyyyz, g_z_0_xxxyyz_xyyzz, g_z_0_xxxyyz_xyzzz, g_z_0_xxxyyz_xzzzz, g_z_0_xxxyyz_yyyyy, g_z_0_xxxyyz_yyyyz, g_z_0_xxxyyz_yyyzz, g_z_0_xxxyyz_yyzzz, g_z_0_xxxyyz_yzzzz, g_z_0_xxxyyz_zzzzz, g_z_0_xxyyz_xxxxx, g_z_0_xxyyz_xxxxxx, g_z_0_xxyyz_xxxxxy, g_z_0_xxyyz_xxxxxz, g_z_0_xxyyz_xxxxy, g_z_0_xxyyz_xxxxyy, g_z_0_xxyyz_xxxxyz, g_z_0_xxyyz_xxxxz, g_z_0_xxyyz_xxxxzz, g_z_0_xxyyz_xxxyy, g_z_0_xxyyz_xxxyyy, g_z_0_xxyyz_xxxyyz, g_z_0_xxyyz_xxxyz, g_z_0_xxyyz_xxxyzz, g_z_0_xxyyz_xxxzz, g_z_0_xxyyz_xxxzzz, g_z_0_xxyyz_xxyyy, g_z_0_xxyyz_xxyyyy, g_z_0_xxyyz_xxyyyz, g_z_0_xxyyz_xxyyz, g_z_0_xxyyz_xxyyzz, g_z_0_xxyyz_xxyzz, g_z_0_xxyyz_xxyzzz, g_z_0_xxyyz_xxzzz, g_z_0_xxyyz_xxzzzz, g_z_0_xxyyz_xyyyy, g_z_0_xxyyz_xyyyyy, g_z_0_xxyyz_xyyyyz, g_z_0_xxyyz_xyyyz, g_z_0_xxyyz_xyyyzz, g_z_0_xxyyz_xyyzz, g_z_0_xxyyz_xyyzzz, g_z_0_xxyyz_xyzzz, g_z_0_xxyyz_xyzzzz, g_z_0_xxyyz_xzzzz, g_z_0_xxyyz_xzzzzz, g_z_0_xxyyz_yyyyy, g_z_0_xxyyz_yyyyz, g_z_0_xxyyz_yyyzz, g_z_0_xxyyz_yyzzz, g_z_0_xxyyz_yzzzz, g_z_0_xxyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_xxxxx[k] = -g_z_0_xxyyz_xxxxx[k] * cd_x[k] + g_z_0_xxyyz_xxxxxx[k];

                g_z_0_xxxyyz_xxxxy[k] = -g_z_0_xxyyz_xxxxy[k] * cd_x[k] + g_z_0_xxyyz_xxxxxy[k];

                g_z_0_xxxyyz_xxxxz[k] = -g_z_0_xxyyz_xxxxz[k] * cd_x[k] + g_z_0_xxyyz_xxxxxz[k];

                g_z_0_xxxyyz_xxxyy[k] = -g_z_0_xxyyz_xxxyy[k] * cd_x[k] + g_z_0_xxyyz_xxxxyy[k];

                g_z_0_xxxyyz_xxxyz[k] = -g_z_0_xxyyz_xxxyz[k] * cd_x[k] + g_z_0_xxyyz_xxxxyz[k];

                g_z_0_xxxyyz_xxxzz[k] = -g_z_0_xxyyz_xxxzz[k] * cd_x[k] + g_z_0_xxyyz_xxxxzz[k];

                g_z_0_xxxyyz_xxyyy[k] = -g_z_0_xxyyz_xxyyy[k] * cd_x[k] + g_z_0_xxyyz_xxxyyy[k];

                g_z_0_xxxyyz_xxyyz[k] = -g_z_0_xxyyz_xxyyz[k] * cd_x[k] + g_z_0_xxyyz_xxxyyz[k];

                g_z_0_xxxyyz_xxyzz[k] = -g_z_0_xxyyz_xxyzz[k] * cd_x[k] + g_z_0_xxyyz_xxxyzz[k];

                g_z_0_xxxyyz_xxzzz[k] = -g_z_0_xxyyz_xxzzz[k] * cd_x[k] + g_z_0_xxyyz_xxxzzz[k];

                g_z_0_xxxyyz_xyyyy[k] = -g_z_0_xxyyz_xyyyy[k] * cd_x[k] + g_z_0_xxyyz_xxyyyy[k];

                g_z_0_xxxyyz_xyyyz[k] = -g_z_0_xxyyz_xyyyz[k] * cd_x[k] + g_z_0_xxyyz_xxyyyz[k];

                g_z_0_xxxyyz_xyyzz[k] = -g_z_0_xxyyz_xyyzz[k] * cd_x[k] + g_z_0_xxyyz_xxyyzz[k];

                g_z_0_xxxyyz_xyzzz[k] = -g_z_0_xxyyz_xyzzz[k] * cd_x[k] + g_z_0_xxyyz_xxyzzz[k];

                g_z_0_xxxyyz_xzzzz[k] = -g_z_0_xxyyz_xzzzz[k] * cd_x[k] + g_z_0_xxyyz_xxzzzz[k];

                g_z_0_xxxyyz_yyyyy[k] = -g_z_0_xxyyz_yyyyy[k] * cd_x[k] + g_z_0_xxyyz_xyyyyy[k];

                g_z_0_xxxyyz_yyyyz[k] = -g_z_0_xxyyz_yyyyz[k] * cd_x[k] + g_z_0_xxyyz_xyyyyz[k];

                g_z_0_xxxyyz_yyyzz[k] = -g_z_0_xxyyz_yyyzz[k] * cd_x[k] + g_z_0_xxyyz_xyyyzz[k];

                g_z_0_xxxyyz_yyzzz[k] = -g_z_0_xxyyz_yyzzz[k] * cd_x[k] + g_z_0_xxyyz_xyyzzz[k];

                g_z_0_xxxyyz_yzzzz[k] = -g_z_0_xxyyz_yzzzz[k] * cd_x[k] + g_z_0_xxyyz_xyzzzz[k];

                g_z_0_xxxyyz_zzzzz[k] = -g_z_0_xxyyz_zzzzz[k] * cd_x[k] + g_z_0_xxyyz_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 168);

            auto g_z_0_xxxyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 169);

            auto g_z_0_xxxyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 170);

            auto g_z_0_xxxyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 171);

            auto g_z_0_xxxyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 172);

            auto g_z_0_xxxyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 173);

            auto g_z_0_xxxyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 174);

            auto g_z_0_xxxyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 175);

            auto g_z_0_xxxyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 176);

            auto g_z_0_xxxyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 177);

            auto g_z_0_xxxyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 178);

            auto g_z_0_xxxyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 179);

            auto g_z_0_xxxyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 180);

            auto g_z_0_xxxyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 181);

            auto g_z_0_xxxyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 182);

            auto g_z_0_xxxyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 183);

            auto g_z_0_xxxyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 184);

            auto g_z_0_xxxyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 185);

            auto g_z_0_xxxyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 186);

            auto g_z_0_xxxyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 187);

            auto g_z_0_xxxyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyzz_xxxxx, g_z_0_xxxyzz_xxxxy, g_z_0_xxxyzz_xxxxz, g_z_0_xxxyzz_xxxyy, g_z_0_xxxyzz_xxxyz, g_z_0_xxxyzz_xxxzz, g_z_0_xxxyzz_xxyyy, g_z_0_xxxyzz_xxyyz, g_z_0_xxxyzz_xxyzz, g_z_0_xxxyzz_xxzzz, g_z_0_xxxyzz_xyyyy, g_z_0_xxxyzz_xyyyz, g_z_0_xxxyzz_xyyzz, g_z_0_xxxyzz_xyzzz, g_z_0_xxxyzz_xzzzz, g_z_0_xxxyzz_yyyyy, g_z_0_xxxyzz_yyyyz, g_z_0_xxxyzz_yyyzz, g_z_0_xxxyzz_yyzzz, g_z_0_xxxyzz_yzzzz, g_z_0_xxxyzz_zzzzz, g_z_0_xxyzz_xxxxx, g_z_0_xxyzz_xxxxxx, g_z_0_xxyzz_xxxxxy, g_z_0_xxyzz_xxxxxz, g_z_0_xxyzz_xxxxy, g_z_0_xxyzz_xxxxyy, g_z_0_xxyzz_xxxxyz, g_z_0_xxyzz_xxxxz, g_z_0_xxyzz_xxxxzz, g_z_0_xxyzz_xxxyy, g_z_0_xxyzz_xxxyyy, g_z_0_xxyzz_xxxyyz, g_z_0_xxyzz_xxxyz, g_z_0_xxyzz_xxxyzz, g_z_0_xxyzz_xxxzz, g_z_0_xxyzz_xxxzzz, g_z_0_xxyzz_xxyyy, g_z_0_xxyzz_xxyyyy, g_z_0_xxyzz_xxyyyz, g_z_0_xxyzz_xxyyz, g_z_0_xxyzz_xxyyzz, g_z_0_xxyzz_xxyzz, g_z_0_xxyzz_xxyzzz, g_z_0_xxyzz_xxzzz, g_z_0_xxyzz_xxzzzz, g_z_0_xxyzz_xyyyy, g_z_0_xxyzz_xyyyyy, g_z_0_xxyzz_xyyyyz, g_z_0_xxyzz_xyyyz, g_z_0_xxyzz_xyyyzz, g_z_0_xxyzz_xyyzz, g_z_0_xxyzz_xyyzzz, g_z_0_xxyzz_xyzzz, g_z_0_xxyzz_xyzzzz, g_z_0_xxyzz_xzzzz, g_z_0_xxyzz_xzzzzz, g_z_0_xxyzz_yyyyy, g_z_0_xxyzz_yyyyz, g_z_0_xxyzz_yyyzz, g_z_0_xxyzz_yyzzz, g_z_0_xxyzz_yzzzz, g_z_0_xxyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_xxxxx[k] = -g_z_0_xxyzz_xxxxx[k] * cd_x[k] + g_z_0_xxyzz_xxxxxx[k];

                g_z_0_xxxyzz_xxxxy[k] = -g_z_0_xxyzz_xxxxy[k] * cd_x[k] + g_z_0_xxyzz_xxxxxy[k];

                g_z_0_xxxyzz_xxxxz[k] = -g_z_0_xxyzz_xxxxz[k] * cd_x[k] + g_z_0_xxyzz_xxxxxz[k];

                g_z_0_xxxyzz_xxxyy[k] = -g_z_0_xxyzz_xxxyy[k] * cd_x[k] + g_z_0_xxyzz_xxxxyy[k];

                g_z_0_xxxyzz_xxxyz[k] = -g_z_0_xxyzz_xxxyz[k] * cd_x[k] + g_z_0_xxyzz_xxxxyz[k];

                g_z_0_xxxyzz_xxxzz[k] = -g_z_0_xxyzz_xxxzz[k] * cd_x[k] + g_z_0_xxyzz_xxxxzz[k];

                g_z_0_xxxyzz_xxyyy[k] = -g_z_0_xxyzz_xxyyy[k] * cd_x[k] + g_z_0_xxyzz_xxxyyy[k];

                g_z_0_xxxyzz_xxyyz[k] = -g_z_0_xxyzz_xxyyz[k] * cd_x[k] + g_z_0_xxyzz_xxxyyz[k];

                g_z_0_xxxyzz_xxyzz[k] = -g_z_0_xxyzz_xxyzz[k] * cd_x[k] + g_z_0_xxyzz_xxxyzz[k];

                g_z_0_xxxyzz_xxzzz[k] = -g_z_0_xxyzz_xxzzz[k] * cd_x[k] + g_z_0_xxyzz_xxxzzz[k];

                g_z_0_xxxyzz_xyyyy[k] = -g_z_0_xxyzz_xyyyy[k] * cd_x[k] + g_z_0_xxyzz_xxyyyy[k];

                g_z_0_xxxyzz_xyyyz[k] = -g_z_0_xxyzz_xyyyz[k] * cd_x[k] + g_z_0_xxyzz_xxyyyz[k];

                g_z_0_xxxyzz_xyyzz[k] = -g_z_0_xxyzz_xyyzz[k] * cd_x[k] + g_z_0_xxyzz_xxyyzz[k];

                g_z_0_xxxyzz_xyzzz[k] = -g_z_0_xxyzz_xyzzz[k] * cd_x[k] + g_z_0_xxyzz_xxyzzz[k];

                g_z_0_xxxyzz_xzzzz[k] = -g_z_0_xxyzz_xzzzz[k] * cd_x[k] + g_z_0_xxyzz_xxzzzz[k];

                g_z_0_xxxyzz_yyyyy[k] = -g_z_0_xxyzz_yyyyy[k] * cd_x[k] + g_z_0_xxyzz_xyyyyy[k];

                g_z_0_xxxyzz_yyyyz[k] = -g_z_0_xxyzz_yyyyz[k] * cd_x[k] + g_z_0_xxyzz_xyyyyz[k];

                g_z_0_xxxyzz_yyyzz[k] = -g_z_0_xxyzz_yyyzz[k] * cd_x[k] + g_z_0_xxyzz_xyyyzz[k];

                g_z_0_xxxyzz_yyzzz[k] = -g_z_0_xxyzz_yyzzz[k] * cd_x[k] + g_z_0_xxyzz_xyyzzz[k];

                g_z_0_xxxyzz_yzzzz[k] = -g_z_0_xxyzz_yzzzz[k] * cd_x[k] + g_z_0_xxyzz_xyzzzz[k];

                g_z_0_xxxyzz_zzzzz[k] = -g_z_0_xxyzz_zzzzz[k] * cd_x[k] + g_z_0_xxyzz_xzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 189);

            auto g_z_0_xxxzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 190);

            auto g_z_0_xxxzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 191);

            auto g_z_0_xxxzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 192);

            auto g_z_0_xxxzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 193);

            auto g_z_0_xxxzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 194);

            auto g_z_0_xxxzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 195);

            auto g_z_0_xxxzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 196);

            auto g_z_0_xxxzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 197);

            auto g_z_0_xxxzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 198);

            auto g_z_0_xxxzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 199);

            auto g_z_0_xxxzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 200);

            auto g_z_0_xxxzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 201);

            auto g_z_0_xxxzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 202);

            auto g_z_0_xxxzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 203);

            auto g_z_0_xxxzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 204);

            auto g_z_0_xxxzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 205);

            auto g_z_0_xxxzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 206);

            auto g_z_0_xxxzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 207);

            auto g_z_0_xxxzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 208);

            auto g_z_0_xxxzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzzz_xxxxx, g_z_0_xxxzzz_xxxxy, g_z_0_xxxzzz_xxxxz, g_z_0_xxxzzz_xxxyy, g_z_0_xxxzzz_xxxyz, g_z_0_xxxzzz_xxxzz, g_z_0_xxxzzz_xxyyy, g_z_0_xxxzzz_xxyyz, g_z_0_xxxzzz_xxyzz, g_z_0_xxxzzz_xxzzz, g_z_0_xxxzzz_xyyyy, g_z_0_xxxzzz_xyyyz, g_z_0_xxxzzz_xyyzz, g_z_0_xxxzzz_xyzzz, g_z_0_xxxzzz_xzzzz, g_z_0_xxxzzz_yyyyy, g_z_0_xxxzzz_yyyyz, g_z_0_xxxzzz_yyyzz, g_z_0_xxxzzz_yyzzz, g_z_0_xxxzzz_yzzzz, g_z_0_xxxzzz_zzzzz, g_z_0_xxzzz_xxxxx, g_z_0_xxzzz_xxxxxx, g_z_0_xxzzz_xxxxxy, g_z_0_xxzzz_xxxxxz, g_z_0_xxzzz_xxxxy, g_z_0_xxzzz_xxxxyy, g_z_0_xxzzz_xxxxyz, g_z_0_xxzzz_xxxxz, g_z_0_xxzzz_xxxxzz, g_z_0_xxzzz_xxxyy, g_z_0_xxzzz_xxxyyy, g_z_0_xxzzz_xxxyyz, g_z_0_xxzzz_xxxyz, g_z_0_xxzzz_xxxyzz, g_z_0_xxzzz_xxxzz, g_z_0_xxzzz_xxxzzz, g_z_0_xxzzz_xxyyy, g_z_0_xxzzz_xxyyyy, g_z_0_xxzzz_xxyyyz, g_z_0_xxzzz_xxyyz, g_z_0_xxzzz_xxyyzz, g_z_0_xxzzz_xxyzz, g_z_0_xxzzz_xxyzzz, g_z_0_xxzzz_xxzzz, g_z_0_xxzzz_xxzzzz, g_z_0_xxzzz_xyyyy, g_z_0_xxzzz_xyyyyy, g_z_0_xxzzz_xyyyyz, g_z_0_xxzzz_xyyyz, g_z_0_xxzzz_xyyyzz, g_z_0_xxzzz_xyyzz, g_z_0_xxzzz_xyyzzz, g_z_0_xxzzz_xyzzz, g_z_0_xxzzz_xyzzzz, g_z_0_xxzzz_xzzzz, g_z_0_xxzzz_xzzzzz, g_z_0_xxzzz_yyyyy, g_z_0_xxzzz_yyyyz, g_z_0_xxzzz_yyyzz, g_z_0_xxzzz_yyzzz, g_z_0_xxzzz_yzzzz, g_z_0_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_xxxxx[k] = -g_z_0_xxzzz_xxxxx[k] * cd_x[k] + g_z_0_xxzzz_xxxxxx[k];

                g_z_0_xxxzzz_xxxxy[k] = -g_z_0_xxzzz_xxxxy[k] * cd_x[k] + g_z_0_xxzzz_xxxxxy[k];

                g_z_0_xxxzzz_xxxxz[k] = -g_z_0_xxzzz_xxxxz[k] * cd_x[k] + g_z_0_xxzzz_xxxxxz[k];

                g_z_0_xxxzzz_xxxyy[k] = -g_z_0_xxzzz_xxxyy[k] * cd_x[k] + g_z_0_xxzzz_xxxxyy[k];

                g_z_0_xxxzzz_xxxyz[k] = -g_z_0_xxzzz_xxxyz[k] * cd_x[k] + g_z_0_xxzzz_xxxxyz[k];

                g_z_0_xxxzzz_xxxzz[k] = -g_z_0_xxzzz_xxxzz[k] * cd_x[k] + g_z_0_xxzzz_xxxxzz[k];

                g_z_0_xxxzzz_xxyyy[k] = -g_z_0_xxzzz_xxyyy[k] * cd_x[k] + g_z_0_xxzzz_xxxyyy[k];

                g_z_0_xxxzzz_xxyyz[k] = -g_z_0_xxzzz_xxyyz[k] * cd_x[k] + g_z_0_xxzzz_xxxyyz[k];

                g_z_0_xxxzzz_xxyzz[k] = -g_z_0_xxzzz_xxyzz[k] * cd_x[k] + g_z_0_xxzzz_xxxyzz[k];

                g_z_0_xxxzzz_xxzzz[k] = -g_z_0_xxzzz_xxzzz[k] * cd_x[k] + g_z_0_xxzzz_xxxzzz[k];

                g_z_0_xxxzzz_xyyyy[k] = -g_z_0_xxzzz_xyyyy[k] * cd_x[k] + g_z_0_xxzzz_xxyyyy[k];

                g_z_0_xxxzzz_xyyyz[k] = -g_z_0_xxzzz_xyyyz[k] * cd_x[k] + g_z_0_xxzzz_xxyyyz[k];

                g_z_0_xxxzzz_xyyzz[k] = -g_z_0_xxzzz_xyyzz[k] * cd_x[k] + g_z_0_xxzzz_xxyyzz[k];

                g_z_0_xxxzzz_xyzzz[k] = -g_z_0_xxzzz_xyzzz[k] * cd_x[k] + g_z_0_xxzzz_xxyzzz[k];

                g_z_0_xxxzzz_xzzzz[k] = -g_z_0_xxzzz_xzzzz[k] * cd_x[k] + g_z_0_xxzzz_xxzzzz[k];

                g_z_0_xxxzzz_yyyyy[k] = -g_z_0_xxzzz_yyyyy[k] * cd_x[k] + g_z_0_xxzzz_xyyyyy[k];

                g_z_0_xxxzzz_yyyyz[k] = -g_z_0_xxzzz_yyyyz[k] * cd_x[k] + g_z_0_xxzzz_xyyyyz[k];

                g_z_0_xxxzzz_yyyzz[k] = -g_z_0_xxzzz_yyyzz[k] * cd_x[k] + g_z_0_xxzzz_xyyyzz[k];

                g_z_0_xxxzzz_yyzzz[k] = -g_z_0_xxzzz_yyzzz[k] * cd_x[k] + g_z_0_xxzzz_xyyzzz[k];

                g_z_0_xxxzzz_yzzzz[k] = -g_z_0_xxzzz_yzzzz[k] * cd_x[k] + g_z_0_xxzzz_xyzzzz[k];

                g_z_0_xxxzzz_zzzzz[k] = -g_z_0_xxzzz_zzzzz[k] * cd_x[k] + g_z_0_xxzzz_xzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 210);

            auto g_z_0_xxyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 211);

            auto g_z_0_xxyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 212);

            auto g_z_0_xxyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 213);

            auto g_z_0_xxyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 214);

            auto g_z_0_xxyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 215);

            auto g_z_0_xxyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 216);

            auto g_z_0_xxyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 217);

            auto g_z_0_xxyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 218);

            auto g_z_0_xxyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 219);

            auto g_z_0_xxyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 220);

            auto g_z_0_xxyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 221);

            auto g_z_0_xxyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 222);

            auto g_z_0_xxyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 223);

            auto g_z_0_xxyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 224);

            auto g_z_0_xxyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 225);

            auto g_z_0_xxyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 226);

            auto g_z_0_xxyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 227);

            auto g_z_0_xxyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 228);

            auto g_z_0_xxyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 229);

            auto g_z_0_xxyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyy_xxxxx, g_z_0_xxyyyy_xxxxy, g_z_0_xxyyyy_xxxxz, g_z_0_xxyyyy_xxxyy, g_z_0_xxyyyy_xxxyz, g_z_0_xxyyyy_xxxzz, g_z_0_xxyyyy_xxyyy, g_z_0_xxyyyy_xxyyz, g_z_0_xxyyyy_xxyzz, g_z_0_xxyyyy_xxzzz, g_z_0_xxyyyy_xyyyy, g_z_0_xxyyyy_xyyyz, g_z_0_xxyyyy_xyyzz, g_z_0_xxyyyy_xyzzz, g_z_0_xxyyyy_xzzzz, g_z_0_xxyyyy_yyyyy, g_z_0_xxyyyy_yyyyz, g_z_0_xxyyyy_yyyzz, g_z_0_xxyyyy_yyzzz, g_z_0_xxyyyy_yzzzz, g_z_0_xxyyyy_zzzzz, g_z_0_xyyyy_xxxxx, g_z_0_xyyyy_xxxxxx, g_z_0_xyyyy_xxxxxy, g_z_0_xyyyy_xxxxxz, g_z_0_xyyyy_xxxxy, g_z_0_xyyyy_xxxxyy, g_z_0_xyyyy_xxxxyz, g_z_0_xyyyy_xxxxz, g_z_0_xyyyy_xxxxzz, g_z_0_xyyyy_xxxyy, g_z_0_xyyyy_xxxyyy, g_z_0_xyyyy_xxxyyz, g_z_0_xyyyy_xxxyz, g_z_0_xyyyy_xxxyzz, g_z_0_xyyyy_xxxzz, g_z_0_xyyyy_xxxzzz, g_z_0_xyyyy_xxyyy, g_z_0_xyyyy_xxyyyy, g_z_0_xyyyy_xxyyyz, g_z_0_xyyyy_xxyyz, g_z_0_xyyyy_xxyyzz, g_z_0_xyyyy_xxyzz, g_z_0_xyyyy_xxyzzz, g_z_0_xyyyy_xxzzz, g_z_0_xyyyy_xxzzzz, g_z_0_xyyyy_xyyyy, g_z_0_xyyyy_xyyyyy, g_z_0_xyyyy_xyyyyz, g_z_0_xyyyy_xyyyz, g_z_0_xyyyy_xyyyzz, g_z_0_xyyyy_xyyzz, g_z_0_xyyyy_xyyzzz, g_z_0_xyyyy_xyzzz, g_z_0_xyyyy_xyzzzz, g_z_0_xyyyy_xzzzz, g_z_0_xyyyy_xzzzzz, g_z_0_xyyyy_yyyyy, g_z_0_xyyyy_yyyyz, g_z_0_xyyyy_yyyzz, g_z_0_xyyyy_yyzzz, g_z_0_xyyyy_yzzzz, g_z_0_xyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_xxxxx[k] = -g_z_0_xyyyy_xxxxx[k] * cd_x[k] + g_z_0_xyyyy_xxxxxx[k];

                g_z_0_xxyyyy_xxxxy[k] = -g_z_0_xyyyy_xxxxy[k] * cd_x[k] + g_z_0_xyyyy_xxxxxy[k];

                g_z_0_xxyyyy_xxxxz[k] = -g_z_0_xyyyy_xxxxz[k] * cd_x[k] + g_z_0_xyyyy_xxxxxz[k];

                g_z_0_xxyyyy_xxxyy[k] = -g_z_0_xyyyy_xxxyy[k] * cd_x[k] + g_z_0_xyyyy_xxxxyy[k];

                g_z_0_xxyyyy_xxxyz[k] = -g_z_0_xyyyy_xxxyz[k] * cd_x[k] + g_z_0_xyyyy_xxxxyz[k];

                g_z_0_xxyyyy_xxxzz[k] = -g_z_0_xyyyy_xxxzz[k] * cd_x[k] + g_z_0_xyyyy_xxxxzz[k];

                g_z_0_xxyyyy_xxyyy[k] = -g_z_0_xyyyy_xxyyy[k] * cd_x[k] + g_z_0_xyyyy_xxxyyy[k];

                g_z_0_xxyyyy_xxyyz[k] = -g_z_0_xyyyy_xxyyz[k] * cd_x[k] + g_z_0_xyyyy_xxxyyz[k];

                g_z_0_xxyyyy_xxyzz[k] = -g_z_0_xyyyy_xxyzz[k] * cd_x[k] + g_z_0_xyyyy_xxxyzz[k];

                g_z_0_xxyyyy_xxzzz[k] = -g_z_0_xyyyy_xxzzz[k] * cd_x[k] + g_z_0_xyyyy_xxxzzz[k];

                g_z_0_xxyyyy_xyyyy[k] = -g_z_0_xyyyy_xyyyy[k] * cd_x[k] + g_z_0_xyyyy_xxyyyy[k];

                g_z_0_xxyyyy_xyyyz[k] = -g_z_0_xyyyy_xyyyz[k] * cd_x[k] + g_z_0_xyyyy_xxyyyz[k];

                g_z_0_xxyyyy_xyyzz[k] = -g_z_0_xyyyy_xyyzz[k] * cd_x[k] + g_z_0_xyyyy_xxyyzz[k];

                g_z_0_xxyyyy_xyzzz[k] = -g_z_0_xyyyy_xyzzz[k] * cd_x[k] + g_z_0_xyyyy_xxyzzz[k];

                g_z_0_xxyyyy_xzzzz[k] = -g_z_0_xyyyy_xzzzz[k] * cd_x[k] + g_z_0_xyyyy_xxzzzz[k];

                g_z_0_xxyyyy_yyyyy[k] = -g_z_0_xyyyy_yyyyy[k] * cd_x[k] + g_z_0_xyyyy_xyyyyy[k];

                g_z_0_xxyyyy_yyyyz[k] = -g_z_0_xyyyy_yyyyz[k] * cd_x[k] + g_z_0_xyyyy_xyyyyz[k];

                g_z_0_xxyyyy_yyyzz[k] = -g_z_0_xyyyy_yyyzz[k] * cd_x[k] + g_z_0_xyyyy_xyyyzz[k];

                g_z_0_xxyyyy_yyzzz[k] = -g_z_0_xyyyy_yyzzz[k] * cd_x[k] + g_z_0_xyyyy_xyyzzz[k];

                g_z_0_xxyyyy_yzzzz[k] = -g_z_0_xyyyy_yzzzz[k] * cd_x[k] + g_z_0_xyyyy_xyzzzz[k];

                g_z_0_xxyyyy_zzzzz[k] = -g_z_0_xyyyy_zzzzz[k] * cd_x[k] + g_z_0_xyyyy_xzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 231);

            auto g_z_0_xxyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 232);

            auto g_z_0_xxyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 233);

            auto g_z_0_xxyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 234);

            auto g_z_0_xxyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 235);

            auto g_z_0_xxyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 236);

            auto g_z_0_xxyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 237);

            auto g_z_0_xxyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 238);

            auto g_z_0_xxyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 239);

            auto g_z_0_xxyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 240);

            auto g_z_0_xxyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 241);

            auto g_z_0_xxyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 242);

            auto g_z_0_xxyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 243);

            auto g_z_0_xxyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 244);

            auto g_z_0_xxyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 245);

            auto g_z_0_xxyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 246);

            auto g_z_0_xxyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 247);

            auto g_z_0_xxyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 248);

            auto g_z_0_xxyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 249);

            auto g_z_0_xxyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 250);

            auto g_z_0_xxyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyz_xxxxx, g_z_0_xxyyyz_xxxxy, g_z_0_xxyyyz_xxxxz, g_z_0_xxyyyz_xxxyy, g_z_0_xxyyyz_xxxyz, g_z_0_xxyyyz_xxxzz, g_z_0_xxyyyz_xxyyy, g_z_0_xxyyyz_xxyyz, g_z_0_xxyyyz_xxyzz, g_z_0_xxyyyz_xxzzz, g_z_0_xxyyyz_xyyyy, g_z_0_xxyyyz_xyyyz, g_z_0_xxyyyz_xyyzz, g_z_0_xxyyyz_xyzzz, g_z_0_xxyyyz_xzzzz, g_z_0_xxyyyz_yyyyy, g_z_0_xxyyyz_yyyyz, g_z_0_xxyyyz_yyyzz, g_z_0_xxyyyz_yyzzz, g_z_0_xxyyyz_yzzzz, g_z_0_xxyyyz_zzzzz, g_z_0_xyyyz_xxxxx, g_z_0_xyyyz_xxxxxx, g_z_0_xyyyz_xxxxxy, g_z_0_xyyyz_xxxxxz, g_z_0_xyyyz_xxxxy, g_z_0_xyyyz_xxxxyy, g_z_0_xyyyz_xxxxyz, g_z_0_xyyyz_xxxxz, g_z_0_xyyyz_xxxxzz, g_z_0_xyyyz_xxxyy, g_z_0_xyyyz_xxxyyy, g_z_0_xyyyz_xxxyyz, g_z_0_xyyyz_xxxyz, g_z_0_xyyyz_xxxyzz, g_z_0_xyyyz_xxxzz, g_z_0_xyyyz_xxxzzz, g_z_0_xyyyz_xxyyy, g_z_0_xyyyz_xxyyyy, g_z_0_xyyyz_xxyyyz, g_z_0_xyyyz_xxyyz, g_z_0_xyyyz_xxyyzz, g_z_0_xyyyz_xxyzz, g_z_0_xyyyz_xxyzzz, g_z_0_xyyyz_xxzzz, g_z_0_xyyyz_xxzzzz, g_z_0_xyyyz_xyyyy, g_z_0_xyyyz_xyyyyy, g_z_0_xyyyz_xyyyyz, g_z_0_xyyyz_xyyyz, g_z_0_xyyyz_xyyyzz, g_z_0_xyyyz_xyyzz, g_z_0_xyyyz_xyyzzz, g_z_0_xyyyz_xyzzz, g_z_0_xyyyz_xyzzzz, g_z_0_xyyyz_xzzzz, g_z_0_xyyyz_xzzzzz, g_z_0_xyyyz_yyyyy, g_z_0_xyyyz_yyyyz, g_z_0_xyyyz_yyyzz, g_z_0_xyyyz_yyzzz, g_z_0_xyyyz_yzzzz, g_z_0_xyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_xxxxx[k] = -g_z_0_xyyyz_xxxxx[k] * cd_x[k] + g_z_0_xyyyz_xxxxxx[k];

                g_z_0_xxyyyz_xxxxy[k] = -g_z_0_xyyyz_xxxxy[k] * cd_x[k] + g_z_0_xyyyz_xxxxxy[k];

                g_z_0_xxyyyz_xxxxz[k] = -g_z_0_xyyyz_xxxxz[k] * cd_x[k] + g_z_0_xyyyz_xxxxxz[k];

                g_z_0_xxyyyz_xxxyy[k] = -g_z_0_xyyyz_xxxyy[k] * cd_x[k] + g_z_0_xyyyz_xxxxyy[k];

                g_z_0_xxyyyz_xxxyz[k] = -g_z_0_xyyyz_xxxyz[k] * cd_x[k] + g_z_0_xyyyz_xxxxyz[k];

                g_z_0_xxyyyz_xxxzz[k] = -g_z_0_xyyyz_xxxzz[k] * cd_x[k] + g_z_0_xyyyz_xxxxzz[k];

                g_z_0_xxyyyz_xxyyy[k] = -g_z_0_xyyyz_xxyyy[k] * cd_x[k] + g_z_0_xyyyz_xxxyyy[k];

                g_z_0_xxyyyz_xxyyz[k] = -g_z_0_xyyyz_xxyyz[k] * cd_x[k] + g_z_0_xyyyz_xxxyyz[k];

                g_z_0_xxyyyz_xxyzz[k] = -g_z_0_xyyyz_xxyzz[k] * cd_x[k] + g_z_0_xyyyz_xxxyzz[k];

                g_z_0_xxyyyz_xxzzz[k] = -g_z_0_xyyyz_xxzzz[k] * cd_x[k] + g_z_0_xyyyz_xxxzzz[k];

                g_z_0_xxyyyz_xyyyy[k] = -g_z_0_xyyyz_xyyyy[k] * cd_x[k] + g_z_0_xyyyz_xxyyyy[k];

                g_z_0_xxyyyz_xyyyz[k] = -g_z_0_xyyyz_xyyyz[k] * cd_x[k] + g_z_0_xyyyz_xxyyyz[k];

                g_z_0_xxyyyz_xyyzz[k] = -g_z_0_xyyyz_xyyzz[k] * cd_x[k] + g_z_0_xyyyz_xxyyzz[k];

                g_z_0_xxyyyz_xyzzz[k] = -g_z_0_xyyyz_xyzzz[k] * cd_x[k] + g_z_0_xyyyz_xxyzzz[k];

                g_z_0_xxyyyz_xzzzz[k] = -g_z_0_xyyyz_xzzzz[k] * cd_x[k] + g_z_0_xyyyz_xxzzzz[k];

                g_z_0_xxyyyz_yyyyy[k] = -g_z_0_xyyyz_yyyyy[k] * cd_x[k] + g_z_0_xyyyz_xyyyyy[k];

                g_z_0_xxyyyz_yyyyz[k] = -g_z_0_xyyyz_yyyyz[k] * cd_x[k] + g_z_0_xyyyz_xyyyyz[k];

                g_z_0_xxyyyz_yyyzz[k] = -g_z_0_xyyyz_yyyzz[k] * cd_x[k] + g_z_0_xyyyz_xyyyzz[k];

                g_z_0_xxyyyz_yyzzz[k] = -g_z_0_xyyyz_yyzzz[k] * cd_x[k] + g_z_0_xyyyz_xyyzzz[k];

                g_z_0_xxyyyz_yzzzz[k] = -g_z_0_xyyyz_yzzzz[k] * cd_x[k] + g_z_0_xyyyz_xyzzzz[k];

                g_z_0_xxyyyz_zzzzz[k] = -g_z_0_xyyyz_zzzzz[k] * cd_x[k] + g_z_0_xyyyz_xzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 252);

            auto g_z_0_xxyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 253);

            auto g_z_0_xxyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 254);

            auto g_z_0_xxyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 255);

            auto g_z_0_xxyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 256);

            auto g_z_0_xxyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 257);

            auto g_z_0_xxyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 258);

            auto g_z_0_xxyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 259);

            auto g_z_0_xxyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 260);

            auto g_z_0_xxyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 261);

            auto g_z_0_xxyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 262);

            auto g_z_0_xxyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 263);

            auto g_z_0_xxyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 264);

            auto g_z_0_xxyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 265);

            auto g_z_0_xxyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 266);

            auto g_z_0_xxyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 267);

            auto g_z_0_xxyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 268);

            auto g_z_0_xxyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 269);

            auto g_z_0_xxyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 270);

            auto g_z_0_xxyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 271);

            auto g_z_0_xxyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyzz_xxxxx, g_z_0_xxyyzz_xxxxy, g_z_0_xxyyzz_xxxxz, g_z_0_xxyyzz_xxxyy, g_z_0_xxyyzz_xxxyz, g_z_0_xxyyzz_xxxzz, g_z_0_xxyyzz_xxyyy, g_z_0_xxyyzz_xxyyz, g_z_0_xxyyzz_xxyzz, g_z_0_xxyyzz_xxzzz, g_z_0_xxyyzz_xyyyy, g_z_0_xxyyzz_xyyyz, g_z_0_xxyyzz_xyyzz, g_z_0_xxyyzz_xyzzz, g_z_0_xxyyzz_xzzzz, g_z_0_xxyyzz_yyyyy, g_z_0_xxyyzz_yyyyz, g_z_0_xxyyzz_yyyzz, g_z_0_xxyyzz_yyzzz, g_z_0_xxyyzz_yzzzz, g_z_0_xxyyzz_zzzzz, g_z_0_xyyzz_xxxxx, g_z_0_xyyzz_xxxxxx, g_z_0_xyyzz_xxxxxy, g_z_0_xyyzz_xxxxxz, g_z_0_xyyzz_xxxxy, g_z_0_xyyzz_xxxxyy, g_z_0_xyyzz_xxxxyz, g_z_0_xyyzz_xxxxz, g_z_0_xyyzz_xxxxzz, g_z_0_xyyzz_xxxyy, g_z_0_xyyzz_xxxyyy, g_z_0_xyyzz_xxxyyz, g_z_0_xyyzz_xxxyz, g_z_0_xyyzz_xxxyzz, g_z_0_xyyzz_xxxzz, g_z_0_xyyzz_xxxzzz, g_z_0_xyyzz_xxyyy, g_z_0_xyyzz_xxyyyy, g_z_0_xyyzz_xxyyyz, g_z_0_xyyzz_xxyyz, g_z_0_xyyzz_xxyyzz, g_z_0_xyyzz_xxyzz, g_z_0_xyyzz_xxyzzz, g_z_0_xyyzz_xxzzz, g_z_0_xyyzz_xxzzzz, g_z_0_xyyzz_xyyyy, g_z_0_xyyzz_xyyyyy, g_z_0_xyyzz_xyyyyz, g_z_0_xyyzz_xyyyz, g_z_0_xyyzz_xyyyzz, g_z_0_xyyzz_xyyzz, g_z_0_xyyzz_xyyzzz, g_z_0_xyyzz_xyzzz, g_z_0_xyyzz_xyzzzz, g_z_0_xyyzz_xzzzz, g_z_0_xyyzz_xzzzzz, g_z_0_xyyzz_yyyyy, g_z_0_xyyzz_yyyyz, g_z_0_xyyzz_yyyzz, g_z_0_xyyzz_yyzzz, g_z_0_xyyzz_yzzzz, g_z_0_xyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_xxxxx[k] = -g_z_0_xyyzz_xxxxx[k] * cd_x[k] + g_z_0_xyyzz_xxxxxx[k];

                g_z_0_xxyyzz_xxxxy[k] = -g_z_0_xyyzz_xxxxy[k] * cd_x[k] + g_z_0_xyyzz_xxxxxy[k];

                g_z_0_xxyyzz_xxxxz[k] = -g_z_0_xyyzz_xxxxz[k] * cd_x[k] + g_z_0_xyyzz_xxxxxz[k];

                g_z_0_xxyyzz_xxxyy[k] = -g_z_0_xyyzz_xxxyy[k] * cd_x[k] + g_z_0_xyyzz_xxxxyy[k];

                g_z_0_xxyyzz_xxxyz[k] = -g_z_0_xyyzz_xxxyz[k] * cd_x[k] + g_z_0_xyyzz_xxxxyz[k];

                g_z_0_xxyyzz_xxxzz[k] = -g_z_0_xyyzz_xxxzz[k] * cd_x[k] + g_z_0_xyyzz_xxxxzz[k];

                g_z_0_xxyyzz_xxyyy[k] = -g_z_0_xyyzz_xxyyy[k] * cd_x[k] + g_z_0_xyyzz_xxxyyy[k];

                g_z_0_xxyyzz_xxyyz[k] = -g_z_0_xyyzz_xxyyz[k] * cd_x[k] + g_z_0_xyyzz_xxxyyz[k];

                g_z_0_xxyyzz_xxyzz[k] = -g_z_0_xyyzz_xxyzz[k] * cd_x[k] + g_z_0_xyyzz_xxxyzz[k];

                g_z_0_xxyyzz_xxzzz[k] = -g_z_0_xyyzz_xxzzz[k] * cd_x[k] + g_z_0_xyyzz_xxxzzz[k];

                g_z_0_xxyyzz_xyyyy[k] = -g_z_0_xyyzz_xyyyy[k] * cd_x[k] + g_z_0_xyyzz_xxyyyy[k];

                g_z_0_xxyyzz_xyyyz[k] = -g_z_0_xyyzz_xyyyz[k] * cd_x[k] + g_z_0_xyyzz_xxyyyz[k];

                g_z_0_xxyyzz_xyyzz[k] = -g_z_0_xyyzz_xyyzz[k] * cd_x[k] + g_z_0_xyyzz_xxyyzz[k];

                g_z_0_xxyyzz_xyzzz[k] = -g_z_0_xyyzz_xyzzz[k] * cd_x[k] + g_z_0_xyyzz_xxyzzz[k];

                g_z_0_xxyyzz_xzzzz[k] = -g_z_0_xyyzz_xzzzz[k] * cd_x[k] + g_z_0_xyyzz_xxzzzz[k];

                g_z_0_xxyyzz_yyyyy[k] = -g_z_0_xyyzz_yyyyy[k] * cd_x[k] + g_z_0_xyyzz_xyyyyy[k];

                g_z_0_xxyyzz_yyyyz[k] = -g_z_0_xyyzz_yyyyz[k] * cd_x[k] + g_z_0_xyyzz_xyyyyz[k];

                g_z_0_xxyyzz_yyyzz[k] = -g_z_0_xyyzz_yyyzz[k] * cd_x[k] + g_z_0_xyyzz_xyyyzz[k];

                g_z_0_xxyyzz_yyzzz[k] = -g_z_0_xyyzz_yyzzz[k] * cd_x[k] + g_z_0_xyyzz_xyyzzz[k];

                g_z_0_xxyyzz_yzzzz[k] = -g_z_0_xyyzz_yzzzz[k] * cd_x[k] + g_z_0_xyyzz_xyzzzz[k];

                g_z_0_xxyyzz_zzzzz[k] = -g_z_0_xyyzz_zzzzz[k] * cd_x[k] + g_z_0_xyyzz_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 273);

            auto g_z_0_xxyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 274);

            auto g_z_0_xxyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 275);

            auto g_z_0_xxyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 276);

            auto g_z_0_xxyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 277);

            auto g_z_0_xxyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 278);

            auto g_z_0_xxyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 279);

            auto g_z_0_xxyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 280);

            auto g_z_0_xxyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 281);

            auto g_z_0_xxyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 282);

            auto g_z_0_xxyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 283);

            auto g_z_0_xxyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 284);

            auto g_z_0_xxyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 285);

            auto g_z_0_xxyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 286);

            auto g_z_0_xxyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 287);

            auto g_z_0_xxyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 288);

            auto g_z_0_xxyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 289);

            auto g_z_0_xxyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 290);

            auto g_z_0_xxyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 291);

            auto g_z_0_xxyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 292);

            auto g_z_0_xxyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzzz_xxxxx, g_z_0_xxyzzz_xxxxy, g_z_0_xxyzzz_xxxxz, g_z_0_xxyzzz_xxxyy, g_z_0_xxyzzz_xxxyz, g_z_0_xxyzzz_xxxzz, g_z_0_xxyzzz_xxyyy, g_z_0_xxyzzz_xxyyz, g_z_0_xxyzzz_xxyzz, g_z_0_xxyzzz_xxzzz, g_z_0_xxyzzz_xyyyy, g_z_0_xxyzzz_xyyyz, g_z_0_xxyzzz_xyyzz, g_z_0_xxyzzz_xyzzz, g_z_0_xxyzzz_xzzzz, g_z_0_xxyzzz_yyyyy, g_z_0_xxyzzz_yyyyz, g_z_0_xxyzzz_yyyzz, g_z_0_xxyzzz_yyzzz, g_z_0_xxyzzz_yzzzz, g_z_0_xxyzzz_zzzzz, g_z_0_xyzzz_xxxxx, g_z_0_xyzzz_xxxxxx, g_z_0_xyzzz_xxxxxy, g_z_0_xyzzz_xxxxxz, g_z_0_xyzzz_xxxxy, g_z_0_xyzzz_xxxxyy, g_z_0_xyzzz_xxxxyz, g_z_0_xyzzz_xxxxz, g_z_0_xyzzz_xxxxzz, g_z_0_xyzzz_xxxyy, g_z_0_xyzzz_xxxyyy, g_z_0_xyzzz_xxxyyz, g_z_0_xyzzz_xxxyz, g_z_0_xyzzz_xxxyzz, g_z_0_xyzzz_xxxzz, g_z_0_xyzzz_xxxzzz, g_z_0_xyzzz_xxyyy, g_z_0_xyzzz_xxyyyy, g_z_0_xyzzz_xxyyyz, g_z_0_xyzzz_xxyyz, g_z_0_xyzzz_xxyyzz, g_z_0_xyzzz_xxyzz, g_z_0_xyzzz_xxyzzz, g_z_0_xyzzz_xxzzz, g_z_0_xyzzz_xxzzzz, g_z_0_xyzzz_xyyyy, g_z_0_xyzzz_xyyyyy, g_z_0_xyzzz_xyyyyz, g_z_0_xyzzz_xyyyz, g_z_0_xyzzz_xyyyzz, g_z_0_xyzzz_xyyzz, g_z_0_xyzzz_xyyzzz, g_z_0_xyzzz_xyzzz, g_z_0_xyzzz_xyzzzz, g_z_0_xyzzz_xzzzz, g_z_0_xyzzz_xzzzzz, g_z_0_xyzzz_yyyyy, g_z_0_xyzzz_yyyyz, g_z_0_xyzzz_yyyzz, g_z_0_xyzzz_yyzzz, g_z_0_xyzzz_yzzzz, g_z_0_xyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_xxxxx[k] = -g_z_0_xyzzz_xxxxx[k] * cd_x[k] + g_z_0_xyzzz_xxxxxx[k];

                g_z_0_xxyzzz_xxxxy[k] = -g_z_0_xyzzz_xxxxy[k] * cd_x[k] + g_z_0_xyzzz_xxxxxy[k];

                g_z_0_xxyzzz_xxxxz[k] = -g_z_0_xyzzz_xxxxz[k] * cd_x[k] + g_z_0_xyzzz_xxxxxz[k];

                g_z_0_xxyzzz_xxxyy[k] = -g_z_0_xyzzz_xxxyy[k] * cd_x[k] + g_z_0_xyzzz_xxxxyy[k];

                g_z_0_xxyzzz_xxxyz[k] = -g_z_0_xyzzz_xxxyz[k] * cd_x[k] + g_z_0_xyzzz_xxxxyz[k];

                g_z_0_xxyzzz_xxxzz[k] = -g_z_0_xyzzz_xxxzz[k] * cd_x[k] + g_z_0_xyzzz_xxxxzz[k];

                g_z_0_xxyzzz_xxyyy[k] = -g_z_0_xyzzz_xxyyy[k] * cd_x[k] + g_z_0_xyzzz_xxxyyy[k];

                g_z_0_xxyzzz_xxyyz[k] = -g_z_0_xyzzz_xxyyz[k] * cd_x[k] + g_z_0_xyzzz_xxxyyz[k];

                g_z_0_xxyzzz_xxyzz[k] = -g_z_0_xyzzz_xxyzz[k] * cd_x[k] + g_z_0_xyzzz_xxxyzz[k];

                g_z_0_xxyzzz_xxzzz[k] = -g_z_0_xyzzz_xxzzz[k] * cd_x[k] + g_z_0_xyzzz_xxxzzz[k];

                g_z_0_xxyzzz_xyyyy[k] = -g_z_0_xyzzz_xyyyy[k] * cd_x[k] + g_z_0_xyzzz_xxyyyy[k];

                g_z_0_xxyzzz_xyyyz[k] = -g_z_0_xyzzz_xyyyz[k] * cd_x[k] + g_z_0_xyzzz_xxyyyz[k];

                g_z_0_xxyzzz_xyyzz[k] = -g_z_0_xyzzz_xyyzz[k] * cd_x[k] + g_z_0_xyzzz_xxyyzz[k];

                g_z_0_xxyzzz_xyzzz[k] = -g_z_0_xyzzz_xyzzz[k] * cd_x[k] + g_z_0_xyzzz_xxyzzz[k];

                g_z_0_xxyzzz_xzzzz[k] = -g_z_0_xyzzz_xzzzz[k] * cd_x[k] + g_z_0_xyzzz_xxzzzz[k];

                g_z_0_xxyzzz_yyyyy[k] = -g_z_0_xyzzz_yyyyy[k] * cd_x[k] + g_z_0_xyzzz_xyyyyy[k];

                g_z_0_xxyzzz_yyyyz[k] = -g_z_0_xyzzz_yyyyz[k] * cd_x[k] + g_z_0_xyzzz_xyyyyz[k];

                g_z_0_xxyzzz_yyyzz[k] = -g_z_0_xyzzz_yyyzz[k] * cd_x[k] + g_z_0_xyzzz_xyyyzz[k];

                g_z_0_xxyzzz_yyzzz[k] = -g_z_0_xyzzz_yyzzz[k] * cd_x[k] + g_z_0_xyzzz_xyyzzz[k];

                g_z_0_xxyzzz_yzzzz[k] = -g_z_0_xyzzz_yzzzz[k] * cd_x[k] + g_z_0_xyzzz_xyzzzz[k];

                g_z_0_xxyzzz_zzzzz[k] = -g_z_0_xyzzz_zzzzz[k] * cd_x[k] + g_z_0_xyzzz_xzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 294);

            auto g_z_0_xxzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 295);

            auto g_z_0_xxzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 296);

            auto g_z_0_xxzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 297);

            auto g_z_0_xxzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 298);

            auto g_z_0_xxzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 299);

            auto g_z_0_xxzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 300);

            auto g_z_0_xxzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 301);

            auto g_z_0_xxzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 302);

            auto g_z_0_xxzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 303);

            auto g_z_0_xxzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 304);

            auto g_z_0_xxzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 305);

            auto g_z_0_xxzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 306);

            auto g_z_0_xxzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 307);

            auto g_z_0_xxzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 308);

            auto g_z_0_xxzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 309);

            auto g_z_0_xxzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 310);

            auto g_z_0_xxzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 311);

            auto g_z_0_xxzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 312);

            auto g_z_0_xxzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 313);

            auto g_z_0_xxzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzzz_xxxxx, g_z_0_xxzzzz_xxxxy, g_z_0_xxzzzz_xxxxz, g_z_0_xxzzzz_xxxyy, g_z_0_xxzzzz_xxxyz, g_z_0_xxzzzz_xxxzz, g_z_0_xxzzzz_xxyyy, g_z_0_xxzzzz_xxyyz, g_z_0_xxzzzz_xxyzz, g_z_0_xxzzzz_xxzzz, g_z_0_xxzzzz_xyyyy, g_z_0_xxzzzz_xyyyz, g_z_0_xxzzzz_xyyzz, g_z_0_xxzzzz_xyzzz, g_z_0_xxzzzz_xzzzz, g_z_0_xxzzzz_yyyyy, g_z_0_xxzzzz_yyyyz, g_z_0_xxzzzz_yyyzz, g_z_0_xxzzzz_yyzzz, g_z_0_xxzzzz_yzzzz, g_z_0_xxzzzz_zzzzz, g_z_0_xzzzz_xxxxx, g_z_0_xzzzz_xxxxxx, g_z_0_xzzzz_xxxxxy, g_z_0_xzzzz_xxxxxz, g_z_0_xzzzz_xxxxy, g_z_0_xzzzz_xxxxyy, g_z_0_xzzzz_xxxxyz, g_z_0_xzzzz_xxxxz, g_z_0_xzzzz_xxxxzz, g_z_0_xzzzz_xxxyy, g_z_0_xzzzz_xxxyyy, g_z_0_xzzzz_xxxyyz, g_z_0_xzzzz_xxxyz, g_z_0_xzzzz_xxxyzz, g_z_0_xzzzz_xxxzz, g_z_0_xzzzz_xxxzzz, g_z_0_xzzzz_xxyyy, g_z_0_xzzzz_xxyyyy, g_z_0_xzzzz_xxyyyz, g_z_0_xzzzz_xxyyz, g_z_0_xzzzz_xxyyzz, g_z_0_xzzzz_xxyzz, g_z_0_xzzzz_xxyzzz, g_z_0_xzzzz_xxzzz, g_z_0_xzzzz_xxzzzz, g_z_0_xzzzz_xyyyy, g_z_0_xzzzz_xyyyyy, g_z_0_xzzzz_xyyyyz, g_z_0_xzzzz_xyyyz, g_z_0_xzzzz_xyyyzz, g_z_0_xzzzz_xyyzz, g_z_0_xzzzz_xyyzzz, g_z_0_xzzzz_xyzzz, g_z_0_xzzzz_xyzzzz, g_z_0_xzzzz_xzzzz, g_z_0_xzzzz_xzzzzz, g_z_0_xzzzz_yyyyy, g_z_0_xzzzz_yyyyz, g_z_0_xzzzz_yyyzz, g_z_0_xzzzz_yyzzz, g_z_0_xzzzz_yzzzz, g_z_0_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_xxxxx[k] = -g_z_0_xzzzz_xxxxx[k] * cd_x[k] + g_z_0_xzzzz_xxxxxx[k];

                g_z_0_xxzzzz_xxxxy[k] = -g_z_0_xzzzz_xxxxy[k] * cd_x[k] + g_z_0_xzzzz_xxxxxy[k];

                g_z_0_xxzzzz_xxxxz[k] = -g_z_0_xzzzz_xxxxz[k] * cd_x[k] + g_z_0_xzzzz_xxxxxz[k];

                g_z_0_xxzzzz_xxxyy[k] = -g_z_0_xzzzz_xxxyy[k] * cd_x[k] + g_z_0_xzzzz_xxxxyy[k];

                g_z_0_xxzzzz_xxxyz[k] = -g_z_0_xzzzz_xxxyz[k] * cd_x[k] + g_z_0_xzzzz_xxxxyz[k];

                g_z_0_xxzzzz_xxxzz[k] = -g_z_0_xzzzz_xxxzz[k] * cd_x[k] + g_z_0_xzzzz_xxxxzz[k];

                g_z_0_xxzzzz_xxyyy[k] = -g_z_0_xzzzz_xxyyy[k] * cd_x[k] + g_z_0_xzzzz_xxxyyy[k];

                g_z_0_xxzzzz_xxyyz[k] = -g_z_0_xzzzz_xxyyz[k] * cd_x[k] + g_z_0_xzzzz_xxxyyz[k];

                g_z_0_xxzzzz_xxyzz[k] = -g_z_0_xzzzz_xxyzz[k] * cd_x[k] + g_z_0_xzzzz_xxxyzz[k];

                g_z_0_xxzzzz_xxzzz[k] = -g_z_0_xzzzz_xxzzz[k] * cd_x[k] + g_z_0_xzzzz_xxxzzz[k];

                g_z_0_xxzzzz_xyyyy[k] = -g_z_0_xzzzz_xyyyy[k] * cd_x[k] + g_z_0_xzzzz_xxyyyy[k];

                g_z_0_xxzzzz_xyyyz[k] = -g_z_0_xzzzz_xyyyz[k] * cd_x[k] + g_z_0_xzzzz_xxyyyz[k];

                g_z_0_xxzzzz_xyyzz[k] = -g_z_0_xzzzz_xyyzz[k] * cd_x[k] + g_z_0_xzzzz_xxyyzz[k];

                g_z_0_xxzzzz_xyzzz[k] = -g_z_0_xzzzz_xyzzz[k] * cd_x[k] + g_z_0_xzzzz_xxyzzz[k];

                g_z_0_xxzzzz_xzzzz[k] = -g_z_0_xzzzz_xzzzz[k] * cd_x[k] + g_z_0_xzzzz_xxzzzz[k];

                g_z_0_xxzzzz_yyyyy[k] = -g_z_0_xzzzz_yyyyy[k] * cd_x[k] + g_z_0_xzzzz_xyyyyy[k];

                g_z_0_xxzzzz_yyyyz[k] = -g_z_0_xzzzz_yyyyz[k] * cd_x[k] + g_z_0_xzzzz_xyyyyz[k];

                g_z_0_xxzzzz_yyyzz[k] = -g_z_0_xzzzz_yyyzz[k] * cd_x[k] + g_z_0_xzzzz_xyyyzz[k];

                g_z_0_xxzzzz_yyzzz[k] = -g_z_0_xzzzz_yyzzz[k] * cd_x[k] + g_z_0_xzzzz_xyyzzz[k];

                g_z_0_xxzzzz_yzzzz[k] = -g_z_0_xzzzz_yzzzz[k] * cd_x[k] + g_z_0_xzzzz_xyzzzz[k];

                g_z_0_xxzzzz_zzzzz[k] = -g_z_0_xzzzz_zzzzz[k] * cd_x[k] + g_z_0_xzzzz_xzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 315);

            auto g_z_0_xyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 316);

            auto g_z_0_xyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 317);

            auto g_z_0_xyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 318);

            auto g_z_0_xyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 319);

            auto g_z_0_xyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 320);

            auto g_z_0_xyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 321);

            auto g_z_0_xyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 322);

            auto g_z_0_xyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 323);

            auto g_z_0_xyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 324);

            auto g_z_0_xyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 325);

            auto g_z_0_xyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 326);

            auto g_z_0_xyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 327);

            auto g_z_0_xyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 328);

            auto g_z_0_xyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 329);

            auto g_z_0_xyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 330);

            auto g_z_0_xyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 331);

            auto g_z_0_xyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 332);

            auto g_z_0_xyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 333);

            auto g_z_0_xyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 334);

            auto g_z_0_xyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyy_xxxxx, g_z_0_xyyyyy_xxxxy, g_z_0_xyyyyy_xxxxz, g_z_0_xyyyyy_xxxyy, g_z_0_xyyyyy_xxxyz, g_z_0_xyyyyy_xxxzz, g_z_0_xyyyyy_xxyyy, g_z_0_xyyyyy_xxyyz, g_z_0_xyyyyy_xxyzz, g_z_0_xyyyyy_xxzzz, g_z_0_xyyyyy_xyyyy, g_z_0_xyyyyy_xyyyz, g_z_0_xyyyyy_xyyzz, g_z_0_xyyyyy_xyzzz, g_z_0_xyyyyy_xzzzz, g_z_0_xyyyyy_yyyyy, g_z_0_xyyyyy_yyyyz, g_z_0_xyyyyy_yyyzz, g_z_0_xyyyyy_yyzzz, g_z_0_xyyyyy_yzzzz, g_z_0_xyyyyy_zzzzz, g_z_0_yyyyy_xxxxx, g_z_0_yyyyy_xxxxxx, g_z_0_yyyyy_xxxxxy, g_z_0_yyyyy_xxxxxz, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxxyy, g_z_0_yyyyy_xxxxyz, g_z_0_yyyyy_xxxxz, g_z_0_yyyyy_xxxxzz, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyyy, g_z_0_yyyyy_xxxyyz, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxyzz, g_z_0_yyyyy_xxxzz, g_z_0_yyyyy_xxxzzz, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyyy, g_z_0_yyyyy_xxyyyz, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyyzz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxyzzz, g_z_0_yyyyy_xxzzz, g_z_0_yyyyy_xxzzzz, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyyy, g_z_0_yyyyy_xyyyyz, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyyzz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyyzzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xyzzzz, g_z_0_yyyyy_xzzzz, g_z_0_yyyyy_xzzzzz, g_z_0_yyyyy_yyyyy, g_z_0_yyyyy_yyyyz, g_z_0_yyyyy_yyyzz, g_z_0_yyyyy_yyzzz, g_z_0_yyyyy_yzzzz, g_z_0_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_xxxxx[k] = -g_z_0_yyyyy_xxxxx[k] * cd_x[k] + g_z_0_yyyyy_xxxxxx[k];

                g_z_0_xyyyyy_xxxxy[k] = -g_z_0_yyyyy_xxxxy[k] * cd_x[k] + g_z_0_yyyyy_xxxxxy[k];

                g_z_0_xyyyyy_xxxxz[k] = -g_z_0_yyyyy_xxxxz[k] * cd_x[k] + g_z_0_yyyyy_xxxxxz[k];

                g_z_0_xyyyyy_xxxyy[k] = -g_z_0_yyyyy_xxxyy[k] * cd_x[k] + g_z_0_yyyyy_xxxxyy[k];

                g_z_0_xyyyyy_xxxyz[k] = -g_z_0_yyyyy_xxxyz[k] * cd_x[k] + g_z_0_yyyyy_xxxxyz[k];

                g_z_0_xyyyyy_xxxzz[k] = -g_z_0_yyyyy_xxxzz[k] * cd_x[k] + g_z_0_yyyyy_xxxxzz[k];

                g_z_0_xyyyyy_xxyyy[k] = -g_z_0_yyyyy_xxyyy[k] * cd_x[k] + g_z_0_yyyyy_xxxyyy[k];

                g_z_0_xyyyyy_xxyyz[k] = -g_z_0_yyyyy_xxyyz[k] * cd_x[k] + g_z_0_yyyyy_xxxyyz[k];

                g_z_0_xyyyyy_xxyzz[k] = -g_z_0_yyyyy_xxyzz[k] * cd_x[k] + g_z_0_yyyyy_xxxyzz[k];

                g_z_0_xyyyyy_xxzzz[k] = -g_z_0_yyyyy_xxzzz[k] * cd_x[k] + g_z_0_yyyyy_xxxzzz[k];

                g_z_0_xyyyyy_xyyyy[k] = -g_z_0_yyyyy_xyyyy[k] * cd_x[k] + g_z_0_yyyyy_xxyyyy[k];

                g_z_0_xyyyyy_xyyyz[k] = -g_z_0_yyyyy_xyyyz[k] * cd_x[k] + g_z_0_yyyyy_xxyyyz[k];

                g_z_0_xyyyyy_xyyzz[k] = -g_z_0_yyyyy_xyyzz[k] * cd_x[k] + g_z_0_yyyyy_xxyyzz[k];

                g_z_0_xyyyyy_xyzzz[k] = -g_z_0_yyyyy_xyzzz[k] * cd_x[k] + g_z_0_yyyyy_xxyzzz[k];

                g_z_0_xyyyyy_xzzzz[k] = -g_z_0_yyyyy_xzzzz[k] * cd_x[k] + g_z_0_yyyyy_xxzzzz[k];

                g_z_0_xyyyyy_yyyyy[k] = -g_z_0_yyyyy_yyyyy[k] * cd_x[k] + g_z_0_yyyyy_xyyyyy[k];

                g_z_0_xyyyyy_yyyyz[k] = -g_z_0_yyyyy_yyyyz[k] * cd_x[k] + g_z_0_yyyyy_xyyyyz[k];

                g_z_0_xyyyyy_yyyzz[k] = -g_z_0_yyyyy_yyyzz[k] * cd_x[k] + g_z_0_yyyyy_xyyyzz[k];

                g_z_0_xyyyyy_yyzzz[k] = -g_z_0_yyyyy_yyzzz[k] * cd_x[k] + g_z_0_yyyyy_xyyzzz[k];

                g_z_0_xyyyyy_yzzzz[k] = -g_z_0_yyyyy_yzzzz[k] * cd_x[k] + g_z_0_yyyyy_xyzzzz[k];

                g_z_0_xyyyyy_zzzzz[k] = -g_z_0_yyyyy_zzzzz[k] * cd_x[k] + g_z_0_yyyyy_xzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 336);

            auto g_z_0_xyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 337);

            auto g_z_0_xyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 338);

            auto g_z_0_xyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 339);

            auto g_z_0_xyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 340);

            auto g_z_0_xyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 341);

            auto g_z_0_xyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 342);

            auto g_z_0_xyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 343);

            auto g_z_0_xyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 344);

            auto g_z_0_xyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 345);

            auto g_z_0_xyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 346);

            auto g_z_0_xyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 347);

            auto g_z_0_xyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 348);

            auto g_z_0_xyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 349);

            auto g_z_0_xyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 350);

            auto g_z_0_xyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 351);

            auto g_z_0_xyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 352);

            auto g_z_0_xyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 353);

            auto g_z_0_xyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 354);

            auto g_z_0_xyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 355);

            auto g_z_0_xyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 356);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyz_xxxxx, g_z_0_xyyyyz_xxxxy, g_z_0_xyyyyz_xxxxz, g_z_0_xyyyyz_xxxyy, g_z_0_xyyyyz_xxxyz, g_z_0_xyyyyz_xxxzz, g_z_0_xyyyyz_xxyyy, g_z_0_xyyyyz_xxyyz, g_z_0_xyyyyz_xxyzz, g_z_0_xyyyyz_xxzzz, g_z_0_xyyyyz_xyyyy, g_z_0_xyyyyz_xyyyz, g_z_0_xyyyyz_xyyzz, g_z_0_xyyyyz_xyzzz, g_z_0_xyyyyz_xzzzz, g_z_0_xyyyyz_yyyyy, g_z_0_xyyyyz_yyyyz, g_z_0_xyyyyz_yyyzz, g_z_0_xyyyyz_yyzzz, g_z_0_xyyyyz_yzzzz, g_z_0_xyyyyz_zzzzz, g_z_0_yyyyz_xxxxx, g_z_0_yyyyz_xxxxxx, g_z_0_yyyyz_xxxxxy, g_z_0_yyyyz_xxxxxz, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxxyy, g_z_0_yyyyz_xxxxyz, g_z_0_yyyyz_xxxxz, g_z_0_yyyyz_xxxxzz, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyyy, g_z_0_yyyyz_xxxyyz, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxyzz, g_z_0_yyyyz_xxxzz, g_z_0_yyyyz_xxxzzz, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyyy, g_z_0_yyyyz_xxyyyz, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyyzz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxyzzz, g_z_0_yyyyz_xxzzz, g_z_0_yyyyz_xxzzzz, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyyy, g_z_0_yyyyz_xyyyyz, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyyzz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyyzzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xyzzzz, g_z_0_yyyyz_xzzzz, g_z_0_yyyyz_xzzzzz, g_z_0_yyyyz_yyyyy, g_z_0_yyyyz_yyyyz, g_z_0_yyyyz_yyyzz, g_z_0_yyyyz_yyzzz, g_z_0_yyyyz_yzzzz, g_z_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_xxxxx[k] = -g_z_0_yyyyz_xxxxx[k] * cd_x[k] + g_z_0_yyyyz_xxxxxx[k];

                g_z_0_xyyyyz_xxxxy[k] = -g_z_0_yyyyz_xxxxy[k] * cd_x[k] + g_z_0_yyyyz_xxxxxy[k];

                g_z_0_xyyyyz_xxxxz[k] = -g_z_0_yyyyz_xxxxz[k] * cd_x[k] + g_z_0_yyyyz_xxxxxz[k];

                g_z_0_xyyyyz_xxxyy[k] = -g_z_0_yyyyz_xxxyy[k] * cd_x[k] + g_z_0_yyyyz_xxxxyy[k];

                g_z_0_xyyyyz_xxxyz[k] = -g_z_0_yyyyz_xxxyz[k] * cd_x[k] + g_z_0_yyyyz_xxxxyz[k];

                g_z_0_xyyyyz_xxxzz[k] = -g_z_0_yyyyz_xxxzz[k] * cd_x[k] + g_z_0_yyyyz_xxxxzz[k];

                g_z_0_xyyyyz_xxyyy[k] = -g_z_0_yyyyz_xxyyy[k] * cd_x[k] + g_z_0_yyyyz_xxxyyy[k];

                g_z_0_xyyyyz_xxyyz[k] = -g_z_0_yyyyz_xxyyz[k] * cd_x[k] + g_z_0_yyyyz_xxxyyz[k];

                g_z_0_xyyyyz_xxyzz[k] = -g_z_0_yyyyz_xxyzz[k] * cd_x[k] + g_z_0_yyyyz_xxxyzz[k];

                g_z_0_xyyyyz_xxzzz[k] = -g_z_0_yyyyz_xxzzz[k] * cd_x[k] + g_z_0_yyyyz_xxxzzz[k];

                g_z_0_xyyyyz_xyyyy[k] = -g_z_0_yyyyz_xyyyy[k] * cd_x[k] + g_z_0_yyyyz_xxyyyy[k];

                g_z_0_xyyyyz_xyyyz[k] = -g_z_0_yyyyz_xyyyz[k] * cd_x[k] + g_z_0_yyyyz_xxyyyz[k];

                g_z_0_xyyyyz_xyyzz[k] = -g_z_0_yyyyz_xyyzz[k] * cd_x[k] + g_z_0_yyyyz_xxyyzz[k];

                g_z_0_xyyyyz_xyzzz[k] = -g_z_0_yyyyz_xyzzz[k] * cd_x[k] + g_z_0_yyyyz_xxyzzz[k];

                g_z_0_xyyyyz_xzzzz[k] = -g_z_0_yyyyz_xzzzz[k] * cd_x[k] + g_z_0_yyyyz_xxzzzz[k];

                g_z_0_xyyyyz_yyyyy[k] = -g_z_0_yyyyz_yyyyy[k] * cd_x[k] + g_z_0_yyyyz_xyyyyy[k];

                g_z_0_xyyyyz_yyyyz[k] = -g_z_0_yyyyz_yyyyz[k] * cd_x[k] + g_z_0_yyyyz_xyyyyz[k];

                g_z_0_xyyyyz_yyyzz[k] = -g_z_0_yyyyz_yyyzz[k] * cd_x[k] + g_z_0_yyyyz_xyyyzz[k];

                g_z_0_xyyyyz_yyzzz[k] = -g_z_0_yyyyz_yyzzz[k] * cd_x[k] + g_z_0_yyyyz_xyyzzz[k];

                g_z_0_xyyyyz_yzzzz[k] = -g_z_0_yyyyz_yzzzz[k] * cd_x[k] + g_z_0_yyyyz_xyzzzz[k];

                g_z_0_xyyyyz_zzzzz[k] = -g_z_0_yyyyz_zzzzz[k] * cd_x[k] + g_z_0_yyyyz_xzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 357);

            auto g_z_0_xyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 358);

            auto g_z_0_xyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 359);

            auto g_z_0_xyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 360);

            auto g_z_0_xyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 361);

            auto g_z_0_xyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 362);

            auto g_z_0_xyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 363);

            auto g_z_0_xyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 364);

            auto g_z_0_xyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 365);

            auto g_z_0_xyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 366);

            auto g_z_0_xyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 367);

            auto g_z_0_xyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 368);

            auto g_z_0_xyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 369);

            auto g_z_0_xyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 370);

            auto g_z_0_xyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 371);

            auto g_z_0_xyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 372);

            auto g_z_0_xyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 373);

            auto g_z_0_xyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 374);

            auto g_z_0_xyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 375);

            auto g_z_0_xyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 376);

            auto g_z_0_xyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 377);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyzz_xxxxx, g_z_0_xyyyzz_xxxxy, g_z_0_xyyyzz_xxxxz, g_z_0_xyyyzz_xxxyy, g_z_0_xyyyzz_xxxyz, g_z_0_xyyyzz_xxxzz, g_z_0_xyyyzz_xxyyy, g_z_0_xyyyzz_xxyyz, g_z_0_xyyyzz_xxyzz, g_z_0_xyyyzz_xxzzz, g_z_0_xyyyzz_xyyyy, g_z_0_xyyyzz_xyyyz, g_z_0_xyyyzz_xyyzz, g_z_0_xyyyzz_xyzzz, g_z_0_xyyyzz_xzzzz, g_z_0_xyyyzz_yyyyy, g_z_0_xyyyzz_yyyyz, g_z_0_xyyyzz_yyyzz, g_z_0_xyyyzz_yyzzz, g_z_0_xyyyzz_yzzzz, g_z_0_xyyyzz_zzzzz, g_z_0_yyyzz_xxxxx, g_z_0_yyyzz_xxxxxx, g_z_0_yyyzz_xxxxxy, g_z_0_yyyzz_xxxxxz, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxxyy, g_z_0_yyyzz_xxxxyz, g_z_0_yyyzz_xxxxz, g_z_0_yyyzz_xxxxzz, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyyy, g_z_0_yyyzz_xxxyyz, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxyzz, g_z_0_yyyzz_xxxzz, g_z_0_yyyzz_xxxzzz, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyyy, g_z_0_yyyzz_xxyyyz, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyyzz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxyzzz, g_z_0_yyyzz_xxzzz, g_z_0_yyyzz_xxzzzz, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyyy, g_z_0_yyyzz_xyyyyz, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyyzz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyyzzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xyzzzz, g_z_0_yyyzz_xzzzz, g_z_0_yyyzz_xzzzzz, g_z_0_yyyzz_yyyyy, g_z_0_yyyzz_yyyyz, g_z_0_yyyzz_yyyzz, g_z_0_yyyzz_yyzzz, g_z_0_yyyzz_yzzzz, g_z_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_xxxxx[k] = -g_z_0_yyyzz_xxxxx[k] * cd_x[k] + g_z_0_yyyzz_xxxxxx[k];

                g_z_0_xyyyzz_xxxxy[k] = -g_z_0_yyyzz_xxxxy[k] * cd_x[k] + g_z_0_yyyzz_xxxxxy[k];

                g_z_0_xyyyzz_xxxxz[k] = -g_z_0_yyyzz_xxxxz[k] * cd_x[k] + g_z_0_yyyzz_xxxxxz[k];

                g_z_0_xyyyzz_xxxyy[k] = -g_z_0_yyyzz_xxxyy[k] * cd_x[k] + g_z_0_yyyzz_xxxxyy[k];

                g_z_0_xyyyzz_xxxyz[k] = -g_z_0_yyyzz_xxxyz[k] * cd_x[k] + g_z_0_yyyzz_xxxxyz[k];

                g_z_0_xyyyzz_xxxzz[k] = -g_z_0_yyyzz_xxxzz[k] * cd_x[k] + g_z_0_yyyzz_xxxxzz[k];

                g_z_0_xyyyzz_xxyyy[k] = -g_z_0_yyyzz_xxyyy[k] * cd_x[k] + g_z_0_yyyzz_xxxyyy[k];

                g_z_0_xyyyzz_xxyyz[k] = -g_z_0_yyyzz_xxyyz[k] * cd_x[k] + g_z_0_yyyzz_xxxyyz[k];

                g_z_0_xyyyzz_xxyzz[k] = -g_z_0_yyyzz_xxyzz[k] * cd_x[k] + g_z_0_yyyzz_xxxyzz[k];

                g_z_0_xyyyzz_xxzzz[k] = -g_z_0_yyyzz_xxzzz[k] * cd_x[k] + g_z_0_yyyzz_xxxzzz[k];

                g_z_0_xyyyzz_xyyyy[k] = -g_z_0_yyyzz_xyyyy[k] * cd_x[k] + g_z_0_yyyzz_xxyyyy[k];

                g_z_0_xyyyzz_xyyyz[k] = -g_z_0_yyyzz_xyyyz[k] * cd_x[k] + g_z_0_yyyzz_xxyyyz[k];

                g_z_0_xyyyzz_xyyzz[k] = -g_z_0_yyyzz_xyyzz[k] * cd_x[k] + g_z_0_yyyzz_xxyyzz[k];

                g_z_0_xyyyzz_xyzzz[k] = -g_z_0_yyyzz_xyzzz[k] * cd_x[k] + g_z_0_yyyzz_xxyzzz[k];

                g_z_0_xyyyzz_xzzzz[k] = -g_z_0_yyyzz_xzzzz[k] * cd_x[k] + g_z_0_yyyzz_xxzzzz[k];

                g_z_0_xyyyzz_yyyyy[k] = -g_z_0_yyyzz_yyyyy[k] * cd_x[k] + g_z_0_yyyzz_xyyyyy[k];

                g_z_0_xyyyzz_yyyyz[k] = -g_z_0_yyyzz_yyyyz[k] * cd_x[k] + g_z_0_yyyzz_xyyyyz[k];

                g_z_0_xyyyzz_yyyzz[k] = -g_z_0_yyyzz_yyyzz[k] * cd_x[k] + g_z_0_yyyzz_xyyyzz[k];

                g_z_0_xyyyzz_yyzzz[k] = -g_z_0_yyyzz_yyzzz[k] * cd_x[k] + g_z_0_yyyzz_xyyzzz[k];

                g_z_0_xyyyzz_yzzzz[k] = -g_z_0_yyyzz_yzzzz[k] * cd_x[k] + g_z_0_yyyzz_xyzzzz[k];

                g_z_0_xyyyzz_zzzzz[k] = -g_z_0_yyyzz_zzzzz[k] * cd_x[k] + g_z_0_yyyzz_xzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 378);

            auto g_z_0_xyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 379);

            auto g_z_0_xyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 380);

            auto g_z_0_xyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 381);

            auto g_z_0_xyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 382);

            auto g_z_0_xyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 383);

            auto g_z_0_xyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 384);

            auto g_z_0_xyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 385);

            auto g_z_0_xyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 386);

            auto g_z_0_xyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 387);

            auto g_z_0_xyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 388);

            auto g_z_0_xyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 389);

            auto g_z_0_xyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 390);

            auto g_z_0_xyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 391);

            auto g_z_0_xyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 392);

            auto g_z_0_xyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 393);

            auto g_z_0_xyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 394);

            auto g_z_0_xyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 395);

            auto g_z_0_xyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 396);

            auto g_z_0_xyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 397);

            auto g_z_0_xyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 398);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzzz_xxxxx, g_z_0_xyyzzz_xxxxy, g_z_0_xyyzzz_xxxxz, g_z_0_xyyzzz_xxxyy, g_z_0_xyyzzz_xxxyz, g_z_0_xyyzzz_xxxzz, g_z_0_xyyzzz_xxyyy, g_z_0_xyyzzz_xxyyz, g_z_0_xyyzzz_xxyzz, g_z_0_xyyzzz_xxzzz, g_z_0_xyyzzz_xyyyy, g_z_0_xyyzzz_xyyyz, g_z_0_xyyzzz_xyyzz, g_z_0_xyyzzz_xyzzz, g_z_0_xyyzzz_xzzzz, g_z_0_xyyzzz_yyyyy, g_z_0_xyyzzz_yyyyz, g_z_0_xyyzzz_yyyzz, g_z_0_xyyzzz_yyzzz, g_z_0_xyyzzz_yzzzz, g_z_0_xyyzzz_zzzzz, g_z_0_yyzzz_xxxxx, g_z_0_yyzzz_xxxxxx, g_z_0_yyzzz_xxxxxy, g_z_0_yyzzz_xxxxxz, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxxyy, g_z_0_yyzzz_xxxxyz, g_z_0_yyzzz_xxxxz, g_z_0_yyzzz_xxxxzz, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyyy, g_z_0_yyzzz_xxxyyz, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxyzz, g_z_0_yyzzz_xxxzz, g_z_0_yyzzz_xxxzzz, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyyy, g_z_0_yyzzz_xxyyyz, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyyzz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxyzzz, g_z_0_yyzzz_xxzzz, g_z_0_yyzzz_xxzzzz, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyyy, g_z_0_yyzzz_xyyyyz, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyyzz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyyzzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xyzzzz, g_z_0_yyzzz_xzzzz, g_z_0_yyzzz_xzzzzz, g_z_0_yyzzz_yyyyy, g_z_0_yyzzz_yyyyz, g_z_0_yyzzz_yyyzz, g_z_0_yyzzz_yyzzz, g_z_0_yyzzz_yzzzz, g_z_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_xxxxx[k] = -g_z_0_yyzzz_xxxxx[k] * cd_x[k] + g_z_0_yyzzz_xxxxxx[k];

                g_z_0_xyyzzz_xxxxy[k] = -g_z_0_yyzzz_xxxxy[k] * cd_x[k] + g_z_0_yyzzz_xxxxxy[k];

                g_z_0_xyyzzz_xxxxz[k] = -g_z_0_yyzzz_xxxxz[k] * cd_x[k] + g_z_0_yyzzz_xxxxxz[k];

                g_z_0_xyyzzz_xxxyy[k] = -g_z_0_yyzzz_xxxyy[k] * cd_x[k] + g_z_0_yyzzz_xxxxyy[k];

                g_z_0_xyyzzz_xxxyz[k] = -g_z_0_yyzzz_xxxyz[k] * cd_x[k] + g_z_0_yyzzz_xxxxyz[k];

                g_z_0_xyyzzz_xxxzz[k] = -g_z_0_yyzzz_xxxzz[k] * cd_x[k] + g_z_0_yyzzz_xxxxzz[k];

                g_z_0_xyyzzz_xxyyy[k] = -g_z_0_yyzzz_xxyyy[k] * cd_x[k] + g_z_0_yyzzz_xxxyyy[k];

                g_z_0_xyyzzz_xxyyz[k] = -g_z_0_yyzzz_xxyyz[k] * cd_x[k] + g_z_0_yyzzz_xxxyyz[k];

                g_z_0_xyyzzz_xxyzz[k] = -g_z_0_yyzzz_xxyzz[k] * cd_x[k] + g_z_0_yyzzz_xxxyzz[k];

                g_z_0_xyyzzz_xxzzz[k] = -g_z_0_yyzzz_xxzzz[k] * cd_x[k] + g_z_0_yyzzz_xxxzzz[k];

                g_z_0_xyyzzz_xyyyy[k] = -g_z_0_yyzzz_xyyyy[k] * cd_x[k] + g_z_0_yyzzz_xxyyyy[k];

                g_z_0_xyyzzz_xyyyz[k] = -g_z_0_yyzzz_xyyyz[k] * cd_x[k] + g_z_0_yyzzz_xxyyyz[k];

                g_z_0_xyyzzz_xyyzz[k] = -g_z_0_yyzzz_xyyzz[k] * cd_x[k] + g_z_0_yyzzz_xxyyzz[k];

                g_z_0_xyyzzz_xyzzz[k] = -g_z_0_yyzzz_xyzzz[k] * cd_x[k] + g_z_0_yyzzz_xxyzzz[k];

                g_z_0_xyyzzz_xzzzz[k] = -g_z_0_yyzzz_xzzzz[k] * cd_x[k] + g_z_0_yyzzz_xxzzzz[k];

                g_z_0_xyyzzz_yyyyy[k] = -g_z_0_yyzzz_yyyyy[k] * cd_x[k] + g_z_0_yyzzz_xyyyyy[k];

                g_z_0_xyyzzz_yyyyz[k] = -g_z_0_yyzzz_yyyyz[k] * cd_x[k] + g_z_0_yyzzz_xyyyyz[k];

                g_z_0_xyyzzz_yyyzz[k] = -g_z_0_yyzzz_yyyzz[k] * cd_x[k] + g_z_0_yyzzz_xyyyzz[k];

                g_z_0_xyyzzz_yyzzz[k] = -g_z_0_yyzzz_yyzzz[k] * cd_x[k] + g_z_0_yyzzz_xyyzzz[k];

                g_z_0_xyyzzz_yzzzz[k] = -g_z_0_yyzzz_yzzzz[k] * cd_x[k] + g_z_0_yyzzz_xyzzzz[k];

                g_z_0_xyyzzz_zzzzz[k] = -g_z_0_yyzzz_zzzzz[k] * cd_x[k] + g_z_0_yyzzz_xzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 399);

            auto g_z_0_xyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 400);

            auto g_z_0_xyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 401);

            auto g_z_0_xyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 402);

            auto g_z_0_xyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 403);

            auto g_z_0_xyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 404);

            auto g_z_0_xyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 405);

            auto g_z_0_xyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 406);

            auto g_z_0_xyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 407);

            auto g_z_0_xyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 408);

            auto g_z_0_xyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 409);

            auto g_z_0_xyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 410);

            auto g_z_0_xyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 411);

            auto g_z_0_xyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 412);

            auto g_z_0_xyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 413);

            auto g_z_0_xyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 414);

            auto g_z_0_xyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 415);

            auto g_z_0_xyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 416);

            auto g_z_0_xyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 417);

            auto g_z_0_xyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 418);

            auto g_z_0_xyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzzz_xxxxx, g_z_0_xyzzzz_xxxxy, g_z_0_xyzzzz_xxxxz, g_z_0_xyzzzz_xxxyy, g_z_0_xyzzzz_xxxyz, g_z_0_xyzzzz_xxxzz, g_z_0_xyzzzz_xxyyy, g_z_0_xyzzzz_xxyyz, g_z_0_xyzzzz_xxyzz, g_z_0_xyzzzz_xxzzz, g_z_0_xyzzzz_xyyyy, g_z_0_xyzzzz_xyyyz, g_z_0_xyzzzz_xyyzz, g_z_0_xyzzzz_xyzzz, g_z_0_xyzzzz_xzzzz, g_z_0_xyzzzz_yyyyy, g_z_0_xyzzzz_yyyyz, g_z_0_xyzzzz_yyyzz, g_z_0_xyzzzz_yyzzz, g_z_0_xyzzzz_yzzzz, g_z_0_xyzzzz_zzzzz, g_z_0_yzzzz_xxxxx, g_z_0_yzzzz_xxxxxx, g_z_0_yzzzz_xxxxxy, g_z_0_yzzzz_xxxxxz, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxxyy, g_z_0_yzzzz_xxxxyz, g_z_0_yzzzz_xxxxz, g_z_0_yzzzz_xxxxzz, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyyy, g_z_0_yzzzz_xxxyyz, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxyzz, g_z_0_yzzzz_xxxzz, g_z_0_yzzzz_xxxzzz, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyyy, g_z_0_yzzzz_xxyyyz, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyyzz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxyzzz, g_z_0_yzzzz_xxzzz, g_z_0_yzzzz_xxzzzz, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyyy, g_z_0_yzzzz_xyyyyz, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyyzz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyyzzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xyzzzz, g_z_0_yzzzz_xzzzz, g_z_0_yzzzz_xzzzzz, g_z_0_yzzzz_yyyyy, g_z_0_yzzzz_yyyyz, g_z_0_yzzzz_yyyzz, g_z_0_yzzzz_yyzzz, g_z_0_yzzzz_yzzzz, g_z_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_xxxxx[k] = -g_z_0_yzzzz_xxxxx[k] * cd_x[k] + g_z_0_yzzzz_xxxxxx[k];

                g_z_0_xyzzzz_xxxxy[k] = -g_z_0_yzzzz_xxxxy[k] * cd_x[k] + g_z_0_yzzzz_xxxxxy[k];

                g_z_0_xyzzzz_xxxxz[k] = -g_z_0_yzzzz_xxxxz[k] * cd_x[k] + g_z_0_yzzzz_xxxxxz[k];

                g_z_0_xyzzzz_xxxyy[k] = -g_z_0_yzzzz_xxxyy[k] * cd_x[k] + g_z_0_yzzzz_xxxxyy[k];

                g_z_0_xyzzzz_xxxyz[k] = -g_z_0_yzzzz_xxxyz[k] * cd_x[k] + g_z_0_yzzzz_xxxxyz[k];

                g_z_0_xyzzzz_xxxzz[k] = -g_z_0_yzzzz_xxxzz[k] * cd_x[k] + g_z_0_yzzzz_xxxxzz[k];

                g_z_0_xyzzzz_xxyyy[k] = -g_z_0_yzzzz_xxyyy[k] * cd_x[k] + g_z_0_yzzzz_xxxyyy[k];

                g_z_0_xyzzzz_xxyyz[k] = -g_z_0_yzzzz_xxyyz[k] * cd_x[k] + g_z_0_yzzzz_xxxyyz[k];

                g_z_0_xyzzzz_xxyzz[k] = -g_z_0_yzzzz_xxyzz[k] * cd_x[k] + g_z_0_yzzzz_xxxyzz[k];

                g_z_0_xyzzzz_xxzzz[k] = -g_z_0_yzzzz_xxzzz[k] * cd_x[k] + g_z_0_yzzzz_xxxzzz[k];

                g_z_0_xyzzzz_xyyyy[k] = -g_z_0_yzzzz_xyyyy[k] * cd_x[k] + g_z_0_yzzzz_xxyyyy[k];

                g_z_0_xyzzzz_xyyyz[k] = -g_z_0_yzzzz_xyyyz[k] * cd_x[k] + g_z_0_yzzzz_xxyyyz[k];

                g_z_0_xyzzzz_xyyzz[k] = -g_z_0_yzzzz_xyyzz[k] * cd_x[k] + g_z_0_yzzzz_xxyyzz[k];

                g_z_0_xyzzzz_xyzzz[k] = -g_z_0_yzzzz_xyzzz[k] * cd_x[k] + g_z_0_yzzzz_xxyzzz[k];

                g_z_0_xyzzzz_xzzzz[k] = -g_z_0_yzzzz_xzzzz[k] * cd_x[k] + g_z_0_yzzzz_xxzzzz[k];

                g_z_0_xyzzzz_yyyyy[k] = -g_z_0_yzzzz_yyyyy[k] * cd_x[k] + g_z_0_yzzzz_xyyyyy[k];

                g_z_0_xyzzzz_yyyyz[k] = -g_z_0_yzzzz_yyyyz[k] * cd_x[k] + g_z_0_yzzzz_xyyyyz[k];

                g_z_0_xyzzzz_yyyzz[k] = -g_z_0_yzzzz_yyyzz[k] * cd_x[k] + g_z_0_yzzzz_xyyyzz[k];

                g_z_0_xyzzzz_yyzzz[k] = -g_z_0_yzzzz_yyzzz[k] * cd_x[k] + g_z_0_yzzzz_xyyzzz[k];

                g_z_0_xyzzzz_yzzzz[k] = -g_z_0_yzzzz_yzzzz[k] * cd_x[k] + g_z_0_yzzzz_xyzzzz[k];

                g_z_0_xyzzzz_zzzzz[k] = -g_z_0_yzzzz_zzzzz[k] * cd_x[k] + g_z_0_yzzzz_xzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 420);

            auto g_z_0_xzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 421);

            auto g_z_0_xzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 422);

            auto g_z_0_xzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 423);

            auto g_z_0_xzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 424);

            auto g_z_0_xzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 425);

            auto g_z_0_xzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 426);

            auto g_z_0_xzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 427);

            auto g_z_0_xzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 428);

            auto g_z_0_xzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 429);

            auto g_z_0_xzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 430);

            auto g_z_0_xzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 431);

            auto g_z_0_xzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 432);

            auto g_z_0_xzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 433);

            auto g_z_0_xzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 434);

            auto g_z_0_xzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 435);

            auto g_z_0_xzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 436);

            auto g_z_0_xzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 437);

            auto g_z_0_xzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 438);

            auto g_z_0_xzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 439);

            auto g_z_0_xzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 440);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzzz_xxxxx, g_z_0_xzzzzz_xxxxy, g_z_0_xzzzzz_xxxxz, g_z_0_xzzzzz_xxxyy, g_z_0_xzzzzz_xxxyz, g_z_0_xzzzzz_xxxzz, g_z_0_xzzzzz_xxyyy, g_z_0_xzzzzz_xxyyz, g_z_0_xzzzzz_xxyzz, g_z_0_xzzzzz_xxzzz, g_z_0_xzzzzz_xyyyy, g_z_0_xzzzzz_xyyyz, g_z_0_xzzzzz_xyyzz, g_z_0_xzzzzz_xyzzz, g_z_0_xzzzzz_xzzzz, g_z_0_xzzzzz_yyyyy, g_z_0_xzzzzz_yyyyz, g_z_0_xzzzzz_yyyzz, g_z_0_xzzzzz_yyzzz, g_z_0_xzzzzz_yzzzz, g_z_0_xzzzzz_zzzzz, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxxx, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_xxxxx[k] = -g_z_0_zzzzz_xxxxx[k] * cd_x[k] + g_z_0_zzzzz_xxxxxx[k];

                g_z_0_xzzzzz_xxxxy[k] = -g_z_0_zzzzz_xxxxy[k] * cd_x[k] + g_z_0_zzzzz_xxxxxy[k];

                g_z_0_xzzzzz_xxxxz[k] = -g_z_0_zzzzz_xxxxz[k] * cd_x[k] + g_z_0_zzzzz_xxxxxz[k];

                g_z_0_xzzzzz_xxxyy[k] = -g_z_0_zzzzz_xxxyy[k] * cd_x[k] + g_z_0_zzzzz_xxxxyy[k];

                g_z_0_xzzzzz_xxxyz[k] = -g_z_0_zzzzz_xxxyz[k] * cd_x[k] + g_z_0_zzzzz_xxxxyz[k];

                g_z_0_xzzzzz_xxxzz[k] = -g_z_0_zzzzz_xxxzz[k] * cd_x[k] + g_z_0_zzzzz_xxxxzz[k];

                g_z_0_xzzzzz_xxyyy[k] = -g_z_0_zzzzz_xxyyy[k] * cd_x[k] + g_z_0_zzzzz_xxxyyy[k];

                g_z_0_xzzzzz_xxyyz[k] = -g_z_0_zzzzz_xxyyz[k] * cd_x[k] + g_z_0_zzzzz_xxxyyz[k];

                g_z_0_xzzzzz_xxyzz[k] = -g_z_0_zzzzz_xxyzz[k] * cd_x[k] + g_z_0_zzzzz_xxxyzz[k];

                g_z_0_xzzzzz_xxzzz[k] = -g_z_0_zzzzz_xxzzz[k] * cd_x[k] + g_z_0_zzzzz_xxxzzz[k];

                g_z_0_xzzzzz_xyyyy[k] = -g_z_0_zzzzz_xyyyy[k] * cd_x[k] + g_z_0_zzzzz_xxyyyy[k];

                g_z_0_xzzzzz_xyyyz[k] = -g_z_0_zzzzz_xyyyz[k] * cd_x[k] + g_z_0_zzzzz_xxyyyz[k];

                g_z_0_xzzzzz_xyyzz[k] = -g_z_0_zzzzz_xyyzz[k] * cd_x[k] + g_z_0_zzzzz_xxyyzz[k];

                g_z_0_xzzzzz_xyzzz[k] = -g_z_0_zzzzz_xyzzz[k] * cd_x[k] + g_z_0_zzzzz_xxyzzz[k];

                g_z_0_xzzzzz_xzzzz[k] = -g_z_0_zzzzz_xzzzz[k] * cd_x[k] + g_z_0_zzzzz_xxzzzz[k];

                g_z_0_xzzzzz_yyyyy[k] = -g_z_0_zzzzz_yyyyy[k] * cd_x[k] + g_z_0_zzzzz_xyyyyy[k];

                g_z_0_xzzzzz_yyyyz[k] = -g_z_0_zzzzz_yyyyz[k] * cd_x[k] + g_z_0_zzzzz_xyyyyz[k];

                g_z_0_xzzzzz_yyyzz[k] = -g_z_0_zzzzz_yyyzz[k] * cd_x[k] + g_z_0_zzzzz_xyyyzz[k];

                g_z_0_xzzzzz_yyzzz[k] = -g_z_0_zzzzz_yyzzz[k] * cd_x[k] + g_z_0_zzzzz_xyyzzz[k];

                g_z_0_xzzzzz_yzzzz[k] = -g_z_0_zzzzz_yzzzz[k] * cd_x[k] + g_z_0_zzzzz_xyzzzz[k];

                g_z_0_xzzzzz_zzzzz[k] = -g_z_0_zzzzz_zzzzz[k] * cd_x[k] + g_z_0_zzzzz_xzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 441);

            auto g_z_0_yyyyyy_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 442);

            auto g_z_0_yyyyyy_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 443);

            auto g_z_0_yyyyyy_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 444);

            auto g_z_0_yyyyyy_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 445);

            auto g_z_0_yyyyyy_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 446);

            auto g_z_0_yyyyyy_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 447);

            auto g_z_0_yyyyyy_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 448);

            auto g_z_0_yyyyyy_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 449);

            auto g_z_0_yyyyyy_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 450);

            auto g_z_0_yyyyyy_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 451);

            auto g_z_0_yyyyyy_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 452);

            auto g_z_0_yyyyyy_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 453);

            auto g_z_0_yyyyyy_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 454);

            auto g_z_0_yyyyyy_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 455);

            auto g_z_0_yyyyyy_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 456);

            auto g_z_0_yyyyyy_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 457);

            auto g_z_0_yyyyyy_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 458);

            auto g_z_0_yyyyyy_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 459);

            auto g_z_0_yyyyyy_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 460);

            auto g_z_0_yyyyyy_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 461);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyy_xxxxx, g_z_0_yyyyy_xxxxxy, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxxyy, g_z_0_yyyyy_xxxxyz, g_z_0_yyyyy_xxxxz, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyyy, g_z_0_yyyyy_xxxyyz, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxyzz, g_z_0_yyyyy_xxxzz, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyyy, g_z_0_yyyyy_xxyyyz, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyyzz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxyzzz, g_z_0_yyyyy_xxzzz, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyyy, g_z_0_yyyyy_xyyyyz, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyyzz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyyzzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xyzzzz, g_z_0_yyyyy_xzzzz, g_z_0_yyyyy_yyyyy, g_z_0_yyyyy_yyyyyy, g_z_0_yyyyy_yyyyyz, g_z_0_yyyyy_yyyyz, g_z_0_yyyyy_yyyyzz, g_z_0_yyyyy_yyyzz, g_z_0_yyyyy_yyyzzz, g_z_0_yyyyy_yyzzz, g_z_0_yyyyy_yyzzzz, g_z_0_yyyyy_yzzzz, g_z_0_yyyyy_yzzzzz, g_z_0_yyyyy_zzzzz, g_z_0_yyyyyy_xxxxx, g_z_0_yyyyyy_xxxxy, g_z_0_yyyyyy_xxxxz, g_z_0_yyyyyy_xxxyy, g_z_0_yyyyyy_xxxyz, g_z_0_yyyyyy_xxxzz, g_z_0_yyyyyy_xxyyy, g_z_0_yyyyyy_xxyyz, g_z_0_yyyyyy_xxyzz, g_z_0_yyyyyy_xxzzz, g_z_0_yyyyyy_xyyyy, g_z_0_yyyyyy_xyyyz, g_z_0_yyyyyy_xyyzz, g_z_0_yyyyyy_xyzzz, g_z_0_yyyyyy_xzzzz, g_z_0_yyyyyy_yyyyy, g_z_0_yyyyyy_yyyyz, g_z_0_yyyyyy_yyyzz, g_z_0_yyyyyy_yyzzz, g_z_0_yyyyyy_yzzzz, g_z_0_yyyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_xxxxx[k] = -g_z_0_yyyyy_xxxxx[k] * cd_y[k] + g_z_0_yyyyy_xxxxxy[k];

                g_z_0_yyyyyy_xxxxy[k] = -g_z_0_yyyyy_xxxxy[k] * cd_y[k] + g_z_0_yyyyy_xxxxyy[k];

                g_z_0_yyyyyy_xxxxz[k] = -g_z_0_yyyyy_xxxxz[k] * cd_y[k] + g_z_0_yyyyy_xxxxyz[k];

                g_z_0_yyyyyy_xxxyy[k] = -g_z_0_yyyyy_xxxyy[k] * cd_y[k] + g_z_0_yyyyy_xxxyyy[k];

                g_z_0_yyyyyy_xxxyz[k] = -g_z_0_yyyyy_xxxyz[k] * cd_y[k] + g_z_0_yyyyy_xxxyyz[k];

                g_z_0_yyyyyy_xxxzz[k] = -g_z_0_yyyyy_xxxzz[k] * cd_y[k] + g_z_0_yyyyy_xxxyzz[k];

                g_z_0_yyyyyy_xxyyy[k] = -g_z_0_yyyyy_xxyyy[k] * cd_y[k] + g_z_0_yyyyy_xxyyyy[k];

                g_z_0_yyyyyy_xxyyz[k] = -g_z_0_yyyyy_xxyyz[k] * cd_y[k] + g_z_0_yyyyy_xxyyyz[k];

                g_z_0_yyyyyy_xxyzz[k] = -g_z_0_yyyyy_xxyzz[k] * cd_y[k] + g_z_0_yyyyy_xxyyzz[k];

                g_z_0_yyyyyy_xxzzz[k] = -g_z_0_yyyyy_xxzzz[k] * cd_y[k] + g_z_0_yyyyy_xxyzzz[k];

                g_z_0_yyyyyy_xyyyy[k] = -g_z_0_yyyyy_xyyyy[k] * cd_y[k] + g_z_0_yyyyy_xyyyyy[k];

                g_z_0_yyyyyy_xyyyz[k] = -g_z_0_yyyyy_xyyyz[k] * cd_y[k] + g_z_0_yyyyy_xyyyyz[k];

                g_z_0_yyyyyy_xyyzz[k] = -g_z_0_yyyyy_xyyzz[k] * cd_y[k] + g_z_0_yyyyy_xyyyzz[k];

                g_z_0_yyyyyy_xyzzz[k] = -g_z_0_yyyyy_xyzzz[k] * cd_y[k] + g_z_0_yyyyy_xyyzzz[k];

                g_z_0_yyyyyy_xzzzz[k] = -g_z_0_yyyyy_xzzzz[k] * cd_y[k] + g_z_0_yyyyy_xyzzzz[k];

                g_z_0_yyyyyy_yyyyy[k] = -g_z_0_yyyyy_yyyyy[k] * cd_y[k] + g_z_0_yyyyy_yyyyyy[k];

                g_z_0_yyyyyy_yyyyz[k] = -g_z_0_yyyyy_yyyyz[k] * cd_y[k] + g_z_0_yyyyy_yyyyyz[k];

                g_z_0_yyyyyy_yyyzz[k] = -g_z_0_yyyyy_yyyzz[k] * cd_y[k] + g_z_0_yyyyy_yyyyzz[k];

                g_z_0_yyyyyy_yyzzz[k] = -g_z_0_yyyyy_yyzzz[k] * cd_y[k] + g_z_0_yyyyy_yyyzzz[k];

                g_z_0_yyyyyy_yzzzz[k] = -g_z_0_yyyyy_yzzzz[k] * cd_y[k] + g_z_0_yyyyy_yyzzzz[k];

                g_z_0_yyyyyy_zzzzz[k] = -g_z_0_yyyyy_zzzzz[k] * cd_y[k] + g_z_0_yyyyy_yzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 462);

            auto g_z_0_yyyyyz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 463);

            auto g_z_0_yyyyyz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 464);

            auto g_z_0_yyyyyz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 465);

            auto g_z_0_yyyyyz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 466);

            auto g_z_0_yyyyyz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 467);

            auto g_z_0_yyyyyz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 468);

            auto g_z_0_yyyyyz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 469);

            auto g_z_0_yyyyyz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 470);

            auto g_z_0_yyyyyz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 471);

            auto g_z_0_yyyyyz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 472);

            auto g_z_0_yyyyyz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 473);

            auto g_z_0_yyyyyz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 474);

            auto g_z_0_yyyyyz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 475);

            auto g_z_0_yyyyyz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 476);

            auto g_z_0_yyyyyz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 477);

            auto g_z_0_yyyyyz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 478);

            auto g_z_0_yyyyyz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 479);

            auto g_z_0_yyyyyz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 480);

            auto g_z_0_yyyyyz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 481);

            auto g_z_0_yyyyyz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 482);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyyz_xxxxx, g_z_0_yyyyyz_xxxxy, g_z_0_yyyyyz_xxxxz, g_z_0_yyyyyz_xxxyy, g_z_0_yyyyyz_xxxyz, g_z_0_yyyyyz_xxxzz, g_z_0_yyyyyz_xxyyy, g_z_0_yyyyyz_xxyyz, g_z_0_yyyyyz_xxyzz, g_z_0_yyyyyz_xxzzz, g_z_0_yyyyyz_xyyyy, g_z_0_yyyyyz_xyyyz, g_z_0_yyyyyz_xyyzz, g_z_0_yyyyyz_xyzzz, g_z_0_yyyyyz_xzzzz, g_z_0_yyyyyz_yyyyy, g_z_0_yyyyyz_yyyyz, g_z_0_yyyyyz_yyyzz, g_z_0_yyyyyz_yyzzz, g_z_0_yyyyyz_yzzzz, g_z_0_yyyyyz_zzzzz, g_z_0_yyyyz_xxxxx, g_z_0_yyyyz_xxxxxy, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxxyy, g_z_0_yyyyz_xxxxyz, g_z_0_yyyyz_xxxxz, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyyy, g_z_0_yyyyz_xxxyyz, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxyzz, g_z_0_yyyyz_xxxzz, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyyy, g_z_0_yyyyz_xxyyyz, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyyzz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxyzzz, g_z_0_yyyyz_xxzzz, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyyy, g_z_0_yyyyz_xyyyyz, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyyzz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyyzzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xyzzzz, g_z_0_yyyyz_xzzzz, g_z_0_yyyyz_yyyyy, g_z_0_yyyyz_yyyyyy, g_z_0_yyyyz_yyyyyz, g_z_0_yyyyz_yyyyz, g_z_0_yyyyz_yyyyzz, g_z_0_yyyyz_yyyzz, g_z_0_yyyyz_yyyzzz, g_z_0_yyyyz_yyzzz, g_z_0_yyyyz_yyzzzz, g_z_0_yyyyz_yzzzz, g_z_0_yyyyz_yzzzzz, g_z_0_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_xxxxx[k] = -g_z_0_yyyyz_xxxxx[k] * cd_y[k] + g_z_0_yyyyz_xxxxxy[k];

                g_z_0_yyyyyz_xxxxy[k] = -g_z_0_yyyyz_xxxxy[k] * cd_y[k] + g_z_0_yyyyz_xxxxyy[k];

                g_z_0_yyyyyz_xxxxz[k] = -g_z_0_yyyyz_xxxxz[k] * cd_y[k] + g_z_0_yyyyz_xxxxyz[k];

                g_z_0_yyyyyz_xxxyy[k] = -g_z_0_yyyyz_xxxyy[k] * cd_y[k] + g_z_0_yyyyz_xxxyyy[k];

                g_z_0_yyyyyz_xxxyz[k] = -g_z_0_yyyyz_xxxyz[k] * cd_y[k] + g_z_0_yyyyz_xxxyyz[k];

                g_z_0_yyyyyz_xxxzz[k] = -g_z_0_yyyyz_xxxzz[k] * cd_y[k] + g_z_0_yyyyz_xxxyzz[k];

                g_z_0_yyyyyz_xxyyy[k] = -g_z_0_yyyyz_xxyyy[k] * cd_y[k] + g_z_0_yyyyz_xxyyyy[k];

                g_z_0_yyyyyz_xxyyz[k] = -g_z_0_yyyyz_xxyyz[k] * cd_y[k] + g_z_0_yyyyz_xxyyyz[k];

                g_z_0_yyyyyz_xxyzz[k] = -g_z_0_yyyyz_xxyzz[k] * cd_y[k] + g_z_0_yyyyz_xxyyzz[k];

                g_z_0_yyyyyz_xxzzz[k] = -g_z_0_yyyyz_xxzzz[k] * cd_y[k] + g_z_0_yyyyz_xxyzzz[k];

                g_z_0_yyyyyz_xyyyy[k] = -g_z_0_yyyyz_xyyyy[k] * cd_y[k] + g_z_0_yyyyz_xyyyyy[k];

                g_z_0_yyyyyz_xyyyz[k] = -g_z_0_yyyyz_xyyyz[k] * cd_y[k] + g_z_0_yyyyz_xyyyyz[k];

                g_z_0_yyyyyz_xyyzz[k] = -g_z_0_yyyyz_xyyzz[k] * cd_y[k] + g_z_0_yyyyz_xyyyzz[k];

                g_z_0_yyyyyz_xyzzz[k] = -g_z_0_yyyyz_xyzzz[k] * cd_y[k] + g_z_0_yyyyz_xyyzzz[k];

                g_z_0_yyyyyz_xzzzz[k] = -g_z_0_yyyyz_xzzzz[k] * cd_y[k] + g_z_0_yyyyz_xyzzzz[k];

                g_z_0_yyyyyz_yyyyy[k] = -g_z_0_yyyyz_yyyyy[k] * cd_y[k] + g_z_0_yyyyz_yyyyyy[k];

                g_z_0_yyyyyz_yyyyz[k] = -g_z_0_yyyyz_yyyyz[k] * cd_y[k] + g_z_0_yyyyz_yyyyyz[k];

                g_z_0_yyyyyz_yyyzz[k] = -g_z_0_yyyyz_yyyzz[k] * cd_y[k] + g_z_0_yyyyz_yyyyzz[k];

                g_z_0_yyyyyz_yyzzz[k] = -g_z_0_yyyyz_yyzzz[k] * cd_y[k] + g_z_0_yyyyz_yyyzzz[k];

                g_z_0_yyyyyz_yzzzz[k] = -g_z_0_yyyyz_yzzzz[k] * cd_y[k] + g_z_0_yyyyz_yyzzzz[k];

                g_z_0_yyyyyz_zzzzz[k] = -g_z_0_yyyyz_zzzzz[k] * cd_y[k] + g_z_0_yyyyz_yzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 483);

            auto g_z_0_yyyyzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 484);

            auto g_z_0_yyyyzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 485);

            auto g_z_0_yyyyzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 486);

            auto g_z_0_yyyyzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 487);

            auto g_z_0_yyyyzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 488);

            auto g_z_0_yyyyzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 489);

            auto g_z_0_yyyyzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 490);

            auto g_z_0_yyyyzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 491);

            auto g_z_0_yyyyzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 492);

            auto g_z_0_yyyyzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 493);

            auto g_z_0_yyyyzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 494);

            auto g_z_0_yyyyzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 495);

            auto g_z_0_yyyyzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 496);

            auto g_z_0_yyyyzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 497);

            auto g_z_0_yyyyzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 498);

            auto g_z_0_yyyyzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 499);

            auto g_z_0_yyyyzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 500);

            auto g_z_0_yyyyzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 501);

            auto g_z_0_yyyyzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 502);

            auto g_z_0_yyyyzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 503);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyzz_xxxxx, g_z_0_yyyyzz_xxxxy, g_z_0_yyyyzz_xxxxz, g_z_0_yyyyzz_xxxyy, g_z_0_yyyyzz_xxxyz, g_z_0_yyyyzz_xxxzz, g_z_0_yyyyzz_xxyyy, g_z_0_yyyyzz_xxyyz, g_z_0_yyyyzz_xxyzz, g_z_0_yyyyzz_xxzzz, g_z_0_yyyyzz_xyyyy, g_z_0_yyyyzz_xyyyz, g_z_0_yyyyzz_xyyzz, g_z_0_yyyyzz_xyzzz, g_z_0_yyyyzz_xzzzz, g_z_0_yyyyzz_yyyyy, g_z_0_yyyyzz_yyyyz, g_z_0_yyyyzz_yyyzz, g_z_0_yyyyzz_yyzzz, g_z_0_yyyyzz_yzzzz, g_z_0_yyyyzz_zzzzz, g_z_0_yyyzz_xxxxx, g_z_0_yyyzz_xxxxxy, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxxyy, g_z_0_yyyzz_xxxxyz, g_z_0_yyyzz_xxxxz, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyyy, g_z_0_yyyzz_xxxyyz, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxyzz, g_z_0_yyyzz_xxxzz, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyyy, g_z_0_yyyzz_xxyyyz, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyyzz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxyzzz, g_z_0_yyyzz_xxzzz, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyyy, g_z_0_yyyzz_xyyyyz, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyyzz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyyzzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xyzzzz, g_z_0_yyyzz_xzzzz, g_z_0_yyyzz_yyyyy, g_z_0_yyyzz_yyyyyy, g_z_0_yyyzz_yyyyyz, g_z_0_yyyzz_yyyyz, g_z_0_yyyzz_yyyyzz, g_z_0_yyyzz_yyyzz, g_z_0_yyyzz_yyyzzz, g_z_0_yyyzz_yyzzz, g_z_0_yyyzz_yyzzzz, g_z_0_yyyzz_yzzzz, g_z_0_yyyzz_yzzzzz, g_z_0_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_xxxxx[k] = -g_z_0_yyyzz_xxxxx[k] * cd_y[k] + g_z_0_yyyzz_xxxxxy[k];

                g_z_0_yyyyzz_xxxxy[k] = -g_z_0_yyyzz_xxxxy[k] * cd_y[k] + g_z_0_yyyzz_xxxxyy[k];

                g_z_0_yyyyzz_xxxxz[k] = -g_z_0_yyyzz_xxxxz[k] * cd_y[k] + g_z_0_yyyzz_xxxxyz[k];

                g_z_0_yyyyzz_xxxyy[k] = -g_z_0_yyyzz_xxxyy[k] * cd_y[k] + g_z_0_yyyzz_xxxyyy[k];

                g_z_0_yyyyzz_xxxyz[k] = -g_z_0_yyyzz_xxxyz[k] * cd_y[k] + g_z_0_yyyzz_xxxyyz[k];

                g_z_0_yyyyzz_xxxzz[k] = -g_z_0_yyyzz_xxxzz[k] * cd_y[k] + g_z_0_yyyzz_xxxyzz[k];

                g_z_0_yyyyzz_xxyyy[k] = -g_z_0_yyyzz_xxyyy[k] * cd_y[k] + g_z_0_yyyzz_xxyyyy[k];

                g_z_0_yyyyzz_xxyyz[k] = -g_z_0_yyyzz_xxyyz[k] * cd_y[k] + g_z_0_yyyzz_xxyyyz[k];

                g_z_0_yyyyzz_xxyzz[k] = -g_z_0_yyyzz_xxyzz[k] * cd_y[k] + g_z_0_yyyzz_xxyyzz[k];

                g_z_0_yyyyzz_xxzzz[k] = -g_z_0_yyyzz_xxzzz[k] * cd_y[k] + g_z_0_yyyzz_xxyzzz[k];

                g_z_0_yyyyzz_xyyyy[k] = -g_z_0_yyyzz_xyyyy[k] * cd_y[k] + g_z_0_yyyzz_xyyyyy[k];

                g_z_0_yyyyzz_xyyyz[k] = -g_z_0_yyyzz_xyyyz[k] * cd_y[k] + g_z_0_yyyzz_xyyyyz[k];

                g_z_0_yyyyzz_xyyzz[k] = -g_z_0_yyyzz_xyyzz[k] * cd_y[k] + g_z_0_yyyzz_xyyyzz[k];

                g_z_0_yyyyzz_xyzzz[k] = -g_z_0_yyyzz_xyzzz[k] * cd_y[k] + g_z_0_yyyzz_xyyzzz[k];

                g_z_0_yyyyzz_xzzzz[k] = -g_z_0_yyyzz_xzzzz[k] * cd_y[k] + g_z_0_yyyzz_xyzzzz[k];

                g_z_0_yyyyzz_yyyyy[k] = -g_z_0_yyyzz_yyyyy[k] * cd_y[k] + g_z_0_yyyzz_yyyyyy[k];

                g_z_0_yyyyzz_yyyyz[k] = -g_z_0_yyyzz_yyyyz[k] * cd_y[k] + g_z_0_yyyzz_yyyyyz[k];

                g_z_0_yyyyzz_yyyzz[k] = -g_z_0_yyyzz_yyyzz[k] * cd_y[k] + g_z_0_yyyzz_yyyyzz[k];

                g_z_0_yyyyzz_yyzzz[k] = -g_z_0_yyyzz_yyzzz[k] * cd_y[k] + g_z_0_yyyzz_yyyzzz[k];

                g_z_0_yyyyzz_yzzzz[k] = -g_z_0_yyyzz_yzzzz[k] * cd_y[k] + g_z_0_yyyzz_yyzzzz[k];

                g_z_0_yyyyzz_zzzzz[k] = -g_z_0_yyyzz_zzzzz[k] * cd_y[k] + g_z_0_yyyzz_yzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 504);

            auto g_z_0_yyyzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 505);

            auto g_z_0_yyyzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 506);

            auto g_z_0_yyyzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 507);

            auto g_z_0_yyyzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 508);

            auto g_z_0_yyyzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 509);

            auto g_z_0_yyyzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 510);

            auto g_z_0_yyyzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 511);

            auto g_z_0_yyyzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 512);

            auto g_z_0_yyyzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 513);

            auto g_z_0_yyyzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 514);

            auto g_z_0_yyyzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 515);

            auto g_z_0_yyyzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 516);

            auto g_z_0_yyyzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 517);

            auto g_z_0_yyyzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 518);

            auto g_z_0_yyyzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 519);

            auto g_z_0_yyyzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 520);

            auto g_z_0_yyyzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 521);

            auto g_z_0_yyyzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 522);

            auto g_z_0_yyyzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 523);

            auto g_z_0_yyyzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 524);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzzz_xxxxx, g_z_0_yyyzzz_xxxxy, g_z_0_yyyzzz_xxxxz, g_z_0_yyyzzz_xxxyy, g_z_0_yyyzzz_xxxyz, g_z_0_yyyzzz_xxxzz, g_z_0_yyyzzz_xxyyy, g_z_0_yyyzzz_xxyyz, g_z_0_yyyzzz_xxyzz, g_z_0_yyyzzz_xxzzz, g_z_0_yyyzzz_xyyyy, g_z_0_yyyzzz_xyyyz, g_z_0_yyyzzz_xyyzz, g_z_0_yyyzzz_xyzzz, g_z_0_yyyzzz_xzzzz, g_z_0_yyyzzz_yyyyy, g_z_0_yyyzzz_yyyyz, g_z_0_yyyzzz_yyyzz, g_z_0_yyyzzz_yyzzz, g_z_0_yyyzzz_yzzzz, g_z_0_yyyzzz_zzzzz, g_z_0_yyzzz_xxxxx, g_z_0_yyzzz_xxxxxy, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxxyy, g_z_0_yyzzz_xxxxyz, g_z_0_yyzzz_xxxxz, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyyy, g_z_0_yyzzz_xxxyyz, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxyzz, g_z_0_yyzzz_xxxzz, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyyy, g_z_0_yyzzz_xxyyyz, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyyzz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxyzzz, g_z_0_yyzzz_xxzzz, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyyy, g_z_0_yyzzz_xyyyyz, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyyzz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyyzzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xyzzzz, g_z_0_yyzzz_xzzzz, g_z_0_yyzzz_yyyyy, g_z_0_yyzzz_yyyyyy, g_z_0_yyzzz_yyyyyz, g_z_0_yyzzz_yyyyz, g_z_0_yyzzz_yyyyzz, g_z_0_yyzzz_yyyzz, g_z_0_yyzzz_yyyzzz, g_z_0_yyzzz_yyzzz, g_z_0_yyzzz_yyzzzz, g_z_0_yyzzz_yzzzz, g_z_0_yyzzz_yzzzzz, g_z_0_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_xxxxx[k] = -g_z_0_yyzzz_xxxxx[k] * cd_y[k] + g_z_0_yyzzz_xxxxxy[k];

                g_z_0_yyyzzz_xxxxy[k] = -g_z_0_yyzzz_xxxxy[k] * cd_y[k] + g_z_0_yyzzz_xxxxyy[k];

                g_z_0_yyyzzz_xxxxz[k] = -g_z_0_yyzzz_xxxxz[k] * cd_y[k] + g_z_0_yyzzz_xxxxyz[k];

                g_z_0_yyyzzz_xxxyy[k] = -g_z_0_yyzzz_xxxyy[k] * cd_y[k] + g_z_0_yyzzz_xxxyyy[k];

                g_z_0_yyyzzz_xxxyz[k] = -g_z_0_yyzzz_xxxyz[k] * cd_y[k] + g_z_0_yyzzz_xxxyyz[k];

                g_z_0_yyyzzz_xxxzz[k] = -g_z_0_yyzzz_xxxzz[k] * cd_y[k] + g_z_0_yyzzz_xxxyzz[k];

                g_z_0_yyyzzz_xxyyy[k] = -g_z_0_yyzzz_xxyyy[k] * cd_y[k] + g_z_0_yyzzz_xxyyyy[k];

                g_z_0_yyyzzz_xxyyz[k] = -g_z_0_yyzzz_xxyyz[k] * cd_y[k] + g_z_0_yyzzz_xxyyyz[k];

                g_z_0_yyyzzz_xxyzz[k] = -g_z_0_yyzzz_xxyzz[k] * cd_y[k] + g_z_0_yyzzz_xxyyzz[k];

                g_z_0_yyyzzz_xxzzz[k] = -g_z_0_yyzzz_xxzzz[k] * cd_y[k] + g_z_0_yyzzz_xxyzzz[k];

                g_z_0_yyyzzz_xyyyy[k] = -g_z_0_yyzzz_xyyyy[k] * cd_y[k] + g_z_0_yyzzz_xyyyyy[k];

                g_z_0_yyyzzz_xyyyz[k] = -g_z_0_yyzzz_xyyyz[k] * cd_y[k] + g_z_0_yyzzz_xyyyyz[k];

                g_z_0_yyyzzz_xyyzz[k] = -g_z_0_yyzzz_xyyzz[k] * cd_y[k] + g_z_0_yyzzz_xyyyzz[k];

                g_z_0_yyyzzz_xyzzz[k] = -g_z_0_yyzzz_xyzzz[k] * cd_y[k] + g_z_0_yyzzz_xyyzzz[k];

                g_z_0_yyyzzz_xzzzz[k] = -g_z_0_yyzzz_xzzzz[k] * cd_y[k] + g_z_0_yyzzz_xyzzzz[k];

                g_z_0_yyyzzz_yyyyy[k] = -g_z_0_yyzzz_yyyyy[k] * cd_y[k] + g_z_0_yyzzz_yyyyyy[k];

                g_z_0_yyyzzz_yyyyz[k] = -g_z_0_yyzzz_yyyyz[k] * cd_y[k] + g_z_0_yyzzz_yyyyyz[k];

                g_z_0_yyyzzz_yyyzz[k] = -g_z_0_yyzzz_yyyzz[k] * cd_y[k] + g_z_0_yyzzz_yyyyzz[k];

                g_z_0_yyyzzz_yyzzz[k] = -g_z_0_yyzzz_yyzzz[k] * cd_y[k] + g_z_0_yyzzz_yyyzzz[k];

                g_z_0_yyyzzz_yzzzz[k] = -g_z_0_yyzzz_yzzzz[k] * cd_y[k] + g_z_0_yyzzz_yyzzzz[k];

                g_z_0_yyyzzz_zzzzz[k] = -g_z_0_yyzzz_zzzzz[k] * cd_y[k] + g_z_0_yyzzz_yzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 525);

            auto g_z_0_yyzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 526);

            auto g_z_0_yyzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 527);

            auto g_z_0_yyzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 528);

            auto g_z_0_yyzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 529);

            auto g_z_0_yyzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 530);

            auto g_z_0_yyzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 531);

            auto g_z_0_yyzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 532);

            auto g_z_0_yyzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 533);

            auto g_z_0_yyzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 534);

            auto g_z_0_yyzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 535);

            auto g_z_0_yyzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 536);

            auto g_z_0_yyzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 537);

            auto g_z_0_yyzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 538);

            auto g_z_0_yyzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 539);

            auto g_z_0_yyzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 540);

            auto g_z_0_yyzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 541);

            auto g_z_0_yyzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 542);

            auto g_z_0_yyzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 543);

            auto g_z_0_yyzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 544);

            auto g_z_0_yyzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 545);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzzz_xxxxx, g_z_0_yyzzzz_xxxxy, g_z_0_yyzzzz_xxxxz, g_z_0_yyzzzz_xxxyy, g_z_0_yyzzzz_xxxyz, g_z_0_yyzzzz_xxxzz, g_z_0_yyzzzz_xxyyy, g_z_0_yyzzzz_xxyyz, g_z_0_yyzzzz_xxyzz, g_z_0_yyzzzz_xxzzz, g_z_0_yyzzzz_xyyyy, g_z_0_yyzzzz_xyyyz, g_z_0_yyzzzz_xyyzz, g_z_0_yyzzzz_xyzzz, g_z_0_yyzzzz_xzzzz, g_z_0_yyzzzz_yyyyy, g_z_0_yyzzzz_yyyyz, g_z_0_yyzzzz_yyyzz, g_z_0_yyzzzz_yyzzz, g_z_0_yyzzzz_yzzzz, g_z_0_yyzzzz_zzzzz, g_z_0_yzzzz_xxxxx, g_z_0_yzzzz_xxxxxy, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxxyy, g_z_0_yzzzz_xxxxyz, g_z_0_yzzzz_xxxxz, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyyy, g_z_0_yzzzz_xxxyyz, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxyzz, g_z_0_yzzzz_xxxzz, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyyy, g_z_0_yzzzz_xxyyyz, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyyzz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxyzzz, g_z_0_yzzzz_xxzzz, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyyy, g_z_0_yzzzz_xyyyyz, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyyzz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyyzzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xyzzzz, g_z_0_yzzzz_xzzzz, g_z_0_yzzzz_yyyyy, g_z_0_yzzzz_yyyyyy, g_z_0_yzzzz_yyyyyz, g_z_0_yzzzz_yyyyz, g_z_0_yzzzz_yyyyzz, g_z_0_yzzzz_yyyzz, g_z_0_yzzzz_yyyzzz, g_z_0_yzzzz_yyzzz, g_z_0_yzzzz_yyzzzz, g_z_0_yzzzz_yzzzz, g_z_0_yzzzz_yzzzzz, g_z_0_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_xxxxx[k] = -g_z_0_yzzzz_xxxxx[k] * cd_y[k] + g_z_0_yzzzz_xxxxxy[k];

                g_z_0_yyzzzz_xxxxy[k] = -g_z_0_yzzzz_xxxxy[k] * cd_y[k] + g_z_0_yzzzz_xxxxyy[k];

                g_z_0_yyzzzz_xxxxz[k] = -g_z_0_yzzzz_xxxxz[k] * cd_y[k] + g_z_0_yzzzz_xxxxyz[k];

                g_z_0_yyzzzz_xxxyy[k] = -g_z_0_yzzzz_xxxyy[k] * cd_y[k] + g_z_0_yzzzz_xxxyyy[k];

                g_z_0_yyzzzz_xxxyz[k] = -g_z_0_yzzzz_xxxyz[k] * cd_y[k] + g_z_0_yzzzz_xxxyyz[k];

                g_z_0_yyzzzz_xxxzz[k] = -g_z_0_yzzzz_xxxzz[k] * cd_y[k] + g_z_0_yzzzz_xxxyzz[k];

                g_z_0_yyzzzz_xxyyy[k] = -g_z_0_yzzzz_xxyyy[k] * cd_y[k] + g_z_0_yzzzz_xxyyyy[k];

                g_z_0_yyzzzz_xxyyz[k] = -g_z_0_yzzzz_xxyyz[k] * cd_y[k] + g_z_0_yzzzz_xxyyyz[k];

                g_z_0_yyzzzz_xxyzz[k] = -g_z_0_yzzzz_xxyzz[k] * cd_y[k] + g_z_0_yzzzz_xxyyzz[k];

                g_z_0_yyzzzz_xxzzz[k] = -g_z_0_yzzzz_xxzzz[k] * cd_y[k] + g_z_0_yzzzz_xxyzzz[k];

                g_z_0_yyzzzz_xyyyy[k] = -g_z_0_yzzzz_xyyyy[k] * cd_y[k] + g_z_0_yzzzz_xyyyyy[k];

                g_z_0_yyzzzz_xyyyz[k] = -g_z_0_yzzzz_xyyyz[k] * cd_y[k] + g_z_0_yzzzz_xyyyyz[k];

                g_z_0_yyzzzz_xyyzz[k] = -g_z_0_yzzzz_xyyzz[k] * cd_y[k] + g_z_0_yzzzz_xyyyzz[k];

                g_z_0_yyzzzz_xyzzz[k] = -g_z_0_yzzzz_xyzzz[k] * cd_y[k] + g_z_0_yzzzz_xyyzzz[k];

                g_z_0_yyzzzz_xzzzz[k] = -g_z_0_yzzzz_xzzzz[k] * cd_y[k] + g_z_0_yzzzz_xyzzzz[k];

                g_z_0_yyzzzz_yyyyy[k] = -g_z_0_yzzzz_yyyyy[k] * cd_y[k] + g_z_0_yzzzz_yyyyyy[k];

                g_z_0_yyzzzz_yyyyz[k] = -g_z_0_yzzzz_yyyyz[k] * cd_y[k] + g_z_0_yzzzz_yyyyyz[k];

                g_z_0_yyzzzz_yyyzz[k] = -g_z_0_yzzzz_yyyzz[k] * cd_y[k] + g_z_0_yzzzz_yyyyzz[k];

                g_z_0_yyzzzz_yyzzz[k] = -g_z_0_yzzzz_yyzzz[k] * cd_y[k] + g_z_0_yzzzz_yyyzzz[k];

                g_z_0_yyzzzz_yzzzz[k] = -g_z_0_yzzzz_yzzzz[k] * cd_y[k] + g_z_0_yzzzz_yyzzzz[k];

                g_z_0_yyzzzz_zzzzz[k] = -g_z_0_yzzzz_zzzzz[k] * cd_y[k] + g_z_0_yzzzz_yzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 546);

            auto g_z_0_yzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 547);

            auto g_z_0_yzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 548);

            auto g_z_0_yzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 549);

            auto g_z_0_yzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 550);

            auto g_z_0_yzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 551);

            auto g_z_0_yzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 552);

            auto g_z_0_yzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 553);

            auto g_z_0_yzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 554);

            auto g_z_0_yzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 555);

            auto g_z_0_yzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 556);

            auto g_z_0_yzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 557);

            auto g_z_0_yzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 558);

            auto g_z_0_yzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 559);

            auto g_z_0_yzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 560);

            auto g_z_0_yzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 561);

            auto g_z_0_yzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 562);

            auto g_z_0_yzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 563);

            auto g_z_0_yzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 564);

            auto g_z_0_yzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 565);

            auto g_z_0_yzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 566);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzzz_xxxxx, g_z_0_yzzzzz_xxxxy, g_z_0_yzzzzz_xxxxz, g_z_0_yzzzzz_xxxyy, g_z_0_yzzzzz_xxxyz, g_z_0_yzzzzz_xxxzz, g_z_0_yzzzzz_xxyyy, g_z_0_yzzzzz_xxyyz, g_z_0_yzzzzz_xxyzz, g_z_0_yzzzzz_xxzzz, g_z_0_yzzzzz_xyyyy, g_z_0_yzzzzz_xyyyz, g_z_0_yzzzzz_xyyzz, g_z_0_yzzzzz_xyzzz, g_z_0_yzzzzz_xzzzz, g_z_0_yzzzzz_yyyyy, g_z_0_yzzzzz_yyyyz, g_z_0_yzzzzz_yyyzz, g_z_0_yzzzzz_yyzzz, g_z_0_yzzzzz_yzzzz, g_z_0_yzzzzz_zzzzz, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyyy, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_xxxxx[k] = -g_z_0_zzzzz_xxxxx[k] * cd_y[k] + g_z_0_zzzzz_xxxxxy[k];

                g_z_0_yzzzzz_xxxxy[k] = -g_z_0_zzzzz_xxxxy[k] * cd_y[k] + g_z_0_zzzzz_xxxxyy[k];

                g_z_0_yzzzzz_xxxxz[k] = -g_z_0_zzzzz_xxxxz[k] * cd_y[k] + g_z_0_zzzzz_xxxxyz[k];

                g_z_0_yzzzzz_xxxyy[k] = -g_z_0_zzzzz_xxxyy[k] * cd_y[k] + g_z_0_zzzzz_xxxyyy[k];

                g_z_0_yzzzzz_xxxyz[k] = -g_z_0_zzzzz_xxxyz[k] * cd_y[k] + g_z_0_zzzzz_xxxyyz[k];

                g_z_0_yzzzzz_xxxzz[k] = -g_z_0_zzzzz_xxxzz[k] * cd_y[k] + g_z_0_zzzzz_xxxyzz[k];

                g_z_0_yzzzzz_xxyyy[k] = -g_z_0_zzzzz_xxyyy[k] * cd_y[k] + g_z_0_zzzzz_xxyyyy[k];

                g_z_0_yzzzzz_xxyyz[k] = -g_z_0_zzzzz_xxyyz[k] * cd_y[k] + g_z_0_zzzzz_xxyyyz[k];

                g_z_0_yzzzzz_xxyzz[k] = -g_z_0_zzzzz_xxyzz[k] * cd_y[k] + g_z_0_zzzzz_xxyyzz[k];

                g_z_0_yzzzzz_xxzzz[k] = -g_z_0_zzzzz_xxzzz[k] * cd_y[k] + g_z_0_zzzzz_xxyzzz[k];

                g_z_0_yzzzzz_xyyyy[k] = -g_z_0_zzzzz_xyyyy[k] * cd_y[k] + g_z_0_zzzzz_xyyyyy[k];

                g_z_0_yzzzzz_xyyyz[k] = -g_z_0_zzzzz_xyyyz[k] * cd_y[k] + g_z_0_zzzzz_xyyyyz[k];

                g_z_0_yzzzzz_xyyzz[k] = -g_z_0_zzzzz_xyyzz[k] * cd_y[k] + g_z_0_zzzzz_xyyyzz[k];

                g_z_0_yzzzzz_xyzzz[k] = -g_z_0_zzzzz_xyzzz[k] * cd_y[k] + g_z_0_zzzzz_xyyzzz[k];

                g_z_0_yzzzzz_xzzzz[k] = -g_z_0_zzzzz_xzzzz[k] * cd_y[k] + g_z_0_zzzzz_xyzzzz[k];

                g_z_0_yzzzzz_yyyyy[k] = -g_z_0_zzzzz_yyyyy[k] * cd_y[k] + g_z_0_zzzzz_yyyyyy[k];

                g_z_0_yzzzzz_yyyyz[k] = -g_z_0_zzzzz_yyyyz[k] * cd_y[k] + g_z_0_zzzzz_yyyyyz[k];

                g_z_0_yzzzzz_yyyzz[k] = -g_z_0_zzzzz_yyyzz[k] * cd_y[k] + g_z_0_zzzzz_yyyyzz[k];

                g_z_0_yzzzzz_yyzzz[k] = -g_z_0_zzzzz_yyzzz[k] * cd_y[k] + g_z_0_zzzzz_yyyzzz[k];

                g_z_0_yzzzzz_yzzzz[k] = -g_z_0_zzzzz_yzzzz[k] * cd_y[k] + g_z_0_zzzzz_yyzzzz[k];

                g_z_0_yzzzzz_zzzzz[k] = -g_z_0_zzzzz_zzzzz[k] * cd_y[k] + g_z_0_zzzzz_yzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_xxxxx = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 567);

            auto g_z_0_zzzzzz_xxxxy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 568);

            auto g_z_0_zzzzzz_xxxxz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 569);

            auto g_z_0_zzzzzz_xxxyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 570);

            auto g_z_0_zzzzzz_xxxyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 571);

            auto g_z_0_zzzzzz_xxxzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 572);

            auto g_z_0_zzzzzz_xxyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 573);

            auto g_z_0_zzzzzz_xxyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 574);

            auto g_z_0_zzzzzz_xxyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 575);

            auto g_z_0_zzzzzz_xxzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 576);

            auto g_z_0_zzzzzz_xyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 577);

            auto g_z_0_zzzzzz_xyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 578);

            auto g_z_0_zzzzzz_xyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 579);

            auto g_z_0_zzzzzz_xyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 580);

            auto g_z_0_zzzzzz_xzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 581);

            auto g_z_0_zzzzzz_yyyyy = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 582);

            auto g_z_0_zzzzzz_yyyyz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 583);

            auto g_z_0_zzzzzz_yyyzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 584);

            auto g_z_0_zzzzzz_yyzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 585);

            auto g_z_0_zzzzzz_yzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 586);

            auto g_z_0_zzzzzz_zzzzz = cbuffer.data(ih_geom_10_off + 1176 * acomps * bcomps + 587);

            #pragma omp simd aligned(cd_z, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_zzzzz, g_z_0_zzzzz_zzzzzz, g_z_0_zzzzzz_xxxxx, g_z_0_zzzzzz_xxxxy, g_z_0_zzzzzz_xxxxz, g_z_0_zzzzzz_xxxyy, g_z_0_zzzzzz_xxxyz, g_z_0_zzzzzz_xxxzz, g_z_0_zzzzzz_xxyyy, g_z_0_zzzzzz_xxyyz, g_z_0_zzzzzz_xxyzz, g_z_0_zzzzzz_xxzzz, g_z_0_zzzzzz_xyyyy, g_z_0_zzzzzz_xyyyz, g_z_0_zzzzzz_xyyzz, g_z_0_zzzzzz_xyzzz, g_z_0_zzzzzz_xzzzz, g_z_0_zzzzzz_yyyyy, g_z_0_zzzzzz_yyyyz, g_z_0_zzzzzz_yyyzz, g_z_0_zzzzzz_yyzzz, g_z_0_zzzzzz_yzzzz, g_z_0_zzzzzz_zzzzz, g_zzzzz_xxxxx, g_zzzzz_xxxxy, g_zzzzz_xxxxz, g_zzzzz_xxxyy, g_zzzzz_xxxyz, g_zzzzz_xxxzz, g_zzzzz_xxyyy, g_zzzzz_xxyyz, g_zzzzz_xxyzz, g_zzzzz_xxzzz, g_zzzzz_xyyyy, g_zzzzz_xyyyz, g_zzzzz_xyyzz, g_zzzzz_xyzzz, g_zzzzz_xzzzz, g_zzzzz_yyyyy, g_zzzzz_yyyyz, g_zzzzz_yyyzz, g_zzzzz_yyzzz, g_zzzzz_yzzzz, g_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_xxxxx[k] = -g_zzzzz_xxxxx[k] - g_z_0_zzzzz_xxxxx[k] * cd_z[k] + g_z_0_zzzzz_xxxxxz[k];

                g_z_0_zzzzzz_xxxxy[k] = -g_zzzzz_xxxxy[k] - g_z_0_zzzzz_xxxxy[k] * cd_z[k] + g_z_0_zzzzz_xxxxyz[k];

                g_z_0_zzzzzz_xxxxz[k] = -g_zzzzz_xxxxz[k] - g_z_0_zzzzz_xxxxz[k] * cd_z[k] + g_z_0_zzzzz_xxxxzz[k];

                g_z_0_zzzzzz_xxxyy[k] = -g_zzzzz_xxxyy[k] - g_z_0_zzzzz_xxxyy[k] * cd_z[k] + g_z_0_zzzzz_xxxyyz[k];

                g_z_0_zzzzzz_xxxyz[k] = -g_zzzzz_xxxyz[k] - g_z_0_zzzzz_xxxyz[k] * cd_z[k] + g_z_0_zzzzz_xxxyzz[k];

                g_z_0_zzzzzz_xxxzz[k] = -g_zzzzz_xxxzz[k] - g_z_0_zzzzz_xxxzz[k] * cd_z[k] + g_z_0_zzzzz_xxxzzz[k];

                g_z_0_zzzzzz_xxyyy[k] = -g_zzzzz_xxyyy[k] - g_z_0_zzzzz_xxyyy[k] * cd_z[k] + g_z_0_zzzzz_xxyyyz[k];

                g_z_0_zzzzzz_xxyyz[k] = -g_zzzzz_xxyyz[k] - g_z_0_zzzzz_xxyyz[k] * cd_z[k] + g_z_0_zzzzz_xxyyzz[k];

                g_z_0_zzzzzz_xxyzz[k] = -g_zzzzz_xxyzz[k] - g_z_0_zzzzz_xxyzz[k] * cd_z[k] + g_z_0_zzzzz_xxyzzz[k];

                g_z_0_zzzzzz_xxzzz[k] = -g_zzzzz_xxzzz[k] - g_z_0_zzzzz_xxzzz[k] * cd_z[k] + g_z_0_zzzzz_xxzzzz[k];

                g_z_0_zzzzzz_xyyyy[k] = -g_zzzzz_xyyyy[k] - g_z_0_zzzzz_xyyyy[k] * cd_z[k] + g_z_0_zzzzz_xyyyyz[k];

                g_z_0_zzzzzz_xyyyz[k] = -g_zzzzz_xyyyz[k] - g_z_0_zzzzz_xyyyz[k] * cd_z[k] + g_z_0_zzzzz_xyyyzz[k];

                g_z_0_zzzzzz_xyyzz[k] = -g_zzzzz_xyyzz[k] - g_z_0_zzzzz_xyyzz[k] * cd_z[k] + g_z_0_zzzzz_xyyzzz[k];

                g_z_0_zzzzzz_xyzzz[k] = -g_zzzzz_xyzzz[k] - g_z_0_zzzzz_xyzzz[k] * cd_z[k] + g_z_0_zzzzz_xyzzzz[k];

                g_z_0_zzzzzz_xzzzz[k] = -g_zzzzz_xzzzz[k] - g_z_0_zzzzz_xzzzz[k] * cd_z[k] + g_z_0_zzzzz_xzzzzz[k];

                g_z_0_zzzzzz_yyyyy[k] = -g_zzzzz_yyyyy[k] - g_z_0_zzzzz_yyyyy[k] * cd_z[k] + g_z_0_zzzzz_yyyyyz[k];

                g_z_0_zzzzzz_yyyyz[k] = -g_zzzzz_yyyyz[k] - g_z_0_zzzzz_yyyyz[k] * cd_z[k] + g_z_0_zzzzz_yyyyzz[k];

                g_z_0_zzzzzz_yyyzz[k] = -g_zzzzz_yyyzz[k] - g_z_0_zzzzz_yyyzz[k] * cd_z[k] + g_z_0_zzzzz_yyyzzz[k];

                g_z_0_zzzzzz_yyzzz[k] = -g_zzzzz_yyzzz[k] - g_z_0_zzzzz_yyzzz[k] * cd_z[k] + g_z_0_zzzzz_yyzzzz[k];

                g_z_0_zzzzzz_yzzzz[k] = -g_zzzzz_yzzzz[k] - g_z_0_zzzzz_yzzzz[k] * cd_z[k] + g_z_0_zzzzz_yzzzzz[k];

                g_z_0_zzzzzz_zzzzz[k] = -g_zzzzz_zzzzz[k] - g_z_0_zzzzz_zzzzz[k] * cd_z[k] + g_z_0_zzzzz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

