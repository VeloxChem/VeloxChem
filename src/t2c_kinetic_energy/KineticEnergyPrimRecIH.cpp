#include "KineticEnergyPrimRecIH.hpp"

namespace kinrec { // kinrec namespace

auto
comp_prim_kinetic_energy_ih(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_ih,
                            const size_t idx_ovl_gh,
                            const size_t idx_kin_gh,
                            const size_t idx_kin_hg,
                            const size_t idx_kin_hh,
                            const size_t idx_ovl_ih,
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

    auto ts_xxxy_xxxxx = pbuffer.data(idx_ovl_gh + 21);

    auto ts_xxxy_xxxxz = pbuffer.data(idx_ovl_gh + 23);

    auto ts_xxxy_xxxzz = pbuffer.data(idx_ovl_gh + 26);

    auto ts_xxxy_xxzzz = pbuffer.data(idx_ovl_gh + 30);

    auto ts_xxxy_xzzzz = pbuffer.data(idx_ovl_gh + 35);

    auto ts_xxxz_xxxxx = pbuffer.data(idx_ovl_gh + 42);

    auto ts_xxxz_xxxxy = pbuffer.data(idx_ovl_gh + 43);

    auto ts_xxxz_xxxyy = pbuffer.data(idx_ovl_gh + 45);

    auto ts_xxxz_xxyyy = pbuffer.data(idx_ovl_gh + 48);

    auto ts_xxxz_xyyyy = pbuffer.data(idx_ovl_gh + 52);

    auto ts_xxyy_xxxxx = pbuffer.data(idx_ovl_gh + 63);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_ovl_gh + 64);

    auto ts_xxyy_xxxxz = pbuffer.data(idx_ovl_gh + 65);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_ovl_gh + 66);

    auto ts_xxyy_xxxyz = pbuffer.data(idx_ovl_gh + 67);

    auto ts_xxyy_xxxzz = pbuffer.data(idx_ovl_gh + 68);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_ovl_gh + 69);

    auto ts_xxyy_xxyyz = pbuffer.data(idx_ovl_gh + 70);

    auto ts_xxyy_xxyzz = pbuffer.data(idx_ovl_gh + 71);

    auto ts_xxyy_xxzzz = pbuffer.data(idx_ovl_gh + 72);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_ovl_gh + 73);

    auto ts_xxyy_xyyyz = pbuffer.data(idx_ovl_gh + 74);

    auto ts_xxyy_xyyzz = pbuffer.data(idx_ovl_gh + 75);

    auto ts_xxyy_xyzzz = pbuffer.data(idx_ovl_gh + 76);

    auto ts_xxyy_xzzzz = pbuffer.data(idx_ovl_gh + 77);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_ovl_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_ovl_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_ovl_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_ovl_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_ovl_gh + 82);

    auto ts_xxyy_zzzzz = pbuffer.data(idx_ovl_gh + 83);

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

    auto ts_xxzz_yyyyy = pbuffer.data(idx_ovl_gh + 120);

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

    auto ts_xyyy_zzzzz = pbuffer.data(idx_ovl_gh + 146);

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

    auto ts_xzzz_yyyyy = pbuffer.data(idx_ovl_gh + 204);

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

    auto ts_yyyz_xxxxy = pbuffer.data(idx_ovl_gh + 232);

    auto ts_yyyz_xxxyy = pbuffer.data(idx_ovl_gh + 234);

    auto ts_yyyz_xxyyy = pbuffer.data(idx_ovl_gh + 237);

    auto ts_yyyz_xyyyy = pbuffer.data(idx_ovl_gh + 241);

    auto ts_yyyz_yyyyy = pbuffer.data(idx_ovl_gh + 246);

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

    auto ts_yzzz_xxxxx = pbuffer.data(idx_ovl_gh + 273);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_ovl_gh + 275);

    auto ts_yzzz_xxxyz = pbuffer.data(idx_ovl_gh + 277);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_ovl_gh + 278);

    auto ts_yzzz_xxyyz = pbuffer.data(idx_ovl_gh + 280);

    auto ts_yzzz_xxyzz = pbuffer.data(idx_ovl_gh + 281);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_ovl_gh + 282);

    auto ts_yzzz_xyyyz = pbuffer.data(idx_ovl_gh + 284);

    auto ts_yzzz_xyyzz = pbuffer.data(idx_ovl_gh + 285);

    auto ts_yzzz_xyzzz = pbuffer.data(idx_ovl_gh + 286);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_ovl_gh + 287);

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

    // Set up components of auxiliary buffer : GH

    auto tk_xxxx_xxxxx = pbuffer.data(idx_kin_gh);

    auto tk_xxxx_xxxxy = pbuffer.data(idx_kin_gh + 1);

    auto tk_xxxx_xxxxz = pbuffer.data(idx_kin_gh + 2);

    auto tk_xxxx_xxxyy = pbuffer.data(idx_kin_gh + 3);

    auto tk_xxxx_xxxyz = pbuffer.data(idx_kin_gh + 4);

    auto tk_xxxx_xxxzz = pbuffer.data(idx_kin_gh + 5);

    auto tk_xxxx_xxyyy = pbuffer.data(idx_kin_gh + 6);

    auto tk_xxxx_xxyyz = pbuffer.data(idx_kin_gh + 7);

    auto tk_xxxx_xxyzz = pbuffer.data(idx_kin_gh + 8);

    auto tk_xxxx_xxzzz = pbuffer.data(idx_kin_gh + 9);

    auto tk_xxxx_xyyyy = pbuffer.data(idx_kin_gh + 10);

    auto tk_xxxx_xyyyz = pbuffer.data(idx_kin_gh + 11);

    auto tk_xxxx_xyyzz = pbuffer.data(idx_kin_gh + 12);

    auto tk_xxxx_xyzzz = pbuffer.data(idx_kin_gh + 13);

    auto tk_xxxx_xzzzz = pbuffer.data(idx_kin_gh + 14);

    auto tk_xxxx_yyyyy = pbuffer.data(idx_kin_gh + 15);

    auto tk_xxxx_yyyyz = pbuffer.data(idx_kin_gh + 16);

    auto tk_xxxx_yyyzz = pbuffer.data(idx_kin_gh + 17);

    auto tk_xxxx_yyzzz = pbuffer.data(idx_kin_gh + 18);

    auto tk_xxxx_yzzzz = pbuffer.data(idx_kin_gh + 19);

    auto tk_xxxx_zzzzz = pbuffer.data(idx_kin_gh + 20);

    auto tk_xxxy_xxxxx = pbuffer.data(idx_kin_gh + 21);

    auto tk_xxxy_xxxxz = pbuffer.data(idx_kin_gh + 23);

    auto tk_xxxy_xxxzz = pbuffer.data(idx_kin_gh + 26);

    auto tk_xxxy_xxzzz = pbuffer.data(idx_kin_gh + 30);

    auto tk_xxxy_xzzzz = pbuffer.data(idx_kin_gh + 35);

    auto tk_xxxz_xxxxx = pbuffer.data(idx_kin_gh + 42);

    auto tk_xxxz_xxxxy = pbuffer.data(idx_kin_gh + 43);

    auto tk_xxxz_xxxyy = pbuffer.data(idx_kin_gh + 45);

    auto tk_xxxz_xxyyy = pbuffer.data(idx_kin_gh + 48);

    auto tk_xxxz_xyyyy = pbuffer.data(idx_kin_gh + 52);

    auto tk_xxyy_xxxxx = pbuffer.data(idx_kin_gh + 63);

    auto tk_xxyy_xxxxy = pbuffer.data(idx_kin_gh + 64);

    auto tk_xxyy_xxxxz = pbuffer.data(idx_kin_gh + 65);

    auto tk_xxyy_xxxyy = pbuffer.data(idx_kin_gh + 66);

    auto tk_xxyy_xxxyz = pbuffer.data(idx_kin_gh + 67);

    auto tk_xxyy_xxxzz = pbuffer.data(idx_kin_gh + 68);

    auto tk_xxyy_xxyyy = pbuffer.data(idx_kin_gh + 69);

    auto tk_xxyy_xxyyz = pbuffer.data(idx_kin_gh + 70);

    auto tk_xxyy_xxyzz = pbuffer.data(idx_kin_gh + 71);

    auto tk_xxyy_xxzzz = pbuffer.data(idx_kin_gh + 72);

    auto tk_xxyy_xyyyy = pbuffer.data(idx_kin_gh + 73);

    auto tk_xxyy_xyyyz = pbuffer.data(idx_kin_gh + 74);

    auto tk_xxyy_xyyzz = pbuffer.data(idx_kin_gh + 75);

    auto tk_xxyy_xyzzz = pbuffer.data(idx_kin_gh + 76);

    auto tk_xxyy_xzzzz = pbuffer.data(idx_kin_gh + 77);

    auto tk_xxyy_yyyyy = pbuffer.data(idx_kin_gh + 78);

    auto tk_xxyy_yyyyz = pbuffer.data(idx_kin_gh + 79);

    auto tk_xxyy_yyyzz = pbuffer.data(idx_kin_gh + 80);

    auto tk_xxyy_yyzzz = pbuffer.data(idx_kin_gh + 81);

    auto tk_xxyy_yzzzz = pbuffer.data(idx_kin_gh + 82);

    auto tk_xxyy_zzzzz = pbuffer.data(idx_kin_gh + 83);

    auto tk_xxzz_xxxxx = pbuffer.data(idx_kin_gh + 105);

    auto tk_xxzz_xxxxy = pbuffer.data(idx_kin_gh + 106);

    auto tk_xxzz_xxxxz = pbuffer.data(idx_kin_gh + 107);

    auto tk_xxzz_xxxyy = pbuffer.data(idx_kin_gh + 108);

    auto tk_xxzz_xxxyz = pbuffer.data(idx_kin_gh + 109);

    auto tk_xxzz_xxxzz = pbuffer.data(idx_kin_gh + 110);

    auto tk_xxzz_xxyyy = pbuffer.data(idx_kin_gh + 111);

    auto tk_xxzz_xxyyz = pbuffer.data(idx_kin_gh + 112);

    auto tk_xxzz_xxyzz = pbuffer.data(idx_kin_gh + 113);

    auto tk_xxzz_xxzzz = pbuffer.data(idx_kin_gh + 114);

    auto tk_xxzz_xyyyy = pbuffer.data(idx_kin_gh + 115);

    auto tk_xxzz_xyyyz = pbuffer.data(idx_kin_gh + 116);

    auto tk_xxzz_xyyzz = pbuffer.data(idx_kin_gh + 117);

    auto tk_xxzz_xyzzz = pbuffer.data(idx_kin_gh + 118);

    auto tk_xxzz_xzzzz = pbuffer.data(idx_kin_gh + 119);

    auto tk_xxzz_yyyyy = pbuffer.data(idx_kin_gh + 120);

    auto tk_xxzz_yyyyz = pbuffer.data(idx_kin_gh + 121);

    auto tk_xxzz_yyyzz = pbuffer.data(idx_kin_gh + 122);

    auto tk_xxzz_yyzzz = pbuffer.data(idx_kin_gh + 123);

    auto tk_xxzz_yzzzz = pbuffer.data(idx_kin_gh + 124);

    auto tk_xxzz_zzzzz = pbuffer.data(idx_kin_gh + 125);

    auto tk_xyyy_xxxxy = pbuffer.data(idx_kin_gh + 127);

    auto tk_xyyy_xxxyy = pbuffer.data(idx_kin_gh + 129);

    auto tk_xyyy_xxxyz = pbuffer.data(idx_kin_gh + 130);

    auto tk_xyyy_xxyyy = pbuffer.data(idx_kin_gh + 132);

    auto tk_xyyy_xxyyz = pbuffer.data(idx_kin_gh + 133);

    auto tk_xyyy_xxyzz = pbuffer.data(idx_kin_gh + 134);

    auto tk_xyyy_xyyyy = pbuffer.data(idx_kin_gh + 136);

    auto tk_xyyy_xyyyz = pbuffer.data(idx_kin_gh + 137);

    auto tk_xyyy_xyyzz = pbuffer.data(idx_kin_gh + 138);

    auto tk_xyyy_xyzzz = pbuffer.data(idx_kin_gh + 139);

    auto tk_xyyy_yyyyy = pbuffer.data(idx_kin_gh + 141);

    auto tk_xyyy_yyyyz = pbuffer.data(idx_kin_gh + 142);

    auto tk_xyyy_yyyzz = pbuffer.data(idx_kin_gh + 143);

    auto tk_xyyy_yyzzz = pbuffer.data(idx_kin_gh + 144);

    auto tk_xyyy_yzzzz = pbuffer.data(idx_kin_gh + 145);

    auto tk_xyyy_zzzzz = pbuffer.data(idx_kin_gh + 146);

    auto tk_xzzz_xxxxz = pbuffer.data(idx_kin_gh + 191);

    auto tk_xzzz_xxxyz = pbuffer.data(idx_kin_gh + 193);

    auto tk_xzzz_xxxzz = pbuffer.data(idx_kin_gh + 194);

    auto tk_xzzz_xxyyz = pbuffer.data(idx_kin_gh + 196);

    auto tk_xzzz_xxyzz = pbuffer.data(idx_kin_gh + 197);

    auto tk_xzzz_xxzzz = pbuffer.data(idx_kin_gh + 198);

    auto tk_xzzz_xyyyz = pbuffer.data(idx_kin_gh + 200);

    auto tk_xzzz_xyyzz = pbuffer.data(idx_kin_gh + 201);

    auto tk_xzzz_xyzzz = pbuffer.data(idx_kin_gh + 202);

    auto tk_xzzz_xzzzz = pbuffer.data(idx_kin_gh + 203);

    auto tk_xzzz_yyyyy = pbuffer.data(idx_kin_gh + 204);

    auto tk_xzzz_yyyyz = pbuffer.data(idx_kin_gh + 205);

    auto tk_xzzz_yyyzz = pbuffer.data(idx_kin_gh + 206);

    auto tk_xzzz_yyzzz = pbuffer.data(idx_kin_gh + 207);

    auto tk_xzzz_yzzzz = pbuffer.data(idx_kin_gh + 208);

    auto tk_xzzz_zzzzz = pbuffer.data(idx_kin_gh + 209);

    auto tk_yyyy_xxxxx = pbuffer.data(idx_kin_gh + 210);

    auto tk_yyyy_xxxxy = pbuffer.data(idx_kin_gh + 211);

    auto tk_yyyy_xxxxz = pbuffer.data(idx_kin_gh + 212);

    auto tk_yyyy_xxxyy = pbuffer.data(idx_kin_gh + 213);

    auto tk_yyyy_xxxyz = pbuffer.data(idx_kin_gh + 214);

    auto tk_yyyy_xxxzz = pbuffer.data(idx_kin_gh + 215);

    auto tk_yyyy_xxyyy = pbuffer.data(idx_kin_gh + 216);

    auto tk_yyyy_xxyyz = pbuffer.data(idx_kin_gh + 217);

    auto tk_yyyy_xxyzz = pbuffer.data(idx_kin_gh + 218);

    auto tk_yyyy_xxzzz = pbuffer.data(idx_kin_gh + 219);

    auto tk_yyyy_xyyyy = pbuffer.data(idx_kin_gh + 220);

    auto tk_yyyy_xyyyz = pbuffer.data(idx_kin_gh + 221);

    auto tk_yyyy_xyyzz = pbuffer.data(idx_kin_gh + 222);

    auto tk_yyyy_xyzzz = pbuffer.data(idx_kin_gh + 223);

    auto tk_yyyy_xzzzz = pbuffer.data(idx_kin_gh + 224);

    auto tk_yyyy_yyyyy = pbuffer.data(idx_kin_gh + 225);

    auto tk_yyyy_yyyyz = pbuffer.data(idx_kin_gh + 226);

    auto tk_yyyy_yyyzz = pbuffer.data(idx_kin_gh + 227);

    auto tk_yyyy_yyzzz = pbuffer.data(idx_kin_gh + 228);

    auto tk_yyyy_yzzzz = pbuffer.data(idx_kin_gh + 229);

    auto tk_yyyy_zzzzz = pbuffer.data(idx_kin_gh + 230);

    auto tk_yyyz_xxxxy = pbuffer.data(idx_kin_gh + 232);

    auto tk_yyyz_xxxyy = pbuffer.data(idx_kin_gh + 234);

    auto tk_yyyz_xxyyy = pbuffer.data(idx_kin_gh + 237);

    auto tk_yyyz_xyyyy = pbuffer.data(idx_kin_gh + 241);

    auto tk_yyyz_yyyyy = pbuffer.data(idx_kin_gh + 246);

    auto tk_yyzz_xxxxx = pbuffer.data(idx_kin_gh + 252);

    auto tk_yyzz_xxxxy = pbuffer.data(idx_kin_gh + 253);

    auto tk_yyzz_xxxxz = pbuffer.data(idx_kin_gh + 254);

    auto tk_yyzz_xxxyy = pbuffer.data(idx_kin_gh + 255);

    auto tk_yyzz_xxxyz = pbuffer.data(idx_kin_gh + 256);

    auto tk_yyzz_xxxzz = pbuffer.data(idx_kin_gh + 257);

    auto tk_yyzz_xxyyy = pbuffer.data(idx_kin_gh + 258);

    auto tk_yyzz_xxyyz = pbuffer.data(idx_kin_gh + 259);

    auto tk_yyzz_xxyzz = pbuffer.data(idx_kin_gh + 260);

    auto tk_yyzz_xxzzz = pbuffer.data(idx_kin_gh + 261);

    auto tk_yyzz_xyyyy = pbuffer.data(idx_kin_gh + 262);

    auto tk_yyzz_xyyyz = pbuffer.data(idx_kin_gh + 263);

    auto tk_yyzz_xyyzz = pbuffer.data(idx_kin_gh + 264);

    auto tk_yyzz_xyzzz = pbuffer.data(idx_kin_gh + 265);

    auto tk_yyzz_xzzzz = pbuffer.data(idx_kin_gh + 266);

    auto tk_yyzz_yyyyy = pbuffer.data(idx_kin_gh + 267);

    auto tk_yyzz_yyyyz = pbuffer.data(idx_kin_gh + 268);

    auto tk_yyzz_yyyzz = pbuffer.data(idx_kin_gh + 269);

    auto tk_yyzz_yyzzz = pbuffer.data(idx_kin_gh + 270);

    auto tk_yyzz_yzzzz = pbuffer.data(idx_kin_gh + 271);

    auto tk_yyzz_zzzzz = pbuffer.data(idx_kin_gh + 272);

    auto tk_yzzz_xxxxx = pbuffer.data(idx_kin_gh + 273);

    auto tk_yzzz_xxxxz = pbuffer.data(idx_kin_gh + 275);

    auto tk_yzzz_xxxyz = pbuffer.data(idx_kin_gh + 277);

    auto tk_yzzz_xxxzz = pbuffer.data(idx_kin_gh + 278);

    auto tk_yzzz_xxyyz = pbuffer.data(idx_kin_gh + 280);

    auto tk_yzzz_xxyzz = pbuffer.data(idx_kin_gh + 281);

    auto tk_yzzz_xxzzz = pbuffer.data(idx_kin_gh + 282);

    auto tk_yzzz_xyyyz = pbuffer.data(idx_kin_gh + 284);

    auto tk_yzzz_xyyzz = pbuffer.data(idx_kin_gh + 285);

    auto tk_yzzz_xyzzz = pbuffer.data(idx_kin_gh + 286);

    auto tk_yzzz_xzzzz = pbuffer.data(idx_kin_gh + 287);

    auto tk_yzzz_yyyyz = pbuffer.data(idx_kin_gh + 289);

    auto tk_yzzz_yyyzz = pbuffer.data(idx_kin_gh + 290);

    auto tk_yzzz_yyzzz = pbuffer.data(idx_kin_gh + 291);

    auto tk_yzzz_yzzzz = pbuffer.data(idx_kin_gh + 292);

    auto tk_yzzz_zzzzz = pbuffer.data(idx_kin_gh + 293);

    auto tk_zzzz_xxxxx = pbuffer.data(idx_kin_gh + 294);

    auto tk_zzzz_xxxxy = pbuffer.data(idx_kin_gh + 295);

    auto tk_zzzz_xxxxz = pbuffer.data(idx_kin_gh + 296);

    auto tk_zzzz_xxxyy = pbuffer.data(idx_kin_gh + 297);

    auto tk_zzzz_xxxyz = pbuffer.data(idx_kin_gh + 298);

    auto tk_zzzz_xxxzz = pbuffer.data(idx_kin_gh + 299);

    auto tk_zzzz_xxyyy = pbuffer.data(idx_kin_gh + 300);

    auto tk_zzzz_xxyyz = pbuffer.data(idx_kin_gh + 301);

    auto tk_zzzz_xxyzz = pbuffer.data(idx_kin_gh + 302);

    auto tk_zzzz_xxzzz = pbuffer.data(idx_kin_gh + 303);

    auto tk_zzzz_xyyyy = pbuffer.data(idx_kin_gh + 304);

    auto tk_zzzz_xyyyz = pbuffer.data(idx_kin_gh + 305);

    auto tk_zzzz_xyyzz = pbuffer.data(idx_kin_gh + 306);

    auto tk_zzzz_xyzzz = pbuffer.data(idx_kin_gh + 307);

    auto tk_zzzz_xzzzz = pbuffer.data(idx_kin_gh + 308);

    auto tk_zzzz_yyyyy = pbuffer.data(idx_kin_gh + 309);

    auto tk_zzzz_yyyyz = pbuffer.data(idx_kin_gh + 310);

    auto tk_zzzz_yyyzz = pbuffer.data(idx_kin_gh + 311);

    auto tk_zzzz_yyzzz = pbuffer.data(idx_kin_gh + 312);

    auto tk_zzzz_yzzzz = pbuffer.data(idx_kin_gh + 313);

    auto tk_zzzz_zzzzz = pbuffer.data(idx_kin_gh + 314);

    // Set up components of auxiliary buffer : HG

    auto tk_xxxxx_xxxx = pbuffer.data(idx_kin_hg);

    auto tk_xxxxx_xxxy = pbuffer.data(idx_kin_hg + 1);

    auto tk_xxxxx_xxxz = pbuffer.data(idx_kin_hg + 2);

    auto tk_xxxxx_xxyy = pbuffer.data(idx_kin_hg + 3);

    auto tk_xxxxx_xxyz = pbuffer.data(idx_kin_hg + 4);

    auto tk_xxxxx_xxzz = pbuffer.data(idx_kin_hg + 5);

    auto tk_xxxxx_xyyy = pbuffer.data(idx_kin_hg + 6);

    auto tk_xxxxx_xyyz = pbuffer.data(idx_kin_hg + 7);

    auto tk_xxxxx_xyzz = pbuffer.data(idx_kin_hg + 8);

    auto tk_xxxxx_xzzz = pbuffer.data(idx_kin_hg + 9);

    auto tk_xxxxx_yyyy = pbuffer.data(idx_kin_hg + 10);

    auto tk_xxxxx_yyyz = pbuffer.data(idx_kin_hg + 11);

    auto tk_xxxxx_yyzz = pbuffer.data(idx_kin_hg + 12);

    auto tk_xxxxx_yzzz = pbuffer.data(idx_kin_hg + 13);

    auto tk_xxxxx_zzzz = pbuffer.data(idx_kin_hg + 14);

    auto tk_xxxxz_xxxz = pbuffer.data(idx_kin_hg + 32);

    auto tk_xxxxz_xxyz = pbuffer.data(idx_kin_hg + 34);

    auto tk_xxxxz_xxzz = pbuffer.data(idx_kin_hg + 35);

    auto tk_xxxxz_xyyz = pbuffer.data(idx_kin_hg + 37);

    auto tk_xxxxz_xyzz = pbuffer.data(idx_kin_hg + 38);

    auto tk_xxxxz_xzzz = pbuffer.data(idx_kin_hg + 39);

    auto tk_xxxxz_yyyz = pbuffer.data(idx_kin_hg + 41);

    auto tk_xxxxz_yyzz = pbuffer.data(idx_kin_hg + 42);

    auto tk_xxxxz_yzzz = pbuffer.data(idx_kin_hg + 43);

    auto tk_xxxxz_zzzz = pbuffer.data(idx_kin_hg + 44);

    auto tk_xxxyy_xxxx = pbuffer.data(idx_kin_hg + 45);

    auto tk_xxxyy_xxxy = pbuffer.data(idx_kin_hg + 46);

    auto tk_xxxyy_xxxz = pbuffer.data(idx_kin_hg + 47);

    auto tk_xxxyy_xxyy = pbuffer.data(idx_kin_hg + 48);

    auto tk_xxxyy_xxyz = pbuffer.data(idx_kin_hg + 49);

    auto tk_xxxyy_xxzz = pbuffer.data(idx_kin_hg + 50);

    auto tk_xxxyy_xyyy = pbuffer.data(idx_kin_hg + 51);

    auto tk_xxxyy_xyyz = pbuffer.data(idx_kin_hg + 52);

    auto tk_xxxyy_xyzz = pbuffer.data(idx_kin_hg + 53);

    auto tk_xxxyy_xzzz = pbuffer.data(idx_kin_hg + 54);

    auto tk_xxxyy_yyyy = pbuffer.data(idx_kin_hg + 55);

    auto tk_xxxyy_yyyz = pbuffer.data(idx_kin_hg + 56);

    auto tk_xxxyy_yyzz = pbuffer.data(idx_kin_hg + 57);

    auto tk_xxxyy_yzzz = pbuffer.data(idx_kin_hg + 58);

    auto tk_xxxyy_zzzz = pbuffer.data(idx_kin_hg + 59);

    auto tk_xxxzz_xxxx = pbuffer.data(idx_kin_hg + 75);

    auto tk_xxxzz_xxxy = pbuffer.data(idx_kin_hg + 76);

    auto tk_xxxzz_xxxz = pbuffer.data(idx_kin_hg + 77);

    auto tk_xxxzz_xxyy = pbuffer.data(idx_kin_hg + 78);

    auto tk_xxxzz_xxyz = pbuffer.data(idx_kin_hg + 79);

    auto tk_xxxzz_xxzz = pbuffer.data(idx_kin_hg + 80);

    auto tk_xxxzz_xyyy = pbuffer.data(idx_kin_hg + 81);

    auto tk_xxxzz_xyyz = pbuffer.data(idx_kin_hg + 82);

    auto tk_xxxzz_xyzz = pbuffer.data(idx_kin_hg + 83);

    auto tk_xxxzz_xzzz = pbuffer.data(idx_kin_hg + 84);

    auto tk_xxxzz_yyyy = pbuffer.data(idx_kin_hg + 85);

    auto tk_xxxzz_yyyz = pbuffer.data(idx_kin_hg + 86);

    auto tk_xxxzz_yyzz = pbuffer.data(idx_kin_hg + 87);

    auto tk_xxxzz_yzzz = pbuffer.data(idx_kin_hg + 88);

    auto tk_xxxzz_zzzz = pbuffer.data(idx_kin_hg + 89);

    auto tk_xxyyy_xxxx = pbuffer.data(idx_kin_hg + 90);

    auto tk_xxyyy_xxxy = pbuffer.data(idx_kin_hg + 91);

    auto tk_xxyyy_xxxz = pbuffer.data(idx_kin_hg + 92);

    auto tk_xxyyy_xxyy = pbuffer.data(idx_kin_hg + 93);

    auto tk_xxyyy_xxyz = pbuffer.data(idx_kin_hg + 94);

    auto tk_xxyyy_xxzz = pbuffer.data(idx_kin_hg + 95);

    auto tk_xxyyy_xyyy = pbuffer.data(idx_kin_hg + 96);

    auto tk_xxyyy_xyyz = pbuffer.data(idx_kin_hg + 97);

    auto tk_xxyyy_xyzz = pbuffer.data(idx_kin_hg + 98);

    auto tk_xxyyy_xzzz = pbuffer.data(idx_kin_hg + 99);

    auto tk_xxyyy_yyyy = pbuffer.data(idx_kin_hg + 100);

    auto tk_xxyyy_yyyz = pbuffer.data(idx_kin_hg + 101);

    auto tk_xxyyy_yyzz = pbuffer.data(idx_kin_hg + 102);

    auto tk_xxyyy_yzzz = pbuffer.data(idx_kin_hg + 103);

    auto tk_xxyyy_zzzz = pbuffer.data(idx_kin_hg + 104);

    auto tk_xxzzz_xxxx = pbuffer.data(idx_kin_hg + 135);

    auto tk_xxzzz_xxxy = pbuffer.data(idx_kin_hg + 136);

    auto tk_xxzzz_xxxz = pbuffer.data(idx_kin_hg + 137);

    auto tk_xxzzz_xxyy = pbuffer.data(idx_kin_hg + 138);

    auto tk_xxzzz_xxyz = pbuffer.data(idx_kin_hg + 139);

    auto tk_xxzzz_xxzz = pbuffer.data(idx_kin_hg + 140);

    auto tk_xxzzz_xyyy = pbuffer.data(idx_kin_hg + 141);

    auto tk_xxzzz_xyyz = pbuffer.data(idx_kin_hg + 142);

    auto tk_xxzzz_xyzz = pbuffer.data(idx_kin_hg + 143);

    auto tk_xxzzz_xzzz = pbuffer.data(idx_kin_hg + 144);

    auto tk_xxzzz_yyyy = pbuffer.data(idx_kin_hg + 145);

    auto tk_xxzzz_yyyz = pbuffer.data(idx_kin_hg + 146);

    auto tk_xxzzz_yyzz = pbuffer.data(idx_kin_hg + 147);

    auto tk_xxzzz_yzzz = pbuffer.data(idx_kin_hg + 148);

    auto tk_xxzzz_zzzz = pbuffer.data(idx_kin_hg + 149);

    auto tk_xyyyy_xxxy = pbuffer.data(idx_kin_hg + 151);

    auto tk_xyyyy_xxyy = pbuffer.data(idx_kin_hg + 153);

    auto tk_xyyyy_xxyz = pbuffer.data(idx_kin_hg + 154);

    auto tk_xyyyy_xyyy = pbuffer.data(idx_kin_hg + 156);

    auto tk_xyyyy_xyyz = pbuffer.data(idx_kin_hg + 157);

    auto tk_xyyyy_xyzz = pbuffer.data(idx_kin_hg + 158);

    auto tk_xyyyy_yyyy = pbuffer.data(idx_kin_hg + 160);

    auto tk_xyyyy_yyyz = pbuffer.data(idx_kin_hg + 161);

    auto tk_xyyyy_yyzz = pbuffer.data(idx_kin_hg + 162);

    auto tk_xyyyy_yzzz = pbuffer.data(idx_kin_hg + 163);

    auto tk_xyyzz_xxyz = pbuffer.data(idx_kin_hg + 184);

    auto tk_xyyzz_xyyz = pbuffer.data(idx_kin_hg + 187);

    auto tk_xyyzz_xyzz = pbuffer.data(idx_kin_hg + 188);

    auto tk_xyyzz_yyyz = pbuffer.data(idx_kin_hg + 191);

    auto tk_xyyzz_yyzz = pbuffer.data(idx_kin_hg + 192);

    auto tk_xyyzz_yzzz = pbuffer.data(idx_kin_hg + 193);

    auto tk_xzzzz_xxxz = pbuffer.data(idx_kin_hg + 212);

    auto tk_xzzzz_xxyz = pbuffer.data(idx_kin_hg + 214);

    auto tk_xzzzz_xxzz = pbuffer.data(idx_kin_hg + 215);

    auto tk_xzzzz_xyyz = pbuffer.data(idx_kin_hg + 217);

    auto tk_xzzzz_xyzz = pbuffer.data(idx_kin_hg + 218);

    auto tk_xzzzz_xzzz = pbuffer.data(idx_kin_hg + 219);

    auto tk_xzzzz_yyyz = pbuffer.data(idx_kin_hg + 221);

    auto tk_xzzzz_yyzz = pbuffer.data(idx_kin_hg + 222);

    auto tk_xzzzz_yzzz = pbuffer.data(idx_kin_hg + 223);

    auto tk_xzzzz_zzzz = pbuffer.data(idx_kin_hg + 224);

    auto tk_yyyyy_xxxx = pbuffer.data(idx_kin_hg + 225);

    auto tk_yyyyy_xxxy = pbuffer.data(idx_kin_hg + 226);

    auto tk_yyyyy_xxxz = pbuffer.data(idx_kin_hg + 227);

    auto tk_yyyyy_xxyy = pbuffer.data(idx_kin_hg + 228);

    auto tk_yyyyy_xxyz = pbuffer.data(idx_kin_hg + 229);

    auto tk_yyyyy_xxzz = pbuffer.data(idx_kin_hg + 230);

    auto tk_yyyyy_xyyy = pbuffer.data(idx_kin_hg + 231);

    auto tk_yyyyy_xyyz = pbuffer.data(idx_kin_hg + 232);

    auto tk_yyyyy_xyzz = pbuffer.data(idx_kin_hg + 233);

    auto tk_yyyyy_xzzz = pbuffer.data(idx_kin_hg + 234);

    auto tk_yyyyy_yyyy = pbuffer.data(idx_kin_hg + 235);

    auto tk_yyyyy_yyyz = pbuffer.data(idx_kin_hg + 236);

    auto tk_yyyyy_yyzz = pbuffer.data(idx_kin_hg + 237);

    auto tk_yyyyy_yzzz = pbuffer.data(idx_kin_hg + 238);

    auto tk_yyyyy_zzzz = pbuffer.data(idx_kin_hg + 239);

    auto tk_yyyyz_xxxz = pbuffer.data(idx_kin_hg + 242);

    auto tk_yyyyz_xxyz = pbuffer.data(idx_kin_hg + 244);

    auto tk_yyyyz_xxzz = pbuffer.data(idx_kin_hg + 245);

    auto tk_yyyyz_xyyz = pbuffer.data(idx_kin_hg + 247);

    auto tk_yyyyz_xyzz = pbuffer.data(idx_kin_hg + 248);

    auto tk_yyyyz_xzzz = pbuffer.data(idx_kin_hg + 249);

    auto tk_yyyyz_yyyz = pbuffer.data(idx_kin_hg + 251);

    auto tk_yyyyz_yyzz = pbuffer.data(idx_kin_hg + 252);

    auto tk_yyyyz_yzzz = pbuffer.data(idx_kin_hg + 253);

    auto tk_yyyyz_zzzz = pbuffer.data(idx_kin_hg + 254);

    auto tk_yyyzz_xxxx = pbuffer.data(idx_kin_hg + 255);

    auto tk_yyyzz_xxxy = pbuffer.data(idx_kin_hg + 256);

    auto tk_yyyzz_xxxz = pbuffer.data(idx_kin_hg + 257);

    auto tk_yyyzz_xxyy = pbuffer.data(idx_kin_hg + 258);

    auto tk_yyyzz_xxyz = pbuffer.data(idx_kin_hg + 259);

    auto tk_yyyzz_xxzz = pbuffer.data(idx_kin_hg + 260);

    auto tk_yyyzz_xyyy = pbuffer.data(idx_kin_hg + 261);

    auto tk_yyyzz_xyyz = pbuffer.data(idx_kin_hg + 262);

    auto tk_yyyzz_xyzz = pbuffer.data(idx_kin_hg + 263);

    auto tk_yyyzz_xzzz = pbuffer.data(idx_kin_hg + 264);

    auto tk_yyyzz_yyyy = pbuffer.data(idx_kin_hg + 265);

    auto tk_yyyzz_yyyz = pbuffer.data(idx_kin_hg + 266);

    auto tk_yyyzz_yyzz = pbuffer.data(idx_kin_hg + 267);

    auto tk_yyyzz_yzzz = pbuffer.data(idx_kin_hg + 268);

    auto tk_yyyzz_zzzz = pbuffer.data(idx_kin_hg + 269);

    auto tk_yyzzz_xxxx = pbuffer.data(idx_kin_hg + 270);

    auto tk_yyzzz_xxxy = pbuffer.data(idx_kin_hg + 271);

    auto tk_yyzzz_xxxz = pbuffer.data(idx_kin_hg + 272);

    auto tk_yyzzz_xxyy = pbuffer.data(idx_kin_hg + 273);

    auto tk_yyzzz_xxyz = pbuffer.data(idx_kin_hg + 274);

    auto tk_yyzzz_xxzz = pbuffer.data(idx_kin_hg + 275);

    auto tk_yyzzz_xyyy = pbuffer.data(idx_kin_hg + 276);

    auto tk_yyzzz_xyyz = pbuffer.data(idx_kin_hg + 277);

    auto tk_yyzzz_xyzz = pbuffer.data(idx_kin_hg + 278);

    auto tk_yyzzz_xzzz = pbuffer.data(idx_kin_hg + 279);

    auto tk_yyzzz_yyyy = pbuffer.data(idx_kin_hg + 280);

    auto tk_yyzzz_yyyz = pbuffer.data(idx_kin_hg + 281);

    auto tk_yyzzz_yyzz = pbuffer.data(idx_kin_hg + 282);

    auto tk_yyzzz_yzzz = pbuffer.data(idx_kin_hg + 283);

    auto tk_yyzzz_zzzz = pbuffer.data(idx_kin_hg + 284);

    auto tk_yzzzz_xxxy = pbuffer.data(idx_kin_hg + 286);

    auto tk_yzzzz_xxxz = pbuffer.data(idx_kin_hg + 287);

    auto tk_yzzzz_xxyy = pbuffer.data(idx_kin_hg + 288);

    auto tk_yzzzz_xxyz = pbuffer.data(idx_kin_hg + 289);

    auto tk_yzzzz_xxzz = pbuffer.data(idx_kin_hg + 290);

    auto tk_yzzzz_xyyy = pbuffer.data(idx_kin_hg + 291);

    auto tk_yzzzz_xyyz = pbuffer.data(idx_kin_hg + 292);

    auto tk_yzzzz_xyzz = pbuffer.data(idx_kin_hg + 293);

    auto tk_yzzzz_xzzz = pbuffer.data(idx_kin_hg + 294);

    auto tk_yzzzz_yyyy = pbuffer.data(idx_kin_hg + 295);

    auto tk_yzzzz_yyyz = pbuffer.data(idx_kin_hg + 296);

    auto tk_yzzzz_yyzz = pbuffer.data(idx_kin_hg + 297);

    auto tk_yzzzz_yzzz = pbuffer.data(idx_kin_hg + 298);

    auto tk_yzzzz_zzzz = pbuffer.data(idx_kin_hg + 299);

    auto tk_zzzzz_xxxx = pbuffer.data(idx_kin_hg + 300);

    auto tk_zzzzz_xxxy = pbuffer.data(idx_kin_hg + 301);

    auto tk_zzzzz_xxxz = pbuffer.data(idx_kin_hg + 302);

    auto tk_zzzzz_xxyy = pbuffer.data(idx_kin_hg + 303);

    auto tk_zzzzz_xxyz = pbuffer.data(idx_kin_hg + 304);

    auto tk_zzzzz_xxzz = pbuffer.data(idx_kin_hg + 305);

    auto tk_zzzzz_xyyy = pbuffer.data(idx_kin_hg + 306);

    auto tk_zzzzz_xyyz = pbuffer.data(idx_kin_hg + 307);

    auto tk_zzzzz_xyzz = pbuffer.data(idx_kin_hg + 308);

    auto tk_zzzzz_xzzz = pbuffer.data(idx_kin_hg + 309);

    auto tk_zzzzz_yyyy = pbuffer.data(idx_kin_hg + 310);

    auto tk_zzzzz_yyyz = pbuffer.data(idx_kin_hg + 311);

    auto tk_zzzzz_yyzz = pbuffer.data(idx_kin_hg + 312);

    auto tk_zzzzz_yzzz = pbuffer.data(idx_kin_hg + 313);

    auto tk_zzzzz_zzzz = pbuffer.data(idx_kin_hg + 314);

    // Set up components of auxiliary buffer : HH

    auto tk_xxxxx_xxxxx = pbuffer.data(idx_kin_hh);

    auto tk_xxxxx_xxxxy = pbuffer.data(idx_kin_hh + 1);

    auto tk_xxxxx_xxxxz = pbuffer.data(idx_kin_hh + 2);

    auto tk_xxxxx_xxxyy = pbuffer.data(idx_kin_hh + 3);

    auto tk_xxxxx_xxxyz = pbuffer.data(idx_kin_hh + 4);

    auto tk_xxxxx_xxxzz = pbuffer.data(idx_kin_hh + 5);

    auto tk_xxxxx_xxyyy = pbuffer.data(idx_kin_hh + 6);

    auto tk_xxxxx_xxyyz = pbuffer.data(idx_kin_hh + 7);

    auto tk_xxxxx_xxyzz = pbuffer.data(idx_kin_hh + 8);

    auto tk_xxxxx_xxzzz = pbuffer.data(idx_kin_hh + 9);

    auto tk_xxxxx_xyyyy = pbuffer.data(idx_kin_hh + 10);

    auto tk_xxxxx_xyyyz = pbuffer.data(idx_kin_hh + 11);

    auto tk_xxxxx_xyyzz = pbuffer.data(idx_kin_hh + 12);

    auto tk_xxxxx_xyzzz = pbuffer.data(idx_kin_hh + 13);

    auto tk_xxxxx_xzzzz = pbuffer.data(idx_kin_hh + 14);

    auto tk_xxxxx_yyyyy = pbuffer.data(idx_kin_hh + 15);

    auto tk_xxxxx_yyyyz = pbuffer.data(idx_kin_hh + 16);

    auto tk_xxxxx_yyyzz = pbuffer.data(idx_kin_hh + 17);

    auto tk_xxxxx_yyzzz = pbuffer.data(idx_kin_hh + 18);

    auto tk_xxxxx_yzzzz = pbuffer.data(idx_kin_hh + 19);

    auto tk_xxxxx_zzzzz = pbuffer.data(idx_kin_hh + 20);

    auto tk_xxxxy_xxxxx = pbuffer.data(idx_kin_hh + 21);

    auto tk_xxxxy_xxxxy = pbuffer.data(idx_kin_hh + 22);

    auto tk_xxxxy_xxxxz = pbuffer.data(idx_kin_hh + 23);

    auto tk_xxxxy_xxxyy = pbuffer.data(idx_kin_hh + 24);

    auto tk_xxxxy_xxxzz = pbuffer.data(idx_kin_hh + 26);

    auto tk_xxxxy_xxyyy = pbuffer.data(idx_kin_hh + 27);

    auto tk_xxxxy_xxzzz = pbuffer.data(idx_kin_hh + 30);

    auto tk_xxxxy_xyyyy = pbuffer.data(idx_kin_hh + 31);

    auto tk_xxxxy_xzzzz = pbuffer.data(idx_kin_hh + 35);

    auto tk_xxxxy_yyyyy = pbuffer.data(idx_kin_hh + 36);

    auto tk_xxxxz_xxxxx = pbuffer.data(idx_kin_hh + 42);

    auto tk_xxxxz_xxxxy = pbuffer.data(idx_kin_hh + 43);

    auto tk_xxxxz_xxxxz = pbuffer.data(idx_kin_hh + 44);

    auto tk_xxxxz_xxxyy = pbuffer.data(idx_kin_hh + 45);

    auto tk_xxxxz_xxxyz = pbuffer.data(idx_kin_hh + 46);

    auto tk_xxxxz_xxxzz = pbuffer.data(idx_kin_hh + 47);

    auto tk_xxxxz_xxyyy = pbuffer.data(idx_kin_hh + 48);

    auto tk_xxxxz_xxyyz = pbuffer.data(idx_kin_hh + 49);

    auto tk_xxxxz_xxyzz = pbuffer.data(idx_kin_hh + 50);

    auto tk_xxxxz_xxzzz = pbuffer.data(idx_kin_hh + 51);

    auto tk_xxxxz_xyyyy = pbuffer.data(idx_kin_hh + 52);

    auto tk_xxxxz_xyyyz = pbuffer.data(idx_kin_hh + 53);

    auto tk_xxxxz_xyyzz = pbuffer.data(idx_kin_hh + 54);

    auto tk_xxxxz_xyzzz = pbuffer.data(idx_kin_hh + 55);

    auto tk_xxxxz_xzzzz = pbuffer.data(idx_kin_hh + 56);

    auto tk_xxxxz_yyyyz = pbuffer.data(idx_kin_hh + 58);

    auto tk_xxxxz_yyyzz = pbuffer.data(idx_kin_hh + 59);

    auto tk_xxxxz_yyzzz = pbuffer.data(idx_kin_hh + 60);

    auto tk_xxxxz_yzzzz = pbuffer.data(idx_kin_hh + 61);

    auto tk_xxxxz_zzzzz = pbuffer.data(idx_kin_hh + 62);

    auto tk_xxxyy_xxxxx = pbuffer.data(idx_kin_hh + 63);

    auto tk_xxxyy_xxxxy = pbuffer.data(idx_kin_hh + 64);

    auto tk_xxxyy_xxxxz = pbuffer.data(idx_kin_hh + 65);

    auto tk_xxxyy_xxxyy = pbuffer.data(idx_kin_hh + 66);

    auto tk_xxxyy_xxxyz = pbuffer.data(idx_kin_hh + 67);

    auto tk_xxxyy_xxxzz = pbuffer.data(idx_kin_hh + 68);

    auto tk_xxxyy_xxyyy = pbuffer.data(idx_kin_hh + 69);

    auto tk_xxxyy_xxyyz = pbuffer.data(idx_kin_hh + 70);

    auto tk_xxxyy_xxyzz = pbuffer.data(idx_kin_hh + 71);

    auto tk_xxxyy_xxzzz = pbuffer.data(idx_kin_hh + 72);

    auto tk_xxxyy_xyyyy = pbuffer.data(idx_kin_hh + 73);

    auto tk_xxxyy_xyyyz = pbuffer.data(idx_kin_hh + 74);

    auto tk_xxxyy_xyyzz = pbuffer.data(idx_kin_hh + 75);

    auto tk_xxxyy_xyzzz = pbuffer.data(idx_kin_hh + 76);

    auto tk_xxxyy_xzzzz = pbuffer.data(idx_kin_hh + 77);

    auto tk_xxxyy_yyyyy = pbuffer.data(idx_kin_hh + 78);

    auto tk_xxxyy_yyyyz = pbuffer.data(idx_kin_hh + 79);

    auto tk_xxxyy_yyyzz = pbuffer.data(idx_kin_hh + 80);

    auto tk_xxxyy_yyzzz = pbuffer.data(idx_kin_hh + 81);

    auto tk_xxxyy_yzzzz = pbuffer.data(idx_kin_hh + 82);

    auto tk_xxxyy_zzzzz = pbuffer.data(idx_kin_hh + 83);

    auto tk_xxxzz_xxxxx = pbuffer.data(idx_kin_hh + 105);

    auto tk_xxxzz_xxxxy = pbuffer.data(idx_kin_hh + 106);

    auto tk_xxxzz_xxxxz = pbuffer.data(idx_kin_hh + 107);

    auto tk_xxxzz_xxxyy = pbuffer.data(idx_kin_hh + 108);

    auto tk_xxxzz_xxxyz = pbuffer.data(idx_kin_hh + 109);

    auto tk_xxxzz_xxxzz = pbuffer.data(idx_kin_hh + 110);

    auto tk_xxxzz_xxyyy = pbuffer.data(idx_kin_hh + 111);

    auto tk_xxxzz_xxyyz = pbuffer.data(idx_kin_hh + 112);

    auto tk_xxxzz_xxyzz = pbuffer.data(idx_kin_hh + 113);

    auto tk_xxxzz_xxzzz = pbuffer.data(idx_kin_hh + 114);

    auto tk_xxxzz_xyyyy = pbuffer.data(idx_kin_hh + 115);

    auto tk_xxxzz_xyyyz = pbuffer.data(idx_kin_hh + 116);

    auto tk_xxxzz_xyyzz = pbuffer.data(idx_kin_hh + 117);

    auto tk_xxxzz_xyzzz = pbuffer.data(idx_kin_hh + 118);

    auto tk_xxxzz_xzzzz = pbuffer.data(idx_kin_hh + 119);

    auto tk_xxxzz_yyyyy = pbuffer.data(idx_kin_hh + 120);

    auto tk_xxxzz_yyyyz = pbuffer.data(idx_kin_hh + 121);

    auto tk_xxxzz_yyyzz = pbuffer.data(idx_kin_hh + 122);

    auto tk_xxxzz_yyzzz = pbuffer.data(idx_kin_hh + 123);

    auto tk_xxxzz_yzzzz = pbuffer.data(idx_kin_hh + 124);

    auto tk_xxxzz_zzzzz = pbuffer.data(idx_kin_hh + 125);

    auto tk_xxyyy_xxxxx = pbuffer.data(idx_kin_hh + 126);

    auto tk_xxyyy_xxxxy = pbuffer.data(idx_kin_hh + 127);

    auto tk_xxyyy_xxxxz = pbuffer.data(idx_kin_hh + 128);

    auto tk_xxyyy_xxxyy = pbuffer.data(idx_kin_hh + 129);

    auto tk_xxyyy_xxxyz = pbuffer.data(idx_kin_hh + 130);

    auto tk_xxyyy_xxxzz = pbuffer.data(idx_kin_hh + 131);

    auto tk_xxyyy_xxyyy = pbuffer.data(idx_kin_hh + 132);

    auto tk_xxyyy_xxyyz = pbuffer.data(idx_kin_hh + 133);

    auto tk_xxyyy_xxyzz = pbuffer.data(idx_kin_hh + 134);

    auto tk_xxyyy_xxzzz = pbuffer.data(idx_kin_hh + 135);

    auto tk_xxyyy_xyyyy = pbuffer.data(idx_kin_hh + 136);

    auto tk_xxyyy_xyyyz = pbuffer.data(idx_kin_hh + 137);

    auto tk_xxyyy_xyyzz = pbuffer.data(idx_kin_hh + 138);

    auto tk_xxyyy_xyzzz = pbuffer.data(idx_kin_hh + 139);

    auto tk_xxyyy_xzzzz = pbuffer.data(idx_kin_hh + 140);

    auto tk_xxyyy_yyyyy = pbuffer.data(idx_kin_hh + 141);

    auto tk_xxyyy_yyyyz = pbuffer.data(idx_kin_hh + 142);

    auto tk_xxyyy_yyyzz = pbuffer.data(idx_kin_hh + 143);

    auto tk_xxyyy_yyzzz = pbuffer.data(idx_kin_hh + 144);

    auto tk_xxyyy_yzzzz = pbuffer.data(idx_kin_hh + 145);

    auto tk_xxyyy_zzzzz = pbuffer.data(idx_kin_hh + 146);

    auto tk_xxyyz_xxxxy = pbuffer.data(idx_kin_hh + 148);

    auto tk_xxyyz_xxxyy = pbuffer.data(idx_kin_hh + 150);

    auto tk_xxyyz_xxyyy = pbuffer.data(idx_kin_hh + 153);

    auto tk_xxyyz_xyyyy = pbuffer.data(idx_kin_hh + 157);

    auto tk_xxyzz_xxxxx = pbuffer.data(idx_kin_hh + 168);

    auto tk_xxyzz_xxxxz = pbuffer.data(idx_kin_hh + 170);

    auto tk_xxyzz_xxxzz = pbuffer.data(idx_kin_hh + 173);

    auto tk_xxyzz_xxzzz = pbuffer.data(idx_kin_hh + 177);

    auto tk_xxyzz_xzzzz = pbuffer.data(idx_kin_hh + 182);

    auto tk_xxzzz_xxxxx = pbuffer.data(idx_kin_hh + 189);

    auto tk_xxzzz_xxxxy = pbuffer.data(idx_kin_hh + 190);

    auto tk_xxzzz_xxxxz = pbuffer.data(idx_kin_hh + 191);

    auto tk_xxzzz_xxxyy = pbuffer.data(idx_kin_hh + 192);

    auto tk_xxzzz_xxxyz = pbuffer.data(idx_kin_hh + 193);

    auto tk_xxzzz_xxxzz = pbuffer.data(idx_kin_hh + 194);

    auto tk_xxzzz_xxyyy = pbuffer.data(idx_kin_hh + 195);

    auto tk_xxzzz_xxyyz = pbuffer.data(idx_kin_hh + 196);

    auto tk_xxzzz_xxyzz = pbuffer.data(idx_kin_hh + 197);

    auto tk_xxzzz_xxzzz = pbuffer.data(idx_kin_hh + 198);

    auto tk_xxzzz_xyyyy = pbuffer.data(idx_kin_hh + 199);

    auto tk_xxzzz_xyyyz = pbuffer.data(idx_kin_hh + 200);

    auto tk_xxzzz_xyyzz = pbuffer.data(idx_kin_hh + 201);

    auto tk_xxzzz_xyzzz = pbuffer.data(idx_kin_hh + 202);

    auto tk_xxzzz_xzzzz = pbuffer.data(idx_kin_hh + 203);

    auto tk_xxzzz_yyyyy = pbuffer.data(idx_kin_hh + 204);

    auto tk_xxzzz_yyyyz = pbuffer.data(idx_kin_hh + 205);

    auto tk_xxzzz_yyyzz = pbuffer.data(idx_kin_hh + 206);

    auto tk_xxzzz_yyzzz = pbuffer.data(idx_kin_hh + 207);

    auto tk_xxzzz_yzzzz = pbuffer.data(idx_kin_hh + 208);

    auto tk_xxzzz_zzzzz = pbuffer.data(idx_kin_hh + 209);

    auto tk_xyyyy_xxxxx = pbuffer.data(idx_kin_hh + 210);

    auto tk_xyyyy_xxxxy = pbuffer.data(idx_kin_hh + 211);

    auto tk_xyyyy_xxxyy = pbuffer.data(idx_kin_hh + 213);

    auto tk_xyyyy_xxxyz = pbuffer.data(idx_kin_hh + 214);

    auto tk_xyyyy_xxyyy = pbuffer.data(idx_kin_hh + 216);

    auto tk_xyyyy_xxyyz = pbuffer.data(idx_kin_hh + 217);

    auto tk_xyyyy_xxyzz = pbuffer.data(idx_kin_hh + 218);

    auto tk_xyyyy_xyyyy = pbuffer.data(idx_kin_hh + 220);

    auto tk_xyyyy_xyyyz = pbuffer.data(idx_kin_hh + 221);

    auto tk_xyyyy_xyyzz = pbuffer.data(idx_kin_hh + 222);

    auto tk_xyyyy_xyzzz = pbuffer.data(idx_kin_hh + 223);

    auto tk_xyyyy_yyyyy = pbuffer.data(idx_kin_hh + 225);

    auto tk_xyyyy_yyyyz = pbuffer.data(idx_kin_hh + 226);

    auto tk_xyyyy_yyyzz = pbuffer.data(idx_kin_hh + 227);

    auto tk_xyyyy_yyzzz = pbuffer.data(idx_kin_hh + 228);

    auto tk_xyyyy_yzzzz = pbuffer.data(idx_kin_hh + 229);

    auto tk_xyyyy_zzzzz = pbuffer.data(idx_kin_hh + 230);

    auto tk_xyyzz_xxxyz = pbuffer.data(idx_kin_hh + 256);

    auto tk_xyyzz_xxyyz = pbuffer.data(idx_kin_hh + 259);

    auto tk_xyyzz_xxyzz = pbuffer.data(idx_kin_hh + 260);

    auto tk_xyyzz_xyyyz = pbuffer.data(idx_kin_hh + 263);

    auto tk_xyyzz_xyyzz = pbuffer.data(idx_kin_hh + 264);

    auto tk_xyyzz_xyzzz = pbuffer.data(idx_kin_hh + 265);

    auto tk_xyyzz_yyyyy = pbuffer.data(idx_kin_hh + 267);

    auto tk_xyyzz_yyyyz = pbuffer.data(idx_kin_hh + 268);

    auto tk_xyyzz_yyyzz = pbuffer.data(idx_kin_hh + 269);

    auto tk_xyyzz_yyzzz = pbuffer.data(idx_kin_hh + 270);

    auto tk_xyyzz_yzzzz = pbuffer.data(idx_kin_hh + 271);

    auto tk_xyyzz_zzzzz = pbuffer.data(idx_kin_hh + 272);

    auto tk_xzzzz_xxxxx = pbuffer.data(idx_kin_hh + 294);

    auto tk_xzzzz_xxxxz = pbuffer.data(idx_kin_hh + 296);

    auto tk_xzzzz_xxxyz = pbuffer.data(idx_kin_hh + 298);

    auto tk_xzzzz_xxxzz = pbuffer.data(idx_kin_hh + 299);

    auto tk_xzzzz_xxyyz = pbuffer.data(idx_kin_hh + 301);

    auto tk_xzzzz_xxyzz = pbuffer.data(idx_kin_hh + 302);

    auto tk_xzzzz_xxzzz = pbuffer.data(idx_kin_hh + 303);

    auto tk_xzzzz_xyyyz = pbuffer.data(idx_kin_hh + 305);

    auto tk_xzzzz_xyyzz = pbuffer.data(idx_kin_hh + 306);

    auto tk_xzzzz_xyzzz = pbuffer.data(idx_kin_hh + 307);

    auto tk_xzzzz_xzzzz = pbuffer.data(idx_kin_hh + 308);

    auto tk_xzzzz_yyyyy = pbuffer.data(idx_kin_hh + 309);

    auto tk_xzzzz_yyyyz = pbuffer.data(idx_kin_hh + 310);

    auto tk_xzzzz_yyyzz = pbuffer.data(idx_kin_hh + 311);

    auto tk_xzzzz_yyzzz = pbuffer.data(idx_kin_hh + 312);

    auto tk_xzzzz_yzzzz = pbuffer.data(idx_kin_hh + 313);

    auto tk_xzzzz_zzzzz = pbuffer.data(idx_kin_hh + 314);

    auto tk_yyyyy_xxxxx = pbuffer.data(idx_kin_hh + 315);

    auto tk_yyyyy_xxxxy = pbuffer.data(idx_kin_hh + 316);

    auto tk_yyyyy_xxxxz = pbuffer.data(idx_kin_hh + 317);

    auto tk_yyyyy_xxxyy = pbuffer.data(idx_kin_hh + 318);

    auto tk_yyyyy_xxxyz = pbuffer.data(idx_kin_hh + 319);

    auto tk_yyyyy_xxxzz = pbuffer.data(idx_kin_hh + 320);

    auto tk_yyyyy_xxyyy = pbuffer.data(idx_kin_hh + 321);

    auto tk_yyyyy_xxyyz = pbuffer.data(idx_kin_hh + 322);

    auto tk_yyyyy_xxyzz = pbuffer.data(idx_kin_hh + 323);

    auto tk_yyyyy_xxzzz = pbuffer.data(idx_kin_hh + 324);

    auto tk_yyyyy_xyyyy = pbuffer.data(idx_kin_hh + 325);

    auto tk_yyyyy_xyyyz = pbuffer.data(idx_kin_hh + 326);

    auto tk_yyyyy_xyyzz = pbuffer.data(idx_kin_hh + 327);

    auto tk_yyyyy_xyzzz = pbuffer.data(idx_kin_hh + 328);

    auto tk_yyyyy_xzzzz = pbuffer.data(idx_kin_hh + 329);

    auto tk_yyyyy_yyyyy = pbuffer.data(idx_kin_hh + 330);

    auto tk_yyyyy_yyyyz = pbuffer.data(idx_kin_hh + 331);

    auto tk_yyyyy_yyyzz = pbuffer.data(idx_kin_hh + 332);

    auto tk_yyyyy_yyzzz = pbuffer.data(idx_kin_hh + 333);

    auto tk_yyyyy_yzzzz = pbuffer.data(idx_kin_hh + 334);

    auto tk_yyyyy_zzzzz = pbuffer.data(idx_kin_hh + 335);

    auto tk_yyyyz_xxxxy = pbuffer.data(idx_kin_hh + 337);

    auto tk_yyyyz_xxxxz = pbuffer.data(idx_kin_hh + 338);

    auto tk_yyyyz_xxxyy = pbuffer.data(idx_kin_hh + 339);

    auto tk_yyyyz_xxxyz = pbuffer.data(idx_kin_hh + 340);

    auto tk_yyyyz_xxxzz = pbuffer.data(idx_kin_hh + 341);

    auto tk_yyyyz_xxyyy = pbuffer.data(idx_kin_hh + 342);

    auto tk_yyyyz_xxyyz = pbuffer.data(idx_kin_hh + 343);

    auto tk_yyyyz_xxyzz = pbuffer.data(idx_kin_hh + 344);

    auto tk_yyyyz_xxzzz = pbuffer.data(idx_kin_hh + 345);

    auto tk_yyyyz_xyyyy = pbuffer.data(idx_kin_hh + 346);

    auto tk_yyyyz_xyyyz = pbuffer.data(idx_kin_hh + 347);

    auto tk_yyyyz_xyyzz = pbuffer.data(idx_kin_hh + 348);

    auto tk_yyyyz_xyzzz = pbuffer.data(idx_kin_hh + 349);

    auto tk_yyyyz_xzzzz = pbuffer.data(idx_kin_hh + 350);

    auto tk_yyyyz_yyyyy = pbuffer.data(idx_kin_hh + 351);

    auto tk_yyyyz_yyyyz = pbuffer.data(idx_kin_hh + 352);

    auto tk_yyyyz_yyyzz = pbuffer.data(idx_kin_hh + 353);

    auto tk_yyyyz_yyzzz = pbuffer.data(idx_kin_hh + 354);

    auto tk_yyyyz_yzzzz = pbuffer.data(idx_kin_hh + 355);

    auto tk_yyyyz_zzzzz = pbuffer.data(idx_kin_hh + 356);

    auto tk_yyyzz_xxxxx = pbuffer.data(idx_kin_hh + 357);

    auto tk_yyyzz_xxxxy = pbuffer.data(idx_kin_hh + 358);

    auto tk_yyyzz_xxxxz = pbuffer.data(idx_kin_hh + 359);

    auto tk_yyyzz_xxxyy = pbuffer.data(idx_kin_hh + 360);

    auto tk_yyyzz_xxxyz = pbuffer.data(idx_kin_hh + 361);

    auto tk_yyyzz_xxxzz = pbuffer.data(idx_kin_hh + 362);

    auto tk_yyyzz_xxyyy = pbuffer.data(idx_kin_hh + 363);

    auto tk_yyyzz_xxyyz = pbuffer.data(idx_kin_hh + 364);

    auto tk_yyyzz_xxyzz = pbuffer.data(idx_kin_hh + 365);

    auto tk_yyyzz_xxzzz = pbuffer.data(idx_kin_hh + 366);

    auto tk_yyyzz_xyyyy = pbuffer.data(idx_kin_hh + 367);

    auto tk_yyyzz_xyyyz = pbuffer.data(idx_kin_hh + 368);

    auto tk_yyyzz_xyyzz = pbuffer.data(idx_kin_hh + 369);

    auto tk_yyyzz_xyzzz = pbuffer.data(idx_kin_hh + 370);

    auto tk_yyyzz_xzzzz = pbuffer.data(idx_kin_hh + 371);

    auto tk_yyyzz_yyyyy = pbuffer.data(idx_kin_hh + 372);

    auto tk_yyyzz_yyyyz = pbuffer.data(idx_kin_hh + 373);

    auto tk_yyyzz_yyyzz = pbuffer.data(idx_kin_hh + 374);

    auto tk_yyyzz_yyzzz = pbuffer.data(idx_kin_hh + 375);

    auto tk_yyyzz_yzzzz = pbuffer.data(idx_kin_hh + 376);

    auto tk_yyyzz_zzzzz = pbuffer.data(idx_kin_hh + 377);

    auto tk_yyzzz_xxxxx = pbuffer.data(idx_kin_hh + 378);

    auto tk_yyzzz_xxxxy = pbuffer.data(idx_kin_hh + 379);

    auto tk_yyzzz_xxxxz = pbuffer.data(idx_kin_hh + 380);

    auto tk_yyzzz_xxxyy = pbuffer.data(idx_kin_hh + 381);

    auto tk_yyzzz_xxxyz = pbuffer.data(idx_kin_hh + 382);

    auto tk_yyzzz_xxxzz = pbuffer.data(idx_kin_hh + 383);

    auto tk_yyzzz_xxyyy = pbuffer.data(idx_kin_hh + 384);

    auto tk_yyzzz_xxyyz = pbuffer.data(idx_kin_hh + 385);

    auto tk_yyzzz_xxyzz = pbuffer.data(idx_kin_hh + 386);

    auto tk_yyzzz_xxzzz = pbuffer.data(idx_kin_hh + 387);

    auto tk_yyzzz_xyyyy = pbuffer.data(idx_kin_hh + 388);

    auto tk_yyzzz_xyyyz = pbuffer.data(idx_kin_hh + 389);

    auto tk_yyzzz_xyyzz = pbuffer.data(idx_kin_hh + 390);

    auto tk_yyzzz_xyzzz = pbuffer.data(idx_kin_hh + 391);

    auto tk_yyzzz_xzzzz = pbuffer.data(idx_kin_hh + 392);

    auto tk_yyzzz_yyyyy = pbuffer.data(idx_kin_hh + 393);

    auto tk_yyzzz_yyyyz = pbuffer.data(idx_kin_hh + 394);

    auto tk_yyzzz_yyyzz = pbuffer.data(idx_kin_hh + 395);

    auto tk_yyzzz_yyzzz = pbuffer.data(idx_kin_hh + 396);

    auto tk_yyzzz_yzzzz = pbuffer.data(idx_kin_hh + 397);

    auto tk_yyzzz_zzzzz = pbuffer.data(idx_kin_hh + 398);

    auto tk_yzzzz_xxxxx = pbuffer.data(idx_kin_hh + 399);

    auto tk_yzzzz_xxxxy = pbuffer.data(idx_kin_hh + 400);

    auto tk_yzzzz_xxxxz = pbuffer.data(idx_kin_hh + 401);

    auto tk_yzzzz_xxxyy = pbuffer.data(idx_kin_hh + 402);

    auto tk_yzzzz_xxxyz = pbuffer.data(idx_kin_hh + 403);

    auto tk_yzzzz_xxxzz = pbuffer.data(idx_kin_hh + 404);

    auto tk_yzzzz_xxyyy = pbuffer.data(idx_kin_hh + 405);

    auto tk_yzzzz_xxyyz = pbuffer.data(idx_kin_hh + 406);

    auto tk_yzzzz_xxyzz = pbuffer.data(idx_kin_hh + 407);

    auto tk_yzzzz_xxzzz = pbuffer.data(idx_kin_hh + 408);

    auto tk_yzzzz_xyyyy = pbuffer.data(idx_kin_hh + 409);

    auto tk_yzzzz_xyyyz = pbuffer.data(idx_kin_hh + 410);

    auto tk_yzzzz_xyyzz = pbuffer.data(idx_kin_hh + 411);

    auto tk_yzzzz_xyzzz = pbuffer.data(idx_kin_hh + 412);

    auto tk_yzzzz_xzzzz = pbuffer.data(idx_kin_hh + 413);

    auto tk_yzzzz_yyyyy = pbuffer.data(idx_kin_hh + 414);

    auto tk_yzzzz_yyyyz = pbuffer.data(idx_kin_hh + 415);

    auto tk_yzzzz_yyyzz = pbuffer.data(idx_kin_hh + 416);

    auto tk_yzzzz_yyzzz = pbuffer.data(idx_kin_hh + 417);

    auto tk_yzzzz_yzzzz = pbuffer.data(idx_kin_hh + 418);

    auto tk_yzzzz_zzzzz = pbuffer.data(idx_kin_hh + 419);

    auto tk_zzzzz_xxxxx = pbuffer.data(idx_kin_hh + 420);

    auto tk_zzzzz_xxxxy = pbuffer.data(idx_kin_hh + 421);

    auto tk_zzzzz_xxxxz = pbuffer.data(idx_kin_hh + 422);

    auto tk_zzzzz_xxxyy = pbuffer.data(idx_kin_hh + 423);

    auto tk_zzzzz_xxxyz = pbuffer.data(idx_kin_hh + 424);

    auto tk_zzzzz_xxxzz = pbuffer.data(idx_kin_hh + 425);

    auto tk_zzzzz_xxyyy = pbuffer.data(idx_kin_hh + 426);

    auto tk_zzzzz_xxyyz = pbuffer.data(idx_kin_hh + 427);

    auto tk_zzzzz_xxyzz = pbuffer.data(idx_kin_hh + 428);

    auto tk_zzzzz_xxzzz = pbuffer.data(idx_kin_hh + 429);

    auto tk_zzzzz_xyyyy = pbuffer.data(idx_kin_hh + 430);

    auto tk_zzzzz_xyyyz = pbuffer.data(idx_kin_hh + 431);

    auto tk_zzzzz_xyyzz = pbuffer.data(idx_kin_hh + 432);

    auto tk_zzzzz_xyzzz = pbuffer.data(idx_kin_hh + 433);

    auto tk_zzzzz_xzzzz = pbuffer.data(idx_kin_hh + 434);

    auto tk_zzzzz_yyyyy = pbuffer.data(idx_kin_hh + 435);

    auto tk_zzzzz_yyyyz = pbuffer.data(idx_kin_hh + 436);

    auto tk_zzzzz_yyyzz = pbuffer.data(idx_kin_hh + 437);

    auto tk_zzzzz_yyzzz = pbuffer.data(idx_kin_hh + 438);

    auto tk_zzzzz_yzzzz = pbuffer.data(idx_kin_hh + 439);

    auto tk_zzzzz_zzzzz = pbuffer.data(idx_kin_hh + 440);

    // Set up components of auxiliary buffer : IH

    auto ts_xxxxxx_xxxxx = pbuffer.data(idx_ovl_ih);

    auto ts_xxxxxx_xxxxy = pbuffer.data(idx_ovl_ih + 1);

    auto ts_xxxxxx_xxxxz = pbuffer.data(idx_ovl_ih + 2);

    auto ts_xxxxxx_xxxyy = pbuffer.data(idx_ovl_ih + 3);

    auto ts_xxxxxx_xxxyz = pbuffer.data(idx_ovl_ih + 4);

    auto ts_xxxxxx_xxxzz = pbuffer.data(idx_ovl_ih + 5);

    auto ts_xxxxxx_xxyyy = pbuffer.data(idx_ovl_ih + 6);

    auto ts_xxxxxx_xxyyz = pbuffer.data(idx_ovl_ih + 7);

    auto ts_xxxxxx_xxyzz = pbuffer.data(idx_ovl_ih + 8);

    auto ts_xxxxxx_xxzzz = pbuffer.data(idx_ovl_ih + 9);

    auto ts_xxxxxx_xyyyy = pbuffer.data(idx_ovl_ih + 10);

    auto ts_xxxxxx_xyyyz = pbuffer.data(idx_ovl_ih + 11);

    auto ts_xxxxxx_xyyzz = pbuffer.data(idx_ovl_ih + 12);

    auto ts_xxxxxx_xyzzz = pbuffer.data(idx_ovl_ih + 13);

    auto ts_xxxxxx_xzzzz = pbuffer.data(idx_ovl_ih + 14);

    auto ts_xxxxxx_yyyyy = pbuffer.data(idx_ovl_ih + 15);

    auto ts_xxxxxx_yyyyz = pbuffer.data(idx_ovl_ih + 16);

    auto ts_xxxxxx_yyyzz = pbuffer.data(idx_ovl_ih + 17);

    auto ts_xxxxxx_yyzzz = pbuffer.data(idx_ovl_ih + 18);

    auto ts_xxxxxx_yzzzz = pbuffer.data(idx_ovl_ih + 19);

    auto ts_xxxxxx_zzzzz = pbuffer.data(idx_ovl_ih + 20);

    auto ts_xxxxxy_xxxxx = pbuffer.data(idx_ovl_ih + 21);

    auto ts_xxxxxy_xxxxy = pbuffer.data(idx_ovl_ih + 22);

    auto ts_xxxxxy_xxxxz = pbuffer.data(idx_ovl_ih + 23);

    auto ts_xxxxxy_xxxyy = pbuffer.data(idx_ovl_ih + 24);

    auto ts_xxxxxy_xxxyz = pbuffer.data(idx_ovl_ih + 25);

    auto ts_xxxxxy_xxxzz = pbuffer.data(idx_ovl_ih + 26);

    auto ts_xxxxxy_xxyyy = pbuffer.data(idx_ovl_ih + 27);

    auto ts_xxxxxy_xxyyz = pbuffer.data(idx_ovl_ih + 28);

    auto ts_xxxxxy_xxyzz = pbuffer.data(idx_ovl_ih + 29);

    auto ts_xxxxxy_xxzzz = pbuffer.data(idx_ovl_ih + 30);

    auto ts_xxxxxy_xyyyy = pbuffer.data(idx_ovl_ih + 31);

    auto ts_xxxxxy_xyyyz = pbuffer.data(idx_ovl_ih + 32);

    auto ts_xxxxxy_xyyzz = pbuffer.data(idx_ovl_ih + 33);

    auto ts_xxxxxy_xyzzz = pbuffer.data(idx_ovl_ih + 34);

    auto ts_xxxxxy_xzzzz = pbuffer.data(idx_ovl_ih + 35);

    auto ts_xxxxxy_yyyyy = pbuffer.data(idx_ovl_ih + 36);

    auto ts_xxxxxy_yyyyz = pbuffer.data(idx_ovl_ih + 37);

    auto ts_xxxxxy_yyyzz = pbuffer.data(idx_ovl_ih + 38);

    auto ts_xxxxxy_yyzzz = pbuffer.data(idx_ovl_ih + 39);

    auto ts_xxxxxy_yzzzz = pbuffer.data(idx_ovl_ih + 40);

    auto ts_xxxxxy_zzzzz = pbuffer.data(idx_ovl_ih + 41);

    auto ts_xxxxxz_xxxxx = pbuffer.data(idx_ovl_ih + 42);

    auto ts_xxxxxz_xxxxy = pbuffer.data(idx_ovl_ih + 43);

    auto ts_xxxxxz_xxxxz = pbuffer.data(idx_ovl_ih + 44);

    auto ts_xxxxxz_xxxyy = pbuffer.data(idx_ovl_ih + 45);

    auto ts_xxxxxz_xxxyz = pbuffer.data(idx_ovl_ih + 46);

    auto ts_xxxxxz_xxxzz = pbuffer.data(idx_ovl_ih + 47);

    auto ts_xxxxxz_xxyyy = pbuffer.data(idx_ovl_ih + 48);

    auto ts_xxxxxz_xxyyz = pbuffer.data(idx_ovl_ih + 49);

    auto ts_xxxxxz_xxyzz = pbuffer.data(idx_ovl_ih + 50);

    auto ts_xxxxxz_xxzzz = pbuffer.data(idx_ovl_ih + 51);

    auto ts_xxxxxz_xyyyy = pbuffer.data(idx_ovl_ih + 52);

    auto ts_xxxxxz_xyyyz = pbuffer.data(idx_ovl_ih + 53);

    auto ts_xxxxxz_xyyzz = pbuffer.data(idx_ovl_ih + 54);

    auto ts_xxxxxz_xyzzz = pbuffer.data(idx_ovl_ih + 55);

    auto ts_xxxxxz_xzzzz = pbuffer.data(idx_ovl_ih + 56);

    auto ts_xxxxxz_yyyyy = pbuffer.data(idx_ovl_ih + 57);

    auto ts_xxxxxz_yyyyz = pbuffer.data(idx_ovl_ih + 58);

    auto ts_xxxxxz_yyyzz = pbuffer.data(idx_ovl_ih + 59);

    auto ts_xxxxxz_yyzzz = pbuffer.data(idx_ovl_ih + 60);

    auto ts_xxxxxz_yzzzz = pbuffer.data(idx_ovl_ih + 61);

    auto ts_xxxxxz_zzzzz = pbuffer.data(idx_ovl_ih + 62);

    auto ts_xxxxyy_xxxxx = pbuffer.data(idx_ovl_ih + 63);

    auto ts_xxxxyy_xxxxy = pbuffer.data(idx_ovl_ih + 64);

    auto ts_xxxxyy_xxxxz = pbuffer.data(idx_ovl_ih + 65);

    auto ts_xxxxyy_xxxyy = pbuffer.data(idx_ovl_ih + 66);

    auto ts_xxxxyy_xxxyz = pbuffer.data(idx_ovl_ih + 67);

    auto ts_xxxxyy_xxxzz = pbuffer.data(idx_ovl_ih + 68);

    auto ts_xxxxyy_xxyyy = pbuffer.data(idx_ovl_ih + 69);

    auto ts_xxxxyy_xxyyz = pbuffer.data(idx_ovl_ih + 70);

    auto ts_xxxxyy_xxyzz = pbuffer.data(idx_ovl_ih + 71);

    auto ts_xxxxyy_xxzzz = pbuffer.data(idx_ovl_ih + 72);

    auto ts_xxxxyy_xyyyy = pbuffer.data(idx_ovl_ih + 73);

    auto ts_xxxxyy_xyyyz = pbuffer.data(idx_ovl_ih + 74);

    auto ts_xxxxyy_xyyzz = pbuffer.data(idx_ovl_ih + 75);

    auto ts_xxxxyy_xyzzz = pbuffer.data(idx_ovl_ih + 76);

    auto ts_xxxxyy_xzzzz = pbuffer.data(idx_ovl_ih + 77);

    auto ts_xxxxyy_yyyyy = pbuffer.data(idx_ovl_ih + 78);

    auto ts_xxxxyy_yyyyz = pbuffer.data(idx_ovl_ih + 79);

    auto ts_xxxxyy_yyyzz = pbuffer.data(idx_ovl_ih + 80);

    auto ts_xxxxyy_yyzzz = pbuffer.data(idx_ovl_ih + 81);

    auto ts_xxxxyy_yzzzz = pbuffer.data(idx_ovl_ih + 82);

    auto ts_xxxxyy_zzzzz = pbuffer.data(idx_ovl_ih + 83);

    auto ts_xxxxyz_xxxxx = pbuffer.data(idx_ovl_ih + 84);

    auto ts_xxxxyz_xxxxy = pbuffer.data(idx_ovl_ih + 85);

    auto ts_xxxxyz_xxxxz = pbuffer.data(idx_ovl_ih + 86);

    auto ts_xxxxyz_xxxyy = pbuffer.data(idx_ovl_ih + 87);

    auto ts_xxxxyz_xxxyz = pbuffer.data(idx_ovl_ih + 88);

    auto ts_xxxxyz_xxxzz = pbuffer.data(idx_ovl_ih + 89);

    auto ts_xxxxyz_xxyyy = pbuffer.data(idx_ovl_ih + 90);

    auto ts_xxxxyz_xxyyz = pbuffer.data(idx_ovl_ih + 91);

    auto ts_xxxxyz_xxyzz = pbuffer.data(idx_ovl_ih + 92);

    auto ts_xxxxyz_xxzzz = pbuffer.data(idx_ovl_ih + 93);

    auto ts_xxxxyz_xyyyy = pbuffer.data(idx_ovl_ih + 94);

    auto ts_xxxxyz_xyyyz = pbuffer.data(idx_ovl_ih + 95);

    auto ts_xxxxyz_xyyzz = pbuffer.data(idx_ovl_ih + 96);

    auto ts_xxxxyz_xyzzz = pbuffer.data(idx_ovl_ih + 97);

    auto ts_xxxxyz_xzzzz = pbuffer.data(idx_ovl_ih + 98);

    auto ts_xxxxyz_yyyyy = pbuffer.data(idx_ovl_ih + 99);

    auto ts_xxxxyz_yyyyz = pbuffer.data(idx_ovl_ih + 100);

    auto ts_xxxxyz_yyyzz = pbuffer.data(idx_ovl_ih + 101);

    auto ts_xxxxyz_yyzzz = pbuffer.data(idx_ovl_ih + 102);

    auto ts_xxxxyz_yzzzz = pbuffer.data(idx_ovl_ih + 103);

    auto ts_xxxxyz_zzzzz = pbuffer.data(idx_ovl_ih + 104);

    auto ts_xxxxzz_xxxxx = pbuffer.data(idx_ovl_ih + 105);

    auto ts_xxxxzz_xxxxy = pbuffer.data(idx_ovl_ih + 106);

    auto ts_xxxxzz_xxxxz = pbuffer.data(idx_ovl_ih + 107);

    auto ts_xxxxzz_xxxyy = pbuffer.data(idx_ovl_ih + 108);

    auto ts_xxxxzz_xxxyz = pbuffer.data(idx_ovl_ih + 109);

    auto ts_xxxxzz_xxxzz = pbuffer.data(idx_ovl_ih + 110);

    auto ts_xxxxzz_xxyyy = pbuffer.data(idx_ovl_ih + 111);

    auto ts_xxxxzz_xxyyz = pbuffer.data(idx_ovl_ih + 112);

    auto ts_xxxxzz_xxyzz = pbuffer.data(idx_ovl_ih + 113);

    auto ts_xxxxzz_xxzzz = pbuffer.data(idx_ovl_ih + 114);

    auto ts_xxxxzz_xyyyy = pbuffer.data(idx_ovl_ih + 115);

    auto ts_xxxxzz_xyyyz = pbuffer.data(idx_ovl_ih + 116);

    auto ts_xxxxzz_xyyzz = pbuffer.data(idx_ovl_ih + 117);

    auto ts_xxxxzz_xyzzz = pbuffer.data(idx_ovl_ih + 118);

    auto ts_xxxxzz_xzzzz = pbuffer.data(idx_ovl_ih + 119);

    auto ts_xxxxzz_yyyyy = pbuffer.data(idx_ovl_ih + 120);

    auto ts_xxxxzz_yyyyz = pbuffer.data(idx_ovl_ih + 121);

    auto ts_xxxxzz_yyyzz = pbuffer.data(idx_ovl_ih + 122);

    auto ts_xxxxzz_yyzzz = pbuffer.data(idx_ovl_ih + 123);

    auto ts_xxxxzz_yzzzz = pbuffer.data(idx_ovl_ih + 124);

    auto ts_xxxxzz_zzzzz = pbuffer.data(idx_ovl_ih + 125);

    auto ts_xxxyyy_xxxxx = pbuffer.data(idx_ovl_ih + 126);

    auto ts_xxxyyy_xxxxy = pbuffer.data(idx_ovl_ih + 127);

    auto ts_xxxyyy_xxxxz = pbuffer.data(idx_ovl_ih + 128);

    auto ts_xxxyyy_xxxyy = pbuffer.data(idx_ovl_ih + 129);

    auto ts_xxxyyy_xxxyz = pbuffer.data(idx_ovl_ih + 130);

    auto ts_xxxyyy_xxxzz = pbuffer.data(idx_ovl_ih + 131);

    auto ts_xxxyyy_xxyyy = pbuffer.data(idx_ovl_ih + 132);

    auto ts_xxxyyy_xxyyz = pbuffer.data(idx_ovl_ih + 133);

    auto ts_xxxyyy_xxyzz = pbuffer.data(idx_ovl_ih + 134);

    auto ts_xxxyyy_xxzzz = pbuffer.data(idx_ovl_ih + 135);

    auto ts_xxxyyy_xyyyy = pbuffer.data(idx_ovl_ih + 136);

    auto ts_xxxyyy_xyyyz = pbuffer.data(idx_ovl_ih + 137);

    auto ts_xxxyyy_xyyzz = pbuffer.data(idx_ovl_ih + 138);

    auto ts_xxxyyy_xyzzz = pbuffer.data(idx_ovl_ih + 139);

    auto ts_xxxyyy_xzzzz = pbuffer.data(idx_ovl_ih + 140);

    auto ts_xxxyyy_yyyyy = pbuffer.data(idx_ovl_ih + 141);

    auto ts_xxxyyy_yyyyz = pbuffer.data(idx_ovl_ih + 142);

    auto ts_xxxyyy_yyyzz = pbuffer.data(idx_ovl_ih + 143);

    auto ts_xxxyyy_yyzzz = pbuffer.data(idx_ovl_ih + 144);

    auto ts_xxxyyy_yzzzz = pbuffer.data(idx_ovl_ih + 145);

    auto ts_xxxyyy_zzzzz = pbuffer.data(idx_ovl_ih + 146);

    auto ts_xxxyyz_xxxxx = pbuffer.data(idx_ovl_ih + 147);

    auto ts_xxxyyz_xxxxy = pbuffer.data(idx_ovl_ih + 148);

    auto ts_xxxyyz_xxxxz = pbuffer.data(idx_ovl_ih + 149);

    auto ts_xxxyyz_xxxyy = pbuffer.data(idx_ovl_ih + 150);

    auto ts_xxxyyz_xxxyz = pbuffer.data(idx_ovl_ih + 151);

    auto ts_xxxyyz_xxxzz = pbuffer.data(idx_ovl_ih + 152);

    auto ts_xxxyyz_xxyyy = pbuffer.data(idx_ovl_ih + 153);

    auto ts_xxxyyz_xxyyz = pbuffer.data(idx_ovl_ih + 154);

    auto ts_xxxyyz_xxyzz = pbuffer.data(idx_ovl_ih + 155);

    auto ts_xxxyyz_xxzzz = pbuffer.data(idx_ovl_ih + 156);

    auto ts_xxxyyz_xyyyy = pbuffer.data(idx_ovl_ih + 157);

    auto ts_xxxyyz_xyyyz = pbuffer.data(idx_ovl_ih + 158);

    auto ts_xxxyyz_xyyzz = pbuffer.data(idx_ovl_ih + 159);

    auto ts_xxxyyz_xyzzz = pbuffer.data(idx_ovl_ih + 160);

    auto ts_xxxyyz_xzzzz = pbuffer.data(idx_ovl_ih + 161);

    auto ts_xxxyyz_yyyyy = pbuffer.data(idx_ovl_ih + 162);

    auto ts_xxxyyz_yyyyz = pbuffer.data(idx_ovl_ih + 163);

    auto ts_xxxyyz_yyyzz = pbuffer.data(idx_ovl_ih + 164);

    auto ts_xxxyyz_yyzzz = pbuffer.data(idx_ovl_ih + 165);

    auto ts_xxxyyz_yzzzz = pbuffer.data(idx_ovl_ih + 166);

    auto ts_xxxyyz_zzzzz = pbuffer.data(idx_ovl_ih + 167);

    auto ts_xxxyzz_xxxxx = pbuffer.data(idx_ovl_ih + 168);

    auto ts_xxxyzz_xxxxy = pbuffer.data(idx_ovl_ih + 169);

    auto ts_xxxyzz_xxxxz = pbuffer.data(idx_ovl_ih + 170);

    auto ts_xxxyzz_xxxyy = pbuffer.data(idx_ovl_ih + 171);

    auto ts_xxxyzz_xxxyz = pbuffer.data(idx_ovl_ih + 172);

    auto ts_xxxyzz_xxxzz = pbuffer.data(idx_ovl_ih + 173);

    auto ts_xxxyzz_xxyyy = pbuffer.data(idx_ovl_ih + 174);

    auto ts_xxxyzz_xxyyz = pbuffer.data(idx_ovl_ih + 175);

    auto ts_xxxyzz_xxyzz = pbuffer.data(idx_ovl_ih + 176);

    auto ts_xxxyzz_xxzzz = pbuffer.data(idx_ovl_ih + 177);

    auto ts_xxxyzz_xyyyy = pbuffer.data(idx_ovl_ih + 178);

    auto ts_xxxyzz_xyyyz = pbuffer.data(idx_ovl_ih + 179);

    auto ts_xxxyzz_xyyzz = pbuffer.data(idx_ovl_ih + 180);

    auto ts_xxxyzz_xyzzz = pbuffer.data(idx_ovl_ih + 181);

    auto ts_xxxyzz_xzzzz = pbuffer.data(idx_ovl_ih + 182);

    auto ts_xxxyzz_yyyyy = pbuffer.data(idx_ovl_ih + 183);

    auto ts_xxxyzz_yyyyz = pbuffer.data(idx_ovl_ih + 184);

    auto ts_xxxyzz_yyyzz = pbuffer.data(idx_ovl_ih + 185);

    auto ts_xxxyzz_yyzzz = pbuffer.data(idx_ovl_ih + 186);

    auto ts_xxxyzz_yzzzz = pbuffer.data(idx_ovl_ih + 187);

    auto ts_xxxyzz_zzzzz = pbuffer.data(idx_ovl_ih + 188);

    auto ts_xxxzzz_xxxxx = pbuffer.data(idx_ovl_ih + 189);

    auto ts_xxxzzz_xxxxy = pbuffer.data(idx_ovl_ih + 190);

    auto ts_xxxzzz_xxxxz = pbuffer.data(idx_ovl_ih + 191);

    auto ts_xxxzzz_xxxyy = pbuffer.data(idx_ovl_ih + 192);

    auto ts_xxxzzz_xxxyz = pbuffer.data(idx_ovl_ih + 193);

    auto ts_xxxzzz_xxxzz = pbuffer.data(idx_ovl_ih + 194);

    auto ts_xxxzzz_xxyyy = pbuffer.data(idx_ovl_ih + 195);

    auto ts_xxxzzz_xxyyz = pbuffer.data(idx_ovl_ih + 196);

    auto ts_xxxzzz_xxyzz = pbuffer.data(idx_ovl_ih + 197);

    auto ts_xxxzzz_xxzzz = pbuffer.data(idx_ovl_ih + 198);

    auto ts_xxxzzz_xyyyy = pbuffer.data(idx_ovl_ih + 199);

    auto ts_xxxzzz_xyyyz = pbuffer.data(idx_ovl_ih + 200);

    auto ts_xxxzzz_xyyzz = pbuffer.data(idx_ovl_ih + 201);

    auto ts_xxxzzz_xyzzz = pbuffer.data(idx_ovl_ih + 202);

    auto ts_xxxzzz_xzzzz = pbuffer.data(idx_ovl_ih + 203);

    auto ts_xxxzzz_yyyyy = pbuffer.data(idx_ovl_ih + 204);

    auto ts_xxxzzz_yyyyz = pbuffer.data(idx_ovl_ih + 205);

    auto ts_xxxzzz_yyyzz = pbuffer.data(idx_ovl_ih + 206);

    auto ts_xxxzzz_yyzzz = pbuffer.data(idx_ovl_ih + 207);

    auto ts_xxxzzz_yzzzz = pbuffer.data(idx_ovl_ih + 208);

    auto ts_xxxzzz_zzzzz = pbuffer.data(idx_ovl_ih + 209);

    auto ts_xxyyyy_xxxxx = pbuffer.data(idx_ovl_ih + 210);

    auto ts_xxyyyy_xxxxy = pbuffer.data(idx_ovl_ih + 211);

    auto ts_xxyyyy_xxxxz = pbuffer.data(idx_ovl_ih + 212);

    auto ts_xxyyyy_xxxyy = pbuffer.data(idx_ovl_ih + 213);

    auto ts_xxyyyy_xxxyz = pbuffer.data(idx_ovl_ih + 214);

    auto ts_xxyyyy_xxxzz = pbuffer.data(idx_ovl_ih + 215);

    auto ts_xxyyyy_xxyyy = pbuffer.data(idx_ovl_ih + 216);

    auto ts_xxyyyy_xxyyz = pbuffer.data(idx_ovl_ih + 217);

    auto ts_xxyyyy_xxyzz = pbuffer.data(idx_ovl_ih + 218);

    auto ts_xxyyyy_xxzzz = pbuffer.data(idx_ovl_ih + 219);

    auto ts_xxyyyy_xyyyy = pbuffer.data(idx_ovl_ih + 220);

    auto ts_xxyyyy_xyyyz = pbuffer.data(idx_ovl_ih + 221);

    auto ts_xxyyyy_xyyzz = pbuffer.data(idx_ovl_ih + 222);

    auto ts_xxyyyy_xyzzz = pbuffer.data(idx_ovl_ih + 223);

    auto ts_xxyyyy_xzzzz = pbuffer.data(idx_ovl_ih + 224);

    auto ts_xxyyyy_yyyyy = pbuffer.data(idx_ovl_ih + 225);

    auto ts_xxyyyy_yyyyz = pbuffer.data(idx_ovl_ih + 226);

    auto ts_xxyyyy_yyyzz = pbuffer.data(idx_ovl_ih + 227);

    auto ts_xxyyyy_yyzzz = pbuffer.data(idx_ovl_ih + 228);

    auto ts_xxyyyy_yzzzz = pbuffer.data(idx_ovl_ih + 229);

    auto ts_xxyyyy_zzzzz = pbuffer.data(idx_ovl_ih + 230);

    auto ts_xxyyyz_xxxxx = pbuffer.data(idx_ovl_ih + 231);

    auto ts_xxyyyz_xxxxy = pbuffer.data(idx_ovl_ih + 232);

    auto ts_xxyyyz_xxxxz = pbuffer.data(idx_ovl_ih + 233);

    auto ts_xxyyyz_xxxyy = pbuffer.data(idx_ovl_ih + 234);

    auto ts_xxyyyz_xxxyz = pbuffer.data(idx_ovl_ih + 235);

    auto ts_xxyyyz_xxxzz = pbuffer.data(idx_ovl_ih + 236);

    auto ts_xxyyyz_xxyyy = pbuffer.data(idx_ovl_ih + 237);

    auto ts_xxyyyz_xxyyz = pbuffer.data(idx_ovl_ih + 238);

    auto ts_xxyyyz_xxyzz = pbuffer.data(idx_ovl_ih + 239);

    auto ts_xxyyyz_xxzzz = pbuffer.data(idx_ovl_ih + 240);

    auto ts_xxyyyz_xyyyy = pbuffer.data(idx_ovl_ih + 241);

    auto ts_xxyyyz_xyyyz = pbuffer.data(idx_ovl_ih + 242);

    auto ts_xxyyyz_xyyzz = pbuffer.data(idx_ovl_ih + 243);

    auto ts_xxyyyz_xyzzz = pbuffer.data(idx_ovl_ih + 244);

    auto ts_xxyyyz_xzzzz = pbuffer.data(idx_ovl_ih + 245);

    auto ts_xxyyyz_yyyyy = pbuffer.data(idx_ovl_ih + 246);

    auto ts_xxyyyz_yyyyz = pbuffer.data(idx_ovl_ih + 247);

    auto ts_xxyyyz_yyyzz = pbuffer.data(idx_ovl_ih + 248);

    auto ts_xxyyyz_yyzzz = pbuffer.data(idx_ovl_ih + 249);

    auto ts_xxyyyz_yzzzz = pbuffer.data(idx_ovl_ih + 250);

    auto ts_xxyyyz_zzzzz = pbuffer.data(idx_ovl_ih + 251);

    auto ts_xxyyzz_xxxxx = pbuffer.data(idx_ovl_ih + 252);

    auto ts_xxyyzz_xxxxy = pbuffer.data(idx_ovl_ih + 253);

    auto ts_xxyyzz_xxxxz = pbuffer.data(idx_ovl_ih + 254);

    auto ts_xxyyzz_xxxyy = pbuffer.data(idx_ovl_ih + 255);

    auto ts_xxyyzz_xxxyz = pbuffer.data(idx_ovl_ih + 256);

    auto ts_xxyyzz_xxxzz = pbuffer.data(idx_ovl_ih + 257);

    auto ts_xxyyzz_xxyyy = pbuffer.data(idx_ovl_ih + 258);

    auto ts_xxyyzz_xxyyz = pbuffer.data(idx_ovl_ih + 259);

    auto ts_xxyyzz_xxyzz = pbuffer.data(idx_ovl_ih + 260);

    auto ts_xxyyzz_xxzzz = pbuffer.data(idx_ovl_ih + 261);

    auto ts_xxyyzz_xyyyy = pbuffer.data(idx_ovl_ih + 262);

    auto ts_xxyyzz_xyyyz = pbuffer.data(idx_ovl_ih + 263);

    auto ts_xxyyzz_xyyzz = pbuffer.data(idx_ovl_ih + 264);

    auto ts_xxyyzz_xyzzz = pbuffer.data(idx_ovl_ih + 265);

    auto ts_xxyyzz_xzzzz = pbuffer.data(idx_ovl_ih + 266);

    auto ts_xxyyzz_yyyyy = pbuffer.data(idx_ovl_ih + 267);

    auto ts_xxyyzz_yyyyz = pbuffer.data(idx_ovl_ih + 268);

    auto ts_xxyyzz_yyyzz = pbuffer.data(idx_ovl_ih + 269);

    auto ts_xxyyzz_yyzzz = pbuffer.data(idx_ovl_ih + 270);

    auto ts_xxyyzz_yzzzz = pbuffer.data(idx_ovl_ih + 271);

    auto ts_xxyyzz_zzzzz = pbuffer.data(idx_ovl_ih + 272);

    auto ts_xxyzzz_xxxxx = pbuffer.data(idx_ovl_ih + 273);

    auto ts_xxyzzz_xxxxy = pbuffer.data(idx_ovl_ih + 274);

    auto ts_xxyzzz_xxxxz = pbuffer.data(idx_ovl_ih + 275);

    auto ts_xxyzzz_xxxyy = pbuffer.data(idx_ovl_ih + 276);

    auto ts_xxyzzz_xxxyz = pbuffer.data(idx_ovl_ih + 277);

    auto ts_xxyzzz_xxxzz = pbuffer.data(idx_ovl_ih + 278);

    auto ts_xxyzzz_xxyyy = pbuffer.data(idx_ovl_ih + 279);

    auto ts_xxyzzz_xxyyz = pbuffer.data(idx_ovl_ih + 280);

    auto ts_xxyzzz_xxyzz = pbuffer.data(idx_ovl_ih + 281);

    auto ts_xxyzzz_xxzzz = pbuffer.data(idx_ovl_ih + 282);

    auto ts_xxyzzz_xyyyy = pbuffer.data(idx_ovl_ih + 283);

    auto ts_xxyzzz_xyyyz = pbuffer.data(idx_ovl_ih + 284);

    auto ts_xxyzzz_xyyzz = pbuffer.data(idx_ovl_ih + 285);

    auto ts_xxyzzz_xyzzz = pbuffer.data(idx_ovl_ih + 286);

    auto ts_xxyzzz_xzzzz = pbuffer.data(idx_ovl_ih + 287);

    auto ts_xxyzzz_yyyyy = pbuffer.data(idx_ovl_ih + 288);

    auto ts_xxyzzz_yyyyz = pbuffer.data(idx_ovl_ih + 289);

    auto ts_xxyzzz_yyyzz = pbuffer.data(idx_ovl_ih + 290);

    auto ts_xxyzzz_yyzzz = pbuffer.data(idx_ovl_ih + 291);

    auto ts_xxyzzz_yzzzz = pbuffer.data(idx_ovl_ih + 292);

    auto ts_xxyzzz_zzzzz = pbuffer.data(idx_ovl_ih + 293);

    auto ts_xxzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 294);

    auto ts_xxzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 295);

    auto ts_xxzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 296);

    auto ts_xxzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 297);

    auto ts_xxzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 298);

    auto ts_xxzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 299);

    auto ts_xxzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 300);

    auto ts_xxzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 301);

    auto ts_xxzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 302);

    auto ts_xxzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 303);

    auto ts_xxzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 304);

    auto ts_xxzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 305);

    auto ts_xxzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 306);

    auto ts_xxzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 307);

    auto ts_xxzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 308);

    auto ts_xxzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 309);

    auto ts_xxzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 310);

    auto ts_xxzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 311);

    auto ts_xxzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 312);

    auto ts_xxzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 313);

    auto ts_xxzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 314);

    auto ts_xyyyyy_xxxxx = pbuffer.data(idx_ovl_ih + 315);

    auto ts_xyyyyy_xxxxy = pbuffer.data(idx_ovl_ih + 316);

    auto ts_xyyyyy_xxxxz = pbuffer.data(idx_ovl_ih + 317);

    auto ts_xyyyyy_xxxyy = pbuffer.data(idx_ovl_ih + 318);

    auto ts_xyyyyy_xxxyz = pbuffer.data(idx_ovl_ih + 319);

    auto ts_xyyyyy_xxxzz = pbuffer.data(idx_ovl_ih + 320);

    auto ts_xyyyyy_xxyyy = pbuffer.data(idx_ovl_ih + 321);

    auto ts_xyyyyy_xxyyz = pbuffer.data(idx_ovl_ih + 322);

    auto ts_xyyyyy_xxyzz = pbuffer.data(idx_ovl_ih + 323);

    auto ts_xyyyyy_xxzzz = pbuffer.data(idx_ovl_ih + 324);

    auto ts_xyyyyy_xyyyy = pbuffer.data(idx_ovl_ih + 325);

    auto ts_xyyyyy_xyyyz = pbuffer.data(idx_ovl_ih + 326);

    auto ts_xyyyyy_xyyzz = pbuffer.data(idx_ovl_ih + 327);

    auto ts_xyyyyy_xyzzz = pbuffer.data(idx_ovl_ih + 328);

    auto ts_xyyyyy_xzzzz = pbuffer.data(idx_ovl_ih + 329);

    auto ts_xyyyyy_yyyyy = pbuffer.data(idx_ovl_ih + 330);

    auto ts_xyyyyy_yyyyz = pbuffer.data(idx_ovl_ih + 331);

    auto ts_xyyyyy_yyyzz = pbuffer.data(idx_ovl_ih + 332);

    auto ts_xyyyyy_yyzzz = pbuffer.data(idx_ovl_ih + 333);

    auto ts_xyyyyy_yzzzz = pbuffer.data(idx_ovl_ih + 334);

    auto ts_xyyyyy_zzzzz = pbuffer.data(idx_ovl_ih + 335);

    auto ts_xyyyyz_xxxxx = pbuffer.data(idx_ovl_ih + 336);

    auto ts_xyyyyz_xxxxy = pbuffer.data(idx_ovl_ih + 337);

    auto ts_xyyyyz_xxxxz = pbuffer.data(idx_ovl_ih + 338);

    auto ts_xyyyyz_xxxyy = pbuffer.data(idx_ovl_ih + 339);

    auto ts_xyyyyz_xxxyz = pbuffer.data(idx_ovl_ih + 340);

    auto ts_xyyyyz_xxxzz = pbuffer.data(idx_ovl_ih + 341);

    auto ts_xyyyyz_xxyyy = pbuffer.data(idx_ovl_ih + 342);

    auto ts_xyyyyz_xxyyz = pbuffer.data(idx_ovl_ih + 343);

    auto ts_xyyyyz_xxyzz = pbuffer.data(idx_ovl_ih + 344);

    auto ts_xyyyyz_xxzzz = pbuffer.data(idx_ovl_ih + 345);

    auto ts_xyyyyz_xyyyy = pbuffer.data(idx_ovl_ih + 346);

    auto ts_xyyyyz_xyyyz = pbuffer.data(idx_ovl_ih + 347);

    auto ts_xyyyyz_xyyzz = pbuffer.data(idx_ovl_ih + 348);

    auto ts_xyyyyz_xyzzz = pbuffer.data(idx_ovl_ih + 349);

    auto ts_xyyyyz_xzzzz = pbuffer.data(idx_ovl_ih + 350);

    auto ts_xyyyyz_yyyyy = pbuffer.data(idx_ovl_ih + 351);

    auto ts_xyyyyz_yyyyz = pbuffer.data(idx_ovl_ih + 352);

    auto ts_xyyyyz_yyyzz = pbuffer.data(idx_ovl_ih + 353);

    auto ts_xyyyyz_yyzzz = pbuffer.data(idx_ovl_ih + 354);

    auto ts_xyyyyz_yzzzz = pbuffer.data(idx_ovl_ih + 355);

    auto ts_xyyyyz_zzzzz = pbuffer.data(idx_ovl_ih + 356);

    auto ts_xyyyzz_xxxxx = pbuffer.data(idx_ovl_ih + 357);

    auto ts_xyyyzz_xxxxy = pbuffer.data(idx_ovl_ih + 358);

    auto ts_xyyyzz_xxxxz = pbuffer.data(idx_ovl_ih + 359);

    auto ts_xyyyzz_xxxyy = pbuffer.data(idx_ovl_ih + 360);

    auto ts_xyyyzz_xxxyz = pbuffer.data(idx_ovl_ih + 361);

    auto ts_xyyyzz_xxxzz = pbuffer.data(idx_ovl_ih + 362);

    auto ts_xyyyzz_xxyyy = pbuffer.data(idx_ovl_ih + 363);

    auto ts_xyyyzz_xxyyz = pbuffer.data(idx_ovl_ih + 364);

    auto ts_xyyyzz_xxyzz = pbuffer.data(idx_ovl_ih + 365);

    auto ts_xyyyzz_xxzzz = pbuffer.data(idx_ovl_ih + 366);

    auto ts_xyyyzz_xyyyy = pbuffer.data(idx_ovl_ih + 367);

    auto ts_xyyyzz_xyyyz = pbuffer.data(idx_ovl_ih + 368);

    auto ts_xyyyzz_xyyzz = pbuffer.data(idx_ovl_ih + 369);

    auto ts_xyyyzz_xyzzz = pbuffer.data(idx_ovl_ih + 370);

    auto ts_xyyyzz_xzzzz = pbuffer.data(idx_ovl_ih + 371);

    auto ts_xyyyzz_yyyyy = pbuffer.data(idx_ovl_ih + 372);

    auto ts_xyyyzz_yyyyz = pbuffer.data(idx_ovl_ih + 373);

    auto ts_xyyyzz_yyyzz = pbuffer.data(idx_ovl_ih + 374);

    auto ts_xyyyzz_yyzzz = pbuffer.data(idx_ovl_ih + 375);

    auto ts_xyyyzz_yzzzz = pbuffer.data(idx_ovl_ih + 376);

    auto ts_xyyyzz_zzzzz = pbuffer.data(idx_ovl_ih + 377);

    auto ts_xyyzzz_xxxxx = pbuffer.data(idx_ovl_ih + 378);

    auto ts_xyyzzz_xxxxy = pbuffer.data(idx_ovl_ih + 379);

    auto ts_xyyzzz_xxxxz = pbuffer.data(idx_ovl_ih + 380);

    auto ts_xyyzzz_xxxyy = pbuffer.data(idx_ovl_ih + 381);

    auto ts_xyyzzz_xxxyz = pbuffer.data(idx_ovl_ih + 382);

    auto ts_xyyzzz_xxxzz = pbuffer.data(idx_ovl_ih + 383);

    auto ts_xyyzzz_xxyyy = pbuffer.data(idx_ovl_ih + 384);

    auto ts_xyyzzz_xxyyz = pbuffer.data(idx_ovl_ih + 385);

    auto ts_xyyzzz_xxyzz = pbuffer.data(idx_ovl_ih + 386);

    auto ts_xyyzzz_xxzzz = pbuffer.data(idx_ovl_ih + 387);

    auto ts_xyyzzz_xyyyy = pbuffer.data(idx_ovl_ih + 388);

    auto ts_xyyzzz_xyyyz = pbuffer.data(idx_ovl_ih + 389);

    auto ts_xyyzzz_xyyzz = pbuffer.data(idx_ovl_ih + 390);

    auto ts_xyyzzz_xyzzz = pbuffer.data(idx_ovl_ih + 391);

    auto ts_xyyzzz_xzzzz = pbuffer.data(idx_ovl_ih + 392);

    auto ts_xyyzzz_yyyyy = pbuffer.data(idx_ovl_ih + 393);

    auto ts_xyyzzz_yyyyz = pbuffer.data(idx_ovl_ih + 394);

    auto ts_xyyzzz_yyyzz = pbuffer.data(idx_ovl_ih + 395);

    auto ts_xyyzzz_yyzzz = pbuffer.data(idx_ovl_ih + 396);

    auto ts_xyyzzz_yzzzz = pbuffer.data(idx_ovl_ih + 397);

    auto ts_xyyzzz_zzzzz = pbuffer.data(idx_ovl_ih + 398);

    auto ts_xyzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 399);

    auto ts_xyzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 400);

    auto ts_xyzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 401);

    auto ts_xyzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 402);

    auto ts_xyzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 403);

    auto ts_xyzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 404);

    auto ts_xyzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 405);

    auto ts_xyzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 406);

    auto ts_xyzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 407);

    auto ts_xyzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 408);

    auto ts_xyzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 409);

    auto ts_xyzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 410);

    auto ts_xyzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 411);

    auto ts_xyzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 412);

    auto ts_xyzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 413);

    auto ts_xyzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 414);

    auto ts_xyzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 415);

    auto ts_xyzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 416);

    auto ts_xyzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 417);

    auto ts_xyzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 418);

    auto ts_xyzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 419);

    auto ts_xzzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 420);

    auto ts_xzzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 421);

    auto ts_xzzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 422);

    auto ts_xzzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 423);

    auto ts_xzzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 424);

    auto ts_xzzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 425);

    auto ts_xzzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 426);

    auto ts_xzzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 427);

    auto ts_xzzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 428);

    auto ts_xzzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 429);

    auto ts_xzzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 430);

    auto ts_xzzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 431);

    auto ts_xzzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 432);

    auto ts_xzzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 433);

    auto ts_xzzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 434);

    auto ts_xzzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 435);

    auto ts_xzzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 436);

    auto ts_xzzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 437);

    auto ts_xzzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 438);

    auto ts_xzzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 439);

    auto ts_xzzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 440);

    auto ts_yyyyyy_xxxxx = pbuffer.data(idx_ovl_ih + 441);

    auto ts_yyyyyy_xxxxy = pbuffer.data(idx_ovl_ih + 442);

    auto ts_yyyyyy_xxxxz = pbuffer.data(idx_ovl_ih + 443);

    auto ts_yyyyyy_xxxyy = pbuffer.data(idx_ovl_ih + 444);

    auto ts_yyyyyy_xxxyz = pbuffer.data(idx_ovl_ih + 445);

    auto ts_yyyyyy_xxxzz = pbuffer.data(idx_ovl_ih + 446);

    auto ts_yyyyyy_xxyyy = pbuffer.data(idx_ovl_ih + 447);

    auto ts_yyyyyy_xxyyz = pbuffer.data(idx_ovl_ih + 448);

    auto ts_yyyyyy_xxyzz = pbuffer.data(idx_ovl_ih + 449);

    auto ts_yyyyyy_xxzzz = pbuffer.data(idx_ovl_ih + 450);

    auto ts_yyyyyy_xyyyy = pbuffer.data(idx_ovl_ih + 451);

    auto ts_yyyyyy_xyyyz = pbuffer.data(idx_ovl_ih + 452);

    auto ts_yyyyyy_xyyzz = pbuffer.data(idx_ovl_ih + 453);

    auto ts_yyyyyy_xyzzz = pbuffer.data(idx_ovl_ih + 454);

    auto ts_yyyyyy_xzzzz = pbuffer.data(idx_ovl_ih + 455);

    auto ts_yyyyyy_yyyyy = pbuffer.data(idx_ovl_ih + 456);

    auto ts_yyyyyy_yyyyz = pbuffer.data(idx_ovl_ih + 457);

    auto ts_yyyyyy_yyyzz = pbuffer.data(idx_ovl_ih + 458);

    auto ts_yyyyyy_yyzzz = pbuffer.data(idx_ovl_ih + 459);

    auto ts_yyyyyy_yzzzz = pbuffer.data(idx_ovl_ih + 460);

    auto ts_yyyyyy_zzzzz = pbuffer.data(idx_ovl_ih + 461);

    auto ts_yyyyyz_xxxxx = pbuffer.data(idx_ovl_ih + 462);

    auto ts_yyyyyz_xxxxy = pbuffer.data(idx_ovl_ih + 463);

    auto ts_yyyyyz_xxxxz = pbuffer.data(idx_ovl_ih + 464);

    auto ts_yyyyyz_xxxyy = pbuffer.data(idx_ovl_ih + 465);

    auto ts_yyyyyz_xxxyz = pbuffer.data(idx_ovl_ih + 466);

    auto ts_yyyyyz_xxxzz = pbuffer.data(idx_ovl_ih + 467);

    auto ts_yyyyyz_xxyyy = pbuffer.data(idx_ovl_ih + 468);

    auto ts_yyyyyz_xxyyz = pbuffer.data(idx_ovl_ih + 469);

    auto ts_yyyyyz_xxyzz = pbuffer.data(idx_ovl_ih + 470);

    auto ts_yyyyyz_xxzzz = pbuffer.data(idx_ovl_ih + 471);

    auto ts_yyyyyz_xyyyy = pbuffer.data(idx_ovl_ih + 472);

    auto ts_yyyyyz_xyyyz = pbuffer.data(idx_ovl_ih + 473);

    auto ts_yyyyyz_xyyzz = pbuffer.data(idx_ovl_ih + 474);

    auto ts_yyyyyz_xyzzz = pbuffer.data(idx_ovl_ih + 475);

    auto ts_yyyyyz_xzzzz = pbuffer.data(idx_ovl_ih + 476);

    auto ts_yyyyyz_yyyyy = pbuffer.data(idx_ovl_ih + 477);

    auto ts_yyyyyz_yyyyz = pbuffer.data(idx_ovl_ih + 478);

    auto ts_yyyyyz_yyyzz = pbuffer.data(idx_ovl_ih + 479);

    auto ts_yyyyyz_yyzzz = pbuffer.data(idx_ovl_ih + 480);

    auto ts_yyyyyz_yzzzz = pbuffer.data(idx_ovl_ih + 481);

    auto ts_yyyyyz_zzzzz = pbuffer.data(idx_ovl_ih + 482);

    auto ts_yyyyzz_xxxxx = pbuffer.data(idx_ovl_ih + 483);

    auto ts_yyyyzz_xxxxy = pbuffer.data(idx_ovl_ih + 484);

    auto ts_yyyyzz_xxxxz = pbuffer.data(idx_ovl_ih + 485);

    auto ts_yyyyzz_xxxyy = pbuffer.data(idx_ovl_ih + 486);

    auto ts_yyyyzz_xxxyz = pbuffer.data(idx_ovl_ih + 487);

    auto ts_yyyyzz_xxxzz = pbuffer.data(idx_ovl_ih + 488);

    auto ts_yyyyzz_xxyyy = pbuffer.data(idx_ovl_ih + 489);

    auto ts_yyyyzz_xxyyz = pbuffer.data(idx_ovl_ih + 490);

    auto ts_yyyyzz_xxyzz = pbuffer.data(idx_ovl_ih + 491);

    auto ts_yyyyzz_xxzzz = pbuffer.data(idx_ovl_ih + 492);

    auto ts_yyyyzz_xyyyy = pbuffer.data(idx_ovl_ih + 493);

    auto ts_yyyyzz_xyyyz = pbuffer.data(idx_ovl_ih + 494);

    auto ts_yyyyzz_xyyzz = pbuffer.data(idx_ovl_ih + 495);

    auto ts_yyyyzz_xyzzz = pbuffer.data(idx_ovl_ih + 496);

    auto ts_yyyyzz_xzzzz = pbuffer.data(idx_ovl_ih + 497);

    auto ts_yyyyzz_yyyyy = pbuffer.data(idx_ovl_ih + 498);

    auto ts_yyyyzz_yyyyz = pbuffer.data(idx_ovl_ih + 499);

    auto ts_yyyyzz_yyyzz = pbuffer.data(idx_ovl_ih + 500);

    auto ts_yyyyzz_yyzzz = pbuffer.data(idx_ovl_ih + 501);

    auto ts_yyyyzz_yzzzz = pbuffer.data(idx_ovl_ih + 502);

    auto ts_yyyyzz_zzzzz = pbuffer.data(idx_ovl_ih + 503);

    auto ts_yyyzzz_xxxxx = pbuffer.data(idx_ovl_ih + 504);

    auto ts_yyyzzz_xxxxy = pbuffer.data(idx_ovl_ih + 505);

    auto ts_yyyzzz_xxxxz = pbuffer.data(idx_ovl_ih + 506);

    auto ts_yyyzzz_xxxyy = pbuffer.data(idx_ovl_ih + 507);

    auto ts_yyyzzz_xxxyz = pbuffer.data(idx_ovl_ih + 508);

    auto ts_yyyzzz_xxxzz = pbuffer.data(idx_ovl_ih + 509);

    auto ts_yyyzzz_xxyyy = pbuffer.data(idx_ovl_ih + 510);

    auto ts_yyyzzz_xxyyz = pbuffer.data(idx_ovl_ih + 511);

    auto ts_yyyzzz_xxyzz = pbuffer.data(idx_ovl_ih + 512);

    auto ts_yyyzzz_xxzzz = pbuffer.data(idx_ovl_ih + 513);

    auto ts_yyyzzz_xyyyy = pbuffer.data(idx_ovl_ih + 514);

    auto ts_yyyzzz_xyyyz = pbuffer.data(idx_ovl_ih + 515);

    auto ts_yyyzzz_xyyzz = pbuffer.data(idx_ovl_ih + 516);

    auto ts_yyyzzz_xyzzz = pbuffer.data(idx_ovl_ih + 517);

    auto ts_yyyzzz_xzzzz = pbuffer.data(idx_ovl_ih + 518);

    auto ts_yyyzzz_yyyyy = pbuffer.data(idx_ovl_ih + 519);

    auto ts_yyyzzz_yyyyz = pbuffer.data(idx_ovl_ih + 520);

    auto ts_yyyzzz_yyyzz = pbuffer.data(idx_ovl_ih + 521);

    auto ts_yyyzzz_yyzzz = pbuffer.data(idx_ovl_ih + 522);

    auto ts_yyyzzz_yzzzz = pbuffer.data(idx_ovl_ih + 523);

    auto ts_yyyzzz_zzzzz = pbuffer.data(idx_ovl_ih + 524);

    auto ts_yyzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 525);

    auto ts_yyzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 526);

    auto ts_yyzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 527);

    auto ts_yyzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 528);

    auto ts_yyzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 529);

    auto ts_yyzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 530);

    auto ts_yyzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 531);

    auto ts_yyzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 532);

    auto ts_yyzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 533);

    auto ts_yyzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 534);

    auto ts_yyzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 535);

    auto ts_yyzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 536);

    auto ts_yyzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 537);

    auto ts_yyzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 538);

    auto ts_yyzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 539);

    auto ts_yyzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 540);

    auto ts_yyzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 541);

    auto ts_yyzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 542);

    auto ts_yyzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 543);

    auto ts_yyzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 544);

    auto ts_yyzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 545);

    auto ts_yzzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 546);

    auto ts_yzzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 547);

    auto ts_yzzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 548);

    auto ts_yzzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 549);

    auto ts_yzzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 550);

    auto ts_yzzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 551);

    auto ts_yzzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 552);

    auto ts_yzzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 553);

    auto ts_yzzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 554);

    auto ts_yzzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 555);

    auto ts_yzzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 556);

    auto ts_yzzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 557);

    auto ts_yzzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 558);

    auto ts_yzzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 559);

    auto ts_yzzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 560);

    auto ts_yzzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 561);

    auto ts_yzzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 562);

    auto ts_yzzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 563);

    auto ts_yzzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 564);

    auto ts_yzzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 565);

    auto ts_yzzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 566);

    auto ts_zzzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 567);

    auto ts_zzzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 568);

    auto ts_zzzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 569);

    auto ts_zzzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 570);

    auto ts_zzzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 571);

    auto ts_zzzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 572);

    auto ts_zzzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 573);

    auto ts_zzzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 574);

    auto ts_zzzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 575);

    auto ts_zzzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 576);

    auto ts_zzzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 577);

    auto ts_zzzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 578);

    auto ts_zzzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 579);

    auto ts_zzzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 580);

    auto ts_zzzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 581);

    auto ts_zzzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 582);

    auto ts_zzzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 583);

    auto ts_zzzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 584);

    auto ts_zzzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 585);

    auto ts_zzzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 586);

    auto ts_zzzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 587);

    // Set up 0-21 components of targeted buffer : IH

    auto tk_xxxxxx_xxxxx = pbuffer.data(idx_kin_ih);

    auto tk_xxxxxx_xxxxy = pbuffer.data(idx_kin_ih + 1);

    auto tk_xxxxxx_xxxxz = pbuffer.data(idx_kin_ih + 2);

    auto tk_xxxxxx_xxxyy = pbuffer.data(idx_kin_ih + 3);

    auto tk_xxxxxx_xxxyz = pbuffer.data(idx_kin_ih + 4);

    auto tk_xxxxxx_xxxzz = pbuffer.data(idx_kin_ih + 5);

    auto tk_xxxxxx_xxyyy = pbuffer.data(idx_kin_ih + 6);

    auto tk_xxxxxx_xxyyz = pbuffer.data(idx_kin_ih + 7);

    auto tk_xxxxxx_xxyzz = pbuffer.data(idx_kin_ih + 8);

    auto tk_xxxxxx_xxzzz = pbuffer.data(idx_kin_ih + 9);

    auto tk_xxxxxx_xyyyy = pbuffer.data(idx_kin_ih + 10);

    auto tk_xxxxxx_xyyyz = pbuffer.data(idx_kin_ih + 11);

    auto tk_xxxxxx_xyyzz = pbuffer.data(idx_kin_ih + 12);

    auto tk_xxxxxx_xyzzz = pbuffer.data(idx_kin_ih + 13);

    auto tk_xxxxxx_xzzzz = pbuffer.data(idx_kin_ih + 14);

    auto tk_xxxxxx_yyyyy = pbuffer.data(idx_kin_ih + 15);

    auto tk_xxxxxx_yyyyz = pbuffer.data(idx_kin_ih + 16);

    auto tk_xxxxxx_yyyzz = pbuffer.data(idx_kin_ih + 17);

    auto tk_xxxxxx_yyzzz = pbuffer.data(idx_kin_ih + 18);

    auto tk_xxxxxx_yzzzz = pbuffer.data(idx_kin_ih + 19);

    auto tk_xxxxxx_zzzzz = pbuffer.data(idx_kin_ih + 20);

    #pragma omp simd aligned(pa_x, tk_xxxx_xxxxx, tk_xxxx_xxxxy, tk_xxxx_xxxxz, tk_xxxx_xxxyy, tk_xxxx_xxxyz, tk_xxxx_xxxzz, tk_xxxx_xxyyy, tk_xxxx_xxyyz, tk_xxxx_xxyzz, tk_xxxx_xxzzz, tk_xxxx_xyyyy, tk_xxxx_xyyyz, tk_xxxx_xyyzz, tk_xxxx_xyzzz, tk_xxxx_xzzzz, tk_xxxx_yyyyy, tk_xxxx_yyyyz, tk_xxxx_yyyzz, tk_xxxx_yyzzz, tk_xxxx_yzzzz, tk_xxxx_zzzzz, tk_xxxxx_xxxx, tk_xxxxx_xxxxx, tk_xxxxx_xxxxy, tk_xxxxx_xxxxz, tk_xxxxx_xxxy, tk_xxxxx_xxxyy, tk_xxxxx_xxxyz, tk_xxxxx_xxxz, tk_xxxxx_xxxzz, tk_xxxxx_xxyy, tk_xxxxx_xxyyy, tk_xxxxx_xxyyz, tk_xxxxx_xxyz, tk_xxxxx_xxyzz, tk_xxxxx_xxzz, tk_xxxxx_xxzzz, tk_xxxxx_xyyy, tk_xxxxx_xyyyy, tk_xxxxx_xyyyz, tk_xxxxx_xyyz, tk_xxxxx_xyyzz, tk_xxxxx_xyzz, tk_xxxxx_xyzzz, tk_xxxxx_xzzz, tk_xxxxx_xzzzz, tk_xxxxx_yyyy, tk_xxxxx_yyyyy, tk_xxxxx_yyyyz, tk_xxxxx_yyyz, tk_xxxxx_yyyzz, tk_xxxxx_yyzz, tk_xxxxx_yyzzz, tk_xxxxx_yzzz, tk_xxxxx_yzzzz, tk_xxxxx_zzzz, tk_xxxxx_zzzzz, tk_xxxxxx_xxxxx, tk_xxxxxx_xxxxy, tk_xxxxxx_xxxxz, tk_xxxxxx_xxxyy, tk_xxxxxx_xxxyz, tk_xxxxxx_xxxzz, tk_xxxxxx_xxyyy, tk_xxxxxx_xxyyz, tk_xxxxxx_xxyzz, tk_xxxxxx_xxzzz, tk_xxxxxx_xyyyy, tk_xxxxxx_xyyyz, tk_xxxxxx_xyyzz, tk_xxxxxx_xyzzz, tk_xxxxxx_xzzzz, tk_xxxxxx_yyyyy, tk_xxxxxx_yyyyz, tk_xxxxxx_yyyzz, tk_xxxxxx_yyzzz, tk_xxxxxx_yzzzz, tk_xxxxxx_zzzzz, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxzz, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyzz, ts_xxxx_xxzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyzz, ts_xxxx_xyzzz, ts_xxxx_xzzzz, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyzz, ts_xxxx_yyzzz, ts_xxxx_yzzzz, ts_xxxx_zzzzz, ts_xxxxxx_xxxxx, ts_xxxxxx_xxxxy, ts_xxxxxx_xxxxz, ts_xxxxxx_xxxyy, ts_xxxxxx_xxxyz, ts_xxxxxx_xxxzz, ts_xxxxxx_xxyyy, ts_xxxxxx_xxyyz, ts_xxxxxx_xxyzz, ts_xxxxxx_xxzzz, ts_xxxxxx_xyyyy, ts_xxxxxx_xyyyz, ts_xxxxxx_xyyzz, ts_xxxxxx_xyzzz, ts_xxxxxx_xzzzz, ts_xxxxxx_yyyyy, ts_xxxxxx_yyyyz, ts_xxxxxx_yyyzz, ts_xxxxxx_yyzzz, ts_xxxxxx_yzzzz, ts_xxxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxxx_xxxxx[i] = -10.0 * ts_xxxx_xxxxx[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxx[i] * fe_0 + 5.0 * tk_xxxxx_xxxx[i] * fe_0 + tk_xxxxx_xxxxx[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxx[i] * fz_0;

        tk_xxxxxx_xxxxy[i] = -10.0 * ts_xxxx_xxxxy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxy[i] * fe_0 + 4.0 * tk_xxxxx_xxxy[i] * fe_0 + tk_xxxxx_xxxxy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxy[i] * fz_0;

        tk_xxxxxx_xxxxz[i] = -10.0 * ts_xxxx_xxxxz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxxz[i] * fe_0 + 4.0 * tk_xxxxx_xxxz[i] * fe_0 + tk_xxxxx_xxxxz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxxz[i] * fz_0;

        tk_xxxxxx_xxxyy[i] = -10.0 * ts_xxxx_xxxyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxyy[i] * fe_0 + 3.0 * tk_xxxxx_xxyy[i] * fe_0 + tk_xxxxx_xxxyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxyy[i] * fz_0;

        tk_xxxxxx_xxxyz[i] = -10.0 * ts_xxxx_xxxyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxyz[i] * fe_0 + 3.0 * tk_xxxxx_xxyz[i] * fe_0 + tk_xxxxx_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxyz[i] * fz_0;

        tk_xxxxxx_xxxzz[i] = -10.0 * ts_xxxx_xxxzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxxzz[i] * fe_0 + 3.0 * tk_xxxxx_xxzz[i] * fe_0 + tk_xxxxx_xxxzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxxzz[i] * fz_0;

        tk_xxxxxx_xxyyy[i] = -10.0 * ts_xxxx_xxyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyyy[i] * fe_0 + 2.0 * tk_xxxxx_xyyy[i] * fe_0 + tk_xxxxx_xxyyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyyy[i] * fz_0;

        tk_xxxxxx_xxyyz[i] = -10.0 * ts_xxxx_xxyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyyz[i] * fe_0 + 2.0 * tk_xxxxx_xyyz[i] * fe_0 + tk_xxxxx_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyyz[i] * fz_0;

        tk_xxxxxx_xxyzz[i] = -10.0 * ts_xxxx_xxyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxyzz[i] * fe_0 + 2.0 * tk_xxxxx_xyzz[i] * fe_0 + tk_xxxxx_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxyzz[i] * fz_0;

        tk_xxxxxx_xxzzz[i] = -10.0 * ts_xxxx_xxzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xxzzz[i] * fe_0 + 2.0 * tk_xxxxx_xzzz[i] * fe_0 + tk_xxxxx_xxzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xxzzz[i] * fz_0;

        tk_xxxxxx_xyyyy[i] = -10.0 * ts_xxxx_xyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyyy[i] * fe_0 + tk_xxxxx_yyyy[i] * fe_0 + tk_xxxxx_xyyyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyyy[i] * fz_0;

        tk_xxxxxx_xyyyz[i] = -10.0 * ts_xxxx_xyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyyz[i] * fe_0 + tk_xxxxx_yyyz[i] * fe_0 + tk_xxxxx_xyyyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyyz[i] * fz_0;

        tk_xxxxxx_xyyzz[i] = -10.0 * ts_xxxx_xyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyyzz[i] * fe_0 + tk_xxxxx_yyzz[i] * fe_0 + tk_xxxxx_xyyzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyyzz[i] * fz_0;

        tk_xxxxxx_xyzzz[i] = -10.0 * ts_xxxx_xyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xyzzz[i] * fe_0 + tk_xxxxx_yzzz[i] * fe_0 + tk_xxxxx_xyzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xyzzz[i] * fz_0;

        tk_xxxxxx_xzzzz[i] = -10.0 * ts_xxxx_xzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xzzzz[i] * fe_0 + tk_xxxxx_zzzz[i] * fe_0 + tk_xxxxx_xzzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_xzzzz[i] * fz_0;

        tk_xxxxxx_yyyyy[i] = -10.0 * ts_xxxx_yyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyyy[i] * fe_0 + tk_xxxxx_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyyyy[i] * fz_0;

        tk_xxxxxx_yyyyz[i] = -10.0 * ts_xxxx_yyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyyz[i] * fe_0 + tk_xxxxx_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyyyz[i] * fz_0;

        tk_xxxxxx_yyyzz[i] = -10.0 * ts_xxxx_yyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyyzz[i] * fe_0 + tk_xxxxx_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyyzz[i] * fz_0;

        tk_xxxxxx_yyzzz[i] = -10.0 * ts_xxxx_yyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yyzzz[i] * fe_0 + tk_xxxxx_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yyzzz[i] * fz_0;

        tk_xxxxxx_yzzzz[i] = -10.0 * ts_xxxx_yzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yzzzz[i] * fe_0 + tk_xxxxx_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yzzzz[i] * fz_0;

        tk_xxxxxx_zzzzz[i] = -10.0 * ts_xxxx_zzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_zzzzz[i] * fe_0 + tk_xxxxx_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxxxx_zzzzz[i] * fz_0;
    }

    // Set up 21-42 components of targeted buffer : IH

    auto tk_xxxxxy_xxxxx = pbuffer.data(idx_kin_ih + 21);

    auto tk_xxxxxy_xxxxy = pbuffer.data(idx_kin_ih + 22);

    auto tk_xxxxxy_xxxxz = pbuffer.data(idx_kin_ih + 23);

    auto tk_xxxxxy_xxxyy = pbuffer.data(idx_kin_ih + 24);

    auto tk_xxxxxy_xxxyz = pbuffer.data(idx_kin_ih + 25);

    auto tk_xxxxxy_xxxzz = pbuffer.data(idx_kin_ih + 26);

    auto tk_xxxxxy_xxyyy = pbuffer.data(idx_kin_ih + 27);

    auto tk_xxxxxy_xxyyz = pbuffer.data(idx_kin_ih + 28);

    auto tk_xxxxxy_xxyzz = pbuffer.data(idx_kin_ih + 29);

    auto tk_xxxxxy_xxzzz = pbuffer.data(idx_kin_ih + 30);

    auto tk_xxxxxy_xyyyy = pbuffer.data(idx_kin_ih + 31);

    auto tk_xxxxxy_xyyyz = pbuffer.data(idx_kin_ih + 32);

    auto tk_xxxxxy_xyyzz = pbuffer.data(idx_kin_ih + 33);

    auto tk_xxxxxy_xyzzz = pbuffer.data(idx_kin_ih + 34);

    auto tk_xxxxxy_xzzzz = pbuffer.data(idx_kin_ih + 35);

    auto tk_xxxxxy_yyyyy = pbuffer.data(idx_kin_ih + 36);

    auto tk_xxxxxy_yyyyz = pbuffer.data(idx_kin_ih + 37);

    auto tk_xxxxxy_yyyzz = pbuffer.data(idx_kin_ih + 38);

    auto tk_xxxxxy_yyzzz = pbuffer.data(idx_kin_ih + 39);

    auto tk_xxxxxy_yzzzz = pbuffer.data(idx_kin_ih + 40);

    auto tk_xxxxxy_zzzzz = pbuffer.data(idx_kin_ih + 41);

    #pragma omp simd aligned(pa_y, tk_xxxxx_xxxx, tk_xxxxx_xxxxx, tk_xxxxx_xxxxy, tk_xxxxx_xxxxz, tk_xxxxx_xxxy, tk_xxxxx_xxxyy, tk_xxxxx_xxxyz, tk_xxxxx_xxxz, tk_xxxxx_xxxzz, tk_xxxxx_xxyy, tk_xxxxx_xxyyy, tk_xxxxx_xxyyz, tk_xxxxx_xxyz, tk_xxxxx_xxyzz, tk_xxxxx_xxzz, tk_xxxxx_xxzzz, tk_xxxxx_xyyy, tk_xxxxx_xyyyy, tk_xxxxx_xyyyz, tk_xxxxx_xyyz, tk_xxxxx_xyyzz, tk_xxxxx_xyzz, tk_xxxxx_xyzzz, tk_xxxxx_xzzz, tk_xxxxx_xzzzz, tk_xxxxx_yyyy, tk_xxxxx_yyyyy, tk_xxxxx_yyyyz, tk_xxxxx_yyyz, tk_xxxxx_yyyzz, tk_xxxxx_yyzz, tk_xxxxx_yyzzz, tk_xxxxx_yzzz, tk_xxxxx_yzzzz, tk_xxxxx_zzzz, tk_xxxxx_zzzzz, tk_xxxxxy_xxxxx, tk_xxxxxy_xxxxy, tk_xxxxxy_xxxxz, tk_xxxxxy_xxxyy, tk_xxxxxy_xxxyz, tk_xxxxxy_xxxzz, tk_xxxxxy_xxyyy, tk_xxxxxy_xxyyz, tk_xxxxxy_xxyzz, tk_xxxxxy_xxzzz, tk_xxxxxy_xyyyy, tk_xxxxxy_xyyyz, tk_xxxxxy_xyyzz, tk_xxxxxy_xyzzz, tk_xxxxxy_xzzzz, tk_xxxxxy_yyyyy, tk_xxxxxy_yyyyz, tk_xxxxxy_yyyzz, tk_xxxxxy_yyzzz, tk_xxxxxy_yzzzz, tk_xxxxxy_zzzzz, ts_xxxxxy_xxxxx, ts_xxxxxy_xxxxy, ts_xxxxxy_xxxxz, ts_xxxxxy_xxxyy, ts_xxxxxy_xxxyz, ts_xxxxxy_xxxzz, ts_xxxxxy_xxyyy, ts_xxxxxy_xxyyz, ts_xxxxxy_xxyzz, ts_xxxxxy_xxzzz, ts_xxxxxy_xyyyy, ts_xxxxxy_xyyyz, ts_xxxxxy_xyyzz, ts_xxxxxy_xyzzz, ts_xxxxxy_xzzzz, ts_xxxxxy_yyyyy, ts_xxxxxy_yyyyz, ts_xxxxxy_yyyzz, ts_xxxxxy_yyzzz, ts_xxxxxy_yzzzz, ts_xxxxxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxy_xxxxx[i] = tk_xxxxx_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxx[i] * fz_0;

        tk_xxxxxy_xxxxy[i] = tk_xxxxx_xxxx[i] * fe_0 + tk_xxxxx_xxxxy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxy[i] * fz_0;

        tk_xxxxxy_xxxxz[i] = tk_xxxxx_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxxz[i] * fz_0;

        tk_xxxxxy_xxxyy[i] = 2.0 * tk_xxxxx_xxxy[i] * fe_0 + tk_xxxxx_xxxyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxyy[i] * fz_0;

        tk_xxxxxy_xxxyz[i] = tk_xxxxx_xxxz[i] * fe_0 + tk_xxxxx_xxxyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxyz[i] * fz_0;

        tk_xxxxxy_xxxzz[i] = tk_xxxxx_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxxzz[i] * fz_0;

        tk_xxxxxy_xxyyy[i] = 3.0 * tk_xxxxx_xxyy[i] * fe_0 + tk_xxxxx_xxyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyyy[i] * fz_0;

        tk_xxxxxy_xxyyz[i] = 2.0 * tk_xxxxx_xxyz[i] * fe_0 + tk_xxxxx_xxyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyyz[i] * fz_0;

        tk_xxxxxy_xxyzz[i] = tk_xxxxx_xxzz[i] * fe_0 + tk_xxxxx_xxyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxyzz[i] * fz_0;

        tk_xxxxxy_xxzzz[i] = tk_xxxxx_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xxzzz[i] * fz_0;

        tk_xxxxxy_xyyyy[i] = 4.0 * tk_xxxxx_xyyy[i] * fe_0 + tk_xxxxx_xyyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyyy[i] * fz_0;

        tk_xxxxxy_xyyyz[i] = 3.0 * tk_xxxxx_xyyz[i] * fe_0 + tk_xxxxx_xyyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyyz[i] * fz_0;

        tk_xxxxxy_xyyzz[i] = 2.0 * tk_xxxxx_xyzz[i] * fe_0 + tk_xxxxx_xyyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyyzz[i] * fz_0;

        tk_xxxxxy_xyzzz[i] = tk_xxxxx_xzzz[i] * fe_0 + tk_xxxxx_xyzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xyzzz[i] * fz_0;

        tk_xxxxxy_xzzzz[i] = tk_xxxxx_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xzzzz[i] * fz_0;

        tk_xxxxxy_yyyyy[i] = 5.0 * tk_xxxxx_yyyy[i] * fe_0 + tk_xxxxx_yyyyy[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyyy[i] * fz_0;

        tk_xxxxxy_yyyyz[i] = 4.0 * tk_xxxxx_yyyz[i] * fe_0 + tk_xxxxx_yyyyz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyyz[i] * fz_0;

        tk_xxxxxy_yyyzz[i] = 3.0 * tk_xxxxx_yyzz[i] * fe_0 + tk_xxxxx_yyyzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyyzz[i] * fz_0;

        tk_xxxxxy_yyzzz[i] = 2.0 * tk_xxxxx_yzzz[i] * fe_0 + tk_xxxxx_yyzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yyzzz[i] * fz_0;

        tk_xxxxxy_yzzzz[i] = tk_xxxxx_zzzz[i] * fe_0 + tk_xxxxx_yzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yzzzz[i] * fz_0;

        tk_xxxxxy_zzzzz[i] = tk_xxxxx_zzzzz[i] * pa_y[i] + 2.0 * ts_xxxxxy_zzzzz[i] * fz_0;
    }

    // Set up 42-63 components of targeted buffer : IH

    auto tk_xxxxxz_xxxxx = pbuffer.data(idx_kin_ih + 42);

    auto tk_xxxxxz_xxxxy = pbuffer.data(idx_kin_ih + 43);

    auto tk_xxxxxz_xxxxz = pbuffer.data(idx_kin_ih + 44);

    auto tk_xxxxxz_xxxyy = pbuffer.data(idx_kin_ih + 45);

    auto tk_xxxxxz_xxxyz = pbuffer.data(idx_kin_ih + 46);

    auto tk_xxxxxz_xxxzz = pbuffer.data(idx_kin_ih + 47);

    auto tk_xxxxxz_xxyyy = pbuffer.data(idx_kin_ih + 48);

    auto tk_xxxxxz_xxyyz = pbuffer.data(idx_kin_ih + 49);

    auto tk_xxxxxz_xxyzz = pbuffer.data(idx_kin_ih + 50);

    auto tk_xxxxxz_xxzzz = pbuffer.data(idx_kin_ih + 51);

    auto tk_xxxxxz_xyyyy = pbuffer.data(idx_kin_ih + 52);

    auto tk_xxxxxz_xyyyz = pbuffer.data(idx_kin_ih + 53);

    auto tk_xxxxxz_xyyzz = pbuffer.data(idx_kin_ih + 54);

    auto tk_xxxxxz_xyzzz = pbuffer.data(idx_kin_ih + 55);

    auto tk_xxxxxz_xzzzz = pbuffer.data(idx_kin_ih + 56);

    auto tk_xxxxxz_yyyyy = pbuffer.data(idx_kin_ih + 57);

    auto tk_xxxxxz_yyyyz = pbuffer.data(idx_kin_ih + 58);

    auto tk_xxxxxz_yyyzz = pbuffer.data(idx_kin_ih + 59);

    auto tk_xxxxxz_yyzzz = pbuffer.data(idx_kin_ih + 60);

    auto tk_xxxxxz_yzzzz = pbuffer.data(idx_kin_ih + 61);

    auto tk_xxxxxz_zzzzz = pbuffer.data(idx_kin_ih + 62);

    #pragma omp simd aligned(pa_z, tk_xxxxx_xxxx, tk_xxxxx_xxxxx, tk_xxxxx_xxxxy, tk_xxxxx_xxxxz, tk_xxxxx_xxxy, tk_xxxxx_xxxyy, tk_xxxxx_xxxyz, tk_xxxxx_xxxz, tk_xxxxx_xxxzz, tk_xxxxx_xxyy, tk_xxxxx_xxyyy, tk_xxxxx_xxyyz, tk_xxxxx_xxyz, tk_xxxxx_xxyzz, tk_xxxxx_xxzz, tk_xxxxx_xxzzz, tk_xxxxx_xyyy, tk_xxxxx_xyyyy, tk_xxxxx_xyyyz, tk_xxxxx_xyyz, tk_xxxxx_xyyzz, tk_xxxxx_xyzz, tk_xxxxx_xyzzz, tk_xxxxx_xzzz, tk_xxxxx_xzzzz, tk_xxxxx_yyyy, tk_xxxxx_yyyyy, tk_xxxxx_yyyyz, tk_xxxxx_yyyz, tk_xxxxx_yyyzz, tk_xxxxx_yyzz, tk_xxxxx_yyzzz, tk_xxxxx_yzzz, tk_xxxxx_yzzzz, tk_xxxxx_zzzz, tk_xxxxx_zzzzz, tk_xxxxxz_xxxxx, tk_xxxxxz_xxxxy, tk_xxxxxz_xxxxz, tk_xxxxxz_xxxyy, tk_xxxxxz_xxxyz, tk_xxxxxz_xxxzz, tk_xxxxxz_xxyyy, tk_xxxxxz_xxyyz, tk_xxxxxz_xxyzz, tk_xxxxxz_xxzzz, tk_xxxxxz_xyyyy, tk_xxxxxz_xyyyz, tk_xxxxxz_xyyzz, tk_xxxxxz_xyzzz, tk_xxxxxz_xzzzz, tk_xxxxxz_yyyyy, tk_xxxxxz_yyyyz, tk_xxxxxz_yyyzz, tk_xxxxxz_yyzzz, tk_xxxxxz_yzzzz, tk_xxxxxz_zzzzz, ts_xxxxxz_xxxxx, ts_xxxxxz_xxxxy, ts_xxxxxz_xxxxz, ts_xxxxxz_xxxyy, ts_xxxxxz_xxxyz, ts_xxxxxz_xxxzz, ts_xxxxxz_xxyyy, ts_xxxxxz_xxyyz, ts_xxxxxz_xxyzz, ts_xxxxxz_xxzzz, ts_xxxxxz_xyyyy, ts_xxxxxz_xyyyz, ts_xxxxxz_xyyzz, ts_xxxxxz_xyzzz, ts_xxxxxz_xzzzz, ts_xxxxxz_yyyyy, ts_xxxxxz_yyyyz, ts_xxxxxz_yyyzz, ts_xxxxxz_yyzzz, ts_xxxxxz_yzzzz, ts_xxxxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxz_xxxxx[i] = tk_xxxxx_xxxxx[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxx[i] * fz_0;

        tk_xxxxxz_xxxxy[i] = tk_xxxxx_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxy[i] * fz_0;

        tk_xxxxxz_xxxxz[i] = tk_xxxxx_xxxx[i] * fe_0 + tk_xxxxx_xxxxz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxxz[i] * fz_0;

        tk_xxxxxz_xxxyy[i] = tk_xxxxx_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxyy[i] * fz_0;

        tk_xxxxxz_xxxyz[i] = tk_xxxxx_xxxy[i] * fe_0 + tk_xxxxx_xxxyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxyz[i] * fz_0;

        tk_xxxxxz_xxxzz[i] = 2.0 * tk_xxxxx_xxxz[i] * fe_0 + tk_xxxxx_xxxzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxxzz[i] * fz_0;

        tk_xxxxxz_xxyyy[i] = tk_xxxxx_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyyy[i] * fz_0;

        tk_xxxxxz_xxyyz[i] = tk_xxxxx_xxyy[i] * fe_0 + tk_xxxxx_xxyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyyz[i] * fz_0;

        tk_xxxxxz_xxyzz[i] = 2.0 * tk_xxxxx_xxyz[i] * fe_0 + tk_xxxxx_xxyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxyzz[i] * fz_0;

        tk_xxxxxz_xxzzz[i] = 3.0 * tk_xxxxx_xxzz[i] * fe_0 + tk_xxxxx_xxzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xxzzz[i] * fz_0;

        tk_xxxxxz_xyyyy[i] = tk_xxxxx_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyyy[i] * fz_0;

        tk_xxxxxz_xyyyz[i] = tk_xxxxx_xyyy[i] * fe_0 + tk_xxxxx_xyyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyyz[i] * fz_0;

        tk_xxxxxz_xyyzz[i] = 2.0 * tk_xxxxx_xyyz[i] * fe_0 + tk_xxxxx_xyyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyyzz[i] * fz_0;

        tk_xxxxxz_xyzzz[i] = 3.0 * tk_xxxxx_xyzz[i] * fe_0 + tk_xxxxx_xyzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xyzzz[i] * fz_0;

        tk_xxxxxz_xzzzz[i] = 4.0 * tk_xxxxx_xzzz[i] * fe_0 + tk_xxxxx_xzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xzzzz[i] * fz_0;

        tk_xxxxxz_yyyyy[i] = tk_xxxxx_yyyyy[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyyy[i] * fz_0;

        tk_xxxxxz_yyyyz[i] = tk_xxxxx_yyyy[i] * fe_0 + tk_xxxxx_yyyyz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyyz[i] * fz_0;

        tk_xxxxxz_yyyzz[i] = 2.0 * tk_xxxxx_yyyz[i] * fe_0 + tk_xxxxx_yyyzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyyzz[i] * fz_0;

        tk_xxxxxz_yyzzz[i] = 3.0 * tk_xxxxx_yyzz[i] * fe_0 + tk_xxxxx_yyzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yyzzz[i] * fz_0;

        tk_xxxxxz_yzzzz[i] = 4.0 * tk_xxxxx_yzzz[i] * fe_0 + tk_xxxxx_yzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yzzzz[i] * fz_0;

        tk_xxxxxz_zzzzz[i] = 5.0 * tk_xxxxx_zzzz[i] * fe_0 + tk_xxxxx_zzzzz[i] * pa_z[i] + 2.0 * ts_xxxxxz_zzzzz[i] * fz_0;
    }

    // Set up 63-84 components of targeted buffer : IH

    auto tk_xxxxyy_xxxxx = pbuffer.data(idx_kin_ih + 63);

    auto tk_xxxxyy_xxxxy = pbuffer.data(idx_kin_ih + 64);

    auto tk_xxxxyy_xxxxz = pbuffer.data(idx_kin_ih + 65);

    auto tk_xxxxyy_xxxyy = pbuffer.data(idx_kin_ih + 66);

    auto tk_xxxxyy_xxxyz = pbuffer.data(idx_kin_ih + 67);

    auto tk_xxxxyy_xxxzz = pbuffer.data(idx_kin_ih + 68);

    auto tk_xxxxyy_xxyyy = pbuffer.data(idx_kin_ih + 69);

    auto tk_xxxxyy_xxyyz = pbuffer.data(idx_kin_ih + 70);

    auto tk_xxxxyy_xxyzz = pbuffer.data(idx_kin_ih + 71);

    auto tk_xxxxyy_xxzzz = pbuffer.data(idx_kin_ih + 72);

    auto tk_xxxxyy_xyyyy = pbuffer.data(idx_kin_ih + 73);

    auto tk_xxxxyy_xyyyz = pbuffer.data(idx_kin_ih + 74);

    auto tk_xxxxyy_xyyzz = pbuffer.data(idx_kin_ih + 75);

    auto tk_xxxxyy_xyzzz = pbuffer.data(idx_kin_ih + 76);

    auto tk_xxxxyy_xzzzz = pbuffer.data(idx_kin_ih + 77);

    auto tk_xxxxyy_yyyyy = pbuffer.data(idx_kin_ih + 78);

    auto tk_xxxxyy_yyyyz = pbuffer.data(idx_kin_ih + 79);

    auto tk_xxxxyy_yyyzz = pbuffer.data(idx_kin_ih + 80);

    auto tk_xxxxyy_yyzzz = pbuffer.data(idx_kin_ih + 81);

    auto tk_xxxxyy_yzzzz = pbuffer.data(idx_kin_ih + 82);

    auto tk_xxxxyy_zzzzz = pbuffer.data(idx_kin_ih + 83);

    #pragma omp simd aligned(pa_x, pa_y, tk_xxxx_xxxxx, tk_xxxx_xxxxz, tk_xxxx_xxxzz, tk_xxxx_xxzzz, tk_xxxx_xzzzz, tk_xxxxy_xxxxx, tk_xxxxy_xxxxz, tk_xxxxy_xxxzz, tk_xxxxy_xxzzz, tk_xxxxy_xzzzz, tk_xxxxyy_xxxxx, tk_xxxxyy_xxxxy, tk_xxxxyy_xxxxz, tk_xxxxyy_xxxyy, tk_xxxxyy_xxxyz, tk_xxxxyy_xxxzz, tk_xxxxyy_xxyyy, tk_xxxxyy_xxyyz, tk_xxxxyy_xxyzz, tk_xxxxyy_xxzzz, tk_xxxxyy_xyyyy, tk_xxxxyy_xyyyz, tk_xxxxyy_xyyzz, tk_xxxxyy_xyzzz, tk_xxxxyy_xzzzz, tk_xxxxyy_yyyyy, tk_xxxxyy_yyyyz, tk_xxxxyy_yyyzz, tk_xxxxyy_yyzzz, tk_xxxxyy_yzzzz, tk_xxxxyy_zzzzz, tk_xxxyy_xxxxy, tk_xxxyy_xxxy, tk_xxxyy_xxxyy, tk_xxxyy_xxxyz, tk_xxxyy_xxyy, tk_xxxyy_xxyyy, tk_xxxyy_xxyyz, tk_xxxyy_xxyz, tk_xxxyy_xxyzz, tk_xxxyy_xyyy, tk_xxxyy_xyyyy, tk_xxxyy_xyyyz, tk_xxxyy_xyyz, tk_xxxyy_xyyzz, tk_xxxyy_xyzz, tk_xxxyy_xyzzz, tk_xxxyy_yyyy, tk_xxxyy_yyyyy, tk_xxxyy_yyyyz, tk_xxxyy_yyyz, tk_xxxyy_yyyzz, tk_xxxyy_yyzz, tk_xxxyy_yyzzz, tk_xxxyy_yzzz, tk_xxxyy_yzzzz, tk_xxxyy_zzzzz, tk_xxyy_xxxxy, tk_xxyy_xxxyy, tk_xxyy_xxxyz, tk_xxyy_xxyyy, tk_xxyy_xxyyz, tk_xxyy_xxyzz, tk_xxyy_xyyyy, tk_xxyy_xyyyz, tk_xxyy_xyyzz, tk_xxyy_xyzzz, tk_xxyy_yyyyy, tk_xxyy_yyyyz, tk_xxyy_yyyzz, tk_xxyy_yyzzz, tk_xxyy_yzzzz, tk_xxyy_zzzzz, ts_xxxx_xxxxx, ts_xxxx_xxxxz, ts_xxxx_xxxzz, ts_xxxx_xxzzz, ts_xxxx_xzzzz, ts_xxxxyy_xxxxx, ts_xxxxyy_xxxxy, ts_xxxxyy_xxxxz, ts_xxxxyy_xxxyy, ts_xxxxyy_xxxyz, ts_xxxxyy_xxxzz, ts_xxxxyy_xxyyy, ts_xxxxyy_xxyyz, ts_xxxxyy_xxyzz, ts_xxxxyy_xxzzz, ts_xxxxyy_xyyyy, ts_xxxxyy_xyyyz, ts_xxxxyy_xyyzz, ts_xxxxyy_xyzzz, ts_xxxxyy_xzzzz, ts_xxxxyy_yyyyy, ts_xxxxyy_yyyyz, ts_xxxxyy_yyyzz, ts_xxxxyy_yyzzz, ts_xxxxyy_yzzzz, ts_xxxxyy_zzzzz, ts_xxyy_xxxxy, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyzz, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyzz, ts_xxyy_xyzzz, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyzz, ts_xxyy_yyzzz, ts_xxyy_yzzzz, ts_xxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxyy_xxxxx[i] = -2.0 * ts_xxxx_xxxxx[i] * fbe_0 * fz_0 + tk_xxxx_xxxxx[i] * fe_0 + tk_xxxxy_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxxx[i] * fz_0;

        tk_xxxxyy_xxxxy[i] = -6.0 * ts_xxyy_xxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxy[i] * fe_0 + 4.0 * tk_xxxyy_xxxy[i] * fe_0 + tk_xxxyy_xxxxy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxxy[i] * fz_0;

        tk_xxxxyy_xxxxz[i] = -2.0 * ts_xxxx_xxxxz[i] * fbe_0 * fz_0 + tk_xxxx_xxxxz[i] * fe_0 + tk_xxxxy_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxxz[i] * fz_0;

        tk_xxxxyy_xxxyy[i] = -6.0 * ts_xxyy_xxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxyy[i] * fe_0 + 3.0 * tk_xxxyy_xxyy[i] * fe_0 + tk_xxxyy_xxxyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxyy[i] * fz_0;

        tk_xxxxyy_xxxyz[i] = -6.0 * ts_xxyy_xxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxyz[i] * fe_0 + 3.0 * tk_xxxyy_xxyz[i] * fe_0 + tk_xxxyy_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxxyz[i] * fz_0;

        tk_xxxxyy_xxxzz[i] = -2.0 * ts_xxxx_xxxzz[i] * fbe_0 * fz_0 + tk_xxxx_xxxzz[i] * fe_0 + tk_xxxxy_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxxzz[i] * fz_0;

        tk_xxxxyy_xxyyy[i] = -6.0 * ts_xxyy_xxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyyy[i] * fe_0 + 2.0 * tk_xxxyy_xyyy[i] * fe_0 + tk_xxxyy_xxyyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyyy[i] * fz_0;

        tk_xxxxyy_xxyyz[i] = -6.0 * ts_xxyy_xxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyyz[i] * fe_0 + 2.0 * tk_xxxyy_xyyz[i] * fe_0 + tk_xxxyy_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyyz[i] * fz_0;

        tk_xxxxyy_xxyzz[i] = -6.0 * ts_xxyy_xxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxyzz[i] * fe_0 + 2.0 * tk_xxxyy_xyzz[i] * fe_0 + tk_xxxyy_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xxyzz[i] * fz_0;

        tk_xxxxyy_xxzzz[i] = -2.0 * ts_xxxx_xxzzz[i] * fbe_0 * fz_0 + tk_xxxx_xxzzz[i] * fe_0 + tk_xxxxy_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xxzzz[i] * fz_0;

        tk_xxxxyy_xyyyy[i] = -6.0 * ts_xxyy_xyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyyy[i] * fe_0 + tk_xxxyy_yyyy[i] * fe_0 + tk_xxxyy_xyyyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyyy[i] * fz_0;

        tk_xxxxyy_xyyyz[i] = -6.0 * ts_xxyy_xyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyyz[i] * fe_0 + tk_xxxyy_yyyz[i] * fe_0 + tk_xxxyy_xyyyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyyz[i] * fz_0;

        tk_xxxxyy_xyyzz[i] = -6.0 * ts_xxyy_xyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyyzz[i] * fe_0 + tk_xxxyy_yyzz[i] * fe_0 + tk_xxxyy_xyyzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyyzz[i] * fz_0;

        tk_xxxxyy_xyzzz[i] = -6.0 * ts_xxyy_xyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xyzzz[i] * fe_0 + tk_xxxyy_yzzz[i] * fe_0 + tk_xxxyy_xyzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_xyzzz[i] * fz_0;

        tk_xxxxyy_xzzzz[i] = -2.0 * ts_xxxx_xzzzz[i] * fbe_0 * fz_0 + tk_xxxx_xzzzz[i] * fe_0 + tk_xxxxy_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xzzzz[i] * fz_0;

        tk_xxxxyy_yyyyy[i] = -6.0 * ts_xxyy_yyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyyy[i] * fe_0 + tk_xxxyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyyyy[i] * fz_0;

        tk_xxxxyy_yyyyz[i] = -6.0 * ts_xxyy_yyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyyz[i] * fe_0 + tk_xxxyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyyyz[i] * fz_0;

        tk_xxxxyy_yyyzz[i] = -6.0 * ts_xxyy_yyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyyzz[i] * fe_0 + tk_xxxyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyyzz[i] * fz_0;

        tk_xxxxyy_yyzzz[i] = -6.0 * ts_xxyy_yyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yyzzz[i] * fe_0 + tk_xxxyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yyzzz[i] * fz_0;

        tk_xxxxyy_yzzzz[i] = -6.0 * ts_xxyy_yzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yzzzz[i] * fe_0 + tk_xxxyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yzzzz[i] * fz_0;

        tk_xxxxyy_zzzzz[i] = -6.0 * ts_xxyy_zzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_zzzzz[i] * fe_0 + tk_xxxyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxxyy_zzzzz[i] * fz_0;
    }

    // Set up 84-105 components of targeted buffer : IH

    auto tk_xxxxyz_xxxxx = pbuffer.data(idx_kin_ih + 84);

    auto tk_xxxxyz_xxxxy = pbuffer.data(idx_kin_ih + 85);

    auto tk_xxxxyz_xxxxz = pbuffer.data(idx_kin_ih + 86);

    auto tk_xxxxyz_xxxyy = pbuffer.data(idx_kin_ih + 87);

    auto tk_xxxxyz_xxxyz = pbuffer.data(idx_kin_ih + 88);

    auto tk_xxxxyz_xxxzz = pbuffer.data(idx_kin_ih + 89);

    auto tk_xxxxyz_xxyyy = pbuffer.data(idx_kin_ih + 90);

    auto tk_xxxxyz_xxyyz = pbuffer.data(idx_kin_ih + 91);

    auto tk_xxxxyz_xxyzz = pbuffer.data(idx_kin_ih + 92);

    auto tk_xxxxyz_xxzzz = pbuffer.data(idx_kin_ih + 93);

    auto tk_xxxxyz_xyyyy = pbuffer.data(idx_kin_ih + 94);

    auto tk_xxxxyz_xyyyz = pbuffer.data(idx_kin_ih + 95);

    auto tk_xxxxyz_xyyzz = pbuffer.data(idx_kin_ih + 96);

    auto tk_xxxxyz_xyzzz = pbuffer.data(idx_kin_ih + 97);

    auto tk_xxxxyz_xzzzz = pbuffer.data(idx_kin_ih + 98);

    auto tk_xxxxyz_yyyyy = pbuffer.data(idx_kin_ih + 99);

    auto tk_xxxxyz_yyyyz = pbuffer.data(idx_kin_ih + 100);

    auto tk_xxxxyz_yyyzz = pbuffer.data(idx_kin_ih + 101);

    auto tk_xxxxyz_yyzzz = pbuffer.data(idx_kin_ih + 102);

    auto tk_xxxxyz_yzzzz = pbuffer.data(idx_kin_ih + 103);

    auto tk_xxxxyz_zzzzz = pbuffer.data(idx_kin_ih + 104);

    #pragma omp simd aligned(pa_y, pa_z, tk_xxxxy_xxxxy, tk_xxxxy_xxxyy, tk_xxxxy_xxyyy, tk_xxxxy_xyyyy, tk_xxxxy_yyyyy, tk_xxxxyz_xxxxx, tk_xxxxyz_xxxxy, tk_xxxxyz_xxxxz, tk_xxxxyz_xxxyy, tk_xxxxyz_xxxyz, tk_xxxxyz_xxxzz, tk_xxxxyz_xxyyy, tk_xxxxyz_xxyyz, tk_xxxxyz_xxyzz, tk_xxxxyz_xxzzz, tk_xxxxyz_xyyyy, tk_xxxxyz_xyyyz, tk_xxxxyz_xyyzz, tk_xxxxyz_xyzzz, tk_xxxxyz_xzzzz, tk_xxxxyz_yyyyy, tk_xxxxyz_yyyyz, tk_xxxxyz_yyyzz, tk_xxxxyz_yyzzz, tk_xxxxyz_yzzzz, tk_xxxxyz_zzzzz, tk_xxxxz_xxxxx, tk_xxxxz_xxxxz, tk_xxxxz_xxxyz, tk_xxxxz_xxxz, tk_xxxxz_xxxzz, tk_xxxxz_xxyyz, tk_xxxxz_xxyz, tk_xxxxz_xxyzz, tk_xxxxz_xxzz, tk_xxxxz_xxzzz, tk_xxxxz_xyyyz, tk_xxxxz_xyyz, tk_xxxxz_xyyzz, tk_xxxxz_xyzz, tk_xxxxz_xyzzz, tk_xxxxz_xzzz, tk_xxxxz_xzzzz, tk_xxxxz_yyyyz, tk_xxxxz_yyyz, tk_xxxxz_yyyzz, tk_xxxxz_yyzz, tk_xxxxz_yyzzz, tk_xxxxz_yzzz, tk_xxxxz_yzzzz, tk_xxxxz_zzzz, tk_xxxxz_zzzzz, ts_xxxxyz_xxxxx, ts_xxxxyz_xxxxy, ts_xxxxyz_xxxxz, ts_xxxxyz_xxxyy, ts_xxxxyz_xxxyz, ts_xxxxyz_xxxzz, ts_xxxxyz_xxyyy, ts_xxxxyz_xxyyz, ts_xxxxyz_xxyzz, ts_xxxxyz_xxzzz, ts_xxxxyz_xyyyy, ts_xxxxyz_xyyyz, ts_xxxxyz_xyyzz, ts_xxxxyz_xyzzz, ts_xxxxyz_xzzzz, ts_xxxxyz_yyyyy, ts_xxxxyz_yyyyz, ts_xxxxyz_yyyzz, ts_xxxxyz_yyzzz, ts_xxxxyz_yzzzz, ts_xxxxyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxyz_xxxxx[i] = tk_xxxxz_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxxx[i] * fz_0;

        tk_xxxxyz_xxxxy[i] = tk_xxxxy_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxxxy[i] * fz_0;

        tk_xxxxyz_xxxxz[i] = tk_xxxxz_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxxz[i] * fz_0;

        tk_xxxxyz_xxxyy[i] = tk_xxxxy_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxxyy[i] * fz_0;

        tk_xxxxyz_xxxyz[i] = tk_xxxxz_xxxz[i] * fe_0 + tk_xxxxz_xxxyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxyz[i] * fz_0;

        tk_xxxxyz_xxxzz[i] = tk_xxxxz_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxxzz[i] * fz_0;

        tk_xxxxyz_xxyyy[i] = tk_xxxxy_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xxyyy[i] * fz_0;

        tk_xxxxyz_xxyyz[i] = 2.0 * tk_xxxxz_xxyz[i] * fe_0 + tk_xxxxz_xxyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxyyz[i] * fz_0;

        tk_xxxxyz_xxyzz[i] = tk_xxxxz_xxzz[i] * fe_0 + tk_xxxxz_xxyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxyzz[i] * fz_0;

        tk_xxxxyz_xxzzz[i] = tk_xxxxz_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xxzzz[i] * fz_0;

        tk_xxxxyz_xyyyy[i] = tk_xxxxy_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xyyyy[i] * fz_0;

        tk_xxxxyz_xyyyz[i] = 3.0 * tk_xxxxz_xyyz[i] * fe_0 + tk_xxxxz_xyyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyyyz[i] * fz_0;

        tk_xxxxyz_xyyzz[i] = 2.0 * tk_xxxxz_xyzz[i] * fe_0 + tk_xxxxz_xyyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyyzz[i] * fz_0;

        tk_xxxxyz_xyzzz[i] = tk_xxxxz_xzzz[i] * fe_0 + tk_xxxxz_xyzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xyzzz[i] * fz_0;

        tk_xxxxyz_xzzzz[i] = tk_xxxxz_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xzzzz[i] * fz_0;

        tk_xxxxyz_yyyyy[i] = tk_xxxxy_yyyyy[i] * pa_z[i] + 2.0 * ts_xxxxyz_yyyyy[i] * fz_0;

        tk_xxxxyz_yyyyz[i] = 4.0 * tk_xxxxz_yyyz[i] * fe_0 + tk_xxxxz_yyyyz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyyyz[i] * fz_0;

        tk_xxxxyz_yyyzz[i] = 3.0 * tk_xxxxz_yyzz[i] * fe_0 + tk_xxxxz_yyyzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyyzz[i] * fz_0;

        tk_xxxxyz_yyzzz[i] = 2.0 * tk_xxxxz_yzzz[i] * fe_0 + tk_xxxxz_yyzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yyzzz[i] * fz_0;

        tk_xxxxyz_yzzzz[i] = tk_xxxxz_zzzz[i] * fe_0 + tk_xxxxz_yzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yzzzz[i] * fz_0;

        tk_xxxxyz_zzzzz[i] = tk_xxxxz_zzzzz[i] * pa_y[i] + 2.0 * ts_xxxxyz_zzzzz[i] * fz_0;
    }

    // Set up 105-126 components of targeted buffer : IH

    auto tk_xxxxzz_xxxxx = pbuffer.data(idx_kin_ih + 105);

    auto tk_xxxxzz_xxxxy = pbuffer.data(idx_kin_ih + 106);

    auto tk_xxxxzz_xxxxz = pbuffer.data(idx_kin_ih + 107);

    auto tk_xxxxzz_xxxyy = pbuffer.data(idx_kin_ih + 108);

    auto tk_xxxxzz_xxxyz = pbuffer.data(idx_kin_ih + 109);

    auto tk_xxxxzz_xxxzz = pbuffer.data(idx_kin_ih + 110);

    auto tk_xxxxzz_xxyyy = pbuffer.data(idx_kin_ih + 111);

    auto tk_xxxxzz_xxyyz = pbuffer.data(idx_kin_ih + 112);

    auto tk_xxxxzz_xxyzz = pbuffer.data(idx_kin_ih + 113);

    auto tk_xxxxzz_xxzzz = pbuffer.data(idx_kin_ih + 114);

    auto tk_xxxxzz_xyyyy = pbuffer.data(idx_kin_ih + 115);

    auto tk_xxxxzz_xyyyz = pbuffer.data(idx_kin_ih + 116);

    auto tk_xxxxzz_xyyzz = pbuffer.data(idx_kin_ih + 117);

    auto tk_xxxxzz_xyzzz = pbuffer.data(idx_kin_ih + 118);

    auto tk_xxxxzz_xzzzz = pbuffer.data(idx_kin_ih + 119);

    auto tk_xxxxzz_yyyyy = pbuffer.data(idx_kin_ih + 120);

    auto tk_xxxxzz_yyyyz = pbuffer.data(idx_kin_ih + 121);

    auto tk_xxxxzz_yyyzz = pbuffer.data(idx_kin_ih + 122);

    auto tk_xxxxzz_yyzzz = pbuffer.data(idx_kin_ih + 123);

    auto tk_xxxxzz_yzzzz = pbuffer.data(idx_kin_ih + 124);

    auto tk_xxxxzz_zzzzz = pbuffer.data(idx_kin_ih + 125);

    #pragma omp simd aligned(pa_x, pa_z, tk_xxxx_xxxxx, tk_xxxx_xxxxy, tk_xxxx_xxxyy, tk_xxxx_xxyyy, tk_xxxx_xyyyy, tk_xxxxz_xxxxx, tk_xxxxz_xxxxy, tk_xxxxz_xxxyy, tk_xxxxz_xxyyy, tk_xxxxz_xyyyy, tk_xxxxzz_xxxxx, tk_xxxxzz_xxxxy, tk_xxxxzz_xxxxz, tk_xxxxzz_xxxyy, tk_xxxxzz_xxxyz, tk_xxxxzz_xxxzz, tk_xxxxzz_xxyyy, tk_xxxxzz_xxyyz, tk_xxxxzz_xxyzz, tk_xxxxzz_xxzzz, tk_xxxxzz_xyyyy, tk_xxxxzz_xyyyz, tk_xxxxzz_xyyzz, tk_xxxxzz_xyzzz, tk_xxxxzz_xzzzz, tk_xxxxzz_yyyyy, tk_xxxxzz_yyyyz, tk_xxxxzz_yyyzz, tk_xxxxzz_yyzzz, tk_xxxxzz_yzzzz, tk_xxxxzz_zzzzz, tk_xxxzz_xxxxz, tk_xxxzz_xxxyz, tk_xxxzz_xxxz, tk_xxxzz_xxxzz, tk_xxxzz_xxyyz, tk_xxxzz_xxyz, tk_xxxzz_xxyzz, tk_xxxzz_xxzz, tk_xxxzz_xxzzz, tk_xxxzz_xyyyz, tk_xxxzz_xyyz, tk_xxxzz_xyyzz, tk_xxxzz_xyzz, tk_xxxzz_xyzzz, tk_xxxzz_xzzz, tk_xxxzz_xzzzz, tk_xxxzz_yyyyy, tk_xxxzz_yyyyz, tk_xxxzz_yyyz, tk_xxxzz_yyyzz, tk_xxxzz_yyzz, tk_xxxzz_yyzzz, tk_xxxzz_yzzz, tk_xxxzz_yzzzz, tk_xxxzz_zzzz, tk_xxxzz_zzzzz, tk_xxzz_xxxxz, tk_xxzz_xxxyz, tk_xxzz_xxxzz, tk_xxzz_xxyyz, tk_xxzz_xxyzz, tk_xxzz_xxzzz, tk_xxzz_xyyyz, tk_xxzz_xyyzz, tk_xxzz_xyzzz, tk_xxzz_xzzzz, tk_xxzz_yyyyy, tk_xxzz_yyyyz, tk_xxzz_yyyzz, tk_xxzz_yyzzz, tk_xxzz_yzzzz, tk_xxzz_zzzzz, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxyy, ts_xxxx_xxyyy, ts_xxxx_xyyyy, ts_xxxxzz_xxxxx, ts_xxxxzz_xxxxy, ts_xxxxzz_xxxxz, ts_xxxxzz_xxxyy, ts_xxxxzz_xxxyz, ts_xxxxzz_xxxzz, ts_xxxxzz_xxyyy, ts_xxxxzz_xxyyz, ts_xxxxzz_xxyzz, ts_xxxxzz_xxzzz, ts_xxxxzz_xyyyy, ts_xxxxzz_xyyyz, ts_xxxxzz_xyyzz, ts_xxxxzz_xyzzz, ts_xxxxzz_xzzzz, ts_xxxxzz_yyyyy, ts_xxxxzz_yyyyz, ts_xxxxzz_yyyzz, ts_xxxxzz_yyzzz, ts_xxxxzz_yzzzz, ts_xxxxzz_zzzzz, ts_xxzz_xxxxz, ts_xxzz_xxxyz, ts_xxzz_xxxzz, ts_xxzz_xxyyz, ts_xxzz_xxyzz, ts_xxzz_xxzzz, ts_xxzz_xyyyz, ts_xxzz_xyyzz, ts_xxzz_xyzzz, ts_xxzz_xzzzz, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyzz, ts_xxzz_yyzzz, ts_xxzz_yzzzz, ts_xxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxzz_xxxxx[i] = -2.0 * ts_xxxx_xxxxx[i] * fbe_0 * fz_0 + tk_xxxx_xxxxx[i] * fe_0 + tk_xxxxz_xxxxx[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxxx[i] * fz_0;

        tk_xxxxzz_xxxxy[i] = -2.0 * ts_xxxx_xxxxy[i] * fbe_0 * fz_0 + tk_xxxx_xxxxy[i] * fe_0 + tk_xxxxz_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxxy[i] * fz_0;

        tk_xxxxzz_xxxxz[i] = -6.0 * ts_xxzz_xxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxz[i] * fe_0 + 4.0 * tk_xxxzz_xxxz[i] * fe_0 + tk_xxxzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxxz[i] * fz_0;

        tk_xxxxzz_xxxyy[i] = -2.0 * ts_xxxx_xxxyy[i] * fbe_0 * fz_0 + tk_xxxx_xxxyy[i] * fe_0 + tk_xxxxz_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxxyy[i] * fz_0;

        tk_xxxxzz_xxxyz[i] = -6.0 * ts_xxzz_xxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxyz[i] * fe_0 + 3.0 * tk_xxxzz_xxyz[i] * fe_0 + tk_xxxzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxyz[i] * fz_0;

        tk_xxxxzz_xxxzz[i] = -6.0 * ts_xxzz_xxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxzz[i] * fe_0 + 3.0 * tk_xxxzz_xxzz[i] * fe_0 + tk_xxxzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxxzz[i] * fz_0;

        tk_xxxxzz_xxyyy[i] = -2.0 * ts_xxxx_xxyyy[i] * fbe_0 * fz_0 + tk_xxxx_xxyyy[i] * fe_0 + tk_xxxxz_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xxyyy[i] * fz_0;

        tk_xxxxzz_xxyyz[i] = -6.0 * ts_xxzz_xxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyyz[i] * fe_0 + 2.0 * tk_xxxzz_xyyz[i] * fe_0 + tk_xxxzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxyyz[i] * fz_0;

        tk_xxxxzz_xxyzz[i] = -6.0 * ts_xxzz_xxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyzz[i] * fe_0 + 2.0 * tk_xxxzz_xyzz[i] * fe_0 + tk_xxxzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxyzz[i] * fz_0;

        tk_xxxxzz_xxzzz[i] = -6.0 * ts_xxzz_xxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxzzz[i] * fe_0 + 2.0 * tk_xxxzz_xzzz[i] * fe_0 + tk_xxxzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xxzzz[i] * fz_0;

        tk_xxxxzz_xyyyy[i] = -2.0 * ts_xxxx_xyyyy[i] * fbe_0 * fz_0 + tk_xxxx_xyyyy[i] * fe_0 + tk_xxxxz_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xyyyy[i] * fz_0;

        tk_xxxxzz_xyyyz[i] = -6.0 * ts_xxzz_xyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyyz[i] * fe_0 + tk_xxxzz_yyyz[i] * fe_0 + tk_xxxzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyyyz[i] * fz_0;

        tk_xxxxzz_xyyzz[i] = -6.0 * ts_xxzz_xyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyzz[i] * fe_0 + tk_xxxzz_yyzz[i] * fe_0 + tk_xxxzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyyzz[i] * fz_0;

        tk_xxxxzz_xyzzz[i] = -6.0 * ts_xxzz_xyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyzzz[i] * fe_0 + tk_xxxzz_yzzz[i] * fe_0 + tk_xxxzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xyzzz[i] * fz_0;

        tk_xxxxzz_xzzzz[i] = -6.0 * ts_xxzz_xzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xzzzz[i] * fe_0 + tk_xxxzz_zzzz[i] * fe_0 + tk_xxxzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_xzzzz[i] * fz_0;

        tk_xxxxzz_yyyyy[i] = -6.0 * ts_xxzz_yyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyyy[i] * fe_0 + tk_xxxzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyyyy[i] * fz_0;

        tk_xxxxzz_yyyyz[i] = -6.0 * ts_xxzz_yyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyyz[i] * fe_0 + tk_xxxzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyyyz[i] * fz_0;

        tk_xxxxzz_yyyzz[i] = -6.0 * ts_xxzz_yyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyyzz[i] * fe_0 + tk_xxxzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyyzz[i] * fz_0;

        tk_xxxxzz_yyzzz[i] = -6.0 * ts_xxzz_yyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yyzzz[i] * fe_0 + tk_xxxzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yyzzz[i] * fz_0;

        tk_xxxxzz_yzzzz[i] = -6.0 * ts_xxzz_yzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yzzzz[i] * fe_0 + tk_xxxzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yzzzz[i] * fz_0;

        tk_xxxxzz_zzzzz[i] = -6.0 * ts_xxzz_zzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_zzzzz[i] * fe_0 + tk_xxxzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxxzz_zzzzz[i] * fz_0;
    }

    // Set up 126-147 components of targeted buffer : IH

    auto tk_xxxyyy_xxxxx = pbuffer.data(idx_kin_ih + 126);

    auto tk_xxxyyy_xxxxy = pbuffer.data(idx_kin_ih + 127);

    auto tk_xxxyyy_xxxxz = pbuffer.data(idx_kin_ih + 128);

    auto tk_xxxyyy_xxxyy = pbuffer.data(idx_kin_ih + 129);

    auto tk_xxxyyy_xxxyz = pbuffer.data(idx_kin_ih + 130);

    auto tk_xxxyyy_xxxzz = pbuffer.data(idx_kin_ih + 131);

    auto tk_xxxyyy_xxyyy = pbuffer.data(idx_kin_ih + 132);

    auto tk_xxxyyy_xxyyz = pbuffer.data(idx_kin_ih + 133);

    auto tk_xxxyyy_xxyzz = pbuffer.data(idx_kin_ih + 134);

    auto tk_xxxyyy_xxzzz = pbuffer.data(idx_kin_ih + 135);

    auto tk_xxxyyy_xyyyy = pbuffer.data(idx_kin_ih + 136);

    auto tk_xxxyyy_xyyyz = pbuffer.data(idx_kin_ih + 137);

    auto tk_xxxyyy_xyyzz = pbuffer.data(idx_kin_ih + 138);

    auto tk_xxxyyy_xyzzz = pbuffer.data(idx_kin_ih + 139);

    auto tk_xxxyyy_xzzzz = pbuffer.data(idx_kin_ih + 140);

    auto tk_xxxyyy_yyyyy = pbuffer.data(idx_kin_ih + 141);

    auto tk_xxxyyy_yyyyz = pbuffer.data(idx_kin_ih + 142);

    auto tk_xxxyyy_yyyzz = pbuffer.data(idx_kin_ih + 143);

    auto tk_xxxyyy_yyzzz = pbuffer.data(idx_kin_ih + 144);

    auto tk_xxxyyy_yzzzz = pbuffer.data(idx_kin_ih + 145);

    auto tk_xxxyyy_zzzzz = pbuffer.data(idx_kin_ih + 146);

    #pragma omp simd aligned(pa_x, pa_y, tk_xxxy_xxxxx, tk_xxxy_xxxxz, tk_xxxy_xxxzz, tk_xxxy_xxzzz, tk_xxxy_xzzzz, tk_xxxyy_xxxxx, tk_xxxyy_xxxxz, tk_xxxyy_xxxzz, tk_xxxyy_xxzzz, tk_xxxyy_xzzzz, tk_xxxyyy_xxxxx, tk_xxxyyy_xxxxy, tk_xxxyyy_xxxxz, tk_xxxyyy_xxxyy, tk_xxxyyy_xxxyz, tk_xxxyyy_xxxzz, tk_xxxyyy_xxyyy, tk_xxxyyy_xxyyz, tk_xxxyyy_xxyzz, tk_xxxyyy_xxzzz, tk_xxxyyy_xyyyy, tk_xxxyyy_xyyyz, tk_xxxyyy_xyyzz, tk_xxxyyy_xyzzz, tk_xxxyyy_xzzzz, tk_xxxyyy_yyyyy, tk_xxxyyy_yyyyz, tk_xxxyyy_yyyzz, tk_xxxyyy_yyzzz, tk_xxxyyy_yzzzz, tk_xxxyyy_zzzzz, tk_xxyyy_xxxxy, tk_xxyyy_xxxy, tk_xxyyy_xxxyy, tk_xxyyy_xxxyz, tk_xxyyy_xxyy, tk_xxyyy_xxyyy, tk_xxyyy_xxyyz, tk_xxyyy_xxyz, tk_xxyyy_xxyzz, tk_xxyyy_xyyy, tk_xxyyy_xyyyy, tk_xxyyy_xyyyz, tk_xxyyy_xyyz, tk_xxyyy_xyyzz, tk_xxyyy_xyzz, tk_xxyyy_xyzzz, tk_xxyyy_yyyy, tk_xxyyy_yyyyy, tk_xxyyy_yyyyz, tk_xxyyy_yyyz, tk_xxyyy_yyyzz, tk_xxyyy_yyzz, tk_xxyyy_yyzzz, tk_xxyyy_yzzz, tk_xxyyy_yzzzz, tk_xxyyy_zzzzz, tk_xyyy_xxxxy, tk_xyyy_xxxyy, tk_xyyy_xxxyz, tk_xyyy_xxyyy, tk_xyyy_xxyyz, tk_xyyy_xxyzz, tk_xyyy_xyyyy, tk_xyyy_xyyyz, tk_xyyy_xyyzz, tk_xyyy_xyzzz, tk_xyyy_yyyyy, tk_xyyy_yyyyz, tk_xyyy_yyyzz, tk_xyyy_yyzzz, tk_xyyy_yzzzz, tk_xyyy_zzzzz, ts_xxxy_xxxxx, ts_xxxy_xxxxz, ts_xxxy_xxxzz, ts_xxxy_xxzzz, ts_xxxy_xzzzz, ts_xxxyyy_xxxxx, ts_xxxyyy_xxxxy, ts_xxxyyy_xxxxz, ts_xxxyyy_xxxyy, ts_xxxyyy_xxxyz, ts_xxxyyy_xxxzz, ts_xxxyyy_xxyyy, ts_xxxyyy_xxyyz, ts_xxxyyy_xxyzz, ts_xxxyyy_xxzzz, ts_xxxyyy_xyyyy, ts_xxxyyy_xyyyz, ts_xxxyyy_xyyzz, ts_xxxyyy_xyzzz, ts_xxxyyy_xzzzz, ts_xxxyyy_yyyyy, ts_xxxyyy_yyyyz, ts_xxxyyy_yyyzz, ts_xxxyyy_yyzzz, ts_xxxyyy_yzzzz, ts_xxxyyy_zzzzz, ts_xyyy_xxxxy, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyzz, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyzz, ts_xyyy_xyzzz, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyzz, ts_xyyy_yyzzz, ts_xyyy_yzzzz, ts_xyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyyy_xxxxx[i] = -4.0 * ts_xxxy_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxxx[i] * fe_0 + tk_xxxyy_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxxxx[i] * fz_0;

        tk_xxxyyy_xxxxy[i] = -4.0 * ts_xyyy_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxxy[i] * fe_0 + 4.0 * tk_xxyyy_xxxy[i] * fe_0 + tk_xxyyy_xxxxy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxxy[i] * fz_0;

        tk_xxxyyy_xxxxz[i] = -4.0 * ts_xxxy_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxxz[i] * fe_0 + tk_xxxyy_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxxxz[i] * fz_0;

        tk_xxxyyy_xxxyy[i] = -4.0 * ts_xyyy_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxyy[i] * fe_0 + 3.0 * tk_xxyyy_xxyy[i] * fe_0 + tk_xxyyy_xxxyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxyy[i] * fz_0;

        tk_xxxyyy_xxxyz[i] = -4.0 * ts_xyyy_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxxyz[i] * fe_0 + 3.0 * tk_xxyyy_xxyz[i] * fe_0 + tk_xxyyy_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxxyz[i] * fz_0;

        tk_xxxyyy_xxxzz[i] = -4.0 * ts_xxxy_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxxzz[i] * fe_0 + tk_xxxyy_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxxzz[i] * fz_0;

        tk_xxxyyy_xxyyy[i] = -4.0 * ts_xyyy_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyyy[i] * fe_0 + 2.0 * tk_xxyyy_xyyy[i] * fe_0 + tk_xxyyy_xxyyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyyy[i] * fz_0;

        tk_xxxyyy_xxyyz[i] = -4.0 * ts_xyyy_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyyz[i] * fe_0 + 2.0 * tk_xxyyy_xyyz[i] * fe_0 + tk_xxyyy_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyyz[i] * fz_0;

        tk_xxxyyy_xxyzz[i] = -4.0 * ts_xyyy_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xxyzz[i] * fe_0 + 2.0 * tk_xxyyy_xyzz[i] * fe_0 + tk_xxyyy_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xxyzz[i] * fz_0;

        tk_xxxyyy_xxzzz[i] = -4.0 * ts_xxxy_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xxzzz[i] * fe_0 + tk_xxxyy_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xxzzz[i] * fz_0;

        tk_xxxyyy_xyyyy[i] = -4.0 * ts_xyyy_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyyy[i] * fe_0 + tk_xxyyy_yyyy[i] * fe_0 + tk_xxyyy_xyyyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyyy[i] * fz_0;

        tk_xxxyyy_xyyyz[i] = -4.0 * ts_xyyy_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyyz[i] * fe_0 + tk_xxyyy_yyyz[i] * fe_0 + tk_xxyyy_xyyyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyyz[i] * fz_0;

        tk_xxxyyy_xyyzz[i] = -4.0 * ts_xyyy_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyyzz[i] * fe_0 + tk_xxyyy_yyzz[i] * fe_0 + tk_xxyyy_xyyzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyyzz[i] * fz_0;

        tk_xxxyyy_xyzzz[i] = -4.0 * ts_xyyy_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xyzzz[i] * fe_0 + tk_xxyyy_yzzz[i] * fe_0 + tk_xxyyy_xyzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_xyzzz[i] * fz_0;

        tk_xxxyyy_xzzzz[i] = -4.0 * ts_xxxy_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xzzzz[i] * fe_0 + tk_xxxyy_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xzzzz[i] * fz_0;

        tk_xxxyyy_yyyyy[i] = -4.0 * ts_xyyy_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyyy[i] * fe_0 + tk_xxyyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyyyy[i] * fz_0;

        tk_xxxyyy_yyyyz[i] = -4.0 * ts_xyyy_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyyz[i] * fe_0 + tk_xxyyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyyyz[i] * fz_0;

        tk_xxxyyy_yyyzz[i] = -4.0 * ts_xyyy_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyyzz[i] * fe_0 + tk_xxyyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyyzz[i] * fz_0;

        tk_xxxyyy_yyzzz[i] = -4.0 * ts_xyyy_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yyzzz[i] * fe_0 + tk_xxyyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yyzzz[i] * fz_0;

        tk_xxxyyy_yzzzz[i] = -4.0 * ts_xyyy_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yzzzz[i] * fe_0 + tk_xxyyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yzzzz[i] * fz_0;

        tk_xxxyyy_zzzzz[i] = -4.0 * ts_xyyy_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_zzzzz[i] * fe_0 + tk_xxyyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxyyy_zzzzz[i] * fz_0;
    }

    // Set up 147-168 components of targeted buffer : IH

    auto tk_xxxyyz_xxxxx = pbuffer.data(idx_kin_ih + 147);

    auto tk_xxxyyz_xxxxy = pbuffer.data(idx_kin_ih + 148);

    auto tk_xxxyyz_xxxxz = pbuffer.data(idx_kin_ih + 149);

    auto tk_xxxyyz_xxxyy = pbuffer.data(idx_kin_ih + 150);

    auto tk_xxxyyz_xxxyz = pbuffer.data(idx_kin_ih + 151);

    auto tk_xxxyyz_xxxzz = pbuffer.data(idx_kin_ih + 152);

    auto tk_xxxyyz_xxyyy = pbuffer.data(idx_kin_ih + 153);

    auto tk_xxxyyz_xxyyz = pbuffer.data(idx_kin_ih + 154);

    auto tk_xxxyyz_xxyzz = pbuffer.data(idx_kin_ih + 155);

    auto tk_xxxyyz_xxzzz = pbuffer.data(idx_kin_ih + 156);

    auto tk_xxxyyz_xyyyy = pbuffer.data(idx_kin_ih + 157);

    auto tk_xxxyyz_xyyyz = pbuffer.data(idx_kin_ih + 158);

    auto tk_xxxyyz_xyyzz = pbuffer.data(idx_kin_ih + 159);

    auto tk_xxxyyz_xyzzz = pbuffer.data(idx_kin_ih + 160);

    auto tk_xxxyyz_xzzzz = pbuffer.data(idx_kin_ih + 161);

    auto tk_xxxyyz_yyyyy = pbuffer.data(idx_kin_ih + 162);

    auto tk_xxxyyz_yyyyz = pbuffer.data(idx_kin_ih + 163);

    auto tk_xxxyyz_yyyzz = pbuffer.data(idx_kin_ih + 164);

    auto tk_xxxyyz_yyzzz = pbuffer.data(idx_kin_ih + 165);

    auto tk_xxxyyz_yzzzz = pbuffer.data(idx_kin_ih + 166);

    auto tk_xxxyyz_zzzzz = pbuffer.data(idx_kin_ih + 167);

    #pragma omp simd aligned(pa_z, tk_xxxyy_xxxx, tk_xxxyy_xxxxx, tk_xxxyy_xxxxy, tk_xxxyy_xxxxz, tk_xxxyy_xxxy, tk_xxxyy_xxxyy, tk_xxxyy_xxxyz, tk_xxxyy_xxxz, tk_xxxyy_xxxzz, tk_xxxyy_xxyy, tk_xxxyy_xxyyy, tk_xxxyy_xxyyz, tk_xxxyy_xxyz, tk_xxxyy_xxyzz, tk_xxxyy_xxzz, tk_xxxyy_xxzzz, tk_xxxyy_xyyy, tk_xxxyy_xyyyy, tk_xxxyy_xyyyz, tk_xxxyy_xyyz, tk_xxxyy_xyyzz, tk_xxxyy_xyzz, tk_xxxyy_xyzzz, tk_xxxyy_xzzz, tk_xxxyy_xzzzz, tk_xxxyy_yyyy, tk_xxxyy_yyyyy, tk_xxxyy_yyyyz, tk_xxxyy_yyyz, tk_xxxyy_yyyzz, tk_xxxyy_yyzz, tk_xxxyy_yyzzz, tk_xxxyy_yzzz, tk_xxxyy_yzzzz, tk_xxxyy_zzzz, tk_xxxyy_zzzzz, tk_xxxyyz_xxxxx, tk_xxxyyz_xxxxy, tk_xxxyyz_xxxxz, tk_xxxyyz_xxxyy, tk_xxxyyz_xxxyz, tk_xxxyyz_xxxzz, tk_xxxyyz_xxyyy, tk_xxxyyz_xxyyz, tk_xxxyyz_xxyzz, tk_xxxyyz_xxzzz, tk_xxxyyz_xyyyy, tk_xxxyyz_xyyyz, tk_xxxyyz_xyyzz, tk_xxxyyz_xyzzz, tk_xxxyyz_xzzzz, tk_xxxyyz_yyyyy, tk_xxxyyz_yyyyz, tk_xxxyyz_yyyzz, tk_xxxyyz_yyzzz, tk_xxxyyz_yzzzz, tk_xxxyyz_zzzzz, ts_xxxyyz_xxxxx, ts_xxxyyz_xxxxy, ts_xxxyyz_xxxxz, ts_xxxyyz_xxxyy, ts_xxxyyz_xxxyz, ts_xxxyyz_xxxzz, ts_xxxyyz_xxyyy, ts_xxxyyz_xxyyz, ts_xxxyyz_xxyzz, ts_xxxyyz_xxzzz, ts_xxxyyz_xyyyy, ts_xxxyyz_xyyyz, ts_xxxyyz_xyyzz, ts_xxxyyz_xyzzz, ts_xxxyyz_xzzzz, ts_xxxyyz_yyyyy, ts_xxxyyz_yyyyz, ts_xxxyyz_yyyzz, ts_xxxyyz_yyzzz, ts_xxxyyz_yzzzz, ts_xxxyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyyz_xxxxx[i] = tk_xxxyy_xxxxx[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxx[i] * fz_0;

        tk_xxxyyz_xxxxy[i] = tk_xxxyy_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxy[i] * fz_0;

        tk_xxxyyz_xxxxz[i] = tk_xxxyy_xxxx[i] * fe_0 + tk_xxxyy_xxxxz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxxz[i] * fz_0;

        tk_xxxyyz_xxxyy[i] = tk_xxxyy_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxyy[i] * fz_0;

        tk_xxxyyz_xxxyz[i] = tk_xxxyy_xxxy[i] * fe_0 + tk_xxxyy_xxxyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxyz[i] * fz_0;

        tk_xxxyyz_xxxzz[i] = 2.0 * tk_xxxyy_xxxz[i] * fe_0 + tk_xxxyy_xxxzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxxzz[i] * fz_0;

        tk_xxxyyz_xxyyy[i] = tk_xxxyy_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyyy[i] * fz_0;

        tk_xxxyyz_xxyyz[i] = tk_xxxyy_xxyy[i] * fe_0 + tk_xxxyy_xxyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyyz[i] * fz_0;

        tk_xxxyyz_xxyzz[i] = 2.0 * tk_xxxyy_xxyz[i] * fe_0 + tk_xxxyy_xxyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxyzz[i] * fz_0;

        tk_xxxyyz_xxzzz[i] = 3.0 * tk_xxxyy_xxzz[i] * fe_0 + tk_xxxyy_xxzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xxzzz[i] * fz_0;

        tk_xxxyyz_xyyyy[i] = tk_xxxyy_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyyy[i] * fz_0;

        tk_xxxyyz_xyyyz[i] = tk_xxxyy_xyyy[i] * fe_0 + tk_xxxyy_xyyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyyz[i] * fz_0;

        tk_xxxyyz_xyyzz[i] = 2.0 * tk_xxxyy_xyyz[i] * fe_0 + tk_xxxyy_xyyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyyzz[i] * fz_0;

        tk_xxxyyz_xyzzz[i] = 3.0 * tk_xxxyy_xyzz[i] * fe_0 + tk_xxxyy_xyzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xyzzz[i] * fz_0;

        tk_xxxyyz_xzzzz[i] = 4.0 * tk_xxxyy_xzzz[i] * fe_0 + tk_xxxyy_xzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xzzzz[i] * fz_0;

        tk_xxxyyz_yyyyy[i] = tk_xxxyy_yyyyy[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyyy[i] * fz_0;

        tk_xxxyyz_yyyyz[i] = tk_xxxyy_yyyy[i] * fe_0 + tk_xxxyy_yyyyz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyyz[i] * fz_0;

        tk_xxxyyz_yyyzz[i] = 2.0 * tk_xxxyy_yyyz[i] * fe_0 + tk_xxxyy_yyyzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyyzz[i] * fz_0;

        tk_xxxyyz_yyzzz[i] = 3.0 * tk_xxxyy_yyzz[i] * fe_0 + tk_xxxyy_yyzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yyzzz[i] * fz_0;

        tk_xxxyyz_yzzzz[i] = 4.0 * tk_xxxyy_yzzz[i] * fe_0 + tk_xxxyy_yzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yzzzz[i] * fz_0;

        tk_xxxyyz_zzzzz[i] = 5.0 * tk_xxxyy_zzzz[i] * fe_0 + tk_xxxyy_zzzzz[i] * pa_z[i] + 2.0 * ts_xxxyyz_zzzzz[i] * fz_0;
    }

    // Set up 168-189 components of targeted buffer : IH

    auto tk_xxxyzz_xxxxx = pbuffer.data(idx_kin_ih + 168);

    auto tk_xxxyzz_xxxxy = pbuffer.data(idx_kin_ih + 169);

    auto tk_xxxyzz_xxxxz = pbuffer.data(idx_kin_ih + 170);

    auto tk_xxxyzz_xxxyy = pbuffer.data(idx_kin_ih + 171);

    auto tk_xxxyzz_xxxyz = pbuffer.data(idx_kin_ih + 172);

    auto tk_xxxyzz_xxxzz = pbuffer.data(idx_kin_ih + 173);

    auto tk_xxxyzz_xxyyy = pbuffer.data(idx_kin_ih + 174);

    auto tk_xxxyzz_xxyyz = pbuffer.data(idx_kin_ih + 175);

    auto tk_xxxyzz_xxyzz = pbuffer.data(idx_kin_ih + 176);

    auto tk_xxxyzz_xxzzz = pbuffer.data(idx_kin_ih + 177);

    auto tk_xxxyzz_xyyyy = pbuffer.data(idx_kin_ih + 178);

    auto tk_xxxyzz_xyyyz = pbuffer.data(idx_kin_ih + 179);

    auto tk_xxxyzz_xyyzz = pbuffer.data(idx_kin_ih + 180);

    auto tk_xxxyzz_xyzzz = pbuffer.data(idx_kin_ih + 181);

    auto tk_xxxyzz_xzzzz = pbuffer.data(idx_kin_ih + 182);

    auto tk_xxxyzz_yyyyy = pbuffer.data(idx_kin_ih + 183);

    auto tk_xxxyzz_yyyyz = pbuffer.data(idx_kin_ih + 184);

    auto tk_xxxyzz_yyyzz = pbuffer.data(idx_kin_ih + 185);

    auto tk_xxxyzz_yyzzz = pbuffer.data(idx_kin_ih + 186);

    auto tk_xxxyzz_yzzzz = pbuffer.data(idx_kin_ih + 187);

    auto tk_xxxyzz_zzzzz = pbuffer.data(idx_kin_ih + 188);

    #pragma omp simd aligned(pa_y, tk_xxxyzz_xxxxx, tk_xxxyzz_xxxxy, tk_xxxyzz_xxxxz, tk_xxxyzz_xxxyy, tk_xxxyzz_xxxyz, tk_xxxyzz_xxxzz, tk_xxxyzz_xxyyy, tk_xxxyzz_xxyyz, tk_xxxyzz_xxyzz, tk_xxxyzz_xxzzz, tk_xxxyzz_xyyyy, tk_xxxyzz_xyyyz, tk_xxxyzz_xyyzz, tk_xxxyzz_xyzzz, tk_xxxyzz_xzzzz, tk_xxxyzz_yyyyy, tk_xxxyzz_yyyyz, tk_xxxyzz_yyyzz, tk_xxxyzz_yyzzz, tk_xxxyzz_yzzzz, tk_xxxyzz_zzzzz, tk_xxxzz_xxxx, tk_xxxzz_xxxxx, tk_xxxzz_xxxxy, tk_xxxzz_xxxxz, tk_xxxzz_xxxy, tk_xxxzz_xxxyy, tk_xxxzz_xxxyz, tk_xxxzz_xxxz, tk_xxxzz_xxxzz, tk_xxxzz_xxyy, tk_xxxzz_xxyyy, tk_xxxzz_xxyyz, tk_xxxzz_xxyz, tk_xxxzz_xxyzz, tk_xxxzz_xxzz, tk_xxxzz_xxzzz, tk_xxxzz_xyyy, tk_xxxzz_xyyyy, tk_xxxzz_xyyyz, tk_xxxzz_xyyz, tk_xxxzz_xyyzz, tk_xxxzz_xyzz, tk_xxxzz_xyzzz, tk_xxxzz_xzzz, tk_xxxzz_xzzzz, tk_xxxzz_yyyy, tk_xxxzz_yyyyy, tk_xxxzz_yyyyz, tk_xxxzz_yyyz, tk_xxxzz_yyyzz, tk_xxxzz_yyzz, tk_xxxzz_yyzzz, tk_xxxzz_yzzz, tk_xxxzz_yzzzz, tk_xxxzz_zzzz, tk_xxxzz_zzzzz, ts_xxxyzz_xxxxx, ts_xxxyzz_xxxxy, ts_xxxyzz_xxxxz, ts_xxxyzz_xxxyy, ts_xxxyzz_xxxyz, ts_xxxyzz_xxxzz, ts_xxxyzz_xxyyy, ts_xxxyzz_xxyyz, ts_xxxyzz_xxyzz, ts_xxxyzz_xxzzz, ts_xxxyzz_xyyyy, ts_xxxyzz_xyyyz, ts_xxxyzz_xyyzz, ts_xxxyzz_xyzzz, ts_xxxyzz_xzzzz, ts_xxxyzz_yyyyy, ts_xxxyzz_yyyyz, ts_xxxyzz_yyyzz, ts_xxxyzz_yyzzz, ts_xxxyzz_yzzzz, ts_xxxyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyzz_xxxxx[i] = tk_xxxzz_xxxxx[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxx[i] * fz_0;

        tk_xxxyzz_xxxxy[i] = tk_xxxzz_xxxx[i] * fe_0 + tk_xxxzz_xxxxy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxy[i] * fz_0;

        tk_xxxyzz_xxxxz[i] = tk_xxxzz_xxxxz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxxz[i] * fz_0;

        tk_xxxyzz_xxxyy[i] = 2.0 * tk_xxxzz_xxxy[i] * fe_0 + tk_xxxzz_xxxyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxyy[i] * fz_0;

        tk_xxxyzz_xxxyz[i] = tk_xxxzz_xxxz[i] * fe_0 + tk_xxxzz_xxxyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxyz[i] * fz_0;

        tk_xxxyzz_xxxzz[i] = tk_xxxzz_xxxzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxxzz[i] * fz_0;

        tk_xxxyzz_xxyyy[i] = 3.0 * tk_xxxzz_xxyy[i] * fe_0 + tk_xxxzz_xxyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyyy[i] * fz_0;

        tk_xxxyzz_xxyyz[i] = 2.0 * tk_xxxzz_xxyz[i] * fe_0 + tk_xxxzz_xxyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyyz[i] * fz_0;

        tk_xxxyzz_xxyzz[i] = tk_xxxzz_xxzz[i] * fe_0 + tk_xxxzz_xxyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxyzz[i] * fz_0;

        tk_xxxyzz_xxzzz[i] = tk_xxxzz_xxzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xxzzz[i] * fz_0;

        tk_xxxyzz_xyyyy[i] = 4.0 * tk_xxxzz_xyyy[i] * fe_0 + tk_xxxzz_xyyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyyy[i] * fz_0;

        tk_xxxyzz_xyyyz[i] = 3.0 * tk_xxxzz_xyyz[i] * fe_0 + tk_xxxzz_xyyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyyz[i] * fz_0;

        tk_xxxyzz_xyyzz[i] = 2.0 * tk_xxxzz_xyzz[i] * fe_0 + tk_xxxzz_xyyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyyzz[i] * fz_0;

        tk_xxxyzz_xyzzz[i] = tk_xxxzz_xzzz[i] * fe_0 + tk_xxxzz_xyzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xyzzz[i] * fz_0;

        tk_xxxyzz_xzzzz[i] = tk_xxxzz_xzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xzzzz[i] * fz_0;

        tk_xxxyzz_yyyyy[i] = 5.0 * tk_xxxzz_yyyy[i] * fe_0 + tk_xxxzz_yyyyy[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyyy[i] * fz_0;

        tk_xxxyzz_yyyyz[i] = 4.0 * tk_xxxzz_yyyz[i] * fe_0 + tk_xxxzz_yyyyz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyyz[i] * fz_0;

        tk_xxxyzz_yyyzz[i] = 3.0 * tk_xxxzz_yyzz[i] * fe_0 + tk_xxxzz_yyyzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyyzz[i] * fz_0;

        tk_xxxyzz_yyzzz[i] = 2.0 * tk_xxxzz_yzzz[i] * fe_0 + tk_xxxzz_yyzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yyzzz[i] * fz_0;

        tk_xxxyzz_yzzzz[i] = tk_xxxzz_zzzz[i] * fe_0 + tk_xxxzz_yzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yzzzz[i] * fz_0;

        tk_xxxyzz_zzzzz[i] = tk_xxxzz_zzzzz[i] * pa_y[i] + 2.0 * ts_xxxyzz_zzzzz[i] * fz_0;
    }

    // Set up 189-210 components of targeted buffer : IH

    auto tk_xxxzzz_xxxxx = pbuffer.data(idx_kin_ih + 189);

    auto tk_xxxzzz_xxxxy = pbuffer.data(idx_kin_ih + 190);

    auto tk_xxxzzz_xxxxz = pbuffer.data(idx_kin_ih + 191);

    auto tk_xxxzzz_xxxyy = pbuffer.data(idx_kin_ih + 192);

    auto tk_xxxzzz_xxxyz = pbuffer.data(idx_kin_ih + 193);

    auto tk_xxxzzz_xxxzz = pbuffer.data(idx_kin_ih + 194);

    auto tk_xxxzzz_xxyyy = pbuffer.data(idx_kin_ih + 195);

    auto tk_xxxzzz_xxyyz = pbuffer.data(idx_kin_ih + 196);

    auto tk_xxxzzz_xxyzz = pbuffer.data(idx_kin_ih + 197);

    auto tk_xxxzzz_xxzzz = pbuffer.data(idx_kin_ih + 198);

    auto tk_xxxzzz_xyyyy = pbuffer.data(idx_kin_ih + 199);

    auto tk_xxxzzz_xyyyz = pbuffer.data(idx_kin_ih + 200);

    auto tk_xxxzzz_xyyzz = pbuffer.data(idx_kin_ih + 201);

    auto tk_xxxzzz_xyzzz = pbuffer.data(idx_kin_ih + 202);

    auto tk_xxxzzz_xzzzz = pbuffer.data(idx_kin_ih + 203);

    auto tk_xxxzzz_yyyyy = pbuffer.data(idx_kin_ih + 204);

    auto tk_xxxzzz_yyyyz = pbuffer.data(idx_kin_ih + 205);

    auto tk_xxxzzz_yyyzz = pbuffer.data(idx_kin_ih + 206);

    auto tk_xxxzzz_yyzzz = pbuffer.data(idx_kin_ih + 207);

    auto tk_xxxzzz_yzzzz = pbuffer.data(idx_kin_ih + 208);

    auto tk_xxxzzz_zzzzz = pbuffer.data(idx_kin_ih + 209);

    #pragma omp simd aligned(pa_x, pa_z, tk_xxxz_xxxxx, tk_xxxz_xxxxy, tk_xxxz_xxxyy, tk_xxxz_xxyyy, tk_xxxz_xyyyy, tk_xxxzz_xxxxx, tk_xxxzz_xxxxy, tk_xxxzz_xxxyy, tk_xxxzz_xxyyy, tk_xxxzz_xyyyy, tk_xxxzzz_xxxxx, tk_xxxzzz_xxxxy, tk_xxxzzz_xxxxz, tk_xxxzzz_xxxyy, tk_xxxzzz_xxxyz, tk_xxxzzz_xxxzz, tk_xxxzzz_xxyyy, tk_xxxzzz_xxyyz, tk_xxxzzz_xxyzz, tk_xxxzzz_xxzzz, tk_xxxzzz_xyyyy, tk_xxxzzz_xyyyz, tk_xxxzzz_xyyzz, tk_xxxzzz_xyzzz, tk_xxxzzz_xzzzz, tk_xxxzzz_yyyyy, tk_xxxzzz_yyyyz, tk_xxxzzz_yyyzz, tk_xxxzzz_yyzzz, tk_xxxzzz_yzzzz, tk_xxxzzz_zzzzz, tk_xxzzz_xxxxz, tk_xxzzz_xxxyz, tk_xxzzz_xxxz, tk_xxzzz_xxxzz, tk_xxzzz_xxyyz, tk_xxzzz_xxyz, tk_xxzzz_xxyzz, tk_xxzzz_xxzz, tk_xxzzz_xxzzz, tk_xxzzz_xyyyz, tk_xxzzz_xyyz, tk_xxzzz_xyyzz, tk_xxzzz_xyzz, tk_xxzzz_xyzzz, tk_xxzzz_xzzz, tk_xxzzz_xzzzz, tk_xxzzz_yyyyy, tk_xxzzz_yyyyz, tk_xxzzz_yyyz, tk_xxzzz_yyyzz, tk_xxzzz_yyzz, tk_xxzzz_yyzzz, tk_xxzzz_yzzz, tk_xxzzz_yzzzz, tk_xxzzz_zzzz, tk_xxzzz_zzzzz, tk_xzzz_xxxxz, tk_xzzz_xxxyz, tk_xzzz_xxxzz, tk_xzzz_xxyyz, tk_xzzz_xxyzz, tk_xzzz_xxzzz, tk_xzzz_xyyyz, tk_xzzz_xyyzz, tk_xzzz_xyzzz, tk_xzzz_xzzzz, tk_xzzz_yyyyy, tk_xzzz_yyyyz, tk_xzzz_yyyzz, tk_xzzz_yyzzz, tk_xzzz_yzzzz, tk_xzzz_zzzzz, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxyy, ts_xxxz_xxyyy, ts_xxxz_xyyyy, ts_xxxzzz_xxxxx, ts_xxxzzz_xxxxy, ts_xxxzzz_xxxxz, ts_xxxzzz_xxxyy, ts_xxxzzz_xxxyz, ts_xxxzzz_xxxzz, ts_xxxzzz_xxyyy, ts_xxxzzz_xxyyz, ts_xxxzzz_xxyzz, ts_xxxzzz_xxzzz, ts_xxxzzz_xyyyy, ts_xxxzzz_xyyyz, ts_xxxzzz_xyyzz, ts_xxxzzz_xyzzz, ts_xxxzzz_xzzzz, ts_xxxzzz_yyyyy, ts_xxxzzz_yyyyz, ts_xxxzzz_yyyzz, ts_xxxzzz_yyzzz, ts_xxxzzz_yzzzz, ts_xxxzzz_zzzzz, ts_xzzz_xxxxz, ts_xzzz_xxxyz, ts_xzzz_xxxzz, ts_xzzz_xxyyz, ts_xzzz_xxyzz, ts_xzzz_xxzzz, ts_xzzz_xyyyz, ts_xzzz_xyyzz, ts_xzzz_xyzzz, ts_xzzz_xzzzz, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyzz, ts_xzzz_yyzzz, ts_xzzz_yzzzz, ts_xzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzzz_xxxxx[i] = -4.0 * ts_xxxz_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxxx[i] * fe_0 + tk_xxxzz_xxxxx[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxxxx[i] * fz_0;

        tk_xxxzzz_xxxxy[i] = -4.0 * ts_xxxz_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxxy[i] * fe_0 + tk_xxxzz_xxxxy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxxxy[i] * fz_0;

        tk_xxxzzz_xxxxz[i] = -4.0 * ts_xzzz_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxxz[i] * fe_0 + 4.0 * tk_xxzzz_xxxz[i] * fe_0 + tk_xxzzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxxz[i] * fz_0;

        tk_xxxzzz_xxxyy[i] = -4.0 * ts_xxxz_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxxyy[i] * fe_0 + tk_xxxzz_xxxyy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxxyy[i] * fz_0;

        tk_xxxzzz_xxxyz[i] = -4.0 * ts_xzzz_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxyz[i] * fe_0 + 3.0 * tk_xxzzz_xxyz[i] * fe_0 + tk_xxzzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxyz[i] * fz_0;

        tk_xxxzzz_xxxzz[i] = -4.0 * ts_xzzz_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxxzz[i] * fe_0 + 3.0 * tk_xxzzz_xxzz[i] * fe_0 + tk_xxzzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxxzz[i] * fz_0;

        tk_xxxzzz_xxyyy[i] = -4.0 * ts_xxxz_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xxyyy[i] * fe_0 + tk_xxxzz_xxyyy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xxyyy[i] * fz_0;

        tk_xxxzzz_xxyyz[i] = -4.0 * ts_xzzz_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxyyz[i] * fe_0 + 2.0 * tk_xxzzz_xyyz[i] * fe_0 + tk_xxzzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxyyz[i] * fz_0;

        tk_xxxzzz_xxyzz[i] = -4.0 * ts_xzzz_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxyzz[i] * fe_0 + 2.0 * tk_xxzzz_xyzz[i] * fe_0 + tk_xxzzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxyzz[i] * fz_0;

        tk_xxxzzz_xxzzz[i] = -4.0 * ts_xzzz_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xxzzz[i] * fe_0 + 2.0 * tk_xxzzz_xzzz[i] * fe_0 + tk_xxzzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xxzzz[i] * fz_0;

        tk_xxxzzz_xyyyy[i] = -4.0 * ts_xxxz_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xyyyy[i] * fe_0 + tk_xxxzz_xyyyy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xyyyy[i] * fz_0;

        tk_xxxzzz_xyyyz[i] = -4.0 * ts_xzzz_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyyyz[i] * fe_0 + tk_xxzzz_yyyz[i] * fe_0 + tk_xxzzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyyyz[i] * fz_0;

        tk_xxxzzz_xyyzz[i] = -4.0 * ts_xzzz_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyyzz[i] * fe_0 + tk_xxzzz_yyzz[i] * fe_0 + tk_xxzzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyyzz[i] * fz_0;

        tk_xxxzzz_xyzzz[i] = -4.0 * ts_xzzz_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xyzzz[i] * fe_0 + tk_xxzzz_yzzz[i] * fe_0 + tk_xxzzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xyzzz[i] * fz_0;

        tk_xxxzzz_xzzzz[i] = -4.0 * ts_xzzz_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xzzzz[i] * fe_0 + tk_xxzzz_zzzz[i] * fe_0 + tk_xxzzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_xzzzz[i] * fz_0;

        tk_xxxzzz_yyyyy[i] = -4.0 * ts_xzzz_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyyy[i] * fe_0 + tk_xxzzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyyyy[i] * fz_0;

        tk_xxxzzz_yyyyz[i] = -4.0 * ts_xzzz_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyyz[i] * fe_0 + tk_xxzzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyyyz[i] * fz_0;

        tk_xxxzzz_yyyzz[i] = -4.0 * ts_xzzz_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyyzz[i] * fe_0 + tk_xxzzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyyzz[i] * fz_0;

        tk_xxxzzz_yyzzz[i] = -4.0 * ts_xzzz_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yyzzz[i] * fe_0 + tk_xxzzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yyzzz[i] * fz_0;

        tk_xxxzzz_yzzzz[i] = -4.0 * ts_xzzz_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yzzzz[i] * fe_0 + tk_xxzzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yzzzz[i] * fz_0;

        tk_xxxzzz_zzzzz[i] = -4.0 * ts_xzzz_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_zzzzz[i] * fe_0 + tk_xxzzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xxxzzz_zzzzz[i] * fz_0;
    }

    // Set up 210-231 components of targeted buffer : IH

    auto tk_xxyyyy_xxxxx = pbuffer.data(idx_kin_ih + 210);

    auto tk_xxyyyy_xxxxy = pbuffer.data(idx_kin_ih + 211);

    auto tk_xxyyyy_xxxxz = pbuffer.data(idx_kin_ih + 212);

    auto tk_xxyyyy_xxxyy = pbuffer.data(idx_kin_ih + 213);

    auto tk_xxyyyy_xxxyz = pbuffer.data(idx_kin_ih + 214);

    auto tk_xxyyyy_xxxzz = pbuffer.data(idx_kin_ih + 215);

    auto tk_xxyyyy_xxyyy = pbuffer.data(idx_kin_ih + 216);

    auto tk_xxyyyy_xxyyz = pbuffer.data(idx_kin_ih + 217);

    auto tk_xxyyyy_xxyzz = pbuffer.data(idx_kin_ih + 218);

    auto tk_xxyyyy_xxzzz = pbuffer.data(idx_kin_ih + 219);

    auto tk_xxyyyy_xyyyy = pbuffer.data(idx_kin_ih + 220);

    auto tk_xxyyyy_xyyyz = pbuffer.data(idx_kin_ih + 221);

    auto tk_xxyyyy_xyyzz = pbuffer.data(idx_kin_ih + 222);

    auto tk_xxyyyy_xyzzz = pbuffer.data(idx_kin_ih + 223);

    auto tk_xxyyyy_xzzzz = pbuffer.data(idx_kin_ih + 224);

    auto tk_xxyyyy_yyyyy = pbuffer.data(idx_kin_ih + 225);

    auto tk_xxyyyy_yyyyz = pbuffer.data(idx_kin_ih + 226);

    auto tk_xxyyyy_yyyzz = pbuffer.data(idx_kin_ih + 227);

    auto tk_xxyyyy_yyzzz = pbuffer.data(idx_kin_ih + 228);

    auto tk_xxyyyy_yzzzz = pbuffer.data(idx_kin_ih + 229);

    auto tk_xxyyyy_zzzzz = pbuffer.data(idx_kin_ih + 230);

    #pragma omp simd aligned(pa_x, pa_y, tk_xxyy_xxxxx, tk_xxyy_xxxxz, tk_xxyy_xxxzz, tk_xxyy_xxzzz, tk_xxyy_xzzzz, tk_xxyyy_xxxxx, tk_xxyyy_xxxxz, tk_xxyyy_xxxzz, tk_xxyyy_xxzzz, tk_xxyyy_xzzzz, tk_xxyyyy_xxxxx, tk_xxyyyy_xxxxy, tk_xxyyyy_xxxxz, tk_xxyyyy_xxxyy, tk_xxyyyy_xxxyz, tk_xxyyyy_xxxzz, tk_xxyyyy_xxyyy, tk_xxyyyy_xxyyz, tk_xxyyyy_xxyzz, tk_xxyyyy_xxzzz, tk_xxyyyy_xyyyy, tk_xxyyyy_xyyyz, tk_xxyyyy_xyyzz, tk_xxyyyy_xyzzz, tk_xxyyyy_xzzzz, tk_xxyyyy_yyyyy, tk_xxyyyy_yyyyz, tk_xxyyyy_yyyzz, tk_xxyyyy_yyzzz, tk_xxyyyy_yzzzz, tk_xxyyyy_zzzzz, tk_xyyyy_xxxxy, tk_xyyyy_xxxy, tk_xyyyy_xxxyy, tk_xyyyy_xxxyz, tk_xyyyy_xxyy, tk_xyyyy_xxyyy, tk_xyyyy_xxyyz, tk_xyyyy_xxyz, tk_xyyyy_xxyzz, tk_xyyyy_xyyy, tk_xyyyy_xyyyy, tk_xyyyy_xyyyz, tk_xyyyy_xyyz, tk_xyyyy_xyyzz, tk_xyyyy_xyzz, tk_xyyyy_xyzzz, tk_xyyyy_yyyy, tk_xyyyy_yyyyy, tk_xyyyy_yyyyz, tk_xyyyy_yyyz, tk_xyyyy_yyyzz, tk_xyyyy_yyzz, tk_xyyyy_yyzzz, tk_xyyyy_yzzz, tk_xyyyy_yzzzz, tk_xyyyy_zzzzz, tk_yyyy_xxxxy, tk_yyyy_xxxyy, tk_yyyy_xxxyz, tk_yyyy_xxyyy, tk_yyyy_xxyyz, tk_yyyy_xxyzz, tk_yyyy_xyyyy, tk_yyyy_xyyyz, tk_yyyy_xyyzz, tk_yyyy_xyzzz, tk_yyyy_yyyyy, tk_yyyy_yyyyz, tk_yyyy_yyyzz, tk_yyyy_yyzzz, tk_yyyy_yzzzz, tk_yyyy_zzzzz, ts_xxyy_xxxxx, ts_xxyy_xxxxz, ts_xxyy_xxxzz, ts_xxyy_xxzzz, ts_xxyy_xzzzz, ts_xxyyyy_xxxxx, ts_xxyyyy_xxxxy, ts_xxyyyy_xxxxz, ts_xxyyyy_xxxyy, ts_xxyyyy_xxxyz, ts_xxyyyy_xxxzz, ts_xxyyyy_xxyyy, ts_xxyyyy_xxyyz, ts_xxyyyy_xxyzz, ts_xxyyyy_xxzzz, ts_xxyyyy_xyyyy, ts_xxyyyy_xyyyz, ts_xxyyyy_xyyzz, ts_xxyyyy_xyzzz, ts_xxyyyy_xzzzz, ts_xxyyyy_yyyyy, ts_xxyyyy_yyyyz, ts_xxyyyy_yyyzz, ts_xxyyyy_yyzzz, ts_xxyyyy_yzzzz, ts_xxyyyy_zzzzz, ts_yyyy_xxxxy, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyzz, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyzz, ts_yyyy_xyzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyzz, ts_yyyy_yyzzz, ts_yyyy_yzzzz, ts_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyyy_xxxxx[i] = -6.0 * ts_xxyy_xxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxx[i] * fe_0 + tk_xxyyy_xxxxx[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxxxx[i] * fz_0;

        tk_xxyyyy_xxxxy[i] = -2.0 * ts_yyyy_xxxxy[i] * fbe_0 * fz_0 + tk_yyyy_xxxxy[i] * fe_0 + 4.0 * tk_xyyyy_xxxy[i] * fe_0 + tk_xyyyy_xxxxy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxxy[i] * fz_0;

        tk_xxyyyy_xxxxz[i] = -6.0 * ts_xxyy_xxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxxz[i] * fe_0 + tk_xxyyy_xxxxz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxxxz[i] * fz_0;

        tk_xxyyyy_xxxyy[i] = -2.0 * ts_yyyy_xxxyy[i] * fbe_0 * fz_0 + tk_yyyy_xxxyy[i] * fe_0 + 3.0 * tk_xyyyy_xxyy[i] * fe_0 + tk_xyyyy_xxxyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxyy[i] * fz_0;

        tk_xxyyyy_xxxyz[i] = -2.0 * ts_yyyy_xxxyz[i] * fbe_0 * fz_0 + tk_yyyy_xxxyz[i] * fe_0 + 3.0 * tk_xyyyy_xxyz[i] * fe_0 + tk_xyyyy_xxxyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxxyz[i] * fz_0;

        tk_xxyyyy_xxxzz[i] = -6.0 * ts_xxyy_xxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxxzz[i] * fe_0 + tk_xxyyy_xxxzz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxxzz[i] * fz_0;

        tk_xxyyyy_xxyyy[i] = -2.0 * ts_yyyy_xxyyy[i] * fbe_0 * fz_0 + tk_yyyy_xxyyy[i] * fe_0 + 2.0 * tk_xyyyy_xyyy[i] * fe_0 + tk_xyyyy_xxyyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyyy[i] * fz_0;

        tk_xxyyyy_xxyyz[i] = -2.0 * ts_yyyy_xxyyz[i] * fbe_0 * fz_0 + tk_yyyy_xxyyz[i] * fe_0 + 2.0 * tk_xyyyy_xyyz[i] * fe_0 + tk_xyyyy_xxyyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyyz[i] * fz_0;

        tk_xxyyyy_xxyzz[i] = -2.0 * ts_yyyy_xxyzz[i] * fbe_0 * fz_0 + tk_yyyy_xxyzz[i] * fe_0 + 2.0 * tk_xyyyy_xyzz[i] * fe_0 + tk_xyyyy_xxyzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xxyzz[i] * fz_0;

        tk_xxyyyy_xxzzz[i] = -6.0 * ts_xxyy_xxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xxzzz[i] * fe_0 + tk_xxyyy_xxzzz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xxzzz[i] * fz_0;

        tk_xxyyyy_xyyyy[i] = -2.0 * ts_yyyy_xyyyy[i] * fbe_0 * fz_0 + tk_yyyy_xyyyy[i] * fe_0 + tk_xyyyy_yyyy[i] * fe_0 + tk_xyyyy_xyyyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyyyy[i] * fz_0;

        tk_xxyyyy_xyyyz[i] = -2.0 * ts_yyyy_xyyyz[i] * fbe_0 * fz_0 + tk_yyyy_xyyyz[i] * fe_0 + tk_xyyyy_yyyz[i] * fe_0 + tk_xyyyy_xyyyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyyyz[i] * fz_0;

        tk_xxyyyy_xyyzz[i] = -2.0 * ts_yyyy_xyyzz[i] * fbe_0 * fz_0 + tk_yyyy_xyyzz[i] * fe_0 + tk_xyyyy_yyzz[i] * fe_0 + tk_xyyyy_xyyzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyyzz[i] * fz_0;

        tk_xxyyyy_xyzzz[i] = -2.0 * ts_yyyy_xyzzz[i] * fbe_0 * fz_0 + tk_yyyy_xyzzz[i] * fe_0 + tk_xyyyy_yzzz[i] * fe_0 + tk_xyyyy_xyzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_xyzzz[i] * fz_0;

        tk_xxyyyy_xzzzz[i] = -6.0 * ts_xxyy_xzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xzzzz[i] * fe_0 + tk_xxyyy_xzzzz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xzzzz[i] * fz_0;

        tk_xxyyyy_yyyyy[i] = -2.0 * ts_yyyy_yyyyy[i] * fbe_0 * fz_0 + tk_yyyy_yyyyy[i] * fe_0 + tk_xyyyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyyy[i] * fz_0;

        tk_xxyyyy_yyyyz[i] = -2.0 * ts_yyyy_yyyyz[i] * fbe_0 * fz_0 + tk_yyyy_yyyyz[i] * fe_0 + tk_xyyyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyyz[i] * fz_0;

        tk_xxyyyy_yyyzz[i] = -2.0 * ts_yyyy_yyyzz[i] * fbe_0 * fz_0 + tk_yyyy_yyyzz[i] * fe_0 + tk_xyyyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyyzz[i] * fz_0;

        tk_xxyyyy_yyzzz[i] = -2.0 * ts_yyyy_yyzzz[i] * fbe_0 * fz_0 + tk_yyyy_yyzzz[i] * fe_0 + tk_xyyyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yyzzz[i] * fz_0;

        tk_xxyyyy_yzzzz[i] = -2.0 * ts_yyyy_yzzzz[i] * fbe_0 * fz_0 + tk_yyyy_yzzzz[i] * fe_0 + tk_xyyyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yzzzz[i] * fz_0;

        tk_xxyyyy_zzzzz[i] = -2.0 * ts_yyyy_zzzzz[i] * fbe_0 * fz_0 + tk_yyyy_zzzzz[i] * fe_0 + tk_xyyyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xxyyyy_zzzzz[i] * fz_0;
    }

    // Set up 231-252 components of targeted buffer : IH

    auto tk_xxyyyz_xxxxx = pbuffer.data(idx_kin_ih + 231);

    auto tk_xxyyyz_xxxxy = pbuffer.data(idx_kin_ih + 232);

    auto tk_xxyyyz_xxxxz = pbuffer.data(idx_kin_ih + 233);

    auto tk_xxyyyz_xxxyy = pbuffer.data(idx_kin_ih + 234);

    auto tk_xxyyyz_xxxyz = pbuffer.data(idx_kin_ih + 235);

    auto tk_xxyyyz_xxxzz = pbuffer.data(idx_kin_ih + 236);

    auto tk_xxyyyz_xxyyy = pbuffer.data(idx_kin_ih + 237);

    auto tk_xxyyyz_xxyyz = pbuffer.data(idx_kin_ih + 238);

    auto tk_xxyyyz_xxyzz = pbuffer.data(idx_kin_ih + 239);

    auto tk_xxyyyz_xxzzz = pbuffer.data(idx_kin_ih + 240);

    auto tk_xxyyyz_xyyyy = pbuffer.data(idx_kin_ih + 241);

    auto tk_xxyyyz_xyyyz = pbuffer.data(idx_kin_ih + 242);

    auto tk_xxyyyz_xyyzz = pbuffer.data(idx_kin_ih + 243);

    auto tk_xxyyyz_xyzzz = pbuffer.data(idx_kin_ih + 244);

    auto tk_xxyyyz_xzzzz = pbuffer.data(idx_kin_ih + 245);

    auto tk_xxyyyz_yyyyy = pbuffer.data(idx_kin_ih + 246);

    auto tk_xxyyyz_yyyyz = pbuffer.data(idx_kin_ih + 247);

    auto tk_xxyyyz_yyyzz = pbuffer.data(idx_kin_ih + 248);

    auto tk_xxyyyz_yyzzz = pbuffer.data(idx_kin_ih + 249);

    auto tk_xxyyyz_yzzzz = pbuffer.data(idx_kin_ih + 250);

    auto tk_xxyyyz_zzzzz = pbuffer.data(idx_kin_ih + 251);

    #pragma omp simd aligned(pa_z, tk_xxyyy_xxxx, tk_xxyyy_xxxxx, tk_xxyyy_xxxxy, tk_xxyyy_xxxxz, tk_xxyyy_xxxy, tk_xxyyy_xxxyy, tk_xxyyy_xxxyz, tk_xxyyy_xxxz, tk_xxyyy_xxxzz, tk_xxyyy_xxyy, tk_xxyyy_xxyyy, tk_xxyyy_xxyyz, tk_xxyyy_xxyz, tk_xxyyy_xxyzz, tk_xxyyy_xxzz, tk_xxyyy_xxzzz, tk_xxyyy_xyyy, tk_xxyyy_xyyyy, tk_xxyyy_xyyyz, tk_xxyyy_xyyz, tk_xxyyy_xyyzz, tk_xxyyy_xyzz, tk_xxyyy_xyzzz, tk_xxyyy_xzzz, tk_xxyyy_xzzzz, tk_xxyyy_yyyy, tk_xxyyy_yyyyy, tk_xxyyy_yyyyz, tk_xxyyy_yyyz, tk_xxyyy_yyyzz, tk_xxyyy_yyzz, tk_xxyyy_yyzzz, tk_xxyyy_yzzz, tk_xxyyy_yzzzz, tk_xxyyy_zzzz, tk_xxyyy_zzzzz, tk_xxyyyz_xxxxx, tk_xxyyyz_xxxxy, tk_xxyyyz_xxxxz, tk_xxyyyz_xxxyy, tk_xxyyyz_xxxyz, tk_xxyyyz_xxxzz, tk_xxyyyz_xxyyy, tk_xxyyyz_xxyyz, tk_xxyyyz_xxyzz, tk_xxyyyz_xxzzz, tk_xxyyyz_xyyyy, tk_xxyyyz_xyyyz, tk_xxyyyz_xyyzz, tk_xxyyyz_xyzzz, tk_xxyyyz_xzzzz, tk_xxyyyz_yyyyy, tk_xxyyyz_yyyyz, tk_xxyyyz_yyyzz, tk_xxyyyz_yyzzz, tk_xxyyyz_yzzzz, tk_xxyyyz_zzzzz, ts_xxyyyz_xxxxx, ts_xxyyyz_xxxxy, ts_xxyyyz_xxxxz, ts_xxyyyz_xxxyy, ts_xxyyyz_xxxyz, ts_xxyyyz_xxxzz, ts_xxyyyz_xxyyy, ts_xxyyyz_xxyyz, ts_xxyyyz_xxyzz, ts_xxyyyz_xxzzz, ts_xxyyyz_xyyyy, ts_xxyyyz_xyyyz, ts_xxyyyz_xyyzz, ts_xxyyyz_xyzzz, ts_xxyyyz_xzzzz, ts_xxyyyz_yyyyy, ts_xxyyyz_yyyyz, ts_xxyyyz_yyyzz, ts_xxyyyz_yyzzz, ts_xxyyyz_yzzzz, ts_xxyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyyz_xxxxx[i] = tk_xxyyy_xxxxx[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxx[i] * fz_0;

        tk_xxyyyz_xxxxy[i] = tk_xxyyy_xxxxy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxy[i] * fz_0;

        tk_xxyyyz_xxxxz[i] = tk_xxyyy_xxxx[i] * fe_0 + tk_xxyyy_xxxxz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxxz[i] * fz_0;

        tk_xxyyyz_xxxyy[i] = tk_xxyyy_xxxyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxyy[i] * fz_0;

        tk_xxyyyz_xxxyz[i] = tk_xxyyy_xxxy[i] * fe_0 + tk_xxyyy_xxxyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxyz[i] * fz_0;

        tk_xxyyyz_xxxzz[i] = 2.0 * tk_xxyyy_xxxz[i] * fe_0 + tk_xxyyy_xxxzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxxzz[i] * fz_0;

        tk_xxyyyz_xxyyy[i] = tk_xxyyy_xxyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyyy[i] * fz_0;

        tk_xxyyyz_xxyyz[i] = tk_xxyyy_xxyy[i] * fe_0 + tk_xxyyy_xxyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyyz[i] * fz_0;

        tk_xxyyyz_xxyzz[i] = 2.0 * tk_xxyyy_xxyz[i] * fe_0 + tk_xxyyy_xxyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxyzz[i] * fz_0;

        tk_xxyyyz_xxzzz[i] = 3.0 * tk_xxyyy_xxzz[i] * fe_0 + tk_xxyyy_xxzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xxzzz[i] * fz_0;

        tk_xxyyyz_xyyyy[i] = tk_xxyyy_xyyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyyy[i] * fz_0;

        tk_xxyyyz_xyyyz[i] = tk_xxyyy_xyyy[i] * fe_0 + tk_xxyyy_xyyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyyz[i] * fz_0;

        tk_xxyyyz_xyyzz[i] = 2.0 * tk_xxyyy_xyyz[i] * fe_0 + tk_xxyyy_xyyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyyzz[i] * fz_0;

        tk_xxyyyz_xyzzz[i] = 3.0 * tk_xxyyy_xyzz[i] * fe_0 + tk_xxyyy_xyzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xyzzz[i] * fz_0;

        tk_xxyyyz_xzzzz[i] = 4.0 * tk_xxyyy_xzzz[i] * fe_0 + tk_xxyyy_xzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xzzzz[i] * fz_0;

        tk_xxyyyz_yyyyy[i] = tk_xxyyy_yyyyy[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyyy[i] * fz_0;

        tk_xxyyyz_yyyyz[i] = tk_xxyyy_yyyy[i] * fe_0 + tk_xxyyy_yyyyz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyyz[i] * fz_0;

        tk_xxyyyz_yyyzz[i] = 2.0 * tk_xxyyy_yyyz[i] * fe_0 + tk_xxyyy_yyyzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyyzz[i] * fz_0;

        tk_xxyyyz_yyzzz[i] = 3.0 * tk_xxyyy_yyzz[i] * fe_0 + tk_xxyyy_yyzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yyzzz[i] * fz_0;

        tk_xxyyyz_yzzzz[i] = 4.0 * tk_xxyyy_yzzz[i] * fe_0 + tk_xxyyy_yzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yzzzz[i] * fz_0;

        tk_xxyyyz_zzzzz[i] = 5.0 * tk_xxyyy_zzzz[i] * fe_0 + tk_xxyyy_zzzzz[i] * pa_z[i] + 2.0 * ts_xxyyyz_zzzzz[i] * fz_0;
    }

    // Set up 252-273 components of targeted buffer : IH

    auto tk_xxyyzz_xxxxx = pbuffer.data(idx_kin_ih + 252);

    auto tk_xxyyzz_xxxxy = pbuffer.data(idx_kin_ih + 253);

    auto tk_xxyyzz_xxxxz = pbuffer.data(idx_kin_ih + 254);

    auto tk_xxyyzz_xxxyy = pbuffer.data(idx_kin_ih + 255);

    auto tk_xxyyzz_xxxyz = pbuffer.data(idx_kin_ih + 256);

    auto tk_xxyyzz_xxxzz = pbuffer.data(idx_kin_ih + 257);

    auto tk_xxyyzz_xxyyy = pbuffer.data(idx_kin_ih + 258);

    auto tk_xxyyzz_xxyyz = pbuffer.data(idx_kin_ih + 259);

    auto tk_xxyyzz_xxyzz = pbuffer.data(idx_kin_ih + 260);

    auto tk_xxyyzz_xxzzz = pbuffer.data(idx_kin_ih + 261);

    auto tk_xxyyzz_xyyyy = pbuffer.data(idx_kin_ih + 262);

    auto tk_xxyyzz_xyyyz = pbuffer.data(idx_kin_ih + 263);

    auto tk_xxyyzz_xyyzz = pbuffer.data(idx_kin_ih + 264);

    auto tk_xxyyzz_xyzzz = pbuffer.data(idx_kin_ih + 265);

    auto tk_xxyyzz_xzzzz = pbuffer.data(idx_kin_ih + 266);

    auto tk_xxyyzz_yyyyy = pbuffer.data(idx_kin_ih + 267);

    auto tk_xxyyzz_yyyyz = pbuffer.data(idx_kin_ih + 268);

    auto tk_xxyyzz_yyyzz = pbuffer.data(idx_kin_ih + 269);

    auto tk_xxyyzz_yyzzz = pbuffer.data(idx_kin_ih + 270);

    auto tk_xxyyzz_yzzzz = pbuffer.data(idx_kin_ih + 271);

    auto tk_xxyyzz_zzzzz = pbuffer.data(idx_kin_ih + 272);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tk_xxyy_xxxxy, tk_xxyy_xxxyy, tk_xxyy_xxyyy, tk_xxyy_xyyyy, tk_xxyyz_xxxxy, tk_xxyyz_xxxyy, tk_xxyyz_xxyyy, tk_xxyyz_xyyyy, tk_xxyyzz_xxxxx, tk_xxyyzz_xxxxy, tk_xxyyzz_xxxxz, tk_xxyyzz_xxxyy, tk_xxyyzz_xxxyz, tk_xxyyzz_xxxzz, tk_xxyyzz_xxyyy, tk_xxyyzz_xxyyz, tk_xxyyzz_xxyzz, tk_xxyyzz_xxzzz, tk_xxyyzz_xyyyy, tk_xxyyzz_xyyyz, tk_xxyyzz_xyyzz, tk_xxyyzz_xyzzz, tk_xxyyzz_xzzzz, tk_xxyyzz_yyyyy, tk_xxyyzz_yyyyz, tk_xxyyzz_yyyzz, tk_xxyyzz_yyzzz, tk_xxyyzz_yzzzz, tk_xxyyzz_zzzzz, tk_xxyzz_xxxxx, tk_xxyzz_xxxxz, tk_xxyzz_xxxzz, tk_xxyzz_xxzzz, tk_xxyzz_xzzzz, tk_xxzz_xxxxx, tk_xxzz_xxxxz, tk_xxzz_xxxzz, tk_xxzz_xxzzz, tk_xxzz_xzzzz, tk_xyyzz_xxxyz, tk_xyyzz_xxyyz, tk_xyyzz_xxyz, tk_xyyzz_xxyzz, tk_xyyzz_xyyyz, tk_xyyzz_xyyz, tk_xyyzz_xyyzz, tk_xyyzz_xyzz, tk_xyyzz_xyzzz, tk_xyyzz_yyyyy, tk_xyyzz_yyyyz, tk_xyyzz_yyyz, tk_xyyzz_yyyzz, tk_xyyzz_yyzz, tk_xyyzz_yyzzz, tk_xyyzz_yzzz, tk_xyyzz_yzzzz, tk_xyyzz_zzzzz, tk_yyzz_xxxyz, tk_yyzz_xxyyz, tk_yyzz_xxyzz, tk_yyzz_xyyyz, tk_yyzz_xyyzz, tk_yyzz_xyzzz, tk_yyzz_yyyyy, tk_yyzz_yyyyz, tk_yyzz_yyyzz, tk_yyzz_yyzzz, tk_yyzz_yzzzz, tk_yyzz_zzzzz, ts_xxyy_xxxxy, ts_xxyy_xxxyy, ts_xxyy_xxyyy, ts_xxyy_xyyyy, ts_xxyyzz_xxxxx, ts_xxyyzz_xxxxy, ts_xxyyzz_xxxxz, ts_xxyyzz_xxxyy, ts_xxyyzz_xxxyz, ts_xxyyzz_xxxzz, ts_xxyyzz_xxyyy, ts_xxyyzz_xxyyz, ts_xxyyzz_xxyzz, ts_xxyyzz_xxzzz, ts_xxyyzz_xyyyy, ts_xxyyzz_xyyyz, ts_xxyyzz_xyyzz, ts_xxyyzz_xyzzz, ts_xxyyzz_xzzzz, ts_xxyyzz_yyyyy, ts_xxyyzz_yyyyz, ts_xxyyzz_yyyzz, ts_xxyyzz_yyzzz, ts_xxyyzz_yzzzz, ts_xxyyzz_zzzzz, ts_xxzz_xxxxx, ts_xxzz_xxxxz, ts_xxzz_xxxzz, ts_xxzz_xxzzz, ts_xxzz_xzzzz, ts_yyzz_xxxyz, ts_yyzz_xxyyz, ts_yyzz_xxyzz, ts_yyzz_xyyyz, ts_yyzz_xyyzz, ts_yyzz_xyzzz, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyzz, ts_yyzz_yyzzz, ts_yyzz_yzzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyzz_xxxxx[i] = -2.0 * ts_xxzz_xxxxx[i] * fbe_0 * fz_0 + tk_xxzz_xxxxx[i] * fe_0 + tk_xxyzz_xxxxx[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxxx[i] * fz_0;

        tk_xxyyzz_xxxxy[i] = -2.0 * ts_xxyy_xxxxy[i] * fbe_0 * fz_0 + tk_xxyy_xxxxy[i] * fe_0 + tk_xxyyz_xxxxy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxxxy[i] * fz_0;

        tk_xxyyzz_xxxxz[i] = -2.0 * ts_xxzz_xxxxz[i] * fbe_0 * fz_0 + tk_xxzz_xxxxz[i] * fe_0 + tk_xxyzz_xxxxz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxxz[i] * fz_0;

        tk_xxyyzz_xxxyy[i] = -2.0 * ts_xxyy_xxxyy[i] * fbe_0 * fz_0 + tk_xxyy_xxxyy[i] * fe_0 + tk_xxyyz_xxxyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxxyy[i] * fz_0;

        tk_xxyyzz_xxxyz[i] = -2.0 * ts_yyzz_xxxyz[i] * fbe_0 * fz_0 + tk_yyzz_xxxyz[i] * fe_0 + 3.0 * tk_xyyzz_xxyz[i] * fe_0 + tk_xyyzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxxyz[i] * fz_0;

        tk_xxyyzz_xxxzz[i] = -2.0 * ts_xxzz_xxxzz[i] * fbe_0 * fz_0 + tk_xxzz_xxxzz[i] * fe_0 + tk_xxyzz_xxxzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxxzz[i] * fz_0;

        tk_xxyyzz_xxyyy[i] = -2.0 * ts_xxyy_xxyyy[i] * fbe_0 * fz_0 + tk_xxyy_xxyyy[i] * fe_0 + tk_xxyyz_xxyyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xxyyy[i] * fz_0;

        tk_xxyyzz_xxyyz[i] = -2.0 * ts_yyzz_xxyyz[i] * fbe_0 * fz_0 + tk_yyzz_xxyyz[i] * fe_0 + 2.0 * tk_xyyzz_xyyz[i] * fe_0 + tk_xyyzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxyyz[i] * fz_0;

        tk_xxyyzz_xxyzz[i] = -2.0 * ts_yyzz_xxyzz[i] * fbe_0 * fz_0 + tk_yyzz_xxyzz[i] * fe_0 + 2.0 * tk_xyyzz_xyzz[i] * fe_0 + tk_xyyzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xxyzz[i] * fz_0;

        tk_xxyyzz_xxzzz[i] = -2.0 * ts_xxzz_xxzzz[i] * fbe_0 * fz_0 + tk_xxzz_xxzzz[i] * fe_0 + tk_xxyzz_xxzzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xxzzz[i] * fz_0;

        tk_xxyyzz_xyyyy[i] = -2.0 * ts_xxyy_xyyyy[i] * fbe_0 * fz_0 + tk_xxyy_xyyyy[i] * fe_0 + tk_xxyyz_xyyyy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xyyyy[i] * fz_0;

        tk_xxyyzz_xyyyz[i] = -2.0 * ts_yyzz_xyyyz[i] * fbe_0 * fz_0 + tk_yyzz_xyyyz[i] * fe_0 + tk_xyyzz_yyyz[i] * fe_0 + tk_xyyzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xyyyz[i] * fz_0;

        tk_xxyyzz_xyyzz[i] = -2.0 * ts_yyzz_xyyzz[i] * fbe_0 * fz_0 + tk_yyzz_xyyzz[i] * fe_0 + tk_xyyzz_yyzz[i] * fe_0 + tk_xyyzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xyyzz[i] * fz_0;

        tk_xxyyzz_xyzzz[i] = -2.0 * ts_yyzz_xyzzz[i] * fbe_0 * fz_0 + tk_yyzz_xyzzz[i] * fe_0 + tk_xyyzz_yzzz[i] * fe_0 + tk_xyyzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_xyzzz[i] * fz_0;

        tk_xxyyzz_xzzzz[i] = -2.0 * ts_xxzz_xzzzz[i] * fbe_0 * fz_0 + tk_xxzz_xzzzz[i] * fe_0 + tk_xxyzz_xzzzz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xzzzz[i] * fz_0;

        tk_xxyyzz_yyyyy[i] = -2.0 * ts_yyzz_yyyyy[i] * fbe_0 * fz_0 + tk_yyzz_yyyyy[i] * fe_0 + tk_xyyzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyyy[i] * fz_0;

        tk_xxyyzz_yyyyz[i] = -2.0 * ts_yyzz_yyyyz[i] * fbe_0 * fz_0 + tk_yyzz_yyyyz[i] * fe_0 + tk_xyyzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyyz[i] * fz_0;

        tk_xxyyzz_yyyzz[i] = -2.0 * ts_yyzz_yyyzz[i] * fbe_0 * fz_0 + tk_yyzz_yyyzz[i] * fe_0 + tk_xyyzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyyzz[i] * fz_0;

        tk_xxyyzz_yyzzz[i] = -2.0 * ts_yyzz_yyzzz[i] * fbe_0 * fz_0 + tk_yyzz_yyzzz[i] * fe_0 + tk_xyyzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yyzzz[i] * fz_0;

        tk_xxyyzz_yzzzz[i] = -2.0 * ts_yyzz_yzzzz[i] * fbe_0 * fz_0 + tk_yyzz_yzzzz[i] * fe_0 + tk_xyyzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yzzzz[i] * fz_0;

        tk_xxyyzz_zzzzz[i] = -2.0 * ts_yyzz_zzzzz[i] * fbe_0 * fz_0 + tk_yyzz_zzzzz[i] * fe_0 + tk_xyyzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xxyyzz_zzzzz[i] * fz_0;
    }

    // Set up 273-294 components of targeted buffer : IH

    auto tk_xxyzzz_xxxxx = pbuffer.data(idx_kin_ih + 273);

    auto tk_xxyzzz_xxxxy = pbuffer.data(idx_kin_ih + 274);

    auto tk_xxyzzz_xxxxz = pbuffer.data(idx_kin_ih + 275);

    auto tk_xxyzzz_xxxyy = pbuffer.data(idx_kin_ih + 276);

    auto tk_xxyzzz_xxxyz = pbuffer.data(idx_kin_ih + 277);

    auto tk_xxyzzz_xxxzz = pbuffer.data(idx_kin_ih + 278);

    auto tk_xxyzzz_xxyyy = pbuffer.data(idx_kin_ih + 279);

    auto tk_xxyzzz_xxyyz = pbuffer.data(idx_kin_ih + 280);

    auto tk_xxyzzz_xxyzz = pbuffer.data(idx_kin_ih + 281);

    auto tk_xxyzzz_xxzzz = pbuffer.data(idx_kin_ih + 282);

    auto tk_xxyzzz_xyyyy = pbuffer.data(idx_kin_ih + 283);

    auto tk_xxyzzz_xyyyz = pbuffer.data(idx_kin_ih + 284);

    auto tk_xxyzzz_xyyzz = pbuffer.data(idx_kin_ih + 285);

    auto tk_xxyzzz_xyzzz = pbuffer.data(idx_kin_ih + 286);

    auto tk_xxyzzz_xzzzz = pbuffer.data(idx_kin_ih + 287);

    auto tk_xxyzzz_yyyyy = pbuffer.data(idx_kin_ih + 288);

    auto tk_xxyzzz_yyyyz = pbuffer.data(idx_kin_ih + 289);

    auto tk_xxyzzz_yyyzz = pbuffer.data(idx_kin_ih + 290);

    auto tk_xxyzzz_yyzzz = pbuffer.data(idx_kin_ih + 291);

    auto tk_xxyzzz_yzzzz = pbuffer.data(idx_kin_ih + 292);

    auto tk_xxyzzz_zzzzz = pbuffer.data(idx_kin_ih + 293);

    #pragma omp simd aligned(pa_y, tk_xxyzzz_xxxxx, tk_xxyzzz_xxxxy, tk_xxyzzz_xxxxz, tk_xxyzzz_xxxyy, tk_xxyzzz_xxxyz, tk_xxyzzz_xxxzz, tk_xxyzzz_xxyyy, tk_xxyzzz_xxyyz, tk_xxyzzz_xxyzz, tk_xxyzzz_xxzzz, tk_xxyzzz_xyyyy, tk_xxyzzz_xyyyz, tk_xxyzzz_xyyzz, tk_xxyzzz_xyzzz, tk_xxyzzz_xzzzz, tk_xxyzzz_yyyyy, tk_xxyzzz_yyyyz, tk_xxyzzz_yyyzz, tk_xxyzzz_yyzzz, tk_xxyzzz_yzzzz, tk_xxyzzz_zzzzz, tk_xxzzz_xxxx, tk_xxzzz_xxxxx, tk_xxzzz_xxxxy, tk_xxzzz_xxxxz, tk_xxzzz_xxxy, tk_xxzzz_xxxyy, tk_xxzzz_xxxyz, tk_xxzzz_xxxz, tk_xxzzz_xxxzz, tk_xxzzz_xxyy, tk_xxzzz_xxyyy, tk_xxzzz_xxyyz, tk_xxzzz_xxyz, tk_xxzzz_xxyzz, tk_xxzzz_xxzz, tk_xxzzz_xxzzz, tk_xxzzz_xyyy, tk_xxzzz_xyyyy, tk_xxzzz_xyyyz, tk_xxzzz_xyyz, tk_xxzzz_xyyzz, tk_xxzzz_xyzz, tk_xxzzz_xyzzz, tk_xxzzz_xzzz, tk_xxzzz_xzzzz, tk_xxzzz_yyyy, tk_xxzzz_yyyyy, tk_xxzzz_yyyyz, tk_xxzzz_yyyz, tk_xxzzz_yyyzz, tk_xxzzz_yyzz, tk_xxzzz_yyzzz, tk_xxzzz_yzzz, tk_xxzzz_yzzzz, tk_xxzzz_zzzz, tk_xxzzz_zzzzz, ts_xxyzzz_xxxxx, ts_xxyzzz_xxxxy, ts_xxyzzz_xxxxz, ts_xxyzzz_xxxyy, ts_xxyzzz_xxxyz, ts_xxyzzz_xxxzz, ts_xxyzzz_xxyyy, ts_xxyzzz_xxyyz, ts_xxyzzz_xxyzz, ts_xxyzzz_xxzzz, ts_xxyzzz_xyyyy, ts_xxyzzz_xyyyz, ts_xxyzzz_xyyzz, ts_xxyzzz_xyzzz, ts_xxyzzz_xzzzz, ts_xxyzzz_yyyyy, ts_xxyzzz_yyyyz, ts_xxyzzz_yyyzz, ts_xxyzzz_yyzzz, ts_xxyzzz_yzzzz, ts_xxyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzzz_xxxxx[i] = tk_xxzzz_xxxxx[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxx[i] * fz_0;

        tk_xxyzzz_xxxxy[i] = tk_xxzzz_xxxx[i] * fe_0 + tk_xxzzz_xxxxy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxy[i] * fz_0;

        tk_xxyzzz_xxxxz[i] = tk_xxzzz_xxxxz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxxz[i] * fz_0;

        tk_xxyzzz_xxxyy[i] = 2.0 * tk_xxzzz_xxxy[i] * fe_0 + tk_xxzzz_xxxyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxyy[i] * fz_0;

        tk_xxyzzz_xxxyz[i] = tk_xxzzz_xxxz[i] * fe_0 + tk_xxzzz_xxxyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxyz[i] * fz_0;

        tk_xxyzzz_xxxzz[i] = tk_xxzzz_xxxzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxxzz[i] * fz_0;

        tk_xxyzzz_xxyyy[i] = 3.0 * tk_xxzzz_xxyy[i] * fe_0 + tk_xxzzz_xxyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyyy[i] * fz_0;

        tk_xxyzzz_xxyyz[i] = 2.0 * tk_xxzzz_xxyz[i] * fe_0 + tk_xxzzz_xxyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyyz[i] * fz_0;

        tk_xxyzzz_xxyzz[i] = tk_xxzzz_xxzz[i] * fe_0 + tk_xxzzz_xxyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxyzz[i] * fz_0;

        tk_xxyzzz_xxzzz[i] = tk_xxzzz_xxzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xxzzz[i] * fz_0;

        tk_xxyzzz_xyyyy[i] = 4.0 * tk_xxzzz_xyyy[i] * fe_0 + tk_xxzzz_xyyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyyy[i] * fz_0;

        tk_xxyzzz_xyyyz[i] = 3.0 * tk_xxzzz_xyyz[i] * fe_0 + tk_xxzzz_xyyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyyz[i] * fz_0;

        tk_xxyzzz_xyyzz[i] = 2.0 * tk_xxzzz_xyzz[i] * fe_0 + tk_xxzzz_xyyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyyzz[i] * fz_0;

        tk_xxyzzz_xyzzz[i] = tk_xxzzz_xzzz[i] * fe_0 + tk_xxzzz_xyzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xyzzz[i] * fz_0;

        tk_xxyzzz_xzzzz[i] = tk_xxzzz_xzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xzzzz[i] * fz_0;

        tk_xxyzzz_yyyyy[i] = 5.0 * tk_xxzzz_yyyy[i] * fe_0 + tk_xxzzz_yyyyy[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyyy[i] * fz_0;

        tk_xxyzzz_yyyyz[i] = 4.0 * tk_xxzzz_yyyz[i] * fe_0 + tk_xxzzz_yyyyz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyyz[i] * fz_0;

        tk_xxyzzz_yyyzz[i] = 3.0 * tk_xxzzz_yyzz[i] * fe_0 + tk_xxzzz_yyyzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyyzz[i] * fz_0;

        tk_xxyzzz_yyzzz[i] = 2.0 * tk_xxzzz_yzzz[i] * fe_0 + tk_xxzzz_yyzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yyzzz[i] * fz_0;

        tk_xxyzzz_yzzzz[i] = tk_xxzzz_zzzz[i] * fe_0 + tk_xxzzz_yzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yzzzz[i] * fz_0;

        tk_xxyzzz_zzzzz[i] = tk_xxzzz_zzzzz[i] * pa_y[i] + 2.0 * ts_xxyzzz_zzzzz[i] * fz_0;
    }

    // Set up 294-315 components of targeted buffer : IH

    auto tk_xxzzzz_xxxxx = pbuffer.data(idx_kin_ih + 294);

    auto tk_xxzzzz_xxxxy = pbuffer.data(idx_kin_ih + 295);

    auto tk_xxzzzz_xxxxz = pbuffer.data(idx_kin_ih + 296);

    auto tk_xxzzzz_xxxyy = pbuffer.data(idx_kin_ih + 297);

    auto tk_xxzzzz_xxxyz = pbuffer.data(idx_kin_ih + 298);

    auto tk_xxzzzz_xxxzz = pbuffer.data(idx_kin_ih + 299);

    auto tk_xxzzzz_xxyyy = pbuffer.data(idx_kin_ih + 300);

    auto tk_xxzzzz_xxyyz = pbuffer.data(idx_kin_ih + 301);

    auto tk_xxzzzz_xxyzz = pbuffer.data(idx_kin_ih + 302);

    auto tk_xxzzzz_xxzzz = pbuffer.data(idx_kin_ih + 303);

    auto tk_xxzzzz_xyyyy = pbuffer.data(idx_kin_ih + 304);

    auto tk_xxzzzz_xyyyz = pbuffer.data(idx_kin_ih + 305);

    auto tk_xxzzzz_xyyzz = pbuffer.data(idx_kin_ih + 306);

    auto tk_xxzzzz_xyzzz = pbuffer.data(idx_kin_ih + 307);

    auto tk_xxzzzz_xzzzz = pbuffer.data(idx_kin_ih + 308);

    auto tk_xxzzzz_yyyyy = pbuffer.data(idx_kin_ih + 309);

    auto tk_xxzzzz_yyyyz = pbuffer.data(idx_kin_ih + 310);

    auto tk_xxzzzz_yyyzz = pbuffer.data(idx_kin_ih + 311);

    auto tk_xxzzzz_yyzzz = pbuffer.data(idx_kin_ih + 312);

    auto tk_xxzzzz_yzzzz = pbuffer.data(idx_kin_ih + 313);

    auto tk_xxzzzz_zzzzz = pbuffer.data(idx_kin_ih + 314);

    #pragma omp simd aligned(pa_x, pa_z, tk_xxzz_xxxxx, tk_xxzz_xxxxy, tk_xxzz_xxxyy, tk_xxzz_xxyyy, tk_xxzz_xyyyy, tk_xxzzz_xxxxx, tk_xxzzz_xxxxy, tk_xxzzz_xxxyy, tk_xxzzz_xxyyy, tk_xxzzz_xyyyy, tk_xxzzzz_xxxxx, tk_xxzzzz_xxxxy, tk_xxzzzz_xxxxz, tk_xxzzzz_xxxyy, tk_xxzzzz_xxxyz, tk_xxzzzz_xxxzz, tk_xxzzzz_xxyyy, tk_xxzzzz_xxyyz, tk_xxzzzz_xxyzz, tk_xxzzzz_xxzzz, tk_xxzzzz_xyyyy, tk_xxzzzz_xyyyz, tk_xxzzzz_xyyzz, tk_xxzzzz_xyzzz, tk_xxzzzz_xzzzz, tk_xxzzzz_yyyyy, tk_xxzzzz_yyyyz, tk_xxzzzz_yyyzz, tk_xxzzzz_yyzzz, tk_xxzzzz_yzzzz, tk_xxzzzz_zzzzz, tk_xzzzz_xxxxz, tk_xzzzz_xxxyz, tk_xzzzz_xxxz, tk_xzzzz_xxxzz, tk_xzzzz_xxyyz, tk_xzzzz_xxyz, tk_xzzzz_xxyzz, tk_xzzzz_xxzz, tk_xzzzz_xxzzz, tk_xzzzz_xyyyz, tk_xzzzz_xyyz, tk_xzzzz_xyyzz, tk_xzzzz_xyzz, tk_xzzzz_xyzzz, tk_xzzzz_xzzz, tk_xzzzz_xzzzz, tk_xzzzz_yyyyy, tk_xzzzz_yyyyz, tk_xzzzz_yyyz, tk_xzzzz_yyyzz, tk_xzzzz_yyzz, tk_xzzzz_yyzzz, tk_xzzzz_yzzz, tk_xzzzz_yzzzz, tk_xzzzz_zzzz, tk_xzzzz_zzzzz, tk_zzzz_xxxxz, tk_zzzz_xxxyz, tk_zzzz_xxxzz, tk_zzzz_xxyyz, tk_zzzz_xxyzz, tk_zzzz_xxzzz, tk_zzzz_xyyyz, tk_zzzz_xyyzz, tk_zzzz_xyzzz, tk_zzzz_xzzzz, tk_zzzz_yyyyy, tk_zzzz_yyyyz, tk_zzzz_yyyzz, tk_zzzz_yyzzz, tk_zzzz_yzzzz, tk_zzzz_zzzzz, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxyy, ts_xxzz_xxyyy, ts_xxzz_xyyyy, ts_xxzzzz_xxxxx, ts_xxzzzz_xxxxy, ts_xxzzzz_xxxxz, ts_xxzzzz_xxxyy, ts_xxzzzz_xxxyz, ts_xxzzzz_xxxzz, ts_xxzzzz_xxyyy, ts_xxzzzz_xxyyz, ts_xxzzzz_xxyzz, ts_xxzzzz_xxzzz, ts_xxzzzz_xyyyy, ts_xxzzzz_xyyyz, ts_xxzzzz_xyyzz, ts_xxzzzz_xyzzz, ts_xxzzzz_xzzzz, ts_xxzzzz_yyyyy, ts_xxzzzz_yyyyz, ts_xxzzzz_yyyzz, ts_xxzzzz_yyzzz, ts_xxzzzz_yzzzz, ts_xxzzzz_zzzzz, ts_zzzz_xxxxz, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzzz_xxxxx[i] = -6.0 * ts_xxzz_xxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxx[i] * fe_0 + tk_xxzzz_xxxxx[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxxxx[i] * fz_0;

        tk_xxzzzz_xxxxy[i] = -6.0 * ts_xxzz_xxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxxy[i] * fe_0 + tk_xxzzz_xxxxy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxxxy[i] * fz_0;

        tk_xxzzzz_xxxxz[i] = -2.0 * ts_zzzz_xxxxz[i] * fbe_0 * fz_0 + tk_zzzz_xxxxz[i] * fe_0 + 4.0 * tk_xzzzz_xxxz[i] * fe_0 + tk_xzzzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxxz[i] * fz_0;

        tk_xxzzzz_xxxyy[i] = -6.0 * ts_xxzz_xxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxxyy[i] * fe_0 + tk_xxzzz_xxxyy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxxyy[i] * fz_0;

        tk_xxzzzz_xxxyz[i] = -2.0 * ts_zzzz_xxxyz[i] * fbe_0 * fz_0 + tk_zzzz_xxxyz[i] * fe_0 + 3.0 * tk_xzzzz_xxyz[i] * fe_0 + tk_xzzzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxyz[i] * fz_0;

        tk_xxzzzz_xxxzz[i] = -2.0 * ts_zzzz_xxxzz[i] * fbe_0 * fz_0 + tk_zzzz_xxxzz[i] * fe_0 + 3.0 * tk_xzzzz_xxzz[i] * fe_0 + tk_xzzzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxxzz[i] * fz_0;

        tk_xxzzzz_xxyyy[i] = -6.0 * ts_xxzz_xxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xxyyy[i] * fe_0 + tk_xxzzz_xxyyy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xxyyy[i] * fz_0;

        tk_xxzzzz_xxyyz[i] = -2.0 * ts_zzzz_xxyyz[i] * fbe_0 * fz_0 + tk_zzzz_xxyyz[i] * fe_0 + 2.0 * tk_xzzzz_xyyz[i] * fe_0 + tk_xzzzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxyyz[i] * fz_0;

        tk_xxzzzz_xxyzz[i] = -2.0 * ts_zzzz_xxyzz[i] * fbe_0 * fz_0 + tk_zzzz_xxyzz[i] * fe_0 + 2.0 * tk_xzzzz_xyzz[i] * fe_0 + tk_xzzzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxyzz[i] * fz_0;

        tk_xxzzzz_xxzzz[i] = -2.0 * ts_zzzz_xxzzz[i] * fbe_0 * fz_0 + tk_zzzz_xxzzz[i] * fe_0 + 2.0 * tk_xzzzz_xzzz[i] * fe_0 + tk_xzzzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xxzzz[i] * fz_0;

        tk_xxzzzz_xyyyy[i] = -6.0 * ts_xxzz_xyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xyyyy[i] * fe_0 + tk_xxzzz_xyyyy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xyyyy[i] * fz_0;

        tk_xxzzzz_xyyyz[i] = -2.0 * ts_zzzz_xyyyz[i] * fbe_0 * fz_0 + tk_zzzz_xyyyz[i] * fe_0 + tk_xzzzz_yyyz[i] * fe_0 + tk_xzzzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xyyyz[i] * fz_0;

        tk_xxzzzz_xyyzz[i] = -2.0 * ts_zzzz_xyyzz[i] * fbe_0 * fz_0 + tk_zzzz_xyyzz[i] * fe_0 + tk_xzzzz_yyzz[i] * fe_0 + tk_xzzzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xyyzz[i] * fz_0;

        tk_xxzzzz_xyzzz[i] = -2.0 * ts_zzzz_xyzzz[i] * fbe_0 * fz_0 + tk_zzzz_xyzzz[i] * fe_0 + tk_xzzzz_yzzz[i] * fe_0 + tk_xzzzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xyzzz[i] * fz_0;

        tk_xxzzzz_xzzzz[i] = -2.0 * ts_zzzz_xzzzz[i] * fbe_0 * fz_0 + tk_zzzz_xzzzz[i] * fe_0 + tk_xzzzz_zzzz[i] * fe_0 + tk_xzzzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_xzzzz[i] * fz_0;

        tk_xxzzzz_yyyyy[i] = -2.0 * ts_zzzz_yyyyy[i] * fbe_0 * fz_0 + tk_zzzz_yyyyy[i] * fe_0 + tk_xzzzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyyy[i] * fz_0;

        tk_xxzzzz_yyyyz[i] = -2.0 * ts_zzzz_yyyyz[i] * fbe_0 * fz_0 + tk_zzzz_yyyyz[i] * fe_0 + tk_xzzzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyyz[i] * fz_0;

        tk_xxzzzz_yyyzz[i] = -2.0 * ts_zzzz_yyyzz[i] * fbe_0 * fz_0 + tk_zzzz_yyyzz[i] * fe_0 + tk_xzzzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyyzz[i] * fz_0;

        tk_xxzzzz_yyzzz[i] = -2.0 * ts_zzzz_yyzzz[i] * fbe_0 * fz_0 + tk_zzzz_yyzzz[i] * fe_0 + tk_xzzzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yyzzz[i] * fz_0;

        tk_xxzzzz_yzzzz[i] = -2.0 * ts_zzzz_yzzzz[i] * fbe_0 * fz_0 + tk_zzzz_yzzzz[i] * fe_0 + tk_xzzzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yzzzz[i] * fz_0;

        tk_xxzzzz_zzzzz[i] = -2.0 * ts_zzzz_zzzzz[i] * fbe_0 * fz_0 + tk_zzzz_zzzzz[i] * fe_0 + tk_xzzzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xxzzzz_zzzzz[i] * fz_0;
    }

    // Set up 315-336 components of targeted buffer : IH

    auto tk_xyyyyy_xxxxx = pbuffer.data(idx_kin_ih + 315);

    auto tk_xyyyyy_xxxxy = pbuffer.data(idx_kin_ih + 316);

    auto tk_xyyyyy_xxxxz = pbuffer.data(idx_kin_ih + 317);

    auto tk_xyyyyy_xxxyy = pbuffer.data(idx_kin_ih + 318);

    auto tk_xyyyyy_xxxyz = pbuffer.data(idx_kin_ih + 319);

    auto tk_xyyyyy_xxxzz = pbuffer.data(idx_kin_ih + 320);

    auto tk_xyyyyy_xxyyy = pbuffer.data(idx_kin_ih + 321);

    auto tk_xyyyyy_xxyyz = pbuffer.data(idx_kin_ih + 322);

    auto tk_xyyyyy_xxyzz = pbuffer.data(idx_kin_ih + 323);

    auto tk_xyyyyy_xxzzz = pbuffer.data(idx_kin_ih + 324);

    auto tk_xyyyyy_xyyyy = pbuffer.data(idx_kin_ih + 325);

    auto tk_xyyyyy_xyyyz = pbuffer.data(idx_kin_ih + 326);

    auto tk_xyyyyy_xyyzz = pbuffer.data(idx_kin_ih + 327);

    auto tk_xyyyyy_xyzzz = pbuffer.data(idx_kin_ih + 328);

    auto tk_xyyyyy_xzzzz = pbuffer.data(idx_kin_ih + 329);

    auto tk_xyyyyy_yyyyy = pbuffer.data(idx_kin_ih + 330);

    auto tk_xyyyyy_yyyyz = pbuffer.data(idx_kin_ih + 331);

    auto tk_xyyyyy_yyyzz = pbuffer.data(idx_kin_ih + 332);

    auto tk_xyyyyy_yyzzz = pbuffer.data(idx_kin_ih + 333);

    auto tk_xyyyyy_yzzzz = pbuffer.data(idx_kin_ih + 334);

    auto tk_xyyyyy_zzzzz = pbuffer.data(idx_kin_ih + 335);

    #pragma omp simd aligned(pa_x, tk_xyyyyy_xxxxx, tk_xyyyyy_xxxxy, tk_xyyyyy_xxxxz, tk_xyyyyy_xxxyy, tk_xyyyyy_xxxyz, tk_xyyyyy_xxxzz, tk_xyyyyy_xxyyy, tk_xyyyyy_xxyyz, tk_xyyyyy_xxyzz, tk_xyyyyy_xxzzz, tk_xyyyyy_xyyyy, tk_xyyyyy_xyyyz, tk_xyyyyy_xyyzz, tk_xyyyyy_xyzzz, tk_xyyyyy_xzzzz, tk_xyyyyy_yyyyy, tk_xyyyyy_yyyyz, tk_xyyyyy_yyyzz, tk_xyyyyy_yyzzz, tk_xyyyyy_yzzzz, tk_xyyyyy_zzzzz, tk_yyyyy_xxxx, tk_yyyyy_xxxxx, tk_yyyyy_xxxxy, tk_yyyyy_xxxxz, tk_yyyyy_xxxy, tk_yyyyy_xxxyy, tk_yyyyy_xxxyz, tk_yyyyy_xxxz, tk_yyyyy_xxxzz, tk_yyyyy_xxyy, tk_yyyyy_xxyyy, tk_yyyyy_xxyyz, tk_yyyyy_xxyz, tk_yyyyy_xxyzz, tk_yyyyy_xxzz, tk_yyyyy_xxzzz, tk_yyyyy_xyyy, tk_yyyyy_xyyyy, tk_yyyyy_xyyyz, tk_yyyyy_xyyz, tk_yyyyy_xyyzz, tk_yyyyy_xyzz, tk_yyyyy_xyzzz, tk_yyyyy_xzzz, tk_yyyyy_xzzzz, tk_yyyyy_yyyy, tk_yyyyy_yyyyy, tk_yyyyy_yyyyz, tk_yyyyy_yyyz, tk_yyyyy_yyyzz, tk_yyyyy_yyzz, tk_yyyyy_yyzzz, tk_yyyyy_yzzz, tk_yyyyy_yzzzz, tk_yyyyy_zzzz, tk_yyyyy_zzzzz, ts_xyyyyy_xxxxx, ts_xyyyyy_xxxxy, ts_xyyyyy_xxxxz, ts_xyyyyy_xxxyy, ts_xyyyyy_xxxyz, ts_xyyyyy_xxxzz, ts_xyyyyy_xxyyy, ts_xyyyyy_xxyyz, ts_xyyyyy_xxyzz, ts_xyyyyy_xxzzz, ts_xyyyyy_xyyyy, ts_xyyyyy_xyyyz, ts_xyyyyy_xyyzz, ts_xyyyyy_xyzzz, ts_xyyyyy_xzzzz, ts_xyyyyy_yyyyy, ts_xyyyyy_yyyyz, ts_xyyyyy_yyyzz, ts_xyyyyy_yyzzz, ts_xyyyyy_yzzzz, ts_xyyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyy_xxxxx[i] = 5.0 * tk_yyyyy_xxxx[i] * fe_0 + tk_yyyyy_xxxxx[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxx[i] * fz_0;

        tk_xyyyyy_xxxxy[i] = 4.0 * tk_yyyyy_xxxy[i] * fe_0 + tk_yyyyy_xxxxy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxy[i] * fz_0;

        tk_xyyyyy_xxxxz[i] = 4.0 * tk_yyyyy_xxxz[i] * fe_0 + tk_yyyyy_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxxz[i] * fz_0;

        tk_xyyyyy_xxxyy[i] = 3.0 * tk_yyyyy_xxyy[i] * fe_0 + tk_yyyyy_xxxyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxyy[i] * fz_0;

        tk_xyyyyy_xxxyz[i] = 3.0 * tk_yyyyy_xxyz[i] * fe_0 + tk_yyyyy_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxyz[i] * fz_0;

        tk_xyyyyy_xxxzz[i] = 3.0 * tk_yyyyy_xxzz[i] * fe_0 + tk_yyyyy_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxxzz[i] * fz_0;

        tk_xyyyyy_xxyyy[i] = 2.0 * tk_yyyyy_xyyy[i] * fe_0 + tk_yyyyy_xxyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyyy[i] * fz_0;

        tk_xyyyyy_xxyyz[i] = 2.0 * tk_yyyyy_xyyz[i] * fe_0 + tk_yyyyy_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyyz[i] * fz_0;

        tk_xyyyyy_xxyzz[i] = 2.0 * tk_yyyyy_xyzz[i] * fe_0 + tk_yyyyy_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxyzz[i] * fz_0;

        tk_xyyyyy_xxzzz[i] = 2.0 * tk_yyyyy_xzzz[i] * fe_0 + tk_yyyyy_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xxzzz[i] * fz_0;

        tk_xyyyyy_xyyyy[i] = tk_yyyyy_yyyy[i] * fe_0 + tk_yyyyy_xyyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyyy[i] * fz_0;

        tk_xyyyyy_xyyyz[i] = tk_yyyyy_yyyz[i] * fe_0 + tk_yyyyy_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyyz[i] * fz_0;

        tk_xyyyyy_xyyzz[i] = tk_yyyyy_yyzz[i] * fe_0 + tk_yyyyy_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyyzz[i] * fz_0;

        tk_xyyyyy_xyzzz[i] = tk_yyyyy_yzzz[i] * fe_0 + tk_yyyyy_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xyzzz[i] * fz_0;

        tk_xyyyyy_xzzzz[i] = tk_yyyyy_zzzz[i] * fe_0 + tk_yyyyy_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xzzzz[i] * fz_0;

        tk_xyyyyy_yyyyy[i] = tk_yyyyy_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyyy[i] * fz_0;

        tk_xyyyyy_yyyyz[i] = tk_yyyyy_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyyz[i] * fz_0;

        tk_xyyyyy_yyyzz[i] = tk_yyyyy_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyyzz[i] * fz_0;

        tk_xyyyyy_yyzzz[i] = tk_yyyyy_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yyzzz[i] * fz_0;

        tk_xyyyyy_yzzzz[i] = tk_yyyyy_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yzzzz[i] * fz_0;

        tk_xyyyyy_zzzzz[i] = tk_yyyyy_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyy_zzzzz[i] * fz_0;
    }

    // Set up 336-357 components of targeted buffer : IH

    auto tk_xyyyyz_xxxxx = pbuffer.data(idx_kin_ih + 336);

    auto tk_xyyyyz_xxxxy = pbuffer.data(idx_kin_ih + 337);

    auto tk_xyyyyz_xxxxz = pbuffer.data(idx_kin_ih + 338);

    auto tk_xyyyyz_xxxyy = pbuffer.data(idx_kin_ih + 339);

    auto tk_xyyyyz_xxxyz = pbuffer.data(idx_kin_ih + 340);

    auto tk_xyyyyz_xxxzz = pbuffer.data(idx_kin_ih + 341);

    auto tk_xyyyyz_xxyyy = pbuffer.data(idx_kin_ih + 342);

    auto tk_xyyyyz_xxyyz = pbuffer.data(idx_kin_ih + 343);

    auto tk_xyyyyz_xxyzz = pbuffer.data(idx_kin_ih + 344);

    auto tk_xyyyyz_xxzzz = pbuffer.data(idx_kin_ih + 345);

    auto tk_xyyyyz_xyyyy = pbuffer.data(idx_kin_ih + 346);

    auto tk_xyyyyz_xyyyz = pbuffer.data(idx_kin_ih + 347);

    auto tk_xyyyyz_xyyzz = pbuffer.data(idx_kin_ih + 348);

    auto tk_xyyyyz_xyzzz = pbuffer.data(idx_kin_ih + 349);

    auto tk_xyyyyz_xzzzz = pbuffer.data(idx_kin_ih + 350);

    auto tk_xyyyyz_yyyyy = pbuffer.data(idx_kin_ih + 351);

    auto tk_xyyyyz_yyyyz = pbuffer.data(idx_kin_ih + 352);

    auto tk_xyyyyz_yyyzz = pbuffer.data(idx_kin_ih + 353);

    auto tk_xyyyyz_yyzzz = pbuffer.data(idx_kin_ih + 354);

    auto tk_xyyyyz_yzzzz = pbuffer.data(idx_kin_ih + 355);

    auto tk_xyyyyz_zzzzz = pbuffer.data(idx_kin_ih + 356);

    #pragma omp simd aligned(pa_x, pa_z, tk_xyyyy_xxxxx, tk_xyyyy_xxxxy, tk_xyyyy_xxxyy, tk_xyyyy_xxyyy, tk_xyyyy_xyyyy, tk_xyyyyz_xxxxx, tk_xyyyyz_xxxxy, tk_xyyyyz_xxxxz, tk_xyyyyz_xxxyy, tk_xyyyyz_xxxyz, tk_xyyyyz_xxxzz, tk_xyyyyz_xxyyy, tk_xyyyyz_xxyyz, tk_xyyyyz_xxyzz, tk_xyyyyz_xxzzz, tk_xyyyyz_xyyyy, tk_xyyyyz_xyyyz, tk_xyyyyz_xyyzz, tk_xyyyyz_xyzzz, tk_xyyyyz_xzzzz, tk_xyyyyz_yyyyy, tk_xyyyyz_yyyyz, tk_xyyyyz_yyyzz, tk_xyyyyz_yyzzz, tk_xyyyyz_yzzzz, tk_xyyyyz_zzzzz, tk_yyyyz_xxxxz, tk_yyyyz_xxxyz, tk_yyyyz_xxxz, tk_yyyyz_xxxzz, tk_yyyyz_xxyyz, tk_yyyyz_xxyz, tk_yyyyz_xxyzz, tk_yyyyz_xxzz, tk_yyyyz_xxzzz, tk_yyyyz_xyyyz, tk_yyyyz_xyyz, tk_yyyyz_xyyzz, tk_yyyyz_xyzz, tk_yyyyz_xyzzz, tk_yyyyz_xzzz, tk_yyyyz_xzzzz, tk_yyyyz_yyyyy, tk_yyyyz_yyyyz, tk_yyyyz_yyyz, tk_yyyyz_yyyzz, tk_yyyyz_yyzz, tk_yyyyz_yyzzz, tk_yyyyz_yzzz, tk_yyyyz_yzzzz, tk_yyyyz_zzzz, tk_yyyyz_zzzzz, ts_xyyyyz_xxxxx, ts_xyyyyz_xxxxy, ts_xyyyyz_xxxxz, ts_xyyyyz_xxxyy, ts_xyyyyz_xxxyz, ts_xyyyyz_xxxzz, ts_xyyyyz_xxyyy, ts_xyyyyz_xxyyz, ts_xyyyyz_xxyzz, ts_xyyyyz_xxzzz, ts_xyyyyz_xyyyy, ts_xyyyyz_xyyyz, ts_xyyyyz_xyyzz, ts_xyyyyz_xyzzz, ts_xyyyyz_xzzzz, ts_xyyyyz_yyyyy, ts_xyyyyz_yyyyz, ts_xyyyyz_yyyzz, ts_xyyyyz_yyzzz, ts_xyyyyz_yzzzz, ts_xyyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyz_xxxxx[i] = tk_xyyyy_xxxxx[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxxx[i] * fz_0;

        tk_xyyyyz_xxxxy[i] = tk_xyyyy_xxxxy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxxy[i] * fz_0;

        tk_xyyyyz_xxxxz[i] = 4.0 * tk_yyyyz_xxxz[i] * fe_0 + tk_yyyyz_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxxz[i] * fz_0;

        tk_xyyyyz_xxxyy[i] = tk_xyyyy_xxxyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxxyy[i] * fz_0;

        tk_xyyyyz_xxxyz[i] = 3.0 * tk_yyyyz_xxyz[i] * fe_0 + tk_yyyyz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxyz[i] * fz_0;

        tk_xyyyyz_xxxzz[i] = 3.0 * tk_yyyyz_xxzz[i] * fe_0 + tk_yyyyz_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxxzz[i] * fz_0;

        tk_xyyyyz_xxyyy[i] = tk_xyyyy_xxyyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xxyyy[i] * fz_0;

        tk_xyyyyz_xxyyz[i] = 2.0 * tk_yyyyz_xyyz[i] * fe_0 + tk_yyyyz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxyyz[i] * fz_0;

        tk_xyyyyz_xxyzz[i] = 2.0 * tk_yyyyz_xyzz[i] * fe_0 + tk_yyyyz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxyzz[i] * fz_0;

        tk_xyyyyz_xxzzz[i] = 2.0 * tk_yyyyz_xzzz[i] * fe_0 + tk_yyyyz_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xxzzz[i] * fz_0;

        tk_xyyyyz_xyyyy[i] = tk_xyyyy_xyyyy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xyyyy[i] * fz_0;

        tk_xyyyyz_xyyyz[i] = tk_yyyyz_yyyz[i] * fe_0 + tk_yyyyz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyyyz[i] * fz_0;

        tk_xyyyyz_xyyzz[i] = tk_yyyyz_yyzz[i] * fe_0 + tk_yyyyz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyyzz[i] * fz_0;

        tk_xyyyyz_xyzzz[i] = tk_yyyyz_yzzz[i] * fe_0 + tk_yyyyz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xyzzz[i] * fz_0;

        tk_xyyyyz_xzzzz[i] = tk_yyyyz_zzzz[i] * fe_0 + tk_yyyyz_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xzzzz[i] * fz_0;

        tk_xyyyyz_yyyyy[i] = tk_yyyyz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyyy[i] * fz_0;

        tk_xyyyyz_yyyyz[i] = tk_yyyyz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyyz[i] * fz_0;

        tk_xyyyyz_yyyzz[i] = tk_yyyyz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyyzz[i] * fz_0;

        tk_xyyyyz_yyzzz[i] = tk_yyyyz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yyzzz[i] * fz_0;

        tk_xyyyyz_yzzzz[i] = tk_yyyyz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yzzzz[i] * fz_0;

        tk_xyyyyz_zzzzz[i] = tk_yyyyz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyyyz_zzzzz[i] * fz_0;
    }

    // Set up 357-378 components of targeted buffer : IH

    auto tk_xyyyzz_xxxxx = pbuffer.data(idx_kin_ih + 357);

    auto tk_xyyyzz_xxxxy = pbuffer.data(idx_kin_ih + 358);

    auto tk_xyyyzz_xxxxz = pbuffer.data(idx_kin_ih + 359);

    auto tk_xyyyzz_xxxyy = pbuffer.data(idx_kin_ih + 360);

    auto tk_xyyyzz_xxxyz = pbuffer.data(idx_kin_ih + 361);

    auto tk_xyyyzz_xxxzz = pbuffer.data(idx_kin_ih + 362);

    auto tk_xyyyzz_xxyyy = pbuffer.data(idx_kin_ih + 363);

    auto tk_xyyyzz_xxyyz = pbuffer.data(idx_kin_ih + 364);

    auto tk_xyyyzz_xxyzz = pbuffer.data(idx_kin_ih + 365);

    auto tk_xyyyzz_xxzzz = pbuffer.data(idx_kin_ih + 366);

    auto tk_xyyyzz_xyyyy = pbuffer.data(idx_kin_ih + 367);

    auto tk_xyyyzz_xyyyz = pbuffer.data(idx_kin_ih + 368);

    auto tk_xyyyzz_xyyzz = pbuffer.data(idx_kin_ih + 369);

    auto tk_xyyyzz_xyzzz = pbuffer.data(idx_kin_ih + 370);

    auto tk_xyyyzz_xzzzz = pbuffer.data(idx_kin_ih + 371);

    auto tk_xyyyzz_yyyyy = pbuffer.data(idx_kin_ih + 372);

    auto tk_xyyyzz_yyyyz = pbuffer.data(idx_kin_ih + 373);

    auto tk_xyyyzz_yyyzz = pbuffer.data(idx_kin_ih + 374);

    auto tk_xyyyzz_yyzzz = pbuffer.data(idx_kin_ih + 375);

    auto tk_xyyyzz_yzzzz = pbuffer.data(idx_kin_ih + 376);

    auto tk_xyyyzz_zzzzz = pbuffer.data(idx_kin_ih + 377);

    #pragma omp simd aligned(pa_x, tk_xyyyzz_xxxxx, tk_xyyyzz_xxxxy, tk_xyyyzz_xxxxz, tk_xyyyzz_xxxyy, tk_xyyyzz_xxxyz, tk_xyyyzz_xxxzz, tk_xyyyzz_xxyyy, tk_xyyyzz_xxyyz, tk_xyyyzz_xxyzz, tk_xyyyzz_xxzzz, tk_xyyyzz_xyyyy, tk_xyyyzz_xyyyz, tk_xyyyzz_xyyzz, tk_xyyyzz_xyzzz, tk_xyyyzz_xzzzz, tk_xyyyzz_yyyyy, tk_xyyyzz_yyyyz, tk_xyyyzz_yyyzz, tk_xyyyzz_yyzzz, tk_xyyyzz_yzzzz, tk_xyyyzz_zzzzz, tk_yyyzz_xxxx, tk_yyyzz_xxxxx, tk_yyyzz_xxxxy, tk_yyyzz_xxxxz, tk_yyyzz_xxxy, tk_yyyzz_xxxyy, tk_yyyzz_xxxyz, tk_yyyzz_xxxz, tk_yyyzz_xxxzz, tk_yyyzz_xxyy, tk_yyyzz_xxyyy, tk_yyyzz_xxyyz, tk_yyyzz_xxyz, tk_yyyzz_xxyzz, tk_yyyzz_xxzz, tk_yyyzz_xxzzz, tk_yyyzz_xyyy, tk_yyyzz_xyyyy, tk_yyyzz_xyyyz, tk_yyyzz_xyyz, tk_yyyzz_xyyzz, tk_yyyzz_xyzz, tk_yyyzz_xyzzz, tk_yyyzz_xzzz, tk_yyyzz_xzzzz, tk_yyyzz_yyyy, tk_yyyzz_yyyyy, tk_yyyzz_yyyyz, tk_yyyzz_yyyz, tk_yyyzz_yyyzz, tk_yyyzz_yyzz, tk_yyyzz_yyzzz, tk_yyyzz_yzzz, tk_yyyzz_yzzzz, tk_yyyzz_zzzz, tk_yyyzz_zzzzz, ts_xyyyzz_xxxxx, ts_xyyyzz_xxxxy, ts_xyyyzz_xxxxz, ts_xyyyzz_xxxyy, ts_xyyyzz_xxxyz, ts_xyyyzz_xxxzz, ts_xyyyzz_xxyyy, ts_xyyyzz_xxyyz, ts_xyyyzz_xxyzz, ts_xyyyzz_xxzzz, ts_xyyyzz_xyyyy, ts_xyyyzz_xyyyz, ts_xyyyzz_xyyzz, ts_xyyyzz_xyzzz, ts_xyyyzz_xzzzz, ts_xyyyzz_yyyyy, ts_xyyyzz_yyyyz, ts_xyyyzz_yyyzz, ts_xyyyzz_yyzzz, ts_xyyyzz_yzzzz, ts_xyyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyzz_xxxxx[i] = 5.0 * tk_yyyzz_xxxx[i] * fe_0 + tk_yyyzz_xxxxx[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxx[i] * fz_0;

        tk_xyyyzz_xxxxy[i] = 4.0 * tk_yyyzz_xxxy[i] * fe_0 + tk_yyyzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxy[i] * fz_0;

        tk_xyyyzz_xxxxz[i] = 4.0 * tk_yyyzz_xxxz[i] * fe_0 + tk_yyyzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxxz[i] * fz_0;

        tk_xyyyzz_xxxyy[i] = 3.0 * tk_yyyzz_xxyy[i] * fe_0 + tk_yyyzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxyy[i] * fz_0;

        tk_xyyyzz_xxxyz[i] = 3.0 * tk_yyyzz_xxyz[i] * fe_0 + tk_yyyzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxyz[i] * fz_0;

        tk_xyyyzz_xxxzz[i] = 3.0 * tk_yyyzz_xxzz[i] * fe_0 + tk_yyyzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxxzz[i] * fz_0;

        tk_xyyyzz_xxyyy[i] = 2.0 * tk_yyyzz_xyyy[i] * fe_0 + tk_yyyzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyyy[i] * fz_0;

        tk_xyyyzz_xxyyz[i] = 2.0 * tk_yyyzz_xyyz[i] * fe_0 + tk_yyyzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyyz[i] * fz_0;

        tk_xyyyzz_xxyzz[i] = 2.0 * tk_yyyzz_xyzz[i] * fe_0 + tk_yyyzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxyzz[i] * fz_0;

        tk_xyyyzz_xxzzz[i] = 2.0 * tk_yyyzz_xzzz[i] * fe_0 + tk_yyyzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xxzzz[i] * fz_0;

        tk_xyyyzz_xyyyy[i] = tk_yyyzz_yyyy[i] * fe_0 + tk_yyyzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyyy[i] * fz_0;

        tk_xyyyzz_xyyyz[i] = tk_yyyzz_yyyz[i] * fe_0 + tk_yyyzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyyz[i] * fz_0;

        tk_xyyyzz_xyyzz[i] = tk_yyyzz_yyzz[i] * fe_0 + tk_yyyzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyyzz[i] * fz_0;

        tk_xyyyzz_xyzzz[i] = tk_yyyzz_yzzz[i] * fe_0 + tk_yyyzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xyzzz[i] * fz_0;

        tk_xyyyzz_xzzzz[i] = tk_yyyzz_zzzz[i] * fe_0 + tk_yyyzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xzzzz[i] * fz_0;

        tk_xyyyzz_yyyyy[i] = tk_yyyzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyyy[i] * fz_0;

        tk_xyyyzz_yyyyz[i] = tk_yyyzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyyz[i] * fz_0;

        tk_xyyyzz_yyyzz[i] = tk_yyyzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyyzz[i] * fz_0;

        tk_xyyyzz_yyzzz[i] = tk_yyyzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yyzzz[i] * fz_0;

        tk_xyyyzz_yzzzz[i] = tk_yyyzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yzzzz[i] * fz_0;

        tk_xyyyzz_zzzzz[i] = tk_yyyzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyyzz_zzzzz[i] * fz_0;
    }

    // Set up 378-399 components of targeted buffer : IH

    auto tk_xyyzzz_xxxxx = pbuffer.data(idx_kin_ih + 378);

    auto tk_xyyzzz_xxxxy = pbuffer.data(idx_kin_ih + 379);

    auto tk_xyyzzz_xxxxz = pbuffer.data(idx_kin_ih + 380);

    auto tk_xyyzzz_xxxyy = pbuffer.data(idx_kin_ih + 381);

    auto tk_xyyzzz_xxxyz = pbuffer.data(idx_kin_ih + 382);

    auto tk_xyyzzz_xxxzz = pbuffer.data(idx_kin_ih + 383);

    auto tk_xyyzzz_xxyyy = pbuffer.data(idx_kin_ih + 384);

    auto tk_xyyzzz_xxyyz = pbuffer.data(idx_kin_ih + 385);

    auto tk_xyyzzz_xxyzz = pbuffer.data(idx_kin_ih + 386);

    auto tk_xyyzzz_xxzzz = pbuffer.data(idx_kin_ih + 387);

    auto tk_xyyzzz_xyyyy = pbuffer.data(idx_kin_ih + 388);

    auto tk_xyyzzz_xyyyz = pbuffer.data(idx_kin_ih + 389);

    auto tk_xyyzzz_xyyzz = pbuffer.data(idx_kin_ih + 390);

    auto tk_xyyzzz_xyzzz = pbuffer.data(idx_kin_ih + 391);

    auto tk_xyyzzz_xzzzz = pbuffer.data(idx_kin_ih + 392);

    auto tk_xyyzzz_yyyyy = pbuffer.data(idx_kin_ih + 393);

    auto tk_xyyzzz_yyyyz = pbuffer.data(idx_kin_ih + 394);

    auto tk_xyyzzz_yyyzz = pbuffer.data(idx_kin_ih + 395);

    auto tk_xyyzzz_yyzzz = pbuffer.data(idx_kin_ih + 396);

    auto tk_xyyzzz_yzzzz = pbuffer.data(idx_kin_ih + 397);

    auto tk_xyyzzz_zzzzz = pbuffer.data(idx_kin_ih + 398);

    #pragma omp simd aligned(pa_x, tk_xyyzzz_xxxxx, tk_xyyzzz_xxxxy, tk_xyyzzz_xxxxz, tk_xyyzzz_xxxyy, tk_xyyzzz_xxxyz, tk_xyyzzz_xxxzz, tk_xyyzzz_xxyyy, tk_xyyzzz_xxyyz, tk_xyyzzz_xxyzz, tk_xyyzzz_xxzzz, tk_xyyzzz_xyyyy, tk_xyyzzz_xyyyz, tk_xyyzzz_xyyzz, tk_xyyzzz_xyzzz, tk_xyyzzz_xzzzz, tk_xyyzzz_yyyyy, tk_xyyzzz_yyyyz, tk_xyyzzz_yyyzz, tk_xyyzzz_yyzzz, tk_xyyzzz_yzzzz, tk_xyyzzz_zzzzz, tk_yyzzz_xxxx, tk_yyzzz_xxxxx, tk_yyzzz_xxxxy, tk_yyzzz_xxxxz, tk_yyzzz_xxxy, tk_yyzzz_xxxyy, tk_yyzzz_xxxyz, tk_yyzzz_xxxz, tk_yyzzz_xxxzz, tk_yyzzz_xxyy, tk_yyzzz_xxyyy, tk_yyzzz_xxyyz, tk_yyzzz_xxyz, tk_yyzzz_xxyzz, tk_yyzzz_xxzz, tk_yyzzz_xxzzz, tk_yyzzz_xyyy, tk_yyzzz_xyyyy, tk_yyzzz_xyyyz, tk_yyzzz_xyyz, tk_yyzzz_xyyzz, tk_yyzzz_xyzz, tk_yyzzz_xyzzz, tk_yyzzz_xzzz, tk_yyzzz_xzzzz, tk_yyzzz_yyyy, tk_yyzzz_yyyyy, tk_yyzzz_yyyyz, tk_yyzzz_yyyz, tk_yyzzz_yyyzz, tk_yyzzz_yyzz, tk_yyzzz_yyzzz, tk_yyzzz_yzzz, tk_yyzzz_yzzzz, tk_yyzzz_zzzz, tk_yyzzz_zzzzz, ts_xyyzzz_xxxxx, ts_xyyzzz_xxxxy, ts_xyyzzz_xxxxz, ts_xyyzzz_xxxyy, ts_xyyzzz_xxxyz, ts_xyyzzz_xxxzz, ts_xyyzzz_xxyyy, ts_xyyzzz_xxyyz, ts_xyyzzz_xxyzz, ts_xyyzzz_xxzzz, ts_xyyzzz_xyyyy, ts_xyyzzz_xyyyz, ts_xyyzzz_xyyzz, ts_xyyzzz_xyzzz, ts_xyyzzz_xzzzz, ts_xyyzzz_yyyyy, ts_xyyzzz_yyyyz, ts_xyyzzz_yyyzz, ts_xyyzzz_yyzzz, ts_xyyzzz_yzzzz, ts_xyyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzzz_xxxxx[i] = 5.0 * tk_yyzzz_xxxx[i] * fe_0 + tk_yyzzz_xxxxx[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxx[i] * fz_0;

        tk_xyyzzz_xxxxy[i] = 4.0 * tk_yyzzz_xxxy[i] * fe_0 + tk_yyzzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxy[i] * fz_0;

        tk_xyyzzz_xxxxz[i] = 4.0 * tk_yyzzz_xxxz[i] * fe_0 + tk_yyzzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxxz[i] * fz_0;

        tk_xyyzzz_xxxyy[i] = 3.0 * tk_yyzzz_xxyy[i] * fe_0 + tk_yyzzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxyy[i] * fz_0;

        tk_xyyzzz_xxxyz[i] = 3.0 * tk_yyzzz_xxyz[i] * fe_0 + tk_yyzzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxyz[i] * fz_0;

        tk_xyyzzz_xxxzz[i] = 3.0 * tk_yyzzz_xxzz[i] * fe_0 + tk_yyzzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxxzz[i] * fz_0;

        tk_xyyzzz_xxyyy[i] = 2.0 * tk_yyzzz_xyyy[i] * fe_0 + tk_yyzzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyyy[i] * fz_0;

        tk_xyyzzz_xxyyz[i] = 2.0 * tk_yyzzz_xyyz[i] * fe_0 + tk_yyzzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyyz[i] * fz_0;

        tk_xyyzzz_xxyzz[i] = 2.0 * tk_yyzzz_xyzz[i] * fe_0 + tk_yyzzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxyzz[i] * fz_0;

        tk_xyyzzz_xxzzz[i] = 2.0 * tk_yyzzz_xzzz[i] * fe_0 + tk_yyzzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xxzzz[i] * fz_0;

        tk_xyyzzz_xyyyy[i] = tk_yyzzz_yyyy[i] * fe_0 + tk_yyzzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyyy[i] * fz_0;

        tk_xyyzzz_xyyyz[i] = tk_yyzzz_yyyz[i] * fe_0 + tk_yyzzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyyz[i] * fz_0;

        tk_xyyzzz_xyyzz[i] = tk_yyzzz_yyzz[i] * fe_0 + tk_yyzzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyyzz[i] * fz_0;

        tk_xyyzzz_xyzzz[i] = tk_yyzzz_yzzz[i] * fe_0 + tk_yyzzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xyzzz[i] * fz_0;

        tk_xyyzzz_xzzzz[i] = tk_yyzzz_zzzz[i] * fe_0 + tk_yyzzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xzzzz[i] * fz_0;

        tk_xyyzzz_yyyyy[i] = tk_yyzzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyyy[i] * fz_0;

        tk_xyyzzz_yyyyz[i] = tk_yyzzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyyz[i] * fz_0;

        tk_xyyzzz_yyyzz[i] = tk_yyzzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyyzz[i] * fz_0;

        tk_xyyzzz_yyzzz[i] = tk_yyzzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yyzzz[i] * fz_0;

        tk_xyyzzz_yzzzz[i] = tk_yyzzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yzzzz[i] * fz_0;

        tk_xyyzzz_zzzzz[i] = tk_yyzzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyyzzz_zzzzz[i] * fz_0;
    }

    // Set up 399-420 components of targeted buffer : IH

    auto tk_xyzzzz_xxxxx = pbuffer.data(idx_kin_ih + 399);

    auto tk_xyzzzz_xxxxy = pbuffer.data(idx_kin_ih + 400);

    auto tk_xyzzzz_xxxxz = pbuffer.data(idx_kin_ih + 401);

    auto tk_xyzzzz_xxxyy = pbuffer.data(idx_kin_ih + 402);

    auto tk_xyzzzz_xxxyz = pbuffer.data(idx_kin_ih + 403);

    auto tk_xyzzzz_xxxzz = pbuffer.data(idx_kin_ih + 404);

    auto tk_xyzzzz_xxyyy = pbuffer.data(idx_kin_ih + 405);

    auto tk_xyzzzz_xxyyz = pbuffer.data(idx_kin_ih + 406);

    auto tk_xyzzzz_xxyzz = pbuffer.data(idx_kin_ih + 407);

    auto tk_xyzzzz_xxzzz = pbuffer.data(idx_kin_ih + 408);

    auto tk_xyzzzz_xyyyy = pbuffer.data(idx_kin_ih + 409);

    auto tk_xyzzzz_xyyyz = pbuffer.data(idx_kin_ih + 410);

    auto tk_xyzzzz_xyyzz = pbuffer.data(idx_kin_ih + 411);

    auto tk_xyzzzz_xyzzz = pbuffer.data(idx_kin_ih + 412);

    auto tk_xyzzzz_xzzzz = pbuffer.data(idx_kin_ih + 413);

    auto tk_xyzzzz_yyyyy = pbuffer.data(idx_kin_ih + 414);

    auto tk_xyzzzz_yyyyz = pbuffer.data(idx_kin_ih + 415);

    auto tk_xyzzzz_yyyzz = pbuffer.data(idx_kin_ih + 416);

    auto tk_xyzzzz_yyzzz = pbuffer.data(idx_kin_ih + 417);

    auto tk_xyzzzz_yzzzz = pbuffer.data(idx_kin_ih + 418);

    auto tk_xyzzzz_zzzzz = pbuffer.data(idx_kin_ih + 419);

    #pragma omp simd aligned(pa_x, pa_y, tk_xyzzzz_xxxxx, tk_xyzzzz_xxxxy, tk_xyzzzz_xxxxz, tk_xyzzzz_xxxyy, tk_xyzzzz_xxxyz, tk_xyzzzz_xxxzz, tk_xyzzzz_xxyyy, tk_xyzzzz_xxyyz, tk_xyzzzz_xxyzz, tk_xyzzzz_xxzzz, tk_xyzzzz_xyyyy, tk_xyzzzz_xyyyz, tk_xyzzzz_xyyzz, tk_xyzzzz_xyzzz, tk_xyzzzz_xzzzz, tk_xyzzzz_yyyyy, tk_xyzzzz_yyyyz, tk_xyzzzz_yyyzz, tk_xyzzzz_yyzzz, tk_xyzzzz_yzzzz, tk_xyzzzz_zzzzz, tk_xzzzz_xxxxx, tk_xzzzz_xxxxz, tk_xzzzz_xxxzz, tk_xzzzz_xxzzz, tk_xzzzz_xzzzz, tk_yzzzz_xxxxy, tk_yzzzz_xxxy, tk_yzzzz_xxxyy, tk_yzzzz_xxxyz, tk_yzzzz_xxyy, tk_yzzzz_xxyyy, tk_yzzzz_xxyyz, tk_yzzzz_xxyz, tk_yzzzz_xxyzz, tk_yzzzz_xyyy, tk_yzzzz_xyyyy, tk_yzzzz_xyyyz, tk_yzzzz_xyyz, tk_yzzzz_xyyzz, tk_yzzzz_xyzz, tk_yzzzz_xyzzz, tk_yzzzz_yyyy, tk_yzzzz_yyyyy, tk_yzzzz_yyyyz, tk_yzzzz_yyyz, tk_yzzzz_yyyzz, tk_yzzzz_yyzz, tk_yzzzz_yyzzz, tk_yzzzz_yzzz, tk_yzzzz_yzzzz, tk_yzzzz_zzzzz, ts_xyzzzz_xxxxx, ts_xyzzzz_xxxxy, ts_xyzzzz_xxxxz, ts_xyzzzz_xxxyy, ts_xyzzzz_xxxyz, ts_xyzzzz_xxxzz, ts_xyzzzz_xxyyy, ts_xyzzzz_xxyyz, ts_xyzzzz_xxyzz, ts_xyzzzz_xxzzz, ts_xyzzzz_xyyyy, ts_xyzzzz_xyyyz, ts_xyzzzz_xyyzz, ts_xyzzzz_xyzzz, ts_xyzzzz_xzzzz, ts_xyzzzz_yyyyy, ts_xyzzzz_yyyyz, ts_xyzzzz_yyyzz, ts_xyzzzz_yyzzz, ts_xyzzzz_yzzzz, ts_xyzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzzz_xxxxx[i] = tk_xzzzz_xxxxx[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxxx[i] * fz_0;

        tk_xyzzzz_xxxxy[i] = 4.0 * tk_yzzzz_xxxy[i] * fe_0 + tk_yzzzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxxy[i] * fz_0;

        tk_xyzzzz_xxxxz[i] = tk_xzzzz_xxxxz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxxz[i] * fz_0;

        tk_xyzzzz_xxxyy[i] = 3.0 * tk_yzzzz_xxyy[i] * fe_0 + tk_yzzzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxyy[i] * fz_0;

        tk_xyzzzz_xxxyz[i] = 3.0 * tk_yzzzz_xxyz[i] * fe_0 + tk_yzzzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxxyz[i] * fz_0;

        tk_xyzzzz_xxxzz[i] = tk_xzzzz_xxxzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxxzz[i] * fz_0;

        tk_xyzzzz_xxyyy[i] = 2.0 * tk_yzzzz_xyyy[i] * fe_0 + tk_yzzzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyyy[i] * fz_0;

        tk_xyzzzz_xxyyz[i] = 2.0 * tk_yzzzz_xyyz[i] * fe_0 + tk_yzzzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyyz[i] * fz_0;

        tk_xyzzzz_xxyzz[i] = 2.0 * tk_yzzzz_xyzz[i] * fe_0 + tk_yzzzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xxyzz[i] * fz_0;

        tk_xyzzzz_xxzzz[i] = tk_xzzzz_xxzzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xxzzz[i] * fz_0;

        tk_xyzzzz_xyyyy[i] = tk_yzzzz_yyyy[i] * fe_0 + tk_yzzzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyyy[i] * fz_0;

        tk_xyzzzz_xyyyz[i] = tk_yzzzz_yyyz[i] * fe_0 + tk_yzzzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyyz[i] * fz_0;

        tk_xyzzzz_xyyzz[i] = tk_yzzzz_yyzz[i] * fe_0 + tk_yzzzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyyzz[i] * fz_0;

        tk_xyzzzz_xyzzz[i] = tk_yzzzz_yzzz[i] * fe_0 + tk_yzzzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_xyzzz[i] * fz_0;

        tk_xyzzzz_xzzzz[i] = tk_xzzzz_xzzzz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xzzzz[i] * fz_0;

        tk_xyzzzz_yyyyy[i] = tk_yzzzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyyy[i] * fz_0;

        tk_xyzzzz_yyyyz[i] = tk_yzzzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyyz[i] * fz_0;

        tk_xyzzzz_yyyzz[i] = tk_yzzzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyyzz[i] * fz_0;

        tk_xyzzzz_yyzzz[i] = tk_yzzzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yyzzz[i] * fz_0;

        tk_xyzzzz_yzzzz[i] = tk_yzzzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yzzzz[i] * fz_0;

        tk_xyzzzz_zzzzz[i] = tk_yzzzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyzzzz_zzzzz[i] * fz_0;
    }

    // Set up 420-441 components of targeted buffer : IH

    auto tk_xzzzzz_xxxxx = pbuffer.data(idx_kin_ih + 420);

    auto tk_xzzzzz_xxxxy = pbuffer.data(idx_kin_ih + 421);

    auto tk_xzzzzz_xxxxz = pbuffer.data(idx_kin_ih + 422);

    auto tk_xzzzzz_xxxyy = pbuffer.data(idx_kin_ih + 423);

    auto tk_xzzzzz_xxxyz = pbuffer.data(idx_kin_ih + 424);

    auto tk_xzzzzz_xxxzz = pbuffer.data(idx_kin_ih + 425);

    auto tk_xzzzzz_xxyyy = pbuffer.data(idx_kin_ih + 426);

    auto tk_xzzzzz_xxyyz = pbuffer.data(idx_kin_ih + 427);

    auto tk_xzzzzz_xxyzz = pbuffer.data(idx_kin_ih + 428);

    auto tk_xzzzzz_xxzzz = pbuffer.data(idx_kin_ih + 429);

    auto tk_xzzzzz_xyyyy = pbuffer.data(idx_kin_ih + 430);

    auto tk_xzzzzz_xyyyz = pbuffer.data(idx_kin_ih + 431);

    auto tk_xzzzzz_xyyzz = pbuffer.data(idx_kin_ih + 432);

    auto tk_xzzzzz_xyzzz = pbuffer.data(idx_kin_ih + 433);

    auto tk_xzzzzz_xzzzz = pbuffer.data(idx_kin_ih + 434);

    auto tk_xzzzzz_yyyyy = pbuffer.data(idx_kin_ih + 435);

    auto tk_xzzzzz_yyyyz = pbuffer.data(idx_kin_ih + 436);

    auto tk_xzzzzz_yyyzz = pbuffer.data(idx_kin_ih + 437);

    auto tk_xzzzzz_yyzzz = pbuffer.data(idx_kin_ih + 438);

    auto tk_xzzzzz_yzzzz = pbuffer.data(idx_kin_ih + 439);

    auto tk_xzzzzz_zzzzz = pbuffer.data(idx_kin_ih + 440);

    #pragma omp simd aligned(pa_x, tk_xzzzzz_xxxxx, tk_xzzzzz_xxxxy, tk_xzzzzz_xxxxz, tk_xzzzzz_xxxyy, tk_xzzzzz_xxxyz, tk_xzzzzz_xxxzz, tk_xzzzzz_xxyyy, tk_xzzzzz_xxyyz, tk_xzzzzz_xxyzz, tk_xzzzzz_xxzzz, tk_xzzzzz_xyyyy, tk_xzzzzz_xyyyz, tk_xzzzzz_xyyzz, tk_xzzzzz_xyzzz, tk_xzzzzz_xzzzz, tk_xzzzzz_yyyyy, tk_xzzzzz_yyyyz, tk_xzzzzz_yyyzz, tk_xzzzzz_yyzzz, tk_xzzzzz_yzzzz, tk_xzzzzz_zzzzz, tk_zzzzz_xxxx, tk_zzzzz_xxxxx, tk_zzzzz_xxxxy, tk_zzzzz_xxxxz, tk_zzzzz_xxxy, tk_zzzzz_xxxyy, tk_zzzzz_xxxyz, tk_zzzzz_xxxz, tk_zzzzz_xxxzz, tk_zzzzz_xxyy, tk_zzzzz_xxyyy, tk_zzzzz_xxyyz, tk_zzzzz_xxyz, tk_zzzzz_xxyzz, tk_zzzzz_xxzz, tk_zzzzz_xxzzz, tk_zzzzz_xyyy, tk_zzzzz_xyyyy, tk_zzzzz_xyyyz, tk_zzzzz_xyyz, tk_zzzzz_xyyzz, tk_zzzzz_xyzz, tk_zzzzz_xyzzz, tk_zzzzz_xzzz, tk_zzzzz_xzzzz, tk_zzzzz_yyyy, tk_zzzzz_yyyyy, tk_zzzzz_yyyyz, tk_zzzzz_yyyz, tk_zzzzz_yyyzz, tk_zzzzz_yyzz, tk_zzzzz_yyzzz, tk_zzzzz_yzzz, tk_zzzzz_yzzzz, tk_zzzzz_zzzz, tk_zzzzz_zzzzz, ts_xzzzzz_xxxxx, ts_xzzzzz_xxxxy, ts_xzzzzz_xxxxz, ts_xzzzzz_xxxyy, ts_xzzzzz_xxxyz, ts_xzzzzz_xxxzz, ts_xzzzzz_xxyyy, ts_xzzzzz_xxyyz, ts_xzzzzz_xxyzz, ts_xzzzzz_xxzzz, ts_xzzzzz_xyyyy, ts_xzzzzz_xyyyz, ts_xzzzzz_xyyzz, ts_xzzzzz_xyzzz, ts_xzzzzz_xzzzz, ts_xzzzzz_yyyyy, ts_xzzzzz_yyyyz, ts_xzzzzz_yyyzz, ts_xzzzzz_yyzzz, ts_xzzzzz_yzzzz, ts_xzzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzzz_xxxxx[i] = 5.0 * tk_zzzzz_xxxx[i] * fe_0 + tk_zzzzz_xxxxx[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxx[i] * fz_0;

        tk_xzzzzz_xxxxy[i] = 4.0 * tk_zzzzz_xxxy[i] * fe_0 + tk_zzzzz_xxxxy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxy[i] * fz_0;

        tk_xzzzzz_xxxxz[i] = 4.0 * tk_zzzzz_xxxz[i] * fe_0 + tk_zzzzz_xxxxz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxxz[i] * fz_0;

        tk_xzzzzz_xxxyy[i] = 3.0 * tk_zzzzz_xxyy[i] * fe_0 + tk_zzzzz_xxxyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxyy[i] * fz_0;

        tk_xzzzzz_xxxyz[i] = 3.0 * tk_zzzzz_xxyz[i] * fe_0 + tk_zzzzz_xxxyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxyz[i] * fz_0;

        tk_xzzzzz_xxxzz[i] = 3.0 * tk_zzzzz_xxzz[i] * fe_0 + tk_zzzzz_xxxzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxxzz[i] * fz_0;

        tk_xzzzzz_xxyyy[i] = 2.0 * tk_zzzzz_xyyy[i] * fe_0 + tk_zzzzz_xxyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyyy[i] * fz_0;

        tk_xzzzzz_xxyyz[i] = 2.0 * tk_zzzzz_xyyz[i] * fe_0 + tk_zzzzz_xxyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyyz[i] * fz_0;

        tk_xzzzzz_xxyzz[i] = 2.0 * tk_zzzzz_xyzz[i] * fe_0 + tk_zzzzz_xxyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxyzz[i] * fz_0;

        tk_xzzzzz_xxzzz[i] = 2.0 * tk_zzzzz_xzzz[i] * fe_0 + tk_zzzzz_xxzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xxzzz[i] * fz_0;

        tk_xzzzzz_xyyyy[i] = tk_zzzzz_yyyy[i] * fe_0 + tk_zzzzz_xyyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyyy[i] * fz_0;

        tk_xzzzzz_xyyyz[i] = tk_zzzzz_yyyz[i] * fe_0 + tk_zzzzz_xyyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyyz[i] * fz_0;

        tk_xzzzzz_xyyzz[i] = tk_zzzzz_yyzz[i] * fe_0 + tk_zzzzz_xyyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyyzz[i] * fz_0;

        tk_xzzzzz_xyzzz[i] = tk_zzzzz_yzzz[i] * fe_0 + tk_zzzzz_xyzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xyzzz[i] * fz_0;

        tk_xzzzzz_xzzzz[i] = tk_zzzzz_zzzz[i] * fe_0 + tk_zzzzz_xzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xzzzz[i] * fz_0;

        tk_xzzzzz_yyyyy[i] = tk_zzzzz_yyyyy[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyyy[i] * fz_0;

        tk_xzzzzz_yyyyz[i] = tk_zzzzz_yyyyz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyyz[i] * fz_0;

        tk_xzzzzz_yyyzz[i] = tk_zzzzz_yyyzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyyzz[i] * fz_0;

        tk_xzzzzz_yyzzz[i] = tk_zzzzz_yyzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yyzzz[i] * fz_0;

        tk_xzzzzz_yzzzz[i] = tk_zzzzz_yzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yzzzz[i] * fz_0;

        tk_xzzzzz_zzzzz[i] = tk_zzzzz_zzzzz[i] * pa_x[i] + 2.0 * ts_xzzzzz_zzzzz[i] * fz_0;
    }

    // Set up 441-462 components of targeted buffer : IH

    auto tk_yyyyyy_xxxxx = pbuffer.data(idx_kin_ih + 441);

    auto tk_yyyyyy_xxxxy = pbuffer.data(idx_kin_ih + 442);

    auto tk_yyyyyy_xxxxz = pbuffer.data(idx_kin_ih + 443);

    auto tk_yyyyyy_xxxyy = pbuffer.data(idx_kin_ih + 444);

    auto tk_yyyyyy_xxxyz = pbuffer.data(idx_kin_ih + 445);

    auto tk_yyyyyy_xxxzz = pbuffer.data(idx_kin_ih + 446);

    auto tk_yyyyyy_xxyyy = pbuffer.data(idx_kin_ih + 447);

    auto tk_yyyyyy_xxyyz = pbuffer.data(idx_kin_ih + 448);

    auto tk_yyyyyy_xxyzz = pbuffer.data(idx_kin_ih + 449);

    auto tk_yyyyyy_xxzzz = pbuffer.data(idx_kin_ih + 450);

    auto tk_yyyyyy_xyyyy = pbuffer.data(idx_kin_ih + 451);

    auto tk_yyyyyy_xyyyz = pbuffer.data(idx_kin_ih + 452);

    auto tk_yyyyyy_xyyzz = pbuffer.data(idx_kin_ih + 453);

    auto tk_yyyyyy_xyzzz = pbuffer.data(idx_kin_ih + 454);

    auto tk_yyyyyy_xzzzz = pbuffer.data(idx_kin_ih + 455);

    auto tk_yyyyyy_yyyyy = pbuffer.data(idx_kin_ih + 456);

    auto tk_yyyyyy_yyyyz = pbuffer.data(idx_kin_ih + 457);

    auto tk_yyyyyy_yyyzz = pbuffer.data(idx_kin_ih + 458);

    auto tk_yyyyyy_yyzzz = pbuffer.data(idx_kin_ih + 459);

    auto tk_yyyyyy_yzzzz = pbuffer.data(idx_kin_ih + 460);

    auto tk_yyyyyy_zzzzz = pbuffer.data(idx_kin_ih + 461);

    #pragma omp simd aligned(pa_y, tk_yyyy_xxxxx, tk_yyyy_xxxxy, tk_yyyy_xxxxz, tk_yyyy_xxxyy, tk_yyyy_xxxyz, tk_yyyy_xxxzz, tk_yyyy_xxyyy, tk_yyyy_xxyyz, tk_yyyy_xxyzz, tk_yyyy_xxzzz, tk_yyyy_xyyyy, tk_yyyy_xyyyz, tk_yyyy_xyyzz, tk_yyyy_xyzzz, tk_yyyy_xzzzz, tk_yyyy_yyyyy, tk_yyyy_yyyyz, tk_yyyy_yyyzz, tk_yyyy_yyzzz, tk_yyyy_yzzzz, tk_yyyy_zzzzz, tk_yyyyy_xxxx, tk_yyyyy_xxxxx, tk_yyyyy_xxxxy, tk_yyyyy_xxxxz, tk_yyyyy_xxxy, tk_yyyyy_xxxyy, tk_yyyyy_xxxyz, tk_yyyyy_xxxz, tk_yyyyy_xxxzz, tk_yyyyy_xxyy, tk_yyyyy_xxyyy, tk_yyyyy_xxyyz, tk_yyyyy_xxyz, tk_yyyyy_xxyzz, tk_yyyyy_xxzz, tk_yyyyy_xxzzz, tk_yyyyy_xyyy, tk_yyyyy_xyyyy, tk_yyyyy_xyyyz, tk_yyyyy_xyyz, tk_yyyyy_xyyzz, tk_yyyyy_xyzz, tk_yyyyy_xyzzz, tk_yyyyy_xzzz, tk_yyyyy_xzzzz, tk_yyyyy_yyyy, tk_yyyyy_yyyyy, tk_yyyyy_yyyyz, tk_yyyyy_yyyz, tk_yyyyy_yyyzz, tk_yyyyy_yyzz, tk_yyyyy_yyzzz, tk_yyyyy_yzzz, tk_yyyyy_yzzzz, tk_yyyyy_zzzz, tk_yyyyy_zzzzz, tk_yyyyyy_xxxxx, tk_yyyyyy_xxxxy, tk_yyyyyy_xxxxz, tk_yyyyyy_xxxyy, tk_yyyyyy_xxxyz, tk_yyyyyy_xxxzz, tk_yyyyyy_xxyyy, tk_yyyyyy_xxyyz, tk_yyyyyy_xxyzz, tk_yyyyyy_xxzzz, tk_yyyyyy_xyyyy, tk_yyyyyy_xyyyz, tk_yyyyyy_xyyzz, tk_yyyyyy_xyzzz, tk_yyyyyy_xzzzz, tk_yyyyyy_yyyyy, tk_yyyyyy_yyyyz, tk_yyyyyy_yyyzz, tk_yyyyyy_yyzzz, tk_yyyyyy_yzzzz, tk_yyyyyy_zzzzz, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxzz, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyzz, ts_yyyy_xxzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyzz, ts_yyyy_xyzzz, ts_yyyy_xzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyzz, ts_yyyy_yyzzz, ts_yyyy_yzzzz, ts_yyyy_zzzzz, ts_yyyyyy_xxxxx, ts_yyyyyy_xxxxy, ts_yyyyyy_xxxxz, ts_yyyyyy_xxxyy, ts_yyyyyy_xxxyz, ts_yyyyyy_xxxzz, ts_yyyyyy_xxyyy, ts_yyyyyy_xxyyz, ts_yyyyyy_xxyzz, ts_yyyyyy_xxzzz, ts_yyyyyy_xyyyy, ts_yyyyyy_xyyyz, ts_yyyyyy_xyyzz, ts_yyyyyy_xyzzz, ts_yyyyyy_xzzzz, ts_yyyyyy_yyyyy, ts_yyyyyy_yyyyz, ts_yyyyyy_yyyzz, ts_yyyyyy_yyzzz, ts_yyyyyy_yzzzz, ts_yyyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyyy_xxxxx[i] = -10.0 * ts_yyyy_xxxxx[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxx[i] * fe_0 + tk_yyyyy_xxxxx[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxxx[i] * fz_0;

        tk_yyyyyy_xxxxy[i] = -10.0 * ts_yyyy_xxxxy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxy[i] * fe_0 + tk_yyyyy_xxxx[i] * fe_0 + tk_yyyyy_xxxxy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxxy[i] * fz_0;

        tk_yyyyyy_xxxxz[i] = -10.0 * ts_yyyy_xxxxz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxxz[i] * fe_0 + tk_yyyyy_xxxxz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxxz[i] * fz_0;

        tk_yyyyyy_xxxyy[i] = -10.0 * ts_yyyy_xxxyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxyy[i] * fe_0 + 2.0 * tk_yyyyy_xxxy[i] * fe_0 + tk_yyyyy_xxxyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxyy[i] * fz_0;

        tk_yyyyyy_xxxyz[i] = -10.0 * ts_yyyy_xxxyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxyz[i] * fe_0 + tk_yyyyy_xxxz[i] * fe_0 + tk_yyyyy_xxxyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxyz[i] * fz_0;

        tk_yyyyyy_xxxzz[i] = -10.0 * ts_yyyy_xxxzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxxzz[i] * fe_0 + tk_yyyyy_xxxzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxxzz[i] * fz_0;

        tk_yyyyyy_xxyyy[i] = -10.0 * ts_yyyy_xxyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyyy[i] * fe_0 + 3.0 * tk_yyyyy_xxyy[i] * fe_0 + tk_yyyyy_xxyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyyy[i] * fz_0;

        tk_yyyyyy_xxyyz[i] = -10.0 * ts_yyyy_xxyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyyz[i] * fe_0 + 2.0 * tk_yyyyy_xxyz[i] * fe_0 + tk_yyyyy_xxyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyyz[i] * fz_0;

        tk_yyyyyy_xxyzz[i] = -10.0 * ts_yyyy_xxyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxyzz[i] * fe_0 + tk_yyyyy_xxzz[i] * fe_0 + tk_yyyyy_xxyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxyzz[i] * fz_0;

        tk_yyyyyy_xxzzz[i] = -10.0 * ts_yyyy_xxzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xxzzz[i] * fe_0 + tk_yyyyy_xxzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xxzzz[i] * fz_0;

        tk_yyyyyy_xyyyy[i] = -10.0 * ts_yyyy_xyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyyy[i] * fe_0 + 4.0 * tk_yyyyy_xyyy[i] * fe_0 + tk_yyyyy_xyyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyyy[i] * fz_0;

        tk_yyyyyy_xyyyz[i] = -10.0 * ts_yyyy_xyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyyz[i] * fe_0 + 3.0 * tk_yyyyy_xyyz[i] * fe_0 + tk_yyyyy_xyyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyyz[i] * fz_0;

        tk_yyyyyy_xyyzz[i] = -10.0 * ts_yyyy_xyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyyzz[i] * fe_0 + 2.0 * tk_yyyyy_xyzz[i] * fe_0 + tk_yyyyy_xyyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyyzz[i] * fz_0;

        tk_yyyyyy_xyzzz[i] = -10.0 * ts_yyyy_xyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xyzzz[i] * fe_0 + tk_yyyyy_xzzz[i] * fe_0 + tk_yyyyy_xyzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xyzzz[i] * fz_0;

        tk_yyyyyy_xzzzz[i] = -10.0 * ts_yyyy_xzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xzzzz[i] * fe_0 + tk_yyyyy_xzzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xzzzz[i] * fz_0;

        tk_yyyyyy_yyyyy[i] = -10.0 * ts_yyyy_yyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyyy[i] * fe_0 + 5.0 * tk_yyyyy_yyyy[i] * fe_0 + tk_yyyyy_yyyyy[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyyy[i] * fz_0;

        tk_yyyyyy_yyyyz[i] = -10.0 * ts_yyyy_yyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyyz[i] * fe_0 + 4.0 * tk_yyyyy_yyyz[i] * fe_0 + tk_yyyyy_yyyyz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyyz[i] * fz_0;

        tk_yyyyyy_yyyzz[i] = -10.0 * ts_yyyy_yyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyyzz[i] * fe_0 + 3.0 * tk_yyyyy_yyzz[i] * fe_0 + tk_yyyyy_yyyzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyyzz[i] * fz_0;

        tk_yyyyyy_yyzzz[i] = -10.0 * ts_yyyy_yyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yyzzz[i] * fe_0 + 2.0 * tk_yyyyy_yzzz[i] * fe_0 + tk_yyyyy_yyzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yyzzz[i] * fz_0;

        tk_yyyyyy_yzzzz[i] = -10.0 * ts_yyyy_yzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yzzzz[i] * fe_0 + tk_yyyyy_zzzz[i] * fe_0 + tk_yyyyy_yzzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_yzzzz[i] * fz_0;

        tk_yyyyyy_zzzzz[i] = -10.0 * ts_yyyy_zzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_zzzzz[i] * fe_0 + tk_yyyyy_zzzzz[i] * pa_y[i] + 2.0 * ts_yyyyyy_zzzzz[i] * fz_0;
    }

    // Set up 462-483 components of targeted buffer : IH

    auto tk_yyyyyz_xxxxx = pbuffer.data(idx_kin_ih + 462);

    auto tk_yyyyyz_xxxxy = pbuffer.data(idx_kin_ih + 463);

    auto tk_yyyyyz_xxxxz = pbuffer.data(idx_kin_ih + 464);

    auto tk_yyyyyz_xxxyy = pbuffer.data(idx_kin_ih + 465);

    auto tk_yyyyyz_xxxyz = pbuffer.data(idx_kin_ih + 466);

    auto tk_yyyyyz_xxxzz = pbuffer.data(idx_kin_ih + 467);

    auto tk_yyyyyz_xxyyy = pbuffer.data(idx_kin_ih + 468);

    auto tk_yyyyyz_xxyyz = pbuffer.data(idx_kin_ih + 469);

    auto tk_yyyyyz_xxyzz = pbuffer.data(idx_kin_ih + 470);

    auto tk_yyyyyz_xxzzz = pbuffer.data(idx_kin_ih + 471);

    auto tk_yyyyyz_xyyyy = pbuffer.data(idx_kin_ih + 472);

    auto tk_yyyyyz_xyyyz = pbuffer.data(idx_kin_ih + 473);

    auto tk_yyyyyz_xyyzz = pbuffer.data(idx_kin_ih + 474);

    auto tk_yyyyyz_xyzzz = pbuffer.data(idx_kin_ih + 475);

    auto tk_yyyyyz_xzzzz = pbuffer.data(idx_kin_ih + 476);

    auto tk_yyyyyz_yyyyy = pbuffer.data(idx_kin_ih + 477);

    auto tk_yyyyyz_yyyyz = pbuffer.data(idx_kin_ih + 478);

    auto tk_yyyyyz_yyyzz = pbuffer.data(idx_kin_ih + 479);

    auto tk_yyyyyz_yyzzz = pbuffer.data(idx_kin_ih + 480);

    auto tk_yyyyyz_yzzzz = pbuffer.data(idx_kin_ih + 481);

    auto tk_yyyyyz_zzzzz = pbuffer.data(idx_kin_ih + 482);

    #pragma omp simd aligned(pa_z, tk_yyyyy_xxxx, tk_yyyyy_xxxxx, tk_yyyyy_xxxxy, tk_yyyyy_xxxxz, tk_yyyyy_xxxy, tk_yyyyy_xxxyy, tk_yyyyy_xxxyz, tk_yyyyy_xxxz, tk_yyyyy_xxxzz, tk_yyyyy_xxyy, tk_yyyyy_xxyyy, tk_yyyyy_xxyyz, tk_yyyyy_xxyz, tk_yyyyy_xxyzz, tk_yyyyy_xxzz, tk_yyyyy_xxzzz, tk_yyyyy_xyyy, tk_yyyyy_xyyyy, tk_yyyyy_xyyyz, tk_yyyyy_xyyz, tk_yyyyy_xyyzz, tk_yyyyy_xyzz, tk_yyyyy_xyzzz, tk_yyyyy_xzzz, tk_yyyyy_xzzzz, tk_yyyyy_yyyy, tk_yyyyy_yyyyy, tk_yyyyy_yyyyz, tk_yyyyy_yyyz, tk_yyyyy_yyyzz, tk_yyyyy_yyzz, tk_yyyyy_yyzzz, tk_yyyyy_yzzz, tk_yyyyy_yzzzz, tk_yyyyy_zzzz, tk_yyyyy_zzzzz, tk_yyyyyz_xxxxx, tk_yyyyyz_xxxxy, tk_yyyyyz_xxxxz, tk_yyyyyz_xxxyy, tk_yyyyyz_xxxyz, tk_yyyyyz_xxxzz, tk_yyyyyz_xxyyy, tk_yyyyyz_xxyyz, tk_yyyyyz_xxyzz, tk_yyyyyz_xxzzz, tk_yyyyyz_xyyyy, tk_yyyyyz_xyyyz, tk_yyyyyz_xyyzz, tk_yyyyyz_xyzzz, tk_yyyyyz_xzzzz, tk_yyyyyz_yyyyy, tk_yyyyyz_yyyyz, tk_yyyyyz_yyyzz, tk_yyyyyz_yyzzz, tk_yyyyyz_yzzzz, tk_yyyyyz_zzzzz, ts_yyyyyz_xxxxx, ts_yyyyyz_xxxxy, ts_yyyyyz_xxxxz, ts_yyyyyz_xxxyy, ts_yyyyyz_xxxyz, ts_yyyyyz_xxxzz, ts_yyyyyz_xxyyy, ts_yyyyyz_xxyyz, ts_yyyyyz_xxyzz, ts_yyyyyz_xxzzz, ts_yyyyyz_xyyyy, ts_yyyyyz_xyyyz, ts_yyyyyz_xyyzz, ts_yyyyyz_xyzzz, ts_yyyyyz_xzzzz, ts_yyyyyz_yyyyy, ts_yyyyyz_yyyyz, ts_yyyyyz_yyyzz, ts_yyyyyz_yyzzz, ts_yyyyyz_yzzzz, ts_yyyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyyz_xxxxx[i] = tk_yyyyy_xxxxx[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxx[i] * fz_0;

        tk_yyyyyz_xxxxy[i] = tk_yyyyy_xxxxy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxy[i] * fz_0;

        tk_yyyyyz_xxxxz[i] = tk_yyyyy_xxxx[i] * fe_0 + tk_yyyyy_xxxxz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxxz[i] * fz_0;

        tk_yyyyyz_xxxyy[i] = tk_yyyyy_xxxyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxyy[i] * fz_0;

        tk_yyyyyz_xxxyz[i] = tk_yyyyy_xxxy[i] * fe_0 + tk_yyyyy_xxxyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxyz[i] * fz_0;

        tk_yyyyyz_xxxzz[i] = 2.0 * tk_yyyyy_xxxz[i] * fe_0 + tk_yyyyy_xxxzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxxzz[i] * fz_0;

        tk_yyyyyz_xxyyy[i] = tk_yyyyy_xxyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyyy[i] * fz_0;

        tk_yyyyyz_xxyyz[i] = tk_yyyyy_xxyy[i] * fe_0 + tk_yyyyy_xxyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyyz[i] * fz_0;

        tk_yyyyyz_xxyzz[i] = 2.0 * tk_yyyyy_xxyz[i] * fe_0 + tk_yyyyy_xxyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxyzz[i] * fz_0;

        tk_yyyyyz_xxzzz[i] = 3.0 * tk_yyyyy_xxzz[i] * fe_0 + tk_yyyyy_xxzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xxzzz[i] * fz_0;

        tk_yyyyyz_xyyyy[i] = tk_yyyyy_xyyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyyy[i] * fz_0;

        tk_yyyyyz_xyyyz[i] = tk_yyyyy_xyyy[i] * fe_0 + tk_yyyyy_xyyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyyz[i] * fz_0;

        tk_yyyyyz_xyyzz[i] = 2.0 * tk_yyyyy_xyyz[i] * fe_0 + tk_yyyyy_xyyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyyzz[i] * fz_0;

        tk_yyyyyz_xyzzz[i] = 3.0 * tk_yyyyy_xyzz[i] * fe_0 + tk_yyyyy_xyzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xyzzz[i] * fz_0;

        tk_yyyyyz_xzzzz[i] = 4.0 * tk_yyyyy_xzzz[i] * fe_0 + tk_yyyyy_xzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xzzzz[i] * fz_0;

        tk_yyyyyz_yyyyy[i] = tk_yyyyy_yyyyy[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyyy[i] * fz_0;

        tk_yyyyyz_yyyyz[i] = tk_yyyyy_yyyy[i] * fe_0 + tk_yyyyy_yyyyz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyyz[i] * fz_0;

        tk_yyyyyz_yyyzz[i] = 2.0 * tk_yyyyy_yyyz[i] * fe_0 + tk_yyyyy_yyyzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyyzz[i] * fz_0;

        tk_yyyyyz_yyzzz[i] = 3.0 * tk_yyyyy_yyzz[i] * fe_0 + tk_yyyyy_yyzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yyzzz[i] * fz_0;

        tk_yyyyyz_yzzzz[i] = 4.0 * tk_yyyyy_yzzz[i] * fe_0 + tk_yyyyy_yzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yzzzz[i] * fz_0;

        tk_yyyyyz_zzzzz[i] = 5.0 * tk_yyyyy_zzzz[i] * fe_0 + tk_yyyyy_zzzzz[i] * pa_z[i] + 2.0 * ts_yyyyyz_zzzzz[i] * fz_0;
    }

    // Set up 483-504 components of targeted buffer : IH

    auto tk_yyyyzz_xxxxx = pbuffer.data(idx_kin_ih + 483);

    auto tk_yyyyzz_xxxxy = pbuffer.data(idx_kin_ih + 484);

    auto tk_yyyyzz_xxxxz = pbuffer.data(idx_kin_ih + 485);

    auto tk_yyyyzz_xxxyy = pbuffer.data(idx_kin_ih + 486);

    auto tk_yyyyzz_xxxyz = pbuffer.data(idx_kin_ih + 487);

    auto tk_yyyyzz_xxxzz = pbuffer.data(idx_kin_ih + 488);

    auto tk_yyyyzz_xxyyy = pbuffer.data(idx_kin_ih + 489);

    auto tk_yyyyzz_xxyyz = pbuffer.data(idx_kin_ih + 490);

    auto tk_yyyyzz_xxyzz = pbuffer.data(idx_kin_ih + 491);

    auto tk_yyyyzz_xxzzz = pbuffer.data(idx_kin_ih + 492);

    auto tk_yyyyzz_xyyyy = pbuffer.data(idx_kin_ih + 493);

    auto tk_yyyyzz_xyyyz = pbuffer.data(idx_kin_ih + 494);

    auto tk_yyyyzz_xyyzz = pbuffer.data(idx_kin_ih + 495);

    auto tk_yyyyzz_xyzzz = pbuffer.data(idx_kin_ih + 496);

    auto tk_yyyyzz_xzzzz = pbuffer.data(idx_kin_ih + 497);

    auto tk_yyyyzz_yyyyy = pbuffer.data(idx_kin_ih + 498);

    auto tk_yyyyzz_yyyyz = pbuffer.data(idx_kin_ih + 499);

    auto tk_yyyyzz_yyyzz = pbuffer.data(idx_kin_ih + 500);

    auto tk_yyyyzz_yyzzz = pbuffer.data(idx_kin_ih + 501);

    auto tk_yyyyzz_yzzzz = pbuffer.data(idx_kin_ih + 502);

    auto tk_yyyyzz_zzzzz = pbuffer.data(idx_kin_ih + 503);

    #pragma omp simd aligned(pa_y, pa_z, tk_yyyy_xxxxy, tk_yyyy_xxxyy, tk_yyyy_xxyyy, tk_yyyy_xyyyy, tk_yyyy_yyyyy, tk_yyyyz_xxxxy, tk_yyyyz_xxxyy, tk_yyyyz_xxyyy, tk_yyyyz_xyyyy, tk_yyyyz_yyyyy, tk_yyyyzz_xxxxx, tk_yyyyzz_xxxxy, tk_yyyyzz_xxxxz, tk_yyyyzz_xxxyy, tk_yyyyzz_xxxyz, tk_yyyyzz_xxxzz, tk_yyyyzz_xxyyy, tk_yyyyzz_xxyyz, tk_yyyyzz_xxyzz, tk_yyyyzz_xxzzz, tk_yyyyzz_xyyyy, tk_yyyyzz_xyyyz, tk_yyyyzz_xyyzz, tk_yyyyzz_xyzzz, tk_yyyyzz_xzzzz, tk_yyyyzz_yyyyy, tk_yyyyzz_yyyyz, tk_yyyyzz_yyyzz, tk_yyyyzz_yyzzz, tk_yyyyzz_yzzzz, tk_yyyyzz_zzzzz, tk_yyyzz_xxxxx, tk_yyyzz_xxxxz, tk_yyyzz_xxxyz, tk_yyyzz_xxxz, tk_yyyzz_xxxzz, tk_yyyzz_xxyyz, tk_yyyzz_xxyz, tk_yyyzz_xxyzz, tk_yyyzz_xxzz, tk_yyyzz_xxzzz, tk_yyyzz_xyyyz, tk_yyyzz_xyyz, tk_yyyzz_xyyzz, tk_yyyzz_xyzz, tk_yyyzz_xyzzz, tk_yyyzz_xzzz, tk_yyyzz_xzzzz, tk_yyyzz_yyyyz, tk_yyyzz_yyyz, tk_yyyzz_yyyzz, tk_yyyzz_yyzz, tk_yyyzz_yyzzz, tk_yyyzz_yzzz, tk_yyyzz_yzzzz, tk_yyyzz_zzzz, tk_yyyzz_zzzzz, tk_yyzz_xxxxx, tk_yyzz_xxxxz, tk_yyzz_xxxyz, tk_yyzz_xxxzz, tk_yyzz_xxyyz, tk_yyzz_xxyzz, tk_yyzz_xxzzz, tk_yyzz_xyyyz, tk_yyzz_xyyzz, tk_yyzz_xyzzz, tk_yyzz_xzzzz, tk_yyzz_yyyyz, tk_yyzz_yyyzz, tk_yyzz_yyzzz, tk_yyzz_yzzzz, tk_yyzz_zzzzz, ts_yyyy_xxxxy, ts_yyyy_xxxyy, ts_yyyy_xxyyy, ts_yyyy_xyyyy, ts_yyyy_yyyyy, ts_yyyyzz_xxxxx, ts_yyyyzz_xxxxy, ts_yyyyzz_xxxxz, ts_yyyyzz_xxxyy, ts_yyyyzz_xxxyz, ts_yyyyzz_xxxzz, ts_yyyyzz_xxyyy, ts_yyyyzz_xxyyz, ts_yyyyzz_xxyzz, ts_yyyyzz_xxzzz, ts_yyyyzz_xyyyy, ts_yyyyzz_xyyyz, ts_yyyyzz_xyyzz, ts_yyyyzz_xyzzz, ts_yyyyzz_xzzzz, ts_yyyyzz_yyyyy, ts_yyyyzz_yyyyz, ts_yyyyzz_yyyzz, ts_yyyyzz_yyzzz, ts_yyyyzz_yzzzz, ts_yyyyzz_zzzzz, ts_yyzz_xxxxx, ts_yyzz_xxxxz, ts_yyzz_xxxyz, ts_yyzz_xxxzz, ts_yyzz_xxyyz, ts_yyzz_xxyzz, ts_yyzz_xxzzz, ts_yyzz_xyyyz, ts_yyzz_xyyzz, ts_yyzz_xyzzz, ts_yyzz_xzzzz, ts_yyzz_yyyyz, ts_yyzz_yyyzz, ts_yyzz_yyzzz, ts_yyzz_yzzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyzz_xxxxx[i] = -6.0 * ts_yyzz_xxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxx[i] * fe_0 + tk_yyyzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxxx[i] * fz_0;

        tk_yyyyzz_xxxxy[i] = -2.0 * ts_yyyy_xxxxy[i] * fbe_0 * fz_0 + tk_yyyy_xxxxy[i] * fe_0 + tk_yyyyz_xxxxy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxxxy[i] * fz_0;

        tk_yyyyzz_xxxxz[i] = -6.0 * ts_yyzz_xxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxz[i] * fe_0 + tk_yyyzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxxz[i] * fz_0;

        tk_yyyyzz_xxxyy[i] = -2.0 * ts_yyyy_xxxyy[i] * fbe_0 * fz_0 + tk_yyyy_xxxyy[i] * fe_0 + tk_yyyyz_xxxyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxxyy[i] * fz_0;

        tk_yyyyzz_xxxyz[i] = -6.0 * ts_yyzz_xxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxyz[i] * fe_0 + tk_yyyzz_xxxz[i] * fe_0 + tk_yyyzz_xxxyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxyz[i] * fz_0;

        tk_yyyyzz_xxxzz[i] = -6.0 * ts_yyzz_xxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxzz[i] * fe_0 + tk_yyyzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxxzz[i] * fz_0;

        tk_yyyyzz_xxyyy[i] = -2.0 * ts_yyyy_xxyyy[i] * fbe_0 * fz_0 + tk_yyyy_xxyyy[i] * fe_0 + tk_yyyyz_xxyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xxyyy[i] * fz_0;

        tk_yyyyzz_xxyyz[i] = -6.0 * ts_yyzz_xxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyyz[i] * fe_0 + 2.0 * tk_yyyzz_xxyz[i] * fe_0 + tk_yyyzz_xxyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxyyz[i] * fz_0;

        tk_yyyyzz_xxyzz[i] = -6.0 * ts_yyzz_xxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyzz[i] * fe_0 + tk_yyyzz_xxzz[i] * fe_0 + tk_yyyzz_xxyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxyzz[i] * fz_0;

        tk_yyyyzz_xxzzz[i] = -6.0 * ts_yyzz_xxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxzzz[i] * fe_0 + tk_yyyzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xxzzz[i] * fz_0;

        tk_yyyyzz_xyyyy[i] = -2.0 * ts_yyyy_xyyyy[i] * fbe_0 * fz_0 + tk_yyyy_xyyyy[i] * fe_0 + tk_yyyyz_xyyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xyyyy[i] * fz_0;

        tk_yyyyzz_xyyyz[i] = -6.0 * ts_yyzz_xyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyyz[i] * fe_0 + 3.0 * tk_yyyzz_xyyz[i] * fe_0 + tk_yyyzz_xyyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyyyz[i] * fz_0;

        tk_yyyyzz_xyyzz[i] = -6.0 * ts_yyzz_xyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyzz[i] * fe_0 + 2.0 * tk_yyyzz_xyzz[i] * fe_0 + tk_yyyzz_xyyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyyzz[i] * fz_0;

        tk_yyyyzz_xyzzz[i] = -6.0 * ts_yyzz_xyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyzzz[i] * fe_0 + tk_yyyzz_xzzz[i] * fe_0 + tk_yyyzz_xyzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xyzzz[i] * fz_0;

        tk_yyyyzz_xzzzz[i] = -6.0 * ts_yyzz_xzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xzzzz[i] * fe_0 + tk_yyyzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xzzzz[i] * fz_0;

        tk_yyyyzz_yyyyy[i] = -2.0 * ts_yyyy_yyyyy[i] * fbe_0 * fz_0 + tk_yyyy_yyyyy[i] * fe_0 + tk_yyyyz_yyyyy[i] * pa_z[i] + 2.0 * ts_yyyyzz_yyyyy[i] * fz_0;

        tk_yyyyzz_yyyyz[i] = -6.0 * ts_yyzz_yyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyyz[i] * fe_0 + 4.0 * tk_yyyzz_yyyz[i] * fe_0 + tk_yyyzz_yyyyz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyyyz[i] * fz_0;

        tk_yyyyzz_yyyzz[i] = -6.0 * ts_yyzz_yyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyzz[i] * fe_0 + 3.0 * tk_yyyzz_yyzz[i] * fe_0 + tk_yyyzz_yyyzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyyzz[i] * fz_0;

        tk_yyyyzz_yyzzz[i] = -6.0 * ts_yyzz_yyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyzzz[i] * fe_0 + 2.0 * tk_yyyzz_yzzz[i] * fe_0 + tk_yyyzz_yyzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yyzzz[i] * fz_0;

        tk_yyyyzz_yzzzz[i] = -6.0 * ts_yyzz_yzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yzzzz[i] * fe_0 + tk_yyyzz_zzzz[i] * fe_0 + tk_yyyzz_yzzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_yzzzz[i] * fz_0;

        tk_yyyyzz_zzzzz[i] = -6.0 * ts_yyzz_zzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_zzzzz[i] * fe_0 + tk_yyyzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yyyyzz_zzzzz[i] * fz_0;
    }

    // Set up 504-525 components of targeted buffer : IH

    auto tk_yyyzzz_xxxxx = pbuffer.data(idx_kin_ih + 504);

    auto tk_yyyzzz_xxxxy = pbuffer.data(idx_kin_ih + 505);

    auto tk_yyyzzz_xxxxz = pbuffer.data(idx_kin_ih + 506);

    auto tk_yyyzzz_xxxyy = pbuffer.data(idx_kin_ih + 507);

    auto tk_yyyzzz_xxxyz = pbuffer.data(idx_kin_ih + 508);

    auto tk_yyyzzz_xxxzz = pbuffer.data(idx_kin_ih + 509);

    auto tk_yyyzzz_xxyyy = pbuffer.data(idx_kin_ih + 510);

    auto tk_yyyzzz_xxyyz = pbuffer.data(idx_kin_ih + 511);

    auto tk_yyyzzz_xxyzz = pbuffer.data(idx_kin_ih + 512);

    auto tk_yyyzzz_xxzzz = pbuffer.data(idx_kin_ih + 513);

    auto tk_yyyzzz_xyyyy = pbuffer.data(idx_kin_ih + 514);

    auto tk_yyyzzz_xyyyz = pbuffer.data(idx_kin_ih + 515);

    auto tk_yyyzzz_xyyzz = pbuffer.data(idx_kin_ih + 516);

    auto tk_yyyzzz_xyzzz = pbuffer.data(idx_kin_ih + 517);

    auto tk_yyyzzz_xzzzz = pbuffer.data(idx_kin_ih + 518);

    auto tk_yyyzzz_yyyyy = pbuffer.data(idx_kin_ih + 519);

    auto tk_yyyzzz_yyyyz = pbuffer.data(idx_kin_ih + 520);

    auto tk_yyyzzz_yyyzz = pbuffer.data(idx_kin_ih + 521);

    auto tk_yyyzzz_yyzzz = pbuffer.data(idx_kin_ih + 522);

    auto tk_yyyzzz_yzzzz = pbuffer.data(idx_kin_ih + 523);

    auto tk_yyyzzz_zzzzz = pbuffer.data(idx_kin_ih + 524);

    #pragma omp simd aligned(pa_y, pa_z, tk_yyyz_xxxxy, tk_yyyz_xxxyy, tk_yyyz_xxyyy, tk_yyyz_xyyyy, tk_yyyz_yyyyy, tk_yyyzz_xxxxy, tk_yyyzz_xxxyy, tk_yyyzz_xxyyy, tk_yyyzz_xyyyy, tk_yyyzz_yyyyy, tk_yyyzzz_xxxxx, tk_yyyzzz_xxxxy, tk_yyyzzz_xxxxz, tk_yyyzzz_xxxyy, tk_yyyzzz_xxxyz, tk_yyyzzz_xxxzz, tk_yyyzzz_xxyyy, tk_yyyzzz_xxyyz, tk_yyyzzz_xxyzz, tk_yyyzzz_xxzzz, tk_yyyzzz_xyyyy, tk_yyyzzz_xyyyz, tk_yyyzzz_xyyzz, tk_yyyzzz_xyzzz, tk_yyyzzz_xzzzz, tk_yyyzzz_yyyyy, tk_yyyzzz_yyyyz, tk_yyyzzz_yyyzz, tk_yyyzzz_yyzzz, tk_yyyzzz_yzzzz, tk_yyyzzz_zzzzz, tk_yyzzz_xxxxx, tk_yyzzz_xxxxz, tk_yyzzz_xxxyz, tk_yyzzz_xxxz, tk_yyzzz_xxxzz, tk_yyzzz_xxyyz, tk_yyzzz_xxyz, tk_yyzzz_xxyzz, tk_yyzzz_xxzz, tk_yyzzz_xxzzz, tk_yyzzz_xyyyz, tk_yyzzz_xyyz, tk_yyzzz_xyyzz, tk_yyzzz_xyzz, tk_yyzzz_xyzzz, tk_yyzzz_xzzz, tk_yyzzz_xzzzz, tk_yyzzz_yyyyz, tk_yyzzz_yyyz, tk_yyzzz_yyyzz, tk_yyzzz_yyzz, tk_yyzzz_yyzzz, tk_yyzzz_yzzz, tk_yyzzz_yzzzz, tk_yyzzz_zzzz, tk_yyzzz_zzzzz, tk_yzzz_xxxxx, tk_yzzz_xxxxz, tk_yzzz_xxxyz, tk_yzzz_xxxzz, tk_yzzz_xxyyz, tk_yzzz_xxyzz, tk_yzzz_xxzzz, tk_yzzz_xyyyz, tk_yzzz_xyyzz, tk_yzzz_xyzzz, tk_yzzz_xzzzz, tk_yzzz_yyyyz, tk_yzzz_yyyzz, tk_yzzz_yyzzz, tk_yzzz_yzzzz, tk_yzzz_zzzzz, ts_yyyz_xxxxy, ts_yyyz_xxxyy, ts_yyyz_xxyyy, ts_yyyz_xyyyy, ts_yyyz_yyyyy, ts_yyyzzz_xxxxx, ts_yyyzzz_xxxxy, ts_yyyzzz_xxxxz, ts_yyyzzz_xxxyy, ts_yyyzzz_xxxyz, ts_yyyzzz_xxxzz, ts_yyyzzz_xxyyy, ts_yyyzzz_xxyyz, ts_yyyzzz_xxyzz, ts_yyyzzz_xxzzz, ts_yyyzzz_xyyyy, ts_yyyzzz_xyyyz, ts_yyyzzz_xyyzz, ts_yyyzzz_xyzzz, ts_yyyzzz_xzzzz, ts_yyyzzz_yyyyy, ts_yyyzzz_yyyyz, ts_yyyzzz_yyyzz, ts_yyyzzz_yyzzz, ts_yyyzzz_yzzzz, ts_yyyzzz_zzzzz, ts_yzzz_xxxxx, ts_yzzz_xxxxz, ts_yzzz_xxxyz, ts_yzzz_xxxzz, ts_yzzz_xxyyz, ts_yzzz_xxyzz, ts_yzzz_xxzzz, ts_yzzz_xyyyz, ts_yzzz_xyyzz, ts_yzzz_xyzzz, ts_yzzz_xzzzz, ts_yzzz_yyyyz, ts_yzzz_yyyzz, ts_yzzz_yyzzz, ts_yzzz_yzzzz, ts_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzzz_xxxxx[i] = -4.0 * ts_yzzz_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxxx[i] * fe_0 + tk_yyzzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxxx[i] * fz_0;

        tk_yyyzzz_xxxxy[i] = -4.0 * ts_yyyz_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxxxy[i] * fe_0 + tk_yyyzz_xxxxy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xxxxy[i] * fz_0;

        tk_yyyzzz_xxxxz[i] = -4.0 * ts_yzzz_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxxz[i] * fe_0 + tk_yyzzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxxz[i] * fz_0;

        tk_yyyzzz_xxxyy[i] = -4.0 * ts_yyyz_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxxyy[i] * fe_0 + tk_yyyzz_xxxyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xxxyy[i] * fz_0;

        tk_yyyzzz_xxxyz[i] = -4.0 * ts_yzzz_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxyz[i] * fe_0 + tk_yyzzz_xxxz[i] * fe_0 + tk_yyzzz_xxxyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxyz[i] * fz_0;

        tk_yyyzzz_xxxzz[i] = -4.0 * ts_yzzz_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxxzz[i] * fe_0 + tk_yyzzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxxzz[i] * fz_0;

        tk_yyyzzz_xxyyy[i] = -4.0 * ts_yyyz_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xxyyy[i] * fe_0 + tk_yyyzz_xxyyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xxyyy[i] * fz_0;

        tk_yyyzzz_xxyyz[i] = -4.0 * ts_yzzz_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxyyz[i] * fe_0 + 2.0 * tk_yyzzz_xxyz[i] * fe_0 + tk_yyzzz_xxyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxyyz[i] * fz_0;

        tk_yyyzzz_xxyzz[i] = -4.0 * ts_yzzz_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxyzz[i] * fe_0 + tk_yyzzz_xxzz[i] * fe_0 + tk_yyzzz_xxyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxyzz[i] * fz_0;

        tk_yyyzzz_xxzzz[i] = -4.0 * ts_yzzz_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xxzzz[i] * fe_0 + tk_yyzzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xxzzz[i] * fz_0;

        tk_yyyzzz_xyyyy[i] = -4.0 * ts_yyyz_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xyyyy[i] * fe_0 + tk_yyyzz_xyyyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xyyyy[i] * fz_0;

        tk_yyyzzz_xyyyz[i] = -4.0 * ts_yzzz_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyyyz[i] * fe_0 + 3.0 * tk_yyzzz_xyyz[i] * fe_0 + tk_yyzzz_xyyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyyyz[i] * fz_0;

        tk_yyyzzz_xyyzz[i] = -4.0 * ts_yzzz_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyyzz[i] * fe_0 + 2.0 * tk_yyzzz_xyzz[i] * fe_0 + tk_yyzzz_xyyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyyzz[i] * fz_0;

        tk_yyyzzz_xyzzz[i] = -4.0 * ts_yzzz_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xyzzz[i] * fe_0 + tk_yyzzz_xzzz[i] * fe_0 + tk_yyzzz_xyzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xyzzz[i] * fz_0;

        tk_yyyzzz_xzzzz[i] = -4.0 * ts_yzzz_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xzzzz[i] * fe_0 + tk_yyzzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xzzzz[i] * fz_0;

        tk_yyyzzz_yyyyy[i] = -4.0 * ts_yyyz_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_yyyyy[i] * fe_0 + tk_yyyzz_yyyyy[i] * pa_z[i] + 2.0 * ts_yyyzzz_yyyyy[i] * fz_0;

        tk_yyyzzz_yyyyz[i] = -4.0 * ts_yzzz_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyyyz[i] * fe_0 + 4.0 * tk_yyzzz_yyyz[i] * fe_0 + tk_yyzzz_yyyyz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyyyz[i] * fz_0;

        tk_yyyzzz_yyyzz[i] = -4.0 * ts_yzzz_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyyzz[i] * fe_0 + 3.0 * tk_yyzzz_yyzz[i] * fe_0 + tk_yyzzz_yyyzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyyzz[i] * fz_0;

        tk_yyyzzz_yyzzz[i] = -4.0 * ts_yzzz_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yyzzz[i] * fe_0 + 2.0 * tk_yyzzz_yzzz[i] * fe_0 + tk_yyzzz_yyzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yyzzz[i] * fz_0;

        tk_yyyzzz_yzzzz[i] = -4.0 * ts_yzzz_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yzzzz[i] * fe_0 + tk_yyzzz_zzzz[i] * fe_0 + tk_yyzzz_yzzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_yzzzz[i] * fz_0;

        tk_yyyzzz_zzzzz[i] = -4.0 * ts_yzzz_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_zzzzz[i] * fe_0 + tk_yyzzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yyyzzz_zzzzz[i] * fz_0;
    }

    // Set up 525-546 components of targeted buffer : IH

    auto tk_yyzzzz_xxxxx = pbuffer.data(idx_kin_ih + 525);

    auto tk_yyzzzz_xxxxy = pbuffer.data(idx_kin_ih + 526);

    auto tk_yyzzzz_xxxxz = pbuffer.data(idx_kin_ih + 527);

    auto tk_yyzzzz_xxxyy = pbuffer.data(idx_kin_ih + 528);

    auto tk_yyzzzz_xxxyz = pbuffer.data(idx_kin_ih + 529);

    auto tk_yyzzzz_xxxzz = pbuffer.data(idx_kin_ih + 530);

    auto tk_yyzzzz_xxyyy = pbuffer.data(idx_kin_ih + 531);

    auto tk_yyzzzz_xxyyz = pbuffer.data(idx_kin_ih + 532);

    auto tk_yyzzzz_xxyzz = pbuffer.data(idx_kin_ih + 533);

    auto tk_yyzzzz_xxzzz = pbuffer.data(idx_kin_ih + 534);

    auto tk_yyzzzz_xyyyy = pbuffer.data(idx_kin_ih + 535);

    auto tk_yyzzzz_xyyyz = pbuffer.data(idx_kin_ih + 536);

    auto tk_yyzzzz_xyyzz = pbuffer.data(idx_kin_ih + 537);

    auto tk_yyzzzz_xyzzz = pbuffer.data(idx_kin_ih + 538);

    auto tk_yyzzzz_xzzzz = pbuffer.data(idx_kin_ih + 539);

    auto tk_yyzzzz_yyyyy = pbuffer.data(idx_kin_ih + 540);

    auto tk_yyzzzz_yyyyz = pbuffer.data(idx_kin_ih + 541);

    auto tk_yyzzzz_yyyzz = pbuffer.data(idx_kin_ih + 542);

    auto tk_yyzzzz_yyzzz = pbuffer.data(idx_kin_ih + 543);

    auto tk_yyzzzz_yzzzz = pbuffer.data(idx_kin_ih + 544);

    auto tk_yyzzzz_zzzzz = pbuffer.data(idx_kin_ih + 545);

    #pragma omp simd aligned(pa_y, pa_z, tk_yyzz_xxxxy, tk_yyzz_xxxyy, tk_yyzz_xxyyy, tk_yyzz_xyyyy, tk_yyzz_yyyyy, tk_yyzzz_xxxxy, tk_yyzzz_xxxyy, tk_yyzzz_xxyyy, tk_yyzzz_xyyyy, tk_yyzzz_yyyyy, tk_yyzzzz_xxxxx, tk_yyzzzz_xxxxy, tk_yyzzzz_xxxxz, tk_yyzzzz_xxxyy, tk_yyzzzz_xxxyz, tk_yyzzzz_xxxzz, tk_yyzzzz_xxyyy, tk_yyzzzz_xxyyz, tk_yyzzzz_xxyzz, tk_yyzzzz_xxzzz, tk_yyzzzz_xyyyy, tk_yyzzzz_xyyyz, tk_yyzzzz_xyyzz, tk_yyzzzz_xyzzz, tk_yyzzzz_xzzzz, tk_yyzzzz_yyyyy, tk_yyzzzz_yyyyz, tk_yyzzzz_yyyzz, tk_yyzzzz_yyzzz, tk_yyzzzz_yzzzz, tk_yyzzzz_zzzzz, tk_yzzzz_xxxxx, tk_yzzzz_xxxxz, tk_yzzzz_xxxyz, tk_yzzzz_xxxz, tk_yzzzz_xxxzz, tk_yzzzz_xxyyz, tk_yzzzz_xxyz, tk_yzzzz_xxyzz, tk_yzzzz_xxzz, tk_yzzzz_xxzzz, tk_yzzzz_xyyyz, tk_yzzzz_xyyz, tk_yzzzz_xyyzz, tk_yzzzz_xyzz, tk_yzzzz_xyzzz, tk_yzzzz_xzzz, tk_yzzzz_xzzzz, tk_yzzzz_yyyyz, tk_yzzzz_yyyz, tk_yzzzz_yyyzz, tk_yzzzz_yyzz, tk_yzzzz_yyzzz, tk_yzzzz_yzzz, tk_yzzzz_yzzzz, tk_yzzzz_zzzz, tk_yzzzz_zzzzz, tk_zzzz_xxxxx, tk_zzzz_xxxxz, tk_zzzz_xxxyz, tk_zzzz_xxxzz, tk_zzzz_xxyyz, tk_zzzz_xxyzz, tk_zzzz_xxzzz, tk_zzzz_xyyyz, tk_zzzz_xyyzz, tk_zzzz_xyzzz, tk_zzzz_xzzzz, tk_zzzz_yyyyz, tk_zzzz_yyyzz, tk_zzzz_yyzzz, tk_zzzz_yzzzz, tk_zzzz_zzzzz, ts_yyzz_xxxxy, ts_yyzz_xxxyy, ts_yyzz_xxyyy, ts_yyzz_xyyyy, ts_yyzz_yyyyy, ts_yyzzzz_xxxxx, ts_yyzzzz_xxxxy, ts_yyzzzz_xxxxz, ts_yyzzzz_xxxyy, ts_yyzzzz_xxxyz, ts_yyzzzz_xxxzz, ts_yyzzzz_xxyyy, ts_yyzzzz_xxyyz, ts_yyzzzz_xxyzz, ts_yyzzzz_xxzzz, ts_yyzzzz_xyyyy, ts_yyzzzz_xyyyz, ts_yyzzzz_xyyzz, ts_yyzzzz_xyzzz, ts_yyzzzz_xzzzz, ts_yyzzzz_yyyyy, ts_yyzzzz_yyyyz, ts_yyzzzz_yyyzz, ts_yyzzzz_yyzzz, ts_yyzzzz_yzzzz, ts_yyzzzz_zzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxz, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzzz_xxxxx[i] = -2.0 * ts_zzzz_xxxxx[i] * fbe_0 * fz_0 + tk_zzzz_xxxxx[i] * fe_0 + tk_yzzzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxxx[i] * fz_0;

        tk_yyzzzz_xxxxy[i] = -6.0 * ts_yyzz_xxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxxy[i] * fe_0 + tk_yyzzz_xxxxy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xxxxy[i] * fz_0;

        tk_yyzzzz_xxxxz[i] = -2.0 * ts_zzzz_xxxxz[i] * fbe_0 * fz_0 + tk_zzzz_xxxxz[i] * fe_0 + tk_yzzzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxxz[i] * fz_0;

        tk_yyzzzz_xxxyy[i] = -6.0 * ts_yyzz_xxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxxyy[i] * fe_0 + tk_yyzzz_xxxyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xxxyy[i] * fz_0;

        tk_yyzzzz_xxxyz[i] = -2.0 * ts_zzzz_xxxyz[i] * fbe_0 * fz_0 + tk_zzzz_xxxyz[i] * fe_0 + tk_yzzzz_xxxz[i] * fe_0 + tk_yzzzz_xxxyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxyz[i] * fz_0;

        tk_yyzzzz_xxxzz[i] = -2.0 * ts_zzzz_xxxzz[i] * fbe_0 * fz_0 + tk_zzzz_xxxzz[i] * fe_0 + tk_yzzzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxxzz[i] * fz_0;

        tk_yyzzzz_xxyyy[i] = -6.0 * ts_yyzz_xxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xxyyy[i] * fe_0 + tk_yyzzz_xxyyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xxyyy[i] * fz_0;

        tk_yyzzzz_xxyyz[i] = -2.0 * ts_zzzz_xxyyz[i] * fbe_0 * fz_0 + tk_zzzz_xxyyz[i] * fe_0 + 2.0 * tk_yzzzz_xxyz[i] * fe_0 + tk_yzzzz_xxyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxyyz[i] * fz_0;

        tk_yyzzzz_xxyzz[i] = -2.0 * ts_zzzz_xxyzz[i] * fbe_0 * fz_0 + tk_zzzz_xxyzz[i] * fe_0 + tk_yzzzz_xxzz[i] * fe_0 + tk_yzzzz_xxyzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxyzz[i] * fz_0;

        tk_yyzzzz_xxzzz[i] = -2.0 * ts_zzzz_xxzzz[i] * fbe_0 * fz_0 + tk_zzzz_xxzzz[i] * fe_0 + tk_yzzzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xxzzz[i] * fz_0;

        tk_yyzzzz_xyyyy[i] = -6.0 * ts_yyzz_xyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xyyyy[i] * fe_0 + tk_yyzzz_xyyyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xyyyy[i] * fz_0;

        tk_yyzzzz_xyyyz[i] = -2.0 * ts_zzzz_xyyyz[i] * fbe_0 * fz_0 + tk_zzzz_xyyyz[i] * fe_0 + 3.0 * tk_yzzzz_xyyz[i] * fe_0 + tk_yzzzz_xyyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xyyyz[i] * fz_0;

        tk_yyzzzz_xyyzz[i] = -2.0 * ts_zzzz_xyyzz[i] * fbe_0 * fz_0 + tk_zzzz_xyyzz[i] * fe_0 + 2.0 * tk_yzzzz_xyzz[i] * fe_0 + tk_yzzzz_xyyzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xyyzz[i] * fz_0;

        tk_yyzzzz_xyzzz[i] = -2.0 * ts_zzzz_xyzzz[i] * fbe_0 * fz_0 + tk_zzzz_xyzzz[i] * fe_0 + tk_yzzzz_xzzz[i] * fe_0 + tk_yzzzz_xyzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xyzzz[i] * fz_0;

        tk_yyzzzz_xzzzz[i] = -2.0 * ts_zzzz_xzzzz[i] * fbe_0 * fz_0 + tk_zzzz_xzzzz[i] * fe_0 + tk_yzzzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xzzzz[i] * fz_0;

        tk_yyzzzz_yyyyy[i] = -6.0 * ts_yyzz_yyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yyyyy[i] * fe_0 + tk_yyzzz_yyyyy[i] * pa_z[i] + 2.0 * ts_yyzzzz_yyyyy[i] * fz_0;

        tk_yyzzzz_yyyyz[i] = -2.0 * ts_zzzz_yyyyz[i] * fbe_0 * fz_0 + tk_zzzz_yyyyz[i] * fe_0 + 4.0 * tk_yzzzz_yyyz[i] * fe_0 + tk_yzzzz_yyyyz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyyyz[i] * fz_0;

        tk_yyzzzz_yyyzz[i] = -2.0 * ts_zzzz_yyyzz[i] * fbe_0 * fz_0 + tk_zzzz_yyyzz[i] * fe_0 + 3.0 * tk_yzzzz_yyzz[i] * fe_0 + tk_yzzzz_yyyzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyyzz[i] * fz_0;

        tk_yyzzzz_yyzzz[i] = -2.0 * ts_zzzz_yyzzz[i] * fbe_0 * fz_0 + tk_zzzz_yyzzz[i] * fe_0 + 2.0 * tk_yzzzz_yzzz[i] * fe_0 + tk_yzzzz_yyzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yyzzz[i] * fz_0;

        tk_yyzzzz_yzzzz[i] = -2.0 * ts_zzzz_yzzzz[i] * fbe_0 * fz_0 + tk_zzzz_yzzzz[i] * fe_0 + tk_yzzzz_zzzz[i] * fe_0 + tk_yzzzz_yzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_yzzzz[i] * fz_0;

        tk_yyzzzz_zzzzz[i] = -2.0 * ts_zzzz_zzzzz[i] * fbe_0 * fz_0 + tk_zzzz_zzzzz[i] * fe_0 + tk_yzzzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yyzzzz_zzzzz[i] * fz_0;
    }

    // Set up 546-567 components of targeted buffer : IH

    auto tk_yzzzzz_xxxxx = pbuffer.data(idx_kin_ih + 546);

    auto tk_yzzzzz_xxxxy = pbuffer.data(idx_kin_ih + 547);

    auto tk_yzzzzz_xxxxz = pbuffer.data(idx_kin_ih + 548);

    auto tk_yzzzzz_xxxyy = pbuffer.data(idx_kin_ih + 549);

    auto tk_yzzzzz_xxxyz = pbuffer.data(idx_kin_ih + 550);

    auto tk_yzzzzz_xxxzz = pbuffer.data(idx_kin_ih + 551);

    auto tk_yzzzzz_xxyyy = pbuffer.data(idx_kin_ih + 552);

    auto tk_yzzzzz_xxyyz = pbuffer.data(idx_kin_ih + 553);

    auto tk_yzzzzz_xxyzz = pbuffer.data(idx_kin_ih + 554);

    auto tk_yzzzzz_xxzzz = pbuffer.data(idx_kin_ih + 555);

    auto tk_yzzzzz_xyyyy = pbuffer.data(idx_kin_ih + 556);

    auto tk_yzzzzz_xyyyz = pbuffer.data(idx_kin_ih + 557);

    auto tk_yzzzzz_xyyzz = pbuffer.data(idx_kin_ih + 558);

    auto tk_yzzzzz_xyzzz = pbuffer.data(idx_kin_ih + 559);

    auto tk_yzzzzz_xzzzz = pbuffer.data(idx_kin_ih + 560);

    auto tk_yzzzzz_yyyyy = pbuffer.data(idx_kin_ih + 561);

    auto tk_yzzzzz_yyyyz = pbuffer.data(idx_kin_ih + 562);

    auto tk_yzzzzz_yyyzz = pbuffer.data(idx_kin_ih + 563);

    auto tk_yzzzzz_yyzzz = pbuffer.data(idx_kin_ih + 564);

    auto tk_yzzzzz_yzzzz = pbuffer.data(idx_kin_ih + 565);

    auto tk_yzzzzz_zzzzz = pbuffer.data(idx_kin_ih + 566);

    #pragma omp simd aligned(pa_y, tk_yzzzzz_xxxxx, tk_yzzzzz_xxxxy, tk_yzzzzz_xxxxz, tk_yzzzzz_xxxyy, tk_yzzzzz_xxxyz, tk_yzzzzz_xxxzz, tk_yzzzzz_xxyyy, tk_yzzzzz_xxyyz, tk_yzzzzz_xxyzz, tk_yzzzzz_xxzzz, tk_yzzzzz_xyyyy, tk_yzzzzz_xyyyz, tk_yzzzzz_xyyzz, tk_yzzzzz_xyzzz, tk_yzzzzz_xzzzz, tk_yzzzzz_yyyyy, tk_yzzzzz_yyyyz, tk_yzzzzz_yyyzz, tk_yzzzzz_yyzzz, tk_yzzzzz_yzzzz, tk_yzzzzz_zzzzz, tk_zzzzz_xxxx, tk_zzzzz_xxxxx, tk_zzzzz_xxxxy, tk_zzzzz_xxxxz, tk_zzzzz_xxxy, tk_zzzzz_xxxyy, tk_zzzzz_xxxyz, tk_zzzzz_xxxz, tk_zzzzz_xxxzz, tk_zzzzz_xxyy, tk_zzzzz_xxyyy, tk_zzzzz_xxyyz, tk_zzzzz_xxyz, tk_zzzzz_xxyzz, tk_zzzzz_xxzz, tk_zzzzz_xxzzz, tk_zzzzz_xyyy, tk_zzzzz_xyyyy, tk_zzzzz_xyyyz, tk_zzzzz_xyyz, tk_zzzzz_xyyzz, tk_zzzzz_xyzz, tk_zzzzz_xyzzz, tk_zzzzz_xzzz, tk_zzzzz_xzzzz, tk_zzzzz_yyyy, tk_zzzzz_yyyyy, tk_zzzzz_yyyyz, tk_zzzzz_yyyz, tk_zzzzz_yyyzz, tk_zzzzz_yyzz, tk_zzzzz_yyzzz, tk_zzzzz_yzzz, tk_zzzzz_yzzzz, tk_zzzzz_zzzz, tk_zzzzz_zzzzz, ts_yzzzzz_xxxxx, ts_yzzzzz_xxxxy, ts_yzzzzz_xxxxz, ts_yzzzzz_xxxyy, ts_yzzzzz_xxxyz, ts_yzzzzz_xxxzz, ts_yzzzzz_xxyyy, ts_yzzzzz_xxyyz, ts_yzzzzz_xxyzz, ts_yzzzzz_xxzzz, ts_yzzzzz_xyyyy, ts_yzzzzz_xyyyz, ts_yzzzzz_xyyzz, ts_yzzzzz_xyzzz, ts_yzzzzz_xzzzz, ts_yzzzzz_yyyyy, ts_yzzzzz_yyyyz, ts_yzzzzz_yyyzz, ts_yzzzzz_yyzzz, ts_yzzzzz_yzzzz, ts_yzzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzzz_xxxxx[i] = tk_zzzzz_xxxxx[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxx[i] * fz_0;

        tk_yzzzzz_xxxxy[i] = tk_zzzzz_xxxx[i] * fe_0 + tk_zzzzz_xxxxy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxy[i] * fz_0;

        tk_yzzzzz_xxxxz[i] = tk_zzzzz_xxxxz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxxz[i] * fz_0;

        tk_yzzzzz_xxxyy[i] = 2.0 * tk_zzzzz_xxxy[i] * fe_0 + tk_zzzzz_xxxyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxyy[i] * fz_0;

        tk_yzzzzz_xxxyz[i] = tk_zzzzz_xxxz[i] * fe_0 + tk_zzzzz_xxxyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxyz[i] * fz_0;

        tk_yzzzzz_xxxzz[i] = tk_zzzzz_xxxzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxxzz[i] * fz_0;

        tk_yzzzzz_xxyyy[i] = 3.0 * tk_zzzzz_xxyy[i] * fe_0 + tk_zzzzz_xxyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyyy[i] * fz_0;

        tk_yzzzzz_xxyyz[i] = 2.0 * tk_zzzzz_xxyz[i] * fe_0 + tk_zzzzz_xxyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyyz[i] * fz_0;

        tk_yzzzzz_xxyzz[i] = tk_zzzzz_xxzz[i] * fe_0 + tk_zzzzz_xxyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxyzz[i] * fz_0;

        tk_yzzzzz_xxzzz[i] = tk_zzzzz_xxzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xxzzz[i] * fz_0;

        tk_yzzzzz_xyyyy[i] = 4.0 * tk_zzzzz_xyyy[i] * fe_0 + tk_zzzzz_xyyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyyy[i] * fz_0;

        tk_yzzzzz_xyyyz[i] = 3.0 * tk_zzzzz_xyyz[i] * fe_0 + tk_zzzzz_xyyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyyz[i] * fz_0;

        tk_yzzzzz_xyyzz[i] = 2.0 * tk_zzzzz_xyzz[i] * fe_0 + tk_zzzzz_xyyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyyzz[i] * fz_0;

        tk_yzzzzz_xyzzz[i] = tk_zzzzz_xzzz[i] * fe_0 + tk_zzzzz_xyzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xyzzz[i] * fz_0;

        tk_yzzzzz_xzzzz[i] = tk_zzzzz_xzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xzzzz[i] * fz_0;

        tk_yzzzzz_yyyyy[i] = 5.0 * tk_zzzzz_yyyy[i] * fe_0 + tk_zzzzz_yyyyy[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyyy[i] * fz_0;

        tk_yzzzzz_yyyyz[i] = 4.0 * tk_zzzzz_yyyz[i] * fe_0 + tk_zzzzz_yyyyz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyyz[i] * fz_0;

        tk_yzzzzz_yyyzz[i] = 3.0 * tk_zzzzz_yyzz[i] * fe_0 + tk_zzzzz_yyyzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyyzz[i] * fz_0;

        tk_yzzzzz_yyzzz[i] = 2.0 * tk_zzzzz_yzzz[i] * fe_0 + tk_zzzzz_yyzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yyzzz[i] * fz_0;

        tk_yzzzzz_yzzzz[i] = tk_zzzzz_zzzz[i] * fe_0 + tk_zzzzz_yzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yzzzz[i] * fz_0;

        tk_yzzzzz_zzzzz[i] = tk_zzzzz_zzzzz[i] * pa_y[i] + 2.0 * ts_yzzzzz_zzzzz[i] * fz_0;
    }

    // Set up 567-588 components of targeted buffer : IH

    auto tk_zzzzzz_xxxxx = pbuffer.data(idx_kin_ih + 567);

    auto tk_zzzzzz_xxxxy = pbuffer.data(idx_kin_ih + 568);

    auto tk_zzzzzz_xxxxz = pbuffer.data(idx_kin_ih + 569);

    auto tk_zzzzzz_xxxyy = pbuffer.data(idx_kin_ih + 570);

    auto tk_zzzzzz_xxxyz = pbuffer.data(idx_kin_ih + 571);

    auto tk_zzzzzz_xxxzz = pbuffer.data(idx_kin_ih + 572);

    auto tk_zzzzzz_xxyyy = pbuffer.data(idx_kin_ih + 573);

    auto tk_zzzzzz_xxyyz = pbuffer.data(idx_kin_ih + 574);

    auto tk_zzzzzz_xxyzz = pbuffer.data(idx_kin_ih + 575);

    auto tk_zzzzzz_xxzzz = pbuffer.data(idx_kin_ih + 576);

    auto tk_zzzzzz_xyyyy = pbuffer.data(idx_kin_ih + 577);

    auto tk_zzzzzz_xyyyz = pbuffer.data(idx_kin_ih + 578);

    auto tk_zzzzzz_xyyzz = pbuffer.data(idx_kin_ih + 579);

    auto tk_zzzzzz_xyzzz = pbuffer.data(idx_kin_ih + 580);

    auto tk_zzzzzz_xzzzz = pbuffer.data(idx_kin_ih + 581);

    auto tk_zzzzzz_yyyyy = pbuffer.data(idx_kin_ih + 582);

    auto tk_zzzzzz_yyyyz = pbuffer.data(idx_kin_ih + 583);

    auto tk_zzzzzz_yyyzz = pbuffer.data(idx_kin_ih + 584);

    auto tk_zzzzzz_yyzzz = pbuffer.data(idx_kin_ih + 585);

    auto tk_zzzzzz_yzzzz = pbuffer.data(idx_kin_ih + 586);

    auto tk_zzzzzz_zzzzz = pbuffer.data(idx_kin_ih + 587);

    #pragma omp simd aligned(pa_z, tk_zzzz_xxxxx, tk_zzzz_xxxxy, tk_zzzz_xxxxz, tk_zzzz_xxxyy, tk_zzzz_xxxyz, tk_zzzz_xxxzz, tk_zzzz_xxyyy, tk_zzzz_xxyyz, tk_zzzz_xxyzz, tk_zzzz_xxzzz, tk_zzzz_xyyyy, tk_zzzz_xyyyz, tk_zzzz_xyyzz, tk_zzzz_xyzzz, tk_zzzz_xzzzz, tk_zzzz_yyyyy, tk_zzzz_yyyyz, tk_zzzz_yyyzz, tk_zzzz_yyzzz, tk_zzzz_yzzzz, tk_zzzz_zzzzz, tk_zzzzz_xxxx, tk_zzzzz_xxxxx, tk_zzzzz_xxxxy, tk_zzzzz_xxxxz, tk_zzzzz_xxxy, tk_zzzzz_xxxyy, tk_zzzzz_xxxyz, tk_zzzzz_xxxz, tk_zzzzz_xxxzz, tk_zzzzz_xxyy, tk_zzzzz_xxyyy, tk_zzzzz_xxyyz, tk_zzzzz_xxyz, tk_zzzzz_xxyzz, tk_zzzzz_xxzz, tk_zzzzz_xxzzz, tk_zzzzz_xyyy, tk_zzzzz_xyyyy, tk_zzzzz_xyyyz, tk_zzzzz_xyyz, tk_zzzzz_xyyzz, tk_zzzzz_xyzz, tk_zzzzz_xyzzz, tk_zzzzz_xzzz, tk_zzzzz_xzzzz, tk_zzzzz_yyyy, tk_zzzzz_yyyyy, tk_zzzzz_yyyyz, tk_zzzzz_yyyz, tk_zzzzz_yyyzz, tk_zzzzz_yyzz, tk_zzzzz_yyzzz, tk_zzzzz_yzzz, tk_zzzzz_yzzzz, tk_zzzzz_zzzz, tk_zzzzz_zzzzz, tk_zzzzzz_xxxxx, tk_zzzzzz_xxxxy, tk_zzzzzz_xxxxz, tk_zzzzzz_xxxyy, tk_zzzzzz_xxxyz, tk_zzzzzz_xxxzz, tk_zzzzzz_xxyyy, tk_zzzzzz_xxyyz, tk_zzzzzz_xxyzz, tk_zzzzzz_xxzzz, tk_zzzzzz_xyyyy, tk_zzzzzz_xyyyz, tk_zzzzzz_xyyzz, tk_zzzzzz_xyzzz, tk_zzzzzz_xzzzz, tk_zzzzzz_yyyyy, tk_zzzzzz_yyyyz, tk_zzzzzz_yyyzz, tk_zzzzzz_yyzzz, tk_zzzzzz_yzzzz, tk_zzzzzz_zzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, ts_zzzzzz_xxxxx, ts_zzzzzz_xxxxy, ts_zzzzzz_xxxxz, ts_zzzzzz_xxxyy, ts_zzzzzz_xxxyz, ts_zzzzzz_xxxzz, ts_zzzzzz_xxyyy, ts_zzzzzz_xxyyz, ts_zzzzzz_xxyzz, ts_zzzzzz_xxzzz, ts_zzzzzz_xyyyy, ts_zzzzzz_xyyyz, ts_zzzzzz_xyyzz, ts_zzzzzz_xyzzz, ts_zzzzzz_xzzzz, ts_zzzzzz_yyyyy, ts_zzzzzz_yyyyz, ts_zzzzzz_yyyzz, ts_zzzzzz_yyzzz, ts_zzzzzz_yzzzz, ts_zzzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzzz_xxxxx[i] = -10.0 * ts_zzzz_xxxxx[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxx[i] * fe_0 + tk_zzzzz_xxxxx[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxxx[i] * fz_0;

        tk_zzzzzz_xxxxy[i] = -10.0 * ts_zzzz_xxxxy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxy[i] * fe_0 + tk_zzzzz_xxxxy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxxy[i] * fz_0;

        tk_zzzzzz_xxxxz[i] = -10.0 * ts_zzzz_xxxxz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxxz[i] * fe_0 + tk_zzzzz_xxxx[i] * fe_0 + tk_zzzzz_xxxxz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxxz[i] * fz_0;

        tk_zzzzzz_xxxyy[i] = -10.0 * ts_zzzz_xxxyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxyy[i] * fe_0 + tk_zzzzz_xxxyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxyy[i] * fz_0;

        tk_zzzzzz_xxxyz[i] = -10.0 * ts_zzzz_xxxyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxyz[i] * fe_0 + tk_zzzzz_xxxy[i] * fe_0 + tk_zzzzz_xxxyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxyz[i] * fz_0;

        tk_zzzzzz_xxxzz[i] = -10.0 * ts_zzzz_xxxzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxxzz[i] * fe_0 + 2.0 * tk_zzzzz_xxxz[i] * fe_0 + tk_zzzzz_xxxzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxxzz[i] * fz_0;

        tk_zzzzzz_xxyyy[i] = -10.0 * ts_zzzz_xxyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyyy[i] * fe_0 + tk_zzzzz_xxyyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxyyy[i] * fz_0;

        tk_zzzzzz_xxyyz[i] = -10.0 * ts_zzzz_xxyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyyz[i] * fe_0 + tk_zzzzz_xxyy[i] * fe_0 + tk_zzzzz_xxyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxyyz[i] * fz_0;

        tk_zzzzzz_xxyzz[i] = -10.0 * ts_zzzz_xxyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxyzz[i] * fe_0 + 2.0 * tk_zzzzz_xxyz[i] * fe_0 + tk_zzzzz_xxyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxyzz[i] * fz_0;

        tk_zzzzzz_xxzzz[i] = -10.0 * ts_zzzz_xxzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xxzzz[i] * fe_0 + 3.0 * tk_zzzzz_xxzz[i] * fe_0 + tk_zzzzz_xxzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xxzzz[i] * fz_0;

        tk_zzzzzz_xyyyy[i] = -10.0 * ts_zzzz_xyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyyy[i] * fe_0 + tk_zzzzz_xyyyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyyyy[i] * fz_0;

        tk_zzzzzz_xyyyz[i] = -10.0 * ts_zzzz_xyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyyz[i] * fe_0 + tk_zzzzz_xyyy[i] * fe_0 + tk_zzzzz_xyyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyyyz[i] * fz_0;

        tk_zzzzzz_xyyzz[i] = -10.0 * ts_zzzz_xyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyyzz[i] * fe_0 + 2.0 * tk_zzzzz_xyyz[i] * fe_0 + tk_zzzzz_xyyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyyzz[i] * fz_0;

        tk_zzzzzz_xyzzz[i] = -10.0 * ts_zzzz_xyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xyzzz[i] * fe_0 + 3.0 * tk_zzzzz_xyzz[i] * fe_0 + tk_zzzzz_xyzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xyzzz[i] * fz_0;

        tk_zzzzzz_xzzzz[i] = -10.0 * ts_zzzz_xzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xzzzz[i] * fe_0 + 4.0 * tk_zzzzz_xzzz[i] * fe_0 + tk_zzzzz_xzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_xzzzz[i] * fz_0;

        tk_zzzzzz_yyyyy[i] = -10.0 * ts_zzzz_yyyyy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyyy[i] * fe_0 + tk_zzzzz_yyyyy[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyyyy[i] * fz_0;

        tk_zzzzzz_yyyyz[i] = -10.0 * ts_zzzz_yyyyz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyyz[i] * fe_0 + tk_zzzzz_yyyy[i] * fe_0 + tk_zzzzz_yyyyz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyyyz[i] * fz_0;

        tk_zzzzzz_yyyzz[i] = -10.0 * ts_zzzz_yyyzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyyzz[i] * fe_0 + 2.0 * tk_zzzzz_yyyz[i] * fe_0 + tk_zzzzz_yyyzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyyzz[i] * fz_0;

        tk_zzzzzz_yyzzz[i] = -10.0 * ts_zzzz_yyzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yyzzz[i] * fe_0 + 3.0 * tk_zzzzz_yyzz[i] * fe_0 + tk_zzzzz_yyzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yyzzz[i] * fz_0;

        tk_zzzzzz_yzzzz[i] = -10.0 * ts_zzzz_yzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yzzzz[i] * fe_0 + 4.0 * tk_zzzzz_yzzz[i] * fe_0 + tk_zzzzz_yzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_yzzzz[i] * fz_0;

        tk_zzzzzz_zzzzz[i] = -10.0 * ts_zzzz_zzzzz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_zzzzz[i] * fe_0 + 5.0 * tk_zzzzz_zzzz[i] * fe_0 + tk_zzzzz_zzzzz[i] * pa_z[i] + 2.0 * ts_zzzzzz_zzzzz[i] * fz_0;
    }

}

} // kinrec namespace

