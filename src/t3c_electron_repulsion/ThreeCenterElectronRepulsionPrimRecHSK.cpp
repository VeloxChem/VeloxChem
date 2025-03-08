#include "ThreeCenterElectronRepulsionPrimRecHSK.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hsk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsk,
                                 size_t idx_eri_0_fsk,
                                 size_t idx_eri_1_fsk,
                                 size_t idx_eri_1_gsi,
                                 size_t idx_eri_1_gsk,
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

    /// Set up components of auxilary buffer : FSK

    auto g_xxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk);

    auto g_xxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 1);

    auto g_xxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 2);

    auto g_xxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 3);

    auto g_xxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 4);

    auto g_xxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 5);

    auto g_xxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 6);

    auto g_xxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 7);

    auto g_xxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 8);

    auto g_xxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 9);

    auto g_xxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 10);

    auto g_xxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 11);

    auto g_xxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 12);

    auto g_xxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 13);

    auto g_xxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 14);

    auto g_xxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 15);

    auto g_xxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 16);

    auto g_xxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 17);

    auto g_xxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 18);

    auto g_xxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 19);

    auto g_xxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 20);

    auto g_xxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 21);

    auto g_xxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 22);

    auto g_xxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 23);

    auto g_xxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 24);

    auto g_xxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 25);

    auto g_xxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 26);

    auto g_xxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 27);

    auto g_xxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 28);

    auto g_xxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 29);

    auto g_xxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 30);

    auto g_xxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 31);

    auto g_xxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 32);

    auto g_xxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 33);

    auto g_xxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 34);

    auto g_xxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 35);

    auto g_xxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 36);

    auto g_xxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 38);

    auto g_xxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 41);

    auto g_xxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 45);

    auto g_xxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 50);

    auto g_xxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 56);

    auto g_xxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 63);

    auto g_xxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 72);

    auto g_xxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 73);

    auto g_xxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 75);

    auto g_xxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 78);

    auto g_xxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 82);

    auto g_xxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 87);

    auto g_xxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 93);

    auto g_xyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 109);

    auto g_xyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 111);

    auto g_xyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 112);

    auto g_xyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 114);

    auto g_xyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 115);

    auto g_xyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 116);

    auto g_xyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 118);

    auto g_xyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 119);

    auto g_xyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 120);

    auto g_xyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 121);

    auto g_xyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 123);

    auto g_xyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 124);

    auto g_xyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 125);

    auto g_xyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 126);

    auto g_xyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 127);

    auto g_xyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 129);

    auto g_xyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 130);

    auto g_xyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 131);

    auto g_xyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 132);

    auto g_xyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 133);

    auto g_xyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 134);

    auto g_xyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 136);

    auto g_xyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 137);

    auto g_xyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 138);

    auto g_xyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 139);

    auto g_xyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 140);

    auto g_xyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 141);

    auto g_xyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 142);

    auto g_xyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 143);

    auto g_xzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 182);

    auto g_xzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 184);

    auto g_xzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 185);

    auto g_xzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 187);

    auto g_xzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 188);

    auto g_xzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 189);

    auto g_xzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 191);

    auto g_xzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 192);

    auto g_xzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 193);

    auto g_xzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 194);

    auto g_xzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 196);

    auto g_xzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 197);

    auto g_xzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 198);

    auto g_xzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 199);

    auto g_xzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 200);

    auto g_xzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 202);

    auto g_xzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 203);

    auto g_xzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 204);

    auto g_xzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 205);

    auto g_xzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 206);

    auto g_xzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 207);

    auto g_xzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 208);

    auto g_xzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 209);

    auto g_xzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 210);

    auto g_xzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 211);

    auto g_xzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 212);

    auto g_xzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 213);

    auto g_xzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 214);

    auto g_xzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 215);

    auto g_yyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 216);

    auto g_yyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 217);

    auto g_yyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 218);

    auto g_yyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 219);

    auto g_yyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 220);

    auto g_yyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 221);

    auto g_yyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 222);

    auto g_yyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 223);

    auto g_yyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 224);

    auto g_yyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 225);

    auto g_yyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 226);

    auto g_yyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 227);

    auto g_yyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 228);

    auto g_yyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 229);

    auto g_yyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 230);

    auto g_yyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 231);

    auto g_yyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 232);

    auto g_yyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 233);

    auto g_yyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 234);

    auto g_yyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 235);

    auto g_yyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 236);

    auto g_yyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 237);

    auto g_yyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 238);

    auto g_yyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 239);

    auto g_yyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 240);

    auto g_yyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 241);

    auto g_yyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 242);

    auto g_yyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 243);

    auto g_yyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 244);

    auto g_yyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 245);

    auto g_yyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 246);

    auto g_yyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 247);

    auto g_yyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 248);

    auto g_yyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 249);

    auto g_yyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 250);

    auto g_yyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 251);

    auto g_yyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 253);

    auto g_yyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 255);

    auto g_yyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 258);

    auto g_yyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 262);

    auto g_yyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 267);

    auto g_yyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 273);

    auto g_yyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 280);

    auto g_yzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 288);

    auto g_yzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 290);

    auto g_yzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 292);

    auto g_yzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 293);

    auto g_yzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 295);

    auto g_yzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 296);

    auto g_yzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 297);

    auto g_yzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 299);

    auto g_yzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 300);

    auto g_yzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 301);

    auto g_yzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 302);

    auto g_yzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 304);

    auto g_yzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 305);

    auto g_yzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 306);

    auto g_yzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 307);

    auto g_yzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 308);

    auto g_yzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 310);

    auto g_yzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 311);

    auto g_yzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 312);

    auto g_yzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 313);

    auto g_yzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 314);

    auto g_yzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 315);

    auto g_yzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 317);

    auto g_yzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 318);

    auto g_yzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 319);

    auto g_yzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 320);

    auto g_yzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 321);

    auto g_yzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 322);

    auto g_yzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 323);

    auto g_zzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_fsk + 324);

    auto g_zzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_fsk + 325);

    auto g_zzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_fsk + 326);

    auto g_zzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_fsk + 327);

    auto g_zzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_fsk + 328);

    auto g_zzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_fsk + 329);

    auto g_zzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_fsk + 330);

    auto g_zzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_fsk + 331);

    auto g_zzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_fsk + 332);

    auto g_zzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_fsk + 333);

    auto g_zzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_fsk + 334);

    auto g_zzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_fsk + 335);

    auto g_zzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_fsk + 336);

    auto g_zzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_fsk + 337);

    auto g_zzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_fsk + 338);

    auto g_zzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 339);

    auto g_zzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 340);

    auto g_zzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 341);

    auto g_zzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 342);

    auto g_zzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 343);

    auto g_zzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 344);

    auto g_zzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 345);

    auto g_zzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 346);

    auto g_zzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 347);

    auto g_zzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 348);

    auto g_zzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 349);

    auto g_zzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 350);

    auto g_zzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 351);

    auto g_zzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_fsk + 352);

    auto g_zzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_fsk + 353);

    auto g_zzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_fsk + 354);

    auto g_zzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_fsk + 355);

    auto g_zzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_fsk + 356);

    auto g_zzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 357);

    auto g_zzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 358);

    auto g_zzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_fsk + 359);

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

    auto g_xxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 36);

    auto g_xxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 38);

    auto g_xxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 41);

    auto g_xxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 45);

    auto g_xxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 50);

    auto g_xxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 56);

    auto g_xxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 63);

    auto g_xxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 72);

    auto g_xxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 73);

    auto g_xxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 75);

    auto g_xxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 78);

    auto g_xxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 82);

    auto g_xxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 87);

    auto g_xxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 93);

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

    auto g_xyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 143);

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

    auto g_xzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 208);

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

    auto g_yyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 253);

    auto g_yyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 255);

    auto g_yyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 258);

    auto g_yyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 262);

    auto g_yyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 267);

    auto g_yyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 273);

    auto g_yyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 280);

    auto g_yzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 288);

    auto g_yzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 290);

    auto g_yzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 292);

    auto g_yzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 293);

    auto g_yzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 295);

    auto g_yzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 296);

    auto g_yzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 297);

    auto g_yzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 299);

    auto g_yzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 300);

    auto g_yzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 301);

    auto g_yzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 302);

    auto g_yzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 304);

    auto g_yzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 305);

    auto g_yzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 306);

    auto g_yzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 307);

    auto g_yzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 308);

    auto g_yzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 310);

    auto g_yzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 311);

    auto g_yzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 312);

    auto g_yzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 313);

    auto g_yzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 314);

    auto g_yzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 315);

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

    /// Set up components of auxilary buffer : GSI

    auto g_xxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi);

    auto g_xxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 1);

    auto g_xxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 2);

    auto g_xxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 3);

    auto g_xxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 4);

    auto g_xxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 5);

    auto g_xxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 6);

    auto g_xxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 7);

    auto g_xxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 8);

    auto g_xxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 9);

    auto g_xxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 10);

    auto g_xxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 11);

    auto g_xxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 12);

    auto g_xxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 13);

    auto g_xxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 14);

    auto g_xxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 15);

    auto g_xxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 16);

    auto g_xxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 17);

    auto g_xxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 18);

    auto g_xxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 19);

    auto g_xxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 20);

    auto g_xxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 21);

    auto g_xxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 22);

    auto g_xxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 23);

    auto g_xxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 24);

    auto g_xxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 25);

    auto g_xxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 26);

    auto g_xxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 27);

    auto g_xxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 58);

    auto g_xxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 60);

    auto g_xxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 61);

    auto g_xxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 63);

    auto g_xxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 64);

    auto g_xxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 65);

    auto g_xxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 67);

    auto g_xxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 68);

    auto g_xxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 69);

    auto g_xxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 70);

    auto g_xxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 72);

    auto g_xxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 73);

    auto g_xxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 74);

    auto g_xxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 75);

    auto g_xxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 76);

    auto g_xxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 78);

    auto g_xxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 79);

    auto g_xxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 80);

    auto g_xxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 81);

    auto g_xxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 82);

    auto g_xxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 83);

    auto g_xxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 84);

    auto g_xxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 85);

    auto g_xxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 86);

    auto g_xxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 87);

    auto g_xxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 88);

    auto g_xxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 89);

    auto g_xxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 90);

    auto g_xxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 91);

    auto g_xxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 92);

    auto g_xxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 93);

    auto g_xxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 94);

    auto g_xxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 95);

    auto g_xxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 96);

    auto g_xxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 97);

    auto g_xxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 98);

    auto g_xxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 99);

    auto g_xxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 100);

    auto g_xxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 101);

    auto g_xxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 102);

    auto g_xxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 103);

    auto g_xxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 104);

    auto g_xxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 105);

    auto g_xxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 106);

    auto g_xxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 107);

    auto g_xxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 108);

    auto g_xxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 109);

    auto g_xxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 110);

    auto g_xxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 111);

    auto g_xxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 140);

    auto g_xxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 141);

    auto g_xxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 142);

    auto g_xxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 143);

    auto g_xxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 144);

    auto g_xxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 145);

    auto g_xxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 146);

    auto g_xxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 147);

    auto g_xxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 148);

    auto g_xxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 149);

    auto g_xxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 150);

    auto g_xxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 151);

    auto g_xxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 152);

    auto g_xxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 153);

    auto g_xxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 154);

    auto g_xxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 155);

    auto g_xxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 156);

    auto g_xxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 157);

    auto g_xxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 158);

    auto g_xxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 159);

    auto g_xxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 160);

    auto g_xxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 161);

    auto g_xxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 162);

    auto g_xxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 163);

    auto g_xxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 164);

    auto g_xxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 165);

    auto g_xxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 166);

    auto g_xxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 167);

    auto g_xyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 169);

    auto g_xyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 171);

    auto g_xyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 172);

    auto g_xyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 174);

    auto g_xyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 175);

    auto g_xyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 176);

    auto g_xyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 178);

    auto g_xyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 179);

    auto g_xyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 180);

    auto g_xyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 181);

    auto g_xyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 183);

    auto g_xyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 184);

    auto g_xyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 185);

    auto g_xyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 186);

    auto g_xyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 187);

    auto g_xyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 189);

    auto g_xyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 190);

    auto g_xyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 191);

    auto g_xyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 192);

    auto g_xyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 193);

    auto g_xyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 194);

    auto g_xzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 254);

    auto g_xzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 256);

    auto g_xzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 257);

    auto g_xzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 259);

    auto g_xzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 260);

    auto g_xzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 261);

    auto g_xzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 263);

    auto g_xzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 264);

    auto g_xzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 265);

    auto g_xzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 266);

    auto g_xzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 268);

    auto g_xzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 269);

    auto g_xzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 270);

    auto g_xzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 271);

    auto g_xzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 272);

    auto g_xzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 274);

    auto g_xzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 275);

    auto g_xzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 276);

    auto g_xzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 277);

    auto g_xzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 278);

    auto g_xzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 279);

    auto g_yyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 280);

    auto g_yyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 281);

    auto g_yyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 282);

    auto g_yyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 283);

    auto g_yyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 284);

    auto g_yyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 285);

    auto g_yyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 286);

    auto g_yyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 287);

    auto g_yyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 288);

    auto g_yyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 289);

    auto g_yyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 290);

    auto g_yyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 291);

    auto g_yyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 292);

    auto g_yyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 293);

    auto g_yyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 294);

    auto g_yyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 295);

    auto g_yyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 296);

    auto g_yyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 297);

    auto g_yyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 298);

    auto g_yyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 299);

    auto g_yyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 300);

    auto g_yyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 301);

    auto g_yyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 302);

    auto g_yyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 303);

    auto g_yyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 304);

    auto g_yyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 305);

    auto g_yyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 306);

    auto g_yyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 307);

    auto g_yyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 310);

    auto g_yyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 312);

    auto g_yyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 313);

    auto g_yyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 315);

    auto g_yyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 316);

    auto g_yyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 317);

    auto g_yyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 319);

    auto g_yyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 320);

    auto g_yyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 321);

    auto g_yyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 322);

    auto g_yyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 324);

    auto g_yyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 325);

    auto g_yyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 326);

    auto g_yyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 327);

    auto g_yyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 328);

    auto g_yyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 330);

    auto g_yyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 331);

    auto g_yyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 332);

    auto g_yyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 333);

    auto g_yyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 334);

    auto g_yyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 335);

    auto g_yyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 336);

    auto g_yyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 337);

    auto g_yyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 338);

    auto g_yyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 339);

    auto g_yyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 340);

    auto g_yyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 341);

    auto g_yyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 342);

    auto g_yyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 343);

    auto g_yyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 344);

    auto g_yyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 345);

    auto g_yyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 346);

    auto g_yyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 347);

    auto g_yyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 348);

    auto g_yyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 349);

    auto g_yyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 350);

    auto g_yyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 351);

    auto g_yyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 352);

    auto g_yyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 353);

    auto g_yyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 354);

    auto g_yyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 355);

    auto g_yyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 356);

    auto g_yyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 357);

    auto g_yyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 358);

    auto g_yyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 359);

    auto g_yyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 360);

    auto g_yyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 361);

    auto g_yyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 362);

    auto g_yyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 363);

    auto g_yzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 365);

    auto g_yzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 366);

    auto g_yzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 367);

    auto g_yzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 368);

    auto g_yzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 369);

    auto g_yzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 370);

    auto g_yzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 371);

    auto g_yzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 372);

    auto g_yzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 373);

    auto g_yzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 374);

    auto g_yzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 375);

    auto g_yzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 376);

    auto g_yzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 377);

    auto g_yzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 378);

    auto g_yzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 379);

    auto g_yzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 380);

    auto g_yzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 381);

    auto g_yzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 382);

    auto g_yzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 383);

    auto g_yzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 384);

    auto g_yzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 385);

    auto g_yzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 386);

    auto g_yzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 387);

    auto g_yzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 388);

    auto g_yzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 389);

    auto g_yzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 390);

    auto g_yzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 391);

    auto g_zzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 392);

    auto g_zzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 393);

    auto g_zzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 394);

    auto g_zzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 395);

    auto g_zzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 396);

    auto g_zzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 397);

    auto g_zzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 398);

    auto g_zzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 399);

    auto g_zzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 400);

    auto g_zzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 401);

    auto g_zzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 402);

    auto g_zzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 403);

    auto g_zzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 404);

    auto g_zzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 405);

    auto g_zzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 406);

    auto g_zzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 407);

    auto g_zzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 408);

    auto g_zzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 409);

    auto g_zzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 410);

    auto g_zzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 411);

    auto g_zzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 412);

    auto g_zzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 413);

    auto g_zzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 414);

    auto g_zzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 415);

    auto g_zzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 416);

    auto g_zzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 417);

    auto g_zzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 418);

    auto g_zzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 419);

    /// Set up components of auxilary buffer : GSK

    auto g_xxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk);

    auto g_xxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 1);

    auto g_xxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 2);

    auto g_xxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 3);

    auto g_xxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 4);

    auto g_xxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 5);

    auto g_xxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 6);

    auto g_xxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 7);

    auto g_xxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 8);

    auto g_xxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 9);

    auto g_xxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 10);

    auto g_xxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 11);

    auto g_xxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 12);

    auto g_xxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 13);

    auto g_xxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 14);

    auto g_xxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 15);

    auto g_xxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 16);

    auto g_xxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 17);

    auto g_xxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 18);

    auto g_xxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 19);

    auto g_xxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 20);

    auto g_xxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 21);

    auto g_xxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 22);

    auto g_xxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 23);

    auto g_xxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 24);

    auto g_xxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 25);

    auto g_xxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 26);

    auto g_xxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 27);

    auto g_xxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 28);

    auto g_xxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 29);

    auto g_xxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 30);

    auto g_xxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 31);

    auto g_xxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 32);

    auto g_xxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 33);

    auto g_xxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 34);

    auto g_xxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 35);

    auto g_xxxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 36);

    auto g_xxxy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 37);

    auto g_xxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 38);

    auto g_xxxy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 39);

    auto g_xxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 41);

    auto g_xxxy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 42);

    auto g_xxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 45);

    auto g_xxxy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 46);

    auto g_xxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 50);

    auto g_xxxy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 51);

    auto g_xxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 56);

    auto g_xxxy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 57);

    auto g_xxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 63);

    auto g_xxxy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 64);

    auto g_xxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 72);

    auto g_xxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 73);

    auto g_xxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 74);

    auto g_xxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 75);

    auto g_xxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 76);

    auto g_xxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 77);

    auto g_xxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 78);

    auto g_xxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 79);

    auto g_xxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 80);

    auto g_xxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 81);

    auto g_xxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 82);

    auto g_xxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 83);

    auto g_xxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 84);

    auto g_xxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 85);

    auto g_xxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 86);

    auto g_xxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 87);

    auto g_xxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 88);

    auto g_xxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 89);

    auto g_xxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 90);

    auto g_xxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 91);

    auto g_xxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 92);

    auto g_xxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 93);

    auto g_xxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 94);

    auto g_xxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 95);

    auto g_xxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 96);

    auto g_xxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 97);

    auto g_xxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 98);

    auto g_xxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 99);

    auto g_xxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 101);

    auto g_xxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 102);

    auto g_xxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 103);

    auto g_xxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 104);

    auto g_xxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 105);

    auto g_xxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 106);

    auto g_xxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 107);

    auto g_xxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 108);

    auto g_xxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 109);

    auto g_xxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 110);

    auto g_xxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 111);

    auto g_xxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 112);

    auto g_xxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 113);

    auto g_xxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 114);

    auto g_xxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 115);

    auto g_xxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 116);

    auto g_xxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 117);

    auto g_xxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 118);

    auto g_xxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 119);

    auto g_xxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 120);

    auto g_xxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 121);

    auto g_xxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 122);

    auto g_xxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 123);

    auto g_xxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 124);

    auto g_xxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 125);

    auto g_xxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 126);

    auto g_xxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 127);

    auto g_xxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 128);

    auto g_xxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 129);

    auto g_xxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 130);

    auto g_xxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 131);

    auto g_xxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 132);

    auto g_xxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 133);

    auto g_xxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 134);

    auto g_xxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 135);

    auto g_xxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 136);

    auto g_xxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 137);

    auto g_xxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 138);

    auto g_xxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 139);

    auto g_xxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 140);

    auto g_xxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 141);

    auto g_xxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 142);

    auto g_xxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 143);

    auto g_xxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 180);

    auto g_xxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 181);

    auto g_xxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 182);

    auto g_xxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 183);

    auto g_xxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 184);

    auto g_xxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 185);

    auto g_xxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 186);

    auto g_xxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 187);

    auto g_xxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 188);

    auto g_xxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 189);

    auto g_xxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 190);

    auto g_xxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 191);

    auto g_xxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 192);

    auto g_xxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 193);

    auto g_xxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 194);

    auto g_xxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 195);

    auto g_xxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 196);

    auto g_xxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 197);

    auto g_xxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 198);

    auto g_xxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 199);

    auto g_xxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 200);

    auto g_xxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 201);

    auto g_xxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 202);

    auto g_xxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 203);

    auto g_xxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 204);

    auto g_xxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 205);

    auto g_xxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 206);

    auto g_xxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 207);

    auto g_xxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 208);

    auto g_xxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 209);

    auto g_xxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 210);

    auto g_xxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 211);

    auto g_xxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 212);

    auto g_xxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 213);

    auto g_xxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 214);

    auto g_xxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 215);

    auto g_xyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 216);

    auto g_xyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 217);

    auto g_xyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 219);

    auto g_xyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 220);

    auto g_xyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 222);

    auto g_xyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 223);

    auto g_xyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 224);

    auto g_xyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 226);

    auto g_xyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 227);

    auto g_xyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 228);

    auto g_xyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 229);

    auto g_xyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 231);

    auto g_xyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 232);

    auto g_xyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 233);

    auto g_xyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 234);

    auto g_xyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 235);

    auto g_xyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 237);

    auto g_xyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 238);

    auto g_xyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 239);

    auto g_xyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 240);

    auto g_xyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 241);

    auto g_xyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 242);

    auto g_xyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 244);

    auto g_xyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 245);

    auto g_xyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 246);

    auto g_xyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 247);

    auto g_xyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 248);

    auto g_xyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 249);

    auto g_xyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 250);

    auto g_xyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 251);

    auto g_xzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 324);

    auto g_xzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 326);

    auto g_xzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 328);

    auto g_xzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 329);

    auto g_xzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 331);

    auto g_xzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 332);

    auto g_xzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 333);

    auto g_xzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 335);

    auto g_xzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 336);

    auto g_xzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 337);

    auto g_xzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 338);

    auto g_xzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 340);

    auto g_xzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 341);

    auto g_xzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 342);

    auto g_xzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 343);

    auto g_xzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 344);

    auto g_xzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 346);

    auto g_xzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 347);

    auto g_xzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 348);

    auto g_xzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 349);

    auto g_xzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 350);

    auto g_xzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 351);

    auto g_xzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 352);

    auto g_xzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 353);

    auto g_xzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 354);

    auto g_xzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 355);

    auto g_xzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 356);

    auto g_xzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 357);

    auto g_xzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 358);

    auto g_xzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 359);

    auto g_yyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 360);

    auto g_yyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 361);

    auto g_yyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 362);

    auto g_yyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 363);

    auto g_yyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 364);

    auto g_yyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 365);

    auto g_yyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 366);

    auto g_yyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 367);

    auto g_yyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 368);

    auto g_yyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 369);

    auto g_yyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 370);

    auto g_yyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 371);

    auto g_yyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 372);

    auto g_yyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 373);

    auto g_yyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 374);

    auto g_yyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 375);

    auto g_yyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 376);

    auto g_yyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 377);

    auto g_yyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 378);

    auto g_yyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 379);

    auto g_yyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 380);

    auto g_yyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 381);

    auto g_yyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 382);

    auto g_yyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 383);

    auto g_yyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 384);

    auto g_yyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 385);

    auto g_yyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 386);

    auto g_yyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 387);

    auto g_yyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 388);

    auto g_yyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 389);

    auto g_yyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 390);

    auto g_yyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 391);

    auto g_yyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 392);

    auto g_yyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 393);

    auto g_yyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 394);

    auto g_yyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 395);

    auto g_yyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 397);

    auto g_yyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 398);

    auto g_yyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 399);

    auto g_yyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 400);

    auto g_yyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 401);

    auto g_yyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 402);

    auto g_yyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 403);

    auto g_yyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 404);

    auto g_yyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 405);

    auto g_yyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 406);

    auto g_yyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 407);

    auto g_yyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 408);

    auto g_yyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 409);

    auto g_yyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 410);

    auto g_yyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 411);

    auto g_yyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 412);

    auto g_yyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 413);

    auto g_yyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 414);

    auto g_yyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 415);

    auto g_yyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 416);

    auto g_yyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 417);

    auto g_yyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 418);

    auto g_yyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 419);

    auto g_yyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 420);

    auto g_yyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 421);

    auto g_yyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 422);

    auto g_yyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 423);

    auto g_yyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 424);

    auto g_yyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 425);

    auto g_yyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 426);

    auto g_yyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 427);

    auto g_yyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 428);

    auto g_yyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 429);

    auto g_yyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 430);

    auto g_yyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 431);

    auto g_yyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 432);

    auto g_yyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 433);

    auto g_yyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 434);

    auto g_yyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 435);

    auto g_yyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 436);

    auto g_yyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 437);

    auto g_yyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 438);

    auto g_yyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 439);

    auto g_yyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 440);

    auto g_yyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 441);

    auto g_yyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 442);

    auto g_yyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 443);

    auto g_yyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 444);

    auto g_yyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 445);

    auto g_yyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 446);

    auto g_yyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 447);

    auto g_yyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 448);

    auto g_yyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 449);

    auto g_yyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 450);

    auto g_yyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 451);

    auto g_yyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 452);

    auto g_yyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 453);

    auto g_yyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 454);

    auto g_yyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 455);

    auto g_yyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 456);

    auto g_yyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 457);

    auto g_yyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 458);

    auto g_yyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 459);

    auto g_yyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 460);

    auto g_yyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 461);

    auto g_yyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 462);

    auto g_yyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 463);

    auto g_yyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 464);

    auto g_yyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 465);

    auto g_yyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 466);

    auto g_yyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 467);

    auto g_yzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 468);

    auto g_yzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 469);

    auto g_yzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 470);

    auto g_yzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 471);

    auto g_yzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 472);

    auto g_yzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 473);

    auto g_yzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 474);

    auto g_yzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 475);

    auto g_yzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 476);

    auto g_yzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 477);

    auto g_yzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 478);

    auto g_yzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 479);

    auto g_yzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 480);

    auto g_yzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 481);

    auto g_yzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 482);

    auto g_yzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 483);

    auto g_yzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 484);

    auto g_yzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 485);

    auto g_yzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 486);

    auto g_yzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 487);

    auto g_yzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 488);

    auto g_yzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 489);

    auto g_yzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 490);

    auto g_yzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 491);

    auto g_yzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 492);

    auto g_yzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 493);

    auto g_yzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 494);

    auto g_yzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 495);

    auto g_yzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 496);

    auto g_yzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 497);

    auto g_yzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 498);

    auto g_yzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 499);

    auto g_yzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 500);

    auto g_yzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 501);

    auto g_yzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 502);

    auto g_yzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 503);

    auto g_zzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 504);

    auto g_zzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 505);

    auto g_zzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 506);

    auto g_zzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 507);

    auto g_zzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 508);

    auto g_zzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 509);

    auto g_zzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 510);

    auto g_zzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 511);

    auto g_zzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 512);

    auto g_zzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 513);

    auto g_zzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 514);

    auto g_zzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 515);

    auto g_zzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 516);

    auto g_zzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 517);

    auto g_zzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 518);

    auto g_zzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 519);

    auto g_zzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 520);

    auto g_zzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 521);

    auto g_zzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 522);

    auto g_zzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 523);

    auto g_zzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 524);

    auto g_zzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 525);

    auto g_zzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 526);

    auto g_zzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 527);

    auto g_zzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 528);

    auto g_zzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 529);

    auto g_zzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 530);

    auto g_zzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 531);

    auto g_zzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 532);

    auto g_zzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 533);

    auto g_zzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 534);

    auto g_zzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 535);

    auto g_zzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 536);

    auto g_zzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 537);

    auto g_zzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 538);

    auto g_zzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 539);

    /// Set up 0-36 components of targeted buffer : HSK

    auto g_xxxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk);

    auto g_xxxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 1);

    auto g_xxxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 2);

    auto g_xxxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 3);

    auto g_xxxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 4);

    auto g_xxxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 5);

    auto g_xxxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 6);

    auto g_xxxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 7);

    auto g_xxxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 8);

    auto g_xxxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 9);

    auto g_xxxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 10);

    auto g_xxxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 11);

    auto g_xxxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 12);

    auto g_xxxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 13);

    auto g_xxxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 14);

    auto g_xxxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 15);

    auto g_xxxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 16);

    auto g_xxxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 17);

    auto g_xxxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 18);

    auto g_xxxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 19);

    auto g_xxxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 20);

    auto g_xxxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 21);

    auto g_xxxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 22);

    auto g_xxxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 23);

    auto g_xxxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 24);

    auto g_xxxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 25);

    auto g_xxxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 26);

    auto g_xxxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 27);

    auto g_xxxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 28);

    auto g_xxxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 29);

    auto g_xxxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 30);

    auto g_xxxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 31);

    auto g_xxxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 32);

    auto g_xxxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 33);

    auto g_xxxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 34);

    auto g_xxxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 35);

    #pragma omp simd aligned(g_xxx_0_xxxxxxx_0, g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxy_0, g_xxx_0_xxxxxxy_1, g_xxx_0_xxxxxxz_0, g_xxx_0_xxxxxxz_1, g_xxx_0_xxxxxyy_0, g_xxx_0_xxxxxyy_1, g_xxx_0_xxxxxyz_0, g_xxx_0_xxxxxyz_1, g_xxx_0_xxxxxzz_0, g_xxx_0_xxxxxzz_1, g_xxx_0_xxxxyyy_0, g_xxx_0_xxxxyyy_1, g_xxx_0_xxxxyyz_0, g_xxx_0_xxxxyyz_1, g_xxx_0_xxxxyzz_0, g_xxx_0_xxxxyzz_1, g_xxx_0_xxxxzzz_0, g_xxx_0_xxxxzzz_1, g_xxx_0_xxxyyyy_0, g_xxx_0_xxxyyyy_1, g_xxx_0_xxxyyyz_0, g_xxx_0_xxxyyyz_1, g_xxx_0_xxxyyzz_0, g_xxx_0_xxxyyzz_1, g_xxx_0_xxxyzzz_0, g_xxx_0_xxxyzzz_1, g_xxx_0_xxxzzzz_0, g_xxx_0_xxxzzzz_1, g_xxx_0_xxyyyyy_0, g_xxx_0_xxyyyyy_1, g_xxx_0_xxyyyyz_0, g_xxx_0_xxyyyyz_1, g_xxx_0_xxyyyzz_0, g_xxx_0_xxyyyzz_1, g_xxx_0_xxyyzzz_0, g_xxx_0_xxyyzzz_1, g_xxx_0_xxyzzzz_0, g_xxx_0_xxyzzzz_1, g_xxx_0_xxzzzzz_0, g_xxx_0_xxzzzzz_1, g_xxx_0_xyyyyyy_0, g_xxx_0_xyyyyyy_1, g_xxx_0_xyyyyyz_0, g_xxx_0_xyyyyyz_1, g_xxx_0_xyyyyzz_0, g_xxx_0_xyyyyzz_1, g_xxx_0_xyyyzzz_0, g_xxx_0_xyyyzzz_1, g_xxx_0_xyyzzzz_0, g_xxx_0_xyyzzzz_1, g_xxx_0_xyzzzzz_0, g_xxx_0_xyzzzzz_1, g_xxx_0_xzzzzzz_0, g_xxx_0_xzzzzzz_1, g_xxx_0_yyyyyyy_0, g_xxx_0_yyyyyyy_1, g_xxx_0_yyyyyyz_0, g_xxx_0_yyyyyyz_1, g_xxx_0_yyyyyzz_0, g_xxx_0_yyyyyzz_1, g_xxx_0_yyyyzzz_0, g_xxx_0_yyyyzzz_1, g_xxx_0_yyyzzzz_0, g_xxx_0_yyyzzzz_1, g_xxx_0_yyzzzzz_0, g_xxx_0_yyzzzzz_1, g_xxx_0_yzzzzzz_0, g_xxx_0_yzzzzzz_1, g_xxx_0_zzzzzzz_0, g_xxx_0_zzzzzzz_1, g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxy_1, g_xxxx_0_xxxxxxz_1, g_xxxx_0_xxxxxy_1, g_xxxx_0_xxxxxyy_1, g_xxxx_0_xxxxxyz_1, g_xxxx_0_xxxxxz_1, g_xxxx_0_xxxxxzz_1, g_xxxx_0_xxxxyy_1, g_xxxx_0_xxxxyyy_1, g_xxxx_0_xxxxyyz_1, g_xxxx_0_xxxxyz_1, g_xxxx_0_xxxxyzz_1, g_xxxx_0_xxxxzz_1, g_xxxx_0_xxxxzzz_1, g_xxxx_0_xxxyyy_1, g_xxxx_0_xxxyyyy_1, g_xxxx_0_xxxyyyz_1, g_xxxx_0_xxxyyz_1, g_xxxx_0_xxxyyzz_1, g_xxxx_0_xxxyzz_1, g_xxxx_0_xxxyzzz_1, g_xxxx_0_xxxzzz_1, g_xxxx_0_xxxzzzz_1, g_xxxx_0_xxyyyy_1, g_xxxx_0_xxyyyyy_1, g_xxxx_0_xxyyyyz_1, g_xxxx_0_xxyyyz_1, g_xxxx_0_xxyyyzz_1, g_xxxx_0_xxyyzz_1, g_xxxx_0_xxyyzzz_1, g_xxxx_0_xxyzzz_1, g_xxxx_0_xxyzzzz_1, g_xxxx_0_xxzzzz_1, g_xxxx_0_xxzzzzz_1, g_xxxx_0_xyyyyy_1, g_xxxx_0_xyyyyyy_1, g_xxxx_0_xyyyyyz_1, g_xxxx_0_xyyyyz_1, g_xxxx_0_xyyyyzz_1, g_xxxx_0_xyyyzz_1, g_xxxx_0_xyyyzzz_1, g_xxxx_0_xyyzzz_1, g_xxxx_0_xyyzzzz_1, g_xxxx_0_xyzzzz_1, g_xxxx_0_xyzzzzz_1, g_xxxx_0_xzzzzz_1, g_xxxx_0_xzzzzzz_1, g_xxxx_0_yyyyyy_1, g_xxxx_0_yyyyyyy_1, g_xxxx_0_yyyyyyz_1, g_xxxx_0_yyyyyz_1, g_xxxx_0_yyyyyzz_1, g_xxxx_0_yyyyzz_1, g_xxxx_0_yyyyzzz_1, g_xxxx_0_yyyzzz_1, g_xxxx_0_yyyzzzz_1, g_xxxx_0_yyzzzz_1, g_xxxx_0_yyzzzzz_1, g_xxxx_0_yzzzzz_1, g_xxxx_0_yzzzzzz_1, g_xxxx_0_zzzzzz_1, g_xxxx_0_zzzzzzz_1, g_xxxxx_0_xxxxxxx_0, g_xxxxx_0_xxxxxxy_0, g_xxxxx_0_xxxxxxz_0, g_xxxxx_0_xxxxxyy_0, g_xxxxx_0_xxxxxyz_0, g_xxxxx_0_xxxxxzz_0, g_xxxxx_0_xxxxyyy_0, g_xxxxx_0_xxxxyyz_0, g_xxxxx_0_xxxxyzz_0, g_xxxxx_0_xxxxzzz_0, g_xxxxx_0_xxxyyyy_0, g_xxxxx_0_xxxyyyz_0, g_xxxxx_0_xxxyyzz_0, g_xxxxx_0_xxxyzzz_0, g_xxxxx_0_xxxzzzz_0, g_xxxxx_0_xxyyyyy_0, g_xxxxx_0_xxyyyyz_0, g_xxxxx_0_xxyyyzz_0, g_xxxxx_0_xxyyzzz_0, g_xxxxx_0_xxyzzzz_0, g_xxxxx_0_xxzzzzz_0, g_xxxxx_0_xyyyyyy_0, g_xxxxx_0_xyyyyyz_0, g_xxxxx_0_xyyyyzz_0, g_xxxxx_0_xyyyzzz_0, g_xxxxx_0_xyyzzzz_0, g_xxxxx_0_xyzzzzz_0, g_xxxxx_0_xzzzzzz_0, g_xxxxx_0_yyyyyyy_0, g_xxxxx_0_yyyyyyz_0, g_xxxxx_0_yyyyyzz_0, g_xxxxx_0_yyyyzzz_0, g_xxxxx_0_yyyzzzz_0, g_xxxxx_0_yyzzzzz_0, g_xxxxx_0_yzzzzzz_0, g_xxxxx_0_zzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxx_0_xxxxxxx_0[i] = 4.0 * g_xxx_0_xxxxxxx_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxx_1[i] * fz_be_0 + 7.0 * g_xxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxx_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxxy_0[i] = 4.0 * g_xxx_0_xxxxxxy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxxz_0[i] = 4.0 * g_xxx_0_xxxxxxz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxyy_0[i] = 4.0 * g_xxx_0_xxxxxyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxyz_0[i] = 4.0 * g_xxx_0_xxxxxyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxzz_0[i] = 4.0 * g_xxx_0_xxxxxzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyyy_0[i] = 4.0 * g_xxx_0_xxxxyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyyz_0[i] = 4.0 * g_xxx_0_xxxxyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyzz_0[i] = 4.0 * g_xxx_0_xxxxyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxzzz_0[i] = 4.0 * g_xxx_0_xxxxzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyyy_0[i] = 4.0 * g_xxx_0_xxxyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyyz_0[i] = 4.0 * g_xxx_0_xxxyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyzz_0[i] = 4.0 * g_xxx_0_xxxyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyzzz_0[i] = 4.0 * g_xxx_0_xxxyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxzzzz_0[i] = 4.0 * g_xxx_0_xxxzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyyy_0[i] = 4.0 * g_xxx_0_xxyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyyz_0[i] = 4.0 * g_xxx_0_xxyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyzz_0[i] = 4.0 * g_xxx_0_xxyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyzzz_0[i] = 4.0 * g_xxx_0_xxyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyzzzz_0[i] = 4.0 * g_xxx_0_xxyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxzzzzz_0[i] = 4.0 * g_xxx_0_xxzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyyy_0[i] = 4.0 * g_xxx_0_xyyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyyy_1[i] * fz_be_0 + g_xxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyyz_0[i] = 4.0 * g_xxx_0_xyyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyyz_1[i] * fz_be_0 + g_xxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyzz_0[i] = 4.0 * g_xxx_0_xyyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyzz_1[i] * fz_be_0 + g_xxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyzzz_0[i] = 4.0 * g_xxx_0_xyyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyzzz_1[i] * fz_be_0 + g_xxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyzzzz_0[i] = 4.0 * g_xxx_0_xyyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyzzzz_1[i] * fz_be_0 + g_xxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyzzzzz_0[i] = 4.0 * g_xxx_0_xyzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyzzzzz_1[i] * fz_be_0 + g_xxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xzzzzzz_0[i] = 4.0 * g_xxx_0_xzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xzzzzzz_1[i] * fz_be_0 + g_xxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyyy_0[i] = 4.0 * g_xxx_0_yyyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyyy_1[i] * fz_be_0 + g_xxxx_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyyz_0[i] = 4.0 * g_xxx_0_yyyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyyz_1[i] * fz_be_0 + g_xxxx_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyzz_0[i] = 4.0 * g_xxx_0_yyyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyzz_1[i] * fz_be_0 + g_xxxx_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyzzz_0[i] = 4.0 * g_xxx_0_yyyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyzzz_1[i] * fz_be_0 + g_xxxx_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyzzzz_0[i] = 4.0 * g_xxx_0_yyyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyzzzz_1[i] * fz_be_0 + g_xxxx_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyzzzzz_0[i] = 4.0 * g_xxx_0_yyzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyzzzzz_1[i] * fz_be_0 + g_xxxx_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yzzzzzz_0[i] = 4.0 * g_xxx_0_yzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yzzzzzz_1[i] * fz_be_0 + g_xxxx_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_zzzzzzz_0[i] = 4.0 * g_xxx_0_zzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_zzzzzzz_1[i] * fz_be_0 + g_xxxx_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 36-72 components of targeted buffer : HSK

    auto g_xxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 36);

    auto g_xxxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 37);

    auto g_xxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 38);

    auto g_xxxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 39);

    auto g_xxxxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 40);

    auto g_xxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 41);

    auto g_xxxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 42);

    auto g_xxxxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 43);

    auto g_xxxxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 44);

    auto g_xxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 45);

    auto g_xxxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 46);

    auto g_xxxxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 47);

    auto g_xxxxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 48);

    auto g_xxxxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 49);

    auto g_xxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 50);

    auto g_xxxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 51);

    auto g_xxxxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 52);

    auto g_xxxxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 53);

    auto g_xxxxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 54);

    auto g_xxxxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 55);

    auto g_xxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 56);

    auto g_xxxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 57);

    auto g_xxxxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 58);

    auto g_xxxxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 59);

    auto g_xxxxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 60);

    auto g_xxxxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 61);

    auto g_xxxxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 62);

    auto g_xxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 63);

    auto g_xxxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 64);

    auto g_xxxxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 65);

    auto g_xxxxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 66);

    auto g_xxxxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 67);

    auto g_xxxxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 68);

    auto g_xxxxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 69);

    auto g_xxxxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 70);

    auto g_xxxxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 71);

    #pragma omp simd aligned(g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxy_1, g_xxxx_0_xxxxxxz_1, g_xxxx_0_xxxxxy_1, g_xxxx_0_xxxxxyy_1, g_xxxx_0_xxxxxyz_1, g_xxxx_0_xxxxxz_1, g_xxxx_0_xxxxxzz_1, g_xxxx_0_xxxxyy_1, g_xxxx_0_xxxxyyy_1, g_xxxx_0_xxxxyyz_1, g_xxxx_0_xxxxyz_1, g_xxxx_0_xxxxyzz_1, g_xxxx_0_xxxxzz_1, g_xxxx_0_xxxxzzz_1, g_xxxx_0_xxxyyy_1, g_xxxx_0_xxxyyyy_1, g_xxxx_0_xxxyyyz_1, g_xxxx_0_xxxyyz_1, g_xxxx_0_xxxyyzz_1, g_xxxx_0_xxxyzz_1, g_xxxx_0_xxxyzzz_1, g_xxxx_0_xxxzzz_1, g_xxxx_0_xxxzzzz_1, g_xxxx_0_xxyyyy_1, g_xxxx_0_xxyyyyy_1, g_xxxx_0_xxyyyyz_1, g_xxxx_0_xxyyyz_1, g_xxxx_0_xxyyyzz_1, g_xxxx_0_xxyyzz_1, g_xxxx_0_xxyyzzz_1, g_xxxx_0_xxyzzz_1, g_xxxx_0_xxyzzzz_1, g_xxxx_0_xxzzzz_1, g_xxxx_0_xxzzzzz_1, g_xxxx_0_xyyyyy_1, g_xxxx_0_xyyyyyy_1, g_xxxx_0_xyyyyyz_1, g_xxxx_0_xyyyyz_1, g_xxxx_0_xyyyyzz_1, g_xxxx_0_xyyyzz_1, g_xxxx_0_xyyyzzz_1, g_xxxx_0_xyyzzz_1, g_xxxx_0_xyyzzzz_1, g_xxxx_0_xyzzzz_1, g_xxxx_0_xyzzzzz_1, g_xxxx_0_xzzzzz_1, g_xxxx_0_xzzzzzz_1, g_xxxx_0_yyyyyy_1, g_xxxx_0_yyyyyyy_1, g_xxxx_0_yyyyyyz_1, g_xxxx_0_yyyyyz_1, g_xxxx_0_yyyyyzz_1, g_xxxx_0_yyyyzz_1, g_xxxx_0_yyyyzzz_1, g_xxxx_0_yyyzzz_1, g_xxxx_0_yyyzzzz_1, g_xxxx_0_yyzzzz_1, g_xxxx_0_yyzzzzz_1, g_xxxx_0_yzzzzz_1, g_xxxx_0_yzzzzzz_1, g_xxxx_0_zzzzzz_1, g_xxxx_0_zzzzzzz_1, g_xxxxy_0_xxxxxxx_0, g_xxxxy_0_xxxxxxy_0, g_xxxxy_0_xxxxxxz_0, g_xxxxy_0_xxxxxyy_0, g_xxxxy_0_xxxxxyz_0, g_xxxxy_0_xxxxxzz_0, g_xxxxy_0_xxxxyyy_0, g_xxxxy_0_xxxxyyz_0, g_xxxxy_0_xxxxyzz_0, g_xxxxy_0_xxxxzzz_0, g_xxxxy_0_xxxyyyy_0, g_xxxxy_0_xxxyyyz_0, g_xxxxy_0_xxxyyzz_0, g_xxxxy_0_xxxyzzz_0, g_xxxxy_0_xxxzzzz_0, g_xxxxy_0_xxyyyyy_0, g_xxxxy_0_xxyyyyz_0, g_xxxxy_0_xxyyyzz_0, g_xxxxy_0_xxyyzzz_0, g_xxxxy_0_xxyzzzz_0, g_xxxxy_0_xxzzzzz_0, g_xxxxy_0_xyyyyyy_0, g_xxxxy_0_xyyyyyz_0, g_xxxxy_0_xyyyyzz_0, g_xxxxy_0_xyyyzzz_0, g_xxxxy_0_xyyzzzz_0, g_xxxxy_0_xyzzzzz_0, g_xxxxy_0_xzzzzzz_0, g_xxxxy_0_yyyyyyy_0, g_xxxxy_0_yyyyyyz_0, g_xxxxy_0_yyyyyzz_0, g_xxxxy_0_yyyyzzz_0, g_xxxxy_0_yyyzzzz_0, g_xxxxy_0_yyzzzzz_0, g_xxxxy_0_yzzzzzz_0, g_xxxxy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxy_0_xxxxxxx_0[i] = g_xxxx_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxxy_0[i] = g_xxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxxz_0[i] = g_xxxx_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxyy_0[i] = 2.0 * g_xxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxyz_0[i] = g_xxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxzz_0[i] = g_xxxx_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyyy_0[i] = 3.0 * g_xxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyyz_0[i] = 2.0 * g_xxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyzz_0[i] = g_xxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxzzz_0[i] = g_xxxx_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyyy_0[i] = 4.0 * g_xxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyyz_0[i] = 3.0 * g_xxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyzz_0[i] = 2.0 * g_xxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyzzz_0[i] = g_xxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxzzzz_0[i] = g_xxxx_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyyy_0[i] = 5.0 * g_xxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyyz_0[i] = 4.0 * g_xxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyzz_0[i] = 3.0 * g_xxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyzzz_0[i] = 2.0 * g_xxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyzzzz_0[i] = g_xxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxzzzzz_0[i] = g_xxxx_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyyy_0[i] = 6.0 * g_xxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyyz_0[i] = 5.0 * g_xxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyzz_0[i] = 4.0 * g_xxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyzzz_0[i] = 3.0 * g_xxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyzzzz_0[i] = 2.0 * g_xxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyzzzzz_0[i] = g_xxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xzzzzzz_0[i] = g_xxxx_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyyy_0[i] = 7.0 * g_xxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxx_0_yyyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyyz_0[i] = 6.0 * g_xxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyzz_0[i] = 5.0 * g_xxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyzzz_0[i] = 4.0 * g_xxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyzzzz_0[i] = 3.0 * g_xxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyzzzzz_0[i] = 2.0 * g_xxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yzzzzzz_0[i] = g_xxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxx_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_zzzzzzz_0[i] = g_xxxx_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 72-108 components of targeted buffer : HSK

    auto g_xxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 72);

    auto g_xxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 73);

    auto g_xxxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 74);

    auto g_xxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 75);

    auto g_xxxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 76);

    auto g_xxxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 77);

    auto g_xxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 78);

    auto g_xxxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 79);

    auto g_xxxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 80);

    auto g_xxxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 81);

    auto g_xxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 82);

    auto g_xxxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 83);

    auto g_xxxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 84);

    auto g_xxxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 85);

    auto g_xxxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 86);

    auto g_xxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 87);

    auto g_xxxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 88);

    auto g_xxxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 89);

    auto g_xxxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 90);

    auto g_xxxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 91);

    auto g_xxxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 92);

    auto g_xxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 93);

    auto g_xxxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 94);

    auto g_xxxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 95);

    auto g_xxxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 96);

    auto g_xxxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 97);

    auto g_xxxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 98);

    auto g_xxxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 99);

    auto g_xxxxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 100);

    auto g_xxxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 101);

    auto g_xxxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 102);

    auto g_xxxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 103);

    auto g_xxxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 104);

    auto g_xxxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 105);

    auto g_xxxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 106);

    auto g_xxxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 107);

    #pragma omp simd aligned(g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxy_1, g_xxxx_0_xxxxxxz_1, g_xxxx_0_xxxxxy_1, g_xxxx_0_xxxxxyy_1, g_xxxx_0_xxxxxyz_1, g_xxxx_0_xxxxxz_1, g_xxxx_0_xxxxxzz_1, g_xxxx_0_xxxxyy_1, g_xxxx_0_xxxxyyy_1, g_xxxx_0_xxxxyyz_1, g_xxxx_0_xxxxyz_1, g_xxxx_0_xxxxyzz_1, g_xxxx_0_xxxxzz_1, g_xxxx_0_xxxxzzz_1, g_xxxx_0_xxxyyy_1, g_xxxx_0_xxxyyyy_1, g_xxxx_0_xxxyyyz_1, g_xxxx_0_xxxyyz_1, g_xxxx_0_xxxyyzz_1, g_xxxx_0_xxxyzz_1, g_xxxx_0_xxxyzzz_1, g_xxxx_0_xxxzzz_1, g_xxxx_0_xxxzzzz_1, g_xxxx_0_xxyyyy_1, g_xxxx_0_xxyyyyy_1, g_xxxx_0_xxyyyyz_1, g_xxxx_0_xxyyyz_1, g_xxxx_0_xxyyyzz_1, g_xxxx_0_xxyyzz_1, g_xxxx_0_xxyyzzz_1, g_xxxx_0_xxyzzz_1, g_xxxx_0_xxyzzzz_1, g_xxxx_0_xxzzzz_1, g_xxxx_0_xxzzzzz_1, g_xxxx_0_xyyyyy_1, g_xxxx_0_xyyyyyy_1, g_xxxx_0_xyyyyyz_1, g_xxxx_0_xyyyyz_1, g_xxxx_0_xyyyyzz_1, g_xxxx_0_xyyyzz_1, g_xxxx_0_xyyyzzz_1, g_xxxx_0_xyyzzz_1, g_xxxx_0_xyyzzzz_1, g_xxxx_0_xyzzzz_1, g_xxxx_0_xyzzzzz_1, g_xxxx_0_xzzzzz_1, g_xxxx_0_xzzzzzz_1, g_xxxx_0_yyyyyy_1, g_xxxx_0_yyyyyyy_1, g_xxxx_0_yyyyyyz_1, g_xxxx_0_yyyyyz_1, g_xxxx_0_yyyyyzz_1, g_xxxx_0_yyyyzz_1, g_xxxx_0_yyyyzzz_1, g_xxxx_0_yyyzzz_1, g_xxxx_0_yyyzzzz_1, g_xxxx_0_yyzzzz_1, g_xxxx_0_yyzzzzz_1, g_xxxx_0_yzzzzz_1, g_xxxx_0_yzzzzzz_1, g_xxxx_0_zzzzzz_1, g_xxxx_0_zzzzzzz_1, g_xxxxz_0_xxxxxxx_0, g_xxxxz_0_xxxxxxy_0, g_xxxxz_0_xxxxxxz_0, g_xxxxz_0_xxxxxyy_0, g_xxxxz_0_xxxxxyz_0, g_xxxxz_0_xxxxxzz_0, g_xxxxz_0_xxxxyyy_0, g_xxxxz_0_xxxxyyz_0, g_xxxxz_0_xxxxyzz_0, g_xxxxz_0_xxxxzzz_0, g_xxxxz_0_xxxyyyy_0, g_xxxxz_0_xxxyyyz_0, g_xxxxz_0_xxxyyzz_0, g_xxxxz_0_xxxyzzz_0, g_xxxxz_0_xxxzzzz_0, g_xxxxz_0_xxyyyyy_0, g_xxxxz_0_xxyyyyz_0, g_xxxxz_0_xxyyyzz_0, g_xxxxz_0_xxyyzzz_0, g_xxxxz_0_xxyzzzz_0, g_xxxxz_0_xxzzzzz_0, g_xxxxz_0_xyyyyyy_0, g_xxxxz_0_xyyyyyz_0, g_xxxxz_0_xyyyyzz_0, g_xxxxz_0_xyyyzzz_0, g_xxxxz_0_xyyzzzz_0, g_xxxxz_0_xyzzzzz_0, g_xxxxz_0_xzzzzzz_0, g_xxxxz_0_yyyyyyy_0, g_xxxxz_0_yyyyyyz_0, g_xxxxz_0_yyyyyzz_0, g_xxxxz_0_yyyyzzz_0, g_xxxxz_0_yyyzzzz_0, g_xxxxz_0_yyzzzzz_0, g_xxxxz_0_yzzzzzz_0, g_xxxxz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxz_0_xxxxxxx_0[i] = g_xxxx_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxxy_0[i] = g_xxxx_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxxz_0[i] = g_xxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxyy_0[i] = g_xxxx_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxyz_0[i] = g_xxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxzz_0[i] = 2.0 * g_xxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyyy_0[i] = g_xxxx_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyyz_0[i] = g_xxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyzz_0[i] = 2.0 * g_xxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxzzz_0[i] = 3.0 * g_xxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyyy_0[i] = g_xxxx_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyyz_0[i] = g_xxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyzz_0[i] = 2.0 * g_xxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyzzz_0[i] = 3.0 * g_xxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxzzzz_0[i] = 4.0 * g_xxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyyy_0[i] = g_xxxx_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyyz_0[i] = g_xxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyzz_0[i] = 2.0 * g_xxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyzzz_0[i] = 3.0 * g_xxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyzzzz_0[i] = 4.0 * g_xxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxzzzzz_0[i] = 5.0 * g_xxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyyy_0[i] = g_xxxx_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyyz_0[i] = g_xxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyzz_0[i] = 2.0 * g_xxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyzzz_0[i] = 3.0 * g_xxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyzzzz_0[i] = 4.0 * g_xxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyzzzzz_0[i] = 5.0 * g_xxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xzzzzzz_0[i] = 6.0 * g_xxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xzzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyyy_0[i] = g_xxxx_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyyz_0[i] = g_xxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxx_0_yyyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyzz_0[i] = 2.0 * g_xxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyzzz_0[i] = 3.0 * g_xxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyzzzz_0[i] = 4.0 * g_xxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyzzzzz_0[i] = 5.0 * g_xxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yzzzzzz_0[i] = 6.0 * g_xxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxx_0_yzzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_zzzzzzz_0[i] = 7.0 * g_xxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxx_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 108-144 components of targeted buffer : HSK

    auto g_xxxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 108);

    auto g_xxxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 109);

    auto g_xxxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 110);

    auto g_xxxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 111);

    auto g_xxxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 112);

    auto g_xxxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 113);

    auto g_xxxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 114);

    auto g_xxxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 115);

    auto g_xxxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 116);

    auto g_xxxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 117);

    auto g_xxxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 118);

    auto g_xxxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 119);

    auto g_xxxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 120);

    auto g_xxxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 121);

    auto g_xxxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 122);

    auto g_xxxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 123);

    auto g_xxxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 124);

    auto g_xxxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 125);

    auto g_xxxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 126);

    auto g_xxxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 127);

    auto g_xxxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 128);

    auto g_xxxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 129);

    auto g_xxxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 130);

    auto g_xxxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 131);

    auto g_xxxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 132);

    auto g_xxxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 133);

    auto g_xxxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 134);

    auto g_xxxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 135);

    auto g_xxxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 136);

    auto g_xxxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 137);

    auto g_xxxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 138);

    auto g_xxxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 139);

    auto g_xxxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 140);

    auto g_xxxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 141);

    auto g_xxxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 142);

    auto g_xxxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 143);

    #pragma omp simd aligned(g_xxx_0_xxxxxxx_0, g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxz_0, g_xxx_0_xxxxxxz_1, g_xxx_0_xxxxxzz_0, g_xxx_0_xxxxxzz_1, g_xxx_0_xxxxzzz_0, g_xxx_0_xxxxzzz_1, g_xxx_0_xxxzzzz_0, g_xxx_0_xxxzzzz_1, g_xxx_0_xxzzzzz_0, g_xxx_0_xxzzzzz_1, g_xxx_0_xzzzzzz_0, g_xxx_0_xzzzzzz_1, g_xxxy_0_xxxxxxx_1, g_xxxy_0_xxxxxxz_1, g_xxxy_0_xxxxxzz_1, g_xxxy_0_xxxxzzz_1, g_xxxy_0_xxxzzzz_1, g_xxxy_0_xxzzzzz_1, g_xxxy_0_xzzzzzz_1, g_xxxyy_0_xxxxxxx_0, g_xxxyy_0_xxxxxxy_0, g_xxxyy_0_xxxxxxz_0, g_xxxyy_0_xxxxxyy_0, g_xxxyy_0_xxxxxyz_0, g_xxxyy_0_xxxxxzz_0, g_xxxyy_0_xxxxyyy_0, g_xxxyy_0_xxxxyyz_0, g_xxxyy_0_xxxxyzz_0, g_xxxyy_0_xxxxzzz_0, g_xxxyy_0_xxxyyyy_0, g_xxxyy_0_xxxyyyz_0, g_xxxyy_0_xxxyyzz_0, g_xxxyy_0_xxxyzzz_0, g_xxxyy_0_xxxzzzz_0, g_xxxyy_0_xxyyyyy_0, g_xxxyy_0_xxyyyyz_0, g_xxxyy_0_xxyyyzz_0, g_xxxyy_0_xxyyzzz_0, g_xxxyy_0_xxyzzzz_0, g_xxxyy_0_xxzzzzz_0, g_xxxyy_0_xyyyyyy_0, g_xxxyy_0_xyyyyyz_0, g_xxxyy_0_xyyyyzz_0, g_xxxyy_0_xyyyzzz_0, g_xxxyy_0_xyyzzzz_0, g_xxxyy_0_xyzzzzz_0, g_xxxyy_0_xzzzzzz_0, g_xxxyy_0_yyyyyyy_0, g_xxxyy_0_yyyyyyz_0, g_xxxyy_0_yyyyyzz_0, g_xxxyy_0_yyyyzzz_0, g_xxxyy_0_yyyzzzz_0, g_xxxyy_0_yyzzzzz_0, g_xxxyy_0_yzzzzzz_0, g_xxxyy_0_zzzzzzz_0, g_xxyy_0_xxxxxxy_1, g_xxyy_0_xxxxxy_1, g_xxyy_0_xxxxxyy_1, g_xxyy_0_xxxxxyz_1, g_xxyy_0_xxxxyy_1, g_xxyy_0_xxxxyyy_1, g_xxyy_0_xxxxyyz_1, g_xxyy_0_xxxxyz_1, g_xxyy_0_xxxxyzz_1, g_xxyy_0_xxxyyy_1, g_xxyy_0_xxxyyyy_1, g_xxyy_0_xxxyyyz_1, g_xxyy_0_xxxyyz_1, g_xxyy_0_xxxyyzz_1, g_xxyy_0_xxxyzz_1, g_xxyy_0_xxxyzzz_1, g_xxyy_0_xxyyyy_1, g_xxyy_0_xxyyyyy_1, g_xxyy_0_xxyyyyz_1, g_xxyy_0_xxyyyz_1, g_xxyy_0_xxyyyzz_1, g_xxyy_0_xxyyzz_1, g_xxyy_0_xxyyzzz_1, g_xxyy_0_xxyzzz_1, g_xxyy_0_xxyzzzz_1, g_xxyy_0_xyyyyy_1, g_xxyy_0_xyyyyyy_1, g_xxyy_0_xyyyyyz_1, g_xxyy_0_xyyyyz_1, g_xxyy_0_xyyyyzz_1, g_xxyy_0_xyyyzz_1, g_xxyy_0_xyyyzzz_1, g_xxyy_0_xyyzzz_1, g_xxyy_0_xyyzzzz_1, g_xxyy_0_xyzzzz_1, g_xxyy_0_xyzzzzz_1, g_xxyy_0_yyyyyy_1, g_xxyy_0_yyyyyyy_1, g_xxyy_0_yyyyyyz_1, g_xxyy_0_yyyyyz_1, g_xxyy_0_yyyyyzz_1, g_xxyy_0_yyyyzz_1, g_xxyy_0_yyyyzzz_1, g_xxyy_0_yyyzzz_1, g_xxyy_0_yyyzzzz_1, g_xxyy_0_yyzzzz_1, g_xxyy_0_yyzzzzz_1, g_xxyy_0_yzzzzz_1, g_xxyy_0_yzzzzzz_1, g_xxyy_0_zzzzzzz_1, g_xyy_0_xxxxxxy_0, g_xyy_0_xxxxxxy_1, g_xyy_0_xxxxxyy_0, g_xyy_0_xxxxxyy_1, g_xyy_0_xxxxxyz_0, g_xyy_0_xxxxxyz_1, g_xyy_0_xxxxyyy_0, g_xyy_0_xxxxyyy_1, g_xyy_0_xxxxyyz_0, g_xyy_0_xxxxyyz_1, g_xyy_0_xxxxyzz_0, g_xyy_0_xxxxyzz_1, g_xyy_0_xxxyyyy_0, g_xyy_0_xxxyyyy_1, g_xyy_0_xxxyyyz_0, g_xyy_0_xxxyyyz_1, g_xyy_0_xxxyyzz_0, g_xyy_0_xxxyyzz_1, g_xyy_0_xxxyzzz_0, g_xyy_0_xxxyzzz_1, g_xyy_0_xxyyyyy_0, g_xyy_0_xxyyyyy_1, g_xyy_0_xxyyyyz_0, g_xyy_0_xxyyyyz_1, g_xyy_0_xxyyyzz_0, g_xyy_0_xxyyyzz_1, g_xyy_0_xxyyzzz_0, g_xyy_0_xxyyzzz_1, g_xyy_0_xxyzzzz_0, g_xyy_0_xxyzzzz_1, g_xyy_0_xyyyyyy_0, g_xyy_0_xyyyyyy_1, g_xyy_0_xyyyyyz_0, g_xyy_0_xyyyyyz_1, g_xyy_0_xyyyyzz_0, g_xyy_0_xyyyyzz_1, g_xyy_0_xyyyzzz_0, g_xyy_0_xyyyzzz_1, g_xyy_0_xyyzzzz_0, g_xyy_0_xyyzzzz_1, g_xyy_0_xyzzzzz_0, g_xyy_0_xyzzzzz_1, g_xyy_0_yyyyyyy_0, g_xyy_0_yyyyyyy_1, g_xyy_0_yyyyyyz_0, g_xyy_0_yyyyyyz_1, g_xyy_0_yyyyyzz_0, g_xyy_0_yyyyyzz_1, g_xyy_0_yyyyzzz_0, g_xyy_0_yyyyzzz_1, g_xyy_0_yyyzzzz_0, g_xyy_0_yyyzzzz_1, g_xyy_0_yyzzzzz_0, g_xyy_0_yyzzzzz_1, g_xyy_0_yzzzzzz_0, g_xyy_0_yzzzzzz_1, g_xyy_0_zzzzzzz_0, g_xyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyy_0_xxxxxxx_0[i] = g_xxx_0_xxxxxxx_0[i] * fbe_0 - g_xxx_0_xxxxxxx_1[i] * fz_be_0 + g_xxxy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxyy_0_xxxxxxy_0[i] = 2.0 * g_xyy_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxxz_0[i] = g_xxx_0_xxxxxxz_0[i] * fbe_0 - g_xxx_0_xxxxxxz_1[i] * fz_be_0 + g_xxxy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxyy_0_xxxxxyy_0[i] = 2.0 * g_xyy_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxyz_0[i] = 2.0 * g_xyy_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxzz_0[i] = g_xxx_0_xxxxxzz_0[i] * fbe_0 - g_xxx_0_xxxxxzz_1[i] * fz_be_0 + g_xxxy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxyy_0_xxxxyyy_0[i] = 2.0 * g_xyy_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxyyz_0[i] = 2.0 * g_xyy_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxyzz_0[i] = 2.0 * g_xyy_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxzzz_0[i] = g_xxx_0_xxxxzzz_0[i] * fbe_0 - g_xxx_0_xxxxzzz_1[i] * fz_be_0 + g_xxxy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxyy_0_xxxyyyy_0[i] = 2.0 * g_xyy_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxyyyz_0[i] = 2.0 * g_xyy_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxyyzz_0[i] = 2.0 * g_xyy_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxyzzz_0[i] = 2.0 * g_xyy_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxzzzz_0[i] = g_xxx_0_xxxzzzz_0[i] * fbe_0 - g_xxx_0_xxxzzzz_1[i] * fz_be_0 + g_xxxy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxyy_0_xxyyyyy_0[i] = 2.0 * g_xyy_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxyyyyz_0[i] = 2.0 * g_xyy_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxyyyzz_0[i] = 2.0 * g_xyy_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxyyzzz_0[i] = 2.0 * g_xyy_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxyzzzz_0[i] = 2.0 * g_xyy_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxzzzzz_0[i] = g_xxx_0_xxzzzzz_0[i] * fbe_0 - g_xxx_0_xxzzzzz_1[i] * fz_be_0 + g_xxxy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxyy_0_xyyyyyy_0[i] = 2.0 * g_xyy_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xyyyyyz_0[i] = 2.0 * g_xyy_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyyz_1[i] * fz_be_0 + g_xxyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xyyyyzz_0[i] = 2.0 * g_xyy_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyzz_1[i] * fz_be_0 + g_xxyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xyyyzzz_0[i] = 2.0 * g_xyy_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyzzz_1[i] * fz_be_0 + g_xxyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xyyzzzz_0[i] = 2.0 * g_xyy_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyzzzz_1[i] * fz_be_0 + g_xxyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xyzzzzz_0[i] = 2.0 * g_xyy_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyzzzzz_1[i] * fz_be_0 + g_xxyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xzzzzzz_0[i] = g_xxx_0_xzzzzzz_0[i] * fbe_0 - g_xxx_0_xzzzzzz_1[i] * fz_be_0 + g_xxxy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxyy_0_yyyyyyy_0[i] = 2.0 * g_xyy_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyyy_1[i] * fz_be_0 + g_xxyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_yyyyyyz_0[i] = 2.0 * g_xyy_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyyz_1[i] * fz_be_0 + g_xxyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_yyyyyzz_0[i] = 2.0 * g_xyy_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyzz_1[i] * fz_be_0 + g_xxyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_yyyyzzz_0[i] = 2.0 * g_xyy_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyzzz_1[i] * fz_be_0 + g_xxyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_yyyzzzz_0[i] = 2.0 * g_xyy_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyzzzz_1[i] * fz_be_0 + g_xxyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_yyzzzzz_0[i] = 2.0 * g_xyy_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyzzzzz_1[i] * fz_be_0 + g_xxyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_yzzzzzz_0[i] = 2.0 * g_xyy_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yzzzzzz_1[i] * fz_be_0 + g_xxyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_zzzzzzz_0[i] = 2.0 * g_xyy_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_zzzzzzz_1[i] * fz_be_0 + g_xxyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 144-180 components of targeted buffer : HSK

    auto g_xxxyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 144);

    auto g_xxxyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 145);

    auto g_xxxyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 146);

    auto g_xxxyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 147);

    auto g_xxxyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 148);

    auto g_xxxyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 149);

    auto g_xxxyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 150);

    auto g_xxxyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 151);

    auto g_xxxyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 152);

    auto g_xxxyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 153);

    auto g_xxxyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 154);

    auto g_xxxyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 155);

    auto g_xxxyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 156);

    auto g_xxxyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 157);

    auto g_xxxyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 158);

    auto g_xxxyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 159);

    auto g_xxxyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 160);

    auto g_xxxyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 161);

    auto g_xxxyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 162);

    auto g_xxxyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 163);

    auto g_xxxyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 164);

    auto g_xxxyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 165);

    auto g_xxxyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 166);

    auto g_xxxyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 167);

    auto g_xxxyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 168);

    auto g_xxxyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 169);

    auto g_xxxyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 170);

    auto g_xxxyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 171);

    auto g_xxxyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 172);

    auto g_xxxyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 173);

    auto g_xxxyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 174);

    auto g_xxxyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 175);

    auto g_xxxyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 176);

    auto g_xxxyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 177);

    auto g_xxxyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 178);

    auto g_xxxyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 179);

    #pragma omp simd aligned(g_xxxy_0_xxxxxxy_1, g_xxxy_0_xxxxxyy_1, g_xxxy_0_xxxxyyy_1, g_xxxy_0_xxxyyyy_1, g_xxxy_0_xxyyyyy_1, g_xxxy_0_xyyyyyy_1, g_xxxy_0_yyyyyyy_1, g_xxxyz_0_xxxxxxx_0, g_xxxyz_0_xxxxxxy_0, g_xxxyz_0_xxxxxxz_0, g_xxxyz_0_xxxxxyy_0, g_xxxyz_0_xxxxxyz_0, g_xxxyz_0_xxxxxzz_0, g_xxxyz_0_xxxxyyy_0, g_xxxyz_0_xxxxyyz_0, g_xxxyz_0_xxxxyzz_0, g_xxxyz_0_xxxxzzz_0, g_xxxyz_0_xxxyyyy_0, g_xxxyz_0_xxxyyyz_0, g_xxxyz_0_xxxyyzz_0, g_xxxyz_0_xxxyzzz_0, g_xxxyz_0_xxxzzzz_0, g_xxxyz_0_xxyyyyy_0, g_xxxyz_0_xxyyyyz_0, g_xxxyz_0_xxyyyzz_0, g_xxxyz_0_xxyyzzz_0, g_xxxyz_0_xxyzzzz_0, g_xxxyz_0_xxzzzzz_0, g_xxxyz_0_xyyyyyy_0, g_xxxyz_0_xyyyyyz_0, g_xxxyz_0_xyyyyzz_0, g_xxxyz_0_xyyyzzz_0, g_xxxyz_0_xyyzzzz_0, g_xxxyz_0_xyzzzzz_0, g_xxxyz_0_xzzzzzz_0, g_xxxyz_0_yyyyyyy_0, g_xxxyz_0_yyyyyyz_0, g_xxxyz_0_yyyyyzz_0, g_xxxyz_0_yyyyzzz_0, g_xxxyz_0_yyyzzzz_0, g_xxxyz_0_yyzzzzz_0, g_xxxyz_0_yzzzzzz_0, g_xxxyz_0_zzzzzzz_0, g_xxxz_0_xxxxxxx_1, g_xxxz_0_xxxxxxz_1, g_xxxz_0_xxxxxyz_1, g_xxxz_0_xxxxxz_1, g_xxxz_0_xxxxxzz_1, g_xxxz_0_xxxxyyz_1, g_xxxz_0_xxxxyz_1, g_xxxz_0_xxxxyzz_1, g_xxxz_0_xxxxzz_1, g_xxxz_0_xxxxzzz_1, g_xxxz_0_xxxyyyz_1, g_xxxz_0_xxxyyz_1, g_xxxz_0_xxxyyzz_1, g_xxxz_0_xxxyzz_1, g_xxxz_0_xxxyzzz_1, g_xxxz_0_xxxzzz_1, g_xxxz_0_xxxzzzz_1, g_xxxz_0_xxyyyyz_1, g_xxxz_0_xxyyyz_1, g_xxxz_0_xxyyyzz_1, g_xxxz_0_xxyyzz_1, g_xxxz_0_xxyyzzz_1, g_xxxz_0_xxyzzz_1, g_xxxz_0_xxyzzzz_1, g_xxxz_0_xxzzzz_1, g_xxxz_0_xxzzzzz_1, g_xxxz_0_xyyyyyz_1, g_xxxz_0_xyyyyz_1, g_xxxz_0_xyyyyzz_1, g_xxxz_0_xyyyzz_1, g_xxxz_0_xyyyzzz_1, g_xxxz_0_xyyzzz_1, g_xxxz_0_xyyzzzz_1, g_xxxz_0_xyzzzz_1, g_xxxz_0_xyzzzzz_1, g_xxxz_0_xzzzzz_1, g_xxxz_0_xzzzzzz_1, g_xxxz_0_yyyyyyz_1, g_xxxz_0_yyyyyz_1, g_xxxz_0_yyyyyzz_1, g_xxxz_0_yyyyzz_1, g_xxxz_0_yyyyzzz_1, g_xxxz_0_yyyzzz_1, g_xxxz_0_yyyzzzz_1, g_xxxz_0_yyzzzz_1, g_xxxz_0_yyzzzzz_1, g_xxxz_0_yzzzzz_1, g_xxxz_0_yzzzzzz_1, g_xxxz_0_zzzzzz_1, g_xxxz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyz_0_xxxxxxx_0[i] = g_xxxz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxxy_0[i] = g_xxxy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxxxz_0[i] = g_xxxz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxyy_0[i] = g_xxxy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxxyz_0[i] = g_xxxz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxzz_0[i] = g_xxxz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxyyy_0[i] = g_xxxy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxyyz_0[i] = 2.0 * g_xxxz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxyzz_0[i] = g_xxxz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxzzz_0[i] = g_xxxz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyyyy_0[i] = g_xxxy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxyyyz_0[i] = 3.0 * g_xxxz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyyzz_0[i] = 2.0 * g_xxxz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyzzz_0[i] = g_xxxz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxzzzz_0[i] = g_xxxz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyyyy_0[i] = g_xxxy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxyyyyz_0[i] = 4.0 * g_xxxz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyyzz_0[i] = 3.0 * g_xxxz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyzzz_0[i] = 2.0 * g_xxxz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyzzzz_0[i] = g_xxxz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxzzzzz_0[i] = g_xxxz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyyyy_0[i] = g_xxxy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xyyyyyz_0[i] = 5.0 * g_xxxz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyyzz_0[i] = 4.0 * g_xxxz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyzzz_0[i] = 3.0 * g_xxxz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyzzzz_0[i] = 2.0 * g_xxxz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyzzzzz_0[i] = g_xxxz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xzzzzzz_0[i] = g_xxxz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyyyy_0[i] = g_xxxy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_yyyyyyz_0[i] = 6.0 * g_xxxz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyyzz_0[i] = 5.0 * g_xxxz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyzzz_0[i] = 4.0 * g_xxxz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyzzzz_0[i] = 3.0 * g_xxxz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyzzzzz_0[i] = 2.0 * g_xxxz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yzzzzzz_0[i] = g_xxxz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_zzzzzzz_0[i] = g_xxxz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 180-216 components of targeted buffer : HSK

    auto g_xxxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 180);

    auto g_xxxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 181);

    auto g_xxxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 182);

    auto g_xxxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 183);

    auto g_xxxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 184);

    auto g_xxxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 185);

    auto g_xxxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 186);

    auto g_xxxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 187);

    auto g_xxxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 188);

    auto g_xxxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 189);

    auto g_xxxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 190);

    auto g_xxxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 191);

    auto g_xxxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 192);

    auto g_xxxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 193);

    auto g_xxxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 194);

    auto g_xxxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 195);

    auto g_xxxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 196);

    auto g_xxxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 197);

    auto g_xxxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 198);

    auto g_xxxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 199);

    auto g_xxxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 200);

    auto g_xxxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 201);

    auto g_xxxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 202);

    auto g_xxxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 203);

    auto g_xxxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 204);

    auto g_xxxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 205);

    auto g_xxxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 206);

    auto g_xxxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 207);

    auto g_xxxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 208);

    auto g_xxxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 209);

    auto g_xxxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 210);

    auto g_xxxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 211);

    auto g_xxxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 212);

    auto g_xxxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 213);

    auto g_xxxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 214);

    auto g_xxxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 215);

    #pragma omp simd aligned(g_xxx_0_xxxxxxx_0, g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxy_0, g_xxx_0_xxxxxxy_1, g_xxx_0_xxxxxyy_0, g_xxx_0_xxxxxyy_1, g_xxx_0_xxxxyyy_0, g_xxx_0_xxxxyyy_1, g_xxx_0_xxxyyyy_0, g_xxx_0_xxxyyyy_1, g_xxx_0_xxyyyyy_0, g_xxx_0_xxyyyyy_1, g_xxx_0_xyyyyyy_0, g_xxx_0_xyyyyyy_1, g_xxxz_0_xxxxxxx_1, g_xxxz_0_xxxxxxy_1, g_xxxz_0_xxxxxyy_1, g_xxxz_0_xxxxyyy_1, g_xxxz_0_xxxyyyy_1, g_xxxz_0_xxyyyyy_1, g_xxxz_0_xyyyyyy_1, g_xxxzz_0_xxxxxxx_0, g_xxxzz_0_xxxxxxy_0, g_xxxzz_0_xxxxxxz_0, g_xxxzz_0_xxxxxyy_0, g_xxxzz_0_xxxxxyz_0, g_xxxzz_0_xxxxxzz_0, g_xxxzz_0_xxxxyyy_0, g_xxxzz_0_xxxxyyz_0, g_xxxzz_0_xxxxyzz_0, g_xxxzz_0_xxxxzzz_0, g_xxxzz_0_xxxyyyy_0, g_xxxzz_0_xxxyyyz_0, g_xxxzz_0_xxxyyzz_0, g_xxxzz_0_xxxyzzz_0, g_xxxzz_0_xxxzzzz_0, g_xxxzz_0_xxyyyyy_0, g_xxxzz_0_xxyyyyz_0, g_xxxzz_0_xxyyyzz_0, g_xxxzz_0_xxyyzzz_0, g_xxxzz_0_xxyzzzz_0, g_xxxzz_0_xxzzzzz_0, g_xxxzz_0_xyyyyyy_0, g_xxxzz_0_xyyyyyz_0, g_xxxzz_0_xyyyyzz_0, g_xxxzz_0_xyyyzzz_0, g_xxxzz_0_xyyzzzz_0, g_xxxzz_0_xyzzzzz_0, g_xxxzz_0_xzzzzzz_0, g_xxxzz_0_yyyyyyy_0, g_xxxzz_0_yyyyyyz_0, g_xxxzz_0_yyyyyzz_0, g_xxxzz_0_yyyyzzz_0, g_xxxzz_0_yyyzzzz_0, g_xxxzz_0_yyzzzzz_0, g_xxxzz_0_yzzzzzz_0, g_xxxzz_0_zzzzzzz_0, g_xxzz_0_xxxxxxz_1, g_xxzz_0_xxxxxyz_1, g_xxzz_0_xxxxxz_1, g_xxzz_0_xxxxxzz_1, g_xxzz_0_xxxxyyz_1, g_xxzz_0_xxxxyz_1, g_xxzz_0_xxxxyzz_1, g_xxzz_0_xxxxzz_1, g_xxzz_0_xxxxzzz_1, g_xxzz_0_xxxyyyz_1, g_xxzz_0_xxxyyz_1, g_xxzz_0_xxxyyzz_1, g_xxzz_0_xxxyzz_1, g_xxzz_0_xxxyzzz_1, g_xxzz_0_xxxzzz_1, g_xxzz_0_xxxzzzz_1, g_xxzz_0_xxyyyyz_1, g_xxzz_0_xxyyyz_1, g_xxzz_0_xxyyyzz_1, g_xxzz_0_xxyyzz_1, g_xxzz_0_xxyyzzz_1, g_xxzz_0_xxyzzz_1, g_xxzz_0_xxyzzzz_1, g_xxzz_0_xxzzzz_1, g_xxzz_0_xxzzzzz_1, g_xxzz_0_xyyyyyz_1, g_xxzz_0_xyyyyz_1, g_xxzz_0_xyyyyzz_1, g_xxzz_0_xyyyzz_1, g_xxzz_0_xyyyzzz_1, g_xxzz_0_xyyzzz_1, g_xxzz_0_xyyzzzz_1, g_xxzz_0_xyzzzz_1, g_xxzz_0_xyzzzzz_1, g_xxzz_0_xzzzzz_1, g_xxzz_0_xzzzzzz_1, g_xxzz_0_yyyyyyy_1, g_xxzz_0_yyyyyyz_1, g_xxzz_0_yyyyyz_1, g_xxzz_0_yyyyyzz_1, g_xxzz_0_yyyyzz_1, g_xxzz_0_yyyyzzz_1, g_xxzz_0_yyyzzz_1, g_xxzz_0_yyyzzzz_1, g_xxzz_0_yyzzzz_1, g_xxzz_0_yyzzzzz_1, g_xxzz_0_yzzzzz_1, g_xxzz_0_yzzzzzz_1, g_xxzz_0_zzzzzz_1, g_xxzz_0_zzzzzzz_1, g_xzz_0_xxxxxxz_0, g_xzz_0_xxxxxxz_1, g_xzz_0_xxxxxyz_0, g_xzz_0_xxxxxyz_1, g_xzz_0_xxxxxzz_0, g_xzz_0_xxxxxzz_1, g_xzz_0_xxxxyyz_0, g_xzz_0_xxxxyyz_1, g_xzz_0_xxxxyzz_0, g_xzz_0_xxxxyzz_1, g_xzz_0_xxxxzzz_0, g_xzz_0_xxxxzzz_1, g_xzz_0_xxxyyyz_0, g_xzz_0_xxxyyyz_1, g_xzz_0_xxxyyzz_0, g_xzz_0_xxxyyzz_1, g_xzz_0_xxxyzzz_0, g_xzz_0_xxxyzzz_1, g_xzz_0_xxxzzzz_0, g_xzz_0_xxxzzzz_1, g_xzz_0_xxyyyyz_0, g_xzz_0_xxyyyyz_1, g_xzz_0_xxyyyzz_0, g_xzz_0_xxyyyzz_1, g_xzz_0_xxyyzzz_0, g_xzz_0_xxyyzzz_1, g_xzz_0_xxyzzzz_0, g_xzz_0_xxyzzzz_1, g_xzz_0_xxzzzzz_0, g_xzz_0_xxzzzzz_1, g_xzz_0_xyyyyyz_0, g_xzz_0_xyyyyyz_1, g_xzz_0_xyyyyzz_0, g_xzz_0_xyyyyzz_1, g_xzz_0_xyyyzzz_0, g_xzz_0_xyyyzzz_1, g_xzz_0_xyyzzzz_0, g_xzz_0_xyyzzzz_1, g_xzz_0_xyzzzzz_0, g_xzz_0_xyzzzzz_1, g_xzz_0_xzzzzzz_0, g_xzz_0_xzzzzzz_1, g_xzz_0_yyyyyyy_0, g_xzz_0_yyyyyyy_1, g_xzz_0_yyyyyyz_0, g_xzz_0_yyyyyyz_1, g_xzz_0_yyyyyzz_0, g_xzz_0_yyyyyzz_1, g_xzz_0_yyyyzzz_0, g_xzz_0_yyyyzzz_1, g_xzz_0_yyyzzzz_0, g_xzz_0_yyyzzzz_1, g_xzz_0_yyzzzzz_0, g_xzz_0_yyzzzzz_1, g_xzz_0_yzzzzzz_0, g_xzz_0_yzzzzzz_1, g_xzz_0_zzzzzzz_0, g_xzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzz_0_xxxxxxx_0[i] = g_xxx_0_xxxxxxx_0[i] * fbe_0 - g_xxx_0_xxxxxxx_1[i] * fz_be_0 + g_xxxz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxxy_0[i] = g_xxx_0_xxxxxxy_0[i] * fbe_0 - g_xxx_0_xxxxxxy_1[i] * fz_be_0 + g_xxxz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxxz_0[i] = 2.0 * g_xzz_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxxyy_0[i] = g_xxx_0_xxxxxyy_0[i] * fbe_0 - g_xxx_0_xxxxxyy_1[i] * fz_be_0 + g_xxxz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxyz_0[i] = 2.0 * g_xzz_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxxzz_0[i] = 2.0 * g_xzz_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxyyy_0[i] = g_xxx_0_xxxxyyy_0[i] * fbe_0 - g_xxx_0_xxxxyyy_1[i] * fz_be_0 + g_xxxz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxyyz_0[i] = 2.0 * g_xzz_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxyzz_0[i] = 2.0 * g_xzz_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxzzz_0[i] = 2.0 * g_xzz_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyyyy_0[i] = g_xxx_0_xxxyyyy_0[i] * fbe_0 - g_xxx_0_xxxyyyy_1[i] * fz_be_0 + g_xxxz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxyyyz_0[i] = 2.0 * g_xzz_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyyzz_0[i] = 2.0 * g_xzz_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyzzz_0[i] = 2.0 * g_xzz_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxzzzz_0[i] = 2.0 * g_xzz_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyyyy_0[i] = g_xxx_0_xxyyyyy_0[i] * fbe_0 - g_xxx_0_xxyyyyy_1[i] * fz_be_0 + g_xxxz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxyyyyz_0[i] = 2.0 * g_xzz_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyyzz_0[i] = 2.0 * g_xzz_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyzzz_0[i] = 2.0 * g_xzz_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyzzzz_0[i] = 2.0 * g_xzz_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxzzzzz_0[i] = 2.0 * g_xzz_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyyyy_0[i] = g_xxx_0_xyyyyyy_0[i] * fbe_0 - g_xxx_0_xyyyyyy_1[i] * fz_be_0 + g_xxxz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xyyyyyz_0[i] = 2.0 * g_xzz_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyyyz_1[i] * fz_be_0 + g_xxzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyyzz_0[i] = 2.0 * g_xzz_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyyzz_1[i] * fz_be_0 + g_xxzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyzzz_0[i] = 2.0 * g_xzz_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyzzz_1[i] * fz_be_0 + g_xxzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyzzzz_0[i] = 2.0 * g_xzz_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyzzzz_1[i] * fz_be_0 + g_xxzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyzzzzz_0[i] = 2.0 * g_xzz_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyzzzzz_1[i] * fz_be_0 + g_xxzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xzzzzzz_0[i] = 2.0 * g_xzz_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyyy_0[i] = 2.0 * g_xzz_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyyy_1[i] * fz_be_0 + g_xxzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyyz_0[i] = 2.0 * g_xzz_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyyz_1[i] * fz_be_0 + g_xxzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyzz_0[i] = 2.0 * g_xzz_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyzz_1[i] * fz_be_0 + g_xxzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyzzz_0[i] = 2.0 * g_xzz_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyzzz_1[i] * fz_be_0 + g_xxzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyzzzz_0[i] = 2.0 * g_xzz_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyzzzz_1[i] * fz_be_0 + g_xxzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyzzzzz_0[i] = 2.0 * g_xzz_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyzzzzz_1[i] * fz_be_0 + g_xxzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yzzzzzz_0[i] = 2.0 * g_xzz_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yzzzzzz_1[i] * fz_be_0 + g_xxzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_zzzzzzz_0[i] = 2.0 * g_xzz_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_zzzzzzz_1[i] * fz_be_0 + g_xxzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 216-252 components of targeted buffer : HSK

    auto g_xxyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 216);

    auto g_xxyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 217);

    auto g_xxyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 218);

    auto g_xxyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 219);

    auto g_xxyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 220);

    auto g_xxyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 221);

    auto g_xxyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 222);

    auto g_xxyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 223);

    auto g_xxyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 224);

    auto g_xxyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 225);

    auto g_xxyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 226);

    auto g_xxyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 227);

    auto g_xxyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 228);

    auto g_xxyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 229);

    auto g_xxyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 230);

    auto g_xxyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 231);

    auto g_xxyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 232);

    auto g_xxyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 233);

    auto g_xxyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 234);

    auto g_xxyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 235);

    auto g_xxyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 236);

    auto g_xxyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 237);

    auto g_xxyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 238);

    auto g_xxyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 239);

    auto g_xxyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 240);

    auto g_xxyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 241);

    auto g_xxyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 242);

    auto g_xxyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 243);

    auto g_xxyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 244);

    auto g_xxyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 245);

    auto g_xxyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 246);

    auto g_xxyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 247);

    auto g_xxyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 248);

    auto g_xxyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 249);

    auto g_xxyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 250);

    auto g_xxyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 251);

    #pragma omp simd aligned(g_xxy_0_xxxxxxx_0, g_xxy_0_xxxxxxx_1, g_xxy_0_xxxxxxz_0, g_xxy_0_xxxxxxz_1, g_xxy_0_xxxxxzz_0, g_xxy_0_xxxxxzz_1, g_xxy_0_xxxxzzz_0, g_xxy_0_xxxxzzz_1, g_xxy_0_xxxzzzz_0, g_xxy_0_xxxzzzz_1, g_xxy_0_xxzzzzz_0, g_xxy_0_xxzzzzz_1, g_xxy_0_xzzzzzz_0, g_xxy_0_xzzzzzz_1, g_xxyy_0_xxxxxxx_1, g_xxyy_0_xxxxxxz_1, g_xxyy_0_xxxxxzz_1, g_xxyy_0_xxxxzzz_1, g_xxyy_0_xxxzzzz_1, g_xxyy_0_xxzzzzz_1, g_xxyy_0_xzzzzzz_1, g_xxyyy_0_xxxxxxx_0, g_xxyyy_0_xxxxxxy_0, g_xxyyy_0_xxxxxxz_0, g_xxyyy_0_xxxxxyy_0, g_xxyyy_0_xxxxxyz_0, g_xxyyy_0_xxxxxzz_0, g_xxyyy_0_xxxxyyy_0, g_xxyyy_0_xxxxyyz_0, g_xxyyy_0_xxxxyzz_0, g_xxyyy_0_xxxxzzz_0, g_xxyyy_0_xxxyyyy_0, g_xxyyy_0_xxxyyyz_0, g_xxyyy_0_xxxyyzz_0, g_xxyyy_0_xxxyzzz_0, g_xxyyy_0_xxxzzzz_0, g_xxyyy_0_xxyyyyy_0, g_xxyyy_0_xxyyyyz_0, g_xxyyy_0_xxyyyzz_0, g_xxyyy_0_xxyyzzz_0, g_xxyyy_0_xxyzzzz_0, g_xxyyy_0_xxzzzzz_0, g_xxyyy_0_xyyyyyy_0, g_xxyyy_0_xyyyyyz_0, g_xxyyy_0_xyyyyzz_0, g_xxyyy_0_xyyyzzz_0, g_xxyyy_0_xyyzzzz_0, g_xxyyy_0_xyzzzzz_0, g_xxyyy_0_xzzzzzz_0, g_xxyyy_0_yyyyyyy_0, g_xxyyy_0_yyyyyyz_0, g_xxyyy_0_yyyyyzz_0, g_xxyyy_0_yyyyzzz_0, g_xxyyy_0_yyyzzzz_0, g_xxyyy_0_yyzzzzz_0, g_xxyyy_0_yzzzzzz_0, g_xxyyy_0_zzzzzzz_0, g_xyyy_0_xxxxxxy_1, g_xyyy_0_xxxxxy_1, g_xyyy_0_xxxxxyy_1, g_xyyy_0_xxxxxyz_1, g_xyyy_0_xxxxyy_1, g_xyyy_0_xxxxyyy_1, g_xyyy_0_xxxxyyz_1, g_xyyy_0_xxxxyz_1, g_xyyy_0_xxxxyzz_1, g_xyyy_0_xxxyyy_1, g_xyyy_0_xxxyyyy_1, g_xyyy_0_xxxyyyz_1, g_xyyy_0_xxxyyz_1, g_xyyy_0_xxxyyzz_1, g_xyyy_0_xxxyzz_1, g_xyyy_0_xxxyzzz_1, g_xyyy_0_xxyyyy_1, g_xyyy_0_xxyyyyy_1, g_xyyy_0_xxyyyyz_1, g_xyyy_0_xxyyyz_1, g_xyyy_0_xxyyyzz_1, g_xyyy_0_xxyyzz_1, g_xyyy_0_xxyyzzz_1, g_xyyy_0_xxyzzz_1, g_xyyy_0_xxyzzzz_1, g_xyyy_0_xyyyyy_1, g_xyyy_0_xyyyyyy_1, g_xyyy_0_xyyyyyz_1, g_xyyy_0_xyyyyz_1, g_xyyy_0_xyyyyzz_1, g_xyyy_0_xyyyzz_1, g_xyyy_0_xyyyzzz_1, g_xyyy_0_xyyzzz_1, g_xyyy_0_xyyzzzz_1, g_xyyy_0_xyzzzz_1, g_xyyy_0_xyzzzzz_1, g_xyyy_0_yyyyyy_1, g_xyyy_0_yyyyyyy_1, g_xyyy_0_yyyyyyz_1, g_xyyy_0_yyyyyz_1, g_xyyy_0_yyyyyzz_1, g_xyyy_0_yyyyzz_1, g_xyyy_0_yyyyzzz_1, g_xyyy_0_yyyzzz_1, g_xyyy_0_yyyzzzz_1, g_xyyy_0_yyzzzz_1, g_xyyy_0_yyzzzzz_1, g_xyyy_0_yzzzzz_1, g_xyyy_0_yzzzzzz_1, g_xyyy_0_zzzzzzz_1, g_yyy_0_xxxxxxy_0, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxyy_0, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxxyz_0, g_yyy_0_xxxxxyz_1, g_yyy_0_xxxxyyy_0, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxxyyz_0, g_yyy_0_xxxxyyz_1, g_yyy_0_xxxxyzz_0, g_yyy_0_xxxxyzz_1, g_yyy_0_xxxyyyy_0, g_yyy_0_xxxyyyy_1, g_yyy_0_xxxyyyz_0, g_yyy_0_xxxyyyz_1, g_yyy_0_xxxyyzz_0, g_yyy_0_xxxyyzz_1, g_yyy_0_xxxyzzz_0, g_yyy_0_xxxyzzz_1, g_yyy_0_xxyyyyy_0, g_yyy_0_xxyyyyy_1, g_yyy_0_xxyyyyz_0, g_yyy_0_xxyyyyz_1, g_yyy_0_xxyyyzz_0, g_yyy_0_xxyyyzz_1, g_yyy_0_xxyyzzz_0, g_yyy_0_xxyyzzz_1, g_yyy_0_xxyzzzz_0, g_yyy_0_xxyzzzz_1, g_yyy_0_xyyyyyy_0, g_yyy_0_xyyyyyy_1, g_yyy_0_xyyyyyz_0, g_yyy_0_xyyyyyz_1, g_yyy_0_xyyyyzz_0, g_yyy_0_xyyyyzz_1, g_yyy_0_xyyyzzz_0, g_yyy_0_xyyyzzz_1, g_yyy_0_xyyzzzz_0, g_yyy_0_xyyzzzz_1, g_yyy_0_xyzzzzz_0, g_yyy_0_xyzzzzz_1, g_yyy_0_yyyyyyy_0, g_yyy_0_yyyyyyy_1, g_yyy_0_yyyyyyz_0, g_yyy_0_yyyyyyz_1, g_yyy_0_yyyyyzz_0, g_yyy_0_yyyyyzz_1, g_yyy_0_yyyyzzz_0, g_yyy_0_yyyyzzz_1, g_yyy_0_yyyzzzz_0, g_yyy_0_yyyzzzz_1, g_yyy_0_yyzzzzz_0, g_yyy_0_yyzzzzz_1, g_yyy_0_yzzzzzz_0, g_yyy_0_yzzzzzz_1, g_yyy_0_zzzzzzz_0, g_yyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyy_0_xxxxxxx_0[i] = 2.0 * g_xxy_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxxx_1[i] * fz_be_0 + g_xxyy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyyy_0_xxxxxxy_0[i] = g_yyy_0_xxxxxxy_0[i] * fbe_0 - g_yyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxxz_0[i] = 2.0 * g_xxy_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxxz_1[i] * fz_be_0 + g_xxyy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyyy_0_xxxxxyy_0[i] = g_yyy_0_xxxxxyy_0[i] * fbe_0 - g_yyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxyz_0[i] = g_yyy_0_xxxxxyz_0[i] * fbe_0 - g_yyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxzz_0[i] = 2.0 * g_xxy_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxzz_1[i] * fz_be_0 + g_xxyy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyyy_0_xxxxyyy_0[i] = g_yyy_0_xxxxyyy_0[i] * fbe_0 - g_yyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxyyz_0[i] = g_yyy_0_xxxxyyz_0[i] * fbe_0 - g_yyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxyzz_0[i] = g_yyy_0_xxxxyzz_0[i] * fbe_0 - g_yyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxzzz_0[i] = 2.0 * g_xxy_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxzzz_1[i] * fz_be_0 + g_xxyy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyyy_0_xxxyyyy_0[i] = g_yyy_0_xxxyyyy_0[i] * fbe_0 - g_yyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxyyyz_0[i] = g_yyy_0_xxxyyyz_0[i] * fbe_0 - g_yyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxyyzz_0[i] = g_yyy_0_xxxyyzz_0[i] * fbe_0 - g_yyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxyzzz_0[i] = g_yyy_0_xxxyzzz_0[i] * fbe_0 - g_yyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxzzzz_0[i] = 2.0 * g_xxy_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxzzzz_1[i] * fz_be_0 + g_xxyy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyyy_0_xxyyyyy_0[i] = g_yyy_0_xxyyyyy_0[i] * fbe_0 - g_yyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxyyyyz_0[i] = g_yyy_0_xxyyyyz_0[i] * fbe_0 - g_yyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxyyyzz_0[i] = g_yyy_0_xxyyyzz_0[i] * fbe_0 - g_yyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxyyzzz_0[i] = g_yyy_0_xxyyzzz_0[i] * fbe_0 - g_yyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxyzzzz_0[i] = g_yyy_0_xxyzzzz_0[i] * fbe_0 - g_yyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxzzzzz_0[i] = 2.0 * g_xxy_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxzzzzz_1[i] * fz_be_0 + g_xxyy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyyy_0_xyyyyyy_0[i] = g_yyy_0_xyyyyyy_0[i] * fbe_0 - g_yyy_0_xyyyyyy_1[i] * fz_be_0 + g_xyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xyyyyyz_0[i] = g_yyy_0_xyyyyyz_0[i] * fbe_0 - g_yyy_0_xyyyyyz_1[i] * fz_be_0 + g_xyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xyyyyzz_0[i] = g_yyy_0_xyyyyzz_0[i] * fbe_0 - g_yyy_0_xyyyyzz_1[i] * fz_be_0 + g_xyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xyyyzzz_0[i] = g_yyy_0_xyyyzzz_0[i] * fbe_0 - g_yyy_0_xyyyzzz_1[i] * fz_be_0 + g_xyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xyyzzzz_0[i] = g_yyy_0_xyyzzzz_0[i] * fbe_0 - g_yyy_0_xyyzzzz_1[i] * fz_be_0 + g_xyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xyzzzzz_0[i] = g_yyy_0_xyzzzzz_0[i] * fbe_0 - g_yyy_0_xyzzzzz_1[i] * fz_be_0 + g_xyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xzzzzzz_0[i] = 2.0 * g_xxy_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xzzzzzz_1[i] * fz_be_0 + g_xxyy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyyy_0_yyyyyyy_0[i] = g_yyy_0_yyyyyyy_0[i] * fbe_0 - g_yyy_0_yyyyyyy_1[i] * fz_be_0 + g_xyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_yyyyyyz_0[i] = g_yyy_0_yyyyyyz_0[i] * fbe_0 - g_yyy_0_yyyyyyz_1[i] * fz_be_0 + g_xyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_yyyyyzz_0[i] = g_yyy_0_yyyyyzz_0[i] * fbe_0 - g_yyy_0_yyyyyzz_1[i] * fz_be_0 + g_xyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_yyyyzzz_0[i] = g_yyy_0_yyyyzzz_0[i] * fbe_0 - g_yyy_0_yyyyzzz_1[i] * fz_be_0 + g_xyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_yyyzzzz_0[i] = g_yyy_0_yyyzzzz_0[i] * fbe_0 - g_yyy_0_yyyzzzz_1[i] * fz_be_0 + g_xyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_yyzzzzz_0[i] = g_yyy_0_yyzzzzz_0[i] * fbe_0 - g_yyy_0_yyzzzzz_1[i] * fz_be_0 + g_xyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_yzzzzzz_0[i] = g_yyy_0_yzzzzzz_0[i] * fbe_0 - g_yyy_0_yzzzzzz_1[i] * fz_be_0 + g_xyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_zzzzzzz_0[i] = g_yyy_0_zzzzzzz_0[i] * fbe_0 - g_yyy_0_zzzzzzz_1[i] * fz_be_0 + g_xyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 252-288 components of targeted buffer : HSK

    auto g_xxyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 252);

    auto g_xxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 253);

    auto g_xxyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 254);

    auto g_xxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 255);

    auto g_xxyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 256);

    auto g_xxyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 257);

    auto g_xxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 258);

    auto g_xxyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 259);

    auto g_xxyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 260);

    auto g_xxyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 261);

    auto g_xxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 262);

    auto g_xxyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 263);

    auto g_xxyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 264);

    auto g_xxyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 265);

    auto g_xxyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 266);

    auto g_xxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 267);

    auto g_xxyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 268);

    auto g_xxyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 269);

    auto g_xxyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 270);

    auto g_xxyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 271);

    auto g_xxyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 272);

    auto g_xxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 273);

    auto g_xxyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 274);

    auto g_xxyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 275);

    auto g_xxyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 276);

    auto g_xxyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 277);

    auto g_xxyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 278);

    auto g_xxyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 279);

    auto g_xxyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 280);

    auto g_xxyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 281);

    auto g_xxyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 282);

    auto g_xxyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 283);

    auto g_xxyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 284);

    auto g_xxyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 285);

    auto g_xxyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 286);

    auto g_xxyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 287);

    #pragma omp simd aligned(g_xxyy_0_xxxxxx_1, g_xxyy_0_xxxxxxx_1, g_xxyy_0_xxxxxxy_1, g_xxyy_0_xxxxxxz_1, g_xxyy_0_xxxxxy_1, g_xxyy_0_xxxxxyy_1, g_xxyy_0_xxxxxyz_1, g_xxyy_0_xxxxxz_1, g_xxyy_0_xxxxxzz_1, g_xxyy_0_xxxxyy_1, g_xxyy_0_xxxxyyy_1, g_xxyy_0_xxxxyyz_1, g_xxyy_0_xxxxyz_1, g_xxyy_0_xxxxyzz_1, g_xxyy_0_xxxxzz_1, g_xxyy_0_xxxxzzz_1, g_xxyy_0_xxxyyy_1, g_xxyy_0_xxxyyyy_1, g_xxyy_0_xxxyyyz_1, g_xxyy_0_xxxyyz_1, g_xxyy_0_xxxyyzz_1, g_xxyy_0_xxxyzz_1, g_xxyy_0_xxxyzzz_1, g_xxyy_0_xxxzzz_1, g_xxyy_0_xxxzzzz_1, g_xxyy_0_xxyyyy_1, g_xxyy_0_xxyyyyy_1, g_xxyy_0_xxyyyyz_1, g_xxyy_0_xxyyyz_1, g_xxyy_0_xxyyyzz_1, g_xxyy_0_xxyyzz_1, g_xxyy_0_xxyyzzz_1, g_xxyy_0_xxyzzz_1, g_xxyy_0_xxyzzzz_1, g_xxyy_0_xxzzzz_1, g_xxyy_0_xxzzzzz_1, g_xxyy_0_xyyyyy_1, g_xxyy_0_xyyyyyy_1, g_xxyy_0_xyyyyyz_1, g_xxyy_0_xyyyyz_1, g_xxyy_0_xyyyyzz_1, g_xxyy_0_xyyyzz_1, g_xxyy_0_xyyyzzz_1, g_xxyy_0_xyyzzz_1, g_xxyy_0_xyyzzzz_1, g_xxyy_0_xyzzzz_1, g_xxyy_0_xyzzzzz_1, g_xxyy_0_xzzzzz_1, g_xxyy_0_xzzzzzz_1, g_xxyy_0_yyyyyy_1, g_xxyy_0_yyyyyyy_1, g_xxyy_0_yyyyyyz_1, g_xxyy_0_yyyyyz_1, g_xxyy_0_yyyyyzz_1, g_xxyy_0_yyyyzz_1, g_xxyy_0_yyyyzzz_1, g_xxyy_0_yyyzzz_1, g_xxyy_0_yyyzzzz_1, g_xxyy_0_yyzzzz_1, g_xxyy_0_yyzzzzz_1, g_xxyy_0_yzzzzz_1, g_xxyy_0_yzzzzzz_1, g_xxyy_0_zzzzzz_1, g_xxyy_0_zzzzzzz_1, g_xxyyz_0_xxxxxxx_0, g_xxyyz_0_xxxxxxy_0, g_xxyyz_0_xxxxxxz_0, g_xxyyz_0_xxxxxyy_0, g_xxyyz_0_xxxxxyz_0, g_xxyyz_0_xxxxxzz_0, g_xxyyz_0_xxxxyyy_0, g_xxyyz_0_xxxxyyz_0, g_xxyyz_0_xxxxyzz_0, g_xxyyz_0_xxxxzzz_0, g_xxyyz_0_xxxyyyy_0, g_xxyyz_0_xxxyyyz_0, g_xxyyz_0_xxxyyzz_0, g_xxyyz_0_xxxyzzz_0, g_xxyyz_0_xxxzzzz_0, g_xxyyz_0_xxyyyyy_0, g_xxyyz_0_xxyyyyz_0, g_xxyyz_0_xxyyyzz_0, g_xxyyz_0_xxyyzzz_0, g_xxyyz_0_xxyzzzz_0, g_xxyyz_0_xxzzzzz_0, g_xxyyz_0_xyyyyyy_0, g_xxyyz_0_xyyyyyz_0, g_xxyyz_0_xyyyyzz_0, g_xxyyz_0_xyyyzzz_0, g_xxyyz_0_xyyzzzz_0, g_xxyyz_0_xyzzzzz_0, g_xxyyz_0_xzzzzzz_0, g_xxyyz_0_yyyyyyy_0, g_xxyyz_0_yyyyyyz_0, g_xxyyz_0_yyyyyzz_0, g_xxyyz_0_yyyyzzz_0, g_xxyyz_0_yyyzzzz_0, g_xxyyz_0_yyzzzzz_0, g_xxyyz_0_yzzzzzz_0, g_xxyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyz_0_xxxxxxx_0[i] = g_xxyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxxy_0[i] = g_xxyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxxz_0[i] = g_xxyy_0_xxxxxx_1[i] * fi_acd_0 + g_xxyy_0_xxxxxxz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxyy_0[i] = g_xxyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxyz_0[i] = g_xxyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxzz_0[i] = 2.0 * g_xxyy_0_xxxxxz_1[i] * fi_acd_0 + g_xxyy_0_xxxxxzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyyy_0[i] = g_xxyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyyz_0[i] = g_xxyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyzz_0[i] = 2.0 * g_xxyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxzzz_0[i] = 3.0 * g_xxyy_0_xxxxzz_1[i] * fi_acd_0 + g_xxyy_0_xxxxzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyyy_0[i] = g_xxyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyyz_0[i] = g_xxyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyzz_0[i] = 2.0 * g_xxyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyzzz_0[i] = 3.0 * g_xxyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxzzzz_0[i] = 4.0 * g_xxyy_0_xxxzzz_1[i] * fi_acd_0 + g_xxyy_0_xxxzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyyy_0[i] = g_xxyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyyz_0[i] = g_xxyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyzz_0[i] = 2.0 * g_xxyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyzzz_0[i] = 3.0 * g_xxyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyzzzz_0[i] = 4.0 * g_xxyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxzzzzz_0[i] = 5.0 * g_xxyy_0_xxzzzz_1[i] * fi_acd_0 + g_xxyy_0_xxzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyyy_0[i] = g_xxyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyyz_0[i] = g_xxyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxyy_0_xyyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyzz_0[i] = 2.0 * g_xxyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyzzz_0[i] = 3.0 * g_xxyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyzzzz_0[i] = 4.0 * g_xxyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyzzzzz_0[i] = 5.0 * g_xxyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xzzzzzz_0[i] = 6.0 * g_xxyy_0_xzzzzz_1[i] * fi_acd_0 + g_xxyy_0_xzzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyyy_0[i] = g_xxyy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyyz_0[i] = g_xxyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxyy_0_yyyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyzz_0[i] = 2.0 * g_xxyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxyy_0_yyyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyzzz_0[i] = 3.0 * g_xxyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxyy_0_yyyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyzzzz_0[i] = 4.0 * g_xxyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxyy_0_yyyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyzzzzz_0[i] = 5.0 * g_xxyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxyy_0_yyzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yzzzzzz_0[i] = 6.0 * g_xxyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxyy_0_yzzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_zzzzzzz_0[i] = 7.0 * g_xxyy_0_zzzzzz_1[i] * fi_acd_0 + g_xxyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 288-324 components of targeted buffer : HSK

    auto g_xxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 288);

    auto g_xxyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 289);

    auto g_xxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 290);

    auto g_xxyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 291);

    auto g_xxyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 292);

    auto g_xxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 293);

    auto g_xxyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 294);

    auto g_xxyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 295);

    auto g_xxyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 296);

    auto g_xxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 297);

    auto g_xxyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 298);

    auto g_xxyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 299);

    auto g_xxyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 300);

    auto g_xxyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 301);

    auto g_xxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 302);

    auto g_xxyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 303);

    auto g_xxyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 304);

    auto g_xxyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 305);

    auto g_xxyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 306);

    auto g_xxyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 307);

    auto g_xxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 308);

    auto g_xxyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 309);

    auto g_xxyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 310);

    auto g_xxyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 311);

    auto g_xxyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 312);

    auto g_xxyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 313);

    auto g_xxyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 314);

    auto g_xxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 315);

    auto g_xxyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 316);

    auto g_xxyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 317);

    auto g_xxyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 318);

    auto g_xxyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 319);

    auto g_xxyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 320);

    auto g_xxyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 321);

    auto g_xxyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 322);

    auto g_xxyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 323);

    #pragma omp simd aligned(g_xxyzz_0_xxxxxxx_0, g_xxyzz_0_xxxxxxy_0, g_xxyzz_0_xxxxxxz_0, g_xxyzz_0_xxxxxyy_0, g_xxyzz_0_xxxxxyz_0, g_xxyzz_0_xxxxxzz_0, g_xxyzz_0_xxxxyyy_0, g_xxyzz_0_xxxxyyz_0, g_xxyzz_0_xxxxyzz_0, g_xxyzz_0_xxxxzzz_0, g_xxyzz_0_xxxyyyy_0, g_xxyzz_0_xxxyyyz_0, g_xxyzz_0_xxxyyzz_0, g_xxyzz_0_xxxyzzz_0, g_xxyzz_0_xxxzzzz_0, g_xxyzz_0_xxyyyyy_0, g_xxyzz_0_xxyyyyz_0, g_xxyzz_0_xxyyyzz_0, g_xxyzz_0_xxyyzzz_0, g_xxyzz_0_xxyzzzz_0, g_xxyzz_0_xxzzzzz_0, g_xxyzz_0_xyyyyyy_0, g_xxyzz_0_xyyyyyz_0, g_xxyzz_0_xyyyyzz_0, g_xxyzz_0_xyyyzzz_0, g_xxyzz_0_xyyzzzz_0, g_xxyzz_0_xyzzzzz_0, g_xxyzz_0_xzzzzzz_0, g_xxyzz_0_yyyyyyy_0, g_xxyzz_0_yyyyyyz_0, g_xxyzz_0_yyyyyzz_0, g_xxyzz_0_yyyyzzz_0, g_xxyzz_0_yyyzzzz_0, g_xxyzz_0_yyzzzzz_0, g_xxyzz_0_yzzzzzz_0, g_xxyzz_0_zzzzzzz_0, g_xxzz_0_xxxxxx_1, g_xxzz_0_xxxxxxx_1, g_xxzz_0_xxxxxxy_1, g_xxzz_0_xxxxxxz_1, g_xxzz_0_xxxxxy_1, g_xxzz_0_xxxxxyy_1, g_xxzz_0_xxxxxyz_1, g_xxzz_0_xxxxxz_1, g_xxzz_0_xxxxxzz_1, g_xxzz_0_xxxxyy_1, g_xxzz_0_xxxxyyy_1, g_xxzz_0_xxxxyyz_1, g_xxzz_0_xxxxyz_1, g_xxzz_0_xxxxyzz_1, g_xxzz_0_xxxxzz_1, g_xxzz_0_xxxxzzz_1, g_xxzz_0_xxxyyy_1, g_xxzz_0_xxxyyyy_1, g_xxzz_0_xxxyyyz_1, g_xxzz_0_xxxyyz_1, g_xxzz_0_xxxyyzz_1, g_xxzz_0_xxxyzz_1, g_xxzz_0_xxxyzzz_1, g_xxzz_0_xxxzzz_1, g_xxzz_0_xxxzzzz_1, g_xxzz_0_xxyyyy_1, g_xxzz_0_xxyyyyy_1, g_xxzz_0_xxyyyyz_1, g_xxzz_0_xxyyyz_1, g_xxzz_0_xxyyyzz_1, g_xxzz_0_xxyyzz_1, g_xxzz_0_xxyyzzz_1, g_xxzz_0_xxyzzz_1, g_xxzz_0_xxyzzzz_1, g_xxzz_0_xxzzzz_1, g_xxzz_0_xxzzzzz_1, g_xxzz_0_xyyyyy_1, g_xxzz_0_xyyyyyy_1, g_xxzz_0_xyyyyyz_1, g_xxzz_0_xyyyyz_1, g_xxzz_0_xyyyyzz_1, g_xxzz_0_xyyyzz_1, g_xxzz_0_xyyyzzz_1, g_xxzz_0_xyyzzz_1, g_xxzz_0_xyyzzzz_1, g_xxzz_0_xyzzzz_1, g_xxzz_0_xyzzzzz_1, g_xxzz_0_xzzzzz_1, g_xxzz_0_xzzzzzz_1, g_xxzz_0_yyyyyy_1, g_xxzz_0_yyyyyyy_1, g_xxzz_0_yyyyyyz_1, g_xxzz_0_yyyyyz_1, g_xxzz_0_yyyyyzz_1, g_xxzz_0_yyyyzz_1, g_xxzz_0_yyyyzzz_1, g_xxzz_0_yyyzzz_1, g_xxzz_0_yyyzzzz_1, g_xxzz_0_yyzzzz_1, g_xxzz_0_yyzzzzz_1, g_xxzz_0_yzzzzz_1, g_xxzz_0_yzzzzzz_1, g_xxzz_0_zzzzzz_1, g_xxzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzz_0_xxxxxxx_0[i] = g_xxzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxxy_0[i] = g_xxzz_0_xxxxxx_1[i] * fi_acd_0 + g_xxzz_0_xxxxxxy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxxz_0[i] = g_xxzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxyy_0[i] = 2.0 * g_xxzz_0_xxxxxy_1[i] * fi_acd_0 + g_xxzz_0_xxxxxyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxyz_0[i] = g_xxzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxzz_0[i] = g_xxzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyyy_0[i] = 3.0 * g_xxzz_0_xxxxyy_1[i] * fi_acd_0 + g_xxzz_0_xxxxyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyyz_0[i] = 2.0 * g_xxzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyzz_0[i] = g_xxzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxzzz_0[i] = g_xxzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyyy_0[i] = 4.0 * g_xxzz_0_xxxyyy_1[i] * fi_acd_0 + g_xxzz_0_xxxyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyyz_0[i] = 3.0 * g_xxzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyzz_0[i] = 2.0 * g_xxzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyzzz_0[i] = g_xxzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxzzzz_0[i] = g_xxzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyyy_0[i] = 5.0 * g_xxzz_0_xxyyyy_1[i] * fi_acd_0 + g_xxzz_0_xxyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyyz_0[i] = 4.0 * g_xxzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyzz_0[i] = 3.0 * g_xxzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyzzz_0[i] = 2.0 * g_xxzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyzzzz_0[i] = g_xxzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxzzzzz_0[i] = g_xxzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyyy_0[i] = 6.0 * g_xxzz_0_xyyyyy_1[i] * fi_acd_0 + g_xxzz_0_xyyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyyz_0[i] = 5.0 * g_xxzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyzz_0[i] = 4.0 * g_xxzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyzzz_0[i] = 3.0 * g_xxzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyzzzz_0[i] = 2.0 * g_xxzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyzzzzz_0[i] = g_xxzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xzzzzzz_0[i] = g_xxzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyyy_0[i] = 7.0 * g_xxzz_0_yyyyyy_1[i] * fi_acd_0 + g_xxzz_0_yyyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyyz_0[i] = 6.0 * g_xxzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxzz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyzz_0[i] = 5.0 * g_xxzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxzz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyzzz_0[i] = 4.0 * g_xxzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxzz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyzzzz_0[i] = 3.0 * g_xxzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxzz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyzzzzz_0[i] = 2.0 * g_xxzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxzz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yzzzzzz_0[i] = g_xxzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxzz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_zzzzzzz_0[i] = g_xxzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 324-360 components of targeted buffer : HSK

    auto g_xxzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 324);

    auto g_xxzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 325);

    auto g_xxzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 326);

    auto g_xxzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 327);

    auto g_xxzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 328);

    auto g_xxzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 329);

    auto g_xxzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 330);

    auto g_xxzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 331);

    auto g_xxzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 332);

    auto g_xxzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 333);

    auto g_xxzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 334);

    auto g_xxzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 335);

    auto g_xxzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 336);

    auto g_xxzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 337);

    auto g_xxzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 338);

    auto g_xxzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 339);

    auto g_xxzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 340);

    auto g_xxzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 341);

    auto g_xxzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 342);

    auto g_xxzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 343);

    auto g_xxzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 344);

    auto g_xxzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 345);

    auto g_xxzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 346);

    auto g_xxzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 347);

    auto g_xxzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 348);

    auto g_xxzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 349);

    auto g_xxzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 350);

    auto g_xxzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 351);

    auto g_xxzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 352);

    auto g_xxzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 353);

    auto g_xxzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 354);

    auto g_xxzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 355);

    auto g_xxzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 356);

    auto g_xxzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 357);

    auto g_xxzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 358);

    auto g_xxzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 359);

    #pragma omp simd aligned(g_xxz_0_xxxxxxx_0, g_xxz_0_xxxxxxx_1, g_xxz_0_xxxxxxy_0, g_xxz_0_xxxxxxy_1, g_xxz_0_xxxxxyy_0, g_xxz_0_xxxxxyy_1, g_xxz_0_xxxxyyy_0, g_xxz_0_xxxxyyy_1, g_xxz_0_xxxyyyy_0, g_xxz_0_xxxyyyy_1, g_xxz_0_xxyyyyy_0, g_xxz_0_xxyyyyy_1, g_xxz_0_xyyyyyy_0, g_xxz_0_xyyyyyy_1, g_xxzz_0_xxxxxxx_1, g_xxzz_0_xxxxxxy_1, g_xxzz_0_xxxxxyy_1, g_xxzz_0_xxxxyyy_1, g_xxzz_0_xxxyyyy_1, g_xxzz_0_xxyyyyy_1, g_xxzz_0_xyyyyyy_1, g_xxzzz_0_xxxxxxx_0, g_xxzzz_0_xxxxxxy_0, g_xxzzz_0_xxxxxxz_0, g_xxzzz_0_xxxxxyy_0, g_xxzzz_0_xxxxxyz_0, g_xxzzz_0_xxxxxzz_0, g_xxzzz_0_xxxxyyy_0, g_xxzzz_0_xxxxyyz_0, g_xxzzz_0_xxxxyzz_0, g_xxzzz_0_xxxxzzz_0, g_xxzzz_0_xxxyyyy_0, g_xxzzz_0_xxxyyyz_0, g_xxzzz_0_xxxyyzz_0, g_xxzzz_0_xxxyzzz_0, g_xxzzz_0_xxxzzzz_0, g_xxzzz_0_xxyyyyy_0, g_xxzzz_0_xxyyyyz_0, g_xxzzz_0_xxyyyzz_0, g_xxzzz_0_xxyyzzz_0, g_xxzzz_0_xxyzzzz_0, g_xxzzz_0_xxzzzzz_0, g_xxzzz_0_xyyyyyy_0, g_xxzzz_0_xyyyyyz_0, g_xxzzz_0_xyyyyzz_0, g_xxzzz_0_xyyyzzz_0, g_xxzzz_0_xyyzzzz_0, g_xxzzz_0_xyzzzzz_0, g_xxzzz_0_xzzzzzz_0, g_xxzzz_0_yyyyyyy_0, g_xxzzz_0_yyyyyyz_0, g_xxzzz_0_yyyyyzz_0, g_xxzzz_0_yyyyzzz_0, g_xxzzz_0_yyyzzzz_0, g_xxzzz_0_yyzzzzz_0, g_xxzzz_0_yzzzzzz_0, g_xxzzz_0_zzzzzzz_0, g_xzzz_0_xxxxxxz_1, g_xzzz_0_xxxxxyz_1, g_xzzz_0_xxxxxz_1, g_xzzz_0_xxxxxzz_1, g_xzzz_0_xxxxyyz_1, g_xzzz_0_xxxxyz_1, g_xzzz_0_xxxxyzz_1, g_xzzz_0_xxxxzz_1, g_xzzz_0_xxxxzzz_1, g_xzzz_0_xxxyyyz_1, g_xzzz_0_xxxyyz_1, g_xzzz_0_xxxyyzz_1, g_xzzz_0_xxxyzz_1, g_xzzz_0_xxxyzzz_1, g_xzzz_0_xxxzzz_1, g_xzzz_0_xxxzzzz_1, g_xzzz_0_xxyyyyz_1, g_xzzz_0_xxyyyz_1, g_xzzz_0_xxyyyzz_1, g_xzzz_0_xxyyzz_1, g_xzzz_0_xxyyzzz_1, g_xzzz_0_xxyzzz_1, g_xzzz_0_xxyzzzz_1, g_xzzz_0_xxzzzz_1, g_xzzz_0_xxzzzzz_1, g_xzzz_0_xyyyyyz_1, g_xzzz_0_xyyyyz_1, g_xzzz_0_xyyyyzz_1, g_xzzz_0_xyyyzz_1, g_xzzz_0_xyyyzzz_1, g_xzzz_0_xyyzzz_1, g_xzzz_0_xyyzzzz_1, g_xzzz_0_xyzzzz_1, g_xzzz_0_xyzzzzz_1, g_xzzz_0_xzzzzz_1, g_xzzz_0_xzzzzzz_1, g_xzzz_0_yyyyyyy_1, g_xzzz_0_yyyyyyz_1, g_xzzz_0_yyyyyz_1, g_xzzz_0_yyyyyzz_1, g_xzzz_0_yyyyzz_1, g_xzzz_0_yyyyzzz_1, g_xzzz_0_yyyzzz_1, g_xzzz_0_yyyzzzz_1, g_xzzz_0_yyzzzz_1, g_xzzz_0_yyzzzzz_1, g_xzzz_0_yzzzzz_1, g_xzzz_0_yzzzzzz_1, g_xzzz_0_zzzzzz_1, g_xzzz_0_zzzzzzz_1, g_zzz_0_xxxxxxz_0, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxyz_0, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxzz_0, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxyyz_0, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyzz_0, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxzzz_0, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxyyyz_0, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyzz_0, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyzzz_0, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxzzzz_0, g_zzz_0_xxxzzzz_1, g_zzz_0_xxyyyyz_0, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyzz_0, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyzzz_0, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyzzzz_0, g_zzz_0_xxyzzzz_1, g_zzz_0_xxzzzzz_0, g_zzz_0_xxzzzzz_1, g_zzz_0_xyyyyyz_0, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyzz_0, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyzzz_0, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyzzzz_0, g_zzz_0_xyyzzzz_1, g_zzz_0_xyzzzzz_0, g_zzz_0_xyzzzzz_1, g_zzz_0_xzzzzzz_0, g_zzz_0_xzzzzzz_1, g_zzz_0_yyyyyyy_0, g_zzz_0_yyyyyyy_1, g_zzz_0_yyyyyyz_0, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyzz_0, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyzzz_0, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyzzzz_0, g_zzz_0_yyyzzzz_1, g_zzz_0_yyzzzzz_0, g_zzz_0_yyzzzzz_1, g_zzz_0_yzzzzzz_0, g_zzz_0_yzzzzzz_1, g_zzz_0_zzzzzzz_0, g_zzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzz_0_xxxxxxx_0[i] = 2.0 * g_xxz_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxxx_1[i] * fz_be_0 + g_xxzz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxxy_0[i] = 2.0 * g_xxz_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxxy_1[i] * fz_be_0 + g_xxzz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxxz_0[i] = g_zzz_0_xxxxxxz_0[i] * fbe_0 - g_zzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxxyy_0[i] = 2.0 * g_xxz_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxyy_1[i] * fz_be_0 + g_xxzz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxyz_0[i] = g_zzz_0_xxxxxyz_0[i] * fbe_0 - g_zzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxxzz_0[i] = g_zzz_0_xxxxxzz_0[i] * fbe_0 - g_zzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxyyy_0[i] = 2.0 * g_xxz_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxyyy_1[i] * fz_be_0 + g_xxzz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxyyz_0[i] = g_zzz_0_xxxxyyz_0[i] * fbe_0 - g_zzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxyzz_0[i] = g_zzz_0_xxxxyzz_0[i] * fbe_0 - g_zzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxzzz_0[i] = g_zzz_0_xxxxzzz_0[i] * fbe_0 - g_zzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyyyy_0[i] = 2.0 * g_xxz_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxyyyy_1[i] * fz_be_0 + g_xxzz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxyyyz_0[i] = g_zzz_0_xxxyyyz_0[i] * fbe_0 - g_zzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyyzz_0[i] = g_zzz_0_xxxyyzz_0[i] * fbe_0 - g_zzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyzzz_0[i] = g_zzz_0_xxxyzzz_0[i] * fbe_0 - g_zzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxzzzz_0[i] = g_zzz_0_xxxzzzz_0[i] * fbe_0 - g_zzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyyyy_0[i] = 2.0 * g_xxz_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxyyyyy_1[i] * fz_be_0 + g_xxzz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxyyyyz_0[i] = g_zzz_0_xxyyyyz_0[i] * fbe_0 - g_zzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyyzz_0[i] = g_zzz_0_xxyyyzz_0[i] * fbe_0 - g_zzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyzzz_0[i] = g_zzz_0_xxyyzzz_0[i] * fbe_0 - g_zzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyzzzz_0[i] = g_zzz_0_xxyzzzz_0[i] * fbe_0 - g_zzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxzzzzz_0[i] = g_zzz_0_xxzzzzz_0[i] * fbe_0 - g_zzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyyyy_0[i] = 2.0 * g_xxz_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xyyyyyy_1[i] * fz_be_0 + g_xxzz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xyyyyyz_0[i] = g_zzz_0_xyyyyyz_0[i] * fbe_0 - g_zzz_0_xyyyyyz_1[i] * fz_be_0 + g_xzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyyzz_0[i] = g_zzz_0_xyyyyzz_0[i] * fbe_0 - g_zzz_0_xyyyyzz_1[i] * fz_be_0 + g_xzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyzzz_0[i] = g_zzz_0_xyyyzzz_0[i] * fbe_0 - g_zzz_0_xyyyzzz_1[i] * fz_be_0 + g_xzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyzzzz_0[i] = g_zzz_0_xyyzzzz_0[i] * fbe_0 - g_zzz_0_xyyzzzz_1[i] * fz_be_0 + g_xzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyzzzzz_0[i] = g_zzz_0_xyzzzzz_0[i] * fbe_0 - g_zzz_0_xyzzzzz_1[i] * fz_be_0 + g_xzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xzzzzzz_0[i] = g_zzz_0_xzzzzzz_0[i] * fbe_0 - g_zzz_0_xzzzzzz_1[i] * fz_be_0 + g_xzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyyy_0[i] = g_zzz_0_yyyyyyy_0[i] * fbe_0 - g_zzz_0_yyyyyyy_1[i] * fz_be_0 + g_xzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyyz_0[i] = g_zzz_0_yyyyyyz_0[i] * fbe_0 - g_zzz_0_yyyyyyz_1[i] * fz_be_0 + g_xzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyzz_0[i] = g_zzz_0_yyyyyzz_0[i] * fbe_0 - g_zzz_0_yyyyyzz_1[i] * fz_be_0 + g_xzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyzzz_0[i] = g_zzz_0_yyyyzzz_0[i] * fbe_0 - g_zzz_0_yyyyzzz_1[i] * fz_be_0 + g_xzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyzzzz_0[i] = g_zzz_0_yyyzzzz_0[i] * fbe_0 - g_zzz_0_yyyzzzz_1[i] * fz_be_0 + g_xzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyzzzzz_0[i] = g_zzz_0_yyzzzzz_0[i] * fbe_0 - g_zzz_0_yyzzzzz_1[i] * fz_be_0 + g_xzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yzzzzzz_0[i] = g_zzz_0_yzzzzzz_0[i] * fbe_0 - g_zzz_0_yzzzzzz_1[i] * fz_be_0 + g_xzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_zzzzzzz_0[i] = g_zzz_0_zzzzzzz_0[i] * fbe_0 - g_zzz_0_zzzzzzz_1[i] * fz_be_0 + g_xzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 360-396 components of targeted buffer : HSK

    auto g_xyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 360);

    auto g_xyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 361);

    auto g_xyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 362);

    auto g_xyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 363);

    auto g_xyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 364);

    auto g_xyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 365);

    auto g_xyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 366);

    auto g_xyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 367);

    auto g_xyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 368);

    auto g_xyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 369);

    auto g_xyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 370);

    auto g_xyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 371);

    auto g_xyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 372);

    auto g_xyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 373);

    auto g_xyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 374);

    auto g_xyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 375);

    auto g_xyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 376);

    auto g_xyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 377);

    auto g_xyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 378);

    auto g_xyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 379);

    auto g_xyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 380);

    auto g_xyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 381);

    auto g_xyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 382);

    auto g_xyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 383);

    auto g_xyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 384);

    auto g_xyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 385);

    auto g_xyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 386);

    auto g_xyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 387);

    auto g_xyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 388);

    auto g_xyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 389);

    auto g_xyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 390);

    auto g_xyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 391);

    auto g_xyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 392);

    auto g_xyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 393);

    auto g_xyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 394);

    auto g_xyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 395);

    #pragma omp simd aligned(g_xyyyy_0_xxxxxxx_0, g_xyyyy_0_xxxxxxy_0, g_xyyyy_0_xxxxxxz_0, g_xyyyy_0_xxxxxyy_0, g_xyyyy_0_xxxxxyz_0, g_xyyyy_0_xxxxxzz_0, g_xyyyy_0_xxxxyyy_0, g_xyyyy_0_xxxxyyz_0, g_xyyyy_0_xxxxyzz_0, g_xyyyy_0_xxxxzzz_0, g_xyyyy_0_xxxyyyy_0, g_xyyyy_0_xxxyyyz_0, g_xyyyy_0_xxxyyzz_0, g_xyyyy_0_xxxyzzz_0, g_xyyyy_0_xxxzzzz_0, g_xyyyy_0_xxyyyyy_0, g_xyyyy_0_xxyyyyz_0, g_xyyyy_0_xxyyyzz_0, g_xyyyy_0_xxyyzzz_0, g_xyyyy_0_xxyzzzz_0, g_xyyyy_0_xxzzzzz_0, g_xyyyy_0_xyyyyyy_0, g_xyyyy_0_xyyyyyz_0, g_xyyyy_0_xyyyyzz_0, g_xyyyy_0_xyyyzzz_0, g_xyyyy_0_xyyzzzz_0, g_xyyyy_0_xyzzzzz_0, g_xyyyy_0_xzzzzzz_0, g_xyyyy_0_yyyyyyy_0, g_xyyyy_0_yyyyyyz_0, g_xyyyy_0_yyyyyzz_0, g_xyyyy_0_yyyyzzz_0, g_xyyyy_0_yyyzzzz_0, g_xyyyy_0_yyzzzzz_0, g_xyyyy_0_yzzzzzz_0, g_xyyyy_0_zzzzzzz_0, g_yyyy_0_xxxxxx_1, g_yyyy_0_xxxxxxx_1, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxxz_1, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxxyz_1, g_yyyy_0_xxxxxz_1, g_yyyy_0_xxxxxzz_1, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxxyyz_1, g_yyyy_0_xxxxyz_1, g_yyyy_0_xxxxyzz_1, g_yyyy_0_xxxxzz_1, g_yyyy_0_xxxxzzz_1, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxxyyyz_1, g_yyyy_0_xxxyyz_1, g_yyyy_0_xxxyyzz_1, g_yyyy_0_xxxyzz_1, g_yyyy_0_xxxyzzz_1, g_yyyy_0_xxxzzz_1, g_yyyy_0_xxxzzzz_1, g_yyyy_0_xxyyyy_1, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xxyyyyz_1, g_yyyy_0_xxyyyz_1, g_yyyy_0_xxyyyzz_1, g_yyyy_0_xxyyzz_1, g_yyyy_0_xxyyzzz_1, g_yyyy_0_xxyzzz_1, g_yyyy_0_xxyzzzz_1, g_yyyy_0_xxzzzz_1, g_yyyy_0_xxzzzzz_1, g_yyyy_0_xyyyyy_1, g_yyyy_0_xyyyyyy_1, g_yyyy_0_xyyyyyz_1, g_yyyy_0_xyyyyz_1, g_yyyy_0_xyyyyzz_1, g_yyyy_0_xyyyzz_1, g_yyyy_0_xyyyzzz_1, g_yyyy_0_xyyzzz_1, g_yyyy_0_xyyzzzz_1, g_yyyy_0_xyzzzz_1, g_yyyy_0_xyzzzzz_1, g_yyyy_0_xzzzzz_1, g_yyyy_0_xzzzzzz_1, g_yyyy_0_yyyyyy_1, g_yyyy_0_yyyyyyy_1, g_yyyy_0_yyyyyyz_1, g_yyyy_0_yyyyyz_1, g_yyyy_0_yyyyyzz_1, g_yyyy_0_yyyyzz_1, g_yyyy_0_yyyyzzz_1, g_yyyy_0_yyyzzz_1, g_yyyy_0_yyyzzzz_1, g_yyyy_0_yyzzzz_1, g_yyyy_0_yyzzzzz_1, g_yyyy_0_yzzzzz_1, g_yyyy_0_yzzzzzz_1, g_yyyy_0_zzzzzz_1, g_yyyy_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyy_0_xxxxxxx_0[i] = 7.0 * g_yyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxxy_0[i] = 6.0 * g_yyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxxz_0[i] = 6.0 * g_yyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxyy_0[i] = 5.0 * g_yyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxyz_0[i] = 5.0 * g_yyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxzz_0[i] = 5.0 * g_yyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyyy_0[i] = 4.0 * g_yyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyyz_0[i] = 4.0 * g_yyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyzz_0[i] = 4.0 * g_yyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxzzz_0[i] = 4.0 * g_yyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyyy_0[i] = 3.0 * g_yyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyyz_0[i] = 3.0 * g_yyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyzz_0[i] = 3.0 * g_yyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyzzz_0[i] = 3.0 * g_yyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxzzzz_0[i] = 3.0 * g_yyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyyy_0[i] = 2.0 * g_yyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyyz_0[i] = 2.0 * g_yyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyzz_0[i] = 2.0 * g_yyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyzzz_0[i] = 2.0 * g_yyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyzzzz_0[i] = 2.0 * g_yyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxzzzzz_0[i] = 2.0 * g_yyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyyy_0[i] = g_yyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyyz_0[i] = g_yyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyzz_0[i] = g_yyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyzzz_0[i] = g_yyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyzzzz_0[i] = g_yyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyzzzzz_0[i] = g_yyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xzzzzzz_0[i] = g_yyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyyy_0[i] = g_yyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyyz_0[i] = g_yyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyzz_0[i] = g_yyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyzzz_0[i] = g_yyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyzzzz_0[i] = g_yyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyzzzzz_0[i] = g_yyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yzzzzzz_0[i] = g_yyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_zzzzzzz_0[i] = g_yyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 396-432 components of targeted buffer : HSK

    auto g_xyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 396);

    auto g_xyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 397);

    auto g_xyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 398);

    auto g_xyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 399);

    auto g_xyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 400);

    auto g_xyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 401);

    auto g_xyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 402);

    auto g_xyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 403);

    auto g_xyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 404);

    auto g_xyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 405);

    auto g_xyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 406);

    auto g_xyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 407);

    auto g_xyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 408);

    auto g_xyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 409);

    auto g_xyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 410);

    auto g_xyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 411);

    auto g_xyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 412);

    auto g_xyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 413);

    auto g_xyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 414);

    auto g_xyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 415);

    auto g_xyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 416);

    auto g_xyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 417);

    auto g_xyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 418);

    auto g_xyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 419);

    auto g_xyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 420);

    auto g_xyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 421);

    auto g_xyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 422);

    auto g_xyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 423);

    auto g_xyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 424);

    auto g_xyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 425);

    auto g_xyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 426);

    auto g_xyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 427);

    auto g_xyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 428);

    auto g_xyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 429);

    auto g_xyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 430);

    auto g_xyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 431);

    #pragma omp simd aligned(g_xyyy_0_xxxxxxx_1, g_xyyy_0_xxxxxxy_1, g_xyyy_0_xxxxxyy_1, g_xyyy_0_xxxxyyy_1, g_xyyy_0_xxxyyyy_1, g_xyyy_0_xxyyyyy_1, g_xyyy_0_xyyyyyy_1, g_xyyyz_0_xxxxxxx_0, g_xyyyz_0_xxxxxxy_0, g_xyyyz_0_xxxxxxz_0, g_xyyyz_0_xxxxxyy_0, g_xyyyz_0_xxxxxyz_0, g_xyyyz_0_xxxxxzz_0, g_xyyyz_0_xxxxyyy_0, g_xyyyz_0_xxxxyyz_0, g_xyyyz_0_xxxxyzz_0, g_xyyyz_0_xxxxzzz_0, g_xyyyz_0_xxxyyyy_0, g_xyyyz_0_xxxyyyz_0, g_xyyyz_0_xxxyyzz_0, g_xyyyz_0_xxxyzzz_0, g_xyyyz_0_xxxzzzz_0, g_xyyyz_0_xxyyyyy_0, g_xyyyz_0_xxyyyyz_0, g_xyyyz_0_xxyyyzz_0, g_xyyyz_0_xxyyzzz_0, g_xyyyz_0_xxyzzzz_0, g_xyyyz_0_xxzzzzz_0, g_xyyyz_0_xyyyyyy_0, g_xyyyz_0_xyyyyyz_0, g_xyyyz_0_xyyyyzz_0, g_xyyyz_0_xyyyzzz_0, g_xyyyz_0_xyyzzzz_0, g_xyyyz_0_xyzzzzz_0, g_xyyyz_0_xzzzzzz_0, g_xyyyz_0_yyyyyyy_0, g_xyyyz_0_yyyyyyz_0, g_xyyyz_0_yyyyyzz_0, g_xyyyz_0_yyyyzzz_0, g_xyyyz_0_yyyzzzz_0, g_xyyyz_0_yyzzzzz_0, g_xyyyz_0_yzzzzzz_0, g_xyyyz_0_zzzzzzz_0, g_yyyz_0_xxxxxxz_1, g_yyyz_0_xxxxxyz_1, g_yyyz_0_xxxxxz_1, g_yyyz_0_xxxxxzz_1, g_yyyz_0_xxxxyyz_1, g_yyyz_0_xxxxyz_1, g_yyyz_0_xxxxyzz_1, g_yyyz_0_xxxxzz_1, g_yyyz_0_xxxxzzz_1, g_yyyz_0_xxxyyyz_1, g_yyyz_0_xxxyyz_1, g_yyyz_0_xxxyyzz_1, g_yyyz_0_xxxyzz_1, g_yyyz_0_xxxyzzz_1, g_yyyz_0_xxxzzz_1, g_yyyz_0_xxxzzzz_1, g_yyyz_0_xxyyyyz_1, g_yyyz_0_xxyyyz_1, g_yyyz_0_xxyyyzz_1, g_yyyz_0_xxyyzz_1, g_yyyz_0_xxyyzzz_1, g_yyyz_0_xxyzzz_1, g_yyyz_0_xxyzzzz_1, g_yyyz_0_xxzzzz_1, g_yyyz_0_xxzzzzz_1, g_yyyz_0_xyyyyyz_1, g_yyyz_0_xyyyyz_1, g_yyyz_0_xyyyyzz_1, g_yyyz_0_xyyyzz_1, g_yyyz_0_xyyyzzz_1, g_yyyz_0_xyyzzz_1, g_yyyz_0_xyyzzzz_1, g_yyyz_0_xyzzzz_1, g_yyyz_0_xyzzzzz_1, g_yyyz_0_xzzzzz_1, g_yyyz_0_xzzzzzz_1, g_yyyz_0_yyyyyyy_1, g_yyyz_0_yyyyyyz_1, g_yyyz_0_yyyyyz_1, g_yyyz_0_yyyyyzz_1, g_yyyz_0_yyyyzz_1, g_yyyz_0_yyyyzzz_1, g_yyyz_0_yyyzzz_1, g_yyyz_0_yyyzzzz_1, g_yyyz_0_yyzzzz_1, g_yyyz_0_yyzzzzz_1, g_yyyz_0_yzzzzz_1, g_yyyz_0_yzzzzzz_1, g_yyyz_0_zzzzzz_1, g_yyyz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyz_0_xxxxxxx_0[i] = g_xyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxxy_0[i] = g_xyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxxz_0[i] = 6.0 * g_yyyz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxxyy_0[i] = g_xyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxyz_0[i] = 5.0 * g_yyyz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxxzz_0[i] = 5.0 * g_yyyz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxyyy_0[i] = g_xyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxyyz_0[i] = 4.0 * g_yyyz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxyzz_0[i] = 4.0 * g_yyyz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxzzz_0[i] = 4.0 * g_yyyz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyyyy_0[i] = g_xyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxyyyz_0[i] = 3.0 * g_yyyz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyyzz_0[i] = 3.0 * g_yyyz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyzzz_0[i] = 3.0 * g_yyyz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxzzzz_0[i] = 3.0 * g_yyyz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyyyy_0[i] = g_xyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxyyyyz_0[i] = 2.0 * g_yyyz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyyzz_0[i] = 2.0 * g_yyyz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyzzz_0[i] = 2.0 * g_yyyz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyzzzz_0[i] = 2.0 * g_yyyz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxzzzzz_0[i] = 2.0 * g_yyyz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyyyy_0[i] = g_xyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xyyyyyz_0[i] = g_yyyz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyyzz_0[i] = g_yyyz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyzzz_0[i] = g_yyyz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyzzzz_0[i] = g_yyyz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyzzzzz_0[i] = g_yyyz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xzzzzzz_0[i] = g_yyyz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyyy_0[i] = g_yyyz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyyz_0[i] = g_yyyz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyzz_0[i] = g_yyyz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyzzz_0[i] = g_yyyz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyzzzz_0[i] = g_yyyz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyzzzzz_0[i] = g_yyyz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yzzzzzz_0[i] = g_yyyz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_zzzzzzz_0[i] = g_yyyz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 432-468 components of targeted buffer : HSK

    auto g_xyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 432);

    auto g_xyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 433);

    auto g_xyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 434);

    auto g_xyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 435);

    auto g_xyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 436);

    auto g_xyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 437);

    auto g_xyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 438);

    auto g_xyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 439);

    auto g_xyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 440);

    auto g_xyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 441);

    auto g_xyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 442);

    auto g_xyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 443);

    auto g_xyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 444);

    auto g_xyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 445);

    auto g_xyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 446);

    auto g_xyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 447);

    auto g_xyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 448);

    auto g_xyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 449);

    auto g_xyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 450);

    auto g_xyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 451);

    auto g_xyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 452);

    auto g_xyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 453);

    auto g_xyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 454);

    auto g_xyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 455);

    auto g_xyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 456);

    auto g_xyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 457);

    auto g_xyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 458);

    auto g_xyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 459);

    auto g_xyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 460);

    auto g_xyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 461);

    auto g_xyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 462);

    auto g_xyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 463);

    auto g_xyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 464);

    auto g_xyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 465);

    auto g_xyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 466);

    auto g_xyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 467);

    #pragma omp simd aligned(g_xyyzz_0_xxxxxxx_0, g_xyyzz_0_xxxxxxy_0, g_xyyzz_0_xxxxxxz_0, g_xyyzz_0_xxxxxyy_0, g_xyyzz_0_xxxxxyz_0, g_xyyzz_0_xxxxxzz_0, g_xyyzz_0_xxxxyyy_0, g_xyyzz_0_xxxxyyz_0, g_xyyzz_0_xxxxyzz_0, g_xyyzz_0_xxxxzzz_0, g_xyyzz_0_xxxyyyy_0, g_xyyzz_0_xxxyyyz_0, g_xyyzz_0_xxxyyzz_0, g_xyyzz_0_xxxyzzz_0, g_xyyzz_0_xxxzzzz_0, g_xyyzz_0_xxyyyyy_0, g_xyyzz_0_xxyyyyz_0, g_xyyzz_0_xxyyyzz_0, g_xyyzz_0_xxyyzzz_0, g_xyyzz_0_xxyzzzz_0, g_xyyzz_0_xxzzzzz_0, g_xyyzz_0_xyyyyyy_0, g_xyyzz_0_xyyyyyz_0, g_xyyzz_0_xyyyyzz_0, g_xyyzz_0_xyyyzzz_0, g_xyyzz_0_xyyzzzz_0, g_xyyzz_0_xyzzzzz_0, g_xyyzz_0_xzzzzzz_0, g_xyyzz_0_yyyyyyy_0, g_xyyzz_0_yyyyyyz_0, g_xyyzz_0_yyyyyzz_0, g_xyyzz_0_yyyyzzz_0, g_xyyzz_0_yyyzzzz_0, g_xyyzz_0_yyzzzzz_0, g_xyyzz_0_yzzzzzz_0, g_xyyzz_0_zzzzzzz_0, g_yyzz_0_xxxxxx_1, g_yyzz_0_xxxxxxx_1, g_yyzz_0_xxxxxxy_1, g_yyzz_0_xxxxxxz_1, g_yyzz_0_xxxxxy_1, g_yyzz_0_xxxxxyy_1, g_yyzz_0_xxxxxyz_1, g_yyzz_0_xxxxxz_1, g_yyzz_0_xxxxxzz_1, g_yyzz_0_xxxxyy_1, g_yyzz_0_xxxxyyy_1, g_yyzz_0_xxxxyyz_1, g_yyzz_0_xxxxyz_1, g_yyzz_0_xxxxyzz_1, g_yyzz_0_xxxxzz_1, g_yyzz_0_xxxxzzz_1, g_yyzz_0_xxxyyy_1, g_yyzz_0_xxxyyyy_1, g_yyzz_0_xxxyyyz_1, g_yyzz_0_xxxyyz_1, g_yyzz_0_xxxyyzz_1, g_yyzz_0_xxxyzz_1, g_yyzz_0_xxxyzzz_1, g_yyzz_0_xxxzzz_1, g_yyzz_0_xxxzzzz_1, g_yyzz_0_xxyyyy_1, g_yyzz_0_xxyyyyy_1, g_yyzz_0_xxyyyyz_1, g_yyzz_0_xxyyyz_1, g_yyzz_0_xxyyyzz_1, g_yyzz_0_xxyyzz_1, g_yyzz_0_xxyyzzz_1, g_yyzz_0_xxyzzz_1, g_yyzz_0_xxyzzzz_1, g_yyzz_0_xxzzzz_1, g_yyzz_0_xxzzzzz_1, g_yyzz_0_xyyyyy_1, g_yyzz_0_xyyyyyy_1, g_yyzz_0_xyyyyyz_1, g_yyzz_0_xyyyyz_1, g_yyzz_0_xyyyyzz_1, g_yyzz_0_xyyyzz_1, g_yyzz_0_xyyyzzz_1, g_yyzz_0_xyyzzz_1, g_yyzz_0_xyyzzzz_1, g_yyzz_0_xyzzzz_1, g_yyzz_0_xyzzzzz_1, g_yyzz_0_xzzzzz_1, g_yyzz_0_xzzzzzz_1, g_yyzz_0_yyyyyy_1, g_yyzz_0_yyyyyyy_1, g_yyzz_0_yyyyyyz_1, g_yyzz_0_yyyyyz_1, g_yyzz_0_yyyyyzz_1, g_yyzz_0_yyyyzz_1, g_yyzz_0_yyyyzzz_1, g_yyzz_0_yyyzzz_1, g_yyzz_0_yyyzzzz_1, g_yyzz_0_yyzzzz_1, g_yyzz_0_yyzzzzz_1, g_yyzz_0_yzzzzz_1, g_yyzz_0_yzzzzzz_1, g_yyzz_0_zzzzzz_1, g_yyzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzz_0_xxxxxxx_0[i] = 7.0 * g_yyzz_0_xxxxxx_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxxy_0[i] = 6.0 * g_yyzz_0_xxxxxy_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxxz_0[i] = 6.0 * g_yyzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxyy_0[i] = 5.0 * g_yyzz_0_xxxxyy_1[i] * fi_acd_0 + g_yyzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxyz_0[i] = 5.0 * g_yyzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxzz_0[i] = 5.0 * g_yyzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyyy_0[i] = 4.0 * g_yyzz_0_xxxyyy_1[i] * fi_acd_0 + g_yyzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyyz_0[i] = 4.0 * g_yyzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyzz_0[i] = 4.0 * g_yyzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxzzz_0[i] = 4.0 * g_yyzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyyy_0[i] = 3.0 * g_yyzz_0_xxyyyy_1[i] * fi_acd_0 + g_yyzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyyz_0[i] = 3.0 * g_yyzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyzz_0[i] = 3.0 * g_yyzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyzzz_0[i] = 3.0 * g_yyzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxzzzz_0[i] = 3.0 * g_yyzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyyy_0[i] = 2.0 * g_yyzz_0_xyyyyy_1[i] * fi_acd_0 + g_yyzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyyz_0[i] = 2.0 * g_yyzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyzz_0[i] = 2.0 * g_yyzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyzzz_0[i] = 2.0 * g_yyzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyzzzz_0[i] = 2.0 * g_yyzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxzzzzz_0[i] = 2.0 * g_yyzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyyy_0[i] = g_yyzz_0_yyyyyy_1[i] * fi_acd_0 + g_yyzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyyz_0[i] = g_yyzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyzz_0[i] = g_yyzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyzzz_0[i] = g_yyzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyzzzz_0[i] = g_yyzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyzzzzz_0[i] = g_yyzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xzzzzzz_0[i] = g_yyzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyyy_0[i] = g_yyzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyyz_0[i] = g_yyzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyzz_0[i] = g_yyzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyzzz_0[i] = g_yyzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyzzzz_0[i] = g_yyzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyzzzzz_0[i] = g_yyzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yzzzzzz_0[i] = g_yyzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_zzzzzzz_0[i] = g_yyzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 468-504 components of targeted buffer : HSK

    auto g_xyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 468);

    auto g_xyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 469);

    auto g_xyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 470);

    auto g_xyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 471);

    auto g_xyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 472);

    auto g_xyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 473);

    auto g_xyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 474);

    auto g_xyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 475);

    auto g_xyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 476);

    auto g_xyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 477);

    auto g_xyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 478);

    auto g_xyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 479);

    auto g_xyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 480);

    auto g_xyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 481);

    auto g_xyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 482);

    auto g_xyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 483);

    auto g_xyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 484);

    auto g_xyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 485);

    auto g_xyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 486);

    auto g_xyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 487);

    auto g_xyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 488);

    auto g_xyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 489);

    auto g_xyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 490);

    auto g_xyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 491);

    auto g_xyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 492);

    auto g_xyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 493);

    auto g_xyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 494);

    auto g_xyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 495);

    auto g_xyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 496);

    auto g_xyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 497);

    auto g_xyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 498);

    auto g_xyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 499);

    auto g_xyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 500);

    auto g_xyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 501);

    auto g_xyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 502);

    auto g_xyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 503);

    #pragma omp simd aligned(g_xyzzz_0_xxxxxxx_0, g_xyzzz_0_xxxxxxy_0, g_xyzzz_0_xxxxxxz_0, g_xyzzz_0_xxxxxyy_0, g_xyzzz_0_xxxxxyz_0, g_xyzzz_0_xxxxxzz_0, g_xyzzz_0_xxxxyyy_0, g_xyzzz_0_xxxxyyz_0, g_xyzzz_0_xxxxyzz_0, g_xyzzz_0_xxxxzzz_0, g_xyzzz_0_xxxyyyy_0, g_xyzzz_0_xxxyyyz_0, g_xyzzz_0_xxxyyzz_0, g_xyzzz_0_xxxyzzz_0, g_xyzzz_0_xxxzzzz_0, g_xyzzz_0_xxyyyyy_0, g_xyzzz_0_xxyyyyz_0, g_xyzzz_0_xxyyyzz_0, g_xyzzz_0_xxyyzzz_0, g_xyzzz_0_xxyzzzz_0, g_xyzzz_0_xxzzzzz_0, g_xyzzz_0_xyyyyyy_0, g_xyzzz_0_xyyyyyz_0, g_xyzzz_0_xyyyyzz_0, g_xyzzz_0_xyyyzzz_0, g_xyzzz_0_xyyzzzz_0, g_xyzzz_0_xyzzzzz_0, g_xyzzz_0_xzzzzzz_0, g_xyzzz_0_yyyyyyy_0, g_xyzzz_0_yyyyyyz_0, g_xyzzz_0_yyyyyzz_0, g_xyzzz_0_yyyyzzz_0, g_xyzzz_0_yyyzzzz_0, g_xyzzz_0_yyzzzzz_0, g_xyzzz_0_yzzzzzz_0, g_xyzzz_0_zzzzzzz_0, g_xzzz_0_xxxxxxx_1, g_xzzz_0_xxxxxxz_1, g_xzzz_0_xxxxxzz_1, g_xzzz_0_xxxxzzz_1, g_xzzz_0_xxxzzzz_1, g_xzzz_0_xxzzzzz_1, g_xzzz_0_xzzzzzz_1, g_yzzz_0_xxxxxxy_1, g_yzzz_0_xxxxxy_1, g_yzzz_0_xxxxxyy_1, g_yzzz_0_xxxxxyz_1, g_yzzz_0_xxxxyy_1, g_yzzz_0_xxxxyyy_1, g_yzzz_0_xxxxyyz_1, g_yzzz_0_xxxxyz_1, g_yzzz_0_xxxxyzz_1, g_yzzz_0_xxxyyy_1, g_yzzz_0_xxxyyyy_1, g_yzzz_0_xxxyyyz_1, g_yzzz_0_xxxyyz_1, g_yzzz_0_xxxyyzz_1, g_yzzz_0_xxxyzz_1, g_yzzz_0_xxxyzzz_1, g_yzzz_0_xxyyyy_1, g_yzzz_0_xxyyyyy_1, g_yzzz_0_xxyyyyz_1, g_yzzz_0_xxyyyz_1, g_yzzz_0_xxyyyzz_1, g_yzzz_0_xxyyzz_1, g_yzzz_0_xxyyzzz_1, g_yzzz_0_xxyzzz_1, g_yzzz_0_xxyzzzz_1, g_yzzz_0_xyyyyy_1, g_yzzz_0_xyyyyyy_1, g_yzzz_0_xyyyyyz_1, g_yzzz_0_xyyyyz_1, g_yzzz_0_xyyyyzz_1, g_yzzz_0_xyyyzz_1, g_yzzz_0_xyyyzzz_1, g_yzzz_0_xyyzzz_1, g_yzzz_0_xyyzzzz_1, g_yzzz_0_xyzzzz_1, g_yzzz_0_xyzzzzz_1, g_yzzz_0_yyyyyy_1, g_yzzz_0_yyyyyyy_1, g_yzzz_0_yyyyyyz_1, g_yzzz_0_yyyyyz_1, g_yzzz_0_yyyyyzz_1, g_yzzz_0_yyyyzz_1, g_yzzz_0_yyyyzzz_1, g_yzzz_0_yyyzzz_1, g_yzzz_0_yyyzzzz_1, g_yzzz_0_yyzzzz_1, g_yzzz_0_yyzzzzz_1, g_yzzz_0_yzzzzz_1, g_yzzz_0_yzzzzzz_1, g_yzzz_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzz_0_xxxxxxx_0[i] = g_xzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xyzzz_0_xxxxxxy_0[i] = 6.0 * g_yzzz_0_xxxxxy_1[i] * fi_acd_0 + g_yzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxxz_0[i] = g_xzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xyzzz_0_xxxxxyy_0[i] = 5.0 * g_yzzz_0_xxxxyy_1[i] * fi_acd_0 + g_yzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxyz_0[i] = 5.0 * g_yzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxzz_0[i] = g_xzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xyzzz_0_xxxxyyy_0[i] = 4.0 * g_yzzz_0_xxxyyy_1[i] * fi_acd_0 + g_yzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxyyz_0[i] = 4.0 * g_yzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxyzz_0[i] = 4.0 * g_yzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxzzz_0[i] = g_xzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xyzzz_0_xxxyyyy_0[i] = 3.0 * g_yzzz_0_xxyyyy_1[i] * fi_acd_0 + g_yzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxyyyz_0[i] = 3.0 * g_yzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxyyzz_0[i] = 3.0 * g_yzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxyzzz_0[i] = 3.0 * g_yzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxzzzz_0[i] = g_xzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xyzzz_0_xxyyyyy_0[i] = 2.0 * g_yzzz_0_xyyyyy_1[i] * fi_acd_0 + g_yzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxyyyyz_0[i] = 2.0 * g_yzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxyyyzz_0[i] = 2.0 * g_yzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxyyzzz_0[i] = 2.0 * g_yzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxyzzzz_0[i] = 2.0 * g_yzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxzzzzz_0[i] = g_xzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xyzzz_0_xyyyyyy_0[i] = g_yzzz_0_yyyyyy_1[i] * fi_acd_0 + g_yzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xyyyyyz_0[i] = g_yzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xyyyyzz_0[i] = g_yzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xyyyzzz_0[i] = g_yzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xyyzzzz_0[i] = g_yzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xyzzzzz_0[i] = g_yzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xzzzzzz_0[i] = g_xzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xyzzz_0_yyyyyyy_0[i] = g_yzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_yyyyyyz_0[i] = g_yzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_yyyyyzz_0[i] = g_yzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_yyyyzzz_0[i] = g_yzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_yyyzzzz_0[i] = g_yzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_yyzzzzz_0[i] = g_yzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_yzzzzzz_0[i] = g_yzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_zzzzzzz_0[i] = g_yzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 504-540 components of targeted buffer : HSK

    auto g_xzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 504);

    auto g_xzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 505);

    auto g_xzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 506);

    auto g_xzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 507);

    auto g_xzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 508);

    auto g_xzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 509);

    auto g_xzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 510);

    auto g_xzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 511);

    auto g_xzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 512);

    auto g_xzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 513);

    auto g_xzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 514);

    auto g_xzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 515);

    auto g_xzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 516);

    auto g_xzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 517);

    auto g_xzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 518);

    auto g_xzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 519);

    auto g_xzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 520);

    auto g_xzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 521);

    auto g_xzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 522);

    auto g_xzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 523);

    auto g_xzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 524);

    auto g_xzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 525);

    auto g_xzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 526);

    auto g_xzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 527);

    auto g_xzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 528);

    auto g_xzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 529);

    auto g_xzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 530);

    auto g_xzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 531);

    auto g_xzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 532);

    auto g_xzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 533);

    auto g_xzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 534);

    auto g_xzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 535);

    auto g_xzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 536);

    auto g_xzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 537);

    auto g_xzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 538);

    auto g_xzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 539);

    #pragma omp simd aligned(g_xzzzz_0_xxxxxxx_0, g_xzzzz_0_xxxxxxy_0, g_xzzzz_0_xxxxxxz_0, g_xzzzz_0_xxxxxyy_0, g_xzzzz_0_xxxxxyz_0, g_xzzzz_0_xxxxxzz_0, g_xzzzz_0_xxxxyyy_0, g_xzzzz_0_xxxxyyz_0, g_xzzzz_0_xxxxyzz_0, g_xzzzz_0_xxxxzzz_0, g_xzzzz_0_xxxyyyy_0, g_xzzzz_0_xxxyyyz_0, g_xzzzz_0_xxxyyzz_0, g_xzzzz_0_xxxyzzz_0, g_xzzzz_0_xxxzzzz_0, g_xzzzz_0_xxyyyyy_0, g_xzzzz_0_xxyyyyz_0, g_xzzzz_0_xxyyyzz_0, g_xzzzz_0_xxyyzzz_0, g_xzzzz_0_xxyzzzz_0, g_xzzzz_0_xxzzzzz_0, g_xzzzz_0_xyyyyyy_0, g_xzzzz_0_xyyyyyz_0, g_xzzzz_0_xyyyyzz_0, g_xzzzz_0_xyyyzzz_0, g_xzzzz_0_xyyzzzz_0, g_xzzzz_0_xyzzzzz_0, g_xzzzz_0_xzzzzzz_0, g_xzzzz_0_yyyyyyy_0, g_xzzzz_0_yyyyyyz_0, g_xzzzz_0_yyyyyzz_0, g_xzzzz_0_yyyyzzz_0, g_xzzzz_0_yyyzzzz_0, g_xzzzz_0_yyzzzzz_0, g_xzzzz_0_yzzzzzz_0, g_xzzzz_0_zzzzzzz_0, g_zzzz_0_xxxxxx_1, g_zzzz_0_xxxxxxx_1, g_zzzz_0_xxxxxxy_1, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxy_1, g_zzzz_0_xxxxxyy_1, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxyy_1, g_zzzz_0_xxxxyyy_1, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxyyy_1, g_zzzz_0_xxxyyyy_1, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxyyyy_1, g_zzzz_0_xxyyyyy_1, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxzzzz_1, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xyyyyy_1, g_zzzz_0_xyyyyyy_1, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyzzzz_1, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xzzzzz_1, g_zzzz_0_xzzzzzz_1, g_zzzz_0_yyyyyy_1, g_zzzz_0_yyyyyyy_1, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyzzzz_1, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yzzzzz_1, g_zzzz_0_yzzzzzz_1, g_zzzz_0_zzzzzz_1, g_zzzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzz_0_xxxxxxx_0[i] = 7.0 * g_zzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxxy_0[i] = 6.0 * g_zzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxxz_0[i] = 6.0 * g_zzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxyy_0[i] = 5.0 * g_zzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxyz_0[i] = 5.0 * g_zzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxzz_0[i] = 5.0 * g_zzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyyy_0[i] = 4.0 * g_zzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyyz_0[i] = 4.0 * g_zzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyzz_0[i] = 4.0 * g_zzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxzzz_0[i] = 4.0 * g_zzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyyy_0[i] = 3.0 * g_zzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyyz_0[i] = 3.0 * g_zzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyzz_0[i] = 3.0 * g_zzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyzzz_0[i] = 3.0 * g_zzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxzzzz_0[i] = 3.0 * g_zzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyyy_0[i] = 2.0 * g_zzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyyz_0[i] = 2.0 * g_zzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyzz_0[i] = 2.0 * g_zzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyzzz_0[i] = 2.0 * g_zzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyzzzz_0[i] = 2.0 * g_zzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxzzzzz_0[i] = 2.0 * g_zzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyyy_0[i] = g_zzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyyz_0[i] = g_zzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyzz_0[i] = g_zzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyzzz_0[i] = g_zzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyzzzz_0[i] = g_zzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyzzzzz_0[i] = g_zzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xzzzzzz_0[i] = g_zzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyyy_0[i] = g_zzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyyz_0[i] = g_zzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyzz_0[i] = g_zzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyzzz_0[i] = g_zzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyzzzz_0[i] = g_zzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyzzzzz_0[i] = g_zzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yzzzzzz_0[i] = g_zzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_zzzzzzz_0[i] = g_zzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 540-576 components of targeted buffer : HSK

    auto g_yyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 540);

    auto g_yyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 541);

    auto g_yyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 542);

    auto g_yyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 543);

    auto g_yyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 544);

    auto g_yyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 545);

    auto g_yyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 546);

    auto g_yyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 547);

    auto g_yyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 548);

    auto g_yyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 549);

    auto g_yyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 550);

    auto g_yyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 551);

    auto g_yyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 552);

    auto g_yyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 553);

    auto g_yyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 554);

    auto g_yyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 555);

    auto g_yyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 556);

    auto g_yyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 557);

    auto g_yyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 558);

    auto g_yyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 559);

    auto g_yyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 560);

    auto g_yyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 561);

    auto g_yyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 562);

    auto g_yyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 563);

    auto g_yyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 564);

    auto g_yyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 565);

    auto g_yyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 566);

    auto g_yyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 567);

    auto g_yyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 568);

    auto g_yyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 569);

    auto g_yyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 570);

    auto g_yyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 571);

    auto g_yyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 572);

    auto g_yyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 573);

    auto g_yyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 574);

    auto g_yyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 575);

    #pragma omp simd aligned(g_yyy_0_xxxxxxx_0, g_yyy_0_xxxxxxx_1, g_yyy_0_xxxxxxy_0, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxxz_0, g_yyy_0_xxxxxxz_1, g_yyy_0_xxxxxyy_0, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxxyz_0, g_yyy_0_xxxxxyz_1, g_yyy_0_xxxxxzz_0, g_yyy_0_xxxxxzz_1, g_yyy_0_xxxxyyy_0, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxxyyz_0, g_yyy_0_xxxxyyz_1, g_yyy_0_xxxxyzz_0, g_yyy_0_xxxxyzz_1, g_yyy_0_xxxxzzz_0, g_yyy_0_xxxxzzz_1, g_yyy_0_xxxyyyy_0, g_yyy_0_xxxyyyy_1, g_yyy_0_xxxyyyz_0, g_yyy_0_xxxyyyz_1, g_yyy_0_xxxyyzz_0, g_yyy_0_xxxyyzz_1, g_yyy_0_xxxyzzz_0, g_yyy_0_xxxyzzz_1, g_yyy_0_xxxzzzz_0, g_yyy_0_xxxzzzz_1, g_yyy_0_xxyyyyy_0, g_yyy_0_xxyyyyy_1, g_yyy_0_xxyyyyz_0, g_yyy_0_xxyyyyz_1, g_yyy_0_xxyyyzz_0, g_yyy_0_xxyyyzz_1, g_yyy_0_xxyyzzz_0, g_yyy_0_xxyyzzz_1, g_yyy_0_xxyzzzz_0, g_yyy_0_xxyzzzz_1, g_yyy_0_xxzzzzz_0, g_yyy_0_xxzzzzz_1, g_yyy_0_xyyyyyy_0, g_yyy_0_xyyyyyy_1, g_yyy_0_xyyyyyz_0, g_yyy_0_xyyyyyz_1, g_yyy_0_xyyyyzz_0, g_yyy_0_xyyyyzz_1, g_yyy_0_xyyyzzz_0, g_yyy_0_xyyyzzz_1, g_yyy_0_xyyzzzz_0, g_yyy_0_xyyzzzz_1, g_yyy_0_xyzzzzz_0, g_yyy_0_xyzzzzz_1, g_yyy_0_xzzzzzz_0, g_yyy_0_xzzzzzz_1, g_yyy_0_yyyyyyy_0, g_yyy_0_yyyyyyy_1, g_yyy_0_yyyyyyz_0, g_yyy_0_yyyyyyz_1, g_yyy_0_yyyyyzz_0, g_yyy_0_yyyyyzz_1, g_yyy_0_yyyyzzz_0, g_yyy_0_yyyyzzz_1, g_yyy_0_yyyzzzz_0, g_yyy_0_yyyzzzz_1, g_yyy_0_yyzzzzz_0, g_yyy_0_yyzzzzz_1, g_yyy_0_yzzzzzz_0, g_yyy_0_yzzzzzz_1, g_yyy_0_zzzzzzz_0, g_yyy_0_zzzzzzz_1, g_yyyy_0_xxxxxx_1, g_yyyy_0_xxxxxxx_1, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxxz_1, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxxyz_1, g_yyyy_0_xxxxxz_1, g_yyyy_0_xxxxxzz_1, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxxyyz_1, g_yyyy_0_xxxxyz_1, g_yyyy_0_xxxxyzz_1, g_yyyy_0_xxxxzz_1, g_yyyy_0_xxxxzzz_1, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxxyyyz_1, g_yyyy_0_xxxyyz_1, g_yyyy_0_xxxyyzz_1, g_yyyy_0_xxxyzz_1, g_yyyy_0_xxxyzzz_1, g_yyyy_0_xxxzzz_1, g_yyyy_0_xxxzzzz_1, g_yyyy_0_xxyyyy_1, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xxyyyyz_1, g_yyyy_0_xxyyyz_1, g_yyyy_0_xxyyyzz_1, g_yyyy_0_xxyyzz_1, g_yyyy_0_xxyyzzz_1, g_yyyy_0_xxyzzz_1, g_yyyy_0_xxyzzzz_1, g_yyyy_0_xxzzzz_1, g_yyyy_0_xxzzzzz_1, g_yyyy_0_xyyyyy_1, g_yyyy_0_xyyyyyy_1, g_yyyy_0_xyyyyyz_1, g_yyyy_0_xyyyyz_1, g_yyyy_0_xyyyyzz_1, g_yyyy_0_xyyyzz_1, g_yyyy_0_xyyyzzz_1, g_yyyy_0_xyyzzz_1, g_yyyy_0_xyyzzzz_1, g_yyyy_0_xyzzzz_1, g_yyyy_0_xyzzzzz_1, g_yyyy_0_xzzzzz_1, g_yyyy_0_xzzzzzz_1, g_yyyy_0_yyyyyy_1, g_yyyy_0_yyyyyyy_1, g_yyyy_0_yyyyyyz_1, g_yyyy_0_yyyyyz_1, g_yyyy_0_yyyyyzz_1, g_yyyy_0_yyyyzz_1, g_yyyy_0_yyyyzzz_1, g_yyyy_0_yyyzzz_1, g_yyyy_0_yyyzzzz_1, g_yyyy_0_yyzzzz_1, g_yyyy_0_yyzzzzz_1, g_yyyy_0_yzzzzz_1, g_yyyy_0_yzzzzzz_1, g_yyyy_0_zzzzzz_1, g_yyyy_0_zzzzzzz_1, g_yyyyy_0_xxxxxxx_0, g_yyyyy_0_xxxxxxy_0, g_yyyyy_0_xxxxxxz_0, g_yyyyy_0_xxxxxyy_0, g_yyyyy_0_xxxxxyz_0, g_yyyyy_0_xxxxxzz_0, g_yyyyy_0_xxxxyyy_0, g_yyyyy_0_xxxxyyz_0, g_yyyyy_0_xxxxyzz_0, g_yyyyy_0_xxxxzzz_0, g_yyyyy_0_xxxyyyy_0, g_yyyyy_0_xxxyyyz_0, g_yyyyy_0_xxxyyzz_0, g_yyyyy_0_xxxyzzz_0, g_yyyyy_0_xxxzzzz_0, g_yyyyy_0_xxyyyyy_0, g_yyyyy_0_xxyyyyz_0, g_yyyyy_0_xxyyyzz_0, g_yyyyy_0_xxyyzzz_0, g_yyyyy_0_xxyzzzz_0, g_yyyyy_0_xxzzzzz_0, g_yyyyy_0_xyyyyyy_0, g_yyyyy_0_xyyyyyz_0, g_yyyyy_0_xyyyyzz_0, g_yyyyy_0_xyyyzzz_0, g_yyyyy_0_xyyzzzz_0, g_yyyyy_0_xyzzzzz_0, g_yyyyy_0_xzzzzzz_0, g_yyyyy_0_yyyyyyy_0, g_yyyyy_0_yyyyyyz_0, g_yyyyy_0_yyyyyzz_0, g_yyyyy_0_yyyyzzz_0, g_yyyyy_0_yyyzzzz_0, g_yyyyy_0_yyzzzzz_0, g_yyyyy_0_yzzzzzz_0, g_yyyyy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyy_0_xxxxxxx_0[i] = 4.0 * g_yyy_0_xxxxxxx_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxx_1[i] * fz_be_0 + g_yyyy_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxxy_0[i] = 4.0 * g_yyy_0_xxxxxxy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxy_1[i] * fz_be_0 + g_yyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxxz_0[i] = 4.0 * g_yyy_0_xxxxxxz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxz_1[i] * fz_be_0 + g_yyyy_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxyy_0[i] = 4.0 * g_yyy_0_xxxxxyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxyz_0[i] = 4.0 * g_yyy_0_xxxxxyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxyz_1[i] * fz_be_0 + g_yyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxzz_0[i] = 4.0 * g_yyy_0_xxxxxzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxzz_1[i] * fz_be_0 + g_yyyy_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyyy_0[i] = 4.0 * g_yyy_0_xxxxyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyyz_0[i] = 4.0 * g_yyy_0_xxxxyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyzz_0[i] = 4.0 * g_yyy_0_xxxxyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyzz_1[i] * fz_be_0 + g_yyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxzzz_0[i] = 4.0 * g_yyy_0_xxxxzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxzzz_1[i] * fz_be_0 + g_yyyy_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyyy_0[i] = 4.0 * g_yyy_0_xxxyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyyz_0[i] = 4.0 * g_yyy_0_xxxyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyzz_0[i] = 4.0 * g_yyy_0_xxxyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyzzz_0[i] = 4.0 * g_yyy_0_xxxyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyzzz_1[i] * fz_be_0 + g_yyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxzzzz_0[i] = 4.0 * g_yyy_0_xxxzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxzzzz_1[i] * fz_be_0 + g_yyyy_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyyy_0[i] = 4.0 * g_yyy_0_xxyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyyz_0[i] = 4.0 * g_yyy_0_xxyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyzz_0[i] = 4.0 * g_yyy_0_xxyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyzzz_0[i] = 4.0 * g_yyy_0_xxyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyzzzz_0[i] = 4.0 * g_yyy_0_xxyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyzzzz_1[i] * fz_be_0 + g_yyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxzzzzz_0[i] = 4.0 * g_yyy_0_xxzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxzzzzz_1[i] * fz_be_0 + g_yyyy_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyyy_0[i] = 4.0 * g_yyy_0_xyyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyyz_0[i] = 4.0 * g_yyy_0_xyyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyzz_0[i] = 4.0 * g_yyy_0_xyyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyzzz_0[i] = 4.0 * g_yyy_0_xyyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyzzzz_0[i] = 4.0 * g_yyy_0_xyyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyzzzzz_0[i] = 4.0 * g_yyy_0_xyzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyzzzzz_1[i] * fz_be_0 + g_yyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xzzzzzz_0[i] = 4.0 * g_yyy_0_xzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xzzzzzz_1[i] * fz_be_0 + g_yyyy_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyyy_0[i] = 4.0 * g_yyy_0_yyyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyy_0_yyyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyyz_0[i] = 4.0 * g_yyy_0_yyyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyzz_0[i] = 4.0 * g_yyy_0_yyyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyzzz_0[i] = 4.0 * g_yyy_0_yyyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyzzzz_0[i] = 4.0 * g_yyy_0_yyyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyzzzzz_0[i] = 4.0 * g_yyy_0_yyzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yzzzzzz_0[i] = 4.0 * g_yyy_0_yzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yzzzzzz_1[i] * fz_be_0 + g_yyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyy_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_zzzzzzz_0[i] = 4.0 * g_yyy_0_zzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_zzzzzzz_1[i] * fz_be_0 + g_yyyy_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 576-612 components of targeted buffer : HSK

    auto g_yyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 576);

    auto g_yyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 577);

    auto g_yyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 578);

    auto g_yyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 579);

    auto g_yyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 580);

    auto g_yyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 581);

    auto g_yyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 582);

    auto g_yyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 583);

    auto g_yyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 584);

    auto g_yyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 585);

    auto g_yyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 586);

    auto g_yyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 587);

    auto g_yyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 588);

    auto g_yyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 589);

    auto g_yyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 590);

    auto g_yyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 591);

    auto g_yyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 592);

    auto g_yyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 593);

    auto g_yyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 594);

    auto g_yyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 595);

    auto g_yyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 596);

    auto g_yyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 597);

    auto g_yyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 598);

    auto g_yyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 599);

    auto g_yyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 600);

    auto g_yyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 601);

    auto g_yyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 602);

    auto g_yyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 603);

    auto g_yyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 604);

    auto g_yyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 605);

    auto g_yyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 606);

    auto g_yyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 607);

    auto g_yyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 608);

    auto g_yyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 609);

    auto g_yyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 610);

    auto g_yyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 611);

    #pragma omp simd aligned(g_yyyy_0_xxxxxx_1, g_yyyy_0_xxxxxxx_1, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxxz_1, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxxyz_1, g_yyyy_0_xxxxxz_1, g_yyyy_0_xxxxxzz_1, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxxyyz_1, g_yyyy_0_xxxxyz_1, g_yyyy_0_xxxxyzz_1, g_yyyy_0_xxxxzz_1, g_yyyy_0_xxxxzzz_1, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxxyyyz_1, g_yyyy_0_xxxyyz_1, g_yyyy_0_xxxyyzz_1, g_yyyy_0_xxxyzz_1, g_yyyy_0_xxxyzzz_1, g_yyyy_0_xxxzzz_1, g_yyyy_0_xxxzzzz_1, g_yyyy_0_xxyyyy_1, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xxyyyyz_1, g_yyyy_0_xxyyyz_1, g_yyyy_0_xxyyyzz_1, g_yyyy_0_xxyyzz_1, g_yyyy_0_xxyyzzz_1, g_yyyy_0_xxyzzz_1, g_yyyy_0_xxyzzzz_1, g_yyyy_0_xxzzzz_1, g_yyyy_0_xxzzzzz_1, g_yyyy_0_xyyyyy_1, g_yyyy_0_xyyyyyy_1, g_yyyy_0_xyyyyyz_1, g_yyyy_0_xyyyyz_1, g_yyyy_0_xyyyyzz_1, g_yyyy_0_xyyyzz_1, g_yyyy_0_xyyyzzz_1, g_yyyy_0_xyyzzz_1, g_yyyy_0_xyyzzzz_1, g_yyyy_0_xyzzzz_1, g_yyyy_0_xyzzzzz_1, g_yyyy_0_xzzzzz_1, g_yyyy_0_xzzzzzz_1, g_yyyy_0_yyyyyy_1, g_yyyy_0_yyyyyyy_1, g_yyyy_0_yyyyyyz_1, g_yyyy_0_yyyyyz_1, g_yyyy_0_yyyyyzz_1, g_yyyy_0_yyyyzz_1, g_yyyy_0_yyyyzzz_1, g_yyyy_0_yyyzzz_1, g_yyyy_0_yyyzzzz_1, g_yyyy_0_yyzzzz_1, g_yyyy_0_yyzzzzz_1, g_yyyy_0_yzzzzz_1, g_yyyy_0_yzzzzzz_1, g_yyyy_0_zzzzzz_1, g_yyyy_0_zzzzzzz_1, g_yyyyz_0_xxxxxxx_0, g_yyyyz_0_xxxxxxy_0, g_yyyyz_0_xxxxxxz_0, g_yyyyz_0_xxxxxyy_0, g_yyyyz_0_xxxxxyz_0, g_yyyyz_0_xxxxxzz_0, g_yyyyz_0_xxxxyyy_0, g_yyyyz_0_xxxxyyz_0, g_yyyyz_0_xxxxyzz_0, g_yyyyz_0_xxxxzzz_0, g_yyyyz_0_xxxyyyy_0, g_yyyyz_0_xxxyyyz_0, g_yyyyz_0_xxxyyzz_0, g_yyyyz_0_xxxyzzz_0, g_yyyyz_0_xxxzzzz_0, g_yyyyz_0_xxyyyyy_0, g_yyyyz_0_xxyyyyz_0, g_yyyyz_0_xxyyyzz_0, g_yyyyz_0_xxyyzzz_0, g_yyyyz_0_xxyzzzz_0, g_yyyyz_0_xxzzzzz_0, g_yyyyz_0_xyyyyyy_0, g_yyyyz_0_xyyyyyz_0, g_yyyyz_0_xyyyyzz_0, g_yyyyz_0_xyyyzzz_0, g_yyyyz_0_xyyzzzz_0, g_yyyyz_0_xyzzzzz_0, g_yyyyz_0_xzzzzzz_0, g_yyyyz_0_yyyyyyy_0, g_yyyyz_0_yyyyyyz_0, g_yyyyz_0_yyyyyzz_0, g_yyyyz_0_yyyyzzz_0, g_yyyyz_0_yyyzzzz_0, g_yyyyz_0_yyzzzzz_0, g_yyyyz_0_yzzzzzz_0, g_yyyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyz_0_xxxxxxx_0[i] = g_yyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxxy_0[i] = g_yyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxxz_0[i] = g_yyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxyy_0[i] = g_yyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxyz_0[i] = g_yyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxzz_0[i] = 2.0 * g_yyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyyy_0[i] = g_yyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyyz_0[i] = g_yyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyzz_0[i] = 2.0 * g_yyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxzzz_0[i] = 3.0 * g_yyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyyy_0[i] = g_yyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyyz_0[i] = g_yyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyzz_0[i] = 2.0 * g_yyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyzzz_0[i] = 3.0 * g_yyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxzzzz_0[i] = 4.0 * g_yyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyyy_0[i] = g_yyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyyz_0[i] = g_yyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyzz_0[i] = 2.0 * g_yyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyzzz_0[i] = 3.0 * g_yyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyzzzz_0[i] = 4.0 * g_yyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxzzzzz_0[i] = 5.0 * g_yyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyyy_0[i] = g_yyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyyz_0[i] = g_yyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyzz_0[i] = 2.0 * g_yyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyzzz_0[i] = 3.0 * g_yyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyzzzz_0[i] = 4.0 * g_yyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyzzzzz_0[i] = 5.0 * g_yyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xzzzzzz_0[i] = 6.0 * g_yyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xzzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyyy_0[i] = g_yyyy_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyyz_0[i] = g_yyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyy_0_yyyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyzz_0[i] = 2.0 * g_yyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyzzz_0[i] = 3.0 * g_yyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyzzzz_0[i] = 4.0 * g_yyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyzzzzz_0[i] = 5.0 * g_yyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yzzzzzz_0[i] = 6.0 * g_yyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyy_0_yzzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_zzzzzzz_0[i] = 7.0 * g_yyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 612-648 components of targeted buffer : HSK

    auto g_yyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 612);

    auto g_yyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 613);

    auto g_yyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 614);

    auto g_yyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 615);

    auto g_yyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 616);

    auto g_yyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 617);

    auto g_yyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 618);

    auto g_yyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 619);

    auto g_yyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 620);

    auto g_yyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 621);

    auto g_yyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 622);

    auto g_yyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 623);

    auto g_yyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 624);

    auto g_yyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 625);

    auto g_yyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 626);

    auto g_yyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 627);

    auto g_yyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 628);

    auto g_yyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 629);

    auto g_yyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 630);

    auto g_yyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 631);

    auto g_yyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 632);

    auto g_yyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 633);

    auto g_yyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 634);

    auto g_yyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 635);

    auto g_yyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 636);

    auto g_yyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 637);

    auto g_yyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 638);

    auto g_yyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 639);

    auto g_yyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 640);

    auto g_yyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 641);

    auto g_yyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 642);

    auto g_yyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 643);

    auto g_yyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 644);

    auto g_yyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 645);

    auto g_yyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 646);

    auto g_yyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 647);

    #pragma omp simd aligned(g_yyy_0_xxxxxxy_0, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxyy_0, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxyyy_0, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxyyyy_0, g_yyy_0_xxxyyyy_1, g_yyy_0_xxyyyyy_0, g_yyy_0_xxyyyyy_1, g_yyy_0_xyyyyyy_0, g_yyy_0_xyyyyyy_1, g_yyy_0_yyyyyyy_0, g_yyy_0_yyyyyyy_1, g_yyyz_0_xxxxxxy_1, g_yyyz_0_xxxxxyy_1, g_yyyz_0_xxxxyyy_1, g_yyyz_0_xxxyyyy_1, g_yyyz_0_xxyyyyy_1, g_yyyz_0_xyyyyyy_1, g_yyyz_0_yyyyyyy_1, g_yyyzz_0_xxxxxxx_0, g_yyyzz_0_xxxxxxy_0, g_yyyzz_0_xxxxxxz_0, g_yyyzz_0_xxxxxyy_0, g_yyyzz_0_xxxxxyz_0, g_yyyzz_0_xxxxxzz_0, g_yyyzz_0_xxxxyyy_0, g_yyyzz_0_xxxxyyz_0, g_yyyzz_0_xxxxyzz_0, g_yyyzz_0_xxxxzzz_0, g_yyyzz_0_xxxyyyy_0, g_yyyzz_0_xxxyyyz_0, g_yyyzz_0_xxxyyzz_0, g_yyyzz_0_xxxyzzz_0, g_yyyzz_0_xxxzzzz_0, g_yyyzz_0_xxyyyyy_0, g_yyyzz_0_xxyyyyz_0, g_yyyzz_0_xxyyyzz_0, g_yyyzz_0_xxyyzzz_0, g_yyyzz_0_xxyzzzz_0, g_yyyzz_0_xxzzzzz_0, g_yyyzz_0_xyyyyyy_0, g_yyyzz_0_xyyyyyz_0, g_yyyzz_0_xyyyyzz_0, g_yyyzz_0_xyyyzzz_0, g_yyyzz_0_xyyzzzz_0, g_yyyzz_0_xyzzzzz_0, g_yyyzz_0_xzzzzzz_0, g_yyyzz_0_yyyyyyy_0, g_yyyzz_0_yyyyyyz_0, g_yyyzz_0_yyyyyzz_0, g_yyyzz_0_yyyyzzz_0, g_yyyzz_0_yyyzzzz_0, g_yyyzz_0_yyzzzzz_0, g_yyyzz_0_yzzzzzz_0, g_yyyzz_0_zzzzzzz_0, g_yyzz_0_xxxxxxx_1, g_yyzz_0_xxxxxxz_1, g_yyzz_0_xxxxxyz_1, g_yyzz_0_xxxxxz_1, g_yyzz_0_xxxxxzz_1, g_yyzz_0_xxxxyyz_1, g_yyzz_0_xxxxyz_1, g_yyzz_0_xxxxyzz_1, g_yyzz_0_xxxxzz_1, g_yyzz_0_xxxxzzz_1, g_yyzz_0_xxxyyyz_1, g_yyzz_0_xxxyyz_1, g_yyzz_0_xxxyyzz_1, g_yyzz_0_xxxyzz_1, g_yyzz_0_xxxyzzz_1, g_yyzz_0_xxxzzz_1, g_yyzz_0_xxxzzzz_1, g_yyzz_0_xxyyyyz_1, g_yyzz_0_xxyyyz_1, g_yyzz_0_xxyyyzz_1, g_yyzz_0_xxyyzz_1, g_yyzz_0_xxyyzzz_1, g_yyzz_0_xxyzzz_1, g_yyzz_0_xxyzzzz_1, g_yyzz_0_xxzzzz_1, g_yyzz_0_xxzzzzz_1, g_yyzz_0_xyyyyyz_1, g_yyzz_0_xyyyyz_1, g_yyzz_0_xyyyyzz_1, g_yyzz_0_xyyyzz_1, g_yyzz_0_xyyyzzz_1, g_yyzz_0_xyyzzz_1, g_yyzz_0_xyyzzzz_1, g_yyzz_0_xyzzzz_1, g_yyzz_0_xyzzzzz_1, g_yyzz_0_xzzzzz_1, g_yyzz_0_xzzzzzz_1, g_yyzz_0_yyyyyyz_1, g_yyzz_0_yyyyyz_1, g_yyzz_0_yyyyyzz_1, g_yyzz_0_yyyyzz_1, g_yyzz_0_yyyyzzz_1, g_yyzz_0_yyyzzz_1, g_yyzz_0_yyyzzzz_1, g_yyzz_0_yyzzzz_1, g_yyzz_0_yyzzzzz_1, g_yyzz_0_yzzzzz_1, g_yyzz_0_yzzzzzz_1, g_yyzz_0_zzzzzz_1, g_yyzz_0_zzzzzzz_1, g_yzz_0_xxxxxxx_0, g_yzz_0_xxxxxxx_1, g_yzz_0_xxxxxxz_0, g_yzz_0_xxxxxxz_1, g_yzz_0_xxxxxyz_0, g_yzz_0_xxxxxyz_1, g_yzz_0_xxxxxzz_0, g_yzz_0_xxxxxzz_1, g_yzz_0_xxxxyyz_0, g_yzz_0_xxxxyyz_1, g_yzz_0_xxxxyzz_0, g_yzz_0_xxxxyzz_1, g_yzz_0_xxxxzzz_0, g_yzz_0_xxxxzzz_1, g_yzz_0_xxxyyyz_0, g_yzz_0_xxxyyyz_1, g_yzz_0_xxxyyzz_0, g_yzz_0_xxxyyzz_1, g_yzz_0_xxxyzzz_0, g_yzz_0_xxxyzzz_1, g_yzz_0_xxxzzzz_0, g_yzz_0_xxxzzzz_1, g_yzz_0_xxyyyyz_0, g_yzz_0_xxyyyyz_1, g_yzz_0_xxyyyzz_0, g_yzz_0_xxyyyzz_1, g_yzz_0_xxyyzzz_0, g_yzz_0_xxyyzzz_1, g_yzz_0_xxyzzzz_0, g_yzz_0_xxyzzzz_1, g_yzz_0_xxzzzzz_0, g_yzz_0_xxzzzzz_1, g_yzz_0_xyyyyyz_0, g_yzz_0_xyyyyyz_1, g_yzz_0_xyyyyzz_0, g_yzz_0_xyyyyzz_1, g_yzz_0_xyyyzzz_0, g_yzz_0_xyyyzzz_1, g_yzz_0_xyyzzzz_0, g_yzz_0_xyyzzzz_1, g_yzz_0_xyzzzzz_0, g_yzz_0_xyzzzzz_1, g_yzz_0_xzzzzzz_0, g_yzz_0_xzzzzzz_1, g_yzz_0_yyyyyyz_0, g_yzz_0_yyyyyyz_1, g_yzz_0_yyyyyzz_0, g_yzz_0_yyyyyzz_1, g_yzz_0_yyyyzzz_0, g_yzz_0_yyyyzzz_1, g_yzz_0_yyyzzzz_0, g_yzz_0_yyyzzzz_1, g_yzz_0_yyzzzzz_0, g_yzz_0_yyzzzzz_1, g_yzz_0_yzzzzzz_0, g_yzz_0_yzzzzzz_1, g_yzz_0_zzzzzzz_0, g_yzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzz_0_xxxxxxx_0[i] = 2.0 * g_yzz_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxxx_1[i] * fz_be_0 + g_yyzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxxy_0[i] = g_yyy_0_xxxxxxy_0[i] * fbe_0 - g_yyy_0_xxxxxxy_1[i] * fz_be_0 + g_yyyz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxxxz_0[i] = 2.0 * g_yzz_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxxz_1[i] * fz_be_0 + g_yyzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxyy_0[i] = g_yyy_0_xxxxxyy_0[i] * fbe_0 - g_yyy_0_xxxxxyy_1[i] * fz_be_0 + g_yyyz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxxyz_0[i] = 2.0 * g_yzz_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxyz_1[i] * fz_be_0 + g_yyzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxzz_0[i] = 2.0 * g_yzz_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxzz_1[i] * fz_be_0 + g_yyzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxyyy_0[i] = g_yyy_0_xxxxyyy_0[i] * fbe_0 - g_yyy_0_xxxxyyy_1[i] * fz_be_0 + g_yyyz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxyyz_0[i] = 2.0 * g_yzz_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxyzz_0[i] = 2.0 * g_yzz_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxyzz_1[i] * fz_be_0 + g_yyzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxzzz_0[i] = 2.0 * g_yzz_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxzzz_1[i] * fz_be_0 + g_yyzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyyyy_0[i] = g_yyy_0_xxxyyyy_0[i] * fbe_0 - g_yyy_0_xxxyyyy_1[i] * fz_be_0 + g_yyyz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxyyyz_0[i] = 2.0 * g_yzz_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyyzz_0[i] = 2.0 * g_yzz_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyzzz_0[i] = 2.0 * g_yzz_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyzzz_1[i] * fz_be_0 + g_yyzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxzzzz_0[i] = 2.0 * g_yzz_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxzzzz_1[i] * fz_be_0 + g_yyzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyyyy_0[i] = g_yyy_0_xxyyyyy_0[i] * fbe_0 - g_yyy_0_xxyyyyy_1[i] * fz_be_0 + g_yyyz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxyyyyz_0[i] = 2.0 * g_yzz_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyyzz_0[i] = 2.0 * g_yzz_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyzzz_0[i] = 2.0 * g_yzz_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyzzzz_0[i] = 2.0 * g_yzz_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyzzzz_1[i] * fz_be_0 + g_yyzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxzzzzz_0[i] = 2.0 * g_yzz_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxzzzzz_1[i] * fz_be_0 + g_yyzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyyyy_0[i] = g_yyy_0_xyyyyyy_0[i] * fbe_0 - g_yyy_0_xyyyyyy_1[i] * fz_be_0 + g_yyyz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xyyyyyz_0[i] = 2.0 * g_yzz_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyyzz_0[i] = 2.0 * g_yzz_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyzzz_0[i] = 2.0 * g_yzz_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyzzzz_0[i] = 2.0 * g_yzz_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyzzzzz_0[i] = 2.0 * g_yzz_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyzzzzz_1[i] * fz_be_0 + g_yyzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xzzzzzz_0[i] = 2.0 * g_yzz_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xzzzzzz_1[i] * fz_be_0 + g_yyzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyyyy_0[i] = g_yyy_0_yyyyyyy_0[i] * fbe_0 - g_yyy_0_yyyyyyy_1[i] * fz_be_0 + g_yyyz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_yyyyyyz_0[i] = 2.0 * g_yzz_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyyzz_0[i] = 2.0 * g_yzz_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyzzz_0[i] = 2.0 * g_yzz_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyzzzz_0[i] = 2.0 * g_yzz_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyzzzzz_0[i] = 2.0 * g_yzz_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yzzzzzz_0[i] = 2.0 * g_yzz_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yzzzzzz_1[i] * fz_be_0 + g_yyzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_zzzzzzz_0[i] = 2.0 * g_yzz_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_zzzzzzz_1[i] * fz_be_0 + g_yyzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 648-684 components of targeted buffer : HSK

    auto g_yyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 648);

    auto g_yyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 649);

    auto g_yyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 650);

    auto g_yyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 651);

    auto g_yyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 652);

    auto g_yyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 653);

    auto g_yyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 654);

    auto g_yyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 655);

    auto g_yyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 656);

    auto g_yyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 657);

    auto g_yyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 658);

    auto g_yyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 659);

    auto g_yyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 660);

    auto g_yyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 661);

    auto g_yyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 662);

    auto g_yyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 663);

    auto g_yyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 664);

    auto g_yyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 665);

    auto g_yyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 666);

    auto g_yyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 667);

    auto g_yyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 668);

    auto g_yyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 669);

    auto g_yyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 670);

    auto g_yyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 671);

    auto g_yyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 672);

    auto g_yyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 673);

    auto g_yyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 674);

    auto g_yyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 675);

    auto g_yyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 676);

    auto g_yyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 677);

    auto g_yyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 678);

    auto g_yyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 679);

    auto g_yyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 680);

    auto g_yyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 681);

    auto g_yyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 682);

    auto g_yyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 683);

    #pragma omp simd aligned(g_yyz_0_xxxxxxy_0, g_yyz_0_xxxxxxy_1, g_yyz_0_xxxxxyy_0, g_yyz_0_xxxxxyy_1, g_yyz_0_xxxxyyy_0, g_yyz_0_xxxxyyy_1, g_yyz_0_xxxyyyy_0, g_yyz_0_xxxyyyy_1, g_yyz_0_xxyyyyy_0, g_yyz_0_xxyyyyy_1, g_yyz_0_xyyyyyy_0, g_yyz_0_xyyyyyy_1, g_yyz_0_yyyyyyy_0, g_yyz_0_yyyyyyy_1, g_yyzz_0_xxxxxxy_1, g_yyzz_0_xxxxxyy_1, g_yyzz_0_xxxxyyy_1, g_yyzz_0_xxxyyyy_1, g_yyzz_0_xxyyyyy_1, g_yyzz_0_xyyyyyy_1, g_yyzz_0_yyyyyyy_1, g_yyzzz_0_xxxxxxx_0, g_yyzzz_0_xxxxxxy_0, g_yyzzz_0_xxxxxxz_0, g_yyzzz_0_xxxxxyy_0, g_yyzzz_0_xxxxxyz_0, g_yyzzz_0_xxxxxzz_0, g_yyzzz_0_xxxxyyy_0, g_yyzzz_0_xxxxyyz_0, g_yyzzz_0_xxxxyzz_0, g_yyzzz_0_xxxxzzz_0, g_yyzzz_0_xxxyyyy_0, g_yyzzz_0_xxxyyyz_0, g_yyzzz_0_xxxyyzz_0, g_yyzzz_0_xxxyzzz_0, g_yyzzz_0_xxxzzzz_0, g_yyzzz_0_xxyyyyy_0, g_yyzzz_0_xxyyyyz_0, g_yyzzz_0_xxyyyzz_0, g_yyzzz_0_xxyyzzz_0, g_yyzzz_0_xxyzzzz_0, g_yyzzz_0_xxzzzzz_0, g_yyzzz_0_xyyyyyy_0, g_yyzzz_0_xyyyyyz_0, g_yyzzz_0_xyyyyzz_0, g_yyzzz_0_xyyyzzz_0, g_yyzzz_0_xyyzzzz_0, g_yyzzz_0_xyzzzzz_0, g_yyzzz_0_xzzzzzz_0, g_yyzzz_0_yyyyyyy_0, g_yyzzz_0_yyyyyyz_0, g_yyzzz_0_yyyyyzz_0, g_yyzzz_0_yyyyzzz_0, g_yyzzz_0_yyyzzzz_0, g_yyzzz_0_yyzzzzz_0, g_yyzzz_0_yzzzzzz_0, g_yyzzz_0_zzzzzzz_0, g_yzzz_0_xxxxxxx_1, g_yzzz_0_xxxxxxz_1, g_yzzz_0_xxxxxyz_1, g_yzzz_0_xxxxxz_1, g_yzzz_0_xxxxxzz_1, g_yzzz_0_xxxxyyz_1, g_yzzz_0_xxxxyz_1, g_yzzz_0_xxxxyzz_1, g_yzzz_0_xxxxzz_1, g_yzzz_0_xxxxzzz_1, g_yzzz_0_xxxyyyz_1, g_yzzz_0_xxxyyz_1, g_yzzz_0_xxxyyzz_1, g_yzzz_0_xxxyzz_1, g_yzzz_0_xxxyzzz_1, g_yzzz_0_xxxzzz_1, g_yzzz_0_xxxzzzz_1, g_yzzz_0_xxyyyyz_1, g_yzzz_0_xxyyyz_1, g_yzzz_0_xxyyyzz_1, g_yzzz_0_xxyyzz_1, g_yzzz_0_xxyyzzz_1, g_yzzz_0_xxyzzz_1, g_yzzz_0_xxyzzzz_1, g_yzzz_0_xxzzzz_1, g_yzzz_0_xxzzzzz_1, g_yzzz_0_xyyyyyz_1, g_yzzz_0_xyyyyz_1, g_yzzz_0_xyyyyzz_1, g_yzzz_0_xyyyzz_1, g_yzzz_0_xyyyzzz_1, g_yzzz_0_xyyzzz_1, g_yzzz_0_xyyzzzz_1, g_yzzz_0_xyzzzz_1, g_yzzz_0_xyzzzzz_1, g_yzzz_0_xzzzzz_1, g_yzzz_0_xzzzzzz_1, g_yzzz_0_yyyyyyz_1, g_yzzz_0_yyyyyz_1, g_yzzz_0_yyyyyzz_1, g_yzzz_0_yyyyzz_1, g_yzzz_0_yyyyzzz_1, g_yzzz_0_yyyzzz_1, g_yzzz_0_yyyzzzz_1, g_yzzz_0_yyzzzz_1, g_yzzz_0_yyzzzzz_1, g_yzzz_0_yzzzzz_1, g_yzzz_0_yzzzzzz_1, g_yzzz_0_zzzzzz_1, g_yzzz_0_zzzzzzz_1, g_zzz_0_xxxxxxx_0, g_zzz_0_xxxxxxx_1, g_zzz_0_xxxxxxz_0, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxyz_0, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxzz_0, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxyyz_0, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyzz_0, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxzzz_0, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxyyyz_0, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyzz_0, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyzzz_0, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxzzzz_0, g_zzz_0_xxxzzzz_1, g_zzz_0_xxyyyyz_0, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyzz_0, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyzzz_0, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyzzzz_0, g_zzz_0_xxyzzzz_1, g_zzz_0_xxzzzzz_0, g_zzz_0_xxzzzzz_1, g_zzz_0_xyyyyyz_0, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyzz_0, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyzzz_0, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyzzzz_0, g_zzz_0_xyyzzzz_1, g_zzz_0_xyzzzzz_0, g_zzz_0_xyzzzzz_1, g_zzz_0_xzzzzzz_0, g_zzz_0_xzzzzzz_1, g_zzz_0_yyyyyyz_0, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyzz_0, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyzzz_0, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyzzzz_0, g_zzz_0_yyyzzzz_1, g_zzz_0_yyzzzzz_0, g_zzz_0_yyzzzzz_1, g_zzz_0_yzzzzzz_0, g_zzz_0_yzzzzzz_1, g_zzz_0_zzzzzzz_0, g_zzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzz_0_xxxxxxx_0[i] = g_zzz_0_xxxxxxx_0[i] * fbe_0 - g_zzz_0_xxxxxxx_1[i] * fz_be_0 + g_yzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxxy_0[i] = 2.0 * g_yyz_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxxxy_1[i] * fz_be_0 + g_yyzz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxxxz_0[i] = g_zzz_0_xxxxxxz_0[i] * fbe_0 - g_zzz_0_xxxxxxz_1[i] * fz_be_0 + g_yzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxyy_0[i] = 2.0 * g_yyz_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxxyy_1[i] * fz_be_0 + g_yyzz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxxyz_0[i] = g_zzz_0_xxxxxyz_0[i] * fbe_0 - g_zzz_0_xxxxxyz_1[i] * fz_be_0 + g_yzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxzz_0[i] = g_zzz_0_xxxxxzz_0[i] * fbe_0 - g_zzz_0_xxxxxzz_1[i] * fz_be_0 + g_yzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxyyy_0[i] = 2.0 * g_yyz_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxyyy_1[i] * fz_be_0 + g_yyzz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxyyz_0[i] = g_zzz_0_xxxxyyz_0[i] * fbe_0 - g_zzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxyzz_0[i] = g_zzz_0_xxxxyzz_0[i] * fbe_0 - g_zzz_0_xxxxyzz_1[i] * fz_be_0 + g_yzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxzzz_0[i] = g_zzz_0_xxxxzzz_0[i] * fbe_0 - g_zzz_0_xxxxzzz_1[i] * fz_be_0 + g_yzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyyyy_0[i] = 2.0 * g_yyz_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxyyyy_1[i] * fz_be_0 + g_yyzz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxyyyz_0[i] = g_zzz_0_xxxyyyz_0[i] * fbe_0 - g_zzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyyzz_0[i] = g_zzz_0_xxxyyzz_0[i] * fbe_0 - g_zzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyzzz_0[i] = g_zzz_0_xxxyzzz_0[i] * fbe_0 - g_zzz_0_xxxyzzz_1[i] * fz_be_0 + g_yzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxzzzz_0[i] = g_zzz_0_xxxzzzz_0[i] * fbe_0 - g_zzz_0_xxxzzzz_1[i] * fz_be_0 + g_yzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyyyy_0[i] = 2.0 * g_yyz_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxyyyyy_1[i] * fz_be_0 + g_yyzz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxyyyyz_0[i] = g_zzz_0_xxyyyyz_0[i] * fbe_0 - g_zzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyyzz_0[i] = g_zzz_0_xxyyyzz_0[i] * fbe_0 - g_zzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyzzz_0[i] = g_zzz_0_xxyyzzz_0[i] * fbe_0 - g_zzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyzzzz_0[i] = g_zzz_0_xxyzzzz_0[i] * fbe_0 - g_zzz_0_xxyzzzz_1[i] * fz_be_0 + g_yzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxzzzzz_0[i] = g_zzz_0_xxzzzzz_0[i] * fbe_0 - g_zzz_0_xxzzzzz_1[i] * fz_be_0 + g_yzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyyyy_0[i] = 2.0 * g_yyz_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xyyyyyy_1[i] * fz_be_0 + g_yyzz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xyyyyyz_0[i] = g_zzz_0_xyyyyyz_0[i] * fbe_0 - g_zzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyyzz_0[i] = g_zzz_0_xyyyyzz_0[i] * fbe_0 - g_zzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyzzz_0[i] = g_zzz_0_xyyyzzz_0[i] * fbe_0 - g_zzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyzzzz_0[i] = g_zzz_0_xyyzzzz_0[i] * fbe_0 - g_zzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyzzzzz_0[i] = g_zzz_0_xyzzzzz_0[i] * fbe_0 - g_zzz_0_xyzzzzz_1[i] * fz_be_0 + g_yzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xzzzzzz_0[i] = g_zzz_0_xzzzzzz_0[i] * fbe_0 - g_zzz_0_xzzzzzz_1[i] * fz_be_0 + g_yzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyyyy_0[i] = 2.0 * g_yyz_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_yyyyyyy_1[i] * fz_be_0 + g_yyzz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_yyyyyyz_0[i] = g_zzz_0_yyyyyyz_0[i] * fbe_0 - g_zzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyyzz_0[i] = g_zzz_0_yyyyyzz_0[i] * fbe_0 - g_zzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyzzz_0[i] = g_zzz_0_yyyyzzz_0[i] * fbe_0 - g_zzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyzzzz_0[i] = g_zzz_0_yyyzzzz_0[i] * fbe_0 - g_zzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyzzzzz_0[i] = g_zzz_0_yyzzzzz_0[i] * fbe_0 - g_zzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yzzzzzz_0[i] = g_zzz_0_yzzzzzz_0[i] * fbe_0 - g_zzz_0_yzzzzzz_1[i] * fz_be_0 + g_yzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_zzzzzzz_0[i] = g_zzz_0_zzzzzzz_0[i] * fbe_0 - g_zzz_0_zzzzzzz_1[i] * fz_be_0 + g_yzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 684-720 components of targeted buffer : HSK

    auto g_yzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 684);

    auto g_yzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 685);

    auto g_yzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 686);

    auto g_yzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 687);

    auto g_yzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 688);

    auto g_yzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 689);

    auto g_yzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 690);

    auto g_yzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 691);

    auto g_yzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 692);

    auto g_yzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 693);

    auto g_yzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 694);

    auto g_yzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 695);

    auto g_yzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 696);

    auto g_yzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 697);

    auto g_yzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 698);

    auto g_yzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 699);

    auto g_yzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 700);

    auto g_yzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 701);

    auto g_yzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 702);

    auto g_yzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 703);

    auto g_yzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 704);

    auto g_yzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 705);

    auto g_yzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 706);

    auto g_yzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 707);

    auto g_yzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 708);

    auto g_yzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 709);

    auto g_yzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 710);

    auto g_yzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 711);

    auto g_yzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 712);

    auto g_yzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 713);

    auto g_yzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 714);

    auto g_yzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 715);

    auto g_yzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 716);

    auto g_yzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 717);

    auto g_yzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 718);

    auto g_yzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 719);

    #pragma omp simd aligned(g_yzzzz_0_xxxxxxx_0, g_yzzzz_0_xxxxxxy_0, g_yzzzz_0_xxxxxxz_0, g_yzzzz_0_xxxxxyy_0, g_yzzzz_0_xxxxxyz_0, g_yzzzz_0_xxxxxzz_0, g_yzzzz_0_xxxxyyy_0, g_yzzzz_0_xxxxyyz_0, g_yzzzz_0_xxxxyzz_0, g_yzzzz_0_xxxxzzz_0, g_yzzzz_0_xxxyyyy_0, g_yzzzz_0_xxxyyyz_0, g_yzzzz_0_xxxyyzz_0, g_yzzzz_0_xxxyzzz_0, g_yzzzz_0_xxxzzzz_0, g_yzzzz_0_xxyyyyy_0, g_yzzzz_0_xxyyyyz_0, g_yzzzz_0_xxyyyzz_0, g_yzzzz_0_xxyyzzz_0, g_yzzzz_0_xxyzzzz_0, g_yzzzz_0_xxzzzzz_0, g_yzzzz_0_xyyyyyy_0, g_yzzzz_0_xyyyyyz_0, g_yzzzz_0_xyyyyzz_0, g_yzzzz_0_xyyyzzz_0, g_yzzzz_0_xyyzzzz_0, g_yzzzz_0_xyzzzzz_0, g_yzzzz_0_xzzzzzz_0, g_yzzzz_0_yyyyyyy_0, g_yzzzz_0_yyyyyyz_0, g_yzzzz_0_yyyyyzz_0, g_yzzzz_0_yyyyzzz_0, g_yzzzz_0_yyyzzzz_0, g_yzzzz_0_yyzzzzz_0, g_yzzzz_0_yzzzzzz_0, g_yzzzz_0_zzzzzzz_0, g_zzzz_0_xxxxxx_1, g_zzzz_0_xxxxxxx_1, g_zzzz_0_xxxxxxy_1, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxy_1, g_zzzz_0_xxxxxyy_1, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxyy_1, g_zzzz_0_xxxxyyy_1, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxyyy_1, g_zzzz_0_xxxyyyy_1, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxyyyy_1, g_zzzz_0_xxyyyyy_1, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxzzzz_1, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xyyyyy_1, g_zzzz_0_xyyyyyy_1, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyzzzz_1, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xzzzzz_1, g_zzzz_0_xzzzzzz_1, g_zzzz_0_yyyyyy_1, g_zzzz_0_yyyyyyy_1, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyzzzz_1, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yzzzzz_1, g_zzzz_0_yzzzzzz_1, g_zzzz_0_zzzzzz_1, g_zzzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzz_0_xxxxxxx_0[i] = g_zzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxxy_0[i] = g_zzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxxz_0[i] = g_zzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxyy_0[i] = 2.0 * g_zzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxyz_0[i] = g_zzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxzz_0[i] = g_zzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyyy_0[i] = 3.0 * g_zzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyyz_0[i] = 2.0 * g_zzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyzz_0[i] = g_zzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxzzz_0[i] = g_zzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyyy_0[i] = 4.0 * g_zzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyyz_0[i] = 3.0 * g_zzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyzz_0[i] = 2.0 * g_zzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyzzz_0[i] = g_zzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxzzzz_0[i] = g_zzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyyy_0[i] = 5.0 * g_zzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyyz_0[i] = 4.0 * g_zzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyzz_0[i] = 3.0 * g_zzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyzzz_0[i] = 2.0 * g_zzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyzzzz_0[i] = g_zzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxzzzzz_0[i] = g_zzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyyy_0[i] = 6.0 * g_zzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyyz_0[i] = 5.0 * g_zzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyzz_0[i] = 4.0 * g_zzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyzzz_0[i] = 3.0 * g_zzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyzzzz_0[i] = 2.0 * g_zzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyzzzzz_0[i] = g_zzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xzzzzzz_0[i] = g_zzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyyy_0[i] = 7.0 * g_zzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzz_0_yyyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyyz_0[i] = 6.0 * g_zzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyzz_0[i] = 5.0 * g_zzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyzzz_0[i] = 4.0 * g_zzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyzzzz_0[i] = 3.0 * g_zzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyzzzzz_0[i] = 2.0 * g_zzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yzzzzzz_0[i] = g_zzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_zzzzzzz_0[i] = g_zzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 720-756 components of targeted buffer : HSK

    auto g_zzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 720);

    auto g_zzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 721);

    auto g_zzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 722);

    auto g_zzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 723);

    auto g_zzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 724);

    auto g_zzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 725);

    auto g_zzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 726);

    auto g_zzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 727);

    auto g_zzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 728);

    auto g_zzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 729);

    auto g_zzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 730);

    auto g_zzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 731);

    auto g_zzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 732);

    auto g_zzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 733);

    auto g_zzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 734);

    auto g_zzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 735);

    auto g_zzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 736);

    auto g_zzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 737);

    auto g_zzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 738);

    auto g_zzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 739);

    auto g_zzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 740);

    auto g_zzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 741);

    auto g_zzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 742);

    auto g_zzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 743);

    auto g_zzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 744);

    auto g_zzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 745);

    auto g_zzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 746);

    auto g_zzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 747);

    auto g_zzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 748);

    auto g_zzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 749);

    auto g_zzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 750);

    auto g_zzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 751);

    auto g_zzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 752);

    auto g_zzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 753);

    auto g_zzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 754);

    auto g_zzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 755);

    #pragma omp simd aligned(g_zzz_0_xxxxxxx_0, g_zzz_0_xxxxxxx_1, g_zzz_0_xxxxxxy_0, g_zzz_0_xxxxxxy_1, g_zzz_0_xxxxxxz_0, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxyy_0, g_zzz_0_xxxxxyy_1, g_zzz_0_xxxxxyz_0, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxzz_0, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxyyy_0, g_zzz_0_xxxxyyy_1, g_zzz_0_xxxxyyz_0, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyzz_0, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxzzz_0, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxyyyy_0, g_zzz_0_xxxyyyy_1, g_zzz_0_xxxyyyz_0, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyzz_0, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyzzz_0, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxzzzz_0, g_zzz_0_xxxzzzz_1, g_zzz_0_xxyyyyy_0, g_zzz_0_xxyyyyy_1, g_zzz_0_xxyyyyz_0, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyzz_0, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyzzz_0, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyzzzz_0, g_zzz_0_xxyzzzz_1, g_zzz_0_xxzzzzz_0, g_zzz_0_xxzzzzz_1, g_zzz_0_xyyyyyy_0, g_zzz_0_xyyyyyy_1, g_zzz_0_xyyyyyz_0, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyzz_0, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyzzz_0, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyzzzz_0, g_zzz_0_xyyzzzz_1, g_zzz_0_xyzzzzz_0, g_zzz_0_xyzzzzz_1, g_zzz_0_xzzzzzz_0, g_zzz_0_xzzzzzz_1, g_zzz_0_yyyyyyy_0, g_zzz_0_yyyyyyy_1, g_zzz_0_yyyyyyz_0, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyzz_0, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyzzz_0, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyzzzz_0, g_zzz_0_yyyzzzz_1, g_zzz_0_yyzzzzz_0, g_zzz_0_yyzzzzz_1, g_zzz_0_yzzzzzz_0, g_zzz_0_yzzzzzz_1, g_zzz_0_zzzzzzz_0, g_zzz_0_zzzzzzz_1, g_zzzz_0_xxxxxx_1, g_zzzz_0_xxxxxxx_1, g_zzzz_0_xxxxxxy_1, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxy_1, g_zzzz_0_xxxxxyy_1, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxyy_1, g_zzzz_0_xxxxyyy_1, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxyyy_1, g_zzzz_0_xxxyyyy_1, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxyyyy_1, g_zzzz_0_xxyyyyy_1, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxzzzz_1, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xyyyyy_1, g_zzzz_0_xyyyyyy_1, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyzzzz_1, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xzzzzz_1, g_zzzz_0_xzzzzzz_1, g_zzzz_0_yyyyyy_1, g_zzzz_0_yyyyyyy_1, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyzzzz_1, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yzzzzz_1, g_zzzz_0_yzzzzzz_1, g_zzzz_0_zzzzzz_1, g_zzzz_0_zzzzzzz_1, g_zzzzz_0_xxxxxxx_0, g_zzzzz_0_xxxxxxy_0, g_zzzzz_0_xxxxxxz_0, g_zzzzz_0_xxxxxyy_0, g_zzzzz_0_xxxxxyz_0, g_zzzzz_0_xxxxxzz_0, g_zzzzz_0_xxxxyyy_0, g_zzzzz_0_xxxxyyz_0, g_zzzzz_0_xxxxyzz_0, g_zzzzz_0_xxxxzzz_0, g_zzzzz_0_xxxyyyy_0, g_zzzzz_0_xxxyyyz_0, g_zzzzz_0_xxxyyzz_0, g_zzzzz_0_xxxyzzz_0, g_zzzzz_0_xxxzzzz_0, g_zzzzz_0_xxyyyyy_0, g_zzzzz_0_xxyyyyz_0, g_zzzzz_0_xxyyyzz_0, g_zzzzz_0_xxyyzzz_0, g_zzzzz_0_xxyzzzz_0, g_zzzzz_0_xxzzzzz_0, g_zzzzz_0_xyyyyyy_0, g_zzzzz_0_xyyyyyz_0, g_zzzzz_0_xyyyyzz_0, g_zzzzz_0_xyyyzzz_0, g_zzzzz_0_xyyzzzz_0, g_zzzzz_0_xyzzzzz_0, g_zzzzz_0_xzzzzzz_0, g_zzzzz_0_yyyyyyy_0, g_zzzzz_0_yyyyyyz_0, g_zzzzz_0_yyyyyzz_0, g_zzzzz_0_yyyyzzz_0, g_zzzzz_0_yyyzzzz_0, g_zzzzz_0_yyzzzzz_0, g_zzzzz_0_yzzzzzz_0, g_zzzzz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzz_0_xxxxxxx_0[i] = 4.0 * g_zzz_0_xxxxxxx_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxx_1[i] * fz_be_0 + g_zzzz_0_xxxxxxx_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxxy_0[i] = 4.0 * g_zzz_0_xxxxxxy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxy_1[i] * fz_be_0 + g_zzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxxz_0[i] = 4.0 * g_zzz_0_xxxxxxz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxz_1[i] * fz_be_0 + g_zzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxyy_0[i] = 4.0 * g_zzz_0_xxxxxyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxyy_1[i] * fz_be_0 + g_zzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxyz_0[i] = 4.0 * g_zzz_0_xxxxxyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxyz_1[i] * fz_be_0 + g_zzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxzz_0[i] = 4.0 * g_zzz_0_xxxxxzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyyy_0[i] = 4.0 * g_zzz_0_xxxxyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyyy_1[i] * fz_be_0 + g_zzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyyz_0[i] = 4.0 * g_zzz_0_xxxxyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyyz_1[i] * fz_be_0 + g_zzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyzz_0[i] = 4.0 * g_zzz_0_xxxxyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxzzz_0[i] = 4.0 * g_zzz_0_xxxxzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyyy_0[i] = 4.0 * g_zzz_0_xxxyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyyy_1[i] * fz_be_0 + g_zzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyyz_0[i] = 4.0 * g_zzz_0_xxxyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyyz_1[i] * fz_be_0 + g_zzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyzz_0[i] = 4.0 * g_zzz_0_xxxyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyzzz_0[i] = 4.0 * g_zzz_0_xxxyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxzzzz_0[i] = 4.0 * g_zzz_0_xxxzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyyy_0[i] = 4.0 * g_zzz_0_xxyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyyy_1[i] * fz_be_0 + g_zzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyyz_0[i] = 4.0 * g_zzz_0_xxyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyyz_1[i] * fz_be_0 + g_zzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyzz_0[i] = 4.0 * g_zzz_0_xxyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyzzz_0[i] = 4.0 * g_zzz_0_xxyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyzzzz_0[i] = 4.0 * g_zzz_0_xxyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxzzzzz_0[i] = 4.0 * g_zzz_0_xxzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyyy_0[i] = 4.0 * g_zzz_0_xyyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyyy_1[i] * fz_be_0 + g_zzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyyz_0[i] = 4.0 * g_zzz_0_xyyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyyz_1[i] * fz_be_0 + g_zzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyzz_0[i] = 4.0 * g_zzz_0_xyyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyzzz_0[i] = 4.0 * g_zzz_0_xyyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyzzzz_0[i] = 4.0 * g_zzz_0_xyyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyzzzzz_0[i] = 4.0 * g_zzz_0_xyzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xzzzzzz_0[i] = 4.0 * g_zzz_0_xzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xzzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyyy_0[i] = 4.0 * g_zzz_0_yyyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyyy_1[i] * fz_be_0 + g_zzzz_0_yyyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyyz_0[i] = 4.0 * g_zzz_0_yyyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyyz_1[i] * fz_be_0 + g_zzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzz_0_yyyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyzz_0[i] = 4.0 * g_zzz_0_yyyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyzzz_0[i] = 4.0 * g_zzz_0_yyyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyzzzz_0[i] = 4.0 * g_zzz_0_yyyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyzzzzz_0[i] = 4.0 * g_zzz_0_yyzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yzzzzzz_0[i] = 4.0 * g_zzz_0_yzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzz_0_yzzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_zzzzzzz_0[i] = 4.0 * g_zzz_0_zzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_zzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzz_0_zzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

