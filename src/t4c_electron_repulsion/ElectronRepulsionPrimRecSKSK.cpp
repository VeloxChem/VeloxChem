#include "ElectronRepulsionPrimRecSKSK.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sksk(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sksk,
                                  size_t idx_eri_0_shsk,
                                  size_t idx_eri_1_shsk,
                                  size_t idx_eri_1_sisi,
                                  size_t idx_eri_0_sisk,
                                  size_t idx_eri_1_sisk,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void
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

    /// Set up components of auxilary buffer : SHSK

    auto g_0_xxxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk);

    auto g_0_xxxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 1);

    auto g_0_xxxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 2);

    auto g_0_xxxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 3);

    auto g_0_xxxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 4);

    auto g_0_xxxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 5);

    auto g_0_xxxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 6);

    auto g_0_xxxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 7);

    auto g_0_xxxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 8);

    auto g_0_xxxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 9);

    auto g_0_xxxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 10);

    auto g_0_xxxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 11);

    auto g_0_xxxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 12);

    auto g_0_xxxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 13);

    auto g_0_xxxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 14);

    auto g_0_xxxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 15);

    auto g_0_xxxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 16);

    auto g_0_xxxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 17);

    auto g_0_xxxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 18);

    auto g_0_xxxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 19);

    auto g_0_xxxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 20);

    auto g_0_xxxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 21);

    auto g_0_xxxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 22);

    auto g_0_xxxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 23);

    auto g_0_xxxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 24);

    auto g_0_xxxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 25);

    auto g_0_xxxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 26);

    auto g_0_xxxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 27);

    auto g_0_xxxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 28);

    auto g_0_xxxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 29);

    auto g_0_xxxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 30);

    auto g_0_xxxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 31);

    auto g_0_xxxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 32);

    auto g_0_xxxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 33);

    auto g_0_xxxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 34);

    auto g_0_xxxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 35);

    auto g_0_xxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 36);

    auto g_0_xxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 38);

    auto g_0_xxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 41);

    auto g_0_xxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 45);

    auto g_0_xxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 50);

    auto g_0_xxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 56);

    auto g_0_xxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 63);

    auto g_0_xxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 72);

    auto g_0_xxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 73);

    auto g_0_xxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 75);

    auto g_0_xxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 78);

    auto g_0_xxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 82);

    auto g_0_xxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 87);

    auto g_0_xxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 93);

    auto g_0_xxxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 108);

    auto g_0_xxxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 109);

    auto g_0_xxxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 110);

    auto g_0_xxxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 111);

    auto g_0_xxxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 112);

    auto g_0_xxxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 113);

    auto g_0_xxxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 114);

    auto g_0_xxxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 115);

    auto g_0_xxxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 116);

    auto g_0_xxxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 117);

    auto g_0_xxxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 118);

    auto g_0_xxxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 119);

    auto g_0_xxxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 120);

    auto g_0_xxxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 121);

    auto g_0_xxxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 122);

    auto g_0_xxxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 123);

    auto g_0_xxxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 124);

    auto g_0_xxxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 125);

    auto g_0_xxxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 126);

    auto g_0_xxxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 127);

    auto g_0_xxxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 128);

    auto g_0_xxxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 129);

    auto g_0_xxxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 130);

    auto g_0_xxxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 131);

    auto g_0_xxxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 132);

    auto g_0_xxxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 133);

    auto g_0_xxxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 134);

    auto g_0_xxxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 135);

    auto g_0_xxxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 136);

    auto g_0_xxxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 137);

    auto g_0_xxxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 138);

    auto g_0_xxxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 139);

    auto g_0_xxxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 140);

    auto g_0_xxxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 141);

    auto g_0_xxxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 142);

    auto g_0_xxxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 143);

    auto g_0_xxxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 180);

    auto g_0_xxxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 181);

    auto g_0_xxxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 182);

    auto g_0_xxxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 183);

    auto g_0_xxxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 184);

    auto g_0_xxxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 185);

    auto g_0_xxxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 186);

    auto g_0_xxxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 187);

    auto g_0_xxxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 188);

    auto g_0_xxxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 189);

    auto g_0_xxxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 190);

    auto g_0_xxxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 191);

    auto g_0_xxxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 192);

    auto g_0_xxxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 193);

    auto g_0_xxxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 194);

    auto g_0_xxxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 195);

    auto g_0_xxxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 196);

    auto g_0_xxxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 197);

    auto g_0_xxxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 198);

    auto g_0_xxxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 199);

    auto g_0_xxxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 200);

    auto g_0_xxxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 201);

    auto g_0_xxxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 202);

    auto g_0_xxxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 203);

    auto g_0_xxxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 204);

    auto g_0_xxxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 205);

    auto g_0_xxxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 206);

    auto g_0_xxxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 207);

    auto g_0_xxxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 208);

    auto g_0_xxxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 209);

    auto g_0_xxxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 210);

    auto g_0_xxxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 211);

    auto g_0_xxxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 212);

    auto g_0_xxxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 213);

    auto g_0_xxxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 214);

    auto g_0_xxxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 215);

    auto g_0_xxyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 216);

    auto g_0_xxyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 217);

    auto g_0_xxyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 218);

    auto g_0_xxyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 219);

    auto g_0_xxyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 220);

    auto g_0_xxyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 221);

    auto g_0_xxyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 222);

    auto g_0_xxyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 223);

    auto g_0_xxyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 224);

    auto g_0_xxyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 225);

    auto g_0_xxyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 226);

    auto g_0_xxyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 227);

    auto g_0_xxyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 228);

    auto g_0_xxyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 229);

    auto g_0_xxyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 230);

    auto g_0_xxyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 231);

    auto g_0_xxyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 232);

    auto g_0_xxyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 233);

    auto g_0_xxyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 234);

    auto g_0_xxyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 235);

    auto g_0_xxyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 236);

    auto g_0_xxyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 237);

    auto g_0_xxyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 238);

    auto g_0_xxyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 239);

    auto g_0_xxyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 240);

    auto g_0_xxyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 241);

    auto g_0_xxyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 242);

    auto g_0_xxyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 243);

    auto g_0_xxyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 244);

    auto g_0_xxyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 245);

    auto g_0_xxyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 246);

    auto g_0_xxyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 247);

    auto g_0_xxyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 248);

    auto g_0_xxyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 249);

    auto g_0_xxyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 250);

    auto g_0_xxyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 251);

    auto g_0_xxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 253);

    auto g_0_xxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 255);

    auto g_0_xxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 258);

    auto g_0_xxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 262);

    auto g_0_xxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 267);

    auto g_0_xxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 273);

    auto g_0_xxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 288);

    auto g_0_xxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 290);

    auto g_0_xxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 293);

    auto g_0_xxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 297);

    auto g_0_xxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 302);

    auto g_0_xxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 308);

    auto g_0_xxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 315);

    auto g_0_xxzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 324);

    auto g_0_xxzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 325);

    auto g_0_xxzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 326);

    auto g_0_xxzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 327);

    auto g_0_xxzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 328);

    auto g_0_xxzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 329);

    auto g_0_xxzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 330);

    auto g_0_xxzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 331);

    auto g_0_xxzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 332);

    auto g_0_xxzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 333);

    auto g_0_xxzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 334);

    auto g_0_xxzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 335);

    auto g_0_xxzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 336);

    auto g_0_xxzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 337);

    auto g_0_xxzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 338);

    auto g_0_xxzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 339);

    auto g_0_xxzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 340);

    auto g_0_xxzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 341);

    auto g_0_xxzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 342);

    auto g_0_xxzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 343);

    auto g_0_xxzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 344);

    auto g_0_xxzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 345);

    auto g_0_xxzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 346);

    auto g_0_xxzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 347);

    auto g_0_xxzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 348);

    auto g_0_xxzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 349);

    auto g_0_xxzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 350);

    auto g_0_xxzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 351);

    auto g_0_xxzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 352);

    auto g_0_xxzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 353);

    auto g_0_xxzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 354);

    auto g_0_xxzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 355);

    auto g_0_xxzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 356);

    auto g_0_xxzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 357);

    auto g_0_xxzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 358);

    auto g_0_xxzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 359);

    auto g_0_xyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 361);

    auto g_0_xyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 363);

    auto g_0_xyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 364);

    auto g_0_xyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 366);

    auto g_0_xyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 367);

    auto g_0_xyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 368);

    auto g_0_xyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 370);

    auto g_0_xyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 371);

    auto g_0_xyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 372);

    auto g_0_xyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 373);

    auto g_0_xyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 375);

    auto g_0_xyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 376);

    auto g_0_xyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 377);

    auto g_0_xyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 378);

    auto g_0_xyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 379);

    auto g_0_xyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 381);

    auto g_0_xyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 382);

    auto g_0_xyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 383);

    auto g_0_xyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 384);

    auto g_0_xyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 385);

    auto g_0_xyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 386);

    auto g_0_xyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 388);

    auto g_0_xyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 389);

    auto g_0_xyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 390);

    auto g_0_xyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 391);

    auto g_0_xyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 392);

    auto g_0_xyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 393);

    auto g_0_xyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 394);

    auto g_0_xyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 395);

    auto g_0_xyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 436);

    auto g_0_xyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 439);

    auto g_0_xyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 440);

    auto g_0_xyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 443);

    auto g_0_xyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 444);

    auto g_0_xyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 445);

    auto g_0_xyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 448);

    auto g_0_xyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 449);

    auto g_0_xyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 450);

    auto g_0_xyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 451);

    auto g_0_xyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 454);

    auto g_0_xyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 455);

    auto g_0_xyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 456);

    auto g_0_xyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 457);

    auto g_0_xyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 458);

    auto g_0_xyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 460);

    auto g_0_xyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 461);

    auto g_0_xyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 462);

    auto g_0_xyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 463);

    auto g_0_xyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 464);

    auto g_0_xyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 465);

    auto g_0_xyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 466);

    auto g_0_xyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 467);

    auto g_0_xzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 506);

    auto g_0_xzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 508);

    auto g_0_xzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 509);

    auto g_0_xzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 511);

    auto g_0_xzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 512);

    auto g_0_xzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 513);

    auto g_0_xzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 515);

    auto g_0_xzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 516);

    auto g_0_xzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 517);

    auto g_0_xzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 518);

    auto g_0_xzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 520);

    auto g_0_xzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 521);

    auto g_0_xzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 522);

    auto g_0_xzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 523);

    auto g_0_xzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 524);

    auto g_0_xzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 526);

    auto g_0_xzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 527);

    auto g_0_xzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 528);

    auto g_0_xzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 529);

    auto g_0_xzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 530);

    auto g_0_xzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 531);

    auto g_0_xzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 532);

    auto g_0_xzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 533);

    auto g_0_xzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 534);

    auto g_0_xzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 535);

    auto g_0_xzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 536);

    auto g_0_xzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 537);

    auto g_0_xzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 538);

    auto g_0_xzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 539);

    auto g_0_yyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 540);

    auto g_0_yyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 541);

    auto g_0_yyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 542);

    auto g_0_yyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 543);

    auto g_0_yyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 544);

    auto g_0_yyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 545);

    auto g_0_yyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 546);

    auto g_0_yyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 547);

    auto g_0_yyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 548);

    auto g_0_yyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 549);

    auto g_0_yyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 550);

    auto g_0_yyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 551);

    auto g_0_yyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 552);

    auto g_0_yyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 553);

    auto g_0_yyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 554);

    auto g_0_yyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 555);

    auto g_0_yyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 556);

    auto g_0_yyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 557);

    auto g_0_yyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 558);

    auto g_0_yyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 559);

    auto g_0_yyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 560);

    auto g_0_yyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 561);

    auto g_0_yyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 562);

    auto g_0_yyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 563);

    auto g_0_yyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 564);

    auto g_0_yyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 565);

    auto g_0_yyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 566);

    auto g_0_yyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 567);

    auto g_0_yyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 568);

    auto g_0_yyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 569);

    auto g_0_yyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 570);

    auto g_0_yyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 571);

    auto g_0_yyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 572);

    auto g_0_yyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 573);

    auto g_0_yyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 574);

    auto g_0_yyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 575);

    auto g_0_yyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 577);

    auto g_0_yyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 579);

    auto g_0_yyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 582);

    auto g_0_yyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 586);

    auto g_0_yyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 591);

    auto g_0_yyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 597);

    auto g_0_yyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 604);

    auto g_0_yyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 612);

    auto g_0_yyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 613);

    auto g_0_yyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 614);

    auto g_0_yyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 615);

    auto g_0_yyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 616);

    auto g_0_yyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 617);

    auto g_0_yyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 618);

    auto g_0_yyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 619);

    auto g_0_yyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 620);

    auto g_0_yyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 621);

    auto g_0_yyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 622);

    auto g_0_yyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 623);

    auto g_0_yyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 624);

    auto g_0_yyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 625);

    auto g_0_yyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 626);

    auto g_0_yyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 627);

    auto g_0_yyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 628);

    auto g_0_yyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 629);

    auto g_0_yyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 630);

    auto g_0_yyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 631);

    auto g_0_yyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 632);

    auto g_0_yyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 633);

    auto g_0_yyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 634);

    auto g_0_yyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 635);

    auto g_0_yyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 636);

    auto g_0_yyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 637);

    auto g_0_yyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 638);

    auto g_0_yyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 639);

    auto g_0_yyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 640);

    auto g_0_yyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 641);

    auto g_0_yyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 642);

    auto g_0_yyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 643);

    auto g_0_yyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 644);

    auto g_0_yyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 645);

    auto g_0_yyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 646);

    auto g_0_yyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 647);

    auto g_0_yyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 648);

    auto g_0_yyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 649);

    auto g_0_yyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 650);

    auto g_0_yyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 651);

    auto g_0_yyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 652);

    auto g_0_yyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 653);

    auto g_0_yyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 654);

    auto g_0_yyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 655);

    auto g_0_yyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 656);

    auto g_0_yyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 657);

    auto g_0_yyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 658);

    auto g_0_yyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 659);

    auto g_0_yyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 660);

    auto g_0_yyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 661);

    auto g_0_yyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 662);

    auto g_0_yyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 663);

    auto g_0_yyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 664);

    auto g_0_yyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 665);

    auto g_0_yyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 666);

    auto g_0_yyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 667);

    auto g_0_yyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 668);

    auto g_0_yyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 669);

    auto g_0_yyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 670);

    auto g_0_yyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 671);

    auto g_0_yyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 672);

    auto g_0_yyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 673);

    auto g_0_yyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 674);

    auto g_0_yyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 675);

    auto g_0_yyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 676);

    auto g_0_yyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 677);

    auto g_0_yyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 678);

    auto g_0_yyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 679);

    auto g_0_yyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 680);

    auto g_0_yyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 681);

    auto g_0_yyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 682);

    auto g_0_yyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 683);

    auto g_0_yzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 684);

    auto g_0_yzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 686);

    auto g_0_yzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 688);

    auto g_0_yzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 689);

    auto g_0_yzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 691);

    auto g_0_yzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 692);

    auto g_0_yzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 693);

    auto g_0_yzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 695);

    auto g_0_yzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 696);

    auto g_0_yzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 697);

    auto g_0_yzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 698);

    auto g_0_yzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 700);

    auto g_0_yzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 701);

    auto g_0_yzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 702);

    auto g_0_yzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 703);

    auto g_0_yzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 704);

    auto g_0_yzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 706);

    auto g_0_yzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 707);

    auto g_0_yzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 708);

    auto g_0_yzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 709);

    auto g_0_yzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 710);

    auto g_0_yzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 711);

    auto g_0_yzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 713);

    auto g_0_yzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 714);

    auto g_0_yzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 715);

    auto g_0_yzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 716);

    auto g_0_yzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 717);

    auto g_0_yzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 718);

    auto g_0_yzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 719);

    auto g_0_zzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 720);

    auto g_0_zzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 721);

    auto g_0_zzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 722);

    auto g_0_zzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 723);

    auto g_0_zzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 724);

    auto g_0_zzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 725);

    auto g_0_zzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 726);

    auto g_0_zzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 727);

    auto g_0_zzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 728);

    auto g_0_zzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 729);

    auto g_0_zzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 730);

    auto g_0_zzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 731);

    auto g_0_zzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 732);

    auto g_0_zzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 733);

    auto g_0_zzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 734);

    auto g_0_zzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 735);

    auto g_0_zzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 736);

    auto g_0_zzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 737);

    auto g_0_zzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 738);

    auto g_0_zzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 739);

    auto g_0_zzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 740);

    auto g_0_zzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 741);

    auto g_0_zzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 742);

    auto g_0_zzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 743);

    auto g_0_zzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 744);

    auto g_0_zzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 745);

    auto g_0_zzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 746);

    auto g_0_zzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 747);

    auto g_0_zzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 748);

    auto g_0_zzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 749);

    auto g_0_zzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 750);

    auto g_0_zzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 751);

    auto g_0_zzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 752);

    auto g_0_zzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 753);

    auto g_0_zzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 754);

    auto g_0_zzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 755);

    /// Set up components of auxilary buffer : SHSK

    auto g_0_xxxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk);

    auto g_0_xxxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 1);

    auto g_0_xxxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 2);

    auto g_0_xxxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 3);

    auto g_0_xxxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 4);

    auto g_0_xxxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 5);

    auto g_0_xxxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 6);

    auto g_0_xxxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 7);

    auto g_0_xxxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 8);

    auto g_0_xxxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 9);

    auto g_0_xxxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 10);

    auto g_0_xxxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 11);

    auto g_0_xxxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 12);

    auto g_0_xxxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 13);

    auto g_0_xxxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 14);

    auto g_0_xxxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 15);

    auto g_0_xxxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 16);

    auto g_0_xxxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 17);

    auto g_0_xxxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 18);

    auto g_0_xxxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 19);

    auto g_0_xxxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 20);

    auto g_0_xxxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 21);

    auto g_0_xxxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 22);

    auto g_0_xxxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 23);

    auto g_0_xxxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 24);

    auto g_0_xxxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 25);

    auto g_0_xxxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 26);

    auto g_0_xxxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 27);

    auto g_0_xxxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 28);

    auto g_0_xxxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 29);

    auto g_0_xxxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 30);

    auto g_0_xxxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 31);

    auto g_0_xxxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 32);

    auto g_0_xxxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 33);

    auto g_0_xxxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 34);

    auto g_0_xxxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 35);

    auto g_0_xxxxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 36);

    auto g_0_xxxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 38);

    auto g_0_xxxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 41);

    auto g_0_xxxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 45);

    auto g_0_xxxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 50);

    auto g_0_xxxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 56);

    auto g_0_xxxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 63);

    auto g_0_xxxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 72);

    auto g_0_xxxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 73);

    auto g_0_xxxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 75);

    auto g_0_xxxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 78);

    auto g_0_xxxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 82);

    auto g_0_xxxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 87);

    auto g_0_xxxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 93);

    auto g_0_xxxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 108);

    auto g_0_xxxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 109);

    auto g_0_xxxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 110);

    auto g_0_xxxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 111);

    auto g_0_xxxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 112);

    auto g_0_xxxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 113);

    auto g_0_xxxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 114);

    auto g_0_xxxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 115);

    auto g_0_xxxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 116);

    auto g_0_xxxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 117);

    auto g_0_xxxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 118);

    auto g_0_xxxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 119);

    auto g_0_xxxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 120);

    auto g_0_xxxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 121);

    auto g_0_xxxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 122);

    auto g_0_xxxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 123);

    auto g_0_xxxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 124);

    auto g_0_xxxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 125);

    auto g_0_xxxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 126);

    auto g_0_xxxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 127);

    auto g_0_xxxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 128);

    auto g_0_xxxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 129);

    auto g_0_xxxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 130);

    auto g_0_xxxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 131);

    auto g_0_xxxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 132);

    auto g_0_xxxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 133);

    auto g_0_xxxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 134);

    auto g_0_xxxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 135);

    auto g_0_xxxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 136);

    auto g_0_xxxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 137);

    auto g_0_xxxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 138);

    auto g_0_xxxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 139);

    auto g_0_xxxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 140);

    auto g_0_xxxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 141);

    auto g_0_xxxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 142);

    auto g_0_xxxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 143);

    auto g_0_xxxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 180);

    auto g_0_xxxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 181);

    auto g_0_xxxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 182);

    auto g_0_xxxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 183);

    auto g_0_xxxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 184);

    auto g_0_xxxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 185);

    auto g_0_xxxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 186);

    auto g_0_xxxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 187);

    auto g_0_xxxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 188);

    auto g_0_xxxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 189);

    auto g_0_xxxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 190);

    auto g_0_xxxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 191);

    auto g_0_xxxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 192);

    auto g_0_xxxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 193);

    auto g_0_xxxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 194);

    auto g_0_xxxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 195);

    auto g_0_xxxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 196);

    auto g_0_xxxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 197);

    auto g_0_xxxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 198);

    auto g_0_xxxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 199);

    auto g_0_xxxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 200);

    auto g_0_xxxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 201);

    auto g_0_xxxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 202);

    auto g_0_xxxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 203);

    auto g_0_xxxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 204);

    auto g_0_xxxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 205);

    auto g_0_xxxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 206);

    auto g_0_xxxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 207);

    auto g_0_xxxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 208);

    auto g_0_xxxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 209);

    auto g_0_xxxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 210);

    auto g_0_xxxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 211);

    auto g_0_xxxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 212);

    auto g_0_xxxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 213);

    auto g_0_xxxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 214);

    auto g_0_xxxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 215);

    auto g_0_xxyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 216);

    auto g_0_xxyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 217);

    auto g_0_xxyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 218);

    auto g_0_xxyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 219);

    auto g_0_xxyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 220);

    auto g_0_xxyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 221);

    auto g_0_xxyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 222);

    auto g_0_xxyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 223);

    auto g_0_xxyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 224);

    auto g_0_xxyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 225);

    auto g_0_xxyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 226);

    auto g_0_xxyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 227);

    auto g_0_xxyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 228);

    auto g_0_xxyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 229);

    auto g_0_xxyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 230);

    auto g_0_xxyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 231);

    auto g_0_xxyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 232);

    auto g_0_xxyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 233);

    auto g_0_xxyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 234);

    auto g_0_xxyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 235);

    auto g_0_xxyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 236);

    auto g_0_xxyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 237);

    auto g_0_xxyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 238);

    auto g_0_xxyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 239);

    auto g_0_xxyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 240);

    auto g_0_xxyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 241);

    auto g_0_xxyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 242);

    auto g_0_xxyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 243);

    auto g_0_xxyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 244);

    auto g_0_xxyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 245);

    auto g_0_xxyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 246);

    auto g_0_xxyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 247);

    auto g_0_xxyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 248);

    auto g_0_xxyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 249);

    auto g_0_xxyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 250);

    auto g_0_xxyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 251);

    auto g_0_xxyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 253);

    auto g_0_xxyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 255);

    auto g_0_xxyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 258);

    auto g_0_xxyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 262);

    auto g_0_xxyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 267);

    auto g_0_xxyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 273);

    auto g_0_xxyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 288);

    auto g_0_xxyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 290);

    auto g_0_xxyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 293);

    auto g_0_xxyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 297);

    auto g_0_xxyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 302);

    auto g_0_xxyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 308);

    auto g_0_xxyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 315);

    auto g_0_xxzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 324);

    auto g_0_xxzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 325);

    auto g_0_xxzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 326);

    auto g_0_xxzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 327);

    auto g_0_xxzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 328);

    auto g_0_xxzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 329);

    auto g_0_xxzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 330);

    auto g_0_xxzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 331);

    auto g_0_xxzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 332);

    auto g_0_xxzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 333);

    auto g_0_xxzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 334);

    auto g_0_xxzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 335);

    auto g_0_xxzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 336);

    auto g_0_xxzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 337);

    auto g_0_xxzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 338);

    auto g_0_xxzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 339);

    auto g_0_xxzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 340);

    auto g_0_xxzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 341);

    auto g_0_xxzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 342);

    auto g_0_xxzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 343);

    auto g_0_xxzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 344);

    auto g_0_xxzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 345);

    auto g_0_xxzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 346);

    auto g_0_xxzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 347);

    auto g_0_xxzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 348);

    auto g_0_xxzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 349);

    auto g_0_xxzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 350);

    auto g_0_xxzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 351);

    auto g_0_xxzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 352);

    auto g_0_xxzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 353);

    auto g_0_xxzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 354);

    auto g_0_xxzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 355);

    auto g_0_xxzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 356);

    auto g_0_xxzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 357);

    auto g_0_xxzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 358);

    auto g_0_xxzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 359);

    auto g_0_xyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 361);

    auto g_0_xyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 363);

    auto g_0_xyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 364);

    auto g_0_xyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 366);

    auto g_0_xyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 367);

    auto g_0_xyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 368);

    auto g_0_xyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 370);

    auto g_0_xyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 371);

    auto g_0_xyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 372);

    auto g_0_xyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 373);

    auto g_0_xyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 375);

    auto g_0_xyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 376);

    auto g_0_xyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 377);

    auto g_0_xyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 378);

    auto g_0_xyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 379);

    auto g_0_xyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 381);

    auto g_0_xyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 382);

    auto g_0_xyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 383);

    auto g_0_xyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 384);

    auto g_0_xyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 385);

    auto g_0_xyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 386);

    auto g_0_xyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 388);

    auto g_0_xyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 389);

    auto g_0_xyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 390);

    auto g_0_xyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 391);

    auto g_0_xyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 392);

    auto g_0_xyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 393);

    auto g_0_xyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 394);

    auto g_0_xyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 395);

    auto g_0_xyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 436);

    auto g_0_xyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 439);

    auto g_0_xyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 440);

    auto g_0_xyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 443);

    auto g_0_xyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 444);

    auto g_0_xyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 445);

    auto g_0_xyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 448);

    auto g_0_xyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 449);

    auto g_0_xyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 450);

    auto g_0_xyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 451);

    auto g_0_xyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 454);

    auto g_0_xyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 455);

    auto g_0_xyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 456);

    auto g_0_xyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 457);

    auto g_0_xyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 458);

    auto g_0_xyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 460);

    auto g_0_xyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 461);

    auto g_0_xyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 462);

    auto g_0_xyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 463);

    auto g_0_xyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 464);

    auto g_0_xyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 465);

    auto g_0_xyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 466);

    auto g_0_xyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 467);

    auto g_0_xzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 506);

    auto g_0_xzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 508);

    auto g_0_xzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 509);

    auto g_0_xzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 511);

    auto g_0_xzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 512);

    auto g_0_xzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 513);

    auto g_0_xzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 515);

    auto g_0_xzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 516);

    auto g_0_xzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 517);

    auto g_0_xzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 518);

    auto g_0_xzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 520);

    auto g_0_xzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 521);

    auto g_0_xzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 522);

    auto g_0_xzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 523);

    auto g_0_xzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 524);

    auto g_0_xzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 526);

    auto g_0_xzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 527);

    auto g_0_xzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 528);

    auto g_0_xzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 529);

    auto g_0_xzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 530);

    auto g_0_xzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 531);

    auto g_0_xzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 532);

    auto g_0_xzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 533);

    auto g_0_xzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 534);

    auto g_0_xzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 535);

    auto g_0_xzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 536);

    auto g_0_xzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 537);

    auto g_0_xzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 538);

    auto g_0_xzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 539);

    auto g_0_yyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 540);

    auto g_0_yyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 541);

    auto g_0_yyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 542);

    auto g_0_yyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 543);

    auto g_0_yyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 544);

    auto g_0_yyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 545);

    auto g_0_yyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 546);

    auto g_0_yyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 547);

    auto g_0_yyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 548);

    auto g_0_yyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 549);

    auto g_0_yyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 550);

    auto g_0_yyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 551);

    auto g_0_yyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 552);

    auto g_0_yyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 553);

    auto g_0_yyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 554);

    auto g_0_yyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 555);

    auto g_0_yyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 556);

    auto g_0_yyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 557);

    auto g_0_yyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 558);

    auto g_0_yyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 559);

    auto g_0_yyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 560);

    auto g_0_yyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 561);

    auto g_0_yyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 562);

    auto g_0_yyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 563);

    auto g_0_yyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 564);

    auto g_0_yyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 565);

    auto g_0_yyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 566);

    auto g_0_yyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 567);

    auto g_0_yyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 568);

    auto g_0_yyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 569);

    auto g_0_yyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 570);

    auto g_0_yyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 571);

    auto g_0_yyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 572);

    auto g_0_yyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 573);

    auto g_0_yyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 574);

    auto g_0_yyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 575);

    auto g_0_yyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 577);

    auto g_0_yyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 579);

    auto g_0_yyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 582);

    auto g_0_yyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 586);

    auto g_0_yyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 591);

    auto g_0_yyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 597);

    auto g_0_yyyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 604);

    auto g_0_yyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 612);

    auto g_0_yyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 613);

    auto g_0_yyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 614);

    auto g_0_yyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 615);

    auto g_0_yyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 616);

    auto g_0_yyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 617);

    auto g_0_yyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 618);

    auto g_0_yyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 619);

    auto g_0_yyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 620);

    auto g_0_yyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 621);

    auto g_0_yyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 622);

    auto g_0_yyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 623);

    auto g_0_yyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 624);

    auto g_0_yyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 625);

    auto g_0_yyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 626);

    auto g_0_yyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 627);

    auto g_0_yyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 628);

    auto g_0_yyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 629);

    auto g_0_yyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 630);

    auto g_0_yyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 631);

    auto g_0_yyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 632);

    auto g_0_yyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 633);

    auto g_0_yyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 634);

    auto g_0_yyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 635);

    auto g_0_yyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 636);

    auto g_0_yyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 637);

    auto g_0_yyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 638);

    auto g_0_yyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 639);

    auto g_0_yyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 640);

    auto g_0_yyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 641);

    auto g_0_yyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 642);

    auto g_0_yyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 643);

    auto g_0_yyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 644);

    auto g_0_yyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 645);

    auto g_0_yyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 646);

    auto g_0_yyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 647);

    auto g_0_yyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 648);

    auto g_0_yyzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 649);

    auto g_0_yyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 650);

    auto g_0_yyzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 651);

    auto g_0_yyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 652);

    auto g_0_yyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 653);

    auto g_0_yyzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 654);

    auto g_0_yyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 655);

    auto g_0_yyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 656);

    auto g_0_yyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 657);

    auto g_0_yyzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 658);

    auto g_0_yyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 659);

    auto g_0_yyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 660);

    auto g_0_yyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 661);

    auto g_0_yyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 662);

    auto g_0_yyzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 663);

    auto g_0_yyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 664);

    auto g_0_yyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 665);

    auto g_0_yyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 666);

    auto g_0_yyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 667);

    auto g_0_yyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 668);

    auto g_0_yyzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 669);

    auto g_0_yyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 670);

    auto g_0_yyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 671);

    auto g_0_yyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 672);

    auto g_0_yyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 673);

    auto g_0_yyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 674);

    auto g_0_yyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 675);

    auto g_0_yyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 676);

    auto g_0_yyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 677);

    auto g_0_yyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 678);

    auto g_0_yyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 679);

    auto g_0_yyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 680);

    auto g_0_yyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 681);

    auto g_0_yyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 682);

    auto g_0_yyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 683);

    auto g_0_yzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 684);

    auto g_0_yzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 686);

    auto g_0_yzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 688);

    auto g_0_yzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 689);

    auto g_0_yzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 691);

    auto g_0_yzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 692);

    auto g_0_yzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 693);

    auto g_0_yzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 695);

    auto g_0_yzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 696);

    auto g_0_yzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 697);

    auto g_0_yzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 698);

    auto g_0_yzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 700);

    auto g_0_yzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 701);

    auto g_0_yzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 702);

    auto g_0_yzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 703);

    auto g_0_yzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 704);

    auto g_0_yzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 706);

    auto g_0_yzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 707);

    auto g_0_yzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 708);

    auto g_0_yzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 709);

    auto g_0_yzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 710);

    auto g_0_yzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 711);

    auto g_0_yzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 713);

    auto g_0_yzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 714);

    auto g_0_yzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 715);

    auto g_0_yzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 716);

    auto g_0_yzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 717);

    auto g_0_yzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 718);

    auto g_0_yzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 719);

    auto g_0_zzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 720);

    auto g_0_zzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 721);

    auto g_0_zzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 722);

    auto g_0_zzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 723);

    auto g_0_zzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 724);

    auto g_0_zzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 725);

    auto g_0_zzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 726);

    auto g_0_zzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 727);

    auto g_0_zzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 728);

    auto g_0_zzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 729);

    auto g_0_zzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 730);

    auto g_0_zzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 731);

    auto g_0_zzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 732);

    auto g_0_zzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 733);

    auto g_0_zzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 734);

    auto g_0_zzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 735);

    auto g_0_zzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 736);

    auto g_0_zzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 737);

    auto g_0_zzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 738);

    auto g_0_zzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 739);

    auto g_0_zzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 740);

    auto g_0_zzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 741);

    auto g_0_zzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 742);

    auto g_0_zzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 743);

    auto g_0_zzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 744);

    auto g_0_zzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 745);

    auto g_0_zzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 746);

    auto g_0_zzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 747);

    auto g_0_zzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 748);

    auto g_0_zzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 749);

    auto g_0_zzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 750);

    auto g_0_zzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 751);

    auto g_0_zzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 752);

    auto g_0_zzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 753);

    auto g_0_zzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 754);

    auto g_0_zzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 755);

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

    auto g_0_xxxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 58);

    auto g_0_xxxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 60);

    auto g_0_xxxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 61);

    auto g_0_xxxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 63);

    auto g_0_xxxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 64);

    auto g_0_xxxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 65);

    auto g_0_xxxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 67);

    auto g_0_xxxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 68);

    auto g_0_xxxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 69);

    auto g_0_xxxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 70);

    auto g_0_xxxxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 72);

    auto g_0_xxxxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 73);

    auto g_0_xxxxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 74);

    auto g_0_xxxxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 75);

    auto g_0_xxxxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 76);

    auto g_0_xxxxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 78);

    auto g_0_xxxxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 79);

    auto g_0_xxxxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 80);

    auto g_0_xxxxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 81);

    auto g_0_xxxxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 82);

    auto g_0_xxxxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 83);

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

    auto g_0_xxyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 340);

    auto g_0_xxyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 343);

    auto g_0_xxyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 344);

    auto g_0_xxyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 347);

    auto g_0_xxyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 348);

    auto g_0_xxyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 349);

    auto g_0_xxyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 352);

    auto g_0_xxyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 353);

    auto g_0_xxyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 354);

    auto g_0_xxyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 355);

    auto g_0_xxyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 358);

    auto g_0_xxyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 359);

    auto g_0_xxyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 360);

    auto g_0_xxyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 361);

    auto g_0_xxyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 362);

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

    auto g_0_xyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 498);

    auto g_0_xyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 499);

    auto g_0_xyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 500);

    auto g_0_xyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 501);

    auto g_0_xyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 502);

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

    auto g_0_xyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 526);

    auto g_0_xyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 527);

    auto g_0_xyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 528);

    auto g_0_xyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 529);

    auto g_0_xyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 530);

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

    auto g_0_yyyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 618);

    auto g_0_yyyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 620);

    auto g_0_yyyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 621);

    auto g_0_yyyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 623);

    auto g_0_yyyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 624);

    auto g_0_yyyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 625);

    auto g_0_yyyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 627);

    auto g_0_yyyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 628);

    auto g_0_yyyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 629);

    auto g_0_yyyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 630);

    auto g_0_yyyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 632);

    auto g_0_yyyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 633);

    auto g_0_yyyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 634);

    auto g_0_yyyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 635);

    auto g_0_yyyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 636);

    auto g_0_yyyyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 638);

    auto g_0_yyyyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 639);

    auto g_0_yyyyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 640);

    auto g_0_yyyyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 641);

    auto g_0_yyyyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 642);

    auto g_0_yyyyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 643);

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

    auto g_0_yzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 729);

    auto g_0_yzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 730);

    auto g_0_yzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 731);

    auto g_0_yzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 732);

    auto g_0_yzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 733);

    auto g_0_yzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 734);

    auto g_0_yzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 735);

    auto g_0_yzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 736);

    auto g_0_yzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 737);

    auto g_0_yzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 738);

    auto g_0_yzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 739);

    auto g_0_yzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 740);

    auto g_0_yzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 741);

    auto g_0_yzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 742);

    auto g_0_yzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 743);

    auto g_0_yzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 744);

    auto g_0_yzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 745);

    auto g_0_yzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 746);

    auto g_0_yzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 747);

    auto g_0_yzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 748);

    auto g_0_yzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 749);

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

    /// Set up components of auxilary buffer : SISK

    auto g_0_xxxxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk);

    auto g_0_xxxxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 1);

    auto g_0_xxxxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 2);

    auto g_0_xxxxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 3);

    auto g_0_xxxxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 4);

    auto g_0_xxxxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 5);

    auto g_0_xxxxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 6);

    auto g_0_xxxxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 7);

    auto g_0_xxxxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 8);

    auto g_0_xxxxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 9);

    auto g_0_xxxxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 10);

    auto g_0_xxxxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 11);

    auto g_0_xxxxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 12);

    auto g_0_xxxxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 13);

    auto g_0_xxxxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 14);

    auto g_0_xxxxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 15);

    auto g_0_xxxxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 16);

    auto g_0_xxxxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 17);

    auto g_0_xxxxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 18);

    auto g_0_xxxxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 19);

    auto g_0_xxxxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 20);

    auto g_0_xxxxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 21);

    auto g_0_xxxxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 22);

    auto g_0_xxxxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 23);

    auto g_0_xxxxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 24);

    auto g_0_xxxxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 25);

    auto g_0_xxxxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 26);

    auto g_0_xxxxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 27);

    auto g_0_xxxxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 28);

    auto g_0_xxxxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 29);

    auto g_0_xxxxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 30);

    auto g_0_xxxxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 31);

    auto g_0_xxxxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 32);

    auto g_0_xxxxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 33);

    auto g_0_xxxxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 34);

    auto g_0_xxxxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 35);

    auto g_0_xxxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 36);

    auto g_0_xxxxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 37);

    auto g_0_xxxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 38);

    auto g_0_xxxxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 39);

    auto g_0_xxxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 41);

    auto g_0_xxxxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 42);

    auto g_0_xxxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 45);

    auto g_0_xxxxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 46);

    auto g_0_xxxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 50);

    auto g_0_xxxxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 51);

    auto g_0_xxxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 56);

    auto g_0_xxxxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 57);

    auto g_0_xxxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 63);

    auto g_0_xxxxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 64);

    auto g_0_xxxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 72);

    auto g_0_xxxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 73);

    auto g_0_xxxxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 74);

    auto g_0_xxxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 75);

    auto g_0_xxxxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 76);

    auto g_0_xxxxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 77);

    auto g_0_xxxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 78);

    auto g_0_xxxxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 79);

    auto g_0_xxxxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 80);

    auto g_0_xxxxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 81);

    auto g_0_xxxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 82);

    auto g_0_xxxxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 83);

    auto g_0_xxxxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 84);

    auto g_0_xxxxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 85);

    auto g_0_xxxxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 86);

    auto g_0_xxxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 87);

    auto g_0_xxxxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 88);

    auto g_0_xxxxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 89);

    auto g_0_xxxxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 90);

    auto g_0_xxxxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 91);

    auto g_0_xxxxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 92);

    auto g_0_xxxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 93);

    auto g_0_xxxxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 94);

    auto g_0_xxxxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 95);

    auto g_0_xxxxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 96);

    auto g_0_xxxxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 97);

    auto g_0_xxxxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 98);

    auto g_0_xxxxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 99);

    auto g_0_xxxxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 101);

    auto g_0_xxxxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 102);

    auto g_0_xxxxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 103);

    auto g_0_xxxxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 104);

    auto g_0_xxxxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 105);

    auto g_0_xxxxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 106);

    auto g_0_xxxxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 107);

    auto g_0_xxxxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 108);

    auto g_0_xxxxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 109);

    auto g_0_xxxxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 110);

    auto g_0_xxxxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 111);

    auto g_0_xxxxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 112);

    auto g_0_xxxxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 113);

    auto g_0_xxxxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 114);

    auto g_0_xxxxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 115);

    auto g_0_xxxxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 116);

    auto g_0_xxxxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 117);

    auto g_0_xxxxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 118);

    auto g_0_xxxxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 119);

    auto g_0_xxxxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 120);

    auto g_0_xxxxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 121);

    auto g_0_xxxxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 122);

    auto g_0_xxxxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 123);

    auto g_0_xxxxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 124);

    auto g_0_xxxxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 125);

    auto g_0_xxxxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 126);

    auto g_0_xxxxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 127);

    auto g_0_xxxxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 128);

    auto g_0_xxxxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 129);

    auto g_0_xxxxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 130);

    auto g_0_xxxxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 131);

    auto g_0_xxxxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 132);

    auto g_0_xxxxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 133);

    auto g_0_xxxxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 134);

    auto g_0_xxxxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 135);

    auto g_0_xxxxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 136);

    auto g_0_xxxxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 137);

    auto g_0_xxxxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 138);

    auto g_0_xxxxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 139);

    auto g_0_xxxxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 140);

    auto g_0_xxxxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 141);

    auto g_0_xxxxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 142);

    auto g_0_xxxxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 143);

    auto g_0_xxxxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 180);

    auto g_0_xxxxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 181);

    auto g_0_xxxxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 182);

    auto g_0_xxxxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 183);

    auto g_0_xxxxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 184);

    auto g_0_xxxxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 185);

    auto g_0_xxxxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 186);

    auto g_0_xxxxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 187);

    auto g_0_xxxxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 188);

    auto g_0_xxxxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 189);

    auto g_0_xxxxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 190);

    auto g_0_xxxxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 191);

    auto g_0_xxxxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 192);

    auto g_0_xxxxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 193);

    auto g_0_xxxxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 194);

    auto g_0_xxxxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 195);

    auto g_0_xxxxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 196);

    auto g_0_xxxxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 197);

    auto g_0_xxxxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 198);

    auto g_0_xxxxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 199);

    auto g_0_xxxxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 200);

    auto g_0_xxxxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 201);

    auto g_0_xxxxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 202);

    auto g_0_xxxxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 203);

    auto g_0_xxxxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 204);

    auto g_0_xxxxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 205);

    auto g_0_xxxxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 206);

    auto g_0_xxxxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 207);

    auto g_0_xxxxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 208);

    auto g_0_xxxxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 209);

    auto g_0_xxxxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 210);

    auto g_0_xxxxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 211);

    auto g_0_xxxxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 212);

    auto g_0_xxxxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 213);

    auto g_0_xxxxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 214);

    auto g_0_xxxxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 215);

    auto g_0_xxxyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 216);

    auto g_0_xxxyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 217);

    auto g_0_xxxyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 218);

    auto g_0_xxxyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 219);

    auto g_0_xxxyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 220);

    auto g_0_xxxyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 221);

    auto g_0_xxxyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 222);

    auto g_0_xxxyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 223);

    auto g_0_xxxyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 224);

    auto g_0_xxxyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 225);

    auto g_0_xxxyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 226);

    auto g_0_xxxyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 227);

    auto g_0_xxxyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 228);

    auto g_0_xxxyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 229);

    auto g_0_xxxyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 230);

    auto g_0_xxxyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 231);

    auto g_0_xxxyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 232);

    auto g_0_xxxyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 233);

    auto g_0_xxxyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 234);

    auto g_0_xxxyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 235);

    auto g_0_xxxyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 236);

    auto g_0_xxxyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 237);

    auto g_0_xxxyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 238);

    auto g_0_xxxyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 239);

    auto g_0_xxxyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 240);

    auto g_0_xxxyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 241);

    auto g_0_xxxyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 242);

    auto g_0_xxxyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 243);

    auto g_0_xxxyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 244);

    auto g_0_xxxyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 245);

    auto g_0_xxxyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 246);

    auto g_0_xxxyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 247);

    auto g_0_xxxyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 248);

    auto g_0_xxxyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 249);

    auto g_0_xxxyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 250);

    auto g_0_xxxyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 251);

    auto g_0_xxxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 253);

    auto g_0_xxxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 255);

    auto g_0_xxxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 258);

    auto g_0_xxxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 262);

    auto g_0_xxxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 267);

    auto g_0_xxxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 273);

    auto g_0_xxxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 288);

    auto g_0_xxxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 290);

    auto g_0_xxxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 293);

    auto g_0_xxxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 297);

    auto g_0_xxxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 302);

    auto g_0_xxxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 308);

    auto g_0_xxxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 315);

    auto g_0_xxxzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 324);

    auto g_0_xxxzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 325);

    auto g_0_xxxzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 326);

    auto g_0_xxxzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 327);

    auto g_0_xxxzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 328);

    auto g_0_xxxzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 329);

    auto g_0_xxxzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 330);

    auto g_0_xxxzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 331);

    auto g_0_xxxzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 332);

    auto g_0_xxxzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 333);

    auto g_0_xxxzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 334);

    auto g_0_xxxzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 335);

    auto g_0_xxxzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 336);

    auto g_0_xxxzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 337);

    auto g_0_xxxzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 338);

    auto g_0_xxxzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 339);

    auto g_0_xxxzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 340);

    auto g_0_xxxzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 341);

    auto g_0_xxxzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 342);

    auto g_0_xxxzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 343);

    auto g_0_xxxzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 344);

    auto g_0_xxxzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 345);

    auto g_0_xxxzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 346);

    auto g_0_xxxzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 347);

    auto g_0_xxxzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 348);

    auto g_0_xxxzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 349);

    auto g_0_xxxzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 350);

    auto g_0_xxxzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 351);

    auto g_0_xxxzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 352);

    auto g_0_xxxzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 353);

    auto g_0_xxxzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 354);

    auto g_0_xxxzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 355);

    auto g_0_xxxzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 356);

    auto g_0_xxxzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 357);

    auto g_0_xxxzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 358);

    auto g_0_xxxzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 359);

    auto g_0_xxyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 360);

    auto g_0_xxyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 361);

    auto g_0_xxyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 362);

    auto g_0_xxyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 363);

    auto g_0_xxyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 364);

    auto g_0_xxyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 365);

    auto g_0_xxyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 366);

    auto g_0_xxyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 367);

    auto g_0_xxyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 368);

    auto g_0_xxyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 369);

    auto g_0_xxyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 370);

    auto g_0_xxyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 371);

    auto g_0_xxyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 372);

    auto g_0_xxyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 373);

    auto g_0_xxyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 374);

    auto g_0_xxyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 375);

    auto g_0_xxyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 376);

    auto g_0_xxyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 377);

    auto g_0_xxyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 378);

    auto g_0_xxyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 379);

    auto g_0_xxyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 380);

    auto g_0_xxyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 381);

    auto g_0_xxyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 382);

    auto g_0_xxyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 383);

    auto g_0_xxyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 384);

    auto g_0_xxyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 385);

    auto g_0_xxyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 386);

    auto g_0_xxyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 387);

    auto g_0_xxyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 388);

    auto g_0_xxyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 389);

    auto g_0_xxyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 390);

    auto g_0_xxyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 391);

    auto g_0_xxyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 392);

    auto g_0_xxyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 393);

    auto g_0_xxyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 394);

    auto g_0_xxyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 395);

    auto g_0_xxyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 397);

    auto g_0_xxyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 399);

    auto g_0_xxyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 402);

    auto g_0_xxyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 406);

    auto g_0_xxyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 411);

    auto g_0_xxyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 417);

    auto g_0_xxyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 432);

    auto g_0_xxyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 433);

    auto g_0_xxyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 434);

    auto g_0_xxyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 435);

    auto g_0_xxyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 436);

    auto g_0_xxyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 437);

    auto g_0_xxyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 438);

    auto g_0_xxyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 439);

    auto g_0_xxyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 440);

    auto g_0_xxyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 441);

    auto g_0_xxyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 442);

    auto g_0_xxyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 443);

    auto g_0_xxyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 444);

    auto g_0_xxyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 445);

    auto g_0_xxyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 446);

    auto g_0_xxyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 447);

    auto g_0_xxyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 448);

    auto g_0_xxyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 449);

    auto g_0_xxyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 450);

    auto g_0_xxyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 451);

    auto g_0_xxyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 452);

    auto g_0_xxyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 453);

    auto g_0_xxyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 454);

    auto g_0_xxyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 455);

    auto g_0_xxyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 456);

    auto g_0_xxyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 457);

    auto g_0_xxyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 458);

    auto g_0_xxyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 459);

    auto g_0_xxyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 460);

    auto g_0_xxyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 461);

    auto g_0_xxyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 462);

    auto g_0_xxyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 463);

    auto g_0_xxyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 464);

    auto g_0_xxyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 465);

    auto g_0_xxyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 466);

    auto g_0_xxyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 467);

    auto g_0_xxyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 468);

    auto g_0_xxyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 470);

    auto g_0_xxyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 473);

    auto g_0_xxyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 477);

    auto g_0_xxyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 482);

    auto g_0_xxyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 488);

    auto g_0_xxyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 495);

    auto g_0_xxzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 504);

    auto g_0_xxzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 505);

    auto g_0_xxzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 506);

    auto g_0_xxzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 507);

    auto g_0_xxzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 508);

    auto g_0_xxzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 509);

    auto g_0_xxzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 510);

    auto g_0_xxzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 511);

    auto g_0_xxzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 512);

    auto g_0_xxzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 513);

    auto g_0_xxzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 514);

    auto g_0_xxzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 515);

    auto g_0_xxzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 516);

    auto g_0_xxzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 517);

    auto g_0_xxzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 518);

    auto g_0_xxzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 519);

    auto g_0_xxzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 520);

    auto g_0_xxzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 521);

    auto g_0_xxzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 522);

    auto g_0_xxzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 523);

    auto g_0_xxzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 524);

    auto g_0_xxzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 525);

    auto g_0_xxzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 526);

    auto g_0_xxzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 527);

    auto g_0_xxzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 528);

    auto g_0_xxzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 529);

    auto g_0_xxzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 530);

    auto g_0_xxzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 531);

    auto g_0_xxzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 532);

    auto g_0_xxzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 533);

    auto g_0_xxzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 534);

    auto g_0_xxzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 535);

    auto g_0_xxzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 536);

    auto g_0_xxzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 537);

    auto g_0_xxzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 538);

    auto g_0_xxzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 539);

    auto g_0_xyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 540);

    auto g_0_xyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 541);

    auto g_0_xyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 543);

    auto g_0_xyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 544);

    auto g_0_xyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 546);

    auto g_0_xyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 547);

    auto g_0_xyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 548);

    auto g_0_xyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 550);

    auto g_0_xyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 551);

    auto g_0_xyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 552);

    auto g_0_xyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 553);

    auto g_0_xyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 555);

    auto g_0_xyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 556);

    auto g_0_xyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 557);

    auto g_0_xyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 558);

    auto g_0_xyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 559);

    auto g_0_xyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 561);

    auto g_0_xyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 562);

    auto g_0_xyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 563);

    auto g_0_xyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 564);

    auto g_0_xyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 565);

    auto g_0_xyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 566);

    auto g_0_xyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 568);

    auto g_0_xyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 569);

    auto g_0_xyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 570);

    auto g_0_xyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 571);

    auto g_0_xyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 572);

    auto g_0_xyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 573);

    auto g_0_xyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 574);

    auto g_0_xyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 575);

    auto g_0_xyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 616);

    auto g_0_xyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 619);

    auto g_0_xyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 620);

    auto g_0_xyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 623);

    auto g_0_xyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 624);

    auto g_0_xyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 625);

    auto g_0_xyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 628);

    auto g_0_xyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 629);

    auto g_0_xyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 630);

    auto g_0_xyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 631);

    auto g_0_xyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 634);

    auto g_0_xyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 635);

    auto g_0_xyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 636);

    auto g_0_xyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 637);

    auto g_0_xyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 638);

    auto g_0_xyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 640);

    auto g_0_xyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 641);

    auto g_0_xyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 642);

    auto g_0_xyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 643);

    auto g_0_xyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 644);

    auto g_0_xyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 645);

    auto g_0_xyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 646);

    auto g_0_xyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 647);

    auto g_0_xyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 652);

    auto g_0_xyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 655);

    auto g_0_xyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 656);

    auto g_0_xyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 659);

    auto g_0_xyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 660);

    auto g_0_xyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 661);

    auto g_0_xyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 664);

    auto g_0_xyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 665);

    auto g_0_xyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 666);

    auto g_0_xyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 667);

    auto g_0_xyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 670);

    auto g_0_xyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 671);

    auto g_0_xyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 672);

    auto g_0_xyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 673);

    auto g_0_xyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 674);

    auto g_0_xyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 676);

    auto g_0_xyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 677);

    auto g_0_xyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 678);

    auto g_0_xyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 679);

    auto g_0_xyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 680);

    auto g_0_xyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 681);

    auto g_0_xyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 682);

    auto g_0_xyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 683);

    auto g_0_xzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 720);

    auto g_0_xzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 722);

    auto g_0_xzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 724);

    auto g_0_xzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 725);

    auto g_0_xzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 727);

    auto g_0_xzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 728);

    auto g_0_xzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 729);

    auto g_0_xzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 731);

    auto g_0_xzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 732);

    auto g_0_xzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 733);

    auto g_0_xzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 734);

    auto g_0_xzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 736);

    auto g_0_xzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 737);

    auto g_0_xzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 738);

    auto g_0_xzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 739);

    auto g_0_xzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 740);

    auto g_0_xzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 742);

    auto g_0_xzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 743);

    auto g_0_xzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 744);

    auto g_0_xzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 745);

    auto g_0_xzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 746);

    auto g_0_xzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 747);

    auto g_0_xzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 748);

    auto g_0_xzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 749);

    auto g_0_xzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 750);

    auto g_0_xzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 751);

    auto g_0_xzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 752);

    auto g_0_xzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 753);

    auto g_0_xzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 754);

    auto g_0_xzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 755);

    auto g_0_yyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 756);

    auto g_0_yyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 757);

    auto g_0_yyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 758);

    auto g_0_yyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 759);

    auto g_0_yyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 760);

    auto g_0_yyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 761);

    auto g_0_yyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 762);

    auto g_0_yyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 763);

    auto g_0_yyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 764);

    auto g_0_yyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 765);

    auto g_0_yyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 766);

    auto g_0_yyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 767);

    auto g_0_yyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 768);

    auto g_0_yyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 769);

    auto g_0_yyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 770);

    auto g_0_yyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 771);

    auto g_0_yyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 772);

    auto g_0_yyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 773);

    auto g_0_yyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 774);

    auto g_0_yyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 775);

    auto g_0_yyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 776);

    auto g_0_yyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 777);

    auto g_0_yyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 778);

    auto g_0_yyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 779);

    auto g_0_yyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 780);

    auto g_0_yyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 781);

    auto g_0_yyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 782);

    auto g_0_yyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 783);

    auto g_0_yyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 784);

    auto g_0_yyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 785);

    auto g_0_yyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 786);

    auto g_0_yyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 787);

    auto g_0_yyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 788);

    auto g_0_yyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 789);

    auto g_0_yyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 790);

    auto g_0_yyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 791);

    auto g_0_yyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 793);

    auto g_0_yyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 794);

    auto g_0_yyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 795);

    auto g_0_yyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 796);

    auto g_0_yyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 797);

    auto g_0_yyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 798);

    auto g_0_yyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 799);

    auto g_0_yyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 800);

    auto g_0_yyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 801);

    auto g_0_yyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 802);

    auto g_0_yyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 803);

    auto g_0_yyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 804);

    auto g_0_yyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 805);

    auto g_0_yyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 806);

    auto g_0_yyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 807);

    auto g_0_yyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 808);

    auto g_0_yyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 809);

    auto g_0_yyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 810);

    auto g_0_yyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 811);

    auto g_0_yyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 812);

    auto g_0_yyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 813);

    auto g_0_yyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 814);

    auto g_0_yyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 815);

    auto g_0_yyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 816);

    auto g_0_yyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 817);

    auto g_0_yyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 818);

    auto g_0_yyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 819);

    auto g_0_yyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 820);

    auto g_0_yyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 821);

    auto g_0_yyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 822);

    auto g_0_yyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 823);

    auto g_0_yyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 824);

    auto g_0_yyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 825);

    auto g_0_yyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 826);

    auto g_0_yyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 827);

    auto g_0_yyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 828);

    auto g_0_yyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 829);

    auto g_0_yyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 830);

    auto g_0_yyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 831);

    auto g_0_yyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 832);

    auto g_0_yyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 833);

    auto g_0_yyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 834);

    auto g_0_yyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 835);

    auto g_0_yyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 836);

    auto g_0_yyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 837);

    auto g_0_yyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 838);

    auto g_0_yyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 839);

    auto g_0_yyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 840);

    auto g_0_yyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 841);

    auto g_0_yyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 842);

    auto g_0_yyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 843);

    auto g_0_yyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 844);

    auto g_0_yyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 845);

    auto g_0_yyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 846);

    auto g_0_yyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 847);

    auto g_0_yyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 848);

    auto g_0_yyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 849);

    auto g_0_yyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 850);

    auto g_0_yyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 851);

    auto g_0_yyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 852);

    auto g_0_yyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 853);

    auto g_0_yyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 854);

    auto g_0_yyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 855);

    auto g_0_yyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 856);

    auto g_0_yyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 857);

    auto g_0_yyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 858);

    auto g_0_yyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 859);

    auto g_0_yyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 860);

    auto g_0_yyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 861);

    auto g_0_yyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 862);

    auto g_0_yyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 863);

    auto g_0_yyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 864);

    auto g_0_yyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 865);

    auto g_0_yyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 866);

    auto g_0_yyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 867);

    auto g_0_yyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 868);

    auto g_0_yyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 869);

    auto g_0_yyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 870);

    auto g_0_yyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 871);

    auto g_0_yyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 872);

    auto g_0_yyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 873);

    auto g_0_yyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 874);

    auto g_0_yyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 875);

    auto g_0_yyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 876);

    auto g_0_yyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 877);

    auto g_0_yyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 878);

    auto g_0_yyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 879);

    auto g_0_yyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 880);

    auto g_0_yyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 881);

    auto g_0_yyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 882);

    auto g_0_yyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 883);

    auto g_0_yyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 884);

    auto g_0_yyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 885);

    auto g_0_yyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 886);

    auto g_0_yyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 887);

    auto g_0_yyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 888);

    auto g_0_yyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 889);

    auto g_0_yyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 890);

    auto g_0_yyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 891);

    auto g_0_yyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 892);

    auto g_0_yyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 893);

    auto g_0_yyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 894);

    auto g_0_yyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 895);

    auto g_0_yyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 896);

    auto g_0_yyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 897);

    auto g_0_yyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 898);

    auto g_0_yyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 899);

    auto g_0_yyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 900);

    auto g_0_yyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 901);

    auto g_0_yyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 902);

    auto g_0_yyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 903);

    auto g_0_yyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 904);

    auto g_0_yyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 905);

    auto g_0_yyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 906);

    auto g_0_yyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 907);

    auto g_0_yyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 908);

    auto g_0_yyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 909);

    auto g_0_yyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 910);

    auto g_0_yyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 911);

    auto g_0_yyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 912);

    auto g_0_yyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 913);

    auto g_0_yyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 914);

    auto g_0_yyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 915);

    auto g_0_yyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 916);

    auto g_0_yyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 917);

    auto g_0_yyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 918);

    auto g_0_yyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 919);

    auto g_0_yyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 920);

    auto g_0_yyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 921);

    auto g_0_yyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 922);

    auto g_0_yyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 923);

    auto g_0_yyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 924);

    auto g_0_yyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 925);

    auto g_0_yyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 926);

    auto g_0_yyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 927);

    auto g_0_yyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 928);

    auto g_0_yyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 929);

    auto g_0_yyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 930);

    auto g_0_yyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 931);

    auto g_0_yyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 932);

    auto g_0_yyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 933);

    auto g_0_yyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 934);

    auto g_0_yyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 935);

    auto g_0_yzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 936);

    auto g_0_yzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 937);

    auto g_0_yzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 938);

    auto g_0_yzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 939);

    auto g_0_yzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 940);

    auto g_0_yzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 941);

    auto g_0_yzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 942);

    auto g_0_yzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 943);

    auto g_0_yzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 944);

    auto g_0_yzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 945);

    auto g_0_yzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 946);

    auto g_0_yzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 947);

    auto g_0_yzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 948);

    auto g_0_yzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 949);

    auto g_0_yzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 950);

    auto g_0_yzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 951);

    auto g_0_yzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 952);

    auto g_0_yzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 953);

    auto g_0_yzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 954);

    auto g_0_yzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 955);

    auto g_0_yzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 956);

    auto g_0_yzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 957);

    auto g_0_yzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 958);

    auto g_0_yzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 959);

    auto g_0_yzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 960);

    auto g_0_yzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 961);

    auto g_0_yzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 962);

    auto g_0_yzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 963);

    auto g_0_yzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 964);

    auto g_0_yzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 965);

    auto g_0_yzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 966);

    auto g_0_yzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 967);

    auto g_0_yzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 968);

    auto g_0_yzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 969);

    auto g_0_yzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 970);

    auto g_0_yzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 971);

    auto g_0_zzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 972);

    auto g_0_zzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 973);

    auto g_0_zzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 974);

    auto g_0_zzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 975);

    auto g_0_zzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 976);

    auto g_0_zzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 977);

    auto g_0_zzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 978);

    auto g_0_zzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 979);

    auto g_0_zzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 980);

    auto g_0_zzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 981);

    auto g_0_zzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 982);

    auto g_0_zzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 983);

    auto g_0_zzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 984);

    auto g_0_zzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 985);

    auto g_0_zzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 986);

    auto g_0_zzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 987);

    auto g_0_zzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 988);

    auto g_0_zzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 989);

    auto g_0_zzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 990);

    auto g_0_zzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 991);

    auto g_0_zzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 992);

    auto g_0_zzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 993);

    auto g_0_zzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 994);

    auto g_0_zzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 995);

    auto g_0_zzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 996);

    auto g_0_zzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 997);

    auto g_0_zzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 998);

    auto g_0_zzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 999);

    auto g_0_zzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 1000);

    auto g_0_zzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 1001);

    auto g_0_zzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 1002);

    auto g_0_zzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 1003);

    auto g_0_zzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 1004);

    auto g_0_zzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 1005);

    auto g_0_zzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 1006);

    auto g_0_zzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 1007);

    /// Set up components of auxilary buffer : SISK

    auto g_0_xxxxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk);

    auto g_0_xxxxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 1);

    auto g_0_xxxxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 2);

    auto g_0_xxxxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 3);

    auto g_0_xxxxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 4);

    auto g_0_xxxxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 5);

    auto g_0_xxxxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 6);

    auto g_0_xxxxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 7);

    auto g_0_xxxxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 8);

    auto g_0_xxxxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 9);

    auto g_0_xxxxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 10);

    auto g_0_xxxxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 11);

    auto g_0_xxxxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 12);

    auto g_0_xxxxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 13);

    auto g_0_xxxxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 14);

    auto g_0_xxxxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 15);

    auto g_0_xxxxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 16);

    auto g_0_xxxxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 17);

    auto g_0_xxxxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 18);

    auto g_0_xxxxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 19);

    auto g_0_xxxxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 20);

    auto g_0_xxxxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 21);

    auto g_0_xxxxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 22);

    auto g_0_xxxxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 23);

    auto g_0_xxxxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 24);

    auto g_0_xxxxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 25);

    auto g_0_xxxxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 26);

    auto g_0_xxxxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 27);

    auto g_0_xxxxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 28);

    auto g_0_xxxxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 29);

    auto g_0_xxxxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 30);

    auto g_0_xxxxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 31);

    auto g_0_xxxxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 32);

    auto g_0_xxxxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 33);

    auto g_0_xxxxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 34);

    auto g_0_xxxxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 35);

    auto g_0_xxxxxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 36);

    auto g_0_xxxxxy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 37);

    auto g_0_xxxxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 38);

    auto g_0_xxxxxy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 39);

    auto g_0_xxxxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 41);

    auto g_0_xxxxxy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 42);

    auto g_0_xxxxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 45);

    auto g_0_xxxxxy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 46);

    auto g_0_xxxxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 50);

    auto g_0_xxxxxy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 51);

    auto g_0_xxxxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 56);

    auto g_0_xxxxxy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 57);

    auto g_0_xxxxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 63);

    auto g_0_xxxxxy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 64);

    auto g_0_xxxxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 72);

    auto g_0_xxxxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 73);

    auto g_0_xxxxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 74);

    auto g_0_xxxxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 75);

    auto g_0_xxxxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 76);

    auto g_0_xxxxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 77);

    auto g_0_xxxxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 78);

    auto g_0_xxxxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 79);

    auto g_0_xxxxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 80);

    auto g_0_xxxxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 81);

    auto g_0_xxxxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 82);

    auto g_0_xxxxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 83);

    auto g_0_xxxxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 84);

    auto g_0_xxxxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 85);

    auto g_0_xxxxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 86);

    auto g_0_xxxxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 87);

    auto g_0_xxxxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 88);

    auto g_0_xxxxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 89);

    auto g_0_xxxxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 90);

    auto g_0_xxxxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 91);

    auto g_0_xxxxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 92);

    auto g_0_xxxxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 93);

    auto g_0_xxxxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 94);

    auto g_0_xxxxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 95);

    auto g_0_xxxxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 96);

    auto g_0_xxxxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 97);

    auto g_0_xxxxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 98);

    auto g_0_xxxxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 99);

    auto g_0_xxxxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 101);

    auto g_0_xxxxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 102);

    auto g_0_xxxxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 103);

    auto g_0_xxxxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 104);

    auto g_0_xxxxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 105);

    auto g_0_xxxxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 106);

    auto g_0_xxxxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 107);

    auto g_0_xxxxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 108);

    auto g_0_xxxxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 109);

    auto g_0_xxxxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 110);

    auto g_0_xxxxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 111);

    auto g_0_xxxxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 112);

    auto g_0_xxxxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 113);

    auto g_0_xxxxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 114);

    auto g_0_xxxxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 115);

    auto g_0_xxxxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 116);

    auto g_0_xxxxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 117);

    auto g_0_xxxxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 118);

    auto g_0_xxxxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 119);

    auto g_0_xxxxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 120);

    auto g_0_xxxxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 121);

    auto g_0_xxxxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 122);

    auto g_0_xxxxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 123);

    auto g_0_xxxxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 124);

    auto g_0_xxxxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 125);

    auto g_0_xxxxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 126);

    auto g_0_xxxxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 127);

    auto g_0_xxxxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 128);

    auto g_0_xxxxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 129);

    auto g_0_xxxxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 130);

    auto g_0_xxxxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 131);

    auto g_0_xxxxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 132);

    auto g_0_xxxxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 133);

    auto g_0_xxxxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 134);

    auto g_0_xxxxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 135);

    auto g_0_xxxxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 136);

    auto g_0_xxxxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 137);

    auto g_0_xxxxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 138);

    auto g_0_xxxxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 139);

    auto g_0_xxxxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 140);

    auto g_0_xxxxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 141);

    auto g_0_xxxxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 142);

    auto g_0_xxxxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 143);

    auto g_0_xxxxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 180);

    auto g_0_xxxxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 181);

    auto g_0_xxxxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 182);

    auto g_0_xxxxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 183);

    auto g_0_xxxxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 184);

    auto g_0_xxxxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 185);

    auto g_0_xxxxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 186);

    auto g_0_xxxxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 187);

    auto g_0_xxxxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 188);

    auto g_0_xxxxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 189);

    auto g_0_xxxxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 190);

    auto g_0_xxxxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 191);

    auto g_0_xxxxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 192);

    auto g_0_xxxxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 193);

    auto g_0_xxxxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 194);

    auto g_0_xxxxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 195);

    auto g_0_xxxxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 196);

    auto g_0_xxxxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 197);

    auto g_0_xxxxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 198);

    auto g_0_xxxxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 199);

    auto g_0_xxxxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 200);

    auto g_0_xxxxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 201);

    auto g_0_xxxxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 202);

    auto g_0_xxxxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 203);

    auto g_0_xxxxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 204);

    auto g_0_xxxxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 205);

    auto g_0_xxxxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 206);

    auto g_0_xxxxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 207);

    auto g_0_xxxxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 208);

    auto g_0_xxxxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 209);

    auto g_0_xxxxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 210);

    auto g_0_xxxxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 211);

    auto g_0_xxxxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 212);

    auto g_0_xxxxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 213);

    auto g_0_xxxxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 214);

    auto g_0_xxxxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 215);

    auto g_0_xxxyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 216);

    auto g_0_xxxyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 217);

    auto g_0_xxxyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 218);

    auto g_0_xxxyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 219);

    auto g_0_xxxyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 220);

    auto g_0_xxxyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 221);

    auto g_0_xxxyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 222);

    auto g_0_xxxyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 223);

    auto g_0_xxxyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 224);

    auto g_0_xxxyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 225);

    auto g_0_xxxyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 226);

    auto g_0_xxxyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 227);

    auto g_0_xxxyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 228);

    auto g_0_xxxyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 229);

    auto g_0_xxxyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 230);

    auto g_0_xxxyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 231);

    auto g_0_xxxyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 232);

    auto g_0_xxxyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 233);

    auto g_0_xxxyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 234);

    auto g_0_xxxyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 235);

    auto g_0_xxxyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 236);

    auto g_0_xxxyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 237);

    auto g_0_xxxyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 238);

    auto g_0_xxxyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 239);

    auto g_0_xxxyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 240);

    auto g_0_xxxyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 241);

    auto g_0_xxxyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 242);

    auto g_0_xxxyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 243);

    auto g_0_xxxyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 244);

    auto g_0_xxxyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 245);

    auto g_0_xxxyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 246);

    auto g_0_xxxyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 247);

    auto g_0_xxxyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 248);

    auto g_0_xxxyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 249);

    auto g_0_xxxyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 250);

    auto g_0_xxxyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 251);

    auto g_0_xxxyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 253);

    auto g_0_xxxyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 255);

    auto g_0_xxxyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 258);

    auto g_0_xxxyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 262);

    auto g_0_xxxyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 267);

    auto g_0_xxxyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 273);

    auto g_0_xxxyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 288);

    auto g_0_xxxyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 290);

    auto g_0_xxxyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 293);

    auto g_0_xxxyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 297);

    auto g_0_xxxyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 302);

    auto g_0_xxxyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 308);

    auto g_0_xxxyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 315);

    auto g_0_xxxzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 324);

    auto g_0_xxxzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 325);

    auto g_0_xxxzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 326);

    auto g_0_xxxzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 327);

    auto g_0_xxxzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 328);

    auto g_0_xxxzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 329);

    auto g_0_xxxzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 330);

    auto g_0_xxxzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 331);

    auto g_0_xxxzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 332);

    auto g_0_xxxzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 333);

    auto g_0_xxxzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 334);

    auto g_0_xxxzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 335);

    auto g_0_xxxzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 336);

    auto g_0_xxxzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 337);

    auto g_0_xxxzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 338);

    auto g_0_xxxzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 339);

    auto g_0_xxxzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 340);

    auto g_0_xxxzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 341);

    auto g_0_xxxzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 342);

    auto g_0_xxxzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 343);

    auto g_0_xxxzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 344);

    auto g_0_xxxzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 345);

    auto g_0_xxxzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 346);

    auto g_0_xxxzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 347);

    auto g_0_xxxzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 348);

    auto g_0_xxxzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 349);

    auto g_0_xxxzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 350);

    auto g_0_xxxzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 351);

    auto g_0_xxxzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 352);

    auto g_0_xxxzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 353);

    auto g_0_xxxzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 354);

    auto g_0_xxxzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 355);

    auto g_0_xxxzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 356);

    auto g_0_xxxzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 357);

    auto g_0_xxxzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 358);

    auto g_0_xxxzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 359);

    auto g_0_xxyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 360);

    auto g_0_xxyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 361);

    auto g_0_xxyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 362);

    auto g_0_xxyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 363);

    auto g_0_xxyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 364);

    auto g_0_xxyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 365);

    auto g_0_xxyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 366);

    auto g_0_xxyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 367);

    auto g_0_xxyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 368);

    auto g_0_xxyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 369);

    auto g_0_xxyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 370);

    auto g_0_xxyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 371);

    auto g_0_xxyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 372);

    auto g_0_xxyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 373);

    auto g_0_xxyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 374);

    auto g_0_xxyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 375);

    auto g_0_xxyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 376);

    auto g_0_xxyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 377);

    auto g_0_xxyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 378);

    auto g_0_xxyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 379);

    auto g_0_xxyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 380);

    auto g_0_xxyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 381);

    auto g_0_xxyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 382);

    auto g_0_xxyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 383);

    auto g_0_xxyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 384);

    auto g_0_xxyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 385);

    auto g_0_xxyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 386);

    auto g_0_xxyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 387);

    auto g_0_xxyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 388);

    auto g_0_xxyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 389);

    auto g_0_xxyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 390);

    auto g_0_xxyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 391);

    auto g_0_xxyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 392);

    auto g_0_xxyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 393);

    auto g_0_xxyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 394);

    auto g_0_xxyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 395);

    auto g_0_xxyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 397);

    auto g_0_xxyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 399);

    auto g_0_xxyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 402);

    auto g_0_xxyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 406);

    auto g_0_xxyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 411);

    auto g_0_xxyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 417);

    auto g_0_xxyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 432);

    auto g_0_xxyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 433);

    auto g_0_xxyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 434);

    auto g_0_xxyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 435);

    auto g_0_xxyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 436);

    auto g_0_xxyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 437);

    auto g_0_xxyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 438);

    auto g_0_xxyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 439);

    auto g_0_xxyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 440);

    auto g_0_xxyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 441);

    auto g_0_xxyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 442);

    auto g_0_xxyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 443);

    auto g_0_xxyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 444);

    auto g_0_xxyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 445);

    auto g_0_xxyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 446);

    auto g_0_xxyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 447);

    auto g_0_xxyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 448);

    auto g_0_xxyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 449);

    auto g_0_xxyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 450);

    auto g_0_xxyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 451);

    auto g_0_xxyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 452);

    auto g_0_xxyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 453);

    auto g_0_xxyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 454);

    auto g_0_xxyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 455);

    auto g_0_xxyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 456);

    auto g_0_xxyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 457);

    auto g_0_xxyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 458);

    auto g_0_xxyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 459);

    auto g_0_xxyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 460);

    auto g_0_xxyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 461);

    auto g_0_xxyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 462);

    auto g_0_xxyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 463);

    auto g_0_xxyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 464);

    auto g_0_xxyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 465);

    auto g_0_xxyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 466);

    auto g_0_xxyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 467);

    auto g_0_xxyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 468);

    auto g_0_xxyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 470);

    auto g_0_xxyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 473);

    auto g_0_xxyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 477);

    auto g_0_xxyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 482);

    auto g_0_xxyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 488);

    auto g_0_xxyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 495);

    auto g_0_xxzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 504);

    auto g_0_xxzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 505);

    auto g_0_xxzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 506);

    auto g_0_xxzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 507);

    auto g_0_xxzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 508);

    auto g_0_xxzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 509);

    auto g_0_xxzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 510);

    auto g_0_xxzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 511);

    auto g_0_xxzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 512);

    auto g_0_xxzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 513);

    auto g_0_xxzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 514);

    auto g_0_xxzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 515);

    auto g_0_xxzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 516);

    auto g_0_xxzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 517);

    auto g_0_xxzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 518);

    auto g_0_xxzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 519);

    auto g_0_xxzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 520);

    auto g_0_xxzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 521);

    auto g_0_xxzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 522);

    auto g_0_xxzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 523);

    auto g_0_xxzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 524);

    auto g_0_xxzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 525);

    auto g_0_xxzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 526);

    auto g_0_xxzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 527);

    auto g_0_xxzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 528);

    auto g_0_xxzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 529);

    auto g_0_xxzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 530);

    auto g_0_xxzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 531);

    auto g_0_xxzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 532);

    auto g_0_xxzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 533);

    auto g_0_xxzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 534);

    auto g_0_xxzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 535);

    auto g_0_xxzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 536);

    auto g_0_xxzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 537);

    auto g_0_xxzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 538);

    auto g_0_xxzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 539);

    auto g_0_xyyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 540);

    auto g_0_xyyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 541);

    auto g_0_xyyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 543);

    auto g_0_xyyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 544);

    auto g_0_xyyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 546);

    auto g_0_xyyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 547);

    auto g_0_xyyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 548);

    auto g_0_xyyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 550);

    auto g_0_xyyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 551);

    auto g_0_xyyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 552);

    auto g_0_xyyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 553);

    auto g_0_xyyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 555);

    auto g_0_xyyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 556);

    auto g_0_xyyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 557);

    auto g_0_xyyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 558);

    auto g_0_xyyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 559);

    auto g_0_xyyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 561);

    auto g_0_xyyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 562);

    auto g_0_xyyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 563);

    auto g_0_xyyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 564);

    auto g_0_xyyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 565);

    auto g_0_xyyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 566);

    auto g_0_xyyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 568);

    auto g_0_xyyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 569);

    auto g_0_xyyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 570);

    auto g_0_xyyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 571);

    auto g_0_xyyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 572);

    auto g_0_xyyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 573);

    auto g_0_xyyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 574);

    auto g_0_xyyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 575);

    auto g_0_xyyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 616);

    auto g_0_xyyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 619);

    auto g_0_xyyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 620);

    auto g_0_xyyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 623);

    auto g_0_xyyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 624);

    auto g_0_xyyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 625);

    auto g_0_xyyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 628);

    auto g_0_xyyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 629);

    auto g_0_xyyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 630);

    auto g_0_xyyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 631);

    auto g_0_xyyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 634);

    auto g_0_xyyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 635);

    auto g_0_xyyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 636);

    auto g_0_xyyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 637);

    auto g_0_xyyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 638);

    auto g_0_xyyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 640);

    auto g_0_xyyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 641);

    auto g_0_xyyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 642);

    auto g_0_xyyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 643);

    auto g_0_xyyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 644);

    auto g_0_xyyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 645);

    auto g_0_xyyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 646);

    auto g_0_xyyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 647);

    auto g_0_xyyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 652);

    auto g_0_xyyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 655);

    auto g_0_xyyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 656);

    auto g_0_xyyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 659);

    auto g_0_xyyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 660);

    auto g_0_xyyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 661);

    auto g_0_xyyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 664);

    auto g_0_xyyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 665);

    auto g_0_xyyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 666);

    auto g_0_xyyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 667);

    auto g_0_xyyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 670);

    auto g_0_xyyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 671);

    auto g_0_xyyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 672);

    auto g_0_xyyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 673);

    auto g_0_xyyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 674);

    auto g_0_xyyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 676);

    auto g_0_xyyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 677);

    auto g_0_xyyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 678);

    auto g_0_xyyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 679);

    auto g_0_xyyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 680);

    auto g_0_xyyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 681);

    auto g_0_xyyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 682);

    auto g_0_xyyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 683);

    auto g_0_xzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 720);

    auto g_0_xzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 722);

    auto g_0_xzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 724);

    auto g_0_xzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 725);

    auto g_0_xzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 727);

    auto g_0_xzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 728);

    auto g_0_xzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 729);

    auto g_0_xzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 731);

    auto g_0_xzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 732);

    auto g_0_xzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 733);

    auto g_0_xzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 734);

    auto g_0_xzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 736);

    auto g_0_xzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 737);

    auto g_0_xzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 738);

    auto g_0_xzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 739);

    auto g_0_xzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 740);

    auto g_0_xzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 742);

    auto g_0_xzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 743);

    auto g_0_xzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 744);

    auto g_0_xzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 745);

    auto g_0_xzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 746);

    auto g_0_xzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 747);

    auto g_0_xzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 748);

    auto g_0_xzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 749);

    auto g_0_xzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 750);

    auto g_0_xzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 751);

    auto g_0_xzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 752);

    auto g_0_xzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 753);

    auto g_0_xzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 754);

    auto g_0_xzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 755);

    auto g_0_yyyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 756);

    auto g_0_yyyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 757);

    auto g_0_yyyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 758);

    auto g_0_yyyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 759);

    auto g_0_yyyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 760);

    auto g_0_yyyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 761);

    auto g_0_yyyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 762);

    auto g_0_yyyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 763);

    auto g_0_yyyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 764);

    auto g_0_yyyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 765);

    auto g_0_yyyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 766);

    auto g_0_yyyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 767);

    auto g_0_yyyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 768);

    auto g_0_yyyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 769);

    auto g_0_yyyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 770);

    auto g_0_yyyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 771);

    auto g_0_yyyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 772);

    auto g_0_yyyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 773);

    auto g_0_yyyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 774);

    auto g_0_yyyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 775);

    auto g_0_yyyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 776);

    auto g_0_yyyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 777);

    auto g_0_yyyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 778);

    auto g_0_yyyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 779);

    auto g_0_yyyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 780);

    auto g_0_yyyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 781);

    auto g_0_yyyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 782);

    auto g_0_yyyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 783);

    auto g_0_yyyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 784);

    auto g_0_yyyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 785);

    auto g_0_yyyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 786);

    auto g_0_yyyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 787);

    auto g_0_yyyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 788);

    auto g_0_yyyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 789);

    auto g_0_yyyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 790);

    auto g_0_yyyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 791);

    auto g_0_yyyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 793);

    auto g_0_yyyyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 794);

    auto g_0_yyyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 795);

    auto g_0_yyyyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 796);

    auto g_0_yyyyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 797);

    auto g_0_yyyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 798);

    auto g_0_yyyyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 799);

    auto g_0_yyyyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 800);

    auto g_0_yyyyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 801);

    auto g_0_yyyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 802);

    auto g_0_yyyyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 803);

    auto g_0_yyyyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 804);

    auto g_0_yyyyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 805);

    auto g_0_yyyyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 806);

    auto g_0_yyyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 807);

    auto g_0_yyyyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 808);

    auto g_0_yyyyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 809);

    auto g_0_yyyyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 810);

    auto g_0_yyyyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 811);

    auto g_0_yyyyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 812);

    auto g_0_yyyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 813);

    auto g_0_yyyyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 814);

    auto g_0_yyyyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 815);

    auto g_0_yyyyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 816);

    auto g_0_yyyyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 817);

    auto g_0_yyyyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 818);

    auto g_0_yyyyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 819);

    auto g_0_yyyyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 820);

    auto g_0_yyyyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 821);

    auto g_0_yyyyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 822);

    auto g_0_yyyyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 823);

    auto g_0_yyyyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 824);

    auto g_0_yyyyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 825);

    auto g_0_yyyyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 826);

    auto g_0_yyyyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 827);

    auto g_0_yyyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 828);

    auto g_0_yyyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 829);

    auto g_0_yyyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 830);

    auto g_0_yyyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 831);

    auto g_0_yyyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 832);

    auto g_0_yyyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 833);

    auto g_0_yyyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 834);

    auto g_0_yyyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 835);

    auto g_0_yyyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 836);

    auto g_0_yyyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 837);

    auto g_0_yyyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 838);

    auto g_0_yyyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 839);

    auto g_0_yyyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 840);

    auto g_0_yyyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 841);

    auto g_0_yyyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 842);

    auto g_0_yyyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 843);

    auto g_0_yyyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 844);

    auto g_0_yyyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 845);

    auto g_0_yyyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 846);

    auto g_0_yyyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 847);

    auto g_0_yyyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 848);

    auto g_0_yyyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 849);

    auto g_0_yyyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 850);

    auto g_0_yyyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 851);

    auto g_0_yyyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 852);

    auto g_0_yyyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 853);

    auto g_0_yyyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 854);

    auto g_0_yyyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 855);

    auto g_0_yyyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 856);

    auto g_0_yyyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 857);

    auto g_0_yyyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 858);

    auto g_0_yyyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 859);

    auto g_0_yyyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 860);

    auto g_0_yyyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 861);

    auto g_0_yyyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 862);

    auto g_0_yyyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 863);

    auto g_0_yyyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 864);

    auto g_0_yyyzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 865);

    auto g_0_yyyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 866);

    auto g_0_yyyzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 867);

    auto g_0_yyyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 868);

    auto g_0_yyyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 869);

    auto g_0_yyyzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 870);

    auto g_0_yyyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 871);

    auto g_0_yyyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 872);

    auto g_0_yyyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 873);

    auto g_0_yyyzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 874);

    auto g_0_yyyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 875);

    auto g_0_yyyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 876);

    auto g_0_yyyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 877);

    auto g_0_yyyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 878);

    auto g_0_yyyzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 879);

    auto g_0_yyyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 880);

    auto g_0_yyyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 881);

    auto g_0_yyyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 882);

    auto g_0_yyyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 883);

    auto g_0_yyyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 884);

    auto g_0_yyyzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 885);

    auto g_0_yyyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 886);

    auto g_0_yyyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 887);

    auto g_0_yyyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 888);

    auto g_0_yyyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 889);

    auto g_0_yyyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 890);

    auto g_0_yyyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 891);

    auto g_0_yyyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 892);

    auto g_0_yyyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 893);

    auto g_0_yyyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 894);

    auto g_0_yyyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 895);

    auto g_0_yyyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 896);

    auto g_0_yyyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 897);

    auto g_0_yyyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 898);

    auto g_0_yyyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 899);

    auto g_0_yyzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 900);

    auto g_0_yyzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 901);

    auto g_0_yyzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 902);

    auto g_0_yyzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 903);

    auto g_0_yyzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 904);

    auto g_0_yyzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 905);

    auto g_0_yyzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 906);

    auto g_0_yyzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 907);

    auto g_0_yyzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 908);

    auto g_0_yyzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 909);

    auto g_0_yyzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 910);

    auto g_0_yyzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 911);

    auto g_0_yyzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 912);

    auto g_0_yyzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 913);

    auto g_0_yyzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 914);

    auto g_0_yyzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 915);

    auto g_0_yyzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 916);

    auto g_0_yyzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 917);

    auto g_0_yyzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 918);

    auto g_0_yyzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 919);

    auto g_0_yyzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 920);

    auto g_0_yyzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 921);

    auto g_0_yyzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 922);

    auto g_0_yyzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 923);

    auto g_0_yyzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 924);

    auto g_0_yyzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 925);

    auto g_0_yyzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 926);

    auto g_0_yyzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 927);

    auto g_0_yyzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 928);

    auto g_0_yyzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 929);

    auto g_0_yyzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 930);

    auto g_0_yyzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 931);

    auto g_0_yyzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 932);

    auto g_0_yyzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 933);

    auto g_0_yyzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 934);

    auto g_0_yyzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 935);

    auto g_0_yzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 936);

    auto g_0_yzzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 937);

    auto g_0_yzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 938);

    auto g_0_yzzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 939);

    auto g_0_yzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 940);

    auto g_0_yzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 941);

    auto g_0_yzzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 942);

    auto g_0_yzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 943);

    auto g_0_yzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 944);

    auto g_0_yzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 945);

    auto g_0_yzzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 946);

    auto g_0_yzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 947);

    auto g_0_yzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 948);

    auto g_0_yzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 949);

    auto g_0_yzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 950);

    auto g_0_yzzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 951);

    auto g_0_yzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 952);

    auto g_0_yzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 953);

    auto g_0_yzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 954);

    auto g_0_yzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 955);

    auto g_0_yzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 956);

    auto g_0_yzzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 957);

    auto g_0_yzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 958);

    auto g_0_yzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 959);

    auto g_0_yzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 960);

    auto g_0_yzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 961);

    auto g_0_yzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 962);

    auto g_0_yzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 963);

    auto g_0_yzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 964);

    auto g_0_yzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 965);

    auto g_0_yzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 966);

    auto g_0_yzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 967);

    auto g_0_yzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 968);

    auto g_0_yzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 969);

    auto g_0_yzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 970);

    auto g_0_yzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 971);

    auto g_0_zzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 972);

    auto g_0_zzzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 973);

    auto g_0_zzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 974);

    auto g_0_zzzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 975);

    auto g_0_zzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 976);

    auto g_0_zzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 977);

    auto g_0_zzzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 978);

    auto g_0_zzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 979);

    auto g_0_zzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 980);

    auto g_0_zzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 981);

    auto g_0_zzzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 982);

    auto g_0_zzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 983);

    auto g_0_zzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 984);

    auto g_0_zzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 985);

    auto g_0_zzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 986);

    auto g_0_zzzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 987);

    auto g_0_zzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 988);

    auto g_0_zzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 989);

    auto g_0_zzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 990);

    auto g_0_zzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 991);

    auto g_0_zzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 992);

    auto g_0_zzzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 993);

    auto g_0_zzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 994);

    auto g_0_zzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 995);

    auto g_0_zzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 996);

    auto g_0_zzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 997);

    auto g_0_zzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 998);

    auto g_0_zzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 999);

    auto g_0_zzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 1000);

    auto g_0_zzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 1001);

    auto g_0_zzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 1002);

    auto g_0_zzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 1003);

    auto g_0_zzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 1004);

    auto g_0_zzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 1005);

    auto g_0_zzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 1006);

    auto g_0_zzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 1007);

    /// Set up 0-36 components of targeted buffer : SKSK

    auto g_0_xxxxxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk);

    auto g_0_xxxxxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1);

    auto g_0_xxxxxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 2);

    auto g_0_xxxxxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 3);

    auto g_0_xxxxxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 4);

    auto g_0_xxxxxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 5);

    auto g_0_xxxxxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 6);

    auto g_0_xxxxxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 7);

    auto g_0_xxxxxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 8);

    auto g_0_xxxxxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 9);

    auto g_0_xxxxxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 10);

    auto g_0_xxxxxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 11);

    auto g_0_xxxxxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 12);

    auto g_0_xxxxxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 13);

    auto g_0_xxxxxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 14);

    auto g_0_xxxxxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 15);

    auto g_0_xxxxxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 16);

    auto g_0_xxxxxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 17);

    auto g_0_xxxxxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 18);

    auto g_0_xxxxxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 19);

    auto g_0_xxxxxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 20);

    auto g_0_xxxxxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 21);

    auto g_0_xxxxxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 22);

    auto g_0_xxxxxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 23);

    auto g_0_xxxxxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 24);

    auto g_0_xxxxxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 25);

    auto g_0_xxxxxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 26);

    auto g_0_xxxxxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 27);

    auto g_0_xxxxxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 28);

    auto g_0_xxxxxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 29);

    auto g_0_xxxxxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 30);

    auto g_0_xxxxxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 31);

    auto g_0_xxxxxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 32);

    auto g_0_xxxxxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 33);

    auto g_0_xxxxxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 34);

    auto g_0_xxxxxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 35);

    #pragma omp simd aligned(g_0_xxxxx_0_xxxxxxx_0, g_0_xxxxx_0_xxxxxxx_1, g_0_xxxxx_0_xxxxxxy_0, g_0_xxxxx_0_xxxxxxy_1, g_0_xxxxx_0_xxxxxxz_0, g_0_xxxxx_0_xxxxxxz_1, g_0_xxxxx_0_xxxxxyy_0, g_0_xxxxx_0_xxxxxyy_1, g_0_xxxxx_0_xxxxxyz_0, g_0_xxxxx_0_xxxxxyz_1, g_0_xxxxx_0_xxxxxzz_0, g_0_xxxxx_0_xxxxxzz_1, g_0_xxxxx_0_xxxxyyy_0, g_0_xxxxx_0_xxxxyyy_1, g_0_xxxxx_0_xxxxyyz_0, g_0_xxxxx_0_xxxxyyz_1, g_0_xxxxx_0_xxxxyzz_0, g_0_xxxxx_0_xxxxyzz_1, g_0_xxxxx_0_xxxxzzz_0, g_0_xxxxx_0_xxxxzzz_1, g_0_xxxxx_0_xxxyyyy_0, g_0_xxxxx_0_xxxyyyy_1, g_0_xxxxx_0_xxxyyyz_0, g_0_xxxxx_0_xxxyyyz_1, g_0_xxxxx_0_xxxyyzz_0, g_0_xxxxx_0_xxxyyzz_1, g_0_xxxxx_0_xxxyzzz_0, g_0_xxxxx_0_xxxyzzz_1, g_0_xxxxx_0_xxxzzzz_0, g_0_xxxxx_0_xxxzzzz_1, g_0_xxxxx_0_xxyyyyy_0, g_0_xxxxx_0_xxyyyyy_1, g_0_xxxxx_0_xxyyyyz_0, g_0_xxxxx_0_xxyyyyz_1, g_0_xxxxx_0_xxyyyzz_0, g_0_xxxxx_0_xxyyyzz_1, g_0_xxxxx_0_xxyyzzz_0, g_0_xxxxx_0_xxyyzzz_1, g_0_xxxxx_0_xxyzzzz_0, g_0_xxxxx_0_xxyzzzz_1, g_0_xxxxx_0_xxzzzzz_0, g_0_xxxxx_0_xxzzzzz_1, g_0_xxxxx_0_xyyyyyy_0, g_0_xxxxx_0_xyyyyyy_1, g_0_xxxxx_0_xyyyyyz_0, g_0_xxxxx_0_xyyyyyz_1, g_0_xxxxx_0_xyyyyzz_0, g_0_xxxxx_0_xyyyyzz_1, g_0_xxxxx_0_xyyyzzz_0, g_0_xxxxx_0_xyyyzzz_1, g_0_xxxxx_0_xyyzzzz_0, g_0_xxxxx_0_xyyzzzz_1, g_0_xxxxx_0_xyzzzzz_0, g_0_xxxxx_0_xyzzzzz_1, g_0_xxxxx_0_xzzzzzz_0, g_0_xxxxx_0_xzzzzzz_1, g_0_xxxxx_0_yyyyyyy_0, g_0_xxxxx_0_yyyyyyy_1, g_0_xxxxx_0_yyyyyyz_0, g_0_xxxxx_0_yyyyyyz_1, g_0_xxxxx_0_yyyyyzz_0, g_0_xxxxx_0_yyyyyzz_1, g_0_xxxxx_0_yyyyzzz_0, g_0_xxxxx_0_yyyyzzz_1, g_0_xxxxx_0_yyyzzzz_0, g_0_xxxxx_0_yyyzzzz_1, g_0_xxxxx_0_yyzzzzz_0, g_0_xxxxx_0_yyzzzzz_1, g_0_xxxxx_0_yzzzzzz_0, g_0_xxxxx_0_yzzzzzz_1, g_0_xxxxx_0_zzzzzzz_0, g_0_xxxxx_0_zzzzzzz_1, g_0_xxxxxx_0_xxxxxx_1, g_0_xxxxxx_0_xxxxxxx_0, g_0_xxxxxx_0_xxxxxxx_1, g_0_xxxxxx_0_xxxxxxy_0, g_0_xxxxxx_0_xxxxxxy_1, g_0_xxxxxx_0_xxxxxxz_0, g_0_xxxxxx_0_xxxxxxz_1, g_0_xxxxxx_0_xxxxxy_1, g_0_xxxxxx_0_xxxxxyy_0, g_0_xxxxxx_0_xxxxxyy_1, g_0_xxxxxx_0_xxxxxyz_0, g_0_xxxxxx_0_xxxxxyz_1, g_0_xxxxxx_0_xxxxxz_1, g_0_xxxxxx_0_xxxxxzz_0, g_0_xxxxxx_0_xxxxxzz_1, g_0_xxxxxx_0_xxxxyy_1, g_0_xxxxxx_0_xxxxyyy_0, g_0_xxxxxx_0_xxxxyyy_1, g_0_xxxxxx_0_xxxxyyz_0, g_0_xxxxxx_0_xxxxyyz_1, g_0_xxxxxx_0_xxxxyz_1, g_0_xxxxxx_0_xxxxyzz_0, g_0_xxxxxx_0_xxxxyzz_1, g_0_xxxxxx_0_xxxxzz_1, g_0_xxxxxx_0_xxxxzzz_0, g_0_xxxxxx_0_xxxxzzz_1, g_0_xxxxxx_0_xxxyyy_1, g_0_xxxxxx_0_xxxyyyy_0, g_0_xxxxxx_0_xxxyyyy_1, g_0_xxxxxx_0_xxxyyyz_0, g_0_xxxxxx_0_xxxyyyz_1, g_0_xxxxxx_0_xxxyyz_1, g_0_xxxxxx_0_xxxyyzz_0, g_0_xxxxxx_0_xxxyyzz_1, g_0_xxxxxx_0_xxxyzz_1, g_0_xxxxxx_0_xxxyzzz_0, g_0_xxxxxx_0_xxxyzzz_1, g_0_xxxxxx_0_xxxzzz_1, g_0_xxxxxx_0_xxxzzzz_0, g_0_xxxxxx_0_xxxzzzz_1, g_0_xxxxxx_0_xxyyyy_1, g_0_xxxxxx_0_xxyyyyy_0, g_0_xxxxxx_0_xxyyyyy_1, g_0_xxxxxx_0_xxyyyyz_0, g_0_xxxxxx_0_xxyyyyz_1, g_0_xxxxxx_0_xxyyyz_1, g_0_xxxxxx_0_xxyyyzz_0, g_0_xxxxxx_0_xxyyyzz_1, g_0_xxxxxx_0_xxyyzz_1, g_0_xxxxxx_0_xxyyzzz_0, g_0_xxxxxx_0_xxyyzzz_1, g_0_xxxxxx_0_xxyzzz_1, g_0_xxxxxx_0_xxyzzzz_0, g_0_xxxxxx_0_xxyzzzz_1, g_0_xxxxxx_0_xxzzzz_1, g_0_xxxxxx_0_xxzzzzz_0, g_0_xxxxxx_0_xxzzzzz_1, g_0_xxxxxx_0_xyyyyy_1, g_0_xxxxxx_0_xyyyyyy_0, g_0_xxxxxx_0_xyyyyyy_1, g_0_xxxxxx_0_xyyyyyz_0, g_0_xxxxxx_0_xyyyyyz_1, g_0_xxxxxx_0_xyyyyz_1, g_0_xxxxxx_0_xyyyyzz_0, g_0_xxxxxx_0_xyyyyzz_1, g_0_xxxxxx_0_xyyyzz_1, g_0_xxxxxx_0_xyyyzzz_0, g_0_xxxxxx_0_xyyyzzz_1, g_0_xxxxxx_0_xyyzzz_1, g_0_xxxxxx_0_xyyzzzz_0, g_0_xxxxxx_0_xyyzzzz_1, g_0_xxxxxx_0_xyzzzz_1, g_0_xxxxxx_0_xyzzzzz_0, g_0_xxxxxx_0_xyzzzzz_1, g_0_xxxxxx_0_xzzzzz_1, g_0_xxxxxx_0_xzzzzzz_0, g_0_xxxxxx_0_xzzzzzz_1, g_0_xxxxxx_0_yyyyyy_1, g_0_xxxxxx_0_yyyyyyy_0, g_0_xxxxxx_0_yyyyyyy_1, g_0_xxxxxx_0_yyyyyyz_0, g_0_xxxxxx_0_yyyyyyz_1, g_0_xxxxxx_0_yyyyyz_1, g_0_xxxxxx_0_yyyyyzz_0, g_0_xxxxxx_0_yyyyyzz_1, g_0_xxxxxx_0_yyyyzz_1, g_0_xxxxxx_0_yyyyzzz_0, g_0_xxxxxx_0_yyyyzzz_1, g_0_xxxxxx_0_yyyzzz_1, g_0_xxxxxx_0_yyyzzzz_0, g_0_xxxxxx_0_yyyzzzz_1, g_0_xxxxxx_0_yyzzzz_1, g_0_xxxxxx_0_yyzzzzz_0, g_0_xxxxxx_0_yyzzzzz_1, g_0_xxxxxx_0_yzzzzz_1, g_0_xxxxxx_0_yzzzzzz_0, g_0_xxxxxx_0_yzzzzzz_1, g_0_xxxxxx_0_zzzzzz_1, g_0_xxxxxx_0_zzzzzzz_0, g_0_xxxxxx_0_zzzzzzz_1, g_0_xxxxxxx_0_xxxxxxx_0, g_0_xxxxxxx_0_xxxxxxy_0, g_0_xxxxxxx_0_xxxxxxz_0, g_0_xxxxxxx_0_xxxxxyy_0, g_0_xxxxxxx_0_xxxxxyz_0, g_0_xxxxxxx_0_xxxxxzz_0, g_0_xxxxxxx_0_xxxxyyy_0, g_0_xxxxxxx_0_xxxxyyz_0, g_0_xxxxxxx_0_xxxxyzz_0, g_0_xxxxxxx_0_xxxxzzz_0, g_0_xxxxxxx_0_xxxyyyy_0, g_0_xxxxxxx_0_xxxyyyz_0, g_0_xxxxxxx_0_xxxyyzz_0, g_0_xxxxxxx_0_xxxyzzz_0, g_0_xxxxxxx_0_xxxzzzz_0, g_0_xxxxxxx_0_xxyyyyy_0, g_0_xxxxxxx_0_xxyyyyz_0, g_0_xxxxxxx_0_xxyyyzz_0, g_0_xxxxxxx_0_xxyyzzz_0, g_0_xxxxxxx_0_xxyzzzz_0, g_0_xxxxxxx_0_xxzzzzz_0, g_0_xxxxxxx_0_xyyyyyy_0, g_0_xxxxxxx_0_xyyyyyz_0, g_0_xxxxxxx_0_xyyyyzz_0, g_0_xxxxxxx_0_xyyyzzz_0, g_0_xxxxxxx_0_xyyzzzz_0, g_0_xxxxxxx_0_xyzzzzz_0, g_0_xxxxxxx_0_xzzzzzz_0, g_0_xxxxxxx_0_yyyyyyy_0, g_0_xxxxxxx_0_yyyyyyz_0, g_0_xxxxxxx_0_yyyyyzz_0, g_0_xxxxxxx_0_yyyyzzz_0, g_0_xxxxxxx_0_yyyzzzz_0, g_0_xxxxxxx_0_yyzzzzz_0, g_0_xxxxxxx_0_yzzzzzz_0, g_0_xxxxxxx_0_zzzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_xxxxxxx_0[i] = 6.0 * g_0_xxxxx_0_xxxxxxx_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxxx_1[i] * fti_ab_0 + 7.0 * g_0_xxxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxxx_0[i] * pb_x + g_0_xxxxxx_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxxxy_0[i] = 6.0 * g_0_xxxxx_0_xxxxxxy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxxy_0[i] * pb_x + g_0_xxxxxx_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxxxz_0[i] = 6.0 * g_0_xxxxx_0_xxxxxxz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxxz_0[i] * pb_x + g_0_xxxxxx_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxxyy_0[i] = 6.0 * g_0_xxxxx_0_xxxxxyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxyy_0[i] * pb_x + g_0_xxxxxx_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxxyz_0[i] = 6.0 * g_0_xxxxx_0_xxxxxyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxyz_0[i] * pb_x + g_0_xxxxxx_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxxzz_0[i] = 6.0 * g_0_xxxxx_0_xxxxxzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxzz_0[i] * pb_x + g_0_xxxxxx_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxyyy_0[i] = 6.0 * g_0_xxxxx_0_xxxxyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyyy_0[i] * pb_x + g_0_xxxxxx_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxyyz_0[i] = 6.0 * g_0_xxxxx_0_xxxxyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyyz_0[i] * pb_x + g_0_xxxxxx_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxyzz_0[i] = 6.0 * g_0_xxxxx_0_xxxxyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyzz_0[i] * pb_x + g_0_xxxxxx_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxzzz_0[i] = 6.0 * g_0_xxxxx_0_xxxxzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxzzz_0[i] * pb_x + g_0_xxxxxx_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyyyy_0[i] = 6.0 * g_0_xxxxx_0_xxxyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyyy_0[i] * pb_x + g_0_xxxxxx_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyyyz_0[i] = 6.0 * g_0_xxxxx_0_xxxyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyyz_0[i] * pb_x + g_0_xxxxxx_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyyzz_0[i] = 6.0 * g_0_xxxxx_0_xxxyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyzz_0[i] * pb_x + g_0_xxxxxx_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyzzz_0[i] = 6.0 * g_0_xxxxx_0_xxxyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyzzz_0[i] * pb_x + g_0_xxxxxx_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxzzzz_0[i] = 6.0 * g_0_xxxxx_0_xxxzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxzzzz_0[i] * pb_x + g_0_xxxxxx_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyyyy_0[i] = 6.0 * g_0_xxxxx_0_xxyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyyy_0[i] * pb_x + g_0_xxxxxx_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyyyz_0[i] = 6.0 * g_0_xxxxx_0_xxyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyyz_0[i] * pb_x + g_0_xxxxxx_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyyzz_0[i] = 6.0 * g_0_xxxxx_0_xxyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyzz_0[i] * pb_x + g_0_xxxxxx_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyzzz_0[i] = 6.0 * g_0_xxxxx_0_xxyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyzzz_0[i] * pb_x + g_0_xxxxxx_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyzzzz_0[i] = 6.0 * g_0_xxxxx_0_xxyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzzzz_0[i] * pb_x + g_0_xxxxxx_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxzzzzz_0[i] = 6.0 * g_0_xxxxx_0_xxzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxzzzzz_0[i] * pb_x + g_0_xxxxxx_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyyyy_0[i] = 6.0 * g_0_xxxxx_0_xyyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyyy_0[i] * pb_x + g_0_xxxxxx_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyyyz_0[i] = 6.0 * g_0_xxxxx_0_xyyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyyz_0[i] * pb_x + g_0_xxxxxx_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyyzz_0[i] = 6.0 * g_0_xxxxx_0_xyyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyzz_0[i] * pb_x + g_0_xxxxxx_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyzzz_0[i] = 6.0 * g_0_xxxxx_0_xyyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyzzz_0[i] * pb_x + g_0_xxxxxx_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyzzzz_0[i] = 6.0 * g_0_xxxxx_0_xyyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzzzz_0[i] * pb_x + g_0_xxxxxx_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyzzzzz_0[i] = 6.0 * g_0_xxxxx_0_xyzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzzzz_0[i] * pb_x + g_0_xxxxxx_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xzzzzzz_0[i] = 6.0 * g_0_xxxxx_0_xzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xzzzzzz_0[i] * pb_x + g_0_xxxxxx_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyyyy_0[i] = 6.0 * g_0_xxxxx_0_yyyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyyyy_0[i] * pb_x + g_0_xxxxxx_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxx_0_yyyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyyyz_0[i] * pb_x + g_0_xxxxxx_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyyzz_0[i] = 6.0 * g_0_xxxxx_0_yyyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyyzz_0[i] * pb_x + g_0_xxxxxx_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyzzz_0[i] = 6.0 * g_0_xxxxx_0_yyyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyyzzz_0[i] * pb_x + g_0_xxxxxx_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyzzzz_0[i] = 6.0 * g_0_xxxxx_0_yyyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyyzzzz_0[i] * pb_x + g_0_xxxxxx_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyzzzzz_0[i] = 6.0 * g_0_xxxxx_0_yyzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yyzzzzz_0[i] * pb_x + g_0_xxxxxx_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yzzzzzz_0[i] = 6.0 * g_0_xxxxx_0_yzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_yzzzzzz_0[i] * pb_x + g_0_xxxxxx_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_zzzzzzz_0[i] = 6.0 * g_0_xxxxx_0_zzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxx_0_zzzzzzz_0[i] * pb_x + g_0_xxxxxx_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 36-72 components of targeted buffer : SKSK

    auto g_0_xxxxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 36);

    auto g_0_xxxxxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 37);

    auto g_0_xxxxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 38);

    auto g_0_xxxxxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 39);

    auto g_0_xxxxxxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 40);

    auto g_0_xxxxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 41);

    auto g_0_xxxxxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 42);

    auto g_0_xxxxxxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 43);

    auto g_0_xxxxxxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 44);

    auto g_0_xxxxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 45);

    auto g_0_xxxxxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 46);

    auto g_0_xxxxxxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 47);

    auto g_0_xxxxxxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 48);

    auto g_0_xxxxxxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 49);

    auto g_0_xxxxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 50);

    auto g_0_xxxxxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 51);

    auto g_0_xxxxxxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 52);

    auto g_0_xxxxxxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 53);

    auto g_0_xxxxxxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 54);

    auto g_0_xxxxxxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 55);

    auto g_0_xxxxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 56);

    auto g_0_xxxxxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 57);

    auto g_0_xxxxxxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 58);

    auto g_0_xxxxxxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 59);

    auto g_0_xxxxxxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 60);

    auto g_0_xxxxxxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 61);

    auto g_0_xxxxxxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 62);

    auto g_0_xxxxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 63);

    auto g_0_xxxxxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 64);

    auto g_0_xxxxxxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 65);

    auto g_0_xxxxxxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 66);

    auto g_0_xxxxxxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 67);

    auto g_0_xxxxxxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 68);

    auto g_0_xxxxxxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 69);

    auto g_0_xxxxxxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 70);

    auto g_0_xxxxxxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 71);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxxxx_1, g_0_xxxxxx_0_xxxxxxx_0, g_0_xxxxxx_0_xxxxxxx_1, g_0_xxxxxx_0_xxxxxxy_0, g_0_xxxxxx_0_xxxxxxy_1, g_0_xxxxxx_0_xxxxxxz_0, g_0_xxxxxx_0_xxxxxxz_1, g_0_xxxxxx_0_xxxxxy_1, g_0_xxxxxx_0_xxxxxyy_0, g_0_xxxxxx_0_xxxxxyy_1, g_0_xxxxxx_0_xxxxxyz_0, g_0_xxxxxx_0_xxxxxyz_1, g_0_xxxxxx_0_xxxxxz_1, g_0_xxxxxx_0_xxxxxzz_0, g_0_xxxxxx_0_xxxxxzz_1, g_0_xxxxxx_0_xxxxyy_1, g_0_xxxxxx_0_xxxxyyy_0, g_0_xxxxxx_0_xxxxyyy_1, g_0_xxxxxx_0_xxxxyyz_0, g_0_xxxxxx_0_xxxxyyz_1, g_0_xxxxxx_0_xxxxyz_1, g_0_xxxxxx_0_xxxxyzz_0, g_0_xxxxxx_0_xxxxyzz_1, g_0_xxxxxx_0_xxxxzz_1, g_0_xxxxxx_0_xxxxzzz_0, g_0_xxxxxx_0_xxxxzzz_1, g_0_xxxxxx_0_xxxyyy_1, g_0_xxxxxx_0_xxxyyyy_0, g_0_xxxxxx_0_xxxyyyy_1, g_0_xxxxxx_0_xxxyyyz_0, g_0_xxxxxx_0_xxxyyyz_1, g_0_xxxxxx_0_xxxyyz_1, g_0_xxxxxx_0_xxxyyzz_0, g_0_xxxxxx_0_xxxyyzz_1, g_0_xxxxxx_0_xxxyzz_1, g_0_xxxxxx_0_xxxyzzz_0, g_0_xxxxxx_0_xxxyzzz_1, g_0_xxxxxx_0_xxxzzz_1, g_0_xxxxxx_0_xxxzzzz_0, g_0_xxxxxx_0_xxxzzzz_1, g_0_xxxxxx_0_xxyyyy_1, g_0_xxxxxx_0_xxyyyyy_0, g_0_xxxxxx_0_xxyyyyy_1, g_0_xxxxxx_0_xxyyyyz_0, g_0_xxxxxx_0_xxyyyyz_1, g_0_xxxxxx_0_xxyyyz_1, g_0_xxxxxx_0_xxyyyzz_0, g_0_xxxxxx_0_xxyyyzz_1, g_0_xxxxxx_0_xxyyzz_1, g_0_xxxxxx_0_xxyyzzz_0, g_0_xxxxxx_0_xxyyzzz_1, g_0_xxxxxx_0_xxyzzz_1, g_0_xxxxxx_0_xxyzzzz_0, g_0_xxxxxx_0_xxyzzzz_1, g_0_xxxxxx_0_xxzzzz_1, g_0_xxxxxx_0_xxzzzzz_0, g_0_xxxxxx_0_xxzzzzz_1, g_0_xxxxxx_0_xyyyyy_1, g_0_xxxxxx_0_xyyyyyy_0, g_0_xxxxxx_0_xyyyyyy_1, g_0_xxxxxx_0_xyyyyyz_0, g_0_xxxxxx_0_xyyyyyz_1, g_0_xxxxxx_0_xyyyyz_1, g_0_xxxxxx_0_xyyyyzz_0, g_0_xxxxxx_0_xyyyyzz_1, g_0_xxxxxx_0_xyyyzz_1, g_0_xxxxxx_0_xyyyzzz_0, g_0_xxxxxx_0_xyyyzzz_1, g_0_xxxxxx_0_xyyzzz_1, g_0_xxxxxx_0_xyyzzzz_0, g_0_xxxxxx_0_xyyzzzz_1, g_0_xxxxxx_0_xyzzzz_1, g_0_xxxxxx_0_xyzzzzz_0, g_0_xxxxxx_0_xyzzzzz_1, g_0_xxxxxx_0_xzzzzz_1, g_0_xxxxxx_0_xzzzzzz_0, g_0_xxxxxx_0_xzzzzzz_1, g_0_xxxxxx_0_yyyyyy_1, g_0_xxxxxx_0_yyyyyyy_0, g_0_xxxxxx_0_yyyyyyy_1, g_0_xxxxxx_0_yyyyyyz_0, g_0_xxxxxx_0_yyyyyyz_1, g_0_xxxxxx_0_yyyyyz_1, g_0_xxxxxx_0_yyyyyzz_0, g_0_xxxxxx_0_yyyyyzz_1, g_0_xxxxxx_0_yyyyzz_1, g_0_xxxxxx_0_yyyyzzz_0, g_0_xxxxxx_0_yyyyzzz_1, g_0_xxxxxx_0_yyyzzz_1, g_0_xxxxxx_0_yyyzzzz_0, g_0_xxxxxx_0_yyyzzzz_1, g_0_xxxxxx_0_yyzzzz_1, g_0_xxxxxx_0_yyzzzzz_0, g_0_xxxxxx_0_yyzzzzz_1, g_0_xxxxxx_0_yzzzzz_1, g_0_xxxxxx_0_yzzzzzz_0, g_0_xxxxxx_0_yzzzzzz_1, g_0_xxxxxx_0_zzzzzz_1, g_0_xxxxxx_0_zzzzzzz_0, g_0_xxxxxx_0_zzzzzzz_1, g_0_xxxxxxy_0_xxxxxxx_0, g_0_xxxxxxy_0_xxxxxxy_0, g_0_xxxxxxy_0_xxxxxxz_0, g_0_xxxxxxy_0_xxxxxyy_0, g_0_xxxxxxy_0_xxxxxyz_0, g_0_xxxxxxy_0_xxxxxzz_0, g_0_xxxxxxy_0_xxxxyyy_0, g_0_xxxxxxy_0_xxxxyyz_0, g_0_xxxxxxy_0_xxxxyzz_0, g_0_xxxxxxy_0_xxxxzzz_0, g_0_xxxxxxy_0_xxxyyyy_0, g_0_xxxxxxy_0_xxxyyyz_0, g_0_xxxxxxy_0_xxxyyzz_0, g_0_xxxxxxy_0_xxxyzzz_0, g_0_xxxxxxy_0_xxxzzzz_0, g_0_xxxxxxy_0_xxyyyyy_0, g_0_xxxxxxy_0_xxyyyyz_0, g_0_xxxxxxy_0_xxyyyzz_0, g_0_xxxxxxy_0_xxyyzzz_0, g_0_xxxxxxy_0_xxyzzzz_0, g_0_xxxxxxy_0_xxzzzzz_0, g_0_xxxxxxy_0_xyyyyyy_0, g_0_xxxxxxy_0_xyyyyyz_0, g_0_xxxxxxy_0_xyyyyzz_0, g_0_xxxxxxy_0_xyyyzzz_0, g_0_xxxxxxy_0_xyyzzzz_0, g_0_xxxxxxy_0_xyzzzzz_0, g_0_xxxxxxy_0_xzzzzzz_0, g_0_xxxxxxy_0_yyyyyyy_0, g_0_xxxxxxy_0_yyyyyyz_0, g_0_xxxxxxy_0_yyyyyzz_0, g_0_xxxxxxy_0_yyyyzzz_0, g_0_xxxxxxy_0_yyyzzzz_0, g_0_xxxxxxy_0_yyzzzzz_0, g_0_xxxxxxy_0_yzzzzzz_0, g_0_xxxxxxy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxy_0_xxxxxxx_0[i] = g_0_xxxxxx_0_xxxxxxx_0[i] * pb_y + g_0_xxxxxx_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxxxy_0[i] = g_0_xxxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxxy_0[i] * pb_y + g_0_xxxxxx_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxxxz_0[i] = g_0_xxxxxx_0_xxxxxxz_0[i] * pb_y + g_0_xxxxxx_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxxyy_0[i] = 2.0 * g_0_xxxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxyy_0[i] * pb_y + g_0_xxxxxx_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxxyz_0[i] = g_0_xxxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxyz_0[i] * pb_y + g_0_xxxxxx_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxxzz_0[i] = g_0_xxxxxx_0_xxxxxzz_0[i] * pb_y + g_0_xxxxxx_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxyyy_0[i] = 3.0 * g_0_xxxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyyy_0[i] * pb_y + g_0_xxxxxx_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyyz_0[i] * pb_y + g_0_xxxxxx_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxyzz_0[i] = g_0_xxxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyzz_0[i] * pb_y + g_0_xxxxxx_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxzzz_0[i] = g_0_xxxxxx_0_xxxxzzz_0[i] * pb_y + g_0_xxxxxx_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyyyy_0[i] = 4.0 * g_0_xxxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyyy_0[i] * pb_y + g_0_xxxxxx_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyyz_0[i] * pb_y + g_0_xxxxxx_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyzz_0[i] * pb_y + g_0_xxxxxx_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyzzz_0[i] = g_0_xxxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyzzz_0[i] * pb_y + g_0_xxxxxx_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxzzzz_0[i] = g_0_xxxxxx_0_xxxzzzz_0[i] * pb_y + g_0_xxxxxx_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyyyy_0[i] = 5.0 * g_0_xxxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyyy_0[i] * pb_y + g_0_xxxxxx_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyyz_0[i] * pb_y + g_0_xxxxxx_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyzz_0[i] * pb_y + g_0_xxxxxx_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyzzz_0[i] * pb_y + g_0_xxxxxx_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyzzzz_0[i] = g_0_xxxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzzzz_0[i] * pb_y + g_0_xxxxxx_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxzzzzz_0[i] = g_0_xxxxxx_0_xxzzzzz_0[i] * pb_y + g_0_xxxxxx_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyyyy_0[i] = 6.0 * g_0_xxxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyyy_0[i] * pb_y + g_0_xxxxxx_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyyz_0[i] * pb_y + g_0_xxxxxx_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyzz_0[i] * pb_y + g_0_xxxxxx_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyzzz_0[i] * pb_y + g_0_xxxxxx_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzzzz_0[i] * pb_y + g_0_xxxxxx_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyzzzzz_0[i] = g_0_xxxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzzzz_0[i] * pb_y + g_0_xxxxxx_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xzzzzzz_0[i] = g_0_xxxxxx_0_xzzzzzz_0[i] * pb_y + g_0_xxxxxx_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyyyy_0[i] = 7.0 * g_0_xxxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyyyy_0[i] * pb_y + g_0_xxxxxx_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyyyz_0[i] * pb_y + g_0_xxxxxx_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyyzz_0[i] * pb_y + g_0_xxxxxx_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyzzz_0[i] * pb_y + g_0_xxxxxx_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyzzzz_0[i] * pb_y + g_0_xxxxxx_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyzzzzz_0[i] * pb_y + g_0_xxxxxx_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yzzzzzz_0[i] = g_0_xxxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzzzzzz_0[i] * pb_y + g_0_xxxxxx_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_zzzzzzz_0[i] = g_0_xxxxxx_0_zzzzzzz_0[i] * pb_y + g_0_xxxxxx_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 72-108 components of targeted buffer : SKSK

    auto g_0_xxxxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 72);

    auto g_0_xxxxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 73);

    auto g_0_xxxxxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 74);

    auto g_0_xxxxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 75);

    auto g_0_xxxxxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 76);

    auto g_0_xxxxxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 77);

    auto g_0_xxxxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 78);

    auto g_0_xxxxxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 79);

    auto g_0_xxxxxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 80);

    auto g_0_xxxxxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 81);

    auto g_0_xxxxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 82);

    auto g_0_xxxxxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 83);

    auto g_0_xxxxxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 84);

    auto g_0_xxxxxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 85);

    auto g_0_xxxxxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 86);

    auto g_0_xxxxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 87);

    auto g_0_xxxxxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 88);

    auto g_0_xxxxxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 89);

    auto g_0_xxxxxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 90);

    auto g_0_xxxxxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 91);

    auto g_0_xxxxxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 92);

    auto g_0_xxxxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 93);

    auto g_0_xxxxxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 94);

    auto g_0_xxxxxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 95);

    auto g_0_xxxxxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 96);

    auto g_0_xxxxxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 97);

    auto g_0_xxxxxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 98);

    auto g_0_xxxxxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 99);

    auto g_0_xxxxxxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 100);

    auto g_0_xxxxxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 101);

    auto g_0_xxxxxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 102);

    auto g_0_xxxxxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 103);

    auto g_0_xxxxxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 104);

    auto g_0_xxxxxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 105);

    auto g_0_xxxxxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 106);

    auto g_0_xxxxxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 107);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxxxx_1, g_0_xxxxxx_0_xxxxxxx_0, g_0_xxxxxx_0_xxxxxxx_1, g_0_xxxxxx_0_xxxxxxy_0, g_0_xxxxxx_0_xxxxxxy_1, g_0_xxxxxx_0_xxxxxxz_0, g_0_xxxxxx_0_xxxxxxz_1, g_0_xxxxxx_0_xxxxxy_1, g_0_xxxxxx_0_xxxxxyy_0, g_0_xxxxxx_0_xxxxxyy_1, g_0_xxxxxx_0_xxxxxyz_0, g_0_xxxxxx_0_xxxxxyz_1, g_0_xxxxxx_0_xxxxxz_1, g_0_xxxxxx_0_xxxxxzz_0, g_0_xxxxxx_0_xxxxxzz_1, g_0_xxxxxx_0_xxxxyy_1, g_0_xxxxxx_0_xxxxyyy_0, g_0_xxxxxx_0_xxxxyyy_1, g_0_xxxxxx_0_xxxxyyz_0, g_0_xxxxxx_0_xxxxyyz_1, g_0_xxxxxx_0_xxxxyz_1, g_0_xxxxxx_0_xxxxyzz_0, g_0_xxxxxx_0_xxxxyzz_1, g_0_xxxxxx_0_xxxxzz_1, g_0_xxxxxx_0_xxxxzzz_0, g_0_xxxxxx_0_xxxxzzz_1, g_0_xxxxxx_0_xxxyyy_1, g_0_xxxxxx_0_xxxyyyy_0, g_0_xxxxxx_0_xxxyyyy_1, g_0_xxxxxx_0_xxxyyyz_0, g_0_xxxxxx_0_xxxyyyz_1, g_0_xxxxxx_0_xxxyyz_1, g_0_xxxxxx_0_xxxyyzz_0, g_0_xxxxxx_0_xxxyyzz_1, g_0_xxxxxx_0_xxxyzz_1, g_0_xxxxxx_0_xxxyzzz_0, g_0_xxxxxx_0_xxxyzzz_1, g_0_xxxxxx_0_xxxzzz_1, g_0_xxxxxx_0_xxxzzzz_0, g_0_xxxxxx_0_xxxzzzz_1, g_0_xxxxxx_0_xxyyyy_1, g_0_xxxxxx_0_xxyyyyy_0, g_0_xxxxxx_0_xxyyyyy_1, g_0_xxxxxx_0_xxyyyyz_0, g_0_xxxxxx_0_xxyyyyz_1, g_0_xxxxxx_0_xxyyyz_1, g_0_xxxxxx_0_xxyyyzz_0, g_0_xxxxxx_0_xxyyyzz_1, g_0_xxxxxx_0_xxyyzz_1, g_0_xxxxxx_0_xxyyzzz_0, g_0_xxxxxx_0_xxyyzzz_1, g_0_xxxxxx_0_xxyzzz_1, g_0_xxxxxx_0_xxyzzzz_0, g_0_xxxxxx_0_xxyzzzz_1, g_0_xxxxxx_0_xxzzzz_1, g_0_xxxxxx_0_xxzzzzz_0, g_0_xxxxxx_0_xxzzzzz_1, g_0_xxxxxx_0_xyyyyy_1, g_0_xxxxxx_0_xyyyyyy_0, g_0_xxxxxx_0_xyyyyyy_1, g_0_xxxxxx_0_xyyyyyz_0, g_0_xxxxxx_0_xyyyyyz_1, g_0_xxxxxx_0_xyyyyz_1, g_0_xxxxxx_0_xyyyyzz_0, g_0_xxxxxx_0_xyyyyzz_1, g_0_xxxxxx_0_xyyyzz_1, g_0_xxxxxx_0_xyyyzzz_0, g_0_xxxxxx_0_xyyyzzz_1, g_0_xxxxxx_0_xyyzzz_1, g_0_xxxxxx_0_xyyzzzz_0, g_0_xxxxxx_0_xyyzzzz_1, g_0_xxxxxx_0_xyzzzz_1, g_0_xxxxxx_0_xyzzzzz_0, g_0_xxxxxx_0_xyzzzzz_1, g_0_xxxxxx_0_xzzzzz_1, g_0_xxxxxx_0_xzzzzzz_0, g_0_xxxxxx_0_xzzzzzz_1, g_0_xxxxxx_0_yyyyyy_1, g_0_xxxxxx_0_yyyyyyy_0, g_0_xxxxxx_0_yyyyyyy_1, g_0_xxxxxx_0_yyyyyyz_0, g_0_xxxxxx_0_yyyyyyz_1, g_0_xxxxxx_0_yyyyyz_1, g_0_xxxxxx_0_yyyyyzz_0, g_0_xxxxxx_0_yyyyyzz_1, g_0_xxxxxx_0_yyyyzz_1, g_0_xxxxxx_0_yyyyzzz_0, g_0_xxxxxx_0_yyyyzzz_1, g_0_xxxxxx_0_yyyzzz_1, g_0_xxxxxx_0_yyyzzzz_0, g_0_xxxxxx_0_yyyzzzz_1, g_0_xxxxxx_0_yyzzzz_1, g_0_xxxxxx_0_yyzzzzz_0, g_0_xxxxxx_0_yyzzzzz_1, g_0_xxxxxx_0_yzzzzz_1, g_0_xxxxxx_0_yzzzzzz_0, g_0_xxxxxx_0_yzzzzzz_1, g_0_xxxxxx_0_zzzzzz_1, g_0_xxxxxx_0_zzzzzzz_0, g_0_xxxxxx_0_zzzzzzz_1, g_0_xxxxxxz_0_xxxxxxx_0, g_0_xxxxxxz_0_xxxxxxy_0, g_0_xxxxxxz_0_xxxxxxz_0, g_0_xxxxxxz_0_xxxxxyy_0, g_0_xxxxxxz_0_xxxxxyz_0, g_0_xxxxxxz_0_xxxxxzz_0, g_0_xxxxxxz_0_xxxxyyy_0, g_0_xxxxxxz_0_xxxxyyz_0, g_0_xxxxxxz_0_xxxxyzz_0, g_0_xxxxxxz_0_xxxxzzz_0, g_0_xxxxxxz_0_xxxyyyy_0, g_0_xxxxxxz_0_xxxyyyz_0, g_0_xxxxxxz_0_xxxyyzz_0, g_0_xxxxxxz_0_xxxyzzz_0, g_0_xxxxxxz_0_xxxzzzz_0, g_0_xxxxxxz_0_xxyyyyy_0, g_0_xxxxxxz_0_xxyyyyz_0, g_0_xxxxxxz_0_xxyyyzz_0, g_0_xxxxxxz_0_xxyyzzz_0, g_0_xxxxxxz_0_xxyzzzz_0, g_0_xxxxxxz_0_xxzzzzz_0, g_0_xxxxxxz_0_xyyyyyy_0, g_0_xxxxxxz_0_xyyyyyz_0, g_0_xxxxxxz_0_xyyyyzz_0, g_0_xxxxxxz_0_xyyyzzz_0, g_0_xxxxxxz_0_xyyzzzz_0, g_0_xxxxxxz_0_xyzzzzz_0, g_0_xxxxxxz_0_xzzzzzz_0, g_0_xxxxxxz_0_yyyyyyy_0, g_0_xxxxxxz_0_yyyyyyz_0, g_0_xxxxxxz_0_yyyyyzz_0, g_0_xxxxxxz_0_yyyyzzz_0, g_0_xxxxxxz_0_yyyzzzz_0, g_0_xxxxxxz_0_yyzzzzz_0, g_0_xxxxxxz_0_yzzzzzz_0, g_0_xxxxxxz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxz_0_xxxxxxx_0[i] = g_0_xxxxxx_0_xxxxxxx_0[i] * pb_z + g_0_xxxxxx_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxxxy_0[i] = g_0_xxxxxx_0_xxxxxxy_0[i] * pb_z + g_0_xxxxxx_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxxxz_0[i] = g_0_xxxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxxz_0[i] * pb_z + g_0_xxxxxx_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxxyy_0[i] = g_0_xxxxxx_0_xxxxxyy_0[i] * pb_z + g_0_xxxxxx_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxxyz_0[i] = g_0_xxxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxyz_0[i] * pb_z + g_0_xxxxxx_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxzz_0[i] * pb_z + g_0_xxxxxx_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxyyy_0[i] = g_0_xxxxxx_0_xxxxyyy_0[i] * pb_z + g_0_xxxxxx_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxyyz_0[i] = g_0_xxxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyyz_0[i] * pb_z + g_0_xxxxxx_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyzz_0[i] * pb_z + g_0_xxxxxx_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxzzz_0[i] * pb_z + g_0_xxxxxx_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyyyy_0[i] = g_0_xxxxxx_0_xxxyyyy_0[i] * pb_z + g_0_xxxxxx_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyyyz_0[i] = g_0_xxxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyyz_0[i] * pb_z + g_0_xxxxxx_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyzz_0[i] * pb_z + g_0_xxxxxx_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyzzz_0[i] * pb_z + g_0_xxxxxx_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxzzzz_0[i] * pb_z + g_0_xxxxxx_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyyyy_0[i] = g_0_xxxxxx_0_xxyyyyy_0[i] * pb_z + g_0_xxxxxx_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyyyz_0[i] = g_0_xxxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyyz_0[i] * pb_z + g_0_xxxxxx_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyzz_0[i] * pb_z + g_0_xxxxxx_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyzzz_0[i] * pb_z + g_0_xxxxxx_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzzzz_0[i] * pb_z + g_0_xxxxxx_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxzzzzz_0[i] * pb_z + g_0_xxxxxx_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyyyy_0[i] = g_0_xxxxxx_0_xyyyyyy_0[i] * pb_z + g_0_xxxxxx_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyyyz_0[i] = g_0_xxxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyyz_0[i] * pb_z + g_0_xxxxxx_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyzz_0[i] * pb_z + g_0_xxxxxx_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyzzz_0[i] * pb_z + g_0_xxxxxx_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzzzz_0[i] * pb_z + g_0_xxxxxx_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzzzz_0[i] * pb_z + g_0_xxxxxx_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xzzzzzz_0[i] * pb_z + g_0_xxxxxx_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyyyy_0[i] = g_0_xxxxxx_0_yyyyyyy_0[i] * pb_z + g_0_xxxxxx_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyyyz_0[i] = g_0_xxxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyyyz_0[i] * pb_z + g_0_xxxxxx_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyyzz_0[i] * pb_z + g_0_xxxxxx_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyzzz_0[i] * pb_z + g_0_xxxxxx_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyzzzz_0[i] * pb_z + g_0_xxxxxx_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyzzzzz_0[i] * pb_z + g_0_xxxxxx_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzzzzzz_0[i] * pb_z + g_0_xxxxxx_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_zzzzzzz_0[i] * pb_z + g_0_xxxxxx_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 108-144 components of targeted buffer : SKSK

    auto g_0_xxxxxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 108);

    auto g_0_xxxxxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 109);

    auto g_0_xxxxxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 110);

    auto g_0_xxxxxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 111);

    auto g_0_xxxxxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 112);

    auto g_0_xxxxxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 113);

    auto g_0_xxxxxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 114);

    auto g_0_xxxxxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 115);

    auto g_0_xxxxxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 116);

    auto g_0_xxxxxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 117);

    auto g_0_xxxxxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 118);

    auto g_0_xxxxxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 119);

    auto g_0_xxxxxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 120);

    auto g_0_xxxxxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 121);

    auto g_0_xxxxxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 122);

    auto g_0_xxxxxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 123);

    auto g_0_xxxxxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 124);

    auto g_0_xxxxxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 125);

    auto g_0_xxxxxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 126);

    auto g_0_xxxxxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 127);

    auto g_0_xxxxxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 128);

    auto g_0_xxxxxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 129);

    auto g_0_xxxxxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 130);

    auto g_0_xxxxxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 131);

    auto g_0_xxxxxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 132);

    auto g_0_xxxxxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 133);

    auto g_0_xxxxxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 134);

    auto g_0_xxxxxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 135);

    auto g_0_xxxxxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 136);

    auto g_0_xxxxxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 137);

    auto g_0_xxxxxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 138);

    auto g_0_xxxxxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 139);

    auto g_0_xxxxxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 140);

    auto g_0_xxxxxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 141);

    auto g_0_xxxxxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 142);

    auto g_0_xxxxxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 143);

    #pragma omp simd aligned(g_0_xxxxx_0_xxxxxxx_0, g_0_xxxxx_0_xxxxxxx_1, g_0_xxxxx_0_xxxxxxz_0, g_0_xxxxx_0_xxxxxxz_1, g_0_xxxxx_0_xxxxxzz_0, g_0_xxxxx_0_xxxxxzz_1, g_0_xxxxx_0_xxxxzzz_0, g_0_xxxxx_0_xxxxzzz_1, g_0_xxxxx_0_xxxzzzz_0, g_0_xxxxx_0_xxxzzzz_1, g_0_xxxxx_0_xxzzzzz_0, g_0_xxxxx_0_xxzzzzz_1, g_0_xxxxx_0_xzzzzzz_0, g_0_xxxxx_0_xzzzzzz_1, g_0_xxxxxy_0_xxxxxxx_0, g_0_xxxxxy_0_xxxxxxx_1, g_0_xxxxxy_0_xxxxxxz_0, g_0_xxxxxy_0_xxxxxxz_1, g_0_xxxxxy_0_xxxxxzz_0, g_0_xxxxxy_0_xxxxxzz_1, g_0_xxxxxy_0_xxxxzzz_0, g_0_xxxxxy_0_xxxxzzz_1, g_0_xxxxxy_0_xxxzzzz_0, g_0_xxxxxy_0_xxxzzzz_1, g_0_xxxxxy_0_xxzzzzz_0, g_0_xxxxxy_0_xxzzzzz_1, g_0_xxxxxy_0_xzzzzzz_0, g_0_xxxxxy_0_xzzzzzz_1, g_0_xxxxxyy_0_xxxxxxx_0, g_0_xxxxxyy_0_xxxxxxy_0, g_0_xxxxxyy_0_xxxxxxz_0, g_0_xxxxxyy_0_xxxxxyy_0, g_0_xxxxxyy_0_xxxxxyz_0, g_0_xxxxxyy_0_xxxxxzz_0, g_0_xxxxxyy_0_xxxxyyy_0, g_0_xxxxxyy_0_xxxxyyz_0, g_0_xxxxxyy_0_xxxxyzz_0, g_0_xxxxxyy_0_xxxxzzz_0, g_0_xxxxxyy_0_xxxyyyy_0, g_0_xxxxxyy_0_xxxyyyz_0, g_0_xxxxxyy_0_xxxyyzz_0, g_0_xxxxxyy_0_xxxyzzz_0, g_0_xxxxxyy_0_xxxzzzz_0, g_0_xxxxxyy_0_xxyyyyy_0, g_0_xxxxxyy_0_xxyyyyz_0, g_0_xxxxxyy_0_xxyyyzz_0, g_0_xxxxxyy_0_xxyyzzz_0, g_0_xxxxxyy_0_xxyzzzz_0, g_0_xxxxxyy_0_xxzzzzz_0, g_0_xxxxxyy_0_xyyyyyy_0, g_0_xxxxxyy_0_xyyyyyz_0, g_0_xxxxxyy_0_xyyyyzz_0, g_0_xxxxxyy_0_xyyyzzz_0, g_0_xxxxxyy_0_xyyzzzz_0, g_0_xxxxxyy_0_xyzzzzz_0, g_0_xxxxxyy_0_xzzzzzz_0, g_0_xxxxxyy_0_yyyyyyy_0, g_0_xxxxxyy_0_yyyyyyz_0, g_0_xxxxxyy_0_yyyyyzz_0, g_0_xxxxxyy_0_yyyyzzz_0, g_0_xxxxxyy_0_yyyzzzz_0, g_0_xxxxxyy_0_yyzzzzz_0, g_0_xxxxxyy_0_yzzzzzz_0, g_0_xxxxxyy_0_zzzzzzz_0, g_0_xxxxyy_0_xxxxxxy_0, g_0_xxxxyy_0_xxxxxxy_1, g_0_xxxxyy_0_xxxxxy_1, g_0_xxxxyy_0_xxxxxyy_0, g_0_xxxxyy_0_xxxxxyy_1, g_0_xxxxyy_0_xxxxxyz_0, g_0_xxxxyy_0_xxxxxyz_1, g_0_xxxxyy_0_xxxxyy_1, g_0_xxxxyy_0_xxxxyyy_0, g_0_xxxxyy_0_xxxxyyy_1, g_0_xxxxyy_0_xxxxyyz_0, g_0_xxxxyy_0_xxxxyyz_1, g_0_xxxxyy_0_xxxxyz_1, g_0_xxxxyy_0_xxxxyzz_0, g_0_xxxxyy_0_xxxxyzz_1, g_0_xxxxyy_0_xxxyyy_1, g_0_xxxxyy_0_xxxyyyy_0, g_0_xxxxyy_0_xxxyyyy_1, g_0_xxxxyy_0_xxxyyyz_0, g_0_xxxxyy_0_xxxyyyz_1, g_0_xxxxyy_0_xxxyyz_1, g_0_xxxxyy_0_xxxyyzz_0, g_0_xxxxyy_0_xxxyyzz_1, g_0_xxxxyy_0_xxxyzz_1, g_0_xxxxyy_0_xxxyzzz_0, g_0_xxxxyy_0_xxxyzzz_1, g_0_xxxxyy_0_xxyyyy_1, g_0_xxxxyy_0_xxyyyyy_0, g_0_xxxxyy_0_xxyyyyy_1, g_0_xxxxyy_0_xxyyyyz_0, g_0_xxxxyy_0_xxyyyyz_1, g_0_xxxxyy_0_xxyyyz_1, g_0_xxxxyy_0_xxyyyzz_0, g_0_xxxxyy_0_xxyyyzz_1, g_0_xxxxyy_0_xxyyzz_1, g_0_xxxxyy_0_xxyyzzz_0, g_0_xxxxyy_0_xxyyzzz_1, g_0_xxxxyy_0_xxyzzz_1, g_0_xxxxyy_0_xxyzzzz_0, g_0_xxxxyy_0_xxyzzzz_1, g_0_xxxxyy_0_xyyyyy_1, g_0_xxxxyy_0_xyyyyyy_0, g_0_xxxxyy_0_xyyyyyy_1, g_0_xxxxyy_0_xyyyyyz_0, g_0_xxxxyy_0_xyyyyyz_1, g_0_xxxxyy_0_xyyyyz_1, g_0_xxxxyy_0_xyyyyzz_0, g_0_xxxxyy_0_xyyyyzz_1, g_0_xxxxyy_0_xyyyzz_1, g_0_xxxxyy_0_xyyyzzz_0, g_0_xxxxyy_0_xyyyzzz_1, g_0_xxxxyy_0_xyyzzz_1, g_0_xxxxyy_0_xyyzzzz_0, g_0_xxxxyy_0_xyyzzzz_1, g_0_xxxxyy_0_xyzzzz_1, g_0_xxxxyy_0_xyzzzzz_0, g_0_xxxxyy_0_xyzzzzz_1, g_0_xxxxyy_0_yyyyyy_1, g_0_xxxxyy_0_yyyyyyy_0, g_0_xxxxyy_0_yyyyyyy_1, g_0_xxxxyy_0_yyyyyyz_0, g_0_xxxxyy_0_yyyyyyz_1, g_0_xxxxyy_0_yyyyyz_1, g_0_xxxxyy_0_yyyyyzz_0, g_0_xxxxyy_0_yyyyyzz_1, g_0_xxxxyy_0_yyyyzz_1, g_0_xxxxyy_0_yyyyzzz_0, g_0_xxxxyy_0_yyyyzzz_1, g_0_xxxxyy_0_yyyzzz_1, g_0_xxxxyy_0_yyyzzzz_0, g_0_xxxxyy_0_yyyzzzz_1, g_0_xxxxyy_0_yyzzzz_1, g_0_xxxxyy_0_yyzzzzz_0, g_0_xxxxyy_0_yyzzzzz_1, g_0_xxxxyy_0_yzzzzz_1, g_0_xxxxyy_0_yzzzzzz_0, g_0_xxxxyy_0_yzzzzzz_1, g_0_xxxxyy_0_zzzzzzz_0, g_0_xxxxyy_0_zzzzzzz_1, g_0_xxxyy_0_xxxxxxy_0, g_0_xxxyy_0_xxxxxxy_1, g_0_xxxyy_0_xxxxxyy_0, g_0_xxxyy_0_xxxxxyy_1, g_0_xxxyy_0_xxxxxyz_0, g_0_xxxyy_0_xxxxxyz_1, g_0_xxxyy_0_xxxxyyy_0, g_0_xxxyy_0_xxxxyyy_1, g_0_xxxyy_0_xxxxyyz_0, g_0_xxxyy_0_xxxxyyz_1, g_0_xxxyy_0_xxxxyzz_0, g_0_xxxyy_0_xxxxyzz_1, g_0_xxxyy_0_xxxyyyy_0, g_0_xxxyy_0_xxxyyyy_1, g_0_xxxyy_0_xxxyyyz_0, g_0_xxxyy_0_xxxyyyz_1, g_0_xxxyy_0_xxxyyzz_0, g_0_xxxyy_0_xxxyyzz_1, g_0_xxxyy_0_xxxyzzz_0, g_0_xxxyy_0_xxxyzzz_1, g_0_xxxyy_0_xxyyyyy_0, g_0_xxxyy_0_xxyyyyy_1, g_0_xxxyy_0_xxyyyyz_0, g_0_xxxyy_0_xxyyyyz_1, g_0_xxxyy_0_xxyyyzz_0, g_0_xxxyy_0_xxyyyzz_1, g_0_xxxyy_0_xxyyzzz_0, g_0_xxxyy_0_xxyyzzz_1, g_0_xxxyy_0_xxyzzzz_0, g_0_xxxyy_0_xxyzzzz_1, g_0_xxxyy_0_xyyyyyy_0, g_0_xxxyy_0_xyyyyyy_1, g_0_xxxyy_0_xyyyyyz_0, g_0_xxxyy_0_xyyyyyz_1, g_0_xxxyy_0_xyyyyzz_0, g_0_xxxyy_0_xyyyyzz_1, g_0_xxxyy_0_xyyyzzz_0, g_0_xxxyy_0_xyyyzzz_1, g_0_xxxyy_0_xyyzzzz_0, g_0_xxxyy_0_xyyzzzz_1, g_0_xxxyy_0_xyzzzzz_0, g_0_xxxyy_0_xyzzzzz_1, g_0_xxxyy_0_yyyyyyy_0, g_0_xxxyy_0_yyyyyyy_1, g_0_xxxyy_0_yyyyyyz_0, g_0_xxxyy_0_yyyyyyz_1, g_0_xxxyy_0_yyyyyzz_0, g_0_xxxyy_0_yyyyyzz_1, g_0_xxxyy_0_yyyyzzz_0, g_0_xxxyy_0_yyyyzzz_1, g_0_xxxyy_0_yyyzzzz_0, g_0_xxxyy_0_yyyzzzz_1, g_0_xxxyy_0_yyzzzzz_0, g_0_xxxyy_0_yyzzzzz_1, g_0_xxxyy_0_yzzzzzz_0, g_0_xxxyy_0_yzzzzzz_1, g_0_xxxyy_0_zzzzzzz_0, g_0_xxxyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyy_0_xxxxxxx_0[i] = g_0_xxxxx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxxxx_0[i] * pb_y + g_0_xxxxxy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxxxxy_0[i] = 4.0 * g_0_xxxyy_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxxyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxxxy_0[i] * pb_x + g_0_xxxxyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxxxz_0[i] = g_0_xxxxx_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxxxz_0[i] * pb_y + g_0_xxxxxy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxxxyy_0[i] = 4.0 * g_0_xxxyy_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxxyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxxyy_0[i] * pb_x + g_0_xxxxyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxxyz_0[i] = 4.0 * g_0_xxxyy_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxxyz_0[i] * pb_x + g_0_xxxxyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxxzz_0[i] = g_0_xxxxx_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxxzz_0[i] * pb_y + g_0_xxxxxy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxxyyy_0[i] = 4.0 * g_0_xxxyy_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxyyy_0[i] * pb_x + g_0_xxxxyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxyyz_0[i] = 4.0 * g_0_xxxyy_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxyyz_0[i] * pb_x + g_0_xxxxyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxyzz_0[i] = 4.0 * g_0_xxxyy_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxyzz_0[i] * pb_x + g_0_xxxxyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxzzz_0[i] = g_0_xxxxx_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxzzz_0[i] * pb_y + g_0_xxxxxy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxyyyy_0[i] = 4.0 * g_0_xxxyy_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyyyy_0[i] * pb_x + g_0_xxxxyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxyyyz_0[i] = 4.0 * g_0_xxxyy_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyyyz_0[i] * pb_x + g_0_xxxxyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxyyzz_0[i] = 4.0 * g_0_xxxyy_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyyzz_0[i] * pb_x + g_0_xxxxyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxyzzz_0[i] = 4.0 * g_0_xxxyy_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyzzz_0[i] * pb_x + g_0_xxxxyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxzzzz_0[i] = g_0_xxxxx_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxzzzz_0[i] * pb_y + g_0_xxxxxy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxyyyyy_0[i] = 4.0 * g_0_xxxyy_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyyyy_0[i] * pb_x + g_0_xxxxyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyyyyz_0[i] = 4.0 * g_0_xxxyy_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyyyz_0[i] * pb_x + g_0_xxxxyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyyyzz_0[i] = 4.0 * g_0_xxxyy_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyyzz_0[i] * pb_x + g_0_xxxxyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyyzzz_0[i] = 4.0 * g_0_xxxyy_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyzzz_0[i] * pb_x + g_0_xxxxyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyzzzz_0[i] = 4.0 * g_0_xxxyy_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyzzzz_0[i] * pb_x + g_0_xxxxyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxzzzzz_0[i] = g_0_xxxxx_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxzzzzz_0[i] * pb_y + g_0_xxxxxy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xyyyyyy_0[i] = 4.0 * g_0_xxxyy_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyyyy_0[i] * pb_x + g_0_xxxxyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyyyyz_0[i] = 4.0 * g_0_xxxyy_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyyyz_0[i] * pb_x + g_0_xxxxyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyyyzz_0[i] = 4.0 * g_0_xxxyy_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyyzz_0[i] * pb_x + g_0_xxxxyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyyzzz_0[i] = 4.0 * g_0_xxxyy_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyzzz_0[i] * pb_x + g_0_xxxxyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyzzzz_0[i] = 4.0 * g_0_xxxyy_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyzzzz_0[i] * pb_x + g_0_xxxxyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyzzzzz_0[i] = 4.0 * g_0_xxxyy_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyzzzzz_0[i] * pb_x + g_0_xxxxyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xzzzzzz_0[i] = g_0_xxxxx_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xzzzzzz_0[i] * pb_y + g_0_xxxxxy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_yyyyyyy_0[i] = 4.0 * g_0_xxxyy_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyyyy_0[i] * pb_x + g_0_xxxxyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyyyyz_0[i] = 4.0 * g_0_xxxyy_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyyyz_0[i] * pb_x + g_0_xxxxyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyyyzz_0[i] = 4.0 * g_0_xxxyy_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyyzz_0[i] * pb_x + g_0_xxxxyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyyzzz_0[i] = 4.0 * g_0_xxxyy_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyyzzz_0[i] * pb_x + g_0_xxxxyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyzzzz_0[i] = 4.0 * g_0_xxxyy_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyyzzzz_0[i] * pb_x + g_0_xxxxyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyzzzzz_0[i] = 4.0 * g_0_xxxyy_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yyzzzzz_0[i] * pb_x + g_0_xxxxyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yzzzzzz_0[i] = 4.0 * g_0_xxxyy_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_yzzzzzz_0[i] * pb_x + g_0_xxxxyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_zzzzzzz_0[i] = 4.0 * g_0_xxxyy_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_zzzzzzz_0[i] * pb_x + g_0_xxxxyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 144-180 components of targeted buffer : SKSK

    auto g_0_xxxxxyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 144);

    auto g_0_xxxxxyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 145);

    auto g_0_xxxxxyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 146);

    auto g_0_xxxxxyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 147);

    auto g_0_xxxxxyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 148);

    auto g_0_xxxxxyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 149);

    auto g_0_xxxxxyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 150);

    auto g_0_xxxxxyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 151);

    auto g_0_xxxxxyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 152);

    auto g_0_xxxxxyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 153);

    auto g_0_xxxxxyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 154);

    auto g_0_xxxxxyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 155);

    auto g_0_xxxxxyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 156);

    auto g_0_xxxxxyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 157);

    auto g_0_xxxxxyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 158);

    auto g_0_xxxxxyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 159);

    auto g_0_xxxxxyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 160);

    auto g_0_xxxxxyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 161);

    auto g_0_xxxxxyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 162);

    auto g_0_xxxxxyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 163);

    auto g_0_xxxxxyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 164);

    auto g_0_xxxxxyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 165);

    auto g_0_xxxxxyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 166);

    auto g_0_xxxxxyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 167);

    auto g_0_xxxxxyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 168);

    auto g_0_xxxxxyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 169);

    auto g_0_xxxxxyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 170);

    auto g_0_xxxxxyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 171);

    auto g_0_xxxxxyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 172);

    auto g_0_xxxxxyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 173);

    auto g_0_xxxxxyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 174);

    auto g_0_xxxxxyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 175);

    auto g_0_xxxxxyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 176);

    auto g_0_xxxxxyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 177);

    auto g_0_xxxxxyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 178);

    auto g_0_xxxxxyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 179);

    #pragma omp simd aligned(g_0_xxxxxy_0_xxxxxxy_0, g_0_xxxxxy_0_xxxxxxy_1, g_0_xxxxxy_0_xxxxxyy_0, g_0_xxxxxy_0_xxxxxyy_1, g_0_xxxxxy_0_xxxxyyy_0, g_0_xxxxxy_0_xxxxyyy_1, g_0_xxxxxy_0_xxxyyyy_0, g_0_xxxxxy_0_xxxyyyy_1, g_0_xxxxxy_0_xxyyyyy_0, g_0_xxxxxy_0_xxyyyyy_1, g_0_xxxxxy_0_xyyyyyy_0, g_0_xxxxxy_0_xyyyyyy_1, g_0_xxxxxy_0_yyyyyyy_0, g_0_xxxxxy_0_yyyyyyy_1, g_0_xxxxxyz_0_xxxxxxx_0, g_0_xxxxxyz_0_xxxxxxy_0, g_0_xxxxxyz_0_xxxxxxz_0, g_0_xxxxxyz_0_xxxxxyy_0, g_0_xxxxxyz_0_xxxxxyz_0, g_0_xxxxxyz_0_xxxxxzz_0, g_0_xxxxxyz_0_xxxxyyy_0, g_0_xxxxxyz_0_xxxxyyz_0, g_0_xxxxxyz_0_xxxxyzz_0, g_0_xxxxxyz_0_xxxxzzz_0, g_0_xxxxxyz_0_xxxyyyy_0, g_0_xxxxxyz_0_xxxyyyz_0, g_0_xxxxxyz_0_xxxyyzz_0, g_0_xxxxxyz_0_xxxyzzz_0, g_0_xxxxxyz_0_xxxzzzz_0, g_0_xxxxxyz_0_xxyyyyy_0, g_0_xxxxxyz_0_xxyyyyz_0, g_0_xxxxxyz_0_xxyyyzz_0, g_0_xxxxxyz_0_xxyyzzz_0, g_0_xxxxxyz_0_xxyzzzz_0, g_0_xxxxxyz_0_xxzzzzz_0, g_0_xxxxxyz_0_xyyyyyy_0, g_0_xxxxxyz_0_xyyyyyz_0, g_0_xxxxxyz_0_xyyyyzz_0, g_0_xxxxxyz_0_xyyyzzz_0, g_0_xxxxxyz_0_xyyzzzz_0, g_0_xxxxxyz_0_xyzzzzz_0, g_0_xxxxxyz_0_xzzzzzz_0, g_0_xxxxxyz_0_yyyyyyy_0, g_0_xxxxxyz_0_yyyyyyz_0, g_0_xxxxxyz_0_yyyyyzz_0, g_0_xxxxxyz_0_yyyyzzz_0, g_0_xxxxxyz_0_yyyzzzz_0, g_0_xxxxxyz_0_yyzzzzz_0, g_0_xxxxxyz_0_yzzzzzz_0, g_0_xxxxxyz_0_zzzzzzz_0, g_0_xxxxxz_0_xxxxxxx_0, g_0_xxxxxz_0_xxxxxxx_1, g_0_xxxxxz_0_xxxxxxz_0, g_0_xxxxxz_0_xxxxxxz_1, g_0_xxxxxz_0_xxxxxyz_0, g_0_xxxxxz_0_xxxxxyz_1, g_0_xxxxxz_0_xxxxxz_1, g_0_xxxxxz_0_xxxxxzz_0, g_0_xxxxxz_0_xxxxxzz_1, g_0_xxxxxz_0_xxxxyyz_0, g_0_xxxxxz_0_xxxxyyz_1, g_0_xxxxxz_0_xxxxyz_1, g_0_xxxxxz_0_xxxxyzz_0, g_0_xxxxxz_0_xxxxyzz_1, g_0_xxxxxz_0_xxxxzz_1, g_0_xxxxxz_0_xxxxzzz_0, g_0_xxxxxz_0_xxxxzzz_1, g_0_xxxxxz_0_xxxyyyz_0, g_0_xxxxxz_0_xxxyyyz_1, g_0_xxxxxz_0_xxxyyz_1, g_0_xxxxxz_0_xxxyyzz_0, g_0_xxxxxz_0_xxxyyzz_1, g_0_xxxxxz_0_xxxyzz_1, g_0_xxxxxz_0_xxxyzzz_0, g_0_xxxxxz_0_xxxyzzz_1, g_0_xxxxxz_0_xxxzzz_1, g_0_xxxxxz_0_xxxzzzz_0, g_0_xxxxxz_0_xxxzzzz_1, g_0_xxxxxz_0_xxyyyyz_0, g_0_xxxxxz_0_xxyyyyz_1, g_0_xxxxxz_0_xxyyyz_1, g_0_xxxxxz_0_xxyyyzz_0, g_0_xxxxxz_0_xxyyyzz_1, g_0_xxxxxz_0_xxyyzz_1, g_0_xxxxxz_0_xxyyzzz_0, g_0_xxxxxz_0_xxyyzzz_1, g_0_xxxxxz_0_xxyzzz_1, g_0_xxxxxz_0_xxyzzzz_0, g_0_xxxxxz_0_xxyzzzz_1, g_0_xxxxxz_0_xxzzzz_1, g_0_xxxxxz_0_xxzzzzz_0, g_0_xxxxxz_0_xxzzzzz_1, g_0_xxxxxz_0_xyyyyyz_0, g_0_xxxxxz_0_xyyyyyz_1, g_0_xxxxxz_0_xyyyyz_1, g_0_xxxxxz_0_xyyyyzz_0, g_0_xxxxxz_0_xyyyyzz_1, g_0_xxxxxz_0_xyyyzz_1, g_0_xxxxxz_0_xyyyzzz_0, g_0_xxxxxz_0_xyyyzzz_1, g_0_xxxxxz_0_xyyzzz_1, g_0_xxxxxz_0_xyyzzzz_0, g_0_xxxxxz_0_xyyzzzz_1, g_0_xxxxxz_0_xyzzzz_1, g_0_xxxxxz_0_xyzzzzz_0, g_0_xxxxxz_0_xyzzzzz_1, g_0_xxxxxz_0_xzzzzz_1, g_0_xxxxxz_0_xzzzzzz_0, g_0_xxxxxz_0_xzzzzzz_1, g_0_xxxxxz_0_yyyyyyz_0, g_0_xxxxxz_0_yyyyyyz_1, g_0_xxxxxz_0_yyyyyz_1, g_0_xxxxxz_0_yyyyyzz_0, g_0_xxxxxz_0_yyyyyzz_1, g_0_xxxxxz_0_yyyyzz_1, g_0_xxxxxz_0_yyyyzzz_0, g_0_xxxxxz_0_yyyyzzz_1, g_0_xxxxxz_0_yyyzzz_1, g_0_xxxxxz_0_yyyzzzz_0, g_0_xxxxxz_0_yyyzzzz_1, g_0_xxxxxz_0_yyzzzz_1, g_0_xxxxxz_0_yyzzzzz_0, g_0_xxxxxz_0_yyzzzzz_1, g_0_xxxxxz_0_yzzzzz_1, g_0_xxxxxz_0_yzzzzzz_0, g_0_xxxxxz_0_yzzzzzz_1, g_0_xxxxxz_0_zzzzzz_1, g_0_xxxxxz_0_zzzzzzz_0, g_0_xxxxxz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyz_0_xxxxxxx_0[i] = g_0_xxxxxz_0_xxxxxxx_0[i] * pb_y + g_0_xxxxxz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxxxy_0[i] = g_0_xxxxxy_0_xxxxxxy_0[i] * pb_z + g_0_xxxxxy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxxxxz_0[i] = g_0_xxxxxz_0_xxxxxxz_0[i] * pb_y + g_0_xxxxxz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxxyy_0[i] = g_0_xxxxxy_0_xxxxxyy_0[i] * pb_z + g_0_xxxxxy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxxxyz_0[i] = g_0_xxxxxz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxxxyz_0[i] * pb_y + g_0_xxxxxz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxxzz_0[i] = g_0_xxxxxz_0_xxxxxzz_0[i] * pb_y + g_0_xxxxxz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxyyy_0[i] = g_0_xxxxxy_0_xxxxyyy_0[i] * pb_z + g_0_xxxxxy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxxz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxxyyz_0[i] * pb_y + g_0_xxxxxz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxyzz_0[i] = g_0_xxxxxz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxxyzz_0[i] * pb_y + g_0_xxxxxz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxzzz_0[i] = g_0_xxxxxz_0_xxxxzzz_0[i] * pb_y + g_0_xxxxxz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxyyyy_0[i] = g_0_xxxxxy_0_xxxyyyy_0[i] * pb_z + g_0_xxxxxy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxxz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxyyyz_0[i] * pb_y + g_0_xxxxxz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxxz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxyyzz_0[i] * pb_y + g_0_xxxxxz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxyzzz_0[i] = g_0_xxxxxz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxyzzz_0[i] * pb_y + g_0_xxxxxz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxzzzz_0[i] = g_0_xxxxxz_0_xxxzzzz_0[i] * pb_y + g_0_xxxxxz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyyyyy_0[i] = g_0_xxxxxy_0_xxyyyyy_0[i] * pb_z + g_0_xxxxxy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxxz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyyyyz_0[i] * pb_y + g_0_xxxxxz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxxz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyyyzz_0[i] * pb_y + g_0_xxxxxz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxxz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyyzzz_0[i] * pb_y + g_0_xxxxxz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyzzzz_0[i] = g_0_xxxxxz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyzzzz_0[i] * pb_y + g_0_xxxxxz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxzzzzz_0[i] = g_0_xxxxxz_0_xxzzzzz_0[i] * pb_y + g_0_xxxxxz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyyyyy_0[i] = g_0_xxxxxy_0_xyyyyyy_0[i] * pb_z + g_0_xxxxxy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxxz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyyyyz_0[i] * pb_y + g_0_xxxxxz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxxz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyyyzz_0[i] * pb_y + g_0_xxxxxz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxxz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyyzzz_0[i] * pb_y + g_0_xxxxxz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxxz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyzzzz_0[i] * pb_y + g_0_xxxxxz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyzzzzz_0[i] = g_0_xxxxxz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyzzzzz_0[i] * pb_y + g_0_xxxxxz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xzzzzzz_0[i] = g_0_xxxxxz_0_xzzzzzz_0[i] * pb_y + g_0_xxxxxz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyyyyy_0[i] = g_0_xxxxxy_0_yyyyyyy_0[i] * pb_z + g_0_xxxxxy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxxz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyyyyz_0[i] * pb_y + g_0_xxxxxz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxxz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyyyzz_0[i] * pb_y + g_0_xxxxxz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxxz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyyzzz_0[i] * pb_y + g_0_xxxxxz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxxz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyzzzz_0[i] * pb_y + g_0_xxxxxz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxxz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyzzzzz_0[i] * pb_y + g_0_xxxxxz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yzzzzzz_0[i] = g_0_xxxxxz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yzzzzzz_0[i] * pb_y + g_0_xxxxxz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_zzzzzzz_0[i] = g_0_xxxxxz_0_zzzzzzz_0[i] * pb_y + g_0_xxxxxz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 180-216 components of targeted buffer : SKSK

    auto g_0_xxxxxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 180);

    auto g_0_xxxxxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 181);

    auto g_0_xxxxxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 182);

    auto g_0_xxxxxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 183);

    auto g_0_xxxxxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 184);

    auto g_0_xxxxxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 185);

    auto g_0_xxxxxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 186);

    auto g_0_xxxxxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 187);

    auto g_0_xxxxxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 188);

    auto g_0_xxxxxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 189);

    auto g_0_xxxxxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 190);

    auto g_0_xxxxxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 191);

    auto g_0_xxxxxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 192);

    auto g_0_xxxxxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 193);

    auto g_0_xxxxxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 194);

    auto g_0_xxxxxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 195);

    auto g_0_xxxxxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 196);

    auto g_0_xxxxxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 197);

    auto g_0_xxxxxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 198);

    auto g_0_xxxxxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 199);

    auto g_0_xxxxxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 200);

    auto g_0_xxxxxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 201);

    auto g_0_xxxxxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 202);

    auto g_0_xxxxxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 203);

    auto g_0_xxxxxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 204);

    auto g_0_xxxxxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 205);

    auto g_0_xxxxxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 206);

    auto g_0_xxxxxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 207);

    auto g_0_xxxxxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 208);

    auto g_0_xxxxxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 209);

    auto g_0_xxxxxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 210);

    auto g_0_xxxxxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 211);

    auto g_0_xxxxxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 212);

    auto g_0_xxxxxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 213);

    auto g_0_xxxxxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 214);

    auto g_0_xxxxxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 215);

    #pragma omp simd aligned(g_0_xxxxx_0_xxxxxxx_0, g_0_xxxxx_0_xxxxxxx_1, g_0_xxxxx_0_xxxxxxy_0, g_0_xxxxx_0_xxxxxxy_1, g_0_xxxxx_0_xxxxxyy_0, g_0_xxxxx_0_xxxxxyy_1, g_0_xxxxx_0_xxxxyyy_0, g_0_xxxxx_0_xxxxyyy_1, g_0_xxxxx_0_xxxyyyy_0, g_0_xxxxx_0_xxxyyyy_1, g_0_xxxxx_0_xxyyyyy_0, g_0_xxxxx_0_xxyyyyy_1, g_0_xxxxx_0_xyyyyyy_0, g_0_xxxxx_0_xyyyyyy_1, g_0_xxxxxz_0_xxxxxxx_0, g_0_xxxxxz_0_xxxxxxx_1, g_0_xxxxxz_0_xxxxxxy_0, g_0_xxxxxz_0_xxxxxxy_1, g_0_xxxxxz_0_xxxxxyy_0, g_0_xxxxxz_0_xxxxxyy_1, g_0_xxxxxz_0_xxxxyyy_0, g_0_xxxxxz_0_xxxxyyy_1, g_0_xxxxxz_0_xxxyyyy_0, g_0_xxxxxz_0_xxxyyyy_1, g_0_xxxxxz_0_xxyyyyy_0, g_0_xxxxxz_0_xxyyyyy_1, g_0_xxxxxz_0_xyyyyyy_0, g_0_xxxxxz_0_xyyyyyy_1, g_0_xxxxxzz_0_xxxxxxx_0, g_0_xxxxxzz_0_xxxxxxy_0, g_0_xxxxxzz_0_xxxxxxz_0, g_0_xxxxxzz_0_xxxxxyy_0, g_0_xxxxxzz_0_xxxxxyz_0, g_0_xxxxxzz_0_xxxxxzz_0, g_0_xxxxxzz_0_xxxxyyy_0, g_0_xxxxxzz_0_xxxxyyz_0, g_0_xxxxxzz_0_xxxxyzz_0, g_0_xxxxxzz_0_xxxxzzz_0, g_0_xxxxxzz_0_xxxyyyy_0, g_0_xxxxxzz_0_xxxyyyz_0, g_0_xxxxxzz_0_xxxyyzz_0, g_0_xxxxxzz_0_xxxyzzz_0, g_0_xxxxxzz_0_xxxzzzz_0, g_0_xxxxxzz_0_xxyyyyy_0, g_0_xxxxxzz_0_xxyyyyz_0, g_0_xxxxxzz_0_xxyyyzz_0, g_0_xxxxxzz_0_xxyyzzz_0, g_0_xxxxxzz_0_xxyzzzz_0, g_0_xxxxxzz_0_xxzzzzz_0, g_0_xxxxxzz_0_xyyyyyy_0, g_0_xxxxxzz_0_xyyyyyz_0, g_0_xxxxxzz_0_xyyyyzz_0, g_0_xxxxxzz_0_xyyyzzz_0, g_0_xxxxxzz_0_xyyzzzz_0, g_0_xxxxxzz_0_xyzzzzz_0, g_0_xxxxxzz_0_xzzzzzz_0, g_0_xxxxxzz_0_yyyyyyy_0, g_0_xxxxxzz_0_yyyyyyz_0, g_0_xxxxxzz_0_yyyyyzz_0, g_0_xxxxxzz_0_yyyyzzz_0, g_0_xxxxxzz_0_yyyzzzz_0, g_0_xxxxxzz_0_yyzzzzz_0, g_0_xxxxxzz_0_yzzzzzz_0, g_0_xxxxxzz_0_zzzzzzz_0, g_0_xxxxzz_0_xxxxxxz_0, g_0_xxxxzz_0_xxxxxxz_1, g_0_xxxxzz_0_xxxxxyz_0, g_0_xxxxzz_0_xxxxxyz_1, g_0_xxxxzz_0_xxxxxz_1, g_0_xxxxzz_0_xxxxxzz_0, g_0_xxxxzz_0_xxxxxzz_1, g_0_xxxxzz_0_xxxxyyz_0, g_0_xxxxzz_0_xxxxyyz_1, g_0_xxxxzz_0_xxxxyz_1, g_0_xxxxzz_0_xxxxyzz_0, g_0_xxxxzz_0_xxxxyzz_1, g_0_xxxxzz_0_xxxxzz_1, g_0_xxxxzz_0_xxxxzzz_0, g_0_xxxxzz_0_xxxxzzz_1, g_0_xxxxzz_0_xxxyyyz_0, g_0_xxxxzz_0_xxxyyyz_1, g_0_xxxxzz_0_xxxyyz_1, g_0_xxxxzz_0_xxxyyzz_0, g_0_xxxxzz_0_xxxyyzz_1, g_0_xxxxzz_0_xxxyzz_1, g_0_xxxxzz_0_xxxyzzz_0, g_0_xxxxzz_0_xxxyzzz_1, g_0_xxxxzz_0_xxxzzz_1, g_0_xxxxzz_0_xxxzzzz_0, g_0_xxxxzz_0_xxxzzzz_1, g_0_xxxxzz_0_xxyyyyz_0, g_0_xxxxzz_0_xxyyyyz_1, g_0_xxxxzz_0_xxyyyz_1, g_0_xxxxzz_0_xxyyyzz_0, g_0_xxxxzz_0_xxyyyzz_1, g_0_xxxxzz_0_xxyyzz_1, g_0_xxxxzz_0_xxyyzzz_0, g_0_xxxxzz_0_xxyyzzz_1, g_0_xxxxzz_0_xxyzzz_1, g_0_xxxxzz_0_xxyzzzz_0, g_0_xxxxzz_0_xxyzzzz_1, g_0_xxxxzz_0_xxzzzz_1, g_0_xxxxzz_0_xxzzzzz_0, g_0_xxxxzz_0_xxzzzzz_1, g_0_xxxxzz_0_xyyyyyz_0, g_0_xxxxzz_0_xyyyyyz_1, g_0_xxxxzz_0_xyyyyz_1, g_0_xxxxzz_0_xyyyyzz_0, g_0_xxxxzz_0_xyyyyzz_1, g_0_xxxxzz_0_xyyyzz_1, g_0_xxxxzz_0_xyyyzzz_0, g_0_xxxxzz_0_xyyyzzz_1, g_0_xxxxzz_0_xyyzzz_1, g_0_xxxxzz_0_xyyzzzz_0, g_0_xxxxzz_0_xyyzzzz_1, g_0_xxxxzz_0_xyzzzz_1, g_0_xxxxzz_0_xyzzzzz_0, g_0_xxxxzz_0_xyzzzzz_1, g_0_xxxxzz_0_xzzzzz_1, g_0_xxxxzz_0_xzzzzzz_0, g_0_xxxxzz_0_xzzzzzz_1, g_0_xxxxzz_0_yyyyyyy_0, g_0_xxxxzz_0_yyyyyyy_1, g_0_xxxxzz_0_yyyyyyz_0, g_0_xxxxzz_0_yyyyyyz_1, g_0_xxxxzz_0_yyyyyz_1, g_0_xxxxzz_0_yyyyyzz_0, g_0_xxxxzz_0_yyyyyzz_1, g_0_xxxxzz_0_yyyyzz_1, g_0_xxxxzz_0_yyyyzzz_0, g_0_xxxxzz_0_yyyyzzz_1, g_0_xxxxzz_0_yyyzzz_1, g_0_xxxxzz_0_yyyzzzz_0, g_0_xxxxzz_0_yyyzzzz_1, g_0_xxxxzz_0_yyzzzz_1, g_0_xxxxzz_0_yyzzzzz_0, g_0_xxxxzz_0_yyzzzzz_1, g_0_xxxxzz_0_yzzzzz_1, g_0_xxxxzz_0_yzzzzzz_0, g_0_xxxxzz_0_yzzzzzz_1, g_0_xxxxzz_0_zzzzzz_1, g_0_xxxxzz_0_zzzzzzz_0, g_0_xxxxzz_0_zzzzzzz_1, g_0_xxxzz_0_xxxxxxz_0, g_0_xxxzz_0_xxxxxxz_1, g_0_xxxzz_0_xxxxxyz_0, g_0_xxxzz_0_xxxxxyz_1, g_0_xxxzz_0_xxxxxzz_0, g_0_xxxzz_0_xxxxxzz_1, g_0_xxxzz_0_xxxxyyz_0, g_0_xxxzz_0_xxxxyyz_1, g_0_xxxzz_0_xxxxyzz_0, g_0_xxxzz_0_xxxxyzz_1, g_0_xxxzz_0_xxxxzzz_0, g_0_xxxzz_0_xxxxzzz_1, g_0_xxxzz_0_xxxyyyz_0, g_0_xxxzz_0_xxxyyyz_1, g_0_xxxzz_0_xxxyyzz_0, g_0_xxxzz_0_xxxyyzz_1, g_0_xxxzz_0_xxxyzzz_0, g_0_xxxzz_0_xxxyzzz_1, g_0_xxxzz_0_xxxzzzz_0, g_0_xxxzz_0_xxxzzzz_1, g_0_xxxzz_0_xxyyyyz_0, g_0_xxxzz_0_xxyyyyz_1, g_0_xxxzz_0_xxyyyzz_0, g_0_xxxzz_0_xxyyyzz_1, g_0_xxxzz_0_xxyyzzz_0, g_0_xxxzz_0_xxyyzzz_1, g_0_xxxzz_0_xxyzzzz_0, g_0_xxxzz_0_xxyzzzz_1, g_0_xxxzz_0_xxzzzzz_0, g_0_xxxzz_0_xxzzzzz_1, g_0_xxxzz_0_xyyyyyz_0, g_0_xxxzz_0_xyyyyyz_1, g_0_xxxzz_0_xyyyyzz_0, g_0_xxxzz_0_xyyyyzz_1, g_0_xxxzz_0_xyyyzzz_0, g_0_xxxzz_0_xyyyzzz_1, g_0_xxxzz_0_xyyzzzz_0, g_0_xxxzz_0_xyyzzzz_1, g_0_xxxzz_0_xyzzzzz_0, g_0_xxxzz_0_xyzzzzz_1, g_0_xxxzz_0_xzzzzzz_0, g_0_xxxzz_0_xzzzzzz_1, g_0_xxxzz_0_yyyyyyy_0, g_0_xxxzz_0_yyyyyyy_1, g_0_xxxzz_0_yyyyyyz_0, g_0_xxxzz_0_yyyyyyz_1, g_0_xxxzz_0_yyyyyzz_0, g_0_xxxzz_0_yyyyyzz_1, g_0_xxxzz_0_yyyyzzz_0, g_0_xxxzz_0_yyyyzzz_1, g_0_xxxzz_0_yyyzzzz_0, g_0_xxxzz_0_yyyzzzz_1, g_0_xxxzz_0_yyzzzzz_0, g_0_xxxzz_0_yyzzzzz_1, g_0_xxxzz_0_yzzzzzz_0, g_0_xxxzz_0_yzzzzzz_1, g_0_xxxzz_0_zzzzzzz_0, g_0_xxxzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzz_0_xxxxxxx_0[i] = g_0_xxxxx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxxxx_0[i] * pb_z + g_0_xxxxxz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxxxy_0[i] = g_0_xxxxx_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxxxy_0[i] * pb_z + g_0_xxxxxz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxxxz_0[i] = 4.0 * g_0_xxxzz_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxxzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxxxz_0[i] * pb_x + g_0_xxxxzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxxxyy_0[i] = g_0_xxxxx_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxxyy_0[i] * pb_z + g_0_xxxxxz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxxyz_0[i] = 4.0 * g_0_xxxzz_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxxyz_0[i] * pb_x + g_0_xxxxzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxxxzz_0[i] = 4.0 * g_0_xxxzz_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxxzz_0[i] * pb_x + g_0_xxxxzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxxyyy_0[i] = g_0_xxxxx_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxyyy_0[i] * pb_z + g_0_xxxxxz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxyyz_0[i] = 4.0 * g_0_xxxzz_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxyyz_0[i] * pb_x + g_0_xxxxzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxxyzz_0[i] = 4.0 * g_0_xxxzz_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxyzz_0[i] * pb_x + g_0_xxxxzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxxzzz_0[i] = 4.0 * g_0_xxxzz_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxzzz_0[i] * pb_x + g_0_xxxxzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxyyyy_0[i] = g_0_xxxxx_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxyyyy_0[i] * pb_z + g_0_xxxxxz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxyyyz_0[i] = 4.0 * g_0_xxxzz_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyyyz_0[i] * pb_x + g_0_xxxxzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxyyzz_0[i] = 4.0 * g_0_xxxzz_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyyzz_0[i] * pb_x + g_0_xxxxzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxyzzz_0[i] = 4.0 * g_0_xxxzz_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyzzz_0[i] * pb_x + g_0_xxxxzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxzz_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxzzzz_0[i] * pb_x + g_0_xxxxzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyyyyy_0[i] = g_0_xxxxx_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxyyyyy_0[i] * pb_z + g_0_xxxxxz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxzz_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyyyz_0[i] * pb_x + g_0_xxxxzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyyyzz_0[i] = 4.0 * g_0_xxxzz_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyyzz_0[i] * pb_x + g_0_xxxxzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyyzzz_0[i] = 4.0 * g_0_xxxzz_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyzzz_0[i] * pb_x + g_0_xxxxzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxzz_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyzzzz_0[i] * pb_x + g_0_xxxxzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxzzzzz_0[i] = 4.0 * g_0_xxxzz_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxzzzzz_0[i] * pb_x + g_0_xxxxzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyyyyy_0[i] = g_0_xxxxx_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xyyyyyy_0[i] * pb_z + g_0_xxxxxz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xyyyyyz_0[i] = 4.0 * g_0_xxxzz_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyyyz_0[i] * pb_x + g_0_xxxxzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxzz_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyyzz_0[i] * pb_x + g_0_xxxxzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyyzzz_0[i] = 4.0 * g_0_xxxzz_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyzzz_0[i] * pb_x + g_0_xxxxzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxzz_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyzzzz_0[i] * pb_x + g_0_xxxxzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyzzzzz_0[i] = 4.0 * g_0_xxxzz_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyzzzzz_0[i] * pb_x + g_0_xxxxzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xzzzzzz_0[i] = 4.0 * g_0_xxxzz_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xzzzzzz_0[i] * pb_x + g_0_xxxxzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyyyy_0[i] = 4.0 * g_0_xxxzz_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyyyyy_0[i] * pb_x + g_0_xxxxzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyyyz_0[i] = 4.0 * g_0_xxxzz_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyyyyz_0[i] * pb_x + g_0_xxxxzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyyzz_0[i] = 4.0 * g_0_xxxzz_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyyyzz_0[i] * pb_x + g_0_xxxxzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxzz_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyyzzz_0[i] * pb_x + g_0_xxxxzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxzz_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyyzzzz_0[i] * pb_x + g_0_xxxxzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyzzzzz_0[i] = 4.0 * g_0_xxxzz_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yyzzzzz_0[i] * pb_x + g_0_xxxxzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yzzzzzz_0[i] = 4.0 * g_0_xxxzz_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_yzzzzzz_0[i] * pb_x + g_0_xxxxzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_zzzzzzz_0[i] = 4.0 * g_0_xxxzz_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxzz_0_zzzzzzz_0[i] * pb_x + g_0_xxxxzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 216-252 components of targeted buffer : SKSK

    auto g_0_xxxxyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 216);

    auto g_0_xxxxyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 217);

    auto g_0_xxxxyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 218);

    auto g_0_xxxxyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 219);

    auto g_0_xxxxyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 220);

    auto g_0_xxxxyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 221);

    auto g_0_xxxxyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 222);

    auto g_0_xxxxyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 223);

    auto g_0_xxxxyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 224);

    auto g_0_xxxxyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 225);

    auto g_0_xxxxyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 226);

    auto g_0_xxxxyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 227);

    auto g_0_xxxxyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 228);

    auto g_0_xxxxyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 229);

    auto g_0_xxxxyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 230);

    auto g_0_xxxxyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 231);

    auto g_0_xxxxyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 232);

    auto g_0_xxxxyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 233);

    auto g_0_xxxxyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 234);

    auto g_0_xxxxyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 235);

    auto g_0_xxxxyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 236);

    auto g_0_xxxxyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 237);

    auto g_0_xxxxyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 238);

    auto g_0_xxxxyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 239);

    auto g_0_xxxxyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 240);

    auto g_0_xxxxyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 241);

    auto g_0_xxxxyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 242);

    auto g_0_xxxxyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 243);

    auto g_0_xxxxyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 244);

    auto g_0_xxxxyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 245);

    auto g_0_xxxxyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 246);

    auto g_0_xxxxyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 247);

    auto g_0_xxxxyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 248);

    auto g_0_xxxxyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 249);

    auto g_0_xxxxyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 250);

    auto g_0_xxxxyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 251);

    #pragma omp simd aligned(g_0_xxxxy_0_xxxxxxx_0, g_0_xxxxy_0_xxxxxxx_1, g_0_xxxxy_0_xxxxxxz_0, g_0_xxxxy_0_xxxxxxz_1, g_0_xxxxy_0_xxxxxzz_0, g_0_xxxxy_0_xxxxxzz_1, g_0_xxxxy_0_xxxxzzz_0, g_0_xxxxy_0_xxxxzzz_1, g_0_xxxxy_0_xxxzzzz_0, g_0_xxxxy_0_xxxzzzz_1, g_0_xxxxy_0_xxzzzzz_0, g_0_xxxxy_0_xxzzzzz_1, g_0_xxxxy_0_xzzzzzz_0, g_0_xxxxy_0_xzzzzzz_1, g_0_xxxxyy_0_xxxxxxx_0, g_0_xxxxyy_0_xxxxxxx_1, g_0_xxxxyy_0_xxxxxxz_0, g_0_xxxxyy_0_xxxxxxz_1, g_0_xxxxyy_0_xxxxxzz_0, g_0_xxxxyy_0_xxxxxzz_1, g_0_xxxxyy_0_xxxxzzz_0, g_0_xxxxyy_0_xxxxzzz_1, g_0_xxxxyy_0_xxxzzzz_0, g_0_xxxxyy_0_xxxzzzz_1, g_0_xxxxyy_0_xxzzzzz_0, g_0_xxxxyy_0_xxzzzzz_1, g_0_xxxxyy_0_xzzzzzz_0, g_0_xxxxyy_0_xzzzzzz_1, g_0_xxxxyyy_0_xxxxxxx_0, g_0_xxxxyyy_0_xxxxxxy_0, g_0_xxxxyyy_0_xxxxxxz_0, g_0_xxxxyyy_0_xxxxxyy_0, g_0_xxxxyyy_0_xxxxxyz_0, g_0_xxxxyyy_0_xxxxxzz_0, g_0_xxxxyyy_0_xxxxyyy_0, g_0_xxxxyyy_0_xxxxyyz_0, g_0_xxxxyyy_0_xxxxyzz_0, g_0_xxxxyyy_0_xxxxzzz_0, g_0_xxxxyyy_0_xxxyyyy_0, g_0_xxxxyyy_0_xxxyyyz_0, g_0_xxxxyyy_0_xxxyyzz_0, g_0_xxxxyyy_0_xxxyzzz_0, g_0_xxxxyyy_0_xxxzzzz_0, g_0_xxxxyyy_0_xxyyyyy_0, g_0_xxxxyyy_0_xxyyyyz_0, g_0_xxxxyyy_0_xxyyyzz_0, g_0_xxxxyyy_0_xxyyzzz_0, g_0_xxxxyyy_0_xxyzzzz_0, g_0_xxxxyyy_0_xxzzzzz_0, g_0_xxxxyyy_0_xyyyyyy_0, g_0_xxxxyyy_0_xyyyyyz_0, g_0_xxxxyyy_0_xyyyyzz_0, g_0_xxxxyyy_0_xyyyzzz_0, g_0_xxxxyyy_0_xyyzzzz_0, g_0_xxxxyyy_0_xyzzzzz_0, g_0_xxxxyyy_0_xzzzzzz_0, g_0_xxxxyyy_0_yyyyyyy_0, g_0_xxxxyyy_0_yyyyyyz_0, g_0_xxxxyyy_0_yyyyyzz_0, g_0_xxxxyyy_0_yyyyzzz_0, g_0_xxxxyyy_0_yyyzzzz_0, g_0_xxxxyyy_0_yyzzzzz_0, g_0_xxxxyyy_0_yzzzzzz_0, g_0_xxxxyyy_0_zzzzzzz_0, g_0_xxxyyy_0_xxxxxxy_0, g_0_xxxyyy_0_xxxxxxy_1, g_0_xxxyyy_0_xxxxxy_1, g_0_xxxyyy_0_xxxxxyy_0, g_0_xxxyyy_0_xxxxxyy_1, g_0_xxxyyy_0_xxxxxyz_0, g_0_xxxyyy_0_xxxxxyz_1, g_0_xxxyyy_0_xxxxyy_1, g_0_xxxyyy_0_xxxxyyy_0, g_0_xxxyyy_0_xxxxyyy_1, g_0_xxxyyy_0_xxxxyyz_0, g_0_xxxyyy_0_xxxxyyz_1, g_0_xxxyyy_0_xxxxyz_1, g_0_xxxyyy_0_xxxxyzz_0, g_0_xxxyyy_0_xxxxyzz_1, g_0_xxxyyy_0_xxxyyy_1, g_0_xxxyyy_0_xxxyyyy_0, g_0_xxxyyy_0_xxxyyyy_1, g_0_xxxyyy_0_xxxyyyz_0, g_0_xxxyyy_0_xxxyyyz_1, g_0_xxxyyy_0_xxxyyz_1, g_0_xxxyyy_0_xxxyyzz_0, g_0_xxxyyy_0_xxxyyzz_1, g_0_xxxyyy_0_xxxyzz_1, g_0_xxxyyy_0_xxxyzzz_0, g_0_xxxyyy_0_xxxyzzz_1, g_0_xxxyyy_0_xxyyyy_1, g_0_xxxyyy_0_xxyyyyy_0, g_0_xxxyyy_0_xxyyyyy_1, g_0_xxxyyy_0_xxyyyyz_0, g_0_xxxyyy_0_xxyyyyz_1, g_0_xxxyyy_0_xxyyyz_1, g_0_xxxyyy_0_xxyyyzz_0, g_0_xxxyyy_0_xxyyyzz_1, g_0_xxxyyy_0_xxyyzz_1, g_0_xxxyyy_0_xxyyzzz_0, g_0_xxxyyy_0_xxyyzzz_1, g_0_xxxyyy_0_xxyzzz_1, g_0_xxxyyy_0_xxyzzzz_0, g_0_xxxyyy_0_xxyzzzz_1, g_0_xxxyyy_0_xyyyyy_1, g_0_xxxyyy_0_xyyyyyy_0, g_0_xxxyyy_0_xyyyyyy_1, g_0_xxxyyy_0_xyyyyyz_0, g_0_xxxyyy_0_xyyyyyz_1, g_0_xxxyyy_0_xyyyyz_1, g_0_xxxyyy_0_xyyyyzz_0, g_0_xxxyyy_0_xyyyyzz_1, g_0_xxxyyy_0_xyyyzz_1, g_0_xxxyyy_0_xyyyzzz_0, g_0_xxxyyy_0_xyyyzzz_1, g_0_xxxyyy_0_xyyzzz_1, g_0_xxxyyy_0_xyyzzzz_0, g_0_xxxyyy_0_xyyzzzz_1, g_0_xxxyyy_0_xyzzzz_1, g_0_xxxyyy_0_xyzzzzz_0, g_0_xxxyyy_0_xyzzzzz_1, g_0_xxxyyy_0_yyyyyy_1, g_0_xxxyyy_0_yyyyyyy_0, g_0_xxxyyy_0_yyyyyyy_1, g_0_xxxyyy_0_yyyyyyz_0, g_0_xxxyyy_0_yyyyyyz_1, g_0_xxxyyy_0_yyyyyz_1, g_0_xxxyyy_0_yyyyyzz_0, g_0_xxxyyy_0_yyyyyzz_1, g_0_xxxyyy_0_yyyyzz_1, g_0_xxxyyy_0_yyyyzzz_0, g_0_xxxyyy_0_yyyyzzz_1, g_0_xxxyyy_0_yyyzzz_1, g_0_xxxyyy_0_yyyzzzz_0, g_0_xxxyyy_0_yyyzzzz_1, g_0_xxxyyy_0_yyzzzz_1, g_0_xxxyyy_0_yyzzzzz_0, g_0_xxxyyy_0_yyzzzzz_1, g_0_xxxyyy_0_yzzzzz_1, g_0_xxxyyy_0_yzzzzzz_0, g_0_xxxyyy_0_yzzzzzz_1, g_0_xxxyyy_0_zzzzzzz_0, g_0_xxxyyy_0_zzzzzzz_1, g_0_xxyyy_0_xxxxxxy_0, g_0_xxyyy_0_xxxxxxy_1, g_0_xxyyy_0_xxxxxyy_0, g_0_xxyyy_0_xxxxxyy_1, g_0_xxyyy_0_xxxxxyz_0, g_0_xxyyy_0_xxxxxyz_1, g_0_xxyyy_0_xxxxyyy_0, g_0_xxyyy_0_xxxxyyy_1, g_0_xxyyy_0_xxxxyyz_0, g_0_xxyyy_0_xxxxyyz_1, g_0_xxyyy_0_xxxxyzz_0, g_0_xxyyy_0_xxxxyzz_1, g_0_xxyyy_0_xxxyyyy_0, g_0_xxyyy_0_xxxyyyy_1, g_0_xxyyy_0_xxxyyyz_0, g_0_xxyyy_0_xxxyyyz_1, g_0_xxyyy_0_xxxyyzz_0, g_0_xxyyy_0_xxxyyzz_1, g_0_xxyyy_0_xxxyzzz_0, g_0_xxyyy_0_xxxyzzz_1, g_0_xxyyy_0_xxyyyyy_0, g_0_xxyyy_0_xxyyyyy_1, g_0_xxyyy_0_xxyyyyz_0, g_0_xxyyy_0_xxyyyyz_1, g_0_xxyyy_0_xxyyyzz_0, g_0_xxyyy_0_xxyyyzz_1, g_0_xxyyy_0_xxyyzzz_0, g_0_xxyyy_0_xxyyzzz_1, g_0_xxyyy_0_xxyzzzz_0, g_0_xxyyy_0_xxyzzzz_1, g_0_xxyyy_0_xyyyyyy_0, g_0_xxyyy_0_xyyyyyy_1, g_0_xxyyy_0_xyyyyyz_0, g_0_xxyyy_0_xyyyyyz_1, g_0_xxyyy_0_xyyyyzz_0, g_0_xxyyy_0_xyyyyzz_1, g_0_xxyyy_0_xyyyzzz_0, g_0_xxyyy_0_xyyyzzz_1, g_0_xxyyy_0_xyyzzzz_0, g_0_xxyyy_0_xyyzzzz_1, g_0_xxyyy_0_xyzzzzz_0, g_0_xxyyy_0_xyzzzzz_1, g_0_xxyyy_0_yyyyyyy_0, g_0_xxyyy_0_yyyyyyy_1, g_0_xxyyy_0_yyyyyyz_0, g_0_xxyyy_0_yyyyyyz_1, g_0_xxyyy_0_yyyyyzz_0, g_0_xxyyy_0_yyyyyzz_1, g_0_xxyyy_0_yyyyzzz_0, g_0_xxyyy_0_yyyyzzz_1, g_0_xxyyy_0_yyyzzzz_0, g_0_xxyyy_0_yyyzzzz_1, g_0_xxyyy_0_yyzzzzz_0, g_0_xxyyy_0_yyzzzzz_1, g_0_xxyyy_0_yzzzzzz_0, g_0_xxyyy_0_yzzzzzz_1, g_0_xxyyy_0_zzzzzzz_0, g_0_xxyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyy_0_xxxxxxx_0[i] = 2.0 * g_0_xxxxy_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxxxxx_0[i] * pb_y + g_0_xxxxyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxxxxy_0[i] = 3.0 * g_0_xxyyy_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxxxy_0[i] * pb_x + g_0_xxxyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxxxz_0[i] = 2.0 * g_0_xxxxy_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxxxxz_0[i] * pb_y + g_0_xxxxyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxxxyy_0[i] = 3.0 * g_0_xxyyy_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxxyy_0[i] * pb_x + g_0_xxxyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxxyz_0[i] = 3.0 * g_0_xxyyy_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxxyz_0[i] * pb_x + g_0_xxxyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxxzz_0[i] = 2.0 * g_0_xxxxy_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxxxzz_0[i] * pb_y + g_0_xxxxyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxxyyy_0[i] = 3.0 * g_0_xxyyy_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxyyy_0[i] * pb_x + g_0_xxxyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxyyz_0[i] = 3.0 * g_0_xxyyy_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxyyz_0[i] * pb_x + g_0_xxxyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxyzz_0[i] = 3.0 * g_0_xxyyy_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxyzz_0[i] * pb_x + g_0_xxxyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxzzz_0[i] = 2.0 * g_0_xxxxy_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxxzzz_0[i] * pb_y + g_0_xxxxyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxyyyy_0[i] = 3.0 * g_0_xxyyy_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyyyy_0[i] * pb_x + g_0_xxxyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxyyyz_0[i] = 3.0 * g_0_xxyyy_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyyyz_0[i] * pb_x + g_0_xxxyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxyyzz_0[i] = 3.0 * g_0_xxyyy_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyyzz_0[i] * pb_x + g_0_xxxyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxyzzz_0[i] = 3.0 * g_0_xxyyy_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyzzz_0[i] * pb_x + g_0_xxxyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxzzzz_0[i] = 2.0 * g_0_xxxxy_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxxzzzz_0[i] * pb_y + g_0_xxxxyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxyyyyy_0[i] = 3.0 * g_0_xxyyy_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyyyy_0[i] * pb_x + g_0_xxxyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyyyyz_0[i] = 3.0 * g_0_xxyyy_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyyyz_0[i] * pb_x + g_0_xxxyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyyyzz_0[i] = 3.0 * g_0_xxyyy_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyyzz_0[i] * pb_x + g_0_xxxyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyyzzz_0[i] = 3.0 * g_0_xxyyy_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyzzz_0[i] * pb_x + g_0_xxxyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyzzzz_0[i] = 3.0 * g_0_xxyyy_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyzzzz_0[i] * pb_x + g_0_xxxyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxzzzzz_0[i] = 2.0 * g_0_xxxxy_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xxzzzzz_0[i] * pb_y + g_0_xxxxyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xyyyyyy_0[i] = 3.0 * g_0_xxyyy_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyyyy_0[i] * pb_x + g_0_xxxyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyyyyz_0[i] = 3.0 * g_0_xxyyy_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyyyz_0[i] * pb_x + g_0_xxxyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyyyzz_0[i] = 3.0 * g_0_xxyyy_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyyzz_0[i] * pb_x + g_0_xxxyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyyzzz_0[i] = 3.0 * g_0_xxyyy_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyzzz_0[i] * pb_x + g_0_xxxyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyzzzz_0[i] = 3.0 * g_0_xxyyy_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyzzzz_0[i] * pb_x + g_0_xxxyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyzzzzz_0[i] = 3.0 * g_0_xxyyy_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyzzzzz_0[i] * pb_x + g_0_xxxyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xzzzzzz_0[i] = 2.0 * g_0_xxxxy_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxyy_0_xzzzzzz_0[i] * pb_y + g_0_xxxxyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_yyyyyyy_0[i] = 3.0 * g_0_xxyyy_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyyyy_0[i] * pb_x + g_0_xxxyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyyyyz_0[i] = 3.0 * g_0_xxyyy_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyyyz_0[i] * pb_x + g_0_xxxyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyyyzz_0[i] = 3.0 * g_0_xxyyy_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyyzz_0[i] * pb_x + g_0_xxxyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyyzzz_0[i] = 3.0 * g_0_xxyyy_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyyzzz_0[i] * pb_x + g_0_xxxyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyzzzz_0[i] = 3.0 * g_0_xxyyy_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyyzzzz_0[i] * pb_x + g_0_xxxyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyzzzzz_0[i] = 3.0 * g_0_xxyyy_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yyzzzzz_0[i] * pb_x + g_0_xxxyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yzzzzzz_0[i] = 3.0 * g_0_xxyyy_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_yzzzzzz_0[i] * pb_x + g_0_xxxyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_zzzzzzz_0[i] = 3.0 * g_0_xxyyy_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_zzzzzzz_0[i] * pb_x + g_0_xxxyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 252-288 components of targeted buffer : SKSK

    auto g_0_xxxxyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 252);

    auto g_0_xxxxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 253);

    auto g_0_xxxxyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 254);

    auto g_0_xxxxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 255);

    auto g_0_xxxxyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 256);

    auto g_0_xxxxyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 257);

    auto g_0_xxxxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 258);

    auto g_0_xxxxyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 259);

    auto g_0_xxxxyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 260);

    auto g_0_xxxxyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 261);

    auto g_0_xxxxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 262);

    auto g_0_xxxxyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 263);

    auto g_0_xxxxyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 264);

    auto g_0_xxxxyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 265);

    auto g_0_xxxxyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 266);

    auto g_0_xxxxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 267);

    auto g_0_xxxxyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 268);

    auto g_0_xxxxyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 269);

    auto g_0_xxxxyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 270);

    auto g_0_xxxxyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 271);

    auto g_0_xxxxyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 272);

    auto g_0_xxxxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 273);

    auto g_0_xxxxyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 274);

    auto g_0_xxxxyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 275);

    auto g_0_xxxxyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 276);

    auto g_0_xxxxyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 277);

    auto g_0_xxxxyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 278);

    auto g_0_xxxxyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 279);

    auto g_0_xxxxyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 280);

    auto g_0_xxxxyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 281);

    auto g_0_xxxxyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 282);

    auto g_0_xxxxyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 283);

    auto g_0_xxxxyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 284);

    auto g_0_xxxxyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 285);

    auto g_0_xxxxyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 286);

    auto g_0_xxxxyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 287);

    #pragma omp simd aligned(g_0_xxxxyy_0_xxxxxx_1, g_0_xxxxyy_0_xxxxxxx_0, g_0_xxxxyy_0_xxxxxxx_1, g_0_xxxxyy_0_xxxxxxy_0, g_0_xxxxyy_0_xxxxxxy_1, g_0_xxxxyy_0_xxxxxxz_0, g_0_xxxxyy_0_xxxxxxz_1, g_0_xxxxyy_0_xxxxxy_1, g_0_xxxxyy_0_xxxxxyy_0, g_0_xxxxyy_0_xxxxxyy_1, g_0_xxxxyy_0_xxxxxyz_0, g_0_xxxxyy_0_xxxxxyz_1, g_0_xxxxyy_0_xxxxxz_1, g_0_xxxxyy_0_xxxxxzz_0, g_0_xxxxyy_0_xxxxxzz_1, g_0_xxxxyy_0_xxxxyy_1, g_0_xxxxyy_0_xxxxyyy_0, g_0_xxxxyy_0_xxxxyyy_1, g_0_xxxxyy_0_xxxxyyz_0, g_0_xxxxyy_0_xxxxyyz_1, g_0_xxxxyy_0_xxxxyz_1, g_0_xxxxyy_0_xxxxyzz_0, g_0_xxxxyy_0_xxxxyzz_1, g_0_xxxxyy_0_xxxxzz_1, g_0_xxxxyy_0_xxxxzzz_0, g_0_xxxxyy_0_xxxxzzz_1, g_0_xxxxyy_0_xxxyyy_1, g_0_xxxxyy_0_xxxyyyy_0, g_0_xxxxyy_0_xxxyyyy_1, g_0_xxxxyy_0_xxxyyyz_0, g_0_xxxxyy_0_xxxyyyz_1, g_0_xxxxyy_0_xxxyyz_1, g_0_xxxxyy_0_xxxyyzz_0, g_0_xxxxyy_0_xxxyyzz_1, g_0_xxxxyy_0_xxxyzz_1, g_0_xxxxyy_0_xxxyzzz_0, g_0_xxxxyy_0_xxxyzzz_1, g_0_xxxxyy_0_xxxzzz_1, g_0_xxxxyy_0_xxxzzzz_0, g_0_xxxxyy_0_xxxzzzz_1, g_0_xxxxyy_0_xxyyyy_1, g_0_xxxxyy_0_xxyyyyy_0, g_0_xxxxyy_0_xxyyyyy_1, g_0_xxxxyy_0_xxyyyyz_0, g_0_xxxxyy_0_xxyyyyz_1, g_0_xxxxyy_0_xxyyyz_1, g_0_xxxxyy_0_xxyyyzz_0, g_0_xxxxyy_0_xxyyyzz_1, g_0_xxxxyy_0_xxyyzz_1, g_0_xxxxyy_0_xxyyzzz_0, g_0_xxxxyy_0_xxyyzzz_1, g_0_xxxxyy_0_xxyzzz_1, g_0_xxxxyy_0_xxyzzzz_0, g_0_xxxxyy_0_xxyzzzz_1, g_0_xxxxyy_0_xxzzzz_1, g_0_xxxxyy_0_xxzzzzz_0, g_0_xxxxyy_0_xxzzzzz_1, g_0_xxxxyy_0_xyyyyy_1, g_0_xxxxyy_0_xyyyyyy_0, g_0_xxxxyy_0_xyyyyyy_1, g_0_xxxxyy_0_xyyyyyz_0, g_0_xxxxyy_0_xyyyyyz_1, g_0_xxxxyy_0_xyyyyz_1, g_0_xxxxyy_0_xyyyyzz_0, g_0_xxxxyy_0_xyyyyzz_1, g_0_xxxxyy_0_xyyyzz_1, g_0_xxxxyy_0_xyyyzzz_0, g_0_xxxxyy_0_xyyyzzz_1, g_0_xxxxyy_0_xyyzzz_1, g_0_xxxxyy_0_xyyzzzz_0, g_0_xxxxyy_0_xyyzzzz_1, g_0_xxxxyy_0_xyzzzz_1, g_0_xxxxyy_0_xyzzzzz_0, g_0_xxxxyy_0_xyzzzzz_1, g_0_xxxxyy_0_xzzzzz_1, g_0_xxxxyy_0_xzzzzzz_0, g_0_xxxxyy_0_xzzzzzz_1, g_0_xxxxyy_0_yyyyyy_1, g_0_xxxxyy_0_yyyyyyy_0, g_0_xxxxyy_0_yyyyyyy_1, g_0_xxxxyy_0_yyyyyyz_0, g_0_xxxxyy_0_yyyyyyz_1, g_0_xxxxyy_0_yyyyyz_1, g_0_xxxxyy_0_yyyyyzz_0, g_0_xxxxyy_0_yyyyyzz_1, g_0_xxxxyy_0_yyyyzz_1, g_0_xxxxyy_0_yyyyzzz_0, g_0_xxxxyy_0_yyyyzzz_1, g_0_xxxxyy_0_yyyzzz_1, g_0_xxxxyy_0_yyyzzzz_0, g_0_xxxxyy_0_yyyzzzz_1, g_0_xxxxyy_0_yyzzzz_1, g_0_xxxxyy_0_yyzzzzz_0, g_0_xxxxyy_0_yyzzzzz_1, g_0_xxxxyy_0_yzzzzz_1, g_0_xxxxyy_0_yzzzzzz_0, g_0_xxxxyy_0_yzzzzzz_1, g_0_xxxxyy_0_zzzzzz_1, g_0_xxxxyy_0_zzzzzzz_0, g_0_xxxxyy_0_zzzzzzz_1, g_0_xxxxyyz_0_xxxxxxx_0, g_0_xxxxyyz_0_xxxxxxy_0, g_0_xxxxyyz_0_xxxxxxz_0, g_0_xxxxyyz_0_xxxxxyy_0, g_0_xxxxyyz_0_xxxxxyz_0, g_0_xxxxyyz_0_xxxxxzz_0, g_0_xxxxyyz_0_xxxxyyy_0, g_0_xxxxyyz_0_xxxxyyz_0, g_0_xxxxyyz_0_xxxxyzz_0, g_0_xxxxyyz_0_xxxxzzz_0, g_0_xxxxyyz_0_xxxyyyy_0, g_0_xxxxyyz_0_xxxyyyz_0, g_0_xxxxyyz_0_xxxyyzz_0, g_0_xxxxyyz_0_xxxyzzz_0, g_0_xxxxyyz_0_xxxzzzz_0, g_0_xxxxyyz_0_xxyyyyy_0, g_0_xxxxyyz_0_xxyyyyz_0, g_0_xxxxyyz_0_xxyyyzz_0, g_0_xxxxyyz_0_xxyyzzz_0, g_0_xxxxyyz_0_xxyzzzz_0, g_0_xxxxyyz_0_xxzzzzz_0, g_0_xxxxyyz_0_xyyyyyy_0, g_0_xxxxyyz_0_xyyyyyz_0, g_0_xxxxyyz_0_xyyyyzz_0, g_0_xxxxyyz_0_xyyyzzz_0, g_0_xxxxyyz_0_xyyzzzz_0, g_0_xxxxyyz_0_xyzzzzz_0, g_0_xxxxyyz_0_xzzzzzz_0, g_0_xxxxyyz_0_yyyyyyy_0, g_0_xxxxyyz_0_yyyyyyz_0, g_0_xxxxyyz_0_yyyyyzz_0, g_0_xxxxyyz_0_yyyyzzz_0, g_0_xxxxyyz_0_yyyzzzz_0, g_0_xxxxyyz_0_yyzzzzz_0, g_0_xxxxyyz_0_yzzzzzz_0, g_0_xxxxyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyz_0_xxxxxxx_0[i] = g_0_xxxxyy_0_xxxxxxx_0[i] * pb_z + g_0_xxxxyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxxxy_0[i] = g_0_xxxxyy_0_xxxxxxy_0[i] * pb_z + g_0_xxxxyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxxxz_0[i] = g_0_xxxxyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxxxz_0[i] * pb_z + g_0_xxxxyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxxyy_0[i] = g_0_xxxxyy_0_xxxxxyy_0[i] * pb_z + g_0_xxxxyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxxyz_0[i] = g_0_xxxxyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxxyz_0[i] * pb_z + g_0_xxxxyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxxyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxxzz_0[i] * pb_z + g_0_xxxxyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxyyy_0[i] = g_0_xxxxyy_0_xxxxyyy_0[i] * pb_z + g_0_xxxxyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxyyz_0[i] = g_0_xxxxyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxyyz_0[i] * pb_z + g_0_xxxxyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxxyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxyzz_0[i] * pb_z + g_0_xxxxyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxzzz_0[i] * pb_z + g_0_xxxxyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyyyy_0[i] = g_0_xxxxyy_0_xxxyyyy_0[i] * pb_z + g_0_xxxxyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyyyz_0[i] = g_0_xxxxyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyyyz_0[i] * pb_z + g_0_xxxxyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyyzz_0[i] * pb_z + g_0_xxxxyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyzzz_0[i] * pb_z + g_0_xxxxyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxxyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxzzzz_0[i] * pb_z + g_0_xxxxyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyyyy_0[i] = g_0_xxxxyy_0_xxyyyyy_0[i] * pb_z + g_0_xxxxyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyyyz_0[i] = g_0_xxxxyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyyyz_0[i] * pb_z + g_0_xxxxyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxxyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyyzz_0[i] * pb_z + g_0_xxxxyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyzzz_0[i] * pb_z + g_0_xxxxyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxxyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyzzzz_0[i] * pb_z + g_0_xxxxyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxzzzzz_0[i] * pb_z + g_0_xxxxyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyyyy_0[i] = g_0_xxxxyy_0_xyyyyyy_0[i] * pb_z + g_0_xxxxyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyyyz_0[i] = g_0_xxxxyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyyyz_0[i] * pb_z + g_0_xxxxyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxxyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyyzz_0[i] * pb_z + g_0_xxxxyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyzzz_0[i] * pb_z + g_0_xxxxyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxxyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyzzzz_0[i] * pb_z + g_0_xxxxyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyzzzzz_0[i] * pb_z + g_0_xxxxyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxxyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xzzzzzz_0[i] * pb_z + g_0_xxxxyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyyyy_0[i] = g_0_xxxxyy_0_yyyyyyy_0[i] * pb_z + g_0_xxxxyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyyyz_0[i] = g_0_xxxxyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyyyyz_0[i] * pb_z + g_0_xxxxyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxxyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyyyzz_0[i] * pb_z + g_0_xxxxyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxxyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyyzzz_0[i] * pb_z + g_0_xxxxyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxxyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyzzzz_0[i] * pb_z + g_0_xxxxyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyzzzzz_0[i] * pb_z + g_0_xxxxyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxxyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yzzzzzz_0[i] * pb_z + g_0_xxxxyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxxyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_zzzzzzz_0[i] * pb_z + g_0_xxxxyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 288-324 components of targeted buffer : SKSK

    auto g_0_xxxxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 288);

    auto g_0_xxxxyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 289);

    auto g_0_xxxxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 290);

    auto g_0_xxxxyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 291);

    auto g_0_xxxxyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 292);

    auto g_0_xxxxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 293);

    auto g_0_xxxxyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 294);

    auto g_0_xxxxyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 295);

    auto g_0_xxxxyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 296);

    auto g_0_xxxxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 297);

    auto g_0_xxxxyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 298);

    auto g_0_xxxxyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 299);

    auto g_0_xxxxyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 300);

    auto g_0_xxxxyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 301);

    auto g_0_xxxxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 302);

    auto g_0_xxxxyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 303);

    auto g_0_xxxxyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 304);

    auto g_0_xxxxyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 305);

    auto g_0_xxxxyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 306);

    auto g_0_xxxxyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 307);

    auto g_0_xxxxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 308);

    auto g_0_xxxxyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 309);

    auto g_0_xxxxyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 310);

    auto g_0_xxxxyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 311);

    auto g_0_xxxxyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 312);

    auto g_0_xxxxyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 313);

    auto g_0_xxxxyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 314);

    auto g_0_xxxxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 315);

    auto g_0_xxxxyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 316);

    auto g_0_xxxxyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 317);

    auto g_0_xxxxyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 318);

    auto g_0_xxxxyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 319);

    auto g_0_xxxxyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 320);

    auto g_0_xxxxyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 321);

    auto g_0_xxxxyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 322);

    auto g_0_xxxxyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 323);

    #pragma omp simd aligned(g_0_xxxxyzz_0_xxxxxxx_0, g_0_xxxxyzz_0_xxxxxxy_0, g_0_xxxxyzz_0_xxxxxxz_0, g_0_xxxxyzz_0_xxxxxyy_0, g_0_xxxxyzz_0_xxxxxyz_0, g_0_xxxxyzz_0_xxxxxzz_0, g_0_xxxxyzz_0_xxxxyyy_0, g_0_xxxxyzz_0_xxxxyyz_0, g_0_xxxxyzz_0_xxxxyzz_0, g_0_xxxxyzz_0_xxxxzzz_0, g_0_xxxxyzz_0_xxxyyyy_0, g_0_xxxxyzz_0_xxxyyyz_0, g_0_xxxxyzz_0_xxxyyzz_0, g_0_xxxxyzz_0_xxxyzzz_0, g_0_xxxxyzz_0_xxxzzzz_0, g_0_xxxxyzz_0_xxyyyyy_0, g_0_xxxxyzz_0_xxyyyyz_0, g_0_xxxxyzz_0_xxyyyzz_0, g_0_xxxxyzz_0_xxyyzzz_0, g_0_xxxxyzz_0_xxyzzzz_0, g_0_xxxxyzz_0_xxzzzzz_0, g_0_xxxxyzz_0_xyyyyyy_0, g_0_xxxxyzz_0_xyyyyyz_0, g_0_xxxxyzz_0_xyyyyzz_0, g_0_xxxxyzz_0_xyyyzzz_0, g_0_xxxxyzz_0_xyyzzzz_0, g_0_xxxxyzz_0_xyzzzzz_0, g_0_xxxxyzz_0_xzzzzzz_0, g_0_xxxxyzz_0_yyyyyyy_0, g_0_xxxxyzz_0_yyyyyyz_0, g_0_xxxxyzz_0_yyyyyzz_0, g_0_xxxxyzz_0_yyyyzzz_0, g_0_xxxxyzz_0_yyyzzzz_0, g_0_xxxxyzz_0_yyzzzzz_0, g_0_xxxxyzz_0_yzzzzzz_0, g_0_xxxxyzz_0_zzzzzzz_0, g_0_xxxxzz_0_xxxxxx_1, g_0_xxxxzz_0_xxxxxxx_0, g_0_xxxxzz_0_xxxxxxx_1, g_0_xxxxzz_0_xxxxxxy_0, g_0_xxxxzz_0_xxxxxxy_1, g_0_xxxxzz_0_xxxxxxz_0, g_0_xxxxzz_0_xxxxxxz_1, g_0_xxxxzz_0_xxxxxy_1, g_0_xxxxzz_0_xxxxxyy_0, g_0_xxxxzz_0_xxxxxyy_1, g_0_xxxxzz_0_xxxxxyz_0, g_0_xxxxzz_0_xxxxxyz_1, g_0_xxxxzz_0_xxxxxz_1, g_0_xxxxzz_0_xxxxxzz_0, g_0_xxxxzz_0_xxxxxzz_1, g_0_xxxxzz_0_xxxxyy_1, g_0_xxxxzz_0_xxxxyyy_0, g_0_xxxxzz_0_xxxxyyy_1, g_0_xxxxzz_0_xxxxyyz_0, g_0_xxxxzz_0_xxxxyyz_1, g_0_xxxxzz_0_xxxxyz_1, g_0_xxxxzz_0_xxxxyzz_0, g_0_xxxxzz_0_xxxxyzz_1, g_0_xxxxzz_0_xxxxzz_1, g_0_xxxxzz_0_xxxxzzz_0, g_0_xxxxzz_0_xxxxzzz_1, g_0_xxxxzz_0_xxxyyy_1, g_0_xxxxzz_0_xxxyyyy_0, g_0_xxxxzz_0_xxxyyyy_1, g_0_xxxxzz_0_xxxyyyz_0, g_0_xxxxzz_0_xxxyyyz_1, g_0_xxxxzz_0_xxxyyz_1, g_0_xxxxzz_0_xxxyyzz_0, g_0_xxxxzz_0_xxxyyzz_1, g_0_xxxxzz_0_xxxyzz_1, g_0_xxxxzz_0_xxxyzzz_0, g_0_xxxxzz_0_xxxyzzz_1, g_0_xxxxzz_0_xxxzzz_1, g_0_xxxxzz_0_xxxzzzz_0, g_0_xxxxzz_0_xxxzzzz_1, g_0_xxxxzz_0_xxyyyy_1, g_0_xxxxzz_0_xxyyyyy_0, g_0_xxxxzz_0_xxyyyyy_1, g_0_xxxxzz_0_xxyyyyz_0, g_0_xxxxzz_0_xxyyyyz_1, g_0_xxxxzz_0_xxyyyz_1, g_0_xxxxzz_0_xxyyyzz_0, g_0_xxxxzz_0_xxyyyzz_1, g_0_xxxxzz_0_xxyyzz_1, g_0_xxxxzz_0_xxyyzzz_0, g_0_xxxxzz_0_xxyyzzz_1, g_0_xxxxzz_0_xxyzzz_1, g_0_xxxxzz_0_xxyzzzz_0, g_0_xxxxzz_0_xxyzzzz_1, g_0_xxxxzz_0_xxzzzz_1, g_0_xxxxzz_0_xxzzzzz_0, g_0_xxxxzz_0_xxzzzzz_1, g_0_xxxxzz_0_xyyyyy_1, g_0_xxxxzz_0_xyyyyyy_0, g_0_xxxxzz_0_xyyyyyy_1, g_0_xxxxzz_0_xyyyyyz_0, g_0_xxxxzz_0_xyyyyyz_1, g_0_xxxxzz_0_xyyyyz_1, g_0_xxxxzz_0_xyyyyzz_0, g_0_xxxxzz_0_xyyyyzz_1, g_0_xxxxzz_0_xyyyzz_1, g_0_xxxxzz_0_xyyyzzz_0, g_0_xxxxzz_0_xyyyzzz_1, g_0_xxxxzz_0_xyyzzz_1, g_0_xxxxzz_0_xyyzzzz_0, g_0_xxxxzz_0_xyyzzzz_1, g_0_xxxxzz_0_xyzzzz_1, g_0_xxxxzz_0_xyzzzzz_0, g_0_xxxxzz_0_xyzzzzz_1, g_0_xxxxzz_0_xzzzzz_1, g_0_xxxxzz_0_xzzzzzz_0, g_0_xxxxzz_0_xzzzzzz_1, g_0_xxxxzz_0_yyyyyy_1, g_0_xxxxzz_0_yyyyyyy_0, g_0_xxxxzz_0_yyyyyyy_1, g_0_xxxxzz_0_yyyyyyz_0, g_0_xxxxzz_0_yyyyyyz_1, g_0_xxxxzz_0_yyyyyz_1, g_0_xxxxzz_0_yyyyyzz_0, g_0_xxxxzz_0_yyyyyzz_1, g_0_xxxxzz_0_yyyyzz_1, g_0_xxxxzz_0_yyyyzzz_0, g_0_xxxxzz_0_yyyyzzz_1, g_0_xxxxzz_0_yyyzzz_1, g_0_xxxxzz_0_yyyzzzz_0, g_0_xxxxzz_0_yyyzzzz_1, g_0_xxxxzz_0_yyzzzz_1, g_0_xxxxzz_0_yyzzzzz_0, g_0_xxxxzz_0_yyzzzzz_1, g_0_xxxxzz_0_yzzzzz_1, g_0_xxxxzz_0_yzzzzzz_0, g_0_xxxxzz_0_yzzzzzz_1, g_0_xxxxzz_0_zzzzzz_1, g_0_xxxxzz_0_zzzzzzz_0, g_0_xxxxzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzz_0_xxxxxxx_0[i] = g_0_xxxxzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxxzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxxxy_0[i] = g_0_xxxxzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxxxy_0[i] * pb_y + g_0_xxxxzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxxxz_0[i] = g_0_xxxxzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxxzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxxzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxxyy_0[i] * pb_y + g_0_xxxxzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxxyz_0[i] = g_0_xxxxzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxxyz_0[i] * pb_y + g_0_xxxxzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxxzz_0[i] = g_0_xxxxzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxxzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxyyy_0[i] * pb_y + g_0_xxxxzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxyyz_0[i] * pb_y + g_0_xxxxzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxyzz_0[i] = g_0_xxxxzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxyzz_0[i] * pb_y + g_0_xxxxzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxzzz_0[i] = g_0_xxxxzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxxzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxxxzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyyyy_0[i] * pb_y + g_0_xxxxzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyyyz_0[i] * pb_y + g_0_xxxxzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyyzz_0[i] * pb_y + g_0_xxxxzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyzzz_0[i] = g_0_xxxxzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyzzz_0[i] * pb_y + g_0_xxxxzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxzzzz_0[i] = g_0_xxxxzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxxzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxxxzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyyyy_0[i] * pb_y + g_0_xxxxzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyyyz_0[i] * pb_y + g_0_xxxxzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyyzz_0[i] * pb_y + g_0_xxxxzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyzzz_0[i] * pb_y + g_0_xxxxzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyzzzz_0[i] = g_0_xxxxzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyzzzz_0[i] * pb_y + g_0_xxxxzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxzzzzz_0[i] = g_0_xxxxzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxxzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxxxzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyyyy_0[i] * pb_y + g_0_xxxxzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyyyz_0[i] * pb_y + g_0_xxxxzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyyzz_0[i] * pb_y + g_0_xxxxzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyzzz_0[i] * pb_y + g_0_xxxxzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyzzzz_0[i] * pb_y + g_0_xxxxzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyzzzzz_0[i] = g_0_xxxxzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyzzzzz_0[i] * pb_y + g_0_xxxxzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xzzzzzz_0[i] = g_0_xxxxzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxxzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxxxzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyyyy_0[i] * pb_y + g_0_xxxxzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyyyz_0[i] * pb_y + g_0_xxxxzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyyzz_0[i] * pb_y + g_0_xxxxzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyzzz_0[i] * pb_y + g_0_xxxxzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyzzzz_0[i] * pb_y + g_0_xxxxzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyzzzzz_0[i] * pb_y + g_0_xxxxzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yzzzzzz_0[i] = g_0_xxxxzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yzzzzzz_0[i] * pb_y + g_0_xxxxzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_zzzzzzz_0[i] = g_0_xxxxzz_0_zzzzzzz_0[i] * pb_y + g_0_xxxxzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 324-360 components of targeted buffer : SKSK

    auto g_0_xxxxzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 324);

    auto g_0_xxxxzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 325);

    auto g_0_xxxxzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 326);

    auto g_0_xxxxzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 327);

    auto g_0_xxxxzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 328);

    auto g_0_xxxxzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 329);

    auto g_0_xxxxzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 330);

    auto g_0_xxxxzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 331);

    auto g_0_xxxxzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 332);

    auto g_0_xxxxzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 333);

    auto g_0_xxxxzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 334);

    auto g_0_xxxxzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 335);

    auto g_0_xxxxzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 336);

    auto g_0_xxxxzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 337);

    auto g_0_xxxxzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 338);

    auto g_0_xxxxzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 339);

    auto g_0_xxxxzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 340);

    auto g_0_xxxxzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 341);

    auto g_0_xxxxzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 342);

    auto g_0_xxxxzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 343);

    auto g_0_xxxxzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 344);

    auto g_0_xxxxzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 345);

    auto g_0_xxxxzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 346);

    auto g_0_xxxxzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 347);

    auto g_0_xxxxzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 348);

    auto g_0_xxxxzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 349);

    auto g_0_xxxxzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 350);

    auto g_0_xxxxzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 351);

    auto g_0_xxxxzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 352);

    auto g_0_xxxxzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 353);

    auto g_0_xxxxzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 354);

    auto g_0_xxxxzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 355);

    auto g_0_xxxxzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 356);

    auto g_0_xxxxzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 357);

    auto g_0_xxxxzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 358);

    auto g_0_xxxxzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 359);

    #pragma omp simd aligned(g_0_xxxxz_0_xxxxxxx_0, g_0_xxxxz_0_xxxxxxx_1, g_0_xxxxz_0_xxxxxxy_0, g_0_xxxxz_0_xxxxxxy_1, g_0_xxxxz_0_xxxxxyy_0, g_0_xxxxz_0_xxxxxyy_1, g_0_xxxxz_0_xxxxyyy_0, g_0_xxxxz_0_xxxxyyy_1, g_0_xxxxz_0_xxxyyyy_0, g_0_xxxxz_0_xxxyyyy_1, g_0_xxxxz_0_xxyyyyy_0, g_0_xxxxz_0_xxyyyyy_1, g_0_xxxxz_0_xyyyyyy_0, g_0_xxxxz_0_xyyyyyy_1, g_0_xxxxzz_0_xxxxxxx_0, g_0_xxxxzz_0_xxxxxxx_1, g_0_xxxxzz_0_xxxxxxy_0, g_0_xxxxzz_0_xxxxxxy_1, g_0_xxxxzz_0_xxxxxyy_0, g_0_xxxxzz_0_xxxxxyy_1, g_0_xxxxzz_0_xxxxyyy_0, g_0_xxxxzz_0_xxxxyyy_1, g_0_xxxxzz_0_xxxyyyy_0, g_0_xxxxzz_0_xxxyyyy_1, g_0_xxxxzz_0_xxyyyyy_0, g_0_xxxxzz_0_xxyyyyy_1, g_0_xxxxzz_0_xyyyyyy_0, g_0_xxxxzz_0_xyyyyyy_1, g_0_xxxxzzz_0_xxxxxxx_0, g_0_xxxxzzz_0_xxxxxxy_0, g_0_xxxxzzz_0_xxxxxxz_0, g_0_xxxxzzz_0_xxxxxyy_0, g_0_xxxxzzz_0_xxxxxyz_0, g_0_xxxxzzz_0_xxxxxzz_0, g_0_xxxxzzz_0_xxxxyyy_0, g_0_xxxxzzz_0_xxxxyyz_0, g_0_xxxxzzz_0_xxxxyzz_0, g_0_xxxxzzz_0_xxxxzzz_0, g_0_xxxxzzz_0_xxxyyyy_0, g_0_xxxxzzz_0_xxxyyyz_0, g_0_xxxxzzz_0_xxxyyzz_0, g_0_xxxxzzz_0_xxxyzzz_0, g_0_xxxxzzz_0_xxxzzzz_0, g_0_xxxxzzz_0_xxyyyyy_0, g_0_xxxxzzz_0_xxyyyyz_0, g_0_xxxxzzz_0_xxyyyzz_0, g_0_xxxxzzz_0_xxyyzzz_0, g_0_xxxxzzz_0_xxyzzzz_0, g_0_xxxxzzz_0_xxzzzzz_0, g_0_xxxxzzz_0_xyyyyyy_0, g_0_xxxxzzz_0_xyyyyyz_0, g_0_xxxxzzz_0_xyyyyzz_0, g_0_xxxxzzz_0_xyyyzzz_0, g_0_xxxxzzz_0_xyyzzzz_0, g_0_xxxxzzz_0_xyzzzzz_0, g_0_xxxxzzz_0_xzzzzzz_0, g_0_xxxxzzz_0_yyyyyyy_0, g_0_xxxxzzz_0_yyyyyyz_0, g_0_xxxxzzz_0_yyyyyzz_0, g_0_xxxxzzz_0_yyyyzzz_0, g_0_xxxxzzz_0_yyyzzzz_0, g_0_xxxxzzz_0_yyzzzzz_0, g_0_xxxxzzz_0_yzzzzzz_0, g_0_xxxxzzz_0_zzzzzzz_0, g_0_xxxzzz_0_xxxxxxz_0, g_0_xxxzzz_0_xxxxxxz_1, g_0_xxxzzz_0_xxxxxyz_0, g_0_xxxzzz_0_xxxxxyz_1, g_0_xxxzzz_0_xxxxxz_1, g_0_xxxzzz_0_xxxxxzz_0, g_0_xxxzzz_0_xxxxxzz_1, g_0_xxxzzz_0_xxxxyyz_0, g_0_xxxzzz_0_xxxxyyz_1, g_0_xxxzzz_0_xxxxyz_1, g_0_xxxzzz_0_xxxxyzz_0, g_0_xxxzzz_0_xxxxyzz_1, g_0_xxxzzz_0_xxxxzz_1, g_0_xxxzzz_0_xxxxzzz_0, g_0_xxxzzz_0_xxxxzzz_1, g_0_xxxzzz_0_xxxyyyz_0, g_0_xxxzzz_0_xxxyyyz_1, g_0_xxxzzz_0_xxxyyz_1, g_0_xxxzzz_0_xxxyyzz_0, g_0_xxxzzz_0_xxxyyzz_1, g_0_xxxzzz_0_xxxyzz_1, g_0_xxxzzz_0_xxxyzzz_0, g_0_xxxzzz_0_xxxyzzz_1, g_0_xxxzzz_0_xxxzzz_1, g_0_xxxzzz_0_xxxzzzz_0, g_0_xxxzzz_0_xxxzzzz_1, g_0_xxxzzz_0_xxyyyyz_0, g_0_xxxzzz_0_xxyyyyz_1, g_0_xxxzzz_0_xxyyyz_1, g_0_xxxzzz_0_xxyyyzz_0, g_0_xxxzzz_0_xxyyyzz_1, g_0_xxxzzz_0_xxyyzz_1, g_0_xxxzzz_0_xxyyzzz_0, g_0_xxxzzz_0_xxyyzzz_1, g_0_xxxzzz_0_xxyzzz_1, g_0_xxxzzz_0_xxyzzzz_0, g_0_xxxzzz_0_xxyzzzz_1, g_0_xxxzzz_0_xxzzzz_1, g_0_xxxzzz_0_xxzzzzz_0, g_0_xxxzzz_0_xxzzzzz_1, g_0_xxxzzz_0_xyyyyyz_0, g_0_xxxzzz_0_xyyyyyz_1, g_0_xxxzzz_0_xyyyyz_1, g_0_xxxzzz_0_xyyyyzz_0, g_0_xxxzzz_0_xyyyyzz_1, g_0_xxxzzz_0_xyyyzz_1, g_0_xxxzzz_0_xyyyzzz_0, g_0_xxxzzz_0_xyyyzzz_1, g_0_xxxzzz_0_xyyzzz_1, g_0_xxxzzz_0_xyyzzzz_0, g_0_xxxzzz_0_xyyzzzz_1, g_0_xxxzzz_0_xyzzzz_1, g_0_xxxzzz_0_xyzzzzz_0, g_0_xxxzzz_0_xyzzzzz_1, g_0_xxxzzz_0_xzzzzz_1, g_0_xxxzzz_0_xzzzzzz_0, g_0_xxxzzz_0_xzzzzzz_1, g_0_xxxzzz_0_yyyyyyy_0, g_0_xxxzzz_0_yyyyyyy_1, g_0_xxxzzz_0_yyyyyyz_0, g_0_xxxzzz_0_yyyyyyz_1, g_0_xxxzzz_0_yyyyyz_1, g_0_xxxzzz_0_yyyyyzz_0, g_0_xxxzzz_0_yyyyyzz_1, g_0_xxxzzz_0_yyyyzz_1, g_0_xxxzzz_0_yyyyzzz_0, g_0_xxxzzz_0_yyyyzzz_1, g_0_xxxzzz_0_yyyzzz_1, g_0_xxxzzz_0_yyyzzzz_0, g_0_xxxzzz_0_yyyzzzz_1, g_0_xxxzzz_0_yyzzzz_1, g_0_xxxzzz_0_yyzzzzz_0, g_0_xxxzzz_0_yyzzzzz_1, g_0_xxxzzz_0_yzzzzz_1, g_0_xxxzzz_0_yzzzzzz_0, g_0_xxxzzz_0_yzzzzzz_1, g_0_xxxzzz_0_zzzzzz_1, g_0_xxxzzz_0_zzzzzzz_0, g_0_xxxzzz_0_zzzzzzz_1, g_0_xxzzz_0_xxxxxxz_0, g_0_xxzzz_0_xxxxxxz_1, g_0_xxzzz_0_xxxxxyz_0, g_0_xxzzz_0_xxxxxyz_1, g_0_xxzzz_0_xxxxxzz_0, g_0_xxzzz_0_xxxxxzz_1, g_0_xxzzz_0_xxxxyyz_0, g_0_xxzzz_0_xxxxyyz_1, g_0_xxzzz_0_xxxxyzz_0, g_0_xxzzz_0_xxxxyzz_1, g_0_xxzzz_0_xxxxzzz_0, g_0_xxzzz_0_xxxxzzz_1, g_0_xxzzz_0_xxxyyyz_0, g_0_xxzzz_0_xxxyyyz_1, g_0_xxzzz_0_xxxyyzz_0, g_0_xxzzz_0_xxxyyzz_1, g_0_xxzzz_0_xxxyzzz_0, g_0_xxzzz_0_xxxyzzz_1, g_0_xxzzz_0_xxxzzzz_0, g_0_xxzzz_0_xxxzzzz_1, g_0_xxzzz_0_xxyyyyz_0, g_0_xxzzz_0_xxyyyyz_1, g_0_xxzzz_0_xxyyyzz_0, g_0_xxzzz_0_xxyyyzz_1, g_0_xxzzz_0_xxyyzzz_0, g_0_xxzzz_0_xxyyzzz_1, g_0_xxzzz_0_xxyzzzz_0, g_0_xxzzz_0_xxyzzzz_1, g_0_xxzzz_0_xxzzzzz_0, g_0_xxzzz_0_xxzzzzz_1, g_0_xxzzz_0_xyyyyyz_0, g_0_xxzzz_0_xyyyyyz_1, g_0_xxzzz_0_xyyyyzz_0, g_0_xxzzz_0_xyyyyzz_1, g_0_xxzzz_0_xyyyzzz_0, g_0_xxzzz_0_xyyyzzz_1, g_0_xxzzz_0_xyyzzzz_0, g_0_xxzzz_0_xyyzzzz_1, g_0_xxzzz_0_xyzzzzz_0, g_0_xxzzz_0_xyzzzzz_1, g_0_xxzzz_0_xzzzzzz_0, g_0_xxzzz_0_xzzzzzz_1, g_0_xxzzz_0_yyyyyyy_0, g_0_xxzzz_0_yyyyyyy_1, g_0_xxzzz_0_yyyyyyz_0, g_0_xxzzz_0_yyyyyyz_1, g_0_xxzzz_0_yyyyyzz_0, g_0_xxzzz_0_yyyyyzz_1, g_0_xxzzz_0_yyyyzzz_0, g_0_xxzzz_0_yyyyzzz_1, g_0_xxzzz_0_yyyzzzz_0, g_0_xxzzz_0_yyyzzzz_1, g_0_xxzzz_0_yyzzzzz_0, g_0_xxzzz_0_yyzzzzz_1, g_0_xxzzz_0_yzzzzzz_0, g_0_xxzzz_0_yzzzzzz_1, g_0_xxzzz_0_zzzzzzz_0, g_0_xxzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzz_0_xxxxxxx_0[i] = 2.0 * g_0_xxxxz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxxxxx_0[i] * pb_z + g_0_xxxxzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxxxy_0[i] = 2.0 * g_0_xxxxz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxxxxy_0[i] * pb_z + g_0_xxxxzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxxxz_0[i] = 3.0 * g_0_xxzzz_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxxxz_0[i] * pb_x + g_0_xxxzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxxz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxxxyy_0[i] * pb_z + g_0_xxxxzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxxyz_0[i] = 3.0 * g_0_xxzzz_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxxyz_0[i] * pb_x + g_0_xxxzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxxxzz_0[i] = 3.0 * g_0_xxzzz_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxxzz_0[i] * pb_x + g_0_xxxzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxxyyy_0[i] = 2.0 * g_0_xxxxz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxxyyy_0[i] * pb_z + g_0_xxxxzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxyyz_0[i] = 3.0 * g_0_xxzzz_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxyyz_0[i] * pb_x + g_0_xxxzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxxyzz_0[i] = 3.0 * g_0_xxzzz_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxyzz_0[i] * pb_x + g_0_xxxzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxxzzz_0[i] = 3.0 * g_0_xxzzz_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxzzz_0[i] * pb_x + g_0_xxxzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxyyyy_0[i] = 2.0 * g_0_xxxxz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxxyyyy_0[i] * pb_z + g_0_xxxxzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxzzz_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyyyz_0[i] * pb_x + g_0_xxxzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxyyzz_0[i] = 3.0 * g_0_xxzzz_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyyzz_0[i] * pb_x + g_0_xxxzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxyzzz_0[i] = 3.0 * g_0_xxzzz_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyzzz_0[i] * pb_x + g_0_xxxzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxzzzz_0[i] = 3.0 * g_0_xxzzz_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxzzzz_0[i] * pb_x + g_0_xxxzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyyyyy_0[i] = 2.0 * g_0_xxxxz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xxyyyyy_0[i] * pb_z + g_0_xxxxzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxyyyyz_0[i] = 3.0 * g_0_xxzzz_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyyyz_0[i] * pb_x + g_0_xxxzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxzzz_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyyzz_0[i] * pb_x + g_0_xxxzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyyzzz_0[i] = 3.0 * g_0_xxzzz_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyzzz_0[i] * pb_x + g_0_xxxzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyzzzz_0[i] = 3.0 * g_0_xxzzz_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyzzzz_0[i] * pb_x + g_0_xxxzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxzzzzz_0[i] = 3.0 * g_0_xxzzz_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxzzzzz_0[i] * pb_x + g_0_xxxzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyyyyy_0[i] = 2.0 * g_0_xxxxz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxzz_0_xyyyyyy_0[i] * pb_z + g_0_xxxxzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xyyyyyz_0[i] = 3.0 * g_0_xxzzz_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyyyz_0[i] * pb_x + g_0_xxxzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyyyzz_0[i] = 3.0 * g_0_xxzzz_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyyzz_0[i] * pb_x + g_0_xxxzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxzzz_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyzzz_0[i] * pb_x + g_0_xxxzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyzzzz_0[i] = 3.0 * g_0_xxzzz_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyzzzz_0[i] * pb_x + g_0_xxxzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyzzzzz_0[i] = 3.0 * g_0_xxzzz_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyzzzzz_0[i] * pb_x + g_0_xxxzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xzzzzzz_0[i] = 3.0 * g_0_xxzzz_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xzzzzzz_0[i] * pb_x + g_0_xxxzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyyyy_0[i] = 3.0 * g_0_xxzzz_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyyyyy_0[i] * pb_x + g_0_xxxzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyyyz_0[i] = 3.0 * g_0_xxzzz_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyyyyz_0[i] * pb_x + g_0_xxxzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyyzz_0[i] = 3.0 * g_0_xxzzz_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyyyzz_0[i] * pb_x + g_0_xxxzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyzzz_0[i] = 3.0 * g_0_xxzzz_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyyzzz_0[i] * pb_x + g_0_xxxzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxzzz_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyyzzzz_0[i] * pb_x + g_0_xxxzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyzzzzz_0[i] = 3.0 * g_0_xxzzz_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yyzzzzz_0[i] * pb_x + g_0_xxxzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yzzzzzz_0[i] = 3.0 * g_0_xxzzz_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_yzzzzzz_0[i] * pb_x + g_0_xxxzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_zzzzzzz_0[i] = 3.0 * g_0_xxzzz_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxzzz_0_zzzzzzz_0[i] * pb_x + g_0_xxxzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 360-396 components of targeted buffer : SKSK

    auto g_0_xxxyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 360);

    auto g_0_xxxyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 361);

    auto g_0_xxxyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 362);

    auto g_0_xxxyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 363);

    auto g_0_xxxyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 364);

    auto g_0_xxxyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 365);

    auto g_0_xxxyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 366);

    auto g_0_xxxyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 367);

    auto g_0_xxxyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 368);

    auto g_0_xxxyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 369);

    auto g_0_xxxyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 370);

    auto g_0_xxxyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 371);

    auto g_0_xxxyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 372);

    auto g_0_xxxyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 373);

    auto g_0_xxxyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 374);

    auto g_0_xxxyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 375);

    auto g_0_xxxyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 376);

    auto g_0_xxxyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 377);

    auto g_0_xxxyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 378);

    auto g_0_xxxyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 379);

    auto g_0_xxxyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 380);

    auto g_0_xxxyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 381);

    auto g_0_xxxyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 382);

    auto g_0_xxxyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 383);

    auto g_0_xxxyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 384);

    auto g_0_xxxyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 385);

    auto g_0_xxxyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 386);

    auto g_0_xxxyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 387);

    auto g_0_xxxyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 388);

    auto g_0_xxxyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 389);

    auto g_0_xxxyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 390);

    auto g_0_xxxyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 391);

    auto g_0_xxxyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 392);

    auto g_0_xxxyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 393);

    auto g_0_xxxyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 394);

    auto g_0_xxxyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 395);

    #pragma omp simd aligned(g_0_xxxyy_0_xxxxxxx_0, g_0_xxxyy_0_xxxxxxx_1, g_0_xxxyy_0_xxxxxxz_0, g_0_xxxyy_0_xxxxxxz_1, g_0_xxxyy_0_xxxxxzz_0, g_0_xxxyy_0_xxxxxzz_1, g_0_xxxyy_0_xxxxzzz_0, g_0_xxxyy_0_xxxxzzz_1, g_0_xxxyy_0_xxxzzzz_0, g_0_xxxyy_0_xxxzzzz_1, g_0_xxxyy_0_xxzzzzz_0, g_0_xxxyy_0_xxzzzzz_1, g_0_xxxyy_0_xzzzzzz_0, g_0_xxxyy_0_xzzzzzz_1, g_0_xxxyyy_0_xxxxxxx_0, g_0_xxxyyy_0_xxxxxxx_1, g_0_xxxyyy_0_xxxxxxz_0, g_0_xxxyyy_0_xxxxxxz_1, g_0_xxxyyy_0_xxxxxzz_0, g_0_xxxyyy_0_xxxxxzz_1, g_0_xxxyyy_0_xxxxzzz_0, g_0_xxxyyy_0_xxxxzzz_1, g_0_xxxyyy_0_xxxzzzz_0, g_0_xxxyyy_0_xxxzzzz_1, g_0_xxxyyy_0_xxzzzzz_0, g_0_xxxyyy_0_xxzzzzz_1, g_0_xxxyyy_0_xzzzzzz_0, g_0_xxxyyy_0_xzzzzzz_1, g_0_xxxyyyy_0_xxxxxxx_0, g_0_xxxyyyy_0_xxxxxxy_0, g_0_xxxyyyy_0_xxxxxxz_0, g_0_xxxyyyy_0_xxxxxyy_0, g_0_xxxyyyy_0_xxxxxyz_0, g_0_xxxyyyy_0_xxxxxzz_0, g_0_xxxyyyy_0_xxxxyyy_0, g_0_xxxyyyy_0_xxxxyyz_0, g_0_xxxyyyy_0_xxxxyzz_0, g_0_xxxyyyy_0_xxxxzzz_0, g_0_xxxyyyy_0_xxxyyyy_0, g_0_xxxyyyy_0_xxxyyyz_0, g_0_xxxyyyy_0_xxxyyzz_0, g_0_xxxyyyy_0_xxxyzzz_0, g_0_xxxyyyy_0_xxxzzzz_0, g_0_xxxyyyy_0_xxyyyyy_0, g_0_xxxyyyy_0_xxyyyyz_0, g_0_xxxyyyy_0_xxyyyzz_0, g_0_xxxyyyy_0_xxyyzzz_0, g_0_xxxyyyy_0_xxyzzzz_0, g_0_xxxyyyy_0_xxzzzzz_0, g_0_xxxyyyy_0_xyyyyyy_0, g_0_xxxyyyy_0_xyyyyyz_0, g_0_xxxyyyy_0_xyyyyzz_0, g_0_xxxyyyy_0_xyyyzzz_0, g_0_xxxyyyy_0_xyyzzzz_0, g_0_xxxyyyy_0_xyzzzzz_0, g_0_xxxyyyy_0_xzzzzzz_0, g_0_xxxyyyy_0_yyyyyyy_0, g_0_xxxyyyy_0_yyyyyyz_0, g_0_xxxyyyy_0_yyyyyzz_0, g_0_xxxyyyy_0_yyyyzzz_0, g_0_xxxyyyy_0_yyyzzzz_0, g_0_xxxyyyy_0_yyzzzzz_0, g_0_xxxyyyy_0_yzzzzzz_0, g_0_xxxyyyy_0_zzzzzzz_0, g_0_xxyyyy_0_xxxxxxy_0, g_0_xxyyyy_0_xxxxxxy_1, g_0_xxyyyy_0_xxxxxy_1, g_0_xxyyyy_0_xxxxxyy_0, g_0_xxyyyy_0_xxxxxyy_1, g_0_xxyyyy_0_xxxxxyz_0, g_0_xxyyyy_0_xxxxxyz_1, g_0_xxyyyy_0_xxxxyy_1, g_0_xxyyyy_0_xxxxyyy_0, g_0_xxyyyy_0_xxxxyyy_1, g_0_xxyyyy_0_xxxxyyz_0, g_0_xxyyyy_0_xxxxyyz_1, g_0_xxyyyy_0_xxxxyz_1, g_0_xxyyyy_0_xxxxyzz_0, g_0_xxyyyy_0_xxxxyzz_1, g_0_xxyyyy_0_xxxyyy_1, g_0_xxyyyy_0_xxxyyyy_0, g_0_xxyyyy_0_xxxyyyy_1, g_0_xxyyyy_0_xxxyyyz_0, g_0_xxyyyy_0_xxxyyyz_1, g_0_xxyyyy_0_xxxyyz_1, g_0_xxyyyy_0_xxxyyzz_0, g_0_xxyyyy_0_xxxyyzz_1, g_0_xxyyyy_0_xxxyzz_1, g_0_xxyyyy_0_xxxyzzz_0, g_0_xxyyyy_0_xxxyzzz_1, g_0_xxyyyy_0_xxyyyy_1, g_0_xxyyyy_0_xxyyyyy_0, g_0_xxyyyy_0_xxyyyyy_1, g_0_xxyyyy_0_xxyyyyz_0, g_0_xxyyyy_0_xxyyyyz_1, g_0_xxyyyy_0_xxyyyz_1, g_0_xxyyyy_0_xxyyyzz_0, g_0_xxyyyy_0_xxyyyzz_1, g_0_xxyyyy_0_xxyyzz_1, g_0_xxyyyy_0_xxyyzzz_0, g_0_xxyyyy_0_xxyyzzz_1, g_0_xxyyyy_0_xxyzzz_1, g_0_xxyyyy_0_xxyzzzz_0, g_0_xxyyyy_0_xxyzzzz_1, g_0_xxyyyy_0_xyyyyy_1, g_0_xxyyyy_0_xyyyyyy_0, g_0_xxyyyy_0_xyyyyyy_1, g_0_xxyyyy_0_xyyyyyz_0, g_0_xxyyyy_0_xyyyyyz_1, g_0_xxyyyy_0_xyyyyz_1, g_0_xxyyyy_0_xyyyyzz_0, g_0_xxyyyy_0_xyyyyzz_1, g_0_xxyyyy_0_xyyyzz_1, g_0_xxyyyy_0_xyyyzzz_0, g_0_xxyyyy_0_xyyyzzz_1, g_0_xxyyyy_0_xyyzzz_1, g_0_xxyyyy_0_xyyzzzz_0, g_0_xxyyyy_0_xyyzzzz_1, g_0_xxyyyy_0_xyzzzz_1, g_0_xxyyyy_0_xyzzzzz_0, g_0_xxyyyy_0_xyzzzzz_1, g_0_xxyyyy_0_yyyyyy_1, g_0_xxyyyy_0_yyyyyyy_0, g_0_xxyyyy_0_yyyyyyy_1, g_0_xxyyyy_0_yyyyyyz_0, g_0_xxyyyy_0_yyyyyyz_1, g_0_xxyyyy_0_yyyyyz_1, g_0_xxyyyy_0_yyyyyzz_0, g_0_xxyyyy_0_yyyyyzz_1, g_0_xxyyyy_0_yyyyzz_1, g_0_xxyyyy_0_yyyyzzz_0, g_0_xxyyyy_0_yyyyzzz_1, g_0_xxyyyy_0_yyyzzz_1, g_0_xxyyyy_0_yyyzzzz_0, g_0_xxyyyy_0_yyyzzzz_1, g_0_xxyyyy_0_yyzzzz_1, g_0_xxyyyy_0_yyzzzzz_0, g_0_xxyyyy_0_yyzzzzz_1, g_0_xxyyyy_0_yzzzzz_1, g_0_xxyyyy_0_yzzzzzz_0, g_0_xxyyyy_0_yzzzzzz_1, g_0_xxyyyy_0_zzzzzzz_0, g_0_xxyyyy_0_zzzzzzz_1, g_0_xyyyy_0_xxxxxxy_0, g_0_xyyyy_0_xxxxxxy_1, g_0_xyyyy_0_xxxxxyy_0, g_0_xyyyy_0_xxxxxyy_1, g_0_xyyyy_0_xxxxxyz_0, g_0_xyyyy_0_xxxxxyz_1, g_0_xyyyy_0_xxxxyyy_0, g_0_xyyyy_0_xxxxyyy_1, g_0_xyyyy_0_xxxxyyz_0, g_0_xyyyy_0_xxxxyyz_1, g_0_xyyyy_0_xxxxyzz_0, g_0_xyyyy_0_xxxxyzz_1, g_0_xyyyy_0_xxxyyyy_0, g_0_xyyyy_0_xxxyyyy_1, g_0_xyyyy_0_xxxyyyz_0, g_0_xyyyy_0_xxxyyyz_1, g_0_xyyyy_0_xxxyyzz_0, g_0_xyyyy_0_xxxyyzz_1, g_0_xyyyy_0_xxxyzzz_0, g_0_xyyyy_0_xxxyzzz_1, g_0_xyyyy_0_xxyyyyy_0, g_0_xyyyy_0_xxyyyyy_1, g_0_xyyyy_0_xxyyyyz_0, g_0_xyyyy_0_xxyyyyz_1, g_0_xyyyy_0_xxyyyzz_0, g_0_xyyyy_0_xxyyyzz_1, g_0_xyyyy_0_xxyyzzz_0, g_0_xyyyy_0_xxyyzzz_1, g_0_xyyyy_0_xxyzzzz_0, g_0_xyyyy_0_xxyzzzz_1, g_0_xyyyy_0_xyyyyyy_0, g_0_xyyyy_0_xyyyyyy_1, g_0_xyyyy_0_xyyyyyz_0, g_0_xyyyy_0_xyyyyyz_1, g_0_xyyyy_0_xyyyyzz_0, g_0_xyyyy_0_xyyyyzz_1, g_0_xyyyy_0_xyyyzzz_0, g_0_xyyyy_0_xyyyzzz_1, g_0_xyyyy_0_xyyzzzz_0, g_0_xyyyy_0_xyyzzzz_1, g_0_xyyyy_0_xyzzzzz_0, g_0_xyyyy_0_xyzzzzz_1, g_0_xyyyy_0_yyyyyyy_0, g_0_xyyyy_0_yyyyyyy_1, g_0_xyyyy_0_yyyyyyz_0, g_0_xyyyy_0_yyyyyyz_1, g_0_xyyyy_0_yyyyyzz_0, g_0_xyyyy_0_yyyyyzz_1, g_0_xyyyy_0_yyyyzzz_0, g_0_xyyyy_0_yyyyzzz_1, g_0_xyyyy_0_yyyzzzz_0, g_0_xyyyy_0_yyyzzzz_1, g_0_xyyyy_0_yyzzzzz_0, g_0_xyyyy_0_yyzzzzz_1, g_0_xyyyy_0_yzzzzzz_0, g_0_xyyyy_0_yzzzzzz_1, g_0_xyyyy_0_zzzzzzz_0, g_0_xyyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyy_0_xxxxxxx_0[i] = 3.0 * g_0_xxxyy_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxxxxx_0[i] * pb_y + g_0_xxxyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxxxxy_0[i] = 2.0 * g_0_xyyyy_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxxxy_0[i] * pb_x + g_0_xxyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxxxz_0[i] = 3.0 * g_0_xxxyy_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxxxxz_0[i] * pb_y + g_0_xxxyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxxxyy_0[i] = 2.0 * g_0_xyyyy_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxxyy_0[i] * pb_x + g_0_xxyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxxyz_0[i] = 2.0 * g_0_xyyyy_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxxyz_0[i] * pb_x + g_0_xxyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxxzz_0[i] = 3.0 * g_0_xxxyy_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxxxzz_0[i] * pb_y + g_0_xxxyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxxyyy_0[i] = 2.0 * g_0_xyyyy_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxyyy_0[i] * pb_x + g_0_xxyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxyyz_0[i] = 2.0 * g_0_xyyyy_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxyyz_0[i] * pb_x + g_0_xxyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxyzz_0[i] = 2.0 * g_0_xyyyy_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxyzz_0[i] * pb_x + g_0_xxyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxzzz_0[i] = 3.0 * g_0_xxxyy_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxxzzz_0[i] * pb_y + g_0_xxxyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxyyyy_0[i] = 2.0 * g_0_xyyyy_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyyyy_0[i] * pb_x + g_0_xxyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxyyyz_0[i] = 2.0 * g_0_xyyyy_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyyyz_0[i] * pb_x + g_0_xxyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxyyzz_0[i] = 2.0 * g_0_xyyyy_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyyzz_0[i] * pb_x + g_0_xxyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxyzzz_0[i] = 2.0 * g_0_xyyyy_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyzzz_0[i] * pb_x + g_0_xxyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxzzzz_0[i] = 3.0 * g_0_xxxyy_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxxzzzz_0[i] * pb_y + g_0_xxxyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxyyyyy_0[i] = 2.0 * g_0_xyyyy_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyyyy_0[i] * pb_x + g_0_xxyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyyyyz_0[i] = 2.0 * g_0_xyyyy_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyyyz_0[i] * pb_x + g_0_xxyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyyyzz_0[i] = 2.0 * g_0_xyyyy_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyyzz_0[i] * pb_x + g_0_xxyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyyzzz_0[i] = 2.0 * g_0_xyyyy_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyzzz_0[i] * pb_x + g_0_xxyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyzzzz_0[i] = 2.0 * g_0_xyyyy_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyzzzz_0[i] * pb_x + g_0_xxyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxzzzzz_0[i] = 3.0 * g_0_xxxyy_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xxzzzzz_0[i] * pb_y + g_0_xxxyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xyyyyyy_0[i] = 2.0 * g_0_xyyyy_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyyyy_0[i] * pb_x + g_0_xxyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyyyyz_0[i] = 2.0 * g_0_xyyyy_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyyyz_0[i] * pb_x + g_0_xxyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyyyzz_0[i] = 2.0 * g_0_xyyyy_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyyzz_0[i] * pb_x + g_0_xxyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyyzzz_0[i] = 2.0 * g_0_xyyyy_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyzzz_0[i] * pb_x + g_0_xxyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyzzzz_0[i] = 2.0 * g_0_xyyyy_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyzzzz_0[i] * pb_x + g_0_xxyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyzzzzz_0[i] = 2.0 * g_0_xyyyy_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyzzzzz_0[i] * pb_x + g_0_xxyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xzzzzzz_0[i] = 3.0 * g_0_xxxyy_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyy_0_xzzzzzz_0[i] * pb_y + g_0_xxxyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_yyyyyyy_0[i] = 2.0 * g_0_xyyyy_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyyyy_0[i] * pb_x + g_0_xxyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyyyyz_0[i] = 2.0 * g_0_xyyyy_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyyyz_0[i] * pb_x + g_0_xxyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyyyzz_0[i] = 2.0 * g_0_xyyyy_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyyzz_0[i] * pb_x + g_0_xxyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyyzzz_0[i] = 2.0 * g_0_xyyyy_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyyzzz_0[i] * pb_x + g_0_xxyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyzzzz_0[i] = 2.0 * g_0_xyyyy_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyyzzzz_0[i] * pb_x + g_0_xxyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyzzzzz_0[i] = 2.0 * g_0_xyyyy_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yyzzzzz_0[i] * pb_x + g_0_xxyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yzzzzzz_0[i] = 2.0 * g_0_xyyyy_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_yzzzzzz_0[i] * pb_x + g_0_xxyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_zzzzzzz_0[i] = 2.0 * g_0_xyyyy_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_zzzzzzz_0[i] * pb_x + g_0_xxyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 396-432 components of targeted buffer : SKSK

    auto g_0_xxxyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 396);

    auto g_0_xxxyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 397);

    auto g_0_xxxyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 398);

    auto g_0_xxxyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 399);

    auto g_0_xxxyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 400);

    auto g_0_xxxyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 401);

    auto g_0_xxxyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 402);

    auto g_0_xxxyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 403);

    auto g_0_xxxyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 404);

    auto g_0_xxxyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 405);

    auto g_0_xxxyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 406);

    auto g_0_xxxyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 407);

    auto g_0_xxxyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 408);

    auto g_0_xxxyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 409);

    auto g_0_xxxyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 410);

    auto g_0_xxxyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 411);

    auto g_0_xxxyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 412);

    auto g_0_xxxyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 413);

    auto g_0_xxxyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 414);

    auto g_0_xxxyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 415);

    auto g_0_xxxyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 416);

    auto g_0_xxxyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 417);

    auto g_0_xxxyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 418);

    auto g_0_xxxyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 419);

    auto g_0_xxxyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 420);

    auto g_0_xxxyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 421);

    auto g_0_xxxyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 422);

    auto g_0_xxxyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 423);

    auto g_0_xxxyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 424);

    auto g_0_xxxyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 425);

    auto g_0_xxxyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 426);

    auto g_0_xxxyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 427);

    auto g_0_xxxyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 428);

    auto g_0_xxxyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 429);

    auto g_0_xxxyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 430);

    auto g_0_xxxyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 431);

    #pragma omp simd aligned(g_0_xxxyyy_0_xxxxxx_1, g_0_xxxyyy_0_xxxxxxx_0, g_0_xxxyyy_0_xxxxxxx_1, g_0_xxxyyy_0_xxxxxxy_0, g_0_xxxyyy_0_xxxxxxy_1, g_0_xxxyyy_0_xxxxxxz_0, g_0_xxxyyy_0_xxxxxxz_1, g_0_xxxyyy_0_xxxxxy_1, g_0_xxxyyy_0_xxxxxyy_0, g_0_xxxyyy_0_xxxxxyy_1, g_0_xxxyyy_0_xxxxxyz_0, g_0_xxxyyy_0_xxxxxyz_1, g_0_xxxyyy_0_xxxxxz_1, g_0_xxxyyy_0_xxxxxzz_0, g_0_xxxyyy_0_xxxxxzz_1, g_0_xxxyyy_0_xxxxyy_1, g_0_xxxyyy_0_xxxxyyy_0, g_0_xxxyyy_0_xxxxyyy_1, g_0_xxxyyy_0_xxxxyyz_0, g_0_xxxyyy_0_xxxxyyz_1, g_0_xxxyyy_0_xxxxyz_1, g_0_xxxyyy_0_xxxxyzz_0, g_0_xxxyyy_0_xxxxyzz_1, g_0_xxxyyy_0_xxxxzz_1, g_0_xxxyyy_0_xxxxzzz_0, g_0_xxxyyy_0_xxxxzzz_1, g_0_xxxyyy_0_xxxyyy_1, g_0_xxxyyy_0_xxxyyyy_0, g_0_xxxyyy_0_xxxyyyy_1, g_0_xxxyyy_0_xxxyyyz_0, g_0_xxxyyy_0_xxxyyyz_1, g_0_xxxyyy_0_xxxyyz_1, g_0_xxxyyy_0_xxxyyzz_0, g_0_xxxyyy_0_xxxyyzz_1, g_0_xxxyyy_0_xxxyzz_1, g_0_xxxyyy_0_xxxyzzz_0, g_0_xxxyyy_0_xxxyzzz_1, g_0_xxxyyy_0_xxxzzz_1, g_0_xxxyyy_0_xxxzzzz_0, g_0_xxxyyy_0_xxxzzzz_1, g_0_xxxyyy_0_xxyyyy_1, g_0_xxxyyy_0_xxyyyyy_0, g_0_xxxyyy_0_xxyyyyy_1, g_0_xxxyyy_0_xxyyyyz_0, g_0_xxxyyy_0_xxyyyyz_1, g_0_xxxyyy_0_xxyyyz_1, g_0_xxxyyy_0_xxyyyzz_0, g_0_xxxyyy_0_xxyyyzz_1, g_0_xxxyyy_0_xxyyzz_1, g_0_xxxyyy_0_xxyyzzz_0, g_0_xxxyyy_0_xxyyzzz_1, g_0_xxxyyy_0_xxyzzz_1, g_0_xxxyyy_0_xxyzzzz_0, g_0_xxxyyy_0_xxyzzzz_1, g_0_xxxyyy_0_xxzzzz_1, g_0_xxxyyy_0_xxzzzzz_0, g_0_xxxyyy_0_xxzzzzz_1, g_0_xxxyyy_0_xyyyyy_1, g_0_xxxyyy_0_xyyyyyy_0, g_0_xxxyyy_0_xyyyyyy_1, g_0_xxxyyy_0_xyyyyyz_0, g_0_xxxyyy_0_xyyyyyz_1, g_0_xxxyyy_0_xyyyyz_1, g_0_xxxyyy_0_xyyyyzz_0, g_0_xxxyyy_0_xyyyyzz_1, g_0_xxxyyy_0_xyyyzz_1, g_0_xxxyyy_0_xyyyzzz_0, g_0_xxxyyy_0_xyyyzzz_1, g_0_xxxyyy_0_xyyzzz_1, g_0_xxxyyy_0_xyyzzzz_0, g_0_xxxyyy_0_xyyzzzz_1, g_0_xxxyyy_0_xyzzzz_1, g_0_xxxyyy_0_xyzzzzz_0, g_0_xxxyyy_0_xyzzzzz_1, g_0_xxxyyy_0_xzzzzz_1, g_0_xxxyyy_0_xzzzzzz_0, g_0_xxxyyy_0_xzzzzzz_1, g_0_xxxyyy_0_yyyyyy_1, g_0_xxxyyy_0_yyyyyyy_0, g_0_xxxyyy_0_yyyyyyy_1, g_0_xxxyyy_0_yyyyyyz_0, g_0_xxxyyy_0_yyyyyyz_1, g_0_xxxyyy_0_yyyyyz_1, g_0_xxxyyy_0_yyyyyzz_0, g_0_xxxyyy_0_yyyyyzz_1, g_0_xxxyyy_0_yyyyzz_1, g_0_xxxyyy_0_yyyyzzz_0, g_0_xxxyyy_0_yyyyzzz_1, g_0_xxxyyy_0_yyyzzz_1, g_0_xxxyyy_0_yyyzzzz_0, g_0_xxxyyy_0_yyyzzzz_1, g_0_xxxyyy_0_yyzzzz_1, g_0_xxxyyy_0_yyzzzzz_0, g_0_xxxyyy_0_yyzzzzz_1, g_0_xxxyyy_0_yzzzzz_1, g_0_xxxyyy_0_yzzzzzz_0, g_0_xxxyyy_0_yzzzzzz_1, g_0_xxxyyy_0_zzzzzz_1, g_0_xxxyyy_0_zzzzzzz_0, g_0_xxxyyy_0_zzzzzzz_1, g_0_xxxyyyz_0_xxxxxxx_0, g_0_xxxyyyz_0_xxxxxxy_0, g_0_xxxyyyz_0_xxxxxxz_0, g_0_xxxyyyz_0_xxxxxyy_0, g_0_xxxyyyz_0_xxxxxyz_0, g_0_xxxyyyz_0_xxxxxzz_0, g_0_xxxyyyz_0_xxxxyyy_0, g_0_xxxyyyz_0_xxxxyyz_0, g_0_xxxyyyz_0_xxxxyzz_0, g_0_xxxyyyz_0_xxxxzzz_0, g_0_xxxyyyz_0_xxxyyyy_0, g_0_xxxyyyz_0_xxxyyyz_0, g_0_xxxyyyz_0_xxxyyzz_0, g_0_xxxyyyz_0_xxxyzzz_0, g_0_xxxyyyz_0_xxxzzzz_0, g_0_xxxyyyz_0_xxyyyyy_0, g_0_xxxyyyz_0_xxyyyyz_0, g_0_xxxyyyz_0_xxyyyzz_0, g_0_xxxyyyz_0_xxyyzzz_0, g_0_xxxyyyz_0_xxyzzzz_0, g_0_xxxyyyz_0_xxzzzzz_0, g_0_xxxyyyz_0_xyyyyyy_0, g_0_xxxyyyz_0_xyyyyyz_0, g_0_xxxyyyz_0_xyyyyzz_0, g_0_xxxyyyz_0_xyyyzzz_0, g_0_xxxyyyz_0_xyyzzzz_0, g_0_xxxyyyz_0_xyzzzzz_0, g_0_xxxyyyz_0_xzzzzzz_0, g_0_xxxyyyz_0_yyyyyyy_0, g_0_xxxyyyz_0_yyyyyyz_0, g_0_xxxyyyz_0_yyyyyzz_0, g_0_xxxyyyz_0_yyyyzzz_0, g_0_xxxyyyz_0_yyyzzzz_0, g_0_xxxyyyz_0_yyzzzzz_0, g_0_xxxyyyz_0_yzzzzzz_0, g_0_xxxyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyz_0_xxxxxxx_0[i] = g_0_xxxyyy_0_xxxxxxx_0[i] * pb_z + g_0_xxxyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxxxy_0[i] = g_0_xxxyyy_0_xxxxxxy_0[i] * pb_z + g_0_xxxyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxxxz_0[i] = g_0_xxxyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxxxz_0[i] * pb_z + g_0_xxxyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxxyy_0[i] = g_0_xxxyyy_0_xxxxxyy_0[i] * pb_z + g_0_xxxyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxxyz_0[i] = g_0_xxxyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxxyz_0[i] * pb_z + g_0_xxxyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxxzz_0[i] * pb_z + g_0_xxxyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxyyy_0[i] = g_0_xxxyyy_0_xxxxyyy_0[i] * pb_z + g_0_xxxyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxyyz_0[i] = g_0_xxxyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxyyz_0[i] * pb_z + g_0_xxxyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxyzz_0[i] * pb_z + g_0_xxxyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxzzz_0[i] * pb_z + g_0_xxxyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyyyy_0[i] = g_0_xxxyyy_0_xxxyyyy_0[i] * pb_z + g_0_xxxyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyyyz_0[i] = g_0_xxxyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyyyz_0[i] * pb_z + g_0_xxxyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyyzz_0[i] * pb_z + g_0_xxxyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyzzz_0[i] * pb_z + g_0_xxxyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxzzzz_0[i] * pb_z + g_0_xxxyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyyyy_0[i] = g_0_xxxyyy_0_xxyyyyy_0[i] * pb_z + g_0_xxxyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyyyz_0[i] = g_0_xxxyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyyyz_0[i] * pb_z + g_0_xxxyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyyzz_0[i] * pb_z + g_0_xxxyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyzzz_0[i] * pb_z + g_0_xxxyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyzzzz_0[i] * pb_z + g_0_xxxyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxzzzzz_0[i] * pb_z + g_0_xxxyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyyyy_0[i] = g_0_xxxyyy_0_xyyyyyy_0[i] * pb_z + g_0_xxxyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyyyz_0[i] = g_0_xxxyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyyyz_0[i] * pb_z + g_0_xxxyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyyzz_0[i] * pb_z + g_0_xxxyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyzzz_0[i] * pb_z + g_0_xxxyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyzzzz_0[i] * pb_z + g_0_xxxyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyzzzzz_0[i] * pb_z + g_0_xxxyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xzzzzzz_0[i] * pb_z + g_0_xxxyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyyyy_0[i] = g_0_xxxyyy_0_yyyyyyy_0[i] * pb_z + g_0_xxxyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyyyz_0[i] = g_0_xxxyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyyyyz_0[i] * pb_z + g_0_xxxyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyyyzz_0[i] * pb_z + g_0_xxxyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyyzzz_0[i] * pb_z + g_0_xxxyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyzzzz_0[i] * pb_z + g_0_xxxyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyzzzzz_0[i] * pb_z + g_0_xxxyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yzzzzzz_0[i] * pb_z + g_0_xxxyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_zzzzzzz_0[i] * pb_z + g_0_xxxyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 432-468 components of targeted buffer : SKSK

    auto g_0_xxxyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 432);

    auto g_0_xxxyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 433);

    auto g_0_xxxyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 434);

    auto g_0_xxxyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 435);

    auto g_0_xxxyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 436);

    auto g_0_xxxyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 437);

    auto g_0_xxxyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 438);

    auto g_0_xxxyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 439);

    auto g_0_xxxyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 440);

    auto g_0_xxxyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 441);

    auto g_0_xxxyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 442);

    auto g_0_xxxyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 443);

    auto g_0_xxxyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 444);

    auto g_0_xxxyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 445);

    auto g_0_xxxyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 446);

    auto g_0_xxxyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 447);

    auto g_0_xxxyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 448);

    auto g_0_xxxyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 449);

    auto g_0_xxxyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 450);

    auto g_0_xxxyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 451);

    auto g_0_xxxyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 452);

    auto g_0_xxxyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 453);

    auto g_0_xxxyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 454);

    auto g_0_xxxyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 455);

    auto g_0_xxxyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 456);

    auto g_0_xxxyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 457);

    auto g_0_xxxyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 458);

    auto g_0_xxxyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 459);

    auto g_0_xxxyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 460);

    auto g_0_xxxyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 461);

    auto g_0_xxxyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 462);

    auto g_0_xxxyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 463);

    auto g_0_xxxyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 464);

    auto g_0_xxxyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 465);

    auto g_0_xxxyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 466);

    auto g_0_xxxyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 467);

    #pragma omp simd aligned(g_0_xxxyy_0_xxxxxxy_0, g_0_xxxyy_0_xxxxxxy_1, g_0_xxxyy_0_xxxxxyy_0, g_0_xxxyy_0_xxxxxyy_1, g_0_xxxyy_0_xxxxyyy_0, g_0_xxxyy_0_xxxxyyy_1, g_0_xxxyy_0_xxxyyyy_0, g_0_xxxyy_0_xxxyyyy_1, g_0_xxxyy_0_xxyyyyy_0, g_0_xxxyy_0_xxyyyyy_1, g_0_xxxyy_0_xyyyyyy_0, g_0_xxxyy_0_xyyyyyy_1, g_0_xxxyyz_0_xxxxxxy_0, g_0_xxxyyz_0_xxxxxxy_1, g_0_xxxyyz_0_xxxxxyy_0, g_0_xxxyyz_0_xxxxxyy_1, g_0_xxxyyz_0_xxxxyyy_0, g_0_xxxyyz_0_xxxxyyy_1, g_0_xxxyyz_0_xxxyyyy_0, g_0_xxxyyz_0_xxxyyyy_1, g_0_xxxyyz_0_xxyyyyy_0, g_0_xxxyyz_0_xxyyyyy_1, g_0_xxxyyz_0_xyyyyyy_0, g_0_xxxyyz_0_xyyyyyy_1, g_0_xxxyyzz_0_xxxxxxx_0, g_0_xxxyyzz_0_xxxxxxy_0, g_0_xxxyyzz_0_xxxxxxz_0, g_0_xxxyyzz_0_xxxxxyy_0, g_0_xxxyyzz_0_xxxxxyz_0, g_0_xxxyyzz_0_xxxxxzz_0, g_0_xxxyyzz_0_xxxxyyy_0, g_0_xxxyyzz_0_xxxxyyz_0, g_0_xxxyyzz_0_xxxxyzz_0, g_0_xxxyyzz_0_xxxxzzz_0, g_0_xxxyyzz_0_xxxyyyy_0, g_0_xxxyyzz_0_xxxyyyz_0, g_0_xxxyyzz_0_xxxyyzz_0, g_0_xxxyyzz_0_xxxyzzz_0, g_0_xxxyyzz_0_xxxzzzz_0, g_0_xxxyyzz_0_xxyyyyy_0, g_0_xxxyyzz_0_xxyyyyz_0, g_0_xxxyyzz_0_xxyyyzz_0, g_0_xxxyyzz_0_xxyyzzz_0, g_0_xxxyyzz_0_xxyzzzz_0, g_0_xxxyyzz_0_xxzzzzz_0, g_0_xxxyyzz_0_xyyyyyy_0, g_0_xxxyyzz_0_xyyyyyz_0, g_0_xxxyyzz_0_xyyyyzz_0, g_0_xxxyyzz_0_xyyyzzz_0, g_0_xxxyyzz_0_xyyzzzz_0, g_0_xxxyyzz_0_xyzzzzz_0, g_0_xxxyyzz_0_xzzzzzz_0, g_0_xxxyyzz_0_yyyyyyy_0, g_0_xxxyyzz_0_yyyyyyz_0, g_0_xxxyyzz_0_yyyyyzz_0, g_0_xxxyyzz_0_yyyyzzz_0, g_0_xxxyyzz_0_yyyzzzz_0, g_0_xxxyyzz_0_yyzzzzz_0, g_0_xxxyyzz_0_yzzzzzz_0, g_0_xxxyyzz_0_zzzzzzz_0, g_0_xxxyzz_0_xxxxxxx_0, g_0_xxxyzz_0_xxxxxxx_1, g_0_xxxyzz_0_xxxxxxz_0, g_0_xxxyzz_0_xxxxxxz_1, g_0_xxxyzz_0_xxxxxzz_0, g_0_xxxyzz_0_xxxxxzz_1, g_0_xxxyzz_0_xxxxzzz_0, g_0_xxxyzz_0_xxxxzzz_1, g_0_xxxyzz_0_xxxzzzz_0, g_0_xxxyzz_0_xxxzzzz_1, g_0_xxxyzz_0_xxzzzzz_0, g_0_xxxyzz_0_xxzzzzz_1, g_0_xxxyzz_0_xzzzzzz_0, g_0_xxxyzz_0_xzzzzzz_1, g_0_xxxzz_0_xxxxxxx_0, g_0_xxxzz_0_xxxxxxx_1, g_0_xxxzz_0_xxxxxxz_0, g_0_xxxzz_0_xxxxxxz_1, g_0_xxxzz_0_xxxxxzz_0, g_0_xxxzz_0_xxxxxzz_1, g_0_xxxzz_0_xxxxzzz_0, g_0_xxxzz_0_xxxxzzz_1, g_0_xxxzz_0_xxxzzzz_0, g_0_xxxzz_0_xxxzzzz_1, g_0_xxxzz_0_xxzzzzz_0, g_0_xxxzz_0_xxzzzzz_1, g_0_xxxzz_0_xzzzzzz_0, g_0_xxxzz_0_xzzzzzz_1, g_0_xxyyzz_0_xxxxxyz_0, g_0_xxyyzz_0_xxxxxyz_1, g_0_xxyyzz_0_xxxxyyz_0, g_0_xxyyzz_0_xxxxyyz_1, g_0_xxyyzz_0_xxxxyz_1, g_0_xxyyzz_0_xxxxyzz_0, g_0_xxyyzz_0_xxxxyzz_1, g_0_xxyyzz_0_xxxyyyz_0, g_0_xxyyzz_0_xxxyyyz_1, g_0_xxyyzz_0_xxxyyz_1, g_0_xxyyzz_0_xxxyyzz_0, g_0_xxyyzz_0_xxxyyzz_1, g_0_xxyyzz_0_xxxyzz_1, g_0_xxyyzz_0_xxxyzzz_0, g_0_xxyyzz_0_xxxyzzz_1, g_0_xxyyzz_0_xxyyyyz_0, g_0_xxyyzz_0_xxyyyyz_1, g_0_xxyyzz_0_xxyyyz_1, g_0_xxyyzz_0_xxyyyzz_0, g_0_xxyyzz_0_xxyyyzz_1, g_0_xxyyzz_0_xxyyzz_1, g_0_xxyyzz_0_xxyyzzz_0, g_0_xxyyzz_0_xxyyzzz_1, g_0_xxyyzz_0_xxyzzz_1, g_0_xxyyzz_0_xxyzzzz_0, g_0_xxyyzz_0_xxyzzzz_1, g_0_xxyyzz_0_xyyyyyz_0, g_0_xxyyzz_0_xyyyyyz_1, g_0_xxyyzz_0_xyyyyz_1, g_0_xxyyzz_0_xyyyyzz_0, g_0_xxyyzz_0_xyyyyzz_1, g_0_xxyyzz_0_xyyyzz_1, g_0_xxyyzz_0_xyyyzzz_0, g_0_xxyyzz_0_xyyyzzz_1, g_0_xxyyzz_0_xyyzzz_1, g_0_xxyyzz_0_xyyzzzz_0, g_0_xxyyzz_0_xyyzzzz_1, g_0_xxyyzz_0_xyzzzz_1, g_0_xxyyzz_0_xyzzzzz_0, g_0_xxyyzz_0_xyzzzzz_1, g_0_xxyyzz_0_yyyyyyy_0, g_0_xxyyzz_0_yyyyyyy_1, g_0_xxyyzz_0_yyyyyyz_0, g_0_xxyyzz_0_yyyyyyz_1, g_0_xxyyzz_0_yyyyyz_1, g_0_xxyyzz_0_yyyyyzz_0, g_0_xxyyzz_0_yyyyyzz_1, g_0_xxyyzz_0_yyyyzz_1, g_0_xxyyzz_0_yyyyzzz_0, g_0_xxyyzz_0_yyyyzzz_1, g_0_xxyyzz_0_yyyzzz_1, g_0_xxyyzz_0_yyyzzzz_0, g_0_xxyyzz_0_yyyzzzz_1, g_0_xxyyzz_0_yyzzzz_1, g_0_xxyyzz_0_yyzzzzz_0, g_0_xxyyzz_0_yyzzzzz_1, g_0_xxyyzz_0_yzzzzz_1, g_0_xxyyzz_0_yzzzzzz_0, g_0_xxyyzz_0_yzzzzzz_1, g_0_xxyyzz_0_zzzzzzz_0, g_0_xxyyzz_0_zzzzzzz_1, g_0_xyyzz_0_xxxxxyz_0, g_0_xyyzz_0_xxxxxyz_1, g_0_xyyzz_0_xxxxyyz_0, g_0_xyyzz_0_xxxxyyz_1, g_0_xyyzz_0_xxxxyzz_0, g_0_xyyzz_0_xxxxyzz_1, g_0_xyyzz_0_xxxyyyz_0, g_0_xyyzz_0_xxxyyyz_1, g_0_xyyzz_0_xxxyyzz_0, g_0_xyyzz_0_xxxyyzz_1, g_0_xyyzz_0_xxxyzzz_0, g_0_xyyzz_0_xxxyzzz_1, g_0_xyyzz_0_xxyyyyz_0, g_0_xyyzz_0_xxyyyyz_1, g_0_xyyzz_0_xxyyyzz_0, g_0_xyyzz_0_xxyyyzz_1, g_0_xyyzz_0_xxyyzzz_0, g_0_xyyzz_0_xxyyzzz_1, g_0_xyyzz_0_xxyzzzz_0, g_0_xyyzz_0_xxyzzzz_1, g_0_xyyzz_0_xyyyyyz_0, g_0_xyyzz_0_xyyyyyz_1, g_0_xyyzz_0_xyyyyzz_0, g_0_xyyzz_0_xyyyyzz_1, g_0_xyyzz_0_xyyyzzz_0, g_0_xyyzz_0_xyyyzzz_1, g_0_xyyzz_0_xyyzzzz_0, g_0_xyyzz_0_xyyzzzz_1, g_0_xyyzz_0_xyzzzzz_0, g_0_xyyzz_0_xyzzzzz_1, g_0_xyyzz_0_yyyyyyy_0, g_0_xyyzz_0_yyyyyyy_1, g_0_xyyzz_0_yyyyyyz_0, g_0_xyyzz_0_yyyyyyz_1, g_0_xyyzz_0_yyyyyzz_0, g_0_xyyzz_0_yyyyyzz_1, g_0_xyyzz_0_yyyyzzz_0, g_0_xyyzz_0_yyyyzzz_1, g_0_xyyzz_0_yyyzzzz_0, g_0_xyyzz_0_yyyzzzz_1, g_0_xyyzz_0_yyzzzzz_0, g_0_xyyzz_0_yyzzzzz_1, g_0_xyyzz_0_yzzzzzz_0, g_0_xyyzz_0_yzzzzzz_1, g_0_xyyzz_0_zzzzzzz_0, g_0_xyyzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzz_0_xxxxxxx_0[i] = g_0_xxxzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxxxxy_0[i] = g_0_xxxyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxxxxy_0[i] * pb_z + g_0_xxxyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxxxxz_0[i] = g_0_xxxzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxxxyy_0[i] = g_0_xxxyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxxxyy_0[i] * pb_z + g_0_xxxyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxxxyz_0[i] = 2.0 * g_0_xyyzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxxxyz_0[i] * pb_x + g_0_xxyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxxxzz_0[i] = g_0_xxxzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxxyyy_0[i] = g_0_xxxyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxxyyy_0[i] * pb_z + g_0_xxxyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxxyyz_0[i] = 2.0 * g_0_xyyzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxxyyz_0[i] * pb_x + g_0_xxyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxxyzz_0[i] = 2.0 * g_0_xyyzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxxyzz_0[i] * pb_x + g_0_xxyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxxzzz_0[i] = g_0_xxxzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxyyyy_0[i] = g_0_xxxyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxyyyy_0[i] * pb_z + g_0_xxxyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxyyyz_0[i] = 2.0 * g_0_xyyzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxyyyz_0[i] * pb_x + g_0_xxyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxyyzz_0[i] = 2.0 * g_0_xyyzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxyyzz_0[i] * pb_x + g_0_xxyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxyzzz_0[i] = 2.0 * g_0_xyyzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxyzzz_0[i] * pb_x + g_0_xxyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxzzzz_0[i] = g_0_xxxzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxyyyyy_0[i] = g_0_xxxyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxyyyyy_0[i] * pb_z + g_0_xxxyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxyyyyz_0[i] = 2.0 * g_0_xyyzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyyyyz_0[i] * pb_x + g_0_xxyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxyyyzz_0[i] = 2.0 * g_0_xyyzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyyyzz_0[i] * pb_x + g_0_xxyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxyyzzz_0[i] = 2.0 * g_0_xyyzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyyzzz_0[i] * pb_x + g_0_xxyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxyzzzz_0[i] = 2.0 * g_0_xyyzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyzzzz_0[i] * pb_x + g_0_xxyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxzzzzz_0[i] = g_0_xxxzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xyyyyyy_0[i] = g_0_xxxyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xyyyyyy_0[i] * pb_z + g_0_xxxyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xyyyyyz_0[i] = 2.0 * g_0_xyyzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyyyyz_0[i] * pb_x + g_0_xxyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyyyyzz_0[i] = 2.0 * g_0_xyyzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyyyzz_0[i] * pb_x + g_0_xxyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyyyzzz_0[i] = 2.0 * g_0_xyyzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyyzzz_0[i] * pb_x + g_0_xxyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyyzzzz_0[i] = 2.0 * g_0_xyyzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyzzzz_0[i] * pb_x + g_0_xxyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyzzzzz_0[i] = 2.0 * g_0_xyyzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyzzzzz_0[i] * pb_x + g_0_xxyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xzzzzzz_0[i] = g_0_xxxzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_yyyyyyy_0[i] = 2.0 * g_0_xyyzz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyyyyy_0[i] * pb_x + g_0_xxyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyyyyz_0[i] = 2.0 * g_0_xyyzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyyyyz_0[i] * pb_x + g_0_xxyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyyyzz_0[i] = 2.0 * g_0_xyyzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyyyzz_0[i] * pb_x + g_0_xxyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyyzzz_0[i] = 2.0 * g_0_xyyzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyyzzz_0[i] * pb_x + g_0_xxyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyzzzz_0[i] = 2.0 * g_0_xyyzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyyzzzz_0[i] * pb_x + g_0_xxyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyzzzzz_0[i] = 2.0 * g_0_xyyzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yyzzzzz_0[i] * pb_x + g_0_xxyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yzzzzzz_0[i] = 2.0 * g_0_xyyzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_yzzzzzz_0[i] * pb_x + g_0_xxyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_zzzzzzz_0[i] = 2.0 * g_0_xyyzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_zzzzzzz_0[i] * pb_x + g_0_xxyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 468-504 components of targeted buffer : SKSK

    auto g_0_xxxyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 468);

    auto g_0_xxxyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 469);

    auto g_0_xxxyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 470);

    auto g_0_xxxyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 471);

    auto g_0_xxxyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 472);

    auto g_0_xxxyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 473);

    auto g_0_xxxyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 474);

    auto g_0_xxxyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 475);

    auto g_0_xxxyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 476);

    auto g_0_xxxyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 477);

    auto g_0_xxxyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 478);

    auto g_0_xxxyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 479);

    auto g_0_xxxyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 480);

    auto g_0_xxxyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 481);

    auto g_0_xxxyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 482);

    auto g_0_xxxyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 483);

    auto g_0_xxxyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 484);

    auto g_0_xxxyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 485);

    auto g_0_xxxyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 486);

    auto g_0_xxxyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 487);

    auto g_0_xxxyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 488);

    auto g_0_xxxyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 489);

    auto g_0_xxxyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 490);

    auto g_0_xxxyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 491);

    auto g_0_xxxyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 492);

    auto g_0_xxxyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 493);

    auto g_0_xxxyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 494);

    auto g_0_xxxyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 495);

    auto g_0_xxxyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 496);

    auto g_0_xxxyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 497);

    auto g_0_xxxyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 498);

    auto g_0_xxxyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 499);

    auto g_0_xxxyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 500);

    auto g_0_xxxyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 501);

    auto g_0_xxxyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 502);

    auto g_0_xxxyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 503);

    #pragma omp simd aligned(g_0_xxxyzzz_0_xxxxxxx_0, g_0_xxxyzzz_0_xxxxxxy_0, g_0_xxxyzzz_0_xxxxxxz_0, g_0_xxxyzzz_0_xxxxxyy_0, g_0_xxxyzzz_0_xxxxxyz_0, g_0_xxxyzzz_0_xxxxxzz_0, g_0_xxxyzzz_0_xxxxyyy_0, g_0_xxxyzzz_0_xxxxyyz_0, g_0_xxxyzzz_0_xxxxyzz_0, g_0_xxxyzzz_0_xxxxzzz_0, g_0_xxxyzzz_0_xxxyyyy_0, g_0_xxxyzzz_0_xxxyyyz_0, g_0_xxxyzzz_0_xxxyyzz_0, g_0_xxxyzzz_0_xxxyzzz_0, g_0_xxxyzzz_0_xxxzzzz_0, g_0_xxxyzzz_0_xxyyyyy_0, g_0_xxxyzzz_0_xxyyyyz_0, g_0_xxxyzzz_0_xxyyyzz_0, g_0_xxxyzzz_0_xxyyzzz_0, g_0_xxxyzzz_0_xxyzzzz_0, g_0_xxxyzzz_0_xxzzzzz_0, g_0_xxxyzzz_0_xyyyyyy_0, g_0_xxxyzzz_0_xyyyyyz_0, g_0_xxxyzzz_0_xyyyyzz_0, g_0_xxxyzzz_0_xyyyzzz_0, g_0_xxxyzzz_0_xyyzzzz_0, g_0_xxxyzzz_0_xyzzzzz_0, g_0_xxxyzzz_0_xzzzzzz_0, g_0_xxxyzzz_0_yyyyyyy_0, g_0_xxxyzzz_0_yyyyyyz_0, g_0_xxxyzzz_0_yyyyyzz_0, g_0_xxxyzzz_0_yyyyzzz_0, g_0_xxxyzzz_0_yyyzzzz_0, g_0_xxxyzzz_0_yyzzzzz_0, g_0_xxxyzzz_0_yzzzzzz_0, g_0_xxxyzzz_0_zzzzzzz_0, g_0_xxxzzz_0_xxxxxx_1, g_0_xxxzzz_0_xxxxxxx_0, g_0_xxxzzz_0_xxxxxxx_1, g_0_xxxzzz_0_xxxxxxy_0, g_0_xxxzzz_0_xxxxxxy_1, g_0_xxxzzz_0_xxxxxxz_0, g_0_xxxzzz_0_xxxxxxz_1, g_0_xxxzzz_0_xxxxxy_1, g_0_xxxzzz_0_xxxxxyy_0, g_0_xxxzzz_0_xxxxxyy_1, g_0_xxxzzz_0_xxxxxyz_0, g_0_xxxzzz_0_xxxxxyz_1, g_0_xxxzzz_0_xxxxxz_1, g_0_xxxzzz_0_xxxxxzz_0, g_0_xxxzzz_0_xxxxxzz_1, g_0_xxxzzz_0_xxxxyy_1, g_0_xxxzzz_0_xxxxyyy_0, g_0_xxxzzz_0_xxxxyyy_1, g_0_xxxzzz_0_xxxxyyz_0, g_0_xxxzzz_0_xxxxyyz_1, g_0_xxxzzz_0_xxxxyz_1, g_0_xxxzzz_0_xxxxyzz_0, g_0_xxxzzz_0_xxxxyzz_1, g_0_xxxzzz_0_xxxxzz_1, g_0_xxxzzz_0_xxxxzzz_0, g_0_xxxzzz_0_xxxxzzz_1, g_0_xxxzzz_0_xxxyyy_1, g_0_xxxzzz_0_xxxyyyy_0, g_0_xxxzzz_0_xxxyyyy_1, g_0_xxxzzz_0_xxxyyyz_0, g_0_xxxzzz_0_xxxyyyz_1, g_0_xxxzzz_0_xxxyyz_1, g_0_xxxzzz_0_xxxyyzz_0, g_0_xxxzzz_0_xxxyyzz_1, g_0_xxxzzz_0_xxxyzz_1, g_0_xxxzzz_0_xxxyzzz_0, g_0_xxxzzz_0_xxxyzzz_1, g_0_xxxzzz_0_xxxzzz_1, g_0_xxxzzz_0_xxxzzzz_0, g_0_xxxzzz_0_xxxzzzz_1, g_0_xxxzzz_0_xxyyyy_1, g_0_xxxzzz_0_xxyyyyy_0, g_0_xxxzzz_0_xxyyyyy_1, g_0_xxxzzz_0_xxyyyyz_0, g_0_xxxzzz_0_xxyyyyz_1, g_0_xxxzzz_0_xxyyyz_1, g_0_xxxzzz_0_xxyyyzz_0, g_0_xxxzzz_0_xxyyyzz_1, g_0_xxxzzz_0_xxyyzz_1, g_0_xxxzzz_0_xxyyzzz_0, g_0_xxxzzz_0_xxyyzzz_1, g_0_xxxzzz_0_xxyzzz_1, g_0_xxxzzz_0_xxyzzzz_0, g_0_xxxzzz_0_xxyzzzz_1, g_0_xxxzzz_0_xxzzzz_1, g_0_xxxzzz_0_xxzzzzz_0, g_0_xxxzzz_0_xxzzzzz_1, g_0_xxxzzz_0_xyyyyy_1, g_0_xxxzzz_0_xyyyyyy_0, g_0_xxxzzz_0_xyyyyyy_1, g_0_xxxzzz_0_xyyyyyz_0, g_0_xxxzzz_0_xyyyyyz_1, g_0_xxxzzz_0_xyyyyz_1, g_0_xxxzzz_0_xyyyyzz_0, g_0_xxxzzz_0_xyyyyzz_1, g_0_xxxzzz_0_xyyyzz_1, g_0_xxxzzz_0_xyyyzzz_0, g_0_xxxzzz_0_xyyyzzz_1, g_0_xxxzzz_0_xyyzzz_1, g_0_xxxzzz_0_xyyzzzz_0, g_0_xxxzzz_0_xyyzzzz_1, g_0_xxxzzz_0_xyzzzz_1, g_0_xxxzzz_0_xyzzzzz_0, g_0_xxxzzz_0_xyzzzzz_1, g_0_xxxzzz_0_xzzzzz_1, g_0_xxxzzz_0_xzzzzzz_0, g_0_xxxzzz_0_xzzzzzz_1, g_0_xxxzzz_0_yyyyyy_1, g_0_xxxzzz_0_yyyyyyy_0, g_0_xxxzzz_0_yyyyyyy_1, g_0_xxxzzz_0_yyyyyyz_0, g_0_xxxzzz_0_yyyyyyz_1, g_0_xxxzzz_0_yyyyyz_1, g_0_xxxzzz_0_yyyyyzz_0, g_0_xxxzzz_0_yyyyyzz_1, g_0_xxxzzz_0_yyyyzz_1, g_0_xxxzzz_0_yyyyzzz_0, g_0_xxxzzz_0_yyyyzzz_1, g_0_xxxzzz_0_yyyzzz_1, g_0_xxxzzz_0_yyyzzzz_0, g_0_xxxzzz_0_yyyzzzz_1, g_0_xxxzzz_0_yyzzzz_1, g_0_xxxzzz_0_yyzzzzz_0, g_0_xxxzzz_0_yyzzzzz_1, g_0_xxxzzz_0_yzzzzz_1, g_0_xxxzzz_0_yzzzzzz_0, g_0_xxxzzz_0_yzzzzzz_1, g_0_xxxzzz_0_zzzzzz_1, g_0_xxxzzz_0_zzzzzzz_0, g_0_xxxzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzz_0_xxxxxxx_0[i] = g_0_xxxzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxxxy_0[i] = g_0_xxxzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxxxy_0[i] * pb_y + g_0_xxxzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxxxz_0[i] = g_0_xxxzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxxyy_0[i] * pb_y + g_0_xxxzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxxyz_0[i] = g_0_xxxzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxxyz_0[i] * pb_y + g_0_xxxzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxxzz_0[i] = g_0_xxxzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxxzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxyyy_0[i] * pb_y + g_0_xxxzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxyyz_0[i] * pb_y + g_0_xxxzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxyzz_0[i] = g_0_xxxzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxyzz_0[i] * pb_y + g_0_xxxzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxzzz_0[i] = g_0_xxxzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyyyy_0[i] * pb_y + g_0_xxxzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyyyz_0[i] * pb_y + g_0_xxxzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyyzz_0[i] * pb_y + g_0_xxxzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyzzz_0[i] = g_0_xxxzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyzzz_0[i] * pb_y + g_0_xxxzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxzzzz_0[i] = g_0_xxxzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxxzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyyyy_0[i] * pb_y + g_0_xxxzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyyyz_0[i] * pb_y + g_0_xxxzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyyzz_0[i] * pb_y + g_0_xxxzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyzzz_0[i] * pb_y + g_0_xxxzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyzzzz_0[i] = g_0_xxxzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyzzzz_0[i] * pb_y + g_0_xxxzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxzzzzz_0[i] = g_0_xxxzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxxzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyyyy_0[i] * pb_y + g_0_xxxzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyyyz_0[i] * pb_y + g_0_xxxzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyyzz_0[i] * pb_y + g_0_xxxzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyzzz_0[i] * pb_y + g_0_xxxzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyzzzz_0[i] * pb_y + g_0_xxxzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyzzzzz_0[i] = g_0_xxxzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyzzzzz_0[i] * pb_y + g_0_xxxzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xzzzzzz_0[i] = g_0_xxxzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxxzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyyyy_0[i] * pb_y + g_0_xxxzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyyyz_0[i] * pb_y + g_0_xxxzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyyzz_0[i] * pb_y + g_0_xxxzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyzzz_0[i] * pb_y + g_0_xxxzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyzzzz_0[i] * pb_y + g_0_xxxzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyzzzzz_0[i] * pb_y + g_0_xxxzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yzzzzzz_0[i] = g_0_xxxzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yzzzzzz_0[i] * pb_y + g_0_xxxzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_zzzzzzz_0[i] = g_0_xxxzzz_0_zzzzzzz_0[i] * pb_y + g_0_xxxzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 504-540 components of targeted buffer : SKSK

    auto g_0_xxxzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 504);

    auto g_0_xxxzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 505);

    auto g_0_xxxzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 506);

    auto g_0_xxxzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 507);

    auto g_0_xxxzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 508);

    auto g_0_xxxzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 509);

    auto g_0_xxxzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 510);

    auto g_0_xxxzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 511);

    auto g_0_xxxzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 512);

    auto g_0_xxxzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 513);

    auto g_0_xxxzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 514);

    auto g_0_xxxzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 515);

    auto g_0_xxxzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 516);

    auto g_0_xxxzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 517);

    auto g_0_xxxzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 518);

    auto g_0_xxxzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 519);

    auto g_0_xxxzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 520);

    auto g_0_xxxzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 521);

    auto g_0_xxxzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 522);

    auto g_0_xxxzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 523);

    auto g_0_xxxzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 524);

    auto g_0_xxxzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 525);

    auto g_0_xxxzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 526);

    auto g_0_xxxzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 527);

    auto g_0_xxxzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 528);

    auto g_0_xxxzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 529);

    auto g_0_xxxzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 530);

    auto g_0_xxxzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 531);

    auto g_0_xxxzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 532);

    auto g_0_xxxzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 533);

    auto g_0_xxxzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 534);

    auto g_0_xxxzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 535);

    auto g_0_xxxzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 536);

    auto g_0_xxxzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 537);

    auto g_0_xxxzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 538);

    auto g_0_xxxzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 539);

    #pragma omp simd aligned(g_0_xxxzz_0_xxxxxxx_0, g_0_xxxzz_0_xxxxxxx_1, g_0_xxxzz_0_xxxxxxy_0, g_0_xxxzz_0_xxxxxxy_1, g_0_xxxzz_0_xxxxxyy_0, g_0_xxxzz_0_xxxxxyy_1, g_0_xxxzz_0_xxxxyyy_0, g_0_xxxzz_0_xxxxyyy_1, g_0_xxxzz_0_xxxyyyy_0, g_0_xxxzz_0_xxxyyyy_1, g_0_xxxzz_0_xxyyyyy_0, g_0_xxxzz_0_xxyyyyy_1, g_0_xxxzz_0_xyyyyyy_0, g_0_xxxzz_0_xyyyyyy_1, g_0_xxxzzz_0_xxxxxxx_0, g_0_xxxzzz_0_xxxxxxx_1, g_0_xxxzzz_0_xxxxxxy_0, g_0_xxxzzz_0_xxxxxxy_1, g_0_xxxzzz_0_xxxxxyy_0, g_0_xxxzzz_0_xxxxxyy_1, g_0_xxxzzz_0_xxxxyyy_0, g_0_xxxzzz_0_xxxxyyy_1, g_0_xxxzzz_0_xxxyyyy_0, g_0_xxxzzz_0_xxxyyyy_1, g_0_xxxzzz_0_xxyyyyy_0, g_0_xxxzzz_0_xxyyyyy_1, g_0_xxxzzz_0_xyyyyyy_0, g_0_xxxzzz_0_xyyyyyy_1, g_0_xxxzzzz_0_xxxxxxx_0, g_0_xxxzzzz_0_xxxxxxy_0, g_0_xxxzzzz_0_xxxxxxz_0, g_0_xxxzzzz_0_xxxxxyy_0, g_0_xxxzzzz_0_xxxxxyz_0, g_0_xxxzzzz_0_xxxxxzz_0, g_0_xxxzzzz_0_xxxxyyy_0, g_0_xxxzzzz_0_xxxxyyz_0, g_0_xxxzzzz_0_xxxxyzz_0, g_0_xxxzzzz_0_xxxxzzz_0, g_0_xxxzzzz_0_xxxyyyy_0, g_0_xxxzzzz_0_xxxyyyz_0, g_0_xxxzzzz_0_xxxyyzz_0, g_0_xxxzzzz_0_xxxyzzz_0, g_0_xxxzzzz_0_xxxzzzz_0, g_0_xxxzzzz_0_xxyyyyy_0, g_0_xxxzzzz_0_xxyyyyz_0, g_0_xxxzzzz_0_xxyyyzz_0, g_0_xxxzzzz_0_xxyyzzz_0, g_0_xxxzzzz_0_xxyzzzz_0, g_0_xxxzzzz_0_xxzzzzz_0, g_0_xxxzzzz_0_xyyyyyy_0, g_0_xxxzzzz_0_xyyyyyz_0, g_0_xxxzzzz_0_xyyyyzz_0, g_0_xxxzzzz_0_xyyyzzz_0, g_0_xxxzzzz_0_xyyzzzz_0, g_0_xxxzzzz_0_xyzzzzz_0, g_0_xxxzzzz_0_xzzzzzz_0, g_0_xxxzzzz_0_yyyyyyy_0, g_0_xxxzzzz_0_yyyyyyz_0, g_0_xxxzzzz_0_yyyyyzz_0, g_0_xxxzzzz_0_yyyyzzz_0, g_0_xxxzzzz_0_yyyzzzz_0, g_0_xxxzzzz_0_yyzzzzz_0, g_0_xxxzzzz_0_yzzzzzz_0, g_0_xxxzzzz_0_zzzzzzz_0, g_0_xxzzzz_0_xxxxxxz_0, g_0_xxzzzz_0_xxxxxxz_1, g_0_xxzzzz_0_xxxxxyz_0, g_0_xxzzzz_0_xxxxxyz_1, g_0_xxzzzz_0_xxxxxz_1, g_0_xxzzzz_0_xxxxxzz_0, g_0_xxzzzz_0_xxxxxzz_1, g_0_xxzzzz_0_xxxxyyz_0, g_0_xxzzzz_0_xxxxyyz_1, g_0_xxzzzz_0_xxxxyz_1, g_0_xxzzzz_0_xxxxyzz_0, g_0_xxzzzz_0_xxxxyzz_1, g_0_xxzzzz_0_xxxxzz_1, g_0_xxzzzz_0_xxxxzzz_0, g_0_xxzzzz_0_xxxxzzz_1, g_0_xxzzzz_0_xxxyyyz_0, g_0_xxzzzz_0_xxxyyyz_1, g_0_xxzzzz_0_xxxyyz_1, g_0_xxzzzz_0_xxxyyzz_0, g_0_xxzzzz_0_xxxyyzz_1, g_0_xxzzzz_0_xxxyzz_1, g_0_xxzzzz_0_xxxyzzz_0, g_0_xxzzzz_0_xxxyzzz_1, g_0_xxzzzz_0_xxxzzz_1, g_0_xxzzzz_0_xxxzzzz_0, g_0_xxzzzz_0_xxxzzzz_1, g_0_xxzzzz_0_xxyyyyz_0, g_0_xxzzzz_0_xxyyyyz_1, g_0_xxzzzz_0_xxyyyz_1, g_0_xxzzzz_0_xxyyyzz_0, g_0_xxzzzz_0_xxyyyzz_1, g_0_xxzzzz_0_xxyyzz_1, g_0_xxzzzz_0_xxyyzzz_0, g_0_xxzzzz_0_xxyyzzz_1, g_0_xxzzzz_0_xxyzzz_1, g_0_xxzzzz_0_xxyzzzz_0, g_0_xxzzzz_0_xxyzzzz_1, g_0_xxzzzz_0_xxzzzz_1, g_0_xxzzzz_0_xxzzzzz_0, g_0_xxzzzz_0_xxzzzzz_1, g_0_xxzzzz_0_xyyyyyz_0, g_0_xxzzzz_0_xyyyyyz_1, g_0_xxzzzz_0_xyyyyz_1, g_0_xxzzzz_0_xyyyyzz_0, g_0_xxzzzz_0_xyyyyzz_1, g_0_xxzzzz_0_xyyyzz_1, g_0_xxzzzz_0_xyyyzzz_0, g_0_xxzzzz_0_xyyyzzz_1, g_0_xxzzzz_0_xyyzzz_1, g_0_xxzzzz_0_xyyzzzz_0, g_0_xxzzzz_0_xyyzzzz_1, g_0_xxzzzz_0_xyzzzz_1, g_0_xxzzzz_0_xyzzzzz_0, g_0_xxzzzz_0_xyzzzzz_1, g_0_xxzzzz_0_xzzzzz_1, g_0_xxzzzz_0_xzzzzzz_0, g_0_xxzzzz_0_xzzzzzz_1, g_0_xxzzzz_0_yyyyyyy_0, g_0_xxzzzz_0_yyyyyyy_1, g_0_xxzzzz_0_yyyyyyz_0, g_0_xxzzzz_0_yyyyyyz_1, g_0_xxzzzz_0_yyyyyz_1, g_0_xxzzzz_0_yyyyyzz_0, g_0_xxzzzz_0_yyyyyzz_1, g_0_xxzzzz_0_yyyyzz_1, g_0_xxzzzz_0_yyyyzzz_0, g_0_xxzzzz_0_yyyyzzz_1, g_0_xxzzzz_0_yyyzzz_1, g_0_xxzzzz_0_yyyzzzz_0, g_0_xxzzzz_0_yyyzzzz_1, g_0_xxzzzz_0_yyzzzz_1, g_0_xxzzzz_0_yyzzzzz_0, g_0_xxzzzz_0_yyzzzzz_1, g_0_xxzzzz_0_yzzzzz_1, g_0_xxzzzz_0_yzzzzzz_0, g_0_xxzzzz_0_yzzzzzz_1, g_0_xxzzzz_0_zzzzzz_1, g_0_xxzzzz_0_zzzzzzz_0, g_0_xxzzzz_0_zzzzzzz_1, g_0_xzzzz_0_xxxxxxz_0, g_0_xzzzz_0_xxxxxxz_1, g_0_xzzzz_0_xxxxxyz_0, g_0_xzzzz_0_xxxxxyz_1, g_0_xzzzz_0_xxxxxzz_0, g_0_xzzzz_0_xxxxxzz_1, g_0_xzzzz_0_xxxxyyz_0, g_0_xzzzz_0_xxxxyyz_1, g_0_xzzzz_0_xxxxyzz_0, g_0_xzzzz_0_xxxxyzz_1, g_0_xzzzz_0_xxxxzzz_0, g_0_xzzzz_0_xxxxzzz_1, g_0_xzzzz_0_xxxyyyz_0, g_0_xzzzz_0_xxxyyyz_1, g_0_xzzzz_0_xxxyyzz_0, g_0_xzzzz_0_xxxyyzz_1, g_0_xzzzz_0_xxxyzzz_0, g_0_xzzzz_0_xxxyzzz_1, g_0_xzzzz_0_xxxzzzz_0, g_0_xzzzz_0_xxxzzzz_1, g_0_xzzzz_0_xxyyyyz_0, g_0_xzzzz_0_xxyyyyz_1, g_0_xzzzz_0_xxyyyzz_0, g_0_xzzzz_0_xxyyyzz_1, g_0_xzzzz_0_xxyyzzz_0, g_0_xzzzz_0_xxyyzzz_1, g_0_xzzzz_0_xxyzzzz_0, g_0_xzzzz_0_xxyzzzz_1, g_0_xzzzz_0_xxzzzzz_0, g_0_xzzzz_0_xxzzzzz_1, g_0_xzzzz_0_xyyyyyz_0, g_0_xzzzz_0_xyyyyyz_1, g_0_xzzzz_0_xyyyyzz_0, g_0_xzzzz_0_xyyyyzz_1, g_0_xzzzz_0_xyyyzzz_0, g_0_xzzzz_0_xyyyzzz_1, g_0_xzzzz_0_xyyzzzz_0, g_0_xzzzz_0_xyyzzzz_1, g_0_xzzzz_0_xyzzzzz_0, g_0_xzzzz_0_xyzzzzz_1, g_0_xzzzz_0_xzzzzzz_0, g_0_xzzzz_0_xzzzzzz_1, g_0_xzzzz_0_yyyyyyy_0, g_0_xzzzz_0_yyyyyyy_1, g_0_xzzzz_0_yyyyyyz_0, g_0_xzzzz_0_yyyyyyz_1, g_0_xzzzz_0_yyyyyzz_0, g_0_xzzzz_0_yyyyyzz_1, g_0_xzzzz_0_yyyyzzz_0, g_0_xzzzz_0_yyyyzzz_1, g_0_xzzzz_0_yyyzzzz_0, g_0_xzzzz_0_yyyzzzz_1, g_0_xzzzz_0_yyzzzzz_0, g_0_xzzzz_0_yyzzzzz_1, g_0_xzzzz_0_yzzzzzz_0, g_0_xzzzz_0_yzzzzzz_1, g_0_xzzzz_0_zzzzzzz_0, g_0_xzzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzz_0_xxxxxxx_0[i] = 3.0 * g_0_xxxzz_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxxxxx_0[i] * pb_z + g_0_xxxzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxxxy_0[i] = 3.0 * g_0_xxxzz_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxxxxy_0[i] * pb_z + g_0_xxxzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxxxz_0[i] = 2.0 * g_0_xzzzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxxxz_0[i] * pb_x + g_0_xxzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxxxyy_0[i] = 3.0 * g_0_xxxzz_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxxxyy_0[i] * pb_z + g_0_xxxzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxxyz_0[i] = 2.0 * g_0_xzzzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxxyz_0[i] * pb_x + g_0_xxzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxxxzz_0[i] = 2.0 * g_0_xzzzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxxzz_0[i] * pb_x + g_0_xxzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxxzz_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxxyyy_0[i] * pb_z + g_0_xxxzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xzzzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxyyz_0[i] * pb_x + g_0_xxzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxxyzz_0[i] = 2.0 * g_0_xzzzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxyzz_0[i] * pb_x + g_0_xxzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxxzzz_0[i] = 2.0 * g_0_xzzzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxzzz_0[i] * pb_x + g_0_xxzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_xxxzz_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxxyyyy_0[i] * pb_z + g_0_xxxzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxyyyz_0[i] = 2.0 * g_0_xzzzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyyyz_0[i] * pb_x + g_0_xxzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xzzzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyyzz_0[i] * pb_x + g_0_xxzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxyzzz_0[i] = 2.0 * g_0_xzzzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyzzz_0[i] * pb_x + g_0_xxzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxzzzz_0[i] = 2.0 * g_0_xzzzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxzzzz_0[i] * pb_x + g_0_xxzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyyyyy_0[i] = 3.0 * g_0_xxxzz_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xxyyyyy_0[i] * pb_z + g_0_xxxzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_xzzzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyyyz_0[i] * pb_x + g_0_xxzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_xzzzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyyzz_0[i] * pb_x + g_0_xxzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xzzzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyzzz_0[i] * pb_x + g_0_xxzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_xzzzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyzzzz_0[i] * pb_x + g_0_xxzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_xzzzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxzzzzz_0[i] * pb_x + g_0_xxzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyyyyy_0[i] = 3.0 * g_0_xxxzz_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxzzz_0_xyyyyyy_0[i] * pb_z + g_0_xxxzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xyyyyyz_0[i] = 2.0 * g_0_xzzzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyyyz_0[i] * pb_x + g_0_xxzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyyyzz_0[i] = 2.0 * g_0_xzzzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyyzz_0[i] * pb_x + g_0_xxzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyyzzz_0[i] = 2.0 * g_0_xzzzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyzzz_0[i] * pb_x + g_0_xxzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xzzzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyzzzz_0[i] * pb_x + g_0_xxzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyzzzzz_0[i] = 2.0 * g_0_xzzzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyzzzzz_0[i] * pb_x + g_0_xxzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xzzzzzz_0[i] = 2.0 * g_0_xzzzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xzzzzzz_0[i] * pb_x + g_0_xxzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyyyy_0[i] = 2.0 * g_0_xzzzz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyyyyy_0[i] * pb_x + g_0_xxzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyyyz_0[i] = 2.0 * g_0_xzzzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyyyyz_0[i] * pb_x + g_0_xxzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyyzz_0[i] = 2.0 * g_0_xzzzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyyyzz_0[i] * pb_x + g_0_xxzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyzzz_0[i] = 2.0 * g_0_xzzzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyyzzz_0[i] * pb_x + g_0_xxzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyzzzz_0[i] = 2.0 * g_0_xzzzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyyzzzz_0[i] * pb_x + g_0_xxzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xzzzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yyzzzzz_0[i] * pb_x + g_0_xxzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yzzzzzz_0[i] = 2.0 * g_0_xzzzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_yzzzzzz_0[i] * pb_x + g_0_xxzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_zzzzzzz_0[i] = 2.0 * g_0_xzzzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxzzzz_0_zzzzzzz_0[i] * pb_x + g_0_xxzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 540-576 components of targeted buffer : SKSK

    auto g_0_xxyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 540);

    auto g_0_xxyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 541);

    auto g_0_xxyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 542);

    auto g_0_xxyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 543);

    auto g_0_xxyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 544);

    auto g_0_xxyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 545);

    auto g_0_xxyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 546);

    auto g_0_xxyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 547);

    auto g_0_xxyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 548);

    auto g_0_xxyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 549);

    auto g_0_xxyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 550);

    auto g_0_xxyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 551);

    auto g_0_xxyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 552);

    auto g_0_xxyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 553);

    auto g_0_xxyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 554);

    auto g_0_xxyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 555);

    auto g_0_xxyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 556);

    auto g_0_xxyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 557);

    auto g_0_xxyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 558);

    auto g_0_xxyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 559);

    auto g_0_xxyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 560);

    auto g_0_xxyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 561);

    auto g_0_xxyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 562);

    auto g_0_xxyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 563);

    auto g_0_xxyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 564);

    auto g_0_xxyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 565);

    auto g_0_xxyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 566);

    auto g_0_xxyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 567);

    auto g_0_xxyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 568);

    auto g_0_xxyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 569);

    auto g_0_xxyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 570);

    auto g_0_xxyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 571);

    auto g_0_xxyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 572);

    auto g_0_xxyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 573);

    auto g_0_xxyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 574);

    auto g_0_xxyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 575);

    #pragma omp simd aligned(g_0_xxyyy_0_xxxxxxx_0, g_0_xxyyy_0_xxxxxxx_1, g_0_xxyyy_0_xxxxxxz_0, g_0_xxyyy_0_xxxxxxz_1, g_0_xxyyy_0_xxxxxzz_0, g_0_xxyyy_0_xxxxxzz_1, g_0_xxyyy_0_xxxxzzz_0, g_0_xxyyy_0_xxxxzzz_1, g_0_xxyyy_0_xxxzzzz_0, g_0_xxyyy_0_xxxzzzz_1, g_0_xxyyy_0_xxzzzzz_0, g_0_xxyyy_0_xxzzzzz_1, g_0_xxyyy_0_xzzzzzz_0, g_0_xxyyy_0_xzzzzzz_1, g_0_xxyyyy_0_xxxxxxx_0, g_0_xxyyyy_0_xxxxxxx_1, g_0_xxyyyy_0_xxxxxxz_0, g_0_xxyyyy_0_xxxxxxz_1, g_0_xxyyyy_0_xxxxxzz_0, g_0_xxyyyy_0_xxxxxzz_1, g_0_xxyyyy_0_xxxxzzz_0, g_0_xxyyyy_0_xxxxzzz_1, g_0_xxyyyy_0_xxxzzzz_0, g_0_xxyyyy_0_xxxzzzz_1, g_0_xxyyyy_0_xxzzzzz_0, g_0_xxyyyy_0_xxzzzzz_1, g_0_xxyyyy_0_xzzzzzz_0, g_0_xxyyyy_0_xzzzzzz_1, g_0_xxyyyyy_0_xxxxxxx_0, g_0_xxyyyyy_0_xxxxxxy_0, g_0_xxyyyyy_0_xxxxxxz_0, g_0_xxyyyyy_0_xxxxxyy_0, g_0_xxyyyyy_0_xxxxxyz_0, g_0_xxyyyyy_0_xxxxxzz_0, g_0_xxyyyyy_0_xxxxyyy_0, g_0_xxyyyyy_0_xxxxyyz_0, g_0_xxyyyyy_0_xxxxyzz_0, g_0_xxyyyyy_0_xxxxzzz_0, g_0_xxyyyyy_0_xxxyyyy_0, g_0_xxyyyyy_0_xxxyyyz_0, g_0_xxyyyyy_0_xxxyyzz_0, g_0_xxyyyyy_0_xxxyzzz_0, g_0_xxyyyyy_0_xxxzzzz_0, g_0_xxyyyyy_0_xxyyyyy_0, g_0_xxyyyyy_0_xxyyyyz_0, g_0_xxyyyyy_0_xxyyyzz_0, g_0_xxyyyyy_0_xxyyzzz_0, g_0_xxyyyyy_0_xxyzzzz_0, g_0_xxyyyyy_0_xxzzzzz_0, g_0_xxyyyyy_0_xyyyyyy_0, g_0_xxyyyyy_0_xyyyyyz_0, g_0_xxyyyyy_0_xyyyyzz_0, g_0_xxyyyyy_0_xyyyzzz_0, g_0_xxyyyyy_0_xyyzzzz_0, g_0_xxyyyyy_0_xyzzzzz_0, g_0_xxyyyyy_0_xzzzzzz_0, g_0_xxyyyyy_0_yyyyyyy_0, g_0_xxyyyyy_0_yyyyyyz_0, g_0_xxyyyyy_0_yyyyyzz_0, g_0_xxyyyyy_0_yyyyzzz_0, g_0_xxyyyyy_0_yyyzzzz_0, g_0_xxyyyyy_0_yyzzzzz_0, g_0_xxyyyyy_0_yzzzzzz_0, g_0_xxyyyyy_0_zzzzzzz_0, g_0_xyyyyy_0_xxxxxxy_0, g_0_xyyyyy_0_xxxxxxy_1, g_0_xyyyyy_0_xxxxxy_1, g_0_xyyyyy_0_xxxxxyy_0, g_0_xyyyyy_0_xxxxxyy_1, g_0_xyyyyy_0_xxxxxyz_0, g_0_xyyyyy_0_xxxxxyz_1, g_0_xyyyyy_0_xxxxyy_1, g_0_xyyyyy_0_xxxxyyy_0, g_0_xyyyyy_0_xxxxyyy_1, g_0_xyyyyy_0_xxxxyyz_0, g_0_xyyyyy_0_xxxxyyz_1, g_0_xyyyyy_0_xxxxyz_1, g_0_xyyyyy_0_xxxxyzz_0, g_0_xyyyyy_0_xxxxyzz_1, g_0_xyyyyy_0_xxxyyy_1, g_0_xyyyyy_0_xxxyyyy_0, g_0_xyyyyy_0_xxxyyyy_1, g_0_xyyyyy_0_xxxyyyz_0, g_0_xyyyyy_0_xxxyyyz_1, g_0_xyyyyy_0_xxxyyz_1, g_0_xyyyyy_0_xxxyyzz_0, g_0_xyyyyy_0_xxxyyzz_1, g_0_xyyyyy_0_xxxyzz_1, g_0_xyyyyy_0_xxxyzzz_0, g_0_xyyyyy_0_xxxyzzz_1, g_0_xyyyyy_0_xxyyyy_1, g_0_xyyyyy_0_xxyyyyy_0, g_0_xyyyyy_0_xxyyyyy_1, g_0_xyyyyy_0_xxyyyyz_0, g_0_xyyyyy_0_xxyyyyz_1, g_0_xyyyyy_0_xxyyyz_1, g_0_xyyyyy_0_xxyyyzz_0, g_0_xyyyyy_0_xxyyyzz_1, g_0_xyyyyy_0_xxyyzz_1, g_0_xyyyyy_0_xxyyzzz_0, g_0_xyyyyy_0_xxyyzzz_1, g_0_xyyyyy_0_xxyzzz_1, g_0_xyyyyy_0_xxyzzzz_0, g_0_xyyyyy_0_xxyzzzz_1, g_0_xyyyyy_0_xyyyyy_1, g_0_xyyyyy_0_xyyyyyy_0, g_0_xyyyyy_0_xyyyyyy_1, g_0_xyyyyy_0_xyyyyyz_0, g_0_xyyyyy_0_xyyyyyz_1, g_0_xyyyyy_0_xyyyyz_1, g_0_xyyyyy_0_xyyyyzz_0, g_0_xyyyyy_0_xyyyyzz_1, g_0_xyyyyy_0_xyyyzz_1, g_0_xyyyyy_0_xyyyzzz_0, g_0_xyyyyy_0_xyyyzzz_1, g_0_xyyyyy_0_xyyzzz_1, g_0_xyyyyy_0_xyyzzzz_0, g_0_xyyyyy_0_xyyzzzz_1, g_0_xyyyyy_0_xyzzzz_1, g_0_xyyyyy_0_xyzzzzz_0, g_0_xyyyyy_0_xyzzzzz_1, g_0_xyyyyy_0_yyyyyy_1, g_0_xyyyyy_0_yyyyyyy_0, g_0_xyyyyy_0_yyyyyyy_1, g_0_xyyyyy_0_yyyyyyz_0, g_0_xyyyyy_0_yyyyyyz_1, g_0_xyyyyy_0_yyyyyz_1, g_0_xyyyyy_0_yyyyyzz_0, g_0_xyyyyy_0_yyyyyzz_1, g_0_xyyyyy_0_yyyyzz_1, g_0_xyyyyy_0_yyyyzzz_0, g_0_xyyyyy_0_yyyyzzz_1, g_0_xyyyyy_0_yyyzzz_1, g_0_xyyyyy_0_yyyzzzz_0, g_0_xyyyyy_0_yyyzzzz_1, g_0_xyyyyy_0_yyzzzz_1, g_0_xyyyyy_0_yyzzzzz_0, g_0_xyyyyy_0_yyzzzzz_1, g_0_xyyyyy_0_yzzzzz_1, g_0_xyyyyy_0_yzzzzzz_0, g_0_xyyyyy_0_yzzzzzz_1, g_0_xyyyyy_0_zzzzzzz_0, g_0_xyyyyy_0_zzzzzzz_1, g_0_yyyyy_0_xxxxxxy_0, g_0_yyyyy_0_xxxxxxy_1, g_0_yyyyy_0_xxxxxyy_0, g_0_yyyyy_0_xxxxxyy_1, g_0_yyyyy_0_xxxxxyz_0, g_0_yyyyy_0_xxxxxyz_1, g_0_yyyyy_0_xxxxyyy_0, g_0_yyyyy_0_xxxxyyy_1, g_0_yyyyy_0_xxxxyyz_0, g_0_yyyyy_0_xxxxyyz_1, g_0_yyyyy_0_xxxxyzz_0, g_0_yyyyy_0_xxxxyzz_1, g_0_yyyyy_0_xxxyyyy_0, g_0_yyyyy_0_xxxyyyy_1, g_0_yyyyy_0_xxxyyyz_0, g_0_yyyyy_0_xxxyyyz_1, g_0_yyyyy_0_xxxyyzz_0, g_0_yyyyy_0_xxxyyzz_1, g_0_yyyyy_0_xxxyzzz_0, g_0_yyyyy_0_xxxyzzz_1, g_0_yyyyy_0_xxyyyyy_0, g_0_yyyyy_0_xxyyyyy_1, g_0_yyyyy_0_xxyyyyz_0, g_0_yyyyy_0_xxyyyyz_1, g_0_yyyyy_0_xxyyyzz_0, g_0_yyyyy_0_xxyyyzz_1, g_0_yyyyy_0_xxyyzzz_0, g_0_yyyyy_0_xxyyzzz_1, g_0_yyyyy_0_xxyzzzz_0, g_0_yyyyy_0_xxyzzzz_1, g_0_yyyyy_0_xyyyyyy_0, g_0_yyyyy_0_xyyyyyy_1, g_0_yyyyy_0_xyyyyyz_0, g_0_yyyyy_0_xyyyyyz_1, g_0_yyyyy_0_xyyyyzz_0, g_0_yyyyy_0_xyyyyzz_1, g_0_yyyyy_0_xyyyzzz_0, g_0_yyyyy_0_xyyyzzz_1, g_0_yyyyy_0_xyyzzzz_0, g_0_yyyyy_0_xyyzzzz_1, g_0_yyyyy_0_xyzzzzz_0, g_0_yyyyy_0_xyzzzzz_1, g_0_yyyyy_0_yyyyyyy_0, g_0_yyyyy_0_yyyyyyy_1, g_0_yyyyy_0_yyyyyyz_0, g_0_yyyyy_0_yyyyyyz_1, g_0_yyyyy_0_yyyyyzz_0, g_0_yyyyy_0_yyyyyzz_1, g_0_yyyyy_0_yyyyzzz_0, g_0_yyyyy_0_yyyyzzz_1, g_0_yyyyy_0_yyyzzzz_0, g_0_yyyyy_0_yyyzzzz_1, g_0_yyyyy_0_yyzzzzz_0, g_0_yyyyy_0_yyzzzzz_1, g_0_yyyyy_0_yzzzzzz_0, g_0_yyyyy_0_yzzzzzz_1, g_0_yyyyy_0_zzzzzzz_0, g_0_yyyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyy_0_xxxxxxx_0[i] = 4.0 * g_0_xxyyy_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxxxxx_0[i] * pb_y + g_0_xxyyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxxxxy_0[i] = g_0_yyyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxxxy_0[i] * pb_x + g_0_xyyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxxxz_0[i] = 4.0 * g_0_xxyyy_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxxxxz_0[i] * pb_y + g_0_xxyyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxxxyy_0[i] = g_0_yyyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxxyy_0[i] * pb_x + g_0_xyyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxxyz_0[i] = g_0_yyyyy_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxxyz_0[i] * pb_x + g_0_xyyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxxzz_0[i] = 4.0 * g_0_xxyyy_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxxxzz_0[i] * pb_y + g_0_xxyyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxxyyy_0[i] = g_0_yyyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxyyy_0[i] * pb_x + g_0_xyyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxyyz_0[i] = g_0_yyyyy_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxyyz_0[i] * pb_x + g_0_xyyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxyzz_0[i] = g_0_yyyyy_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxyzz_0[i] * pb_x + g_0_xyyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxzzz_0[i] = 4.0 * g_0_xxyyy_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxxzzz_0[i] * pb_y + g_0_xxyyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxyyyy_0[i] = g_0_yyyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyyyy_0[i] * pb_x + g_0_xyyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxyyyz_0[i] = g_0_yyyyy_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyyyz_0[i] * pb_x + g_0_xyyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxyyzz_0[i] = g_0_yyyyy_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyyzz_0[i] * pb_x + g_0_xyyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxyzzz_0[i] = g_0_yyyyy_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyzzz_0[i] * pb_x + g_0_xyyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxzzzz_0[i] = 4.0 * g_0_xxyyy_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxxzzzz_0[i] * pb_y + g_0_xxyyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxyyyyy_0[i] = g_0_yyyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyyyy_0[i] * pb_x + g_0_xyyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyyyyz_0[i] = g_0_yyyyy_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyyyz_0[i] * pb_x + g_0_xyyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyyyzz_0[i] = g_0_yyyyy_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyyzz_0[i] * pb_x + g_0_xyyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyyzzz_0[i] = g_0_yyyyy_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyzzz_0[i] * pb_x + g_0_xyyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyzzzz_0[i] = g_0_yyyyy_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyzzzz_0[i] * pb_x + g_0_xyyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxzzzzz_0[i] = 4.0 * g_0_xxyyy_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xxzzzzz_0[i] * pb_y + g_0_xxyyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xyyyyyy_0[i] = g_0_yyyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyyyyyy_0[i] * pb_x + g_0_xyyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyyyyz_0[i] = g_0_yyyyy_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyyyyyz_0[i] * pb_x + g_0_xyyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyyyzz_0[i] = g_0_yyyyy_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyyyyzz_0[i] * pb_x + g_0_xyyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyyzzz_0[i] = g_0_yyyyy_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyyyzzz_0[i] * pb_x + g_0_xyyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyzzzz_0[i] = g_0_yyyyy_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyyzzzz_0[i] * pb_x + g_0_xyyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyzzzzz_0[i] = g_0_yyyyy_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xyzzzzz_0[i] * pb_x + g_0_xyyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xzzzzzz_0[i] = 4.0 * g_0_xxyyy_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyy_0_xzzzzzz_0[i] * pb_y + g_0_xxyyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_yyyyyyy_0[i] = g_0_yyyyy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyyyy_0[i] * pb_x + g_0_xyyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyyyyz_0[i] = g_0_yyyyy_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyyyz_0[i] * pb_x + g_0_xyyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyyyzz_0[i] = g_0_yyyyy_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyyzz_0[i] * pb_x + g_0_xyyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyyzzz_0[i] = g_0_yyyyy_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyzzz_0[i] * pb_x + g_0_xyyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyzzzz_0[i] = g_0_yyyyy_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyzzzz_0[i] * pb_x + g_0_xyyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyzzzzz_0[i] = g_0_yyyyy_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyzzzzz_0[i] * pb_x + g_0_xyyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yzzzzzz_0[i] = g_0_yyyyy_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzzzzzz_0[i] * pb_x + g_0_xyyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_zzzzzzz_0[i] = g_0_yyyyy_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_zzzzzzz_0[i] * pb_x + g_0_xyyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 576-612 components of targeted buffer : SKSK

    auto g_0_xxyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 576);

    auto g_0_xxyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 577);

    auto g_0_xxyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 578);

    auto g_0_xxyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 579);

    auto g_0_xxyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 580);

    auto g_0_xxyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 581);

    auto g_0_xxyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 582);

    auto g_0_xxyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 583);

    auto g_0_xxyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 584);

    auto g_0_xxyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 585);

    auto g_0_xxyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 586);

    auto g_0_xxyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 587);

    auto g_0_xxyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 588);

    auto g_0_xxyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 589);

    auto g_0_xxyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 590);

    auto g_0_xxyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 591);

    auto g_0_xxyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 592);

    auto g_0_xxyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 593);

    auto g_0_xxyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 594);

    auto g_0_xxyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 595);

    auto g_0_xxyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 596);

    auto g_0_xxyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 597);

    auto g_0_xxyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 598);

    auto g_0_xxyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 599);

    auto g_0_xxyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 600);

    auto g_0_xxyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 601);

    auto g_0_xxyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 602);

    auto g_0_xxyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 603);

    auto g_0_xxyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 604);

    auto g_0_xxyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 605);

    auto g_0_xxyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 606);

    auto g_0_xxyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 607);

    auto g_0_xxyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 608);

    auto g_0_xxyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 609);

    auto g_0_xxyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 610);

    auto g_0_xxyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 611);

    #pragma omp simd aligned(g_0_xxyyyy_0_xxxxxx_1, g_0_xxyyyy_0_xxxxxxx_0, g_0_xxyyyy_0_xxxxxxx_1, g_0_xxyyyy_0_xxxxxxy_0, g_0_xxyyyy_0_xxxxxxy_1, g_0_xxyyyy_0_xxxxxxz_0, g_0_xxyyyy_0_xxxxxxz_1, g_0_xxyyyy_0_xxxxxy_1, g_0_xxyyyy_0_xxxxxyy_0, g_0_xxyyyy_0_xxxxxyy_1, g_0_xxyyyy_0_xxxxxyz_0, g_0_xxyyyy_0_xxxxxyz_1, g_0_xxyyyy_0_xxxxxz_1, g_0_xxyyyy_0_xxxxxzz_0, g_0_xxyyyy_0_xxxxxzz_1, g_0_xxyyyy_0_xxxxyy_1, g_0_xxyyyy_0_xxxxyyy_0, g_0_xxyyyy_0_xxxxyyy_1, g_0_xxyyyy_0_xxxxyyz_0, g_0_xxyyyy_0_xxxxyyz_1, g_0_xxyyyy_0_xxxxyz_1, g_0_xxyyyy_0_xxxxyzz_0, g_0_xxyyyy_0_xxxxyzz_1, g_0_xxyyyy_0_xxxxzz_1, g_0_xxyyyy_0_xxxxzzz_0, g_0_xxyyyy_0_xxxxzzz_1, g_0_xxyyyy_0_xxxyyy_1, g_0_xxyyyy_0_xxxyyyy_0, g_0_xxyyyy_0_xxxyyyy_1, g_0_xxyyyy_0_xxxyyyz_0, g_0_xxyyyy_0_xxxyyyz_1, g_0_xxyyyy_0_xxxyyz_1, g_0_xxyyyy_0_xxxyyzz_0, g_0_xxyyyy_0_xxxyyzz_1, g_0_xxyyyy_0_xxxyzz_1, g_0_xxyyyy_0_xxxyzzz_0, g_0_xxyyyy_0_xxxyzzz_1, g_0_xxyyyy_0_xxxzzz_1, g_0_xxyyyy_0_xxxzzzz_0, g_0_xxyyyy_0_xxxzzzz_1, g_0_xxyyyy_0_xxyyyy_1, g_0_xxyyyy_0_xxyyyyy_0, g_0_xxyyyy_0_xxyyyyy_1, g_0_xxyyyy_0_xxyyyyz_0, g_0_xxyyyy_0_xxyyyyz_1, g_0_xxyyyy_0_xxyyyz_1, g_0_xxyyyy_0_xxyyyzz_0, g_0_xxyyyy_0_xxyyyzz_1, g_0_xxyyyy_0_xxyyzz_1, g_0_xxyyyy_0_xxyyzzz_0, g_0_xxyyyy_0_xxyyzzz_1, g_0_xxyyyy_0_xxyzzz_1, g_0_xxyyyy_0_xxyzzzz_0, g_0_xxyyyy_0_xxyzzzz_1, g_0_xxyyyy_0_xxzzzz_1, g_0_xxyyyy_0_xxzzzzz_0, g_0_xxyyyy_0_xxzzzzz_1, g_0_xxyyyy_0_xyyyyy_1, g_0_xxyyyy_0_xyyyyyy_0, g_0_xxyyyy_0_xyyyyyy_1, g_0_xxyyyy_0_xyyyyyz_0, g_0_xxyyyy_0_xyyyyyz_1, g_0_xxyyyy_0_xyyyyz_1, g_0_xxyyyy_0_xyyyyzz_0, g_0_xxyyyy_0_xyyyyzz_1, g_0_xxyyyy_0_xyyyzz_1, g_0_xxyyyy_0_xyyyzzz_0, g_0_xxyyyy_0_xyyyzzz_1, g_0_xxyyyy_0_xyyzzz_1, g_0_xxyyyy_0_xyyzzzz_0, g_0_xxyyyy_0_xyyzzzz_1, g_0_xxyyyy_0_xyzzzz_1, g_0_xxyyyy_0_xyzzzzz_0, g_0_xxyyyy_0_xyzzzzz_1, g_0_xxyyyy_0_xzzzzz_1, g_0_xxyyyy_0_xzzzzzz_0, g_0_xxyyyy_0_xzzzzzz_1, g_0_xxyyyy_0_yyyyyy_1, g_0_xxyyyy_0_yyyyyyy_0, g_0_xxyyyy_0_yyyyyyy_1, g_0_xxyyyy_0_yyyyyyz_0, g_0_xxyyyy_0_yyyyyyz_1, g_0_xxyyyy_0_yyyyyz_1, g_0_xxyyyy_0_yyyyyzz_0, g_0_xxyyyy_0_yyyyyzz_1, g_0_xxyyyy_0_yyyyzz_1, g_0_xxyyyy_0_yyyyzzz_0, g_0_xxyyyy_0_yyyyzzz_1, g_0_xxyyyy_0_yyyzzz_1, g_0_xxyyyy_0_yyyzzzz_0, g_0_xxyyyy_0_yyyzzzz_1, g_0_xxyyyy_0_yyzzzz_1, g_0_xxyyyy_0_yyzzzzz_0, g_0_xxyyyy_0_yyzzzzz_1, g_0_xxyyyy_0_yzzzzz_1, g_0_xxyyyy_0_yzzzzzz_0, g_0_xxyyyy_0_yzzzzzz_1, g_0_xxyyyy_0_zzzzzz_1, g_0_xxyyyy_0_zzzzzzz_0, g_0_xxyyyy_0_zzzzzzz_1, g_0_xxyyyyz_0_xxxxxxx_0, g_0_xxyyyyz_0_xxxxxxy_0, g_0_xxyyyyz_0_xxxxxxz_0, g_0_xxyyyyz_0_xxxxxyy_0, g_0_xxyyyyz_0_xxxxxyz_0, g_0_xxyyyyz_0_xxxxxzz_0, g_0_xxyyyyz_0_xxxxyyy_0, g_0_xxyyyyz_0_xxxxyyz_0, g_0_xxyyyyz_0_xxxxyzz_0, g_0_xxyyyyz_0_xxxxzzz_0, g_0_xxyyyyz_0_xxxyyyy_0, g_0_xxyyyyz_0_xxxyyyz_0, g_0_xxyyyyz_0_xxxyyzz_0, g_0_xxyyyyz_0_xxxyzzz_0, g_0_xxyyyyz_0_xxxzzzz_0, g_0_xxyyyyz_0_xxyyyyy_0, g_0_xxyyyyz_0_xxyyyyz_0, g_0_xxyyyyz_0_xxyyyzz_0, g_0_xxyyyyz_0_xxyyzzz_0, g_0_xxyyyyz_0_xxyzzzz_0, g_0_xxyyyyz_0_xxzzzzz_0, g_0_xxyyyyz_0_xyyyyyy_0, g_0_xxyyyyz_0_xyyyyyz_0, g_0_xxyyyyz_0_xyyyyzz_0, g_0_xxyyyyz_0_xyyyzzz_0, g_0_xxyyyyz_0_xyyzzzz_0, g_0_xxyyyyz_0_xyzzzzz_0, g_0_xxyyyyz_0_xzzzzzz_0, g_0_xxyyyyz_0_yyyyyyy_0, g_0_xxyyyyz_0_yyyyyyz_0, g_0_xxyyyyz_0_yyyyyzz_0, g_0_xxyyyyz_0_yyyyzzz_0, g_0_xxyyyyz_0_yyyzzzz_0, g_0_xxyyyyz_0_yyzzzzz_0, g_0_xxyyyyz_0_yzzzzzz_0, g_0_xxyyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyz_0_xxxxxxx_0[i] = g_0_xxyyyy_0_xxxxxxx_0[i] * pb_z + g_0_xxyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxxxy_0[i] = g_0_xxyyyy_0_xxxxxxy_0[i] * pb_z + g_0_xxyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxxxz_0[i] = g_0_xxyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxxxz_0[i] * pb_z + g_0_xxyyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxxyy_0[i] = g_0_xxyyyy_0_xxxxxyy_0[i] * pb_z + g_0_xxyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxxyz_0[i] = g_0_xxyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxxyz_0[i] * pb_z + g_0_xxyyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxxzz_0[i] * pb_z + g_0_xxyyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxyyy_0[i] = g_0_xxyyyy_0_xxxxyyy_0[i] * pb_z + g_0_xxyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxyyz_0[i] = g_0_xxyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxyyz_0[i] * pb_z + g_0_xxyyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxyzz_0[i] * pb_z + g_0_xxyyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxzzz_0[i] * pb_z + g_0_xxyyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyyyy_0[i] = g_0_xxyyyy_0_xxxyyyy_0[i] * pb_z + g_0_xxyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyyyz_0[i] = g_0_xxyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyyyz_0[i] * pb_z + g_0_xxyyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyyzz_0[i] * pb_z + g_0_xxyyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyzzz_0[i] * pb_z + g_0_xxyyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxzzzz_0[i] * pb_z + g_0_xxyyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyyyy_0[i] = g_0_xxyyyy_0_xxyyyyy_0[i] * pb_z + g_0_xxyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyyyz_0[i] = g_0_xxyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyyyz_0[i] * pb_z + g_0_xxyyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyyzz_0[i] * pb_z + g_0_xxyyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyzzz_0[i] * pb_z + g_0_xxyyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyzzzz_0[i] * pb_z + g_0_xxyyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxzzzzz_0[i] * pb_z + g_0_xxyyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyyyy_0[i] = g_0_xxyyyy_0_xyyyyyy_0[i] * pb_z + g_0_xxyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyyyz_0[i] = g_0_xxyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyyyz_0[i] * pb_z + g_0_xxyyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyyzz_0[i] * pb_z + g_0_xxyyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyzzz_0[i] * pb_z + g_0_xxyyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyzzzz_0[i] * pb_z + g_0_xxyyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyzzzzz_0[i] * pb_z + g_0_xxyyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xzzzzzz_0[i] * pb_z + g_0_xxyyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyyyy_0[i] = g_0_xxyyyy_0_yyyyyyy_0[i] * pb_z + g_0_xxyyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyyyz_0[i] = g_0_xxyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyyyyz_0[i] * pb_z + g_0_xxyyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyyyzz_0[i] * pb_z + g_0_xxyyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyyzzz_0[i] * pb_z + g_0_xxyyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyzzzz_0[i] * pb_z + g_0_xxyyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyzzzzz_0[i] * pb_z + g_0_xxyyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yzzzzzz_0[i] * pb_z + g_0_xxyyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_zzzzzzz_0[i] * pb_z + g_0_xxyyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 612-648 components of targeted buffer : SKSK

    auto g_0_xxyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 612);

    auto g_0_xxyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 613);

    auto g_0_xxyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 614);

    auto g_0_xxyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 615);

    auto g_0_xxyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 616);

    auto g_0_xxyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 617);

    auto g_0_xxyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 618);

    auto g_0_xxyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 619);

    auto g_0_xxyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 620);

    auto g_0_xxyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 621);

    auto g_0_xxyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 622);

    auto g_0_xxyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 623);

    auto g_0_xxyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 624);

    auto g_0_xxyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 625);

    auto g_0_xxyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 626);

    auto g_0_xxyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 627);

    auto g_0_xxyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 628);

    auto g_0_xxyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 629);

    auto g_0_xxyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 630);

    auto g_0_xxyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 631);

    auto g_0_xxyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 632);

    auto g_0_xxyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 633);

    auto g_0_xxyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 634);

    auto g_0_xxyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 635);

    auto g_0_xxyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 636);

    auto g_0_xxyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 637);

    auto g_0_xxyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 638);

    auto g_0_xxyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 639);

    auto g_0_xxyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 640);

    auto g_0_xxyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 641);

    auto g_0_xxyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 642);

    auto g_0_xxyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 643);

    auto g_0_xxyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 644);

    auto g_0_xxyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 645);

    auto g_0_xxyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 646);

    auto g_0_xxyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 647);

    #pragma omp simd aligned(g_0_xxyyy_0_xxxxxxy_0, g_0_xxyyy_0_xxxxxxy_1, g_0_xxyyy_0_xxxxxyy_0, g_0_xxyyy_0_xxxxxyy_1, g_0_xxyyy_0_xxxxyyy_0, g_0_xxyyy_0_xxxxyyy_1, g_0_xxyyy_0_xxxyyyy_0, g_0_xxyyy_0_xxxyyyy_1, g_0_xxyyy_0_xxyyyyy_0, g_0_xxyyy_0_xxyyyyy_1, g_0_xxyyy_0_xyyyyyy_0, g_0_xxyyy_0_xyyyyyy_1, g_0_xxyyyz_0_xxxxxxy_0, g_0_xxyyyz_0_xxxxxxy_1, g_0_xxyyyz_0_xxxxxyy_0, g_0_xxyyyz_0_xxxxxyy_1, g_0_xxyyyz_0_xxxxyyy_0, g_0_xxyyyz_0_xxxxyyy_1, g_0_xxyyyz_0_xxxyyyy_0, g_0_xxyyyz_0_xxxyyyy_1, g_0_xxyyyz_0_xxyyyyy_0, g_0_xxyyyz_0_xxyyyyy_1, g_0_xxyyyz_0_xyyyyyy_0, g_0_xxyyyz_0_xyyyyyy_1, g_0_xxyyyzz_0_xxxxxxx_0, g_0_xxyyyzz_0_xxxxxxy_0, g_0_xxyyyzz_0_xxxxxxz_0, g_0_xxyyyzz_0_xxxxxyy_0, g_0_xxyyyzz_0_xxxxxyz_0, g_0_xxyyyzz_0_xxxxxzz_0, g_0_xxyyyzz_0_xxxxyyy_0, g_0_xxyyyzz_0_xxxxyyz_0, g_0_xxyyyzz_0_xxxxyzz_0, g_0_xxyyyzz_0_xxxxzzz_0, g_0_xxyyyzz_0_xxxyyyy_0, g_0_xxyyyzz_0_xxxyyyz_0, g_0_xxyyyzz_0_xxxyyzz_0, g_0_xxyyyzz_0_xxxyzzz_0, g_0_xxyyyzz_0_xxxzzzz_0, g_0_xxyyyzz_0_xxyyyyy_0, g_0_xxyyyzz_0_xxyyyyz_0, g_0_xxyyyzz_0_xxyyyzz_0, g_0_xxyyyzz_0_xxyyzzz_0, g_0_xxyyyzz_0_xxyzzzz_0, g_0_xxyyyzz_0_xxzzzzz_0, g_0_xxyyyzz_0_xyyyyyy_0, g_0_xxyyyzz_0_xyyyyyz_0, g_0_xxyyyzz_0_xyyyyzz_0, g_0_xxyyyzz_0_xyyyzzz_0, g_0_xxyyyzz_0_xyyzzzz_0, g_0_xxyyyzz_0_xyzzzzz_0, g_0_xxyyyzz_0_xzzzzzz_0, g_0_xxyyyzz_0_yyyyyyy_0, g_0_xxyyyzz_0_yyyyyyz_0, g_0_xxyyyzz_0_yyyyyzz_0, g_0_xxyyyzz_0_yyyyzzz_0, g_0_xxyyyzz_0_yyyzzzz_0, g_0_xxyyyzz_0_yyzzzzz_0, g_0_xxyyyzz_0_yzzzzzz_0, g_0_xxyyyzz_0_zzzzzzz_0, g_0_xxyyzz_0_xxxxxxx_0, g_0_xxyyzz_0_xxxxxxx_1, g_0_xxyyzz_0_xxxxxxz_0, g_0_xxyyzz_0_xxxxxxz_1, g_0_xxyyzz_0_xxxxxzz_0, g_0_xxyyzz_0_xxxxxzz_1, g_0_xxyyzz_0_xxxxzzz_0, g_0_xxyyzz_0_xxxxzzz_1, g_0_xxyyzz_0_xxxzzzz_0, g_0_xxyyzz_0_xxxzzzz_1, g_0_xxyyzz_0_xxzzzzz_0, g_0_xxyyzz_0_xxzzzzz_1, g_0_xxyyzz_0_xzzzzzz_0, g_0_xxyyzz_0_xzzzzzz_1, g_0_xxyzz_0_xxxxxxx_0, g_0_xxyzz_0_xxxxxxx_1, g_0_xxyzz_0_xxxxxxz_0, g_0_xxyzz_0_xxxxxxz_1, g_0_xxyzz_0_xxxxxzz_0, g_0_xxyzz_0_xxxxxzz_1, g_0_xxyzz_0_xxxxzzz_0, g_0_xxyzz_0_xxxxzzz_1, g_0_xxyzz_0_xxxzzzz_0, g_0_xxyzz_0_xxxzzzz_1, g_0_xxyzz_0_xxzzzzz_0, g_0_xxyzz_0_xxzzzzz_1, g_0_xxyzz_0_xzzzzzz_0, g_0_xxyzz_0_xzzzzzz_1, g_0_xyyyzz_0_xxxxxyz_0, g_0_xyyyzz_0_xxxxxyz_1, g_0_xyyyzz_0_xxxxyyz_0, g_0_xyyyzz_0_xxxxyyz_1, g_0_xyyyzz_0_xxxxyz_1, g_0_xyyyzz_0_xxxxyzz_0, g_0_xyyyzz_0_xxxxyzz_1, g_0_xyyyzz_0_xxxyyyz_0, g_0_xyyyzz_0_xxxyyyz_1, g_0_xyyyzz_0_xxxyyz_1, g_0_xyyyzz_0_xxxyyzz_0, g_0_xyyyzz_0_xxxyyzz_1, g_0_xyyyzz_0_xxxyzz_1, g_0_xyyyzz_0_xxxyzzz_0, g_0_xyyyzz_0_xxxyzzz_1, g_0_xyyyzz_0_xxyyyyz_0, g_0_xyyyzz_0_xxyyyyz_1, g_0_xyyyzz_0_xxyyyz_1, g_0_xyyyzz_0_xxyyyzz_0, g_0_xyyyzz_0_xxyyyzz_1, g_0_xyyyzz_0_xxyyzz_1, g_0_xyyyzz_0_xxyyzzz_0, g_0_xyyyzz_0_xxyyzzz_1, g_0_xyyyzz_0_xxyzzz_1, g_0_xyyyzz_0_xxyzzzz_0, g_0_xyyyzz_0_xxyzzzz_1, g_0_xyyyzz_0_xyyyyyz_0, g_0_xyyyzz_0_xyyyyyz_1, g_0_xyyyzz_0_xyyyyz_1, g_0_xyyyzz_0_xyyyyzz_0, g_0_xyyyzz_0_xyyyyzz_1, g_0_xyyyzz_0_xyyyzz_1, g_0_xyyyzz_0_xyyyzzz_0, g_0_xyyyzz_0_xyyyzzz_1, g_0_xyyyzz_0_xyyzzz_1, g_0_xyyyzz_0_xyyzzzz_0, g_0_xyyyzz_0_xyyzzzz_1, g_0_xyyyzz_0_xyzzzz_1, g_0_xyyyzz_0_xyzzzzz_0, g_0_xyyyzz_0_xyzzzzz_1, g_0_xyyyzz_0_yyyyyyy_0, g_0_xyyyzz_0_yyyyyyy_1, g_0_xyyyzz_0_yyyyyyz_0, g_0_xyyyzz_0_yyyyyyz_1, g_0_xyyyzz_0_yyyyyz_1, g_0_xyyyzz_0_yyyyyzz_0, g_0_xyyyzz_0_yyyyyzz_1, g_0_xyyyzz_0_yyyyzz_1, g_0_xyyyzz_0_yyyyzzz_0, g_0_xyyyzz_0_yyyyzzz_1, g_0_xyyyzz_0_yyyzzz_1, g_0_xyyyzz_0_yyyzzzz_0, g_0_xyyyzz_0_yyyzzzz_1, g_0_xyyyzz_0_yyzzzz_1, g_0_xyyyzz_0_yyzzzzz_0, g_0_xyyyzz_0_yyzzzzz_1, g_0_xyyyzz_0_yzzzzz_1, g_0_xyyyzz_0_yzzzzzz_0, g_0_xyyyzz_0_yzzzzzz_1, g_0_xyyyzz_0_zzzzzzz_0, g_0_xyyyzz_0_zzzzzzz_1, g_0_yyyzz_0_xxxxxyz_0, g_0_yyyzz_0_xxxxxyz_1, g_0_yyyzz_0_xxxxyyz_0, g_0_yyyzz_0_xxxxyyz_1, g_0_yyyzz_0_xxxxyzz_0, g_0_yyyzz_0_xxxxyzz_1, g_0_yyyzz_0_xxxyyyz_0, g_0_yyyzz_0_xxxyyyz_1, g_0_yyyzz_0_xxxyyzz_0, g_0_yyyzz_0_xxxyyzz_1, g_0_yyyzz_0_xxxyzzz_0, g_0_yyyzz_0_xxxyzzz_1, g_0_yyyzz_0_xxyyyyz_0, g_0_yyyzz_0_xxyyyyz_1, g_0_yyyzz_0_xxyyyzz_0, g_0_yyyzz_0_xxyyyzz_1, g_0_yyyzz_0_xxyyzzz_0, g_0_yyyzz_0_xxyyzzz_1, g_0_yyyzz_0_xxyzzzz_0, g_0_yyyzz_0_xxyzzzz_1, g_0_yyyzz_0_xyyyyyz_0, g_0_yyyzz_0_xyyyyyz_1, g_0_yyyzz_0_xyyyyzz_0, g_0_yyyzz_0_xyyyyzz_1, g_0_yyyzz_0_xyyyzzz_0, g_0_yyyzz_0_xyyyzzz_1, g_0_yyyzz_0_xyyzzzz_0, g_0_yyyzz_0_xyyzzzz_1, g_0_yyyzz_0_xyzzzzz_0, g_0_yyyzz_0_xyzzzzz_1, g_0_yyyzz_0_yyyyyyy_0, g_0_yyyzz_0_yyyyyyy_1, g_0_yyyzz_0_yyyyyyz_0, g_0_yyyzz_0_yyyyyyz_1, g_0_yyyzz_0_yyyyyzz_0, g_0_yyyzz_0_yyyyyzz_1, g_0_yyyzz_0_yyyyzzz_0, g_0_yyyzz_0_yyyyzzz_1, g_0_yyyzz_0_yyyzzzz_0, g_0_yyyzz_0_yyyzzzz_1, g_0_yyyzz_0_yyzzzzz_0, g_0_yyyzz_0_yyzzzzz_1, g_0_yyyzz_0_yzzzzzz_0, g_0_yyyzz_0_yzzzzzz_1, g_0_yyyzz_0_zzzzzzz_0, g_0_yyyzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzz_0_xxxxxxx_0[i] = 2.0 * g_0_xxyzz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxxxx_0[i] * pb_y + g_0_xxyyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxxxxy_0[i] = g_0_xxyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxxxxy_0[i] * pb_z + g_0_xxyyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxxxxz_0[i] = 2.0 * g_0_xxyzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxxxz_0[i] * pb_y + g_0_xxyyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxxxyy_0[i] = g_0_xxyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxxxyy_0[i] * pb_z + g_0_xxyyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxxxyz_0[i] = g_0_yyyzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxxxyz_0[i] * pb_x + g_0_xyyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxxxzz_0[i] = 2.0 * g_0_xxyzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxxzz_0[i] * pb_y + g_0_xxyyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxxyyy_0[i] = g_0_xxyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxxyyy_0[i] * pb_z + g_0_xxyyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxxyyz_0[i] = g_0_yyyzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxxyyz_0[i] * pb_x + g_0_xyyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxxyzz_0[i] = g_0_yyyzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxxyzz_0[i] * pb_x + g_0_xyyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxxzzz_0[i] = 2.0 * g_0_xxyzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxzzz_0[i] * pb_y + g_0_xxyyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxyyyy_0[i] = g_0_xxyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxyyyy_0[i] * pb_z + g_0_xxyyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxyyyz_0[i] = g_0_yyyzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxyyyz_0[i] * pb_x + g_0_xyyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxyyzz_0[i] = g_0_yyyzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxyyzz_0[i] * pb_x + g_0_xyyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxyzzz_0[i] = g_0_yyyzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxyzzz_0[i] * pb_x + g_0_xyyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxzzzz_0[i] = 2.0 * g_0_xxyzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxzzzz_0[i] * pb_y + g_0_xxyyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxyyyyy_0[i] = g_0_xxyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxyyyyy_0[i] * pb_z + g_0_xxyyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxyyyyz_0[i] = g_0_yyyzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyyyyz_0[i] * pb_x + g_0_xyyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxyyyzz_0[i] = g_0_yyyzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyyyzz_0[i] * pb_x + g_0_xyyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxyyzzz_0[i] = g_0_yyyzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyyzzz_0[i] * pb_x + g_0_xyyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxyzzzz_0[i] = g_0_yyyzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyzzzz_0[i] * pb_x + g_0_xyyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxzzzzz_0[i] = 2.0 * g_0_xxyzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxzzzzz_0[i] * pb_y + g_0_xxyyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xyyyyyy_0[i] = g_0_xxyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xyyyyyy_0[i] * pb_z + g_0_xxyyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xyyyyyz_0[i] = g_0_yyyzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xyyyyyz_0[i] * pb_x + g_0_xyyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyyyyzz_0[i] = g_0_yyyzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xyyyyzz_0[i] * pb_x + g_0_xyyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyyyzzz_0[i] = g_0_yyyzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xyyyzzz_0[i] * pb_x + g_0_xyyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyyzzzz_0[i] = g_0_yyyzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xyyzzzz_0[i] * pb_x + g_0_xyyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyzzzzz_0[i] = g_0_yyyzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xyzzzzz_0[i] * pb_x + g_0_xyyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xzzzzzz_0[i] = 2.0 * g_0_xxyzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyyzz_0_xzzzzzz_0[i] * pb_y + g_0_xxyyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_yyyyyyy_0[i] = g_0_yyyzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyyyy_0[i] * pb_x + g_0_xyyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyyyyz_0[i] = g_0_yyyzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyyyz_0[i] * pb_x + g_0_xyyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyyyzz_0[i] = g_0_yyyzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyyzz_0[i] * pb_x + g_0_xyyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyyzzz_0[i] = g_0_yyyzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyzzz_0[i] * pb_x + g_0_xyyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyzzzz_0[i] = g_0_yyyzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyzzzz_0[i] * pb_x + g_0_xyyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyzzzzz_0[i] = g_0_yyyzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyzzzzz_0[i] * pb_x + g_0_xyyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yzzzzzz_0[i] = g_0_yyyzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzzzzzz_0[i] * pb_x + g_0_xyyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_zzzzzzz_0[i] = g_0_yyyzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_zzzzzzz_0[i] * pb_x + g_0_xyyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 648-684 components of targeted buffer : SKSK

    auto g_0_xxyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 648);

    auto g_0_xxyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 649);

    auto g_0_xxyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 650);

    auto g_0_xxyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 651);

    auto g_0_xxyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 652);

    auto g_0_xxyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 653);

    auto g_0_xxyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 654);

    auto g_0_xxyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 655);

    auto g_0_xxyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 656);

    auto g_0_xxyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 657);

    auto g_0_xxyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 658);

    auto g_0_xxyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 659);

    auto g_0_xxyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 660);

    auto g_0_xxyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 661);

    auto g_0_xxyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 662);

    auto g_0_xxyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 663);

    auto g_0_xxyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 664);

    auto g_0_xxyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 665);

    auto g_0_xxyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 666);

    auto g_0_xxyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 667);

    auto g_0_xxyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 668);

    auto g_0_xxyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 669);

    auto g_0_xxyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 670);

    auto g_0_xxyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 671);

    auto g_0_xxyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 672);

    auto g_0_xxyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 673);

    auto g_0_xxyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 674);

    auto g_0_xxyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 675);

    auto g_0_xxyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 676);

    auto g_0_xxyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 677);

    auto g_0_xxyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 678);

    auto g_0_xxyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 679);

    auto g_0_xxyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 680);

    auto g_0_xxyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 681);

    auto g_0_xxyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 682);

    auto g_0_xxyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 683);

    #pragma omp simd aligned(g_0_xxyyz_0_xxxxxxy_0, g_0_xxyyz_0_xxxxxxy_1, g_0_xxyyz_0_xxxxxyy_0, g_0_xxyyz_0_xxxxxyy_1, g_0_xxyyz_0_xxxxyyy_0, g_0_xxyyz_0_xxxxyyy_1, g_0_xxyyz_0_xxxyyyy_0, g_0_xxyyz_0_xxxyyyy_1, g_0_xxyyz_0_xxyyyyy_0, g_0_xxyyz_0_xxyyyyy_1, g_0_xxyyz_0_xyyyyyy_0, g_0_xxyyz_0_xyyyyyy_1, g_0_xxyyzz_0_xxxxxxy_0, g_0_xxyyzz_0_xxxxxxy_1, g_0_xxyyzz_0_xxxxxyy_0, g_0_xxyyzz_0_xxxxxyy_1, g_0_xxyyzz_0_xxxxyyy_0, g_0_xxyyzz_0_xxxxyyy_1, g_0_xxyyzz_0_xxxyyyy_0, g_0_xxyyzz_0_xxxyyyy_1, g_0_xxyyzz_0_xxyyyyy_0, g_0_xxyyzz_0_xxyyyyy_1, g_0_xxyyzz_0_xyyyyyy_0, g_0_xxyyzz_0_xyyyyyy_1, g_0_xxyyzzz_0_xxxxxxx_0, g_0_xxyyzzz_0_xxxxxxy_0, g_0_xxyyzzz_0_xxxxxxz_0, g_0_xxyyzzz_0_xxxxxyy_0, g_0_xxyyzzz_0_xxxxxyz_0, g_0_xxyyzzz_0_xxxxxzz_0, g_0_xxyyzzz_0_xxxxyyy_0, g_0_xxyyzzz_0_xxxxyyz_0, g_0_xxyyzzz_0_xxxxyzz_0, g_0_xxyyzzz_0_xxxxzzz_0, g_0_xxyyzzz_0_xxxyyyy_0, g_0_xxyyzzz_0_xxxyyyz_0, g_0_xxyyzzz_0_xxxyyzz_0, g_0_xxyyzzz_0_xxxyzzz_0, g_0_xxyyzzz_0_xxxzzzz_0, g_0_xxyyzzz_0_xxyyyyy_0, g_0_xxyyzzz_0_xxyyyyz_0, g_0_xxyyzzz_0_xxyyyzz_0, g_0_xxyyzzz_0_xxyyzzz_0, g_0_xxyyzzz_0_xxyzzzz_0, g_0_xxyyzzz_0_xxzzzzz_0, g_0_xxyyzzz_0_xyyyyyy_0, g_0_xxyyzzz_0_xyyyyyz_0, g_0_xxyyzzz_0_xyyyyzz_0, g_0_xxyyzzz_0_xyyyzzz_0, g_0_xxyyzzz_0_xyyzzzz_0, g_0_xxyyzzz_0_xyzzzzz_0, g_0_xxyyzzz_0_xzzzzzz_0, g_0_xxyyzzz_0_yyyyyyy_0, g_0_xxyyzzz_0_yyyyyyz_0, g_0_xxyyzzz_0_yyyyyzz_0, g_0_xxyyzzz_0_yyyyzzz_0, g_0_xxyyzzz_0_yyyzzzz_0, g_0_xxyyzzz_0_yyzzzzz_0, g_0_xxyyzzz_0_yzzzzzz_0, g_0_xxyyzzz_0_zzzzzzz_0, g_0_xxyzzz_0_xxxxxxx_0, g_0_xxyzzz_0_xxxxxxx_1, g_0_xxyzzz_0_xxxxxxz_0, g_0_xxyzzz_0_xxxxxxz_1, g_0_xxyzzz_0_xxxxxzz_0, g_0_xxyzzz_0_xxxxxzz_1, g_0_xxyzzz_0_xxxxzzz_0, g_0_xxyzzz_0_xxxxzzz_1, g_0_xxyzzz_0_xxxzzzz_0, g_0_xxyzzz_0_xxxzzzz_1, g_0_xxyzzz_0_xxzzzzz_0, g_0_xxyzzz_0_xxzzzzz_1, g_0_xxyzzz_0_xzzzzzz_0, g_0_xxyzzz_0_xzzzzzz_1, g_0_xxzzz_0_xxxxxxx_0, g_0_xxzzz_0_xxxxxxx_1, g_0_xxzzz_0_xxxxxxz_0, g_0_xxzzz_0_xxxxxxz_1, g_0_xxzzz_0_xxxxxzz_0, g_0_xxzzz_0_xxxxxzz_1, g_0_xxzzz_0_xxxxzzz_0, g_0_xxzzz_0_xxxxzzz_1, g_0_xxzzz_0_xxxzzzz_0, g_0_xxzzz_0_xxxzzzz_1, g_0_xxzzz_0_xxzzzzz_0, g_0_xxzzz_0_xxzzzzz_1, g_0_xxzzz_0_xzzzzzz_0, g_0_xxzzz_0_xzzzzzz_1, g_0_xyyzzz_0_xxxxxyz_0, g_0_xyyzzz_0_xxxxxyz_1, g_0_xyyzzz_0_xxxxyyz_0, g_0_xyyzzz_0_xxxxyyz_1, g_0_xyyzzz_0_xxxxyz_1, g_0_xyyzzz_0_xxxxyzz_0, g_0_xyyzzz_0_xxxxyzz_1, g_0_xyyzzz_0_xxxyyyz_0, g_0_xyyzzz_0_xxxyyyz_1, g_0_xyyzzz_0_xxxyyz_1, g_0_xyyzzz_0_xxxyyzz_0, g_0_xyyzzz_0_xxxyyzz_1, g_0_xyyzzz_0_xxxyzz_1, g_0_xyyzzz_0_xxxyzzz_0, g_0_xyyzzz_0_xxxyzzz_1, g_0_xyyzzz_0_xxyyyyz_0, g_0_xyyzzz_0_xxyyyyz_1, g_0_xyyzzz_0_xxyyyz_1, g_0_xyyzzz_0_xxyyyzz_0, g_0_xyyzzz_0_xxyyyzz_1, g_0_xyyzzz_0_xxyyzz_1, g_0_xyyzzz_0_xxyyzzz_0, g_0_xyyzzz_0_xxyyzzz_1, g_0_xyyzzz_0_xxyzzz_1, g_0_xyyzzz_0_xxyzzzz_0, g_0_xyyzzz_0_xxyzzzz_1, g_0_xyyzzz_0_xyyyyyz_0, g_0_xyyzzz_0_xyyyyyz_1, g_0_xyyzzz_0_xyyyyz_1, g_0_xyyzzz_0_xyyyyzz_0, g_0_xyyzzz_0_xyyyyzz_1, g_0_xyyzzz_0_xyyyzz_1, g_0_xyyzzz_0_xyyyzzz_0, g_0_xyyzzz_0_xyyyzzz_1, g_0_xyyzzz_0_xyyzzz_1, g_0_xyyzzz_0_xyyzzzz_0, g_0_xyyzzz_0_xyyzzzz_1, g_0_xyyzzz_0_xyzzzz_1, g_0_xyyzzz_0_xyzzzzz_0, g_0_xyyzzz_0_xyzzzzz_1, g_0_xyyzzz_0_yyyyyyy_0, g_0_xyyzzz_0_yyyyyyy_1, g_0_xyyzzz_0_yyyyyyz_0, g_0_xyyzzz_0_yyyyyyz_1, g_0_xyyzzz_0_yyyyyz_1, g_0_xyyzzz_0_yyyyyzz_0, g_0_xyyzzz_0_yyyyyzz_1, g_0_xyyzzz_0_yyyyzz_1, g_0_xyyzzz_0_yyyyzzz_0, g_0_xyyzzz_0_yyyyzzz_1, g_0_xyyzzz_0_yyyzzz_1, g_0_xyyzzz_0_yyyzzzz_0, g_0_xyyzzz_0_yyyzzzz_1, g_0_xyyzzz_0_yyzzzz_1, g_0_xyyzzz_0_yyzzzzz_0, g_0_xyyzzz_0_yyzzzzz_1, g_0_xyyzzz_0_yzzzzz_1, g_0_xyyzzz_0_yzzzzzz_0, g_0_xyyzzz_0_yzzzzzz_1, g_0_xyyzzz_0_zzzzzzz_0, g_0_xyyzzz_0_zzzzzzz_1, g_0_yyzzz_0_xxxxxyz_0, g_0_yyzzz_0_xxxxxyz_1, g_0_yyzzz_0_xxxxyyz_0, g_0_yyzzz_0_xxxxyyz_1, g_0_yyzzz_0_xxxxyzz_0, g_0_yyzzz_0_xxxxyzz_1, g_0_yyzzz_0_xxxyyyz_0, g_0_yyzzz_0_xxxyyyz_1, g_0_yyzzz_0_xxxyyzz_0, g_0_yyzzz_0_xxxyyzz_1, g_0_yyzzz_0_xxxyzzz_0, g_0_yyzzz_0_xxxyzzz_1, g_0_yyzzz_0_xxyyyyz_0, g_0_yyzzz_0_xxyyyyz_1, g_0_yyzzz_0_xxyyyzz_0, g_0_yyzzz_0_xxyyyzz_1, g_0_yyzzz_0_xxyyzzz_0, g_0_yyzzz_0_xxyyzzz_1, g_0_yyzzz_0_xxyzzzz_0, g_0_yyzzz_0_xxyzzzz_1, g_0_yyzzz_0_xyyyyyz_0, g_0_yyzzz_0_xyyyyyz_1, g_0_yyzzz_0_xyyyyzz_0, g_0_yyzzz_0_xyyyyzz_1, g_0_yyzzz_0_xyyyzzz_0, g_0_yyzzz_0_xyyyzzz_1, g_0_yyzzz_0_xyyzzzz_0, g_0_yyzzz_0_xyyzzzz_1, g_0_yyzzz_0_xyzzzzz_0, g_0_yyzzz_0_xyzzzzz_1, g_0_yyzzz_0_yyyyyyy_0, g_0_yyzzz_0_yyyyyyy_1, g_0_yyzzz_0_yyyyyyz_0, g_0_yyzzz_0_yyyyyyz_1, g_0_yyzzz_0_yyyyyzz_0, g_0_yyzzz_0_yyyyyzz_1, g_0_yyzzz_0_yyyyzzz_0, g_0_yyzzz_0_yyyyzzz_1, g_0_yyzzz_0_yyyzzzz_0, g_0_yyzzz_0_yyyzzzz_1, g_0_yyzzz_0_yyzzzzz_0, g_0_yyzzz_0_yyzzzzz_1, g_0_yyzzz_0_yzzzzzz_0, g_0_yyzzz_0_yzzzzzz_1, g_0_yyzzz_0_zzzzzzz_0, g_0_yyzzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzz_0_xxxxxxx_0[i] = g_0_xxzzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxyzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxxxxy_0[i] = 2.0 * g_0_xxyyz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxxxy_0[i] * pb_z + g_0_xxyyzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxxxxz_0[i] = g_0_xxzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxyzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxyyz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxxyy_0[i] * pb_z + g_0_xxyyzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxxxyz_0[i] = g_0_yyzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxxxyz_0[i] * pb_x + g_0_xyyzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxxxzz_0[i] = g_0_xxzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxyzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxxyyy_0[i] = 2.0 * g_0_xxyyz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxxyyy_0[i] * pb_z + g_0_xxyyzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxxyyz_0[i] = g_0_yyzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxxyyz_0[i] * pb_x + g_0_xyyzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxxyzz_0[i] = g_0_yyzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxxyzz_0[i] * pb_x + g_0_xyyzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxxzzz_0[i] = g_0_xxzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxyzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxyyyy_0[i] = 2.0 * g_0_xxyyz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxxyyyy_0[i] * pb_z + g_0_xxyyzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxyyyz_0[i] = g_0_yyzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxyyyz_0[i] * pb_x + g_0_xyyzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxyyzz_0[i] = g_0_yyzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxyyzz_0[i] * pb_x + g_0_xyyzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxyzzz_0[i] = g_0_yyzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxyzzz_0[i] * pb_x + g_0_xyyzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxzzzz_0[i] = g_0_xxzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxyzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_xxyyz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xxyyyyy_0[i] * pb_z + g_0_xxyyzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxyyyyz_0[i] = g_0_yyzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyyyyz_0[i] * pb_x + g_0_xyyzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxyyyzz_0[i] = g_0_yyzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyyyzz_0[i] * pb_x + g_0_xyyzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxyyzzz_0[i] = g_0_yyzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyyzzz_0[i] * pb_x + g_0_xyyzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxyzzzz_0[i] = g_0_yyzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyzzzz_0[i] * pb_x + g_0_xyyzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxzzzzz_0[i] = g_0_xxzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxyzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xyyyyyy_0[i] = 2.0 * g_0_xxyyz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyzz_0_xyyyyyy_0[i] * pb_z + g_0_xxyyzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xyyyyyz_0[i] = g_0_yyzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xyyyyyz_0[i] * pb_x + g_0_xyyzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyyyyzz_0[i] = g_0_yyzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xyyyyzz_0[i] * pb_x + g_0_xyyzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyyyzzz_0[i] = g_0_yyzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xyyyzzz_0[i] * pb_x + g_0_xyyzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyyzzzz_0[i] = g_0_yyzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xyyzzzz_0[i] * pb_x + g_0_xyyzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyzzzzz_0[i] = g_0_yyzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xyzzzzz_0[i] * pb_x + g_0_xyyzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xzzzzzz_0[i] = g_0_xxzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxyzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_yyyyyyy_0[i] = g_0_yyzzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyyyy_0[i] * pb_x + g_0_xyyzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyyyyz_0[i] = g_0_yyzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyyyz_0[i] * pb_x + g_0_xyyzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyyyzz_0[i] = g_0_yyzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyyzz_0[i] * pb_x + g_0_xyyzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyyzzz_0[i] = g_0_yyzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyzzz_0[i] * pb_x + g_0_xyyzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyzzzz_0[i] = g_0_yyzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyzzzz_0[i] * pb_x + g_0_xyyzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyzzzzz_0[i] = g_0_yyzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyzzzzz_0[i] * pb_x + g_0_xyyzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yzzzzzz_0[i] = g_0_yyzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzzzzzz_0[i] * pb_x + g_0_xyyzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_zzzzzzz_0[i] = g_0_yyzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_zzzzzzz_0[i] * pb_x + g_0_xyyzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 684-720 components of targeted buffer : SKSK

    auto g_0_xxyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 684);

    auto g_0_xxyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 685);

    auto g_0_xxyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 686);

    auto g_0_xxyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 687);

    auto g_0_xxyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 688);

    auto g_0_xxyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 689);

    auto g_0_xxyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 690);

    auto g_0_xxyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 691);

    auto g_0_xxyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 692);

    auto g_0_xxyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 693);

    auto g_0_xxyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 694);

    auto g_0_xxyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 695);

    auto g_0_xxyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 696);

    auto g_0_xxyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 697);

    auto g_0_xxyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 698);

    auto g_0_xxyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 699);

    auto g_0_xxyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 700);

    auto g_0_xxyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 701);

    auto g_0_xxyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 702);

    auto g_0_xxyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 703);

    auto g_0_xxyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 704);

    auto g_0_xxyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 705);

    auto g_0_xxyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 706);

    auto g_0_xxyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 707);

    auto g_0_xxyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 708);

    auto g_0_xxyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 709);

    auto g_0_xxyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 710);

    auto g_0_xxyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 711);

    auto g_0_xxyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 712);

    auto g_0_xxyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 713);

    auto g_0_xxyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 714);

    auto g_0_xxyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 715);

    auto g_0_xxyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 716);

    auto g_0_xxyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 717);

    auto g_0_xxyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 718);

    auto g_0_xxyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 719);

    #pragma omp simd aligned(g_0_xxyzzzz_0_xxxxxxx_0, g_0_xxyzzzz_0_xxxxxxy_0, g_0_xxyzzzz_0_xxxxxxz_0, g_0_xxyzzzz_0_xxxxxyy_0, g_0_xxyzzzz_0_xxxxxyz_0, g_0_xxyzzzz_0_xxxxxzz_0, g_0_xxyzzzz_0_xxxxyyy_0, g_0_xxyzzzz_0_xxxxyyz_0, g_0_xxyzzzz_0_xxxxyzz_0, g_0_xxyzzzz_0_xxxxzzz_0, g_0_xxyzzzz_0_xxxyyyy_0, g_0_xxyzzzz_0_xxxyyyz_0, g_0_xxyzzzz_0_xxxyyzz_0, g_0_xxyzzzz_0_xxxyzzz_0, g_0_xxyzzzz_0_xxxzzzz_0, g_0_xxyzzzz_0_xxyyyyy_0, g_0_xxyzzzz_0_xxyyyyz_0, g_0_xxyzzzz_0_xxyyyzz_0, g_0_xxyzzzz_0_xxyyzzz_0, g_0_xxyzzzz_0_xxyzzzz_0, g_0_xxyzzzz_0_xxzzzzz_0, g_0_xxyzzzz_0_xyyyyyy_0, g_0_xxyzzzz_0_xyyyyyz_0, g_0_xxyzzzz_0_xyyyyzz_0, g_0_xxyzzzz_0_xyyyzzz_0, g_0_xxyzzzz_0_xyyzzzz_0, g_0_xxyzzzz_0_xyzzzzz_0, g_0_xxyzzzz_0_xzzzzzz_0, g_0_xxyzzzz_0_yyyyyyy_0, g_0_xxyzzzz_0_yyyyyyz_0, g_0_xxyzzzz_0_yyyyyzz_0, g_0_xxyzzzz_0_yyyyzzz_0, g_0_xxyzzzz_0_yyyzzzz_0, g_0_xxyzzzz_0_yyzzzzz_0, g_0_xxyzzzz_0_yzzzzzz_0, g_0_xxyzzzz_0_zzzzzzz_0, g_0_xxzzzz_0_xxxxxx_1, g_0_xxzzzz_0_xxxxxxx_0, g_0_xxzzzz_0_xxxxxxx_1, g_0_xxzzzz_0_xxxxxxy_0, g_0_xxzzzz_0_xxxxxxy_1, g_0_xxzzzz_0_xxxxxxz_0, g_0_xxzzzz_0_xxxxxxz_1, g_0_xxzzzz_0_xxxxxy_1, g_0_xxzzzz_0_xxxxxyy_0, g_0_xxzzzz_0_xxxxxyy_1, g_0_xxzzzz_0_xxxxxyz_0, g_0_xxzzzz_0_xxxxxyz_1, g_0_xxzzzz_0_xxxxxz_1, g_0_xxzzzz_0_xxxxxzz_0, g_0_xxzzzz_0_xxxxxzz_1, g_0_xxzzzz_0_xxxxyy_1, g_0_xxzzzz_0_xxxxyyy_0, g_0_xxzzzz_0_xxxxyyy_1, g_0_xxzzzz_0_xxxxyyz_0, g_0_xxzzzz_0_xxxxyyz_1, g_0_xxzzzz_0_xxxxyz_1, g_0_xxzzzz_0_xxxxyzz_0, g_0_xxzzzz_0_xxxxyzz_1, g_0_xxzzzz_0_xxxxzz_1, g_0_xxzzzz_0_xxxxzzz_0, g_0_xxzzzz_0_xxxxzzz_1, g_0_xxzzzz_0_xxxyyy_1, g_0_xxzzzz_0_xxxyyyy_0, g_0_xxzzzz_0_xxxyyyy_1, g_0_xxzzzz_0_xxxyyyz_0, g_0_xxzzzz_0_xxxyyyz_1, g_0_xxzzzz_0_xxxyyz_1, g_0_xxzzzz_0_xxxyyzz_0, g_0_xxzzzz_0_xxxyyzz_1, g_0_xxzzzz_0_xxxyzz_1, g_0_xxzzzz_0_xxxyzzz_0, g_0_xxzzzz_0_xxxyzzz_1, g_0_xxzzzz_0_xxxzzz_1, g_0_xxzzzz_0_xxxzzzz_0, g_0_xxzzzz_0_xxxzzzz_1, g_0_xxzzzz_0_xxyyyy_1, g_0_xxzzzz_0_xxyyyyy_0, g_0_xxzzzz_0_xxyyyyy_1, g_0_xxzzzz_0_xxyyyyz_0, g_0_xxzzzz_0_xxyyyyz_1, g_0_xxzzzz_0_xxyyyz_1, g_0_xxzzzz_0_xxyyyzz_0, g_0_xxzzzz_0_xxyyyzz_1, g_0_xxzzzz_0_xxyyzz_1, g_0_xxzzzz_0_xxyyzzz_0, g_0_xxzzzz_0_xxyyzzz_1, g_0_xxzzzz_0_xxyzzz_1, g_0_xxzzzz_0_xxyzzzz_0, g_0_xxzzzz_0_xxyzzzz_1, g_0_xxzzzz_0_xxzzzz_1, g_0_xxzzzz_0_xxzzzzz_0, g_0_xxzzzz_0_xxzzzzz_1, g_0_xxzzzz_0_xyyyyy_1, g_0_xxzzzz_0_xyyyyyy_0, g_0_xxzzzz_0_xyyyyyy_1, g_0_xxzzzz_0_xyyyyyz_0, g_0_xxzzzz_0_xyyyyyz_1, g_0_xxzzzz_0_xyyyyz_1, g_0_xxzzzz_0_xyyyyzz_0, g_0_xxzzzz_0_xyyyyzz_1, g_0_xxzzzz_0_xyyyzz_1, g_0_xxzzzz_0_xyyyzzz_0, g_0_xxzzzz_0_xyyyzzz_1, g_0_xxzzzz_0_xyyzzz_1, g_0_xxzzzz_0_xyyzzzz_0, g_0_xxzzzz_0_xyyzzzz_1, g_0_xxzzzz_0_xyzzzz_1, g_0_xxzzzz_0_xyzzzzz_0, g_0_xxzzzz_0_xyzzzzz_1, g_0_xxzzzz_0_xzzzzz_1, g_0_xxzzzz_0_xzzzzzz_0, g_0_xxzzzz_0_xzzzzzz_1, g_0_xxzzzz_0_yyyyyy_1, g_0_xxzzzz_0_yyyyyyy_0, g_0_xxzzzz_0_yyyyyyy_1, g_0_xxzzzz_0_yyyyyyz_0, g_0_xxzzzz_0_yyyyyyz_1, g_0_xxzzzz_0_yyyyyz_1, g_0_xxzzzz_0_yyyyyzz_0, g_0_xxzzzz_0_yyyyyzz_1, g_0_xxzzzz_0_yyyyzz_1, g_0_xxzzzz_0_yyyyzzz_0, g_0_xxzzzz_0_yyyyzzz_1, g_0_xxzzzz_0_yyyzzz_1, g_0_xxzzzz_0_yyyzzzz_0, g_0_xxzzzz_0_yyyzzzz_1, g_0_xxzzzz_0_yyzzzz_1, g_0_xxzzzz_0_yyzzzzz_0, g_0_xxzzzz_0_yyzzzzz_1, g_0_xxzzzz_0_yzzzzz_1, g_0_xxzzzz_0_yzzzzzz_0, g_0_xxzzzz_0_yzzzzzz_1, g_0_xxzzzz_0_zzzzzz_1, g_0_xxzzzz_0_zzzzzzz_0, g_0_xxzzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzz_0_xxxxxxx_0[i] = g_0_xxzzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxxxy_0[i] = g_0_xxzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxxxy_0[i] * pb_y + g_0_xxzzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxxxz_0[i] = g_0_xxzzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxxyy_0[i] * pb_y + g_0_xxzzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxxyz_0[i] = g_0_xxzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxxyz_0[i] * pb_y + g_0_xxzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxxzz_0[i] = g_0_xxzzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxyyy_0[i] * pb_y + g_0_xxzzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxyyz_0[i] * pb_y + g_0_xxzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxyzz_0[i] = g_0_xxzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxyzz_0[i] * pb_y + g_0_xxzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxzzz_0[i] = g_0_xxzzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyyyy_0[i] * pb_y + g_0_xxzzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyyyz_0[i] * pb_y + g_0_xxzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyyzz_0[i] * pb_y + g_0_xxzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyzzz_0[i] = g_0_xxzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyzzz_0[i] * pb_y + g_0_xxzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxzzzz_0[i] = g_0_xxzzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyyyy_0[i] * pb_y + g_0_xxzzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyyyz_0[i] * pb_y + g_0_xxzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyyzz_0[i] * pb_y + g_0_xxzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyzzz_0[i] * pb_y + g_0_xxzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyzzzz_0[i] = g_0_xxzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyzzzz_0[i] * pb_y + g_0_xxzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxzzzzz_0[i] = g_0_xxzzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyyyy_0[i] * pb_y + g_0_xxzzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyyyz_0[i] * pb_y + g_0_xxzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyyzz_0[i] * pb_y + g_0_xxzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyzzz_0[i] * pb_y + g_0_xxzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyzzzz_0[i] * pb_y + g_0_xxzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyzzzzz_0[i] = g_0_xxzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyzzzzz_0[i] * pb_y + g_0_xxzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xzzzzzz_0[i] = g_0_xxzzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyyyy_0[i] * pb_y + g_0_xxzzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyyyz_0[i] * pb_y + g_0_xxzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyyzz_0[i] * pb_y + g_0_xxzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyzzz_0[i] * pb_y + g_0_xxzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyzzzz_0[i] * pb_y + g_0_xxzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyzzzzz_0[i] * pb_y + g_0_xxzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yzzzzzz_0[i] = g_0_xxzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yzzzzzz_0[i] * pb_y + g_0_xxzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_zzzzzzz_0[i] = g_0_xxzzzz_0_zzzzzzz_0[i] * pb_y + g_0_xxzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 720-756 components of targeted buffer : SKSK

    auto g_0_xxzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 720);

    auto g_0_xxzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 721);

    auto g_0_xxzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 722);

    auto g_0_xxzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 723);

    auto g_0_xxzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 724);

    auto g_0_xxzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 725);

    auto g_0_xxzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 726);

    auto g_0_xxzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 727);

    auto g_0_xxzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 728);

    auto g_0_xxzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 729);

    auto g_0_xxzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 730);

    auto g_0_xxzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 731);

    auto g_0_xxzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 732);

    auto g_0_xxzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 733);

    auto g_0_xxzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 734);

    auto g_0_xxzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 735);

    auto g_0_xxzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 736);

    auto g_0_xxzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 737);

    auto g_0_xxzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 738);

    auto g_0_xxzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 739);

    auto g_0_xxzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 740);

    auto g_0_xxzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 741);

    auto g_0_xxzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 742);

    auto g_0_xxzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 743);

    auto g_0_xxzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 744);

    auto g_0_xxzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 745);

    auto g_0_xxzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 746);

    auto g_0_xxzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 747);

    auto g_0_xxzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 748);

    auto g_0_xxzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 749);

    auto g_0_xxzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 750);

    auto g_0_xxzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 751);

    auto g_0_xxzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 752);

    auto g_0_xxzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 753);

    auto g_0_xxzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 754);

    auto g_0_xxzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 755);

    #pragma omp simd aligned(g_0_xxzzz_0_xxxxxxx_0, g_0_xxzzz_0_xxxxxxx_1, g_0_xxzzz_0_xxxxxxy_0, g_0_xxzzz_0_xxxxxxy_1, g_0_xxzzz_0_xxxxxyy_0, g_0_xxzzz_0_xxxxxyy_1, g_0_xxzzz_0_xxxxyyy_0, g_0_xxzzz_0_xxxxyyy_1, g_0_xxzzz_0_xxxyyyy_0, g_0_xxzzz_0_xxxyyyy_1, g_0_xxzzz_0_xxyyyyy_0, g_0_xxzzz_0_xxyyyyy_1, g_0_xxzzz_0_xyyyyyy_0, g_0_xxzzz_0_xyyyyyy_1, g_0_xxzzzz_0_xxxxxxx_0, g_0_xxzzzz_0_xxxxxxx_1, g_0_xxzzzz_0_xxxxxxy_0, g_0_xxzzzz_0_xxxxxxy_1, g_0_xxzzzz_0_xxxxxyy_0, g_0_xxzzzz_0_xxxxxyy_1, g_0_xxzzzz_0_xxxxyyy_0, g_0_xxzzzz_0_xxxxyyy_1, g_0_xxzzzz_0_xxxyyyy_0, g_0_xxzzzz_0_xxxyyyy_1, g_0_xxzzzz_0_xxyyyyy_0, g_0_xxzzzz_0_xxyyyyy_1, g_0_xxzzzz_0_xyyyyyy_0, g_0_xxzzzz_0_xyyyyyy_1, g_0_xxzzzzz_0_xxxxxxx_0, g_0_xxzzzzz_0_xxxxxxy_0, g_0_xxzzzzz_0_xxxxxxz_0, g_0_xxzzzzz_0_xxxxxyy_0, g_0_xxzzzzz_0_xxxxxyz_0, g_0_xxzzzzz_0_xxxxxzz_0, g_0_xxzzzzz_0_xxxxyyy_0, g_0_xxzzzzz_0_xxxxyyz_0, g_0_xxzzzzz_0_xxxxyzz_0, g_0_xxzzzzz_0_xxxxzzz_0, g_0_xxzzzzz_0_xxxyyyy_0, g_0_xxzzzzz_0_xxxyyyz_0, g_0_xxzzzzz_0_xxxyyzz_0, g_0_xxzzzzz_0_xxxyzzz_0, g_0_xxzzzzz_0_xxxzzzz_0, g_0_xxzzzzz_0_xxyyyyy_0, g_0_xxzzzzz_0_xxyyyyz_0, g_0_xxzzzzz_0_xxyyyzz_0, g_0_xxzzzzz_0_xxyyzzz_0, g_0_xxzzzzz_0_xxyzzzz_0, g_0_xxzzzzz_0_xxzzzzz_0, g_0_xxzzzzz_0_xyyyyyy_0, g_0_xxzzzzz_0_xyyyyyz_0, g_0_xxzzzzz_0_xyyyyzz_0, g_0_xxzzzzz_0_xyyyzzz_0, g_0_xxzzzzz_0_xyyzzzz_0, g_0_xxzzzzz_0_xyzzzzz_0, g_0_xxzzzzz_0_xzzzzzz_0, g_0_xxzzzzz_0_yyyyyyy_0, g_0_xxzzzzz_0_yyyyyyz_0, g_0_xxzzzzz_0_yyyyyzz_0, g_0_xxzzzzz_0_yyyyzzz_0, g_0_xxzzzzz_0_yyyzzzz_0, g_0_xxzzzzz_0_yyzzzzz_0, g_0_xxzzzzz_0_yzzzzzz_0, g_0_xxzzzzz_0_zzzzzzz_0, g_0_xzzzzz_0_xxxxxxz_0, g_0_xzzzzz_0_xxxxxxz_1, g_0_xzzzzz_0_xxxxxyz_0, g_0_xzzzzz_0_xxxxxyz_1, g_0_xzzzzz_0_xxxxxz_1, g_0_xzzzzz_0_xxxxxzz_0, g_0_xzzzzz_0_xxxxxzz_1, g_0_xzzzzz_0_xxxxyyz_0, g_0_xzzzzz_0_xxxxyyz_1, g_0_xzzzzz_0_xxxxyz_1, g_0_xzzzzz_0_xxxxyzz_0, g_0_xzzzzz_0_xxxxyzz_1, g_0_xzzzzz_0_xxxxzz_1, g_0_xzzzzz_0_xxxxzzz_0, g_0_xzzzzz_0_xxxxzzz_1, g_0_xzzzzz_0_xxxyyyz_0, g_0_xzzzzz_0_xxxyyyz_1, g_0_xzzzzz_0_xxxyyz_1, g_0_xzzzzz_0_xxxyyzz_0, g_0_xzzzzz_0_xxxyyzz_1, g_0_xzzzzz_0_xxxyzz_1, g_0_xzzzzz_0_xxxyzzz_0, g_0_xzzzzz_0_xxxyzzz_1, g_0_xzzzzz_0_xxxzzz_1, g_0_xzzzzz_0_xxxzzzz_0, g_0_xzzzzz_0_xxxzzzz_1, g_0_xzzzzz_0_xxyyyyz_0, g_0_xzzzzz_0_xxyyyyz_1, g_0_xzzzzz_0_xxyyyz_1, g_0_xzzzzz_0_xxyyyzz_0, g_0_xzzzzz_0_xxyyyzz_1, g_0_xzzzzz_0_xxyyzz_1, g_0_xzzzzz_0_xxyyzzz_0, g_0_xzzzzz_0_xxyyzzz_1, g_0_xzzzzz_0_xxyzzz_1, g_0_xzzzzz_0_xxyzzzz_0, g_0_xzzzzz_0_xxyzzzz_1, g_0_xzzzzz_0_xxzzzz_1, g_0_xzzzzz_0_xxzzzzz_0, g_0_xzzzzz_0_xxzzzzz_1, g_0_xzzzzz_0_xyyyyyz_0, g_0_xzzzzz_0_xyyyyyz_1, g_0_xzzzzz_0_xyyyyz_1, g_0_xzzzzz_0_xyyyyzz_0, g_0_xzzzzz_0_xyyyyzz_1, g_0_xzzzzz_0_xyyyzz_1, g_0_xzzzzz_0_xyyyzzz_0, g_0_xzzzzz_0_xyyyzzz_1, g_0_xzzzzz_0_xyyzzz_1, g_0_xzzzzz_0_xyyzzzz_0, g_0_xzzzzz_0_xyyzzzz_1, g_0_xzzzzz_0_xyzzzz_1, g_0_xzzzzz_0_xyzzzzz_0, g_0_xzzzzz_0_xyzzzzz_1, g_0_xzzzzz_0_xzzzzz_1, g_0_xzzzzz_0_xzzzzzz_0, g_0_xzzzzz_0_xzzzzzz_1, g_0_xzzzzz_0_yyyyyyy_0, g_0_xzzzzz_0_yyyyyyy_1, g_0_xzzzzz_0_yyyyyyz_0, g_0_xzzzzz_0_yyyyyyz_1, g_0_xzzzzz_0_yyyyyz_1, g_0_xzzzzz_0_yyyyyzz_0, g_0_xzzzzz_0_yyyyyzz_1, g_0_xzzzzz_0_yyyyzz_1, g_0_xzzzzz_0_yyyyzzz_0, g_0_xzzzzz_0_yyyyzzz_1, g_0_xzzzzz_0_yyyzzz_1, g_0_xzzzzz_0_yyyzzzz_0, g_0_xzzzzz_0_yyyzzzz_1, g_0_xzzzzz_0_yyzzzz_1, g_0_xzzzzz_0_yyzzzzz_0, g_0_xzzzzz_0_yyzzzzz_1, g_0_xzzzzz_0_yzzzzz_1, g_0_xzzzzz_0_yzzzzzz_0, g_0_xzzzzz_0_yzzzzzz_1, g_0_xzzzzz_0_zzzzzz_1, g_0_xzzzzz_0_zzzzzzz_0, g_0_xzzzzz_0_zzzzzzz_1, g_0_zzzzz_0_xxxxxxz_0, g_0_zzzzz_0_xxxxxxz_1, g_0_zzzzz_0_xxxxxyz_0, g_0_zzzzz_0_xxxxxyz_1, g_0_zzzzz_0_xxxxxzz_0, g_0_zzzzz_0_xxxxxzz_1, g_0_zzzzz_0_xxxxyyz_0, g_0_zzzzz_0_xxxxyyz_1, g_0_zzzzz_0_xxxxyzz_0, g_0_zzzzz_0_xxxxyzz_1, g_0_zzzzz_0_xxxxzzz_0, g_0_zzzzz_0_xxxxzzz_1, g_0_zzzzz_0_xxxyyyz_0, g_0_zzzzz_0_xxxyyyz_1, g_0_zzzzz_0_xxxyyzz_0, g_0_zzzzz_0_xxxyyzz_1, g_0_zzzzz_0_xxxyzzz_0, g_0_zzzzz_0_xxxyzzz_1, g_0_zzzzz_0_xxxzzzz_0, g_0_zzzzz_0_xxxzzzz_1, g_0_zzzzz_0_xxyyyyz_0, g_0_zzzzz_0_xxyyyyz_1, g_0_zzzzz_0_xxyyyzz_0, g_0_zzzzz_0_xxyyyzz_1, g_0_zzzzz_0_xxyyzzz_0, g_0_zzzzz_0_xxyyzzz_1, g_0_zzzzz_0_xxyzzzz_0, g_0_zzzzz_0_xxyzzzz_1, g_0_zzzzz_0_xxzzzzz_0, g_0_zzzzz_0_xxzzzzz_1, g_0_zzzzz_0_xyyyyyz_0, g_0_zzzzz_0_xyyyyyz_1, g_0_zzzzz_0_xyyyyzz_0, g_0_zzzzz_0_xyyyyzz_1, g_0_zzzzz_0_xyyyzzz_0, g_0_zzzzz_0_xyyyzzz_1, g_0_zzzzz_0_xyyzzzz_0, g_0_zzzzz_0_xyyzzzz_1, g_0_zzzzz_0_xyzzzzz_0, g_0_zzzzz_0_xyzzzzz_1, g_0_zzzzz_0_xzzzzzz_0, g_0_zzzzz_0_xzzzzzz_1, g_0_zzzzz_0_yyyyyyy_0, g_0_zzzzz_0_yyyyyyy_1, g_0_zzzzz_0_yyyyyyz_0, g_0_zzzzz_0_yyyyyyz_1, g_0_zzzzz_0_yyyyyzz_0, g_0_zzzzz_0_yyyyyzz_1, g_0_zzzzz_0_yyyyzzz_0, g_0_zzzzz_0_yyyyzzz_1, g_0_zzzzz_0_yyyzzzz_0, g_0_zzzzz_0_yyyzzzz_1, g_0_zzzzz_0_yyzzzzz_0, g_0_zzzzz_0_yyzzzzz_1, g_0_zzzzz_0_yzzzzzz_0, g_0_zzzzz_0_yzzzzzz_1, g_0_zzzzz_0_zzzzzzz_0, g_0_zzzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzz_0_xxxxxxx_0[i] = 4.0 * g_0_xxzzz_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxxxxx_0[i] * pb_z + g_0_xxzzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxxxy_0[i] = 4.0 * g_0_xxzzz_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxxxxy_0[i] * pb_z + g_0_xxzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxxxz_0[i] = g_0_zzzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxxxz_0[i] * pb_x + g_0_xzzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxxxyy_0[i] = 4.0 * g_0_xxzzz_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxxxyy_0[i] * pb_z + g_0_xxzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxxyz_0[i] = g_0_zzzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxxyz_0[i] * pb_x + g_0_xzzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxxxzz_0[i] = g_0_zzzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxxzz_0[i] * pb_x + g_0_xzzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_xxzzz_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxxyyy_0[i] * pb_z + g_0_xxzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxyyz_0[i] = g_0_zzzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxyyz_0[i] * pb_x + g_0_xzzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxxyzz_0[i] = g_0_zzzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxyzz_0[i] * pb_x + g_0_xzzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxxzzz_0[i] = g_0_zzzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxzzz_0[i] * pb_x + g_0_xzzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxzzz_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxxyyyy_0[i] * pb_z + g_0_xxzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxyyyz_0[i] = g_0_zzzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxyyyz_0[i] * pb_x + g_0_xzzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxyyzz_0[i] = g_0_zzzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxyyzz_0[i] * pb_x + g_0_xzzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxyzzz_0[i] = g_0_zzzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxyzzz_0[i] * pb_x + g_0_xzzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxzzzz_0[i] = g_0_zzzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxzzzz_0[i] * pb_x + g_0_xzzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyyyyy_0[i] = 4.0 * g_0_xxzzz_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xxyyyyy_0[i] * pb_z + g_0_xxzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxyyyyz_0[i] = g_0_zzzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyyyyz_0[i] * pb_x + g_0_xzzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyyyzz_0[i] = g_0_zzzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyyyzz_0[i] * pb_x + g_0_xzzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyyzzz_0[i] = g_0_zzzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyyzzz_0[i] * pb_x + g_0_xzzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyzzzz_0[i] = g_0_zzzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyzzzz_0[i] * pb_x + g_0_xzzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxzzzzz_0[i] = g_0_zzzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxzzzzz_0[i] * pb_x + g_0_xzzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyyyyy_0[i] = 4.0 * g_0_xxzzz_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxzzzz_0_xyyyyyy_0[i] * pb_z + g_0_xxzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xyyyyyz_0[i] = g_0_zzzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xyyyyyz_0[i] * pb_x + g_0_xzzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyyyzz_0[i] = g_0_zzzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xyyyyzz_0[i] * pb_x + g_0_xzzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyyzzz_0[i] = g_0_zzzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xyyyzzz_0[i] * pb_x + g_0_xzzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyzzzz_0[i] = g_0_zzzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xyyzzzz_0[i] * pb_x + g_0_xzzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyzzzzz_0[i] = g_0_zzzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xyzzzzz_0[i] * pb_x + g_0_xzzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xzzzzzz_0[i] = g_0_zzzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xzzzzzz_0[i] * pb_x + g_0_xzzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyyyy_0[i] = g_0_zzzzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyyyy_0[i] * pb_x + g_0_xzzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyyyz_0[i] = g_0_zzzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyyyz_0[i] * pb_x + g_0_xzzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyyzz_0[i] = g_0_zzzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyyzz_0[i] * pb_x + g_0_xzzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyzzz_0[i] = g_0_zzzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyzzz_0[i] * pb_x + g_0_xzzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyzzzz_0[i] = g_0_zzzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyzzzz_0[i] * pb_x + g_0_xzzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyzzzzz_0[i] = g_0_zzzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyzzzzz_0[i] * pb_x + g_0_xzzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yzzzzzz_0[i] = g_0_zzzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzzzzzz_0[i] * pb_x + g_0_xzzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_zzzzzzz_0[i] = g_0_zzzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzzzzzz_0[i] * pb_x + g_0_xzzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 756-792 components of targeted buffer : SKSK

    auto g_0_xyyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 756);

    auto g_0_xyyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 757);

    auto g_0_xyyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 758);

    auto g_0_xyyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 759);

    auto g_0_xyyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 760);

    auto g_0_xyyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 761);

    auto g_0_xyyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 762);

    auto g_0_xyyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 763);

    auto g_0_xyyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 764);

    auto g_0_xyyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 765);

    auto g_0_xyyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 766);

    auto g_0_xyyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 767);

    auto g_0_xyyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 768);

    auto g_0_xyyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 769);

    auto g_0_xyyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 770);

    auto g_0_xyyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 771);

    auto g_0_xyyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 772);

    auto g_0_xyyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 773);

    auto g_0_xyyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 774);

    auto g_0_xyyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 775);

    auto g_0_xyyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 776);

    auto g_0_xyyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 777);

    auto g_0_xyyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 778);

    auto g_0_xyyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 779);

    auto g_0_xyyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 780);

    auto g_0_xyyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 781);

    auto g_0_xyyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 782);

    auto g_0_xyyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 783);

    auto g_0_xyyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 784);

    auto g_0_xyyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 785);

    auto g_0_xyyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 786);

    auto g_0_xyyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 787);

    auto g_0_xyyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 788);

    auto g_0_xyyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 789);

    auto g_0_xyyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 790);

    auto g_0_xyyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 791);

    #pragma omp simd aligned(g_0_xyyyyyy_0_xxxxxxx_0, g_0_xyyyyyy_0_xxxxxxy_0, g_0_xyyyyyy_0_xxxxxxz_0, g_0_xyyyyyy_0_xxxxxyy_0, g_0_xyyyyyy_0_xxxxxyz_0, g_0_xyyyyyy_0_xxxxxzz_0, g_0_xyyyyyy_0_xxxxyyy_0, g_0_xyyyyyy_0_xxxxyyz_0, g_0_xyyyyyy_0_xxxxyzz_0, g_0_xyyyyyy_0_xxxxzzz_0, g_0_xyyyyyy_0_xxxyyyy_0, g_0_xyyyyyy_0_xxxyyyz_0, g_0_xyyyyyy_0_xxxyyzz_0, g_0_xyyyyyy_0_xxxyzzz_0, g_0_xyyyyyy_0_xxxzzzz_0, g_0_xyyyyyy_0_xxyyyyy_0, g_0_xyyyyyy_0_xxyyyyz_0, g_0_xyyyyyy_0_xxyyyzz_0, g_0_xyyyyyy_0_xxyyzzz_0, g_0_xyyyyyy_0_xxyzzzz_0, g_0_xyyyyyy_0_xxzzzzz_0, g_0_xyyyyyy_0_xyyyyyy_0, g_0_xyyyyyy_0_xyyyyyz_0, g_0_xyyyyyy_0_xyyyyzz_0, g_0_xyyyyyy_0_xyyyzzz_0, g_0_xyyyyyy_0_xyyzzzz_0, g_0_xyyyyyy_0_xyzzzzz_0, g_0_xyyyyyy_0_xzzzzzz_0, g_0_xyyyyyy_0_yyyyyyy_0, g_0_xyyyyyy_0_yyyyyyz_0, g_0_xyyyyyy_0_yyyyyzz_0, g_0_xyyyyyy_0_yyyyzzz_0, g_0_xyyyyyy_0_yyyzzzz_0, g_0_xyyyyyy_0_yyzzzzz_0, g_0_xyyyyyy_0_yzzzzzz_0, g_0_xyyyyyy_0_zzzzzzz_0, g_0_yyyyyy_0_xxxxxx_1, g_0_yyyyyy_0_xxxxxxx_0, g_0_yyyyyy_0_xxxxxxx_1, g_0_yyyyyy_0_xxxxxxy_0, g_0_yyyyyy_0_xxxxxxy_1, g_0_yyyyyy_0_xxxxxxz_0, g_0_yyyyyy_0_xxxxxxz_1, g_0_yyyyyy_0_xxxxxy_1, g_0_yyyyyy_0_xxxxxyy_0, g_0_yyyyyy_0_xxxxxyy_1, g_0_yyyyyy_0_xxxxxyz_0, g_0_yyyyyy_0_xxxxxyz_1, g_0_yyyyyy_0_xxxxxz_1, g_0_yyyyyy_0_xxxxxzz_0, g_0_yyyyyy_0_xxxxxzz_1, g_0_yyyyyy_0_xxxxyy_1, g_0_yyyyyy_0_xxxxyyy_0, g_0_yyyyyy_0_xxxxyyy_1, g_0_yyyyyy_0_xxxxyyz_0, g_0_yyyyyy_0_xxxxyyz_1, g_0_yyyyyy_0_xxxxyz_1, g_0_yyyyyy_0_xxxxyzz_0, g_0_yyyyyy_0_xxxxyzz_1, g_0_yyyyyy_0_xxxxzz_1, g_0_yyyyyy_0_xxxxzzz_0, g_0_yyyyyy_0_xxxxzzz_1, g_0_yyyyyy_0_xxxyyy_1, g_0_yyyyyy_0_xxxyyyy_0, g_0_yyyyyy_0_xxxyyyy_1, g_0_yyyyyy_0_xxxyyyz_0, g_0_yyyyyy_0_xxxyyyz_1, g_0_yyyyyy_0_xxxyyz_1, g_0_yyyyyy_0_xxxyyzz_0, g_0_yyyyyy_0_xxxyyzz_1, g_0_yyyyyy_0_xxxyzz_1, g_0_yyyyyy_0_xxxyzzz_0, g_0_yyyyyy_0_xxxyzzz_1, g_0_yyyyyy_0_xxxzzz_1, g_0_yyyyyy_0_xxxzzzz_0, g_0_yyyyyy_0_xxxzzzz_1, g_0_yyyyyy_0_xxyyyy_1, g_0_yyyyyy_0_xxyyyyy_0, g_0_yyyyyy_0_xxyyyyy_1, g_0_yyyyyy_0_xxyyyyz_0, g_0_yyyyyy_0_xxyyyyz_1, g_0_yyyyyy_0_xxyyyz_1, g_0_yyyyyy_0_xxyyyzz_0, g_0_yyyyyy_0_xxyyyzz_1, g_0_yyyyyy_0_xxyyzz_1, g_0_yyyyyy_0_xxyyzzz_0, g_0_yyyyyy_0_xxyyzzz_1, g_0_yyyyyy_0_xxyzzz_1, g_0_yyyyyy_0_xxyzzzz_0, g_0_yyyyyy_0_xxyzzzz_1, g_0_yyyyyy_0_xxzzzz_1, g_0_yyyyyy_0_xxzzzzz_0, g_0_yyyyyy_0_xxzzzzz_1, g_0_yyyyyy_0_xyyyyy_1, g_0_yyyyyy_0_xyyyyyy_0, g_0_yyyyyy_0_xyyyyyy_1, g_0_yyyyyy_0_xyyyyyz_0, g_0_yyyyyy_0_xyyyyyz_1, g_0_yyyyyy_0_xyyyyz_1, g_0_yyyyyy_0_xyyyyzz_0, g_0_yyyyyy_0_xyyyyzz_1, g_0_yyyyyy_0_xyyyzz_1, g_0_yyyyyy_0_xyyyzzz_0, g_0_yyyyyy_0_xyyyzzz_1, g_0_yyyyyy_0_xyyzzz_1, g_0_yyyyyy_0_xyyzzzz_0, g_0_yyyyyy_0_xyyzzzz_1, g_0_yyyyyy_0_xyzzzz_1, g_0_yyyyyy_0_xyzzzzz_0, g_0_yyyyyy_0_xyzzzzz_1, g_0_yyyyyy_0_xzzzzz_1, g_0_yyyyyy_0_xzzzzzz_0, g_0_yyyyyy_0_xzzzzzz_1, g_0_yyyyyy_0_yyyyyy_1, g_0_yyyyyy_0_yyyyyyy_0, g_0_yyyyyy_0_yyyyyyy_1, g_0_yyyyyy_0_yyyyyyz_0, g_0_yyyyyy_0_yyyyyyz_1, g_0_yyyyyy_0_yyyyyz_1, g_0_yyyyyy_0_yyyyyzz_0, g_0_yyyyyy_0_yyyyyzz_1, g_0_yyyyyy_0_yyyyzz_1, g_0_yyyyyy_0_yyyyzzz_0, g_0_yyyyyy_0_yyyyzzz_1, g_0_yyyyyy_0_yyyzzz_1, g_0_yyyyyy_0_yyyzzzz_0, g_0_yyyyyy_0_yyyzzzz_1, g_0_yyyyyy_0_yyzzzz_1, g_0_yyyyyy_0_yyzzzzz_0, g_0_yyyyyy_0_yyzzzzz_1, g_0_yyyyyy_0_yzzzzz_1, g_0_yyyyyy_0_yzzzzzz_0, g_0_yyyyyy_0_yzzzzzz_1, g_0_yyyyyy_0_zzzzzz_1, g_0_yyyyyy_0_zzzzzzz_0, g_0_yyyyyy_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyy_0_xxxxxxx_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxxx_0[i] * pb_x + g_0_yyyyyy_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxxxy_0[i] = 6.0 * g_0_yyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxxy_0[i] * pb_x + g_0_yyyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxxz_0[i] * pb_x + g_0_yyyyyy_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxxyy_0[i] = 5.0 * g_0_yyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxyy_0[i] * pb_x + g_0_yyyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxyz_0[i] * pb_x + g_0_yyyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxzz_0[i] * pb_x + g_0_yyyyyy_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxyyy_0[i] = 4.0 * g_0_yyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyyy_0[i] * pb_x + g_0_yyyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyyz_0[i] * pb_x + g_0_yyyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyzz_0[i] * pb_x + g_0_yyyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxzzz_0[i] * pb_x + g_0_yyyyyy_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyyyy_0[i] = 3.0 * g_0_yyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyyy_0[i] * pb_x + g_0_yyyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyyz_0[i] * pb_x + g_0_yyyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyzz_0[i] * pb_x + g_0_yyyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyzzz_0[i] * pb_x + g_0_yyyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxzzzz_0[i] * pb_x + g_0_yyyyyy_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyyyy_0[i] = 2.0 * g_0_yyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyyy_0[i] * pb_x + g_0_yyyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyyz_0[i] * pb_x + g_0_yyyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyzz_0[i] * pb_x + g_0_yyyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyzzz_0[i] * pb_x + g_0_yyyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzzzz_0[i] * pb_x + g_0_yyyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxzzzzz_0[i] * pb_x + g_0_yyyyyy_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyyyy_0[i] = g_0_yyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyyy_0[i] * pb_x + g_0_yyyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyyyz_0[i] = g_0_yyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyyz_0[i] * pb_x + g_0_yyyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyyzz_0[i] = g_0_yyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyzz_0[i] * pb_x + g_0_yyyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyzzz_0[i] = g_0_yyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyzzz_0[i] * pb_x + g_0_yyyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyzzzz_0[i] = g_0_yyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzzzz_0[i] * pb_x + g_0_yyyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyzzzzz_0[i] = g_0_yyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzzzz_0[i] * pb_x + g_0_yyyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xzzzzzz_0[i] = g_0_yyyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzzzzzz_0[i] * pb_x + g_0_yyyyyy_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyyyy_0[i] = g_0_yyyyyy_0_yyyyyyy_0[i] * pb_x + g_0_yyyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyyyz_0[i] = g_0_yyyyyy_0_yyyyyyz_0[i] * pb_x + g_0_yyyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyyzz_0[i] = g_0_yyyyyy_0_yyyyyzz_0[i] * pb_x + g_0_yyyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyzzz_0[i] = g_0_yyyyyy_0_yyyyzzz_0[i] * pb_x + g_0_yyyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyzzzz_0[i] = g_0_yyyyyy_0_yyyzzzz_0[i] * pb_x + g_0_yyyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyzzzzz_0[i] = g_0_yyyyyy_0_yyzzzzz_0[i] * pb_x + g_0_yyyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yzzzzzz_0[i] = g_0_yyyyyy_0_yzzzzzz_0[i] * pb_x + g_0_yyyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_zzzzzzz_0[i] = g_0_yyyyyy_0_zzzzzzz_0[i] * pb_x + g_0_yyyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 792-828 components of targeted buffer : SKSK

    auto g_0_xyyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 792);

    auto g_0_xyyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 793);

    auto g_0_xyyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 794);

    auto g_0_xyyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 795);

    auto g_0_xyyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 796);

    auto g_0_xyyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 797);

    auto g_0_xyyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 798);

    auto g_0_xyyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 799);

    auto g_0_xyyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 800);

    auto g_0_xyyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 801);

    auto g_0_xyyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 802);

    auto g_0_xyyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 803);

    auto g_0_xyyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 804);

    auto g_0_xyyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 805);

    auto g_0_xyyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 806);

    auto g_0_xyyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 807);

    auto g_0_xyyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 808);

    auto g_0_xyyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 809);

    auto g_0_xyyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 810);

    auto g_0_xyyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 811);

    auto g_0_xyyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 812);

    auto g_0_xyyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 813);

    auto g_0_xyyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 814);

    auto g_0_xyyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 815);

    auto g_0_xyyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 816);

    auto g_0_xyyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 817);

    auto g_0_xyyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 818);

    auto g_0_xyyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 819);

    auto g_0_xyyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 820);

    auto g_0_xyyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 821);

    auto g_0_xyyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 822);

    auto g_0_xyyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 823);

    auto g_0_xyyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 824);

    auto g_0_xyyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 825);

    auto g_0_xyyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 826);

    auto g_0_xyyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 827);

    #pragma omp simd aligned(g_0_xyyyyy_0_xxxxxxx_0, g_0_xyyyyy_0_xxxxxxx_1, g_0_xyyyyy_0_xxxxxxy_0, g_0_xyyyyy_0_xxxxxxy_1, g_0_xyyyyy_0_xxxxxyy_0, g_0_xyyyyy_0_xxxxxyy_1, g_0_xyyyyy_0_xxxxyyy_0, g_0_xyyyyy_0_xxxxyyy_1, g_0_xyyyyy_0_xxxyyyy_0, g_0_xyyyyy_0_xxxyyyy_1, g_0_xyyyyy_0_xxyyyyy_0, g_0_xyyyyy_0_xxyyyyy_1, g_0_xyyyyy_0_xyyyyyy_0, g_0_xyyyyy_0_xyyyyyy_1, g_0_xyyyyyz_0_xxxxxxx_0, g_0_xyyyyyz_0_xxxxxxy_0, g_0_xyyyyyz_0_xxxxxxz_0, g_0_xyyyyyz_0_xxxxxyy_0, g_0_xyyyyyz_0_xxxxxyz_0, g_0_xyyyyyz_0_xxxxxzz_0, g_0_xyyyyyz_0_xxxxyyy_0, g_0_xyyyyyz_0_xxxxyyz_0, g_0_xyyyyyz_0_xxxxyzz_0, g_0_xyyyyyz_0_xxxxzzz_0, g_0_xyyyyyz_0_xxxyyyy_0, g_0_xyyyyyz_0_xxxyyyz_0, g_0_xyyyyyz_0_xxxyyzz_0, g_0_xyyyyyz_0_xxxyzzz_0, g_0_xyyyyyz_0_xxxzzzz_0, g_0_xyyyyyz_0_xxyyyyy_0, g_0_xyyyyyz_0_xxyyyyz_0, g_0_xyyyyyz_0_xxyyyzz_0, g_0_xyyyyyz_0_xxyyzzz_0, g_0_xyyyyyz_0_xxyzzzz_0, g_0_xyyyyyz_0_xxzzzzz_0, g_0_xyyyyyz_0_xyyyyyy_0, g_0_xyyyyyz_0_xyyyyyz_0, g_0_xyyyyyz_0_xyyyyzz_0, g_0_xyyyyyz_0_xyyyzzz_0, g_0_xyyyyyz_0_xyyzzzz_0, g_0_xyyyyyz_0_xyzzzzz_0, g_0_xyyyyyz_0_xzzzzzz_0, g_0_xyyyyyz_0_yyyyyyy_0, g_0_xyyyyyz_0_yyyyyyz_0, g_0_xyyyyyz_0_yyyyyzz_0, g_0_xyyyyyz_0_yyyyzzz_0, g_0_xyyyyyz_0_yyyzzzz_0, g_0_xyyyyyz_0_yyzzzzz_0, g_0_xyyyyyz_0_yzzzzzz_0, g_0_xyyyyyz_0_zzzzzzz_0, g_0_yyyyyz_0_xxxxxxz_0, g_0_yyyyyz_0_xxxxxxz_1, g_0_yyyyyz_0_xxxxxyz_0, g_0_yyyyyz_0_xxxxxyz_1, g_0_yyyyyz_0_xxxxxz_1, g_0_yyyyyz_0_xxxxxzz_0, g_0_yyyyyz_0_xxxxxzz_1, g_0_yyyyyz_0_xxxxyyz_0, g_0_yyyyyz_0_xxxxyyz_1, g_0_yyyyyz_0_xxxxyz_1, g_0_yyyyyz_0_xxxxyzz_0, g_0_yyyyyz_0_xxxxyzz_1, g_0_yyyyyz_0_xxxxzz_1, g_0_yyyyyz_0_xxxxzzz_0, g_0_yyyyyz_0_xxxxzzz_1, g_0_yyyyyz_0_xxxyyyz_0, g_0_yyyyyz_0_xxxyyyz_1, g_0_yyyyyz_0_xxxyyz_1, g_0_yyyyyz_0_xxxyyzz_0, g_0_yyyyyz_0_xxxyyzz_1, g_0_yyyyyz_0_xxxyzz_1, g_0_yyyyyz_0_xxxyzzz_0, g_0_yyyyyz_0_xxxyzzz_1, g_0_yyyyyz_0_xxxzzz_1, g_0_yyyyyz_0_xxxzzzz_0, g_0_yyyyyz_0_xxxzzzz_1, g_0_yyyyyz_0_xxyyyyz_0, g_0_yyyyyz_0_xxyyyyz_1, g_0_yyyyyz_0_xxyyyz_1, g_0_yyyyyz_0_xxyyyzz_0, g_0_yyyyyz_0_xxyyyzz_1, g_0_yyyyyz_0_xxyyzz_1, g_0_yyyyyz_0_xxyyzzz_0, g_0_yyyyyz_0_xxyyzzz_1, g_0_yyyyyz_0_xxyzzz_1, g_0_yyyyyz_0_xxyzzzz_0, g_0_yyyyyz_0_xxyzzzz_1, g_0_yyyyyz_0_xxzzzz_1, g_0_yyyyyz_0_xxzzzzz_0, g_0_yyyyyz_0_xxzzzzz_1, g_0_yyyyyz_0_xyyyyyz_0, g_0_yyyyyz_0_xyyyyyz_1, g_0_yyyyyz_0_xyyyyz_1, g_0_yyyyyz_0_xyyyyzz_0, g_0_yyyyyz_0_xyyyyzz_1, g_0_yyyyyz_0_xyyyzz_1, g_0_yyyyyz_0_xyyyzzz_0, g_0_yyyyyz_0_xyyyzzz_1, g_0_yyyyyz_0_xyyzzz_1, g_0_yyyyyz_0_xyyzzzz_0, g_0_yyyyyz_0_xyyzzzz_1, g_0_yyyyyz_0_xyzzzz_1, g_0_yyyyyz_0_xyzzzzz_0, g_0_yyyyyz_0_xyzzzzz_1, g_0_yyyyyz_0_xzzzzz_1, g_0_yyyyyz_0_xzzzzzz_0, g_0_yyyyyz_0_xzzzzzz_1, g_0_yyyyyz_0_yyyyyyy_0, g_0_yyyyyz_0_yyyyyyy_1, g_0_yyyyyz_0_yyyyyyz_0, g_0_yyyyyz_0_yyyyyyz_1, g_0_yyyyyz_0_yyyyyz_1, g_0_yyyyyz_0_yyyyyzz_0, g_0_yyyyyz_0_yyyyyzz_1, g_0_yyyyyz_0_yyyyzz_1, g_0_yyyyyz_0_yyyyzzz_0, g_0_yyyyyz_0_yyyyzzz_1, g_0_yyyyyz_0_yyyzzz_1, g_0_yyyyyz_0_yyyzzzz_0, g_0_yyyyyz_0_yyyzzzz_1, g_0_yyyyyz_0_yyzzzz_1, g_0_yyyyyz_0_yyzzzzz_0, g_0_yyyyyz_0_yyzzzzz_1, g_0_yyyyyz_0_yzzzzz_1, g_0_yyyyyz_0_yzzzzzz_0, g_0_yyyyyz_0_yzzzzzz_1, g_0_yyyyyz_0_zzzzzz_1, g_0_yyyyyz_0_zzzzzzz_0, g_0_yyyyyz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyz_0_xxxxxxx_0[i] = g_0_xyyyyy_0_xxxxxxx_0[i] * pb_z + g_0_xyyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxxxy_0[i] = g_0_xyyyyy_0_xxxxxxy_0[i] * pb_z + g_0_xyyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyyz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxxxz_0[i] * pb_x + g_0_yyyyyz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxxxyy_0[i] = g_0_xyyyyy_0_xxxxxyy_0[i] * pb_z + g_0_xyyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyyz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxxyz_0[i] * pb_x + g_0_yyyyyz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyyz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxxzz_0[i] * pb_x + g_0_yyyyyz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxxyyy_0[i] = g_0_xyyyyy_0_xxxxyyy_0[i] * pb_z + g_0_xyyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyyz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxyyz_0[i] * pb_x + g_0_yyyyyz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyyz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxyzz_0[i] * pb_x + g_0_yyyyyz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyyz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxzzz_0[i] * pb_x + g_0_yyyyyz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxyyyy_0[i] = g_0_xyyyyy_0_xxxyyyy_0[i] * pb_z + g_0_xyyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyyz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxyyyz_0[i] * pb_x + g_0_yyyyyz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyyz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxyyzz_0[i] * pb_x + g_0_yyyyyz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyyz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxyzzz_0[i] * pb_x + g_0_yyyyyz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyyz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxzzzz_0[i] * pb_x + g_0_yyyyyz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyyyyy_0[i] = g_0_xyyyyy_0_xxyyyyy_0[i] * pb_z + g_0_xyyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyyz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyyyyz_0[i] * pb_x + g_0_yyyyyz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyyz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyyyzz_0[i] * pb_x + g_0_yyyyyz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyyz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyyzzz_0[i] * pb_x + g_0_yyyyyz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyyz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyzzzz_0[i] * pb_x + g_0_yyyyyz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyyz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxzzzzz_0[i] * pb_x + g_0_yyyyyz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyyyyy_0[i] = g_0_xyyyyy_0_xyyyyyy_0[i] * pb_z + g_0_xyyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xyyyyyz_0[i] = g_0_yyyyyz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyyyyz_0[i] * pb_x + g_0_yyyyyz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyyyzz_0[i] = g_0_yyyyyz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyyyzz_0[i] * pb_x + g_0_yyyyyz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyyzzz_0[i] = g_0_yyyyyz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyyzzz_0[i] * pb_x + g_0_yyyyyz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyzzzz_0[i] = g_0_yyyyyz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyzzzz_0[i] * pb_x + g_0_yyyyyz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyzzzzz_0[i] = g_0_yyyyyz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyzzzzz_0[i] * pb_x + g_0_yyyyyz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xzzzzzz_0[i] = g_0_yyyyyz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xzzzzzz_0[i] * pb_x + g_0_yyyyyz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyyyy_0[i] = g_0_yyyyyz_0_yyyyyyy_0[i] * pb_x + g_0_yyyyyz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyyyz_0[i] = g_0_yyyyyz_0_yyyyyyz_0[i] * pb_x + g_0_yyyyyz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyyzz_0[i] = g_0_yyyyyz_0_yyyyyzz_0[i] * pb_x + g_0_yyyyyz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyzzz_0[i] = g_0_yyyyyz_0_yyyyzzz_0[i] * pb_x + g_0_yyyyyz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyzzzz_0[i] = g_0_yyyyyz_0_yyyzzzz_0[i] * pb_x + g_0_yyyyyz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyzzzzz_0[i] = g_0_yyyyyz_0_yyzzzzz_0[i] * pb_x + g_0_yyyyyz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yzzzzzz_0[i] = g_0_yyyyyz_0_yzzzzzz_0[i] * pb_x + g_0_yyyyyz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_zzzzzzz_0[i] = g_0_yyyyyz_0_zzzzzzz_0[i] * pb_x + g_0_yyyyyz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 828-864 components of targeted buffer : SKSK

    auto g_0_xyyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 828);

    auto g_0_xyyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 829);

    auto g_0_xyyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 830);

    auto g_0_xyyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 831);

    auto g_0_xyyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 832);

    auto g_0_xyyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 833);

    auto g_0_xyyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 834);

    auto g_0_xyyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 835);

    auto g_0_xyyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 836);

    auto g_0_xyyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 837);

    auto g_0_xyyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 838);

    auto g_0_xyyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 839);

    auto g_0_xyyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 840);

    auto g_0_xyyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 841);

    auto g_0_xyyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 842);

    auto g_0_xyyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 843);

    auto g_0_xyyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 844);

    auto g_0_xyyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 845);

    auto g_0_xyyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 846);

    auto g_0_xyyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 847);

    auto g_0_xyyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 848);

    auto g_0_xyyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 849);

    auto g_0_xyyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 850);

    auto g_0_xyyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 851);

    auto g_0_xyyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 852);

    auto g_0_xyyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 853);

    auto g_0_xyyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 854);

    auto g_0_xyyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 855);

    auto g_0_xyyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 856);

    auto g_0_xyyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 857);

    auto g_0_xyyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 858);

    auto g_0_xyyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 859);

    auto g_0_xyyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 860);

    auto g_0_xyyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 861);

    auto g_0_xyyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 862);

    auto g_0_xyyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 863);

    #pragma omp simd aligned(g_0_xyyyyzz_0_xxxxxxx_0, g_0_xyyyyzz_0_xxxxxxy_0, g_0_xyyyyzz_0_xxxxxxz_0, g_0_xyyyyzz_0_xxxxxyy_0, g_0_xyyyyzz_0_xxxxxyz_0, g_0_xyyyyzz_0_xxxxxzz_0, g_0_xyyyyzz_0_xxxxyyy_0, g_0_xyyyyzz_0_xxxxyyz_0, g_0_xyyyyzz_0_xxxxyzz_0, g_0_xyyyyzz_0_xxxxzzz_0, g_0_xyyyyzz_0_xxxyyyy_0, g_0_xyyyyzz_0_xxxyyyz_0, g_0_xyyyyzz_0_xxxyyzz_0, g_0_xyyyyzz_0_xxxyzzz_0, g_0_xyyyyzz_0_xxxzzzz_0, g_0_xyyyyzz_0_xxyyyyy_0, g_0_xyyyyzz_0_xxyyyyz_0, g_0_xyyyyzz_0_xxyyyzz_0, g_0_xyyyyzz_0_xxyyzzz_0, g_0_xyyyyzz_0_xxyzzzz_0, g_0_xyyyyzz_0_xxzzzzz_0, g_0_xyyyyzz_0_xyyyyyy_0, g_0_xyyyyzz_0_xyyyyyz_0, g_0_xyyyyzz_0_xyyyyzz_0, g_0_xyyyyzz_0_xyyyzzz_0, g_0_xyyyyzz_0_xyyzzzz_0, g_0_xyyyyzz_0_xyzzzzz_0, g_0_xyyyyzz_0_xzzzzzz_0, g_0_xyyyyzz_0_yyyyyyy_0, g_0_xyyyyzz_0_yyyyyyz_0, g_0_xyyyyzz_0_yyyyyzz_0, g_0_xyyyyzz_0_yyyyzzz_0, g_0_xyyyyzz_0_yyyzzzz_0, g_0_xyyyyzz_0_yyzzzzz_0, g_0_xyyyyzz_0_yzzzzzz_0, g_0_xyyyyzz_0_zzzzzzz_0, g_0_yyyyzz_0_xxxxxx_1, g_0_yyyyzz_0_xxxxxxx_0, g_0_yyyyzz_0_xxxxxxx_1, g_0_yyyyzz_0_xxxxxxy_0, g_0_yyyyzz_0_xxxxxxy_1, g_0_yyyyzz_0_xxxxxxz_0, g_0_yyyyzz_0_xxxxxxz_1, g_0_yyyyzz_0_xxxxxy_1, g_0_yyyyzz_0_xxxxxyy_0, g_0_yyyyzz_0_xxxxxyy_1, g_0_yyyyzz_0_xxxxxyz_0, g_0_yyyyzz_0_xxxxxyz_1, g_0_yyyyzz_0_xxxxxz_1, g_0_yyyyzz_0_xxxxxzz_0, g_0_yyyyzz_0_xxxxxzz_1, g_0_yyyyzz_0_xxxxyy_1, g_0_yyyyzz_0_xxxxyyy_0, g_0_yyyyzz_0_xxxxyyy_1, g_0_yyyyzz_0_xxxxyyz_0, g_0_yyyyzz_0_xxxxyyz_1, g_0_yyyyzz_0_xxxxyz_1, g_0_yyyyzz_0_xxxxyzz_0, g_0_yyyyzz_0_xxxxyzz_1, g_0_yyyyzz_0_xxxxzz_1, g_0_yyyyzz_0_xxxxzzz_0, g_0_yyyyzz_0_xxxxzzz_1, g_0_yyyyzz_0_xxxyyy_1, g_0_yyyyzz_0_xxxyyyy_0, g_0_yyyyzz_0_xxxyyyy_1, g_0_yyyyzz_0_xxxyyyz_0, g_0_yyyyzz_0_xxxyyyz_1, g_0_yyyyzz_0_xxxyyz_1, g_0_yyyyzz_0_xxxyyzz_0, g_0_yyyyzz_0_xxxyyzz_1, g_0_yyyyzz_0_xxxyzz_1, g_0_yyyyzz_0_xxxyzzz_0, g_0_yyyyzz_0_xxxyzzz_1, g_0_yyyyzz_0_xxxzzz_1, g_0_yyyyzz_0_xxxzzzz_0, g_0_yyyyzz_0_xxxzzzz_1, g_0_yyyyzz_0_xxyyyy_1, g_0_yyyyzz_0_xxyyyyy_0, g_0_yyyyzz_0_xxyyyyy_1, g_0_yyyyzz_0_xxyyyyz_0, g_0_yyyyzz_0_xxyyyyz_1, g_0_yyyyzz_0_xxyyyz_1, g_0_yyyyzz_0_xxyyyzz_0, g_0_yyyyzz_0_xxyyyzz_1, g_0_yyyyzz_0_xxyyzz_1, g_0_yyyyzz_0_xxyyzzz_0, g_0_yyyyzz_0_xxyyzzz_1, g_0_yyyyzz_0_xxyzzz_1, g_0_yyyyzz_0_xxyzzzz_0, g_0_yyyyzz_0_xxyzzzz_1, g_0_yyyyzz_0_xxzzzz_1, g_0_yyyyzz_0_xxzzzzz_0, g_0_yyyyzz_0_xxzzzzz_1, g_0_yyyyzz_0_xyyyyy_1, g_0_yyyyzz_0_xyyyyyy_0, g_0_yyyyzz_0_xyyyyyy_1, g_0_yyyyzz_0_xyyyyyz_0, g_0_yyyyzz_0_xyyyyyz_1, g_0_yyyyzz_0_xyyyyz_1, g_0_yyyyzz_0_xyyyyzz_0, g_0_yyyyzz_0_xyyyyzz_1, g_0_yyyyzz_0_xyyyzz_1, g_0_yyyyzz_0_xyyyzzz_0, g_0_yyyyzz_0_xyyyzzz_1, g_0_yyyyzz_0_xyyzzz_1, g_0_yyyyzz_0_xyyzzzz_0, g_0_yyyyzz_0_xyyzzzz_1, g_0_yyyyzz_0_xyzzzz_1, g_0_yyyyzz_0_xyzzzzz_0, g_0_yyyyzz_0_xyzzzzz_1, g_0_yyyyzz_0_xzzzzz_1, g_0_yyyyzz_0_xzzzzzz_0, g_0_yyyyzz_0_xzzzzzz_1, g_0_yyyyzz_0_yyyyyy_1, g_0_yyyyzz_0_yyyyyyy_0, g_0_yyyyzz_0_yyyyyyy_1, g_0_yyyyzz_0_yyyyyyz_0, g_0_yyyyzz_0_yyyyyyz_1, g_0_yyyyzz_0_yyyyyz_1, g_0_yyyyzz_0_yyyyyzz_0, g_0_yyyyzz_0_yyyyyzz_1, g_0_yyyyzz_0_yyyyzz_1, g_0_yyyyzz_0_yyyyzzz_0, g_0_yyyyzz_0_yyyyzzz_1, g_0_yyyyzz_0_yyyzzz_1, g_0_yyyyzz_0_yyyzzzz_0, g_0_yyyyzz_0_yyyzzzz_1, g_0_yyyyzz_0_yyzzzz_1, g_0_yyyyzz_0_yyzzzzz_0, g_0_yyyyzz_0_yyzzzzz_1, g_0_yyyyzz_0_yzzzzz_1, g_0_yyyyzz_0_yzzzzzz_0, g_0_yyyyzz_0_yzzzzzz_1, g_0_yyyyzz_0_zzzzzz_1, g_0_yyyyzz_0_zzzzzzz_0, g_0_yyyyzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyyyzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxxx_0[i] * pb_x + g_0_yyyyzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyyyzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxxy_0[i] * pb_x + g_0_yyyyzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxxz_0[i] * pb_x + g_0_yyyyzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyyyzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxyy_0[i] * pb_x + g_0_yyyyzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxyz_0[i] * pb_x + g_0_yyyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxzz_0[i] * pb_x + g_0_yyyyzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyyyzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxyyy_0[i] * pb_x + g_0_yyyyzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxyyz_0[i] * pb_x + g_0_yyyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxyzz_0[i] * pb_x + g_0_yyyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxzzz_0[i] * pb_x + g_0_yyyyzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyyyzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyyyy_0[i] * pb_x + g_0_yyyyzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyyyz_0[i] * pb_x + g_0_yyyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyyzz_0[i] * pb_x + g_0_yyyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyzzz_0[i] * pb_x + g_0_yyyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxzzzz_0[i] * pb_x + g_0_yyyyzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyyzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyyyy_0[i] * pb_x + g_0_yyyyzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyyyz_0[i] * pb_x + g_0_yyyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyyzz_0[i] * pb_x + g_0_yyyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyzzz_0[i] * pb_x + g_0_yyyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyzzzz_0[i] * pb_x + g_0_yyyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxzzzzz_0[i] * pb_x + g_0_yyyyzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyyyy_0[i] = g_0_yyyyzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyyyy_0[i] * pb_x + g_0_yyyyzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyyyz_0[i] = g_0_yyyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyyyz_0[i] * pb_x + g_0_yyyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyyzz_0[i] = g_0_yyyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyyzz_0[i] * pb_x + g_0_yyyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyzzz_0[i] = g_0_yyyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyzzz_0[i] * pb_x + g_0_yyyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyzzzz_0[i] = g_0_yyyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyzzzz_0[i] * pb_x + g_0_yyyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyzzzzz_0[i] = g_0_yyyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyzzzzz_0[i] * pb_x + g_0_yyyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xzzzzzz_0[i] = g_0_yyyyzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xzzzzzz_0[i] * pb_x + g_0_yyyyzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyyyy_0[i] = g_0_yyyyzz_0_yyyyyyy_0[i] * pb_x + g_0_yyyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyyyz_0[i] = g_0_yyyyzz_0_yyyyyyz_0[i] * pb_x + g_0_yyyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyyzz_0[i] = g_0_yyyyzz_0_yyyyyzz_0[i] * pb_x + g_0_yyyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyzzz_0[i] = g_0_yyyyzz_0_yyyyzzz_0[i] * pb_x + g_0_yyyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyzzzz_0[i] = g_0_yyyyzz_0_yyyzzzz_0[i] * pb_x + g_0_yyyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyzzzzz_0[i] = g_0_yyyyzz_0_yyzzzzz_0[i] * pb_x + g_0_yyyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yzzzzzz_0[i] = g_0_yyyyzz_0_yzzzzzz_0[i] * pb_x + g_0_yyyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_zzzzzzz_0[i] = g_0_yyyyzz_0_zzzzzzz_0[i] * pb_x + g_0_yyyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 864-900 components of targeted buffer : SKSK

    auto g_0_xyyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 864);

    auto g_0_xyyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 865);

    auto g_0_xyyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 866);

    auto g_0_xyyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 867);

    auto g_0_xyyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 868);

    auto g_0_xyyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 869);

    auto g_0_xyyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 870);

    auto g_0_xyyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 871);

    auto g_0_xyyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 872);

    auto g_0_xyyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 873);

    auto g_0_xyyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 874);

    auto g_0_xyyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 875);

    auto g_0_xyyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 876);

    auto g_0_xyyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 877);

    auto g_0_xyyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 878);

    auto g_0_xyyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 879);

    auto g_0_xyyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 880);

    auto g_0_xyyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 881);

    auto g_0_xyyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 882);

    auto g_0_xyyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 883);

    auto g_0_xyyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 884);

    auto g_0_xyyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 885);

    auto g_0_xyyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 886);

    auto g_0_xyyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 887);

    auto g_0_xyyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 888);

    auto g_0_xyyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 889);

    auto g_0_xyyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 890);

    auto g_0_xyyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 891);

    auto g_0_xyyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 892);

    auto g_0_xyyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 893);

    auto g_0_xyyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 894);

    auto g_0_xyyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 895);

    auto g_0_xyyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 896);

    auto g_0_xyyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 897);

    auto g_0_xyyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 898);

    auto g_0_xyyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 899);

    #pragma omp simd aligned(g_0_xyyyzzz_0_xxxxxxx_0, g_0_xyyyzzz_0_xxxxxxy_0, g_0_xyyyzzz_0_xxxxxxz_0, g_0_xyyyzzz_0_xxxxxyy_0, g_0_xyyyzzz_0_xxxxxyz_0, g_0_xyyyzzz_0_xxxxxzz_0, g_0_xyyyzzz_0_xxxxyyy_0, g_0_xyyyzzz_0_xxxxyyz_0, g_0_xyyyzzz_0_xxxxyzz_0, g_0_xyyyzzz_0_xxxxzzz_0, g_0_xyyyzzz_0_xxxyyyy_0, g_0_xyyyzzz_0_xxxyyyz_0, g_0_xyyyzzz_0_xxxyyzz_0, g_0_xyyyzzz_0_xxxyzzz_0, g_0_xyyyzzz_0_xxxzzzz_0, g_0_xyyyzzz_0_xxyyyyy_0, g_0_xyyyzzz_0_xxyyyyz_0, g_0_xyyyzzz_0_xxyyyzz_0, g_0_xyyyzzz_0_xxyyzzz_0, g_0_xyyyzzz_0_xxyzzzz_0, g_0_xyyyzzz_0_xxzzzzz_0, g_0_xyyyzzz_0_xyyyyyy_0, g_0_xyyyzzz_0_xyyyyyz_0, g_0_xyyyzzz_0_xyyyyzz_0, g_0_xyyyzzz_0_xyyyzzz_0, g_0_xyyyzzz_0_xyyzzzz_0, g_0_xyyyzzz_0_xyzzzzz_0, g_0_xyyyzzz_0_xzzzzzz_0, g_0_xyyyzzz_0_yyyyyyy_0, g_0_xyyyzzz_0_yyyyyyz_0, g_0_xyyyzzz_0_yyyyyzz_0, g_0_xyyyzzz_0_yyyyzzz_0, g_0_xyyyzzz_0_yyyzzzz_0, g_0_xyyyzzz_0_yyzzzzz_0, g_0_xyyyzzz_0_yzzzzzz_0, g_0_xyyyzzz_0_zzzzzzz_0, g_0_yyyzzz_0_xxxxxx_1, g_0_yyyzzz_0_xxxxxxx_0, g_0_yyyzzz_0_xxxxxxx_1, g_0_yyyzzz_0_xxxxxxy_0, g_0_yyyzzz_0_xxxxxxy_1, g_0_yyyzzz_0_xxxxxxz_0, g_0_yyyzzz_0_xxxxxxz_1, g_0_yyyzzz_0_xxxxxy_1, g_0_yyyzzz_0_xxxxxyy_0, g_0_yyyzzz_0_xxxxxyy_1, g_0_yyyzzz_0_xxxxxyz_0, g_0_yyyzzz_0_xxxxxyz_1, g_0_yyyzzz_0_xxxxxz_1, g_0_yyyzzz_0_xxxxxzz_0, g_0_yyyzzz_0_xxxxxzz_1, g_0_yyyzzz_0_xxxxyy_1, g_0_yyyzzz_0_xxxxyyy_0, g_0_yyyzzz_0_xxxxyyy_1, g_0_yyyzzz_0_xxxxyyz_0, g_0_yyyzzz_0_xxxxyyz_1, g_0_yyyzzz_0_xxxxyz_1, g_0_yyyzzz_0_xxxxyzz_0, g_0_yyyzzz_0_xxxxyzz_1, g_0_yyyzzz_0_xxxxzz_1, g_0_yyyzzz_0_xxxxzzz_0, g_0_yyyzzz_0_xxxxzzz_1, g_0_yyyzzz_0_xxxyyy_1, g_0_yyyzzz_0_xxxyyyy_0, g_0_yyyzzz_0_xxxyyyy_1, g_0_yyyzzz_0_xxxyyyz_0, g_0_yyyzzz_0_xxxyyyz_1, g_0_yyyzzz_0_xxxyyz_1, g_0_yyyzzz_0_xxxyyzz_0, g_0_yyyzzz_0_xxxyyzz_1, g_0_yyyzzz_0_xxxyzz_1, g_0_yyyzzz_0_xxxyzzz_0, g_0_yyyzzz_0_xxxyzzz_1, g_0_yyyzzz_0_xxxzzz_1, g_0_yyyzzz_0_xxxzzzz_0, g_0_yyyzzz_0_xxxzzzz_1, g_0_yyyzzz_0_xxyyyy_1, g_0_yyyzzz_0_xxyyyyy_0, g_0_yyyzzz_0_xxyyyyy_1, g_0_yyyzzz_0_xxyyyyz_0, g_0_yyyzzz_0_xxyyyyz_1, g_0_yyyzzz_0_xxyyyz_1, g_0_yyyzzz_0_xxyyyzz_0, g_0_yyyzzz_0_xxyyyzz_1, g_0_yyyzzz_0_xxyyzz_1, g_0_yyyzzz_0_xxyyzzz_0, g_0_yyyzzz_0_xxyyzzz_1, g_0_yyyzzz_0_xxyzzz_1, g_0_yyyzzz_0_xxyzzzz_0, g_0_yyyzzz_0_xxyzzzz_1, g_0_yyyzzz_0_xxzzzz_1, g_0_yyyzzz_0_xxzzzzz_0, g_0_yyyzzz_0_xxzzzzz_1, g_0_yyyzzz_0_xyyyyy_1, g_0_yyyzzz_0_xyyyyyy_0, g_0_yyyzzz_0_xyyyyyy_1, g_0_yyyzzz_0_xyyyyyz_0, g_0_yyyzzz_0_xyyyyyz_1, g_0_yyyzzz_0_xyyyyz_1, g_0_yyyzzz_0_xyyyyzz_0, g_0_yyyzzz_0_xyyyyzz_1, g_0_yyyzzz_0_xyyyzz_1, g_0_yyyzzz_0_xyyyzzz_0, g_0_yyyzzz_0_xyyyzzz_1, g_0_yyyzzz_0_xyyzzz_1, g_0_yyyzzz_0_xyyzzzz_0, g_0_yyyzzz_0_xyyzzzz_1, g_0_yyyzzz_0_xyzzzz_1, g_0_yyyzzz_0_xyzzzzz_0, g_0_yyyzzz_0_xyzzzzz_1, g_0_yyyzzz_0_xzzzzz_1, g_0_yyyzzz_0_xzzzzzz_0, g_0_yyyzzz_0_xzzzzzz_1, g_0_yyyzzz_0_yyyyyy_1, g_0_yyyzzz_0_yyyyyyy_0, g_0_yyyzzz_0_yyyyyyy_1, g_0_yyyzzz_0_yyyyyyz_0, g_0_yyyzzz_0_yyyyyyz_1, g_0_yyyzzz_0_yyyyyz_1, g_0_yyyzzz_0_yyyyyzz_0, g_0_yyyzzz_0_yyyyyzz_1, g_0_yyyzzz_0_yyyyzz_1, g_0_yyyzzz_0_yyyyzzz_0, g_0_yyyzzz_0_yyyyzzz_1, g_0_yyyzzz_0_yyyzzz_1, g_0_yyyzzz_0_yyyzzzz_0, g_0_yyyzzz_0_yyyzzzz_1, g_0_yyyzzz_0_yyzzzz_1, g_0_yyyzzz_0_yyzzzzz_0, g_0_yyyzzz_0_yyzzzzz_1, g_0_yyyzzz_0_yzzzzz_1, g_0_yyyzzz_0_yzzzzzz_0, g_0_yyyzzz_0_yzzzzzz_1, g_0_yyyzzz_0_zzzzzz_1, g_0_yyyzzz_0_zzzzzzz_0, g_0_yyyzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyyzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxxx_0[i] * pb_x + g_0_yyyzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyyzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxxy_0[i] * pb_x + g_0_yyyzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxxz_0[i] * pb_x + g_0_yyyzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyyzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxyy_0[i] * pb_x + g_0_yyyzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxyz_0[i] * pb_x + g_0_yyyzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxzz_0[i] * pb_x + g_0_yyyzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyyzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxyyy_0[i] * pb_x + g_0_yyyzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxyyz_0[i] * pb_x + g_0_yyyzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxyzz_0[i] * pb_x + g_0_yyyzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxzzz_0[i] * pb_x + g_0_yyyzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyyzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyyyy_0[i] * pb_x + g_0_yyyzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyyyz_0[i] * pb_x + g_0_yyyzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyyzz_0[i] * pb_x + g_0_yyyzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyzzz_0[i] * pb_x + g_0_yyyzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxzzzz_0[i] * pb_x + g_0_yyyzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyyyy_0[i] * pb_x + g_0_yyyzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyyyz_0[i] * pb_x + g_0_yyyzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyyzz_0[i] * pb_x + g_0_yyyzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyzzz_0[i] * pb_x + g_0_yyyzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyzzzz_0[i] * pb_x + g_0_yyyzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxzzzzz_0[i] * pb_x + g_0_yyyzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyyyy_0[i] = g_0_yyyzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyyyy_0[i] * pb_x + g_0_yyyzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyyyz_0[i] = g_0_yyyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyyyz_0[i] * pb_x + g_0_yyyzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyyzz_0[i] = g_0_yyyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyyzz_0[i] * pb_x + g_0_yyyzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyzzz_0[i] = g_0_yyyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyzzz_0[i] * pb_x + g_0_yyyzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyzzzz_0[i] = g_0_yyyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyzzzz_0[i] * pb_x + g_0_yyyzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyzzzzz_0[i] = g_0_yyyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyzzzzz_0[i] * pb_x + g_0_yyyzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xzzzzzz_0[i] = g_0_yyyzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xzzzzzz_0[i] * pb_x + g_0_yyyzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyyyy_0[i] = g_0_yyyzzz_0_yyyyyyy_0[i] * pb_x + g_0_yyyzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyyyz_0[i] = g_0_yyyzzz_0_yyyyyyz_0[i] * pb_x + g_0_yyyzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyyzz_0[i] = g_0_yyyzzz_0_yyyyyzz_0[i] * pb_x + g_0_yyyzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyzzz_0[i] = g_0_yyyzzz_0_yyyyzzz_0[i] * pb_x + g_0_yyyzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyzzzz_0[i] = g_0_yyyzzz_0_yyyzzzz_0[i] * pb_x + g_0_yyyzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyzzzzz_0[i] = g_0_yyyzzz_0_yyzzzzz_0[i] * pb_x + g_0_yyyzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yzzzzzz_0[i] = g_0_yyyzzz_0_yzzzzzz_0[i] * pb_x + g_0_yyyzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_zzzzzzz_0[i] = g_0_yyyzzz_0_zzzzzzz_0[i] * pb_x + g_0_yyyzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 900-936 components of targeted buffer : SKSK

    auto g_0_xyyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 900);

    auto g_0_xyyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 901);

    auto g_0_xyyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 902);

    auto g_0_xyyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 903);

    auto g_0_xyyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 904);

    auto g_0_xyyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 905);

    auto g_0_xyyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 906);

    auto g_0_xyyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 907);

    auto g_0_xyyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 908);

    auto g_0_xyyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 909);

    auto g_0_xyyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 910);

    auto g_0_xyyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 911);

    auto g_0_xyyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 912);

    auto g_0_xyyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 913);

    auto g_0_xyyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 914);

    auto g_0_xyyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 915);

    auto g_0_xyyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 916);

    auto g_0_xyyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 917);

    auto g_0_xyyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 918);

    auto g_0_xyyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 919);

    auto g_0_xyyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 920);

    auto g_0_xyyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 921);

    auto g_0_xyyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 922);

    auto g_0_xyyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 923);

    auto g_0_xyyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 924);

    auto g_0_xyyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 925);

    auto g_0_xyyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 926);

    auto g_0_xyyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 927);

    auto g_0_xyyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 928);

    auto g_0_xyyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 929);

    auto g_0_xyyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 930);

    auto g_0_xyyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 931);

    auto g_0_xyyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 932);

    auto g_0_xyyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 933);

    auto g_0_xyyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 934);

    auto g_0_xyyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 935);

    #pragma omp simd aligned(g_0_xyyzzzz_0_xxxxxxx_0, g_0_xyyzzzz_0_xxxxxxy_0, g_0_xyyzzzz_0_xxxxxxz_0, g_0_xyyzzzz_0_xxxxxyy_0, g_0_xyyzzzz_0_xxxxxyz_0, g_0_xyyzzzz_0_xxxxxzz_0, g_0_xyyzzzz_0_xxxxyyy_0, g_0_xyyzzzz_0_xxxxyyz_0, g_0_xyyzzzz_0_xxxxyzz_0, g_0_xyyzzzz_0_xxxxzzz_0, g_0_xyyzzzz_0_xxxyyyy_0, g_0_xyyzzzz_0_xxxyyyz_0, g_0_xyyzzzz_0_xxxyyzz_0, g_0_xyyzzzz_0_xxxyzzz_0, g_0_xyyzzzz_0_xxxzzzz_0, g_0_xyyzzzz_0_xxyyyyy_0, g_0_xyyzzzz_0_xxyyyyz_0, g_0_xyyzzzz_0_xxyyyzz_0, g_0_xyyzzzz_0_xxyyzzz_0, g_0_xyyzzzz_0_xxyzzzz_0, g_0_xyyzzzz_0_xxzzzzz_0, g_0_xyyzzzz_0_xyyyyyy_0, g_0_xyyzzzz_0_xyyyyyz_0, g_0_xyyzzzz_0_xyyyyzz_0, g_0_xyyzzzz_0_xyyyzzz_0, g_0_xyyzzzz_0_xyyzzzz_0, g_0_xyyzzzz_0_xyzzzzz_0, g_0_xyyzzzz_0_xzzzzzz_0, g_0_xyyzzzz_0_yyyyyyy_0, g_0_xyyzzzz_0_yyyyyyz_0, g_0_xyyzzzz_0_yyyyyzz_0, g_0_xyyzzzz_0_yyyyzzz_0, g_0_xyyzzzz_0_yyyzzzz_0, g_0_xyyzzzz_0_yyzzzzz_0, g_0_xyyzzzz_0_yzzzzzz_0, g_0_xyyzzzz_0_zzzzzzz_0, g_0_yyzzzz_0_xxxxxx_1, g_0_yyzzzz_0_xxxxxxx_0, g_0_yyzzzz_0_xxxxxxx_1, g_0_yyzzzz_0_xxxxxxy_0, g_0_yyzzzz_0_xxxxxxy_1, g_0_yyzzzz_0_xxxxxxz_0, g_0_yyzzzz_0_xxxxxxz_1, g_0_yyzzzz_0_xxxxxy_1, g_0_yyzzzz_0_xxxxxyy_0, g_0_yyzzzz_0_xxxxxyy_1, g_0_yyzzzz_0_xxxxxyz_0, g_0_yyzzzz_0_xxxxxyz_1, g_0_yyzzzz_0_xxxxxz_1, g_0_yyzzzz_0_xxxxxzz_0, g_0_yyzzzz_0_xxxxxzz_1, g_0_yyzzzz_0_xxxxyy_1, g_0_yyzzzz_0_xxxxyyy_0, g_0_yyzzzz_0_xxxxyyy_1, g_0_yyzzzz_0_xxxxyyz_0, g_0_yyzzzz_0_xxxxyyz_1, g_0_yyzzzz_0_xxxxyz_1, g_0_yyzzzz_0_xxxxyzz_0, g_0_yyzzzz_0_xxxxyzz_1, g_0_yyzzzz_0_xxxxzz_1, g_0_yyzzzz_0_xxxxzzz_0, g_0_yyzzzz_0_xxxxzzz_1, g_0_yyzzzz_0_xxxyyy_1, g_0_yyzzzz_0_xxxyyyy_0, g_0_yyzzzz_0_xxxyyyy_1, g_0_yyzzzz_0_xxxyyyz_0, g_0_yyzzzz_0_xxxyyyz_1, g_0_yyzzzz_0_xxxyyz_1, g_0_yyzzzz_0_xxxyyzz_0, g_0_yyzzzz_0_xxxyyzz_1, g_0_yyzzzz_0_xxxyzz_1, g_0_yyzzzz_0_xxxyzzz_0, g_0_yyzzzz_0_xxxyzzz_1, g_0_yyzzzz_0_xxxzzz_1, g_0_yyzzzz_0_xxxzzzz_0, g_0_yyzzzz_0_xxxzzzz_1, g_0_yyzzzz_0_xxyyyy_1, g_0_yyzzzz_0_xxyyyyy_0, g_0_yyzzzz_0_xxyyyyy_1, g_0_yyzzzz_0_xxyyyyz_0, g_0_yyzzzz_0_xxyyyyz_1, g_0_yyzzzz_0_xxyyyz_1, g_0_yyzzzz_0_xxyyyzz_0, g_0_yyzzzz_0_xxyyyzz_1, g_0_yyzzzz_0_xxyyzz_1, g_0_yyzzzz_0_xxyyzzz_0, g_0_yyzzzz_0_xxyyzzz_1, g_0_yyzzzz_0_xxyzzz_1, g_0_yyzzzz_0_xxyzzzz_0, g_0_yyzzzz_0_xxyzzzz_1, g_0_yyzzzz_0_xxzzzz_1, g_0_yyzzzz_0_xxzzzzz_0, g_0_yyzzzz_0_xxzzzzz_1, g_0_yyzzzz_0_xyyyyy_1, g_0_yyzzzz_0_xyyyyyy_0, g_0_yyzzzz_0_xyyyyyy_1, g_0_yyzzzz_0_xyyyyyz_0, g_0_yyzzzz_0_xyyyyyz_1, g_0_yyzzzz_0_xyyyyz_1, g_0_yyzzzz_0_xyyyyzz_0, g_0_yyzzzz_0_xyyyyzz_1, g_0_yyzzzz_0_xyyyzz_1, g_0_yyzzzz_0_xyyyzzz_0, g_0_yyzzzz_0_xyyyzzz_1, g_0_yyzzzz_0_xyyzzz_1, g_0_yyzzzz_0_xyyzzzz_0, g_0_yyzzzz_0_xyyzzzz_1, g_0_yyzzzz_0_xyzzzz_1, g_0_yyzzzz_0_xyzzzzz_0, g_0_yyzzzz_0_xyzzzzz_1, g_0_yyzzzz_0_xzzzzz_1, g_0_yyzzzz_0_xzzzzzz_0, g_0_yyzzzz_0_xzzzzzz_1, g_0_yyzzzz_0_yyyyyy_1, g_0_yyzzzz_0_yyyyyyy_0, g_0_yyzzzz_0_yyyyyyy_1, g_0_yyzzzz_0_yyyyyyz_0, g_0_yyzzzz_0_yyyyyyz_1, g_0_yyzzzz_0_yyyyyz_1, g_0_yyzzzz_0_yyyyyzz_0, g_0_yyzzzz_0_yyyyyzz_1, g_0_yyzzzz_0_yyyyzz_1, g_0_yyzzzz_0_yyyyzzz_0, g_0_yyzzzz_0_yyyyzzz_1, g_0_yyzzzz_0_yyyzzz_1, g_0_yyzzzz_0_yyyzzzz_0, g_0_yyzzzz_0_yyyzzzz_1, g_0_yyzzzz_0_yyzzzz_1, g_0_yyzzzz_0_yyzzzzz_0, g_0_yyzzzz_0_yyzzzzz_1, g_0_yyzzzz_0_yzzzzz_1, g_0_yyzzzz_0_yzzzzzz_0, g_0_yyzzzz_0_yzzzzzz_1, g_0_yyzzzz_0_zzzzzz_1, g_0_yyzzzz_0_zzzzzzz_0, g_0_yyzzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxxx_0[i] * pb_x + g_0_yyzzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxxy_0[i] * pb_x + g_0_yyzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxxz_0[i] * pb_x + g_0_yyzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxyy_0[i] * pb_x + g_0_yyzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxyz_0[i] * pb_x + g_0_yyzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxzz_0[i] * pb_x + g_0_yyzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxyyy_0[i] * pb_x + g_0_yyzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxyyz_0[i] * pb_x + g_0_yyzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxyzz_0[i] * pb_x + g_0_yyzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxzzz_0[i] * pb_x + g_0_yyzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyyyy_0[i] * pb_x + g_0_yyzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyyyz_0[i] * pb_x + g_0_yyzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyyzz_0[i] * pb_x + g_0_yyzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyzzz_0[i] * pb_x + g_0_yyzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxzzzz_0[i] * pb_x + g_0_yyzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyyyy_0[i] * pb_x + g_0_yyzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyyyz_0[i] * pb_x + g_0_yyzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyyzz_0[i] * pb_x + g_0_yyzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyzzz_0[i] * pb_x + g_0_yyzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyzzzz_0[i] * pb_x + g_0_yyzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxzzzzz_0[i] * pb_x + g_0_yyzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyyyy_0[i] = g_0_yyzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyyyy_0[i] * pb_x + g_0_yyzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyyyz_0[i] = g_0_yyzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyyyz_0[i] * pb_x + g_0_yyzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyyzz_0[i] = g_0_yyzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyyzz_0[i] * pb_x + g_0_yyzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyzzz_0[i] = g_0_yyzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyzzz_0[i] * pb_x + g_0_yyzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyzzzz_0[i] = g_0_yyzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyzzzz_0[i] * pb_x + g_0_yyzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyzzzzz_0[i] = g_0_yyzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyzzzzz_0[i] * pb_x + g_0_yyzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xzzzzzz_0[i] = g_0_yyzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xzzzzzz_0[i] * pb_x + g_0_yyzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyyyy_0[i] = g_0_yyzzzz_0_yyyyyyy_0[i] * pb_x + g_0_yyzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyyyz_0[i] = g_0_yyzzzz_0_yyyyyyz_0[i] * pb_x + g_0_yyzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyyzz_0[i] = g_0_yyzzzz_0_yyyyyzz_0[i] * pb_x + g_0_yyzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyzzz_0[i] = g_0_yyzzzz_0_yyyyzzz_0[i] * pb_x + g_0_yyzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyzzzz_0[i] = g_0_yyzzzz_0_yyyzzzz_0[i] * pb_x + g_0_yyzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyzzzzz_0[i] = g_0_yyzzzz_0_yyzzzzz_0[i] * pb_x + g_0_yyzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yzzzzzz_0[i] = g_0_yyzzzz_0_yzzzzzz_0[i] * pb_x + g_0_yyzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_zzzzzzz_0[i] = g_0_yyzzzz_0_zzzzzzz_0[i] * pb_x + g_0_yyzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 936-972 components of targeted buffer : SKSK

    auto g_0_xyzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 936);

    auto g_0_xyzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 937);

    auto g_0_xyzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 938);

    auto g_0_xyzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 939);

    auto g_0_xyzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 940);

    auto g_0_xyzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 941);

    auto g_0_xyzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 942);

    auto g_0_xyzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 943);

    auto g_0_xyzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 944);

    auto g_0_xyzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 945);

    auto g_0_xyzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 946);

    auto g_0_xyzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 947);

    auto g_0_xyzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 948);

    auto g_0_xyzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 949);

    auto g_0_xyzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 950);

    auto g_0_xyzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 951);

    auto g_0_xyzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 952);

    auto g_0_xyzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 953);

    auto g_0_xyzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 954);

    auto g_0_xyzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 955);

    auto g_0_xyzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 956);

    auto g_0_xyzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 957);

    auto g_0_xyzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 958);

    auto g_0_xyzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 959);

    auto g_0_xyzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 960);

    auto g_0_xyzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 961);

    auto g_0_xyzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 962);

    auto g_0_xyzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 963);

    auto g_0_xyzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 964);

    auto g_0_xyzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 965);

    auto g_0_xyzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 966);

    auto g_0_xyzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 967);

    auto g_0_xyzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 968);

    auto g_0_xyzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 969);

    auto g_0_xyzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 970);

    auto g_0_xyzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 971);

    #pragma omp simd aligned(g_0_xyzzzzz_0_xxxxxxx_0, g_0_xyzzzzz_0_xxxxxxy_0, g_0_xyzzzzz_0_xxxxxxz_0, g_0_xyzzzzz_0_xxxxxyy_0, g_0_xyzzzzz_0_xxxxxyz_0, g_0_xyzzzzz_0_xxxxxzz_0, g_0_xyzzzzz_0_xxxxyyy_0, g_0_xyzzzzz_0_xxxxyyz_0, g_0_xyzzzzz_0_xxxxyzz_0, g_0_xyzzzzz_0_xxxxzzz_0, g_0_xyzzzzz_0_xxxyyyy_0, g_0_xyzzzzz_0_xxxyyyz_0, g_0_xyzzzzz_0_xxxyyzz_0, g_0_xyzzzzz_0_xxxyzzz_0, g_0_xyzzzzz_0_xxxzzzz_0, g_0_xyzzzzz_0_xxyyyyy_0, g_0_xyzzzzz_0_xxyyyyz_0, g_0_xyzzzzz_0_xxyyyzz_0, g_0_xyzzzzz_0_xxyyzzz_0, g_0_xyzzzzz_0_xxyzzzz_0, g_0_xyzzzzz_0_xxzzzzz_0, g_0_xyzzzzz_0_xyyyyyy_0, g_0_xyzzzzz_0_xyyyyyz_0, g_0_xyzzzzz_0_xyyyyzz_0, g_0_xyzzzzz_0_xyyyzzz_0, g_0_xyzzzzz_0_xyyzzzz_0, g_0_xyzzzzz_0_xyzzzzz_0, g_0_xyzzzzz_0_xzzzzzz_0, g_0_xyzzzzz_0_yyyyyyy_0, g_0_xyzzzzz_0_yyyyyyz_0, g_0_xyzzzzz_0_yyyyyzz_0, g_0_xyzzzzz_0_yyyyzzz_0, g_0_xyzzzzz_0_yyyzzzz_0, g_0_xyzzzzz_0_yyzzzzz_0, g_0_xyzzzzz_0_yzzzzzz_0, g_0_xyzzzzz_0_zzzzzzz_0, g_0_xzzzzz_0_xxxxxxx_0, g_0_xzzzzz_0_xxxxxxx_1, g_0_xzzzzz_0_xxxxxxz_0, g_0_xzzzzz_0_xxxxxxz_1, g_0_xzzzzz_0_xxxxxzz_0, g_0_xzzzzz_0_xxxxxzz_1, g_0_xzzzzz_0_xxxxzzz_0, g_0_xzzzzz_0_xxxxzzz_1, g_0_xzzzzz_0_xxxzzzz_0, g_0_xzzzzz_0_xxxzzzz_1, g_0_xzzzzz_0_xxzzzzz_0, g_0_xzzzzz_0_xxzzzzz_1, g_0_xzzzzz_0_xzzzzzz_0, g_0_xzzzzz_0_xzzzzzz_1, g_0_yzzzzz_0_xxxxxxy_0, g_0_yzzzzz_0_xxxxxxy_1, g_0_yzzzzz_0_xxxxxy_1, g_0_yzzzzz_0_xxxxxyy_0, g_0_yzzzzz_0_xxxxxyy_1, g_0_yzzzzz_0_xxxxxyz_0, g_0_yzzzzz_0_xxxxxyz_1, g_0_yzzzzz_0_xxxxyy_1, g_0_yzzzzz_0_xxxxyyy_0, g_0_yzzzzz_0_xxxxyyy_1, g_0_yzzzzz_0_xxxxyyz_0, g_0_yzzzzz_0_xxxxyyz_1, g_0_yzzzzz_0_xxxxyz_1, g_0_yzzzzz_0_xxxxyzz_0, g_0_yzzzzz_0_xxxxyzz_1, g_0_yzzzzz_0_xxxyyy_1, g_0_yzzzzz_0_xxxyyyy_0, g_0_yzzzzz_0_xxxyyyy_1, g_0_yzzzzz_0_xxxyyyz_0, g_0_yzzzzz_0_xxxyyyz_1, g_0_yzzzzz_0_xxxyyz_1, g_0_yzzzzz_0_xxxyyzz_0, g_0_yzzzzz_0_xxxyyzz_1, g_0_yzzzzz_0_xxxyzz_1, g_0_yzzzzz_0_xxxyzzz_0, g_0_yzzzzz_0_xxxyzzz_1, g_0_yzzzzz_0_xxyyyy_1, g_0_yzzzzz_0_xxyyyyy_0, g_0_yzzzzz_0_xxyyyyy_1, g_0_yzzzzz_0_xxyyyyz_0, g_0_yzzzzz_0_xxyyyyz_1, g_0_yzzzzz_0_xxyyyz_1, g_0_yzzzzz_0_xxyyyzz_0, g_0_yzzzzz_0_xxyyyzz_1, g_0_yzzzzz_0_xxyyzz_1, g_0_yzzzzz_0_xxyyzzz_0, g_0_yzzzzz_0_xxyyzzz_1, g_0_yzzzzz_0_xxyzzz_1, g_0_yzzzzz_0_xxyzzzz_0, g_0_yzzzzz_0_xxyzzzz_1, g_0_yzzzzz_0_xyyyyy_1, g_0_yzzzzz_0_xyyyyyy_0, g_0_yzzzzz_0_xyyyyyy_1, g_0_yzzzzz_0_xyyyyyz_0, g_0_yzzzzz_0_xyyyyyz_1, g_0_yzzzzz_0_xyyyyz_1, g_0_yzzzzz_0_xyyyyzz_0, g_0_yzzzzz_0_xyyyyzz_1, g_0_yzzzzz_0_xyyyzz_1, g_0_yzzzzz_0_xyyyzzz_0, g_0_yzzzzz_0_xyyyzzz_1, g_0_yzzzzz_0_xyyzzz_1, g_0_yzzzzz_0_xyyzzzz_0, g_0_yzzzzz_0_xyyzzzz_1, g_0_yzzzzz_0_xyzzzz_1, g_0_yzzzzz_0_xyzzzzz_0, g_0_yzzzzz_0_xyzzzzz_1, g_0_yzzzzz_0_yyyyyy_1, g_0_yzzzzz_0_yyyyyyy_0, g_0_yzzzzz_0_yyyyyyy_1, g_0_yzzzzz_0_yyyyyyz_0, g_0_yzzzzz_0_yyyyyyz_1, g_0_yzzzzz_0_yyyyyz_1, g_0_yzzzzz_0_yyyyyzz_0, g_0_yzzzzz_0_yyyyyzz_1, g_0_yzzzzz_0_yyyyzz_1, g_0_yzzzzz_0_yyyyzzz_0, g_0_yzzzzz_0_yyyyzzz_1, g_0_yzzzzz_0_yyyzzz_1, g_0_yzzzzz_0_yyyzzzz_0, g_0_yzzzzz_0_yyyzzzz_1, g_0_yzzzzz_0_yyzzzz_1, g_0_yzzzzz_0_yyzzzzz_0, g_0_yzzzzz_0_yyzzzzz_1, g_0_yzzzzz_0_yzzzzz_1, g_0_yzzzzz_0_yzzzzzz_0, g_0_yzzzzz_0_yzzzzzz_1, g_0_yzzzzz_0_zzzzzzz_0, g_0_yzzzzz_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzz_0_xxxxxxx_0[i] = g_0_xzzzzz_0_xxxxxxx_0[i] * pb_y + g_0_xzzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxxxy_0[i] * pb_x + g_0_yzzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxxxz_0[i] = g_0_xzzzzz_0_xxxxxxz_0[i] * pb_y + g_0_xzzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxxyy_0[i] * pb_x + g_0_yzzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxxyz_0[i] * pb_x + g_0_yzzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxxzz_0[i] = g_0_xzzzzz_0_xxxxxzz_0[i] * pb_y + g_0_xzzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxyyy_0[i] * pb_x + g_0_yzzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxyyz_0[i] * pb_x + g_0_yzzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxyzz_0[i] * pb_x + g_0_yzzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxzzz_0[i] = g_0_xzzzzz_0_xxxxzzz_0[i] * pb_y + g_0_xzzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyyyy_0[i] * pb_x + g_0_yzzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyyyz_0[i] * pb_x + g_0_yzzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyyzz_0[i] * pb_x + g_0_yzzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyzzz_0[i] * pb_x + g_0_yzzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxzzzz_0[i] = g_0_xzzzzz_0_xxxzzzz_0[i] * pb_y + g_0_xzzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyyyy_0[i] * pb_x + g_0_yzzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyyyz_0[i] * pb_x + g_0_yzzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyyzz_0[i] * pb_x + g_0_yzzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyzzz_0[i] * pb_x + g_0_yzzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyzzzz_0[i] * pb_x + g_0_yzzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxzzzzz_0[i] = g_0_xzzzzz_0_xxzzzzz_0[i] * pb_y + g_0_xzzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xyyyyyy_0[i] = g_0_yzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyyyy_0[i] * pb_x + g_0_yzzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyyyyz_0[i] = g_0_yzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyyyz_0[i] * pb_x + g_0_yzzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyyyzz_0[i] = g_0_yzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyyzz_0[i] * pb_x + g_0_yzzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyyzzz_0[i] = g_0_yzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyzzz_0[i] * pb_x + g_0_yzzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyzzzz_0[i] = g_0_yzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyzzzz_0[i] * pb_x + g_0_yzzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyzzzzz_0[i] = g_0_yzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyzzzzz_0[i] * pb_x + g_0_yzzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xzzzzzz_0[i] = g_0_xzzzzz_0_xzzzzzz_0[i] * pb_y + g_0_xzzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_yyyyyyy_0[i] = g_0_yzzzzz_0_yyyyyyy_0[i] * pb_x + g_0_yzzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyyyyz_0[i] = g_0_yzzzzz_0_yyyyyyz_0[i] * pb_x + g_0_yzzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyyyzz_0[i] = g_0_yzzzzz_0_yyyyyzz_0[i] * pb_x + g_0_yzzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyyzzz_0[i] = g_0_yzzzzz_0_yyyyzzz_0[i] * pb_x + g_0_yzzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyzzzz_0[i] = g_0_yzzzzz_0_yyyzzzz_0[i] * pb_x + g_0_yzzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyzzzzz_0[i] = g_0_yzzzzz_0_yyzzzzz_0[i] * pb_x + g_0_yzzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yzzzzzz_0[i] = g_0_yzzzzz_0_yzzzzzz_0[i] * pb_x + g_0_yzzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_zzzzzzz_0[i] = g_0_yzzzzz_0_zzzzzzz_0[i] * pb_x + g_0_yzzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 972-1008 components of targeted buffer : SKSK

    auto g_0_xzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 972);

    auto g_0_xzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 973);

    auto g_0_xzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 974);

    auto g_0_xzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 975);

    auto g_0_xzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 976);

    auto g_0_xzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 977);

    auto g_0_xzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 978);

    auto g_0_xzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 979);

    auto g_0_xzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 980);

    auto g_0_xzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 981);

    auto g_0_xzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 982);

    auto g_0_xzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 983);

    auto g_0_xzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 984);

    auto g_0_xzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 985);

    auto g_0_xzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 986);

    auto g_0_xzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 987);

    auto g_0_xzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 988);

    auto g_0_xzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 989);

    auto g_0_xzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 990);

    auto g_0_xzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 991);

    auto g_0_xzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 992);

    auto g_0_xzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 993);

    auto g_0_xzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 994);

    auto g_0_xzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 995);

    auto g_0_xzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 996);

    auto g_0_xzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 997);

    auto g_0_xzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 998);

    auto g_0_xzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 999);

    auto g_0_xzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1000);

    auto g_0_xzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1001);

    auto g_0_xzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1002);

    auto g_0_xzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1003);

    auto g_0_xzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1004);

    auto g_0_xzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1005);

    auto g_0_xzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1006);

    auto g_0_xzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1007);

    #pragma omp simd aligned(g_0_xzzzzzz_0_xxxxxxx_0, g_0_xzzzzzz_0_xxxxxxy_0, g_0_xzzzzzz_0_xxxxxxz_0, g_0_xzzzzzz_0_xxxxxyy_0, g_0_xzzzzzz_0_xxxxxyz_0, g_0_xzzzzzz_0_xxxxxzz_0, g_0_xzzzzzz_0_xxxxyyy_0, g_0_xzzzzzz_0_xxxxyyz_0, g_0_xzzzzzz_0_xxxxyzz_0, g_0_xzzzzzz_0_xxxxzzz_0, g_0_xzzzzzz_0_xxxyyyy_0, g_0_xzzzzzz_0_xxxyyyz_0, g_0_xzzzzzz_0_xxxyyzz_0, g_0_xzzzzzz_0_xxxyzzz_0, g_0_xzzzzzz_0_xxxzzzz_0, g_0_xzzzzzz_0_xxyyyyy_0, g_0_xzzzzzz_0_xxyyyyz_0, g_0_xzzzzzz_0_xxyyyzz_0, g_0_xzzzzzz_0_xxyyzzz_0, g_0_xzzzzzz_0_xxyzzzz_0, g_0_xzzzzzz_0_xxzzzzz_0, g_0_xzzzzzz_0_xyyyyyy_0, g_0_xzzzzzz_0_xyyyyyz_0, g_0_xzzzzzz_0_xyyyyzz_0, g_0_xzzzzzz_0_xyyyzzz_0, g_0_xzzzzzz_0_xyyzzzz_0, g_0_xzzzzzz_0_xyzzzzz_0, g_0_xzzzzzz_0_xzzzzzz_0, g_0_xzzzzzz_0_yyyyyyy_0, g_0_xzzzzzz_0_yyyyyyz_0, g_0_xzzzzzz_0_yyyyyzz_0, g_0_xzzzzzz_0_yyyyzzz_0, g_0_xzzzzzz_0_yyyzzzz_0, g_0_xzzzzzz_0_yyzzzzz_0, g_0_xzzzzzz_0_yzzzzzz_0, g_0_xzzzzzz_0_zzzzzzz_0, g_0_zzzzzz_0_xxxxxx_1, g_0_zzzzzz_0_xxxxxxx_0, g_0_zzzzzz_0_xxxxxxx_1, g_0_zzzzzz_0_xxxxxxy_0, g_0_zzzzzz_0_xxxxxxy_1, g_0_zzzzzz_0_xxxxxxz_0, g_0_zzzzzz_0_xxxxxxz_1, g_0_zzzzzz_0_xxxxxy_1, g_0_zzzzzz_0_xxxxxyy_0, g_0_zzzzzz_0_xxxxxyy_1, g_0_zzzzzz_0_xxxxxyz_0, g_0_zzzzzz_0_xxxxxyz_1, g_0_zzzzzz_0_xxxxxz_1, g_0_zzzzzz_0_xxxxxzz_0, g_0_zzzzzz_0_xxxxxzz_1, g_0_zzzzzz_0_xxxxyy_1, g_0_zzzzzz_0_xxxxyyy_0, g_0_zzzzzz_0_xxxxyyy_1, g_0_zzzzzz_0_xxxxyyz_0, g_0_zzzzzz_0_xxxxyyz_1, g_0_zzzzzz_0_xxxxyz_1, g_0_zzzzzz_0_xxxxyzz_0, g_0_zzzzzz_0_xxxxyzz_1, g_0_zzzzzz_0_xxxxzz_1, g_0_zzzzzz_0_xxxxzzz_0, g_0_zzzzzz_0_xxxxzzz_1, g_0_zzzzzz_0_xxxyyy_1, g_0_zzzzzz_0_xxxyyyy_0, g_0_zzzzzz_0_xxxyyyy_1, g_0_zzzzzz_0_xxxyyyz_0, g_0_zzzzzz_0_xxxyyyz_1, g_0_zzzzzz_0_xxxyyz_1, g_0_zzzzzz_0_xxxyyzz_0, g_0_zzzzzz_0_xxxyyzz_1, g_0_zzzzzz_0_xxxyzz_1, g_0_zzzzzz_0_xxxyzzz_0, g_0_zzzzzz_0_xxxyzzz_1, g_0_zzzzzz_0_xxxzzz_1, g_0_zzzzzz_0_xxxzzzz_0, g_0_zzzzzz_0_xxxzzzz_1, g_0_zzzzzz_0_xxyyyy_1, g_0_zzzzzz_0_xxyyyyy_0, g_0_zzzzzz_0_xxyyyyy_1, g_0_zzzzzz_0_xxyyyyz_0, g_0_zzzzzz_0_xxyyyyz_1, g_0_zzzzzz_0_xxyyyz_1, g_0_zzzzzz_0_xxyyyzz_0, g_0_zzzzzz_0_xxyyyzz_1, g_0_zzzzzz_0_xxyyzz_1, g_0_zzzzzz_0_xxyyzzz_0, g_0_zzzzzz_0_xxyyzzz_1, g_0_zzzzzz_0_xxyzzz_1, g_0_zzzzzz_0_xxyzzzz_0, g_0_zzzzzz_0_xxyzzzz_1, g_0_zzzzzz_0_xxzzzz_1, g_0_zzzzzz_0_xxzzzzz_0, g_0_zzzzzz_0_xxzzzzz_1, g_0_zzzzzz_0_xyyyyy_1, g_0_zzzzzz_0_xyyyyyy_0, g_0_zzzzzz_0_xyyyyyy_1, g_0_zzzzzz_0_xyyyyyz_0, g_0_zzzzzz_0_xyyyyyz_1, g_0_zzzzzz_0_xyyyyz_1, g_0_zzzzzz_0_xyyyyzz_0, g_0_zzzzzz_0_xyyyyzz_1, g_0_zzzzzz_0_xyyyzz_1, g_0_zzzzzz_0_xyyyzzz_0, g_0_zzzzzz_0_xyyyzzz_1, g_0_zzzzzz_0_xyyzzz_1, g_0_zzzzzz_0_xyyzzzz_0, g_0_zzzzzz_0_xyyzzzz_1, g_0_zzzzzz_0_xyzzzz_1, g_0_zzzzzz_0_xyzzzzz_0, g_0_zzzzzz_0_xyzzzzz_1, g_0_zzzzzz_0_xzzzzz_1, g_0_zzzzzz_0_xzzzzzz_0, g_0_zzzzzz_0_xzzzzzz_1, g_0_zzzzzz_0_yyyyyy_1, g_0_zzzzzz_0_yyyyyyy_0, g_0_zzzzzz_0_yyyyyyy_1, g_0_zzzzzz_0_yyyyyyz_0, g_0_zzzzzz_0_yyyyyyz_1, g_0_zzzzzz_0_yyyyyz_1, g_0_zzzzzz_0_yyyyyzz_0, g_0_zzzzzz_0_yyyyyzz_1, g_0_zzzzzz_0_yyyyzz_1, g_0_zzzzzz_0_yyyyzzz_0, g_0_zzzzzz_0_yyyyzzz_1, g_0_zzzzzz_0_yyyzzz_1, g_0_zzzzzz_0_yyyzzzz_0, g_0_zzzzzz_0_yyyzzzz_1, g_0_zzzzzz_0_yyzzzz_1, g_0_zzzzzz_0_yyzzzzz_0, g_0_zzzzzz_0_yyzzzzz_1, g_0_zzzzzz_0_yzzzzz_1, g_0_zzzzzz_0_yzzzzzz_0, g_0_zzzzzz_0_yzzzzzz_1, g_0_zzzzzz_0_zzzzzz_1, g_0_zzzzzz_0_zzzzzzz_0, g_0_zzzzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzz_0_xxxxxxx_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxxx_0[i] * pb_x + g_0_zzzzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_zzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxxy_0[i] * pb_x + g_0_zzzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxxxz_0[i] = 6.0 * g_0_zzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxxz_0[i] * pb_x + g_0_zzzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_zzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxyy_0[i] * pb_x + g_0_zzzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_zzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxyz_0[i] * pb_x + g_0_zzzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxxzz_0[i] = 5.0 * g_0_zzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxzz_0[i] * pb_x + g_0_zzzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_zzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyyy_0[i] * pb_x + g_0_zzzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_zzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyyz_0[i] * pb_x + g_0_zzzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_zzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyzz_0[i] * pb_x + g_0_zzzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxzzz_0[i] = 4.0 * g_0_zzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxzzz_0[i] * pb_x + g_0_zzzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_zzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyyy_0[i] * pb_x + g_0_zzzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyyz_0[i] * pb_x + g_0_zzzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_zzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyzz_0[i] * pb_x + g_0_zzzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_zzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyzzz_0[i] * pb_x + g_0_zzzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_zzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxzzzz_0[i] * pb_x + g_0_zzzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_zzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyyy_0[i] * pb_x + g_0_zzzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_zzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyyz_0[i] * pb_x + g_0_zzzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_zzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyzz_0[i] * pb_x + g_0_zzzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyzzz_0[i] * pb_x + g_0_zzzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_zzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzzzz_0[i] * pb_x + g_0_zzzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_zzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxzzzzz_0[i] * pb_x + g_0_zzzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyyyy_0[i] = g_0_zzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyyy_0[i] * pb_x + g_0_zzzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyyyz_0[i] = g_0_zzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyyz_0[i] * pb_x + g_0_zzzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyyzz_0[i] = g_0_zzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyzz_0[i] * pb_x + g_0_zzzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyzzz_0[i] = g_0_zzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyzzz_0[i] * pb_x + g_0_zzzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyzzzz_0[i] = g_0_zzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzzzz_0[i] * pb_x + g_0_zzzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyzzzzz_0[i] = g_0_zzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzzzz_0[i] * pb_x + g_0_zzzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xzzzzzz_0[i] = g_0_zzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzzzzzz_0[i] * pb_x + g_0_zzzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyyyy_0[i] = g_0_zzzzzz_0_yyyyyyy_0[i] * pb_x + g_0_zzzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyyyz_0[i] = g_0_zzzzzz_0_yyyyyyz_0[i] * pb_x + g_0_zzzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyyzz_0[i] = g_0_zzzzzz_0_yyyyyzz_0[i] * pb_x + g_0_zzzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyzzz_0[i] = g_0_zzzzzz_0_yyyyzzz_0[i] * pb_x + g_0_zzzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyzzzz_0[i] = g_0_zzzzzz_0_yyyzzzz_0[i] * pb_x + g_0_zzzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyzzzzz_0[i] = g_0_zzzzzz_0_yyzzzzz_0[i] * pb_x + g_0_zzzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yzzzzzz_0[i] = g_0_zzzzzz_0_yzzzzzz_0[i] * pb_x + g_0_zzzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_zzzzzzz_0[i] = g_0_zzzzzz_0_zzzzzzz_0[i] * pb_x + g_0_zzzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1008-1044 components of targeted buffer : SKSK

    auto g_0_yyyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 1008);

    auto g_0_yyyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1009);

    auto g_0_yyyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 1010);

    auto g_0_yyyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 1011);

    auto g_0_yyyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 1012);

    auto g_0_yyyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 1013);

    auto g_0_yyyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 1014);

    auto g_0_yyyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 1015);

    auto g_0_yyyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 1016);

    auto g_0_yyyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 1017);

    auto g_0_yyyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1018);

    auto g_0_yyyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1019);

    auto g_0_yyyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1020);

    auto g_0_yyyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1021);

    auto g_0_yyyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1022);

    auto g_0_yyyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1023);

    auto g_0_yyyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1024);

    auto g_0_yyyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1025);

    auto g_0_yyyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1026);

    auto g_0_yyyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1027);

    auto g_0_yyyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1028);

    auto g_0_yyyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1029);

    auto g_0_yyyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1030);

    auto g_0_yyyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1031);

    auto g_0_yyyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1032);

    auto g_0_yyyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1033);

    auto g_0_yyyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1034);

    auto g_0_yyyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1035);

    auto g_0_yyyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1036);

    auto g_0_yyyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1037);

    auto g_0_yyyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1038);

    auto g_0_yyyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1039);

    auto g_0_yyyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1040);

    auto g_0_yyyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1041);

    auto g_0_yyyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1042);

    auto g_0_yyyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1043);

    #pragma omp simd aligned(g_0_yyyyy_0_xxxxxxx_0, g_0_yyyyy_0_xxxxxxx_1, g_0_yyyyy_0_xxxxxxy_0, g_0_yyyyy_0_xxxxxxy_1, g_0_yyyyy_0_xxxxxxz_0, g_0_yyyyy_0_xxxxxxz_1, g_0_yyyyy_0_xxxxxyy_0, g_0_yyyyy_0_xxxxxyy_1, g_0_yyyyy_0_xxxxxyz_0, g_0_yyyyy_0_xxxxxyz_1, g_0_yyyyy_0_xxxxxzz_0, g_0_yyyyy_0_xxxxxzz_1, g_0_yyyyy_0_xxxxyyy_0, g_0_yyyyy_0_xxxxyyy_1, g_0_yyyyy_0_xxxxyyz_0, g_0_yyyyy_0_xxxxyyz_1, g_0_yyyyy_0_xxxxyzz_0, g_0_yyyyy_0_xxxxyzz_1, g_0_yyyyy_0_xxxxzzz_0, g_0_yyyyy_0_xxxxzzz_1, g_0_yyyyy_0_xxxyyyy_0, g_0_yyyyy_0_xxxyyyy_1, g_0_yyyyy_0_xxxyyyz_0, g_0_yyyyy_0_xxxyyyz_1, g_0_yyyyy_0_xxxyyzz_0, g_0_yyyyy_0_xxxyyzz_1, g_0_yyyyy_0_xxxyzzz_0, g_0_yyyyy_0_xxxyzzz_1, g_0_yyyyy_0_xxxzzzz_0, g_0_yyyyy_0_xxxzzzz_1, g_0_yyyyy_0_xxyyyyy_0, g_0_yyyyy_0_xxyyyyy_1, g_0_yyyyy_0_xxyyyyz_0, g_0_yyyyy_0_xxyyyyz_1, g_0_yyyyy_0_xxyyyzz_0, g_0_yyyyy_0_xxyyyzz_1, g_0_yyyyy_0_xxyyzzz_0, g_0_yyyyy_0_xxyyzzz_1, g_0_yyyyy_0_xxyzzzz_0, g_0_yyyyy_0_xxyzzzz_1, g_0_yyyyy_0_xxzzzzz_0, g_0_yyyyy_0_xxzzzzz_1, g_0_yyyyy_0_xyyyyyy_0, g_0_yyyyy_0_xyyyyyy_1, g_0_yyyyy_0_xyyyyyz_0, g_0_yyyyy_0_xyyyyyz_1, g_0_yyyyy_0_xyyyyzz_0, g_0_yyyyy_0_xyyyyzz_1, g_0_yyyyy_0_xyyyzzz_0, g_0_yyyyy_0_xyyyzzz_1, g_0_yyyyy_0_xyyzzzz_0, g_0_yyyyy_0_xyyzzzz_1, g_0_yyyyy_0_xyzzzzz_0, g_0_yyyyy_0_xyzzzzz_1, g_0_yyyyy_0_xzzzzzz_0, g_0_yyyyy_0_xzzzzzz_1, g_0_yyyyy_0_yyyyyyy_0, g_0_yyyyy_0_yyyyyyy_1, g_0_yyyyy_0_yyyyyyz_0, g_0_yyyyy_0_yyyyyyz_1, g_0_yyyyy_0_yyyyyzz_0, g_0_yyyyy_0_yyyyyzz_1, g_0_yyyyy_0_yyyyzzz_0, g_0_yyyyy_0_yyyyzzz_1, g_0_yyyyy_0_yyyzzzz_0, g_0_yyyyy_0_yyyzzzz_1, g_0_yyyyy_0_yyzzzzz_0, g_0_yyyyy_0_yyzzzzz_1, g_0_yyyyy_0_yzzzzzz_0, g_0_yyyyy_0_yzzzzzz_1, g_0_yyyyy_0_zzzzzzz_0, g_0_yyyyy_0_zzzzzzz_1, g_0_yyyyyy_0_xxxxxx_1, g_0_yyyyyy_0_xxxxxxx_0, g_0_yyyyyy_0_xxxxxxx_1, g_0_yyyyyy_0_xxxxxxy_0, g_0_yyyyyy_0_xxxxxxy_1, g_0_yyyyyy_0_xxxxxxz_0, g_0_yyyyyy_0_xxxxxxz_1, g_0_yyyyyy_0_xxxxxy_1, g_0_yyyyyy_0_xxxxxyy_0, g_0_yyyyyy_0_xxxxxyy_1, g_0_yyyyyy_0_xxxxxyz_0, g_0_yyyyyy_0_xxxxxyz_1, g_0_yyyyyy_0_xxxxxz_1, g_0_yyyyyy_0_xxxxxzz_0, g_0_yyyyyy_0_xxxxxzz_1, g_0_yyyyyy_0_xxxxyy_1, g_0_yyyyyy_0_xxxxyyy_0, g_0_yyyyyy_0_xxxxyyy_1, g_0_yyyyyy_0_xxxxyyz_0, g_0_yyyyyy_0_xxxxyyz_1, g_0_yyyyyy_0_xxxxyz_1, g_0_yyyyyy_0_xxxxyzz_0, g_0_yyyyyy_0_xxxxyzz_1, g_0_yyyyyy_0_xxxxzz_1, g_0_yyyyyy_0_xxxxzzz_0, g_0_yyyyyy_0_xxxxzzz_1, g_0_yyyyyy_0_xxxyyy_1, g_0_yyyyyy_0_xxxyyyy_0, g_0_yyyyyy_0_xxxyyyy_1, g_0_yyyyyy_0_xxxyyyz_0, g_0_yyyyyy_0_xxxyyyz_1, g_0_yyyyyy_0_xxxyyz_1, g_0_yyyyyy_0_xxxyyzz_0, g_0_yyyyyy_0_xxxyyzz_1, g_0_yyyyyy_0_xxxyzz_1, g_0_yyyyyy_0_xxxyzzz_0, g_0_yyyyyy_0_xxxyzzz_1, g_0_yyyyyy_0_xxxzzz_1, g_0_yyyyyy_0_xxxzzzz_0, g_0_yyyyyy_0_xxxzzzz_1, g_0_yyyyyy_0_xxyyyy_1, g_0_yyyyyy_0_xxyyyyy_0, g_0_yyyyyy_0_xxyyyyy_1, g_0_yyyyyy_0_xxyyyyz_0, g_0_yyyyyy_0_xxyyyyz_1, g_0_yyyyyy_0_xxyyyz_1, g_0_yyyyyy_0_xxyyyzz_0, g_0_yyyyyy_0_xxyyyzz_1, g_0_yyyyyy_0_xxyyzz_1, g_0_yyyyyy_0_xxyyzzz_0, g_0_yyyyyy_0_xxyyzzz_1, g_0_yyyyyy_0_xxyzzz_1, g_0_yyyyyy_0_xxyzzzz_0, g_0_yyyyyy_0_xxyzzzz_1, g_0_yyyyyy_0_xxzzzz_1, g_0_yyyyyy_0_xxzzzzz_0, g_0_yyyyyy_0_xxzzzzz_1, g_0_yyyyyy_0_xyyyyy_1, g_0_yyyyyy_0_xyyyyyy_0, g_0_yyyyyy_0_xyyyyyy_1, g_0_yyyyyy_0_xyyyyyz_0, g_0_yyyyyy_0_xyyyyyz_1, g_0_yyyyyy_0_xyyyyz_1, g_0_yyyyyy_0_xyyyyzz_0, g_0_yyyyyy_0_xyyyyzz_1, g_0_yyyyyy_0_xyyyzz_1, g_0_yyyyyy_0_xyyyzzz_0, g_0_yyyyyy_0_xyyyzzz_1, g_0_yyyyyy_0_xyyzzz_1, g_0_yyyyyy_0_xyyzzzz_0, g_0_yyyyyy_0_xyyzzzz_1, g_0_yyyyyy_0_xyzzzz_1, g_0_yyyyyy_0_xyzzzzz_0, g_0_yyyyyy_0_xyzzzzz_1, g_0_yyyyyy_0_xzzzzz_1, g_0_yyyyyy_0_xzzzzzz_0, g_0_yyyyyy_0_xzzzzzz_1, g_0_yyyyyy_0_yyyyyy_1, g_0_yyyyyy_0_yyyyyyy_0, g_0_yyyyyy_0_yyyyyyy_1, g_0_yyyyyy_0_yyyyyyz_0, g_0_yyyyyy_0_yyyyyyz_1, g_0_yyyyyy_0_yyyyyz_1, g_0_yyyyyy_0_yyyyyzz_0, g_0_yyyyyy_0_yyyyyzz_1, g_0_yyyyyy_0_yyyyzz_1, g_0_yyyyyy_0_yyyyzzz_0, g_0_yyyyyy_0_yyyyzzz_1, g_0_yyyyyy_0_yyyzzz_1, g_0_yyyyyy_0_yyyzzzz_0, g_0_yyyyyy_0_yyyzzzz_1, g_0_yyyyyy_0_yyzzzz_1, g_0_yyyyyy_0_yyzzzzz_0, g_0_yyyyyy_0_yyzzzzz_1, g_0_yyyyyy_0_yzzzzz_1, g_0_yyyyyy_0_yzzzzzz_0, g_0_yyyyyy_0_yzzzzzz_1, g_0_yyyyyy_0_zzzzzz_1, g_0_yyyyyy_0_zzzzzzz_0, g_0_yyyyyy_0_zzzzzzz_1, g_0_yyyyyyy_0_xxxxxxx_0, g_0_yyyyyyy_0_xxxxxxy_0, g_0_yyyyyyy_0_xxxxxxz_0, g_0_yyyyyyy_0_xxxxxyy_0, g_0_yyyyyyy_0_xxxxxyz_0, g_0_yyyyyyy_0_xxxxxzz_0, g_0_yyyyyyy_0_xxxxyyy_0, g_0_yyyyyyy_0_xxxxyyz_0, g_0_yyyyyyy_0_xxxxyzz_0, g_0_yyyyyyy_0_xxxxzzz_0, g_0_yyyyyyy_0_xxxyyyy_0, g_0_yyyyyyy_0_xxxyyyz_0, g_0_yyyyyyy_0_xxxyyzz_0, g_0_yyyyyyy_0_xxxyzzz_0, g_0_yyyyyyy_0_xxxzzzz_0, g_0_yyyyyyy_0_xxyyyyy_0, g_0_yyyyyyy_0_xxyyyyz_0, g_0_yyyyyyy_0_xxyyyzz_0, g_0_yyyyyyy_0_xxyyzzz_0, g_0_yyyyyyy_0_xxyzzzz_0, g_0_yyyyyyy_0_xxzzzzz_0, g_0_yyyyyyy_0_xyyyyyy_0, g_0_yyyyyyy_0_xyyyyyz_0, g_0_yyyyyyy_0_xyyyyzz_0, g_0_yyyyyyy_0_xyyyzzz_0, g_0_yyyyyyy_0_xyyzzzz_0, g_0_yyyyyyy_0_xyzzzzz_0, g_0_yyyyyyy_0_xzzzzzz_0, g_0_yyyyyyy_0_yyyyyyy_0, g_0_yyyyyyy_0_yyyyyyz_0, g_0_yyyyyyy_0_yyyyyzz_0, g_0_yyyyyyy_0_yyyyzzz_0, g_0_yyyyyyy_0_yyyzzzz_0, g_0_yyyyyyy_0_yyzzzzz_0, g_0_yyyyyyy_0_yzzzzzz_0, g_0_yyyyyyy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyy_0_xxxxxxx_0[i] = 6.0 * g_0_yyyyy_0_xxxxxxx_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxxxx_0[i] * pb_y + g_0_yyyyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxxxy_0[i] = 6.0 * g_0_yyyyy_0_xxxxxxy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxxy_0[i] * pb_y + g_0_yyyyyy_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyy_0_xxxxxxz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxxxz_0[i] * pb_y + g_0_yyyyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxxyy_0[i] = 6.0 * g_0_yyyyy_0_xxxxxyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxyy_0[i] * pb_y + g_0_yyyyyy_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxxyz_0[i] = 6.0 * g_0_yyyyy_0_xxxxxyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxyz_0[i] * pb_y + g_0_yyyyyy_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxxzz_0[i] = 6.0 * g_0_yyyyy_0_xxxxxzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxxzz_0[i] * pb_y + g_0_yyyyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxyyy_0[i] = 6.0 * g_0_yyyyy_0_xxxxyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyyy_0[i] * pb_y + g_0_yyyyyy_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxyyz_0[i] = 6.0 * g_0_yyyyy_0_xxxxyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyyz_0[i] * pb_y + g_0_yyyyyy_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxyzz_0[i] = 6.0 * g_0_yyyyy_0_xxxxyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyzz_0[i] * pb_y + g_0_yyyyyy_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxzzz_0[i] = 6.0 * g_0_yyyyy_0_xxxxzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxxzzz_0[i] * pb_y + g_0_yyyyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyyyy_0[i] = 6.0 * g_0_yyyyy_0_xxxyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyyy_0[i] * pb_y + g_0_yyyyyy_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyyyz_0[i] = 6.0 * g_0_yyyyy_0_xxxyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyyz_0[i] * pb_y + g_0_yyyyyy_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyyzz_0[i] = 6.0 * g_0_yyyyy_0_xxxyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyzz_0[i] * pb_y + g_0_yyyyyy_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyzzz_0[i] = 6.0 * g_0_yyyyy_0_xxxyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyzzz_0[i] * pb_y + g_0_yyyyyy_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxzzzz_0[i] = 6.0 * g_0_yyyyy_0_xxxzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxxzzzz_0[i] * pb_y + g_0_yyyyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyyyy_0[i] = 6.0 * g_0_yyyyy_0_xxyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyyy_0[i] * pb_y + g_0_yyyyyy_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyyyz_0[i] = 6.0 * g_0_yyyyy_0_xxyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyyz_0[i] * pb_y + g_0_yyyyyy_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyyzz_0[i] = 6.0 * g_0_yyyyy_0_xxyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyzz_0[i] * pb_y + g_0_yyyyyy_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyzzz_0[i] = 6.0 * g_0_yyyyy_0_xxyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyzzz_0[i] * pb_y + g_0_yyyyyy_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyzzzz_0[i] = 6.0 * g_0_yyyyy_0_xxyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzzzz_0[i] * pb_y + g_0_yyyyyy_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxzzzzz_0[i] = 6.0 * g_0_yyyyy_0_xxzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xxzzzzz_0[i] * pb_y + g_0_yyyyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyyyy_0[i] = 6.0 * g_0_yyyyy_0_xyyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_yyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyyy_0[i] * pb_y + g_0_yyyyyy_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyyyz_0[i] = 6.0 * g_0_yyyyy_0_xyyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyyz_0[i] * pb_y + g_0_yyyyyy_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyyzz_0[i] = 6.0 * g_0_yyyyy_0_xyyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyzz_0[i] * pb_y + g_0_yyyyyy_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyzzz_0[i] = 6.0 * g_0_yyyyy_0_xyyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyzzz_0[i] * pb_y + g_0_yyyyyy_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyzzzz_0[i] = 6.0 * g_0_yyyyy_0_xyyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzzzz_0[i] * pb_y + g_0_yyyyyy_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyzzzzz_0[i] = 6.0 * g_0_yyyyy_0_xyzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzzzz_0[i] * pb_y + g_0_yyyyyy_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xzzzzzz_0[i] = 6.0 * g_0_yyyyy_0_xzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_xzzzzzz_0[i] * pb_y + g_0_yyyyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyyyy_0[i] = 6.0 * g_0_yyyyy_0_yyyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyyyy_1[i] * fti_ab_0 + 7.0 * g_0_yyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyyyy_0[i] * pb_y + g_0_yyyyyy_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyyyz_0[i] = 6.0 * g_0_yyyyy_0_yyyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyyyz_0[i] * pb_y + g_0_yyyyyy_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyyzz_0[i] = 6.0 * g_0_yyyyy_0_yyyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyyzz_0[i] * pb_y + g_0_yyyyyy_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyzzz_0[i] = 6.0 * g_0_yyyyy_0_yyyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyzzz_0[i] * pb_y + g_0_yyyyyy_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyzzzz_0[i] = 6.0 * g_0_yyyyy_0_yyyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyzzzz_0[i] * pb_y + g_0_yyyyyy_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyzzzzz_0[i] = 6.0 * g_0_yyyyy_0_yyzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyzzzzz_0[i] * pb_y + g_0_yyyyyy_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yzzzzzz_0[i] = 6.0 * g_0_yyyyy_0_yzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yzzzzzz_0[i] * pb_y + g_0_yyyyyy_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_zzzzzzz_0[i] = 6.0 * g_0_yyyyy_0_zzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyy_0_zzzzzzz_0[i] * pb_y + g_0_yyyyyy_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1044-1080 components of targeted buffer : SKSK

    auto g_0_yyyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 1044);

    auto g_0_yyyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1045);

    auto g_0_yyyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 1046);

    auto g_0_yyyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 1047);

    auto g_0_yyyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 1048);

    auto g_0_yyyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 1049);

    auto g_0_yyyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 1050);

    auto g_0_yyyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 1051);

    auto g_0_yyyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 1052);

    auto g_0_yyyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 1053);

    auto g_0_yyyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1054);

    auto g_0_yyyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1055);

    auto g_0_yyyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1056);

    auto g_0_yyyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1057);

    auto g_0_yyyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1058);

    auto g_0_yyyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1059);

    auto g_0_yyyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1060);

    auto g_0_yyyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1061);

    auto g_0_yyyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1062);

    auto g_0_yyyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1063);

    auto g_0_yyyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1064);

    auto g_0_yyyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1065);

    auto g_0_yyyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1066);

    auto g_0_yyyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1067);

    auto g_0_yyyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1068);

    auto g_0_yyyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1069);

    auto g_0_yyyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1070);

    auto g_0_yyyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1071);

    auto g_0_yyyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1072);

    auto g_0_yyyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1073);

    auto g_0_yyyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1074);

    auto g_0_yyyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1075);

    auto g_0_yyyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1076);

    auto g_0_yyyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1077);

    auto g_0_yyyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1078);

    auto g_0_yyyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1079);

    #pragma omp simd aligned(g_0_yyyyyy_0_xxxxxx_1, g_0_yyyyyy_0_xxxxxxx_0, g_0_yyyyyy_0_xxxxxxx_1, g_0_yyyyyy_0_xxxxxxy_0, g_0_yyyyyy_0_xxxxxxy_1, g_0_yyyyyy_0_xxxxxxz_0, g_0_yyyyyy_0_xxxxxxz_1, g_0_yyyyyy_0_xxxxxy_1, g_0_yyyyyy_0_xxxxxyy_0, g_0_yyyyyy_0_xxxxxyy_1, g_0_yyyyyy_0_xxxxxyz_0, g_0_yyyyyy_0_xxxxxyz_1, g_0_yyyyyy_0_xxxxxz_1, g_0_yyyyyy_0_xxxxxzz_0, g_0_yyyyyy_0_xxxxxzz_1, g_0_yyyyyy_0_xxxxyy_1, g_0_yyyyyy_0_xxxxyyy_0, g_0_yyyyyy_0_xxxxyyy_1, g_0_yyyyyy_0_xxxxyyz_0, g_0_yyyyyy_0_xxxxyyz_1, g_0_yyyyyy_0_xxxxyz_1, g_0_yyyyyy_0_xxxxyzz_0, g_0_yyyyyy_0_xxxxyzz_1, g_0_yyyyyy_0_xxxxzz_1, g_0_yyyyyy_0_xxxxzzz_0, g_0_yyyyyy_0_xxxxzzz_1, g_0_yyyyyy_0_xxxyyy_1, g_0_yyyyyy_0_xxxyyyy_0, g_0_yyyyyy_0_xxxyyyy_1, g_0_yyyyyy_0_xxxyyyz_0, g_0_yyyyyy_0_xxxyyyz_1, g_0_yyyyyy_0_xxxyyz_1, g_0_yyyyyy_0_xxxyyzz_0, g_0_yyyyyy_0_xxxyyzz_1, g_0_yyyyyy_0_xxxyzz_1, g_0_yyyyyy_0_xxxyzzz_0, g_0_yyyyyy_0_xxxyzzz_1, g_0_yyyyyy_0_xxxzzz_1, g_0_yyyyyy_0_xxxzzzz_0, g_0_yyyyyy_0_xxxzzzz_1, g_0_yyyyyy_0_xxyyyy_1, g_0_yyyyyy_0_xxyyyyy_0, g_0_yyyyyy_0_xxyyyyy_1, g_0_yyyyyy_0_xxyyyyz_0, g_0_yyyyyy_0_xxyyyyz_1, g_0_yyyyyy_0_xxyyyz_1, g_0_yyyyyy_0_xxyyyzz_0, g_0_yyyyyy_0_xxyyyzz_1, g_0_yyyyyy_0_xxyyzz_1, g_0_yyyyyy_0_xxyyzzz_0, g_0_yyyyyy_0_xxyyzzz_1, g_0_yyyyyy_0_xxyzzz_1, g_0_yyyyyy_0_xxyzzzz_0, g_0_yyyyyy_0_xxyzzzz_1, g_0_yyyyyy_0_xxzzzz_1, g_0_yyyyyy_0_xxzzzzz_0, g_0_yyyyyy_0_xxzzzzz_1, g_0_yyyyyy_0_xyyyyy_1, g_0_yyyyyy_0_xyyyyyy_0, g_0_yyyyyy_0_xyyyyyy_1, g_0_yyyyyy_0_xyyyyyz_0, g_0_yyyyyy_0_xyyyyyz_1, g_0_yyyyyy_0_xyyyyz_1, g_0_yyyyyy_0_xyyyyzz_0, g_0_yyyyyy_0_xyyyyzz_1, g_0_yyyyyy_0_xyyyzz_1, g_0_yyyyyy_0_xyyyzzz_0, g_0_yyyyyy_0_xyyyzzz_1, g_0_yyyyyy_0_xyyzzz_1, g_0_yyyyyy_0_xyyzzzz_0, g_0_yyyyyy_0_xyyzzzz_1, g_0_yyyyyy_0_xyzzzz_1, g_0_yyyyyy_0_xyzzzzz_0, g_0_yyyyyy_0_xyzzzzz_1, g_0_yyyyyy_0_xzzzzz_1, g_0_yyyyyy_0_xzzzzzz_0, g_0_yyyyyy_0_xzzzzzz_1, g_0_yyyyyy_0_yyyyyy_1, g_0_yyyyyy_0_yyyyyyy_0, g_0_yyyyyy_0_yyyyyyy_1, g_0_yyyyyy_0_yyyyyyz_0, g_0_yyyyyy_0_yyyyyyz_1, g_0_yyyyyy_0_yyyyyz_1, g_0_yyyyyy_0_yyyyyzz_0, g_0_yyyyyy_0_yyyyyzz_1, g_0_yyyyyy_0_yyyyzz_1, g_0_yyyyyy_0_yyyyzzz_0, g_0_yyyyyy_0_yyyyzzz_1, g_0_yyyyyy_0_yyyzzz_1, g_0_yyyyyy_0_yyyzzzz_0, g_0_yyyyyy_0_yyyzzzz_1, g_0_yyyyyy_0_yyzzzz_1, g_0_yyyyyy_0_yyzzzzz_0, g_0_yyyyyy_0_yyzzzzz_1, g_0_yyyyyy_0_yzzzzz_1, g_0_yyyyyy_0_yzzzzzz_0, g_0_yyyyyy_0_yzzzzzz_1, g_0_yyyyyy_0_zzzzzz_1, g_0_yyyyyy_0_zzzzzzz_0, g_0_yyyyyy_0_zzzzzzz_1, g_0_yyyyyyz_0_xxxxxxx_0, g_0_yyyyyyz_0_xxxxxxy_0, g_0_yyyyyyz_0_xxxxxxz_0, g_0_yyyyyyz_0_xxxxxyy_0, g_0_yyyyyyz_0_xxxxxyz_0, g_0_yyyyyyz_0_xxxxxzz_0, g_0_yyyyyyz_0_xxxxyyy_0, g_0_yyyyyyz_0_xxxxyyz_0, g_0_yyyyyyz_0_xxxxyzz_0, g_0_yyyyyyz_0_xxxxzzz_0, g_0_yyyyyyz_0_xxxyyyy_0, g_0_yyyyyyz_0_xxxyyyz_0, g_0_yyyyyyz_0_xxxyyzz_0, g_0_yyyyyyz_0_xxxyzzz_0, g_0_yyyyyyz_0_xxxzzzz_0, g_0_yyyyyyz_0_xxyyyyy_0, g_0_yyyyyyz_0_xxyyyyz_0, g_0_yyyyyyz_0_xxyyyzz_0, g_0_yyyyyyz_0_xxyyzzz_0, g_0_yyyyyyz_0_xxyzzzz_0, g_0_yyyyyyz_0_xxzzzzz_0, g_0_yyyyyyz_0_xyyyyyy_0, g_0_yyyyyyz_0_xyyyyyz_0, g_0_yyyyyyz_0_xyyyyzz_0, g_0_yyyyyyz_0_xyyyzzz_0, g_0_yyyyyyz_0_xyyzzzz_0, g_0_yyyyyyz_0_xyzzzzz_0, g_0_yyyyyyz_0_xzzzzzz_0, g_0_yyyyyyz_0_yyyyyyy_0, g_0_yyyyyyz_0_yyyyyyz_0, g_0_yyyyyyz_0_yyyyyzz_0, g_0_yyyyyyz_0_yyyyzzz_0, g_0_yyyyyyz_0_yyyzzzz_0, g_0_yyyyyyz_0_yyzzzzz_0, g_0_yyyyyyz_0_yzzzzzz_0, g_0_yyyyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyz_0_xxxxxxx_0[i] = g_0_yyyyyy_0_xxxxxxx_0[i] * pb_z + g_0_yyyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxxxy_0[i] = g_0_yyyyyy_0_xxxxxxy_0[i] * pb_z + g_0_yyyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxxxz_0[i] = g_0_yyyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxxz_0[i] * pb_z + g_0_yyyyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxxyy_0[i] = g_0_yyyyyy_0_xxxxxyy_0[i] * pb_z + g_0_yyyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxxyz_0[i] = g_0_yyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxyz_0[i] * pb_z + g_0_yyyyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_yyyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxzz_0[i] * pb_z + g_0_yyyyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxyyy_0[i] = g_0_yyyyyy_0_xxxxyyy_0[i] * pb_z + g_0_yyyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxyyz_0[i] = g_0_yyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyyz_0[i] * pb_z + g_0_yyyyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_yyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyzz_0[i] * pb_z + g_0_yyyyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxzzz_0[i] * pb_z + g_0_yyyyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyyyy_0[i] = g_0_yyyyyy_0_xxxyyyy_0[i] * pb_z + g_0_yyyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyyyz_0[i] = g_0_yyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyyz_0[i] * pb_z + g_0_yyyyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_yyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyzz_0[i] * pb_z + g_0_yyyyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyzzz_0[i] * pb_z + g_0_yyyyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_yyyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxzzzz_0[i] * pb_z + g_0_yyyyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyyyy_0[i] = g_0_yyyyyy_0_xxyyyyy_0[i] * pb_z + g_0_yyyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyyyz_0[i] = g_0_yyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyyz_0[i] * pb_z + g_0_yyyyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyzz_0[i] * pb_z + g_0_yyyyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyzzz_0[i] * pb_z + g_0_yyyyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_yyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzzzz_0[i] * pb_z + g_0_yyyyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_yyyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxzzzzz_0[i] * pb_z + g_0_yyyyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyyyy_0[i] = g_0_yyyyyy_0_xyyyyyy_0[i] * pb_z + g_0_yyyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyyyz_0[i] = g_0_yyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyyz_0[i] * pb_z + g_0_yyyyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_yyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyzz_0[i] * pb_z + g_0_yyyyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_yyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyzzz_0[i] * pb_z + g_0_yyyyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_yyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzzzz_0[i] * pb_z + g_0_yyyyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_yyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzzzz_0[i] * pb_z + g_0_yyyyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_yyyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzzzzzz_0[i] * pb_z + g_0_yyyyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyyyy_0[i] = g_0_yyyyyy_0_yyyyyyy_0[i] * pb_z + g_0_yyyyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyyyz_0[i] = g_0_yyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyyyz_0[i] * pb_z + g_0_yyyyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_yyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyyzz_0[i] * pb_z + g_0_yyyyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_yyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyzzz_0[i] * pb_z + g_0_yyyyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_yyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyzzzz_0[i] * pb_z + g_0_yyyyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_yyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyzzzzz_0[i] * pb_z + g_0_yyyyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_yyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yzzzzzz_0[i] * pb_z + g_0_yyyyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_zzzzzzz_0[i] * pb_z + g_0_yyyyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 1080-1116 components of targeted buffer : SKSK

    auto g_0_yyyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 1080);

    auto g_0_yyyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1081);

    auto g_0_yyyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 1082);

    auto g_0_yyyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 1083);

    auto g_0_yyyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 1084);

    auto g_0_yyyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 1085);

    auto g_0_yyyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 1086);

    auto g_0_yyyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 1087);

    auto g_0_yyyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 1088);

    auto g_0_yyyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 1089);

    auto g_0_yyyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1090);

    auto g_0_yyyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1091);

    auto g_0_yyyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1092);

    auto g_0_yyyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1093);

    auto g_0_yyyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1094);

    auto g_0_yyyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1095);

    auto g_0_yyyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1096);

    auto g_0_yyyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1097);

    auto g_0_yyyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1098);

    auto g_0_yyyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1099);

    auto g_0_yyyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1100);

    auto g_0_yyyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1101);

    auto g_0_yyyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1102);

    auto g_0_yyyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1103);

    auto g_0_yyyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1104);

    auto g_0_yyyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1105);

    auto g_0_yyyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1106);

    auto g_0_yyyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1107);

    auto g_0_yyyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1108);

    auto g_0_yyyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1109);

    auto g_0_yyyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1110);

    auto g_0_yyyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1111);

    auto g_0_yyyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1112);

    auto g_0_yyyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1113);

    auto g_0_yyyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1114);

    auto g_0_yyyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1115);

    #pragma omp simd aligned(g_0_yyyyy_0_xxxxxxy_0, g_0_yyyyy_0_xxxxxxy_1, g_0_yyyyy_0_xxxxxyy_0, g_0_yyyyy_0_xxxxxyy_1, g_0_yyyyy_0_xxxxyyy_0, g_0_yyyyy_0_xxxxyyy_1, g_0_yyyyy_0_xxxyyyy_0, g_0_yyyyy_0_xxxyyyy_1, g_0_yyyyy_0_xxyyyyy_0, g_0_yyyyy_0_xxyyyyy_1, g_0_yyyyy_0_xyyyyyy_0, g_0_yyyyy_0_xyyyyyy_1, g_0_yyyyy_0_yyyyyyy_0, g_0_yyyyy_0_yyyyyyy_1, g_0_yyyyyz_0_xxxxxxy_0, g_0_yyyyyz_0_xxxxxxy_1, g_0_yyyyyz_0_xxxxxyy_0, g_0_yyyyyz_0_xxxxxyy_1, g_0_yyyyyz_0_xxxxyyy_0, g_0_yyyyyz_0_xxxxyyy_1, g_0_yyyyyz_0_xxxyyyy_0, g_0_yyyyyz_0_xxxyyyy_1, g_0_yyyyyz_0_xxyyyyy_0, g_0_yyyyyz_0_xxyyyyy_1, g_0_yyyyyz_0_xyyyyyy_0, g_0_yyyyyz_0_xyyyyyy_1, g_0_yyyyyz_0_yyyyyyy_0, g_0_yyyyyz_0_yyyyyyy_1, g_0_yyyyyzz_0_xxxxxxx_0, g_0_yyyyyzz_0_xxxxxxy_0, g_0_yyyyyzz_0_xxxxxxz_0, g_0_yyyyyzz_0_xxxxxyy_0, g_0_yyyyyzz_0_xxxxxyz_0, g_0_yyyyyzz_0_xxxxxzz_0, g_0_yyyyyzz_0_xxxxyyy_0, g_0_yyyyyzz_0_xxxxyyz_0, g_0_yyyyyzz_0_xxxxyzz_0, g_0_yyyyyzz_0_xxxxzzz_0, g_0_yyyyyzz_0_xxxyyyy_0, g_0_yyyyyzz_0_xxxyyyz_0, g_0_yyyyyzz_0_xxxyyzz_0, g_0_yyyyyzz_0_xxxyzzz_0, g_0_yyyyyzz_0_xxxzzzz_0, g_0_yyyyyzz_0_xxyyyyy_0, g_0_yyyyyzz_0_xxyyyyz_0, g_0_yyyyyzz_0_xxyyyzz_0, g_0_yyyyyzz_0_xxyyzzz_0, g_0_yyyyyzz_0_xxyzzzz_0, g_0_yyyyyzz_0_xxzzzzz_0, g_0_yyyyyzz_0_xyyyyyy_0, g_0_yyyyyzz_0_xyyyyyz_0, g_0_yyyyyzz_0_xyyyyzz_0, g_0_yyyyyzz_0_xyyyzzz_0, g_0_yyyyyzz_0_xyyzzzz_0, g_0_yyyyyzz_0_xyzzzzz_0, g_0_yyyyyzz_0_xzzzzzz_0, g_0_yyyyyzz_0_yyyyyyy_0, g_0_yyyyyzz_0_yyyyyyz_0, g_0_yyyyyzz_0_yyyyyzz_0, g_0_yyyyyzz_0_yyyyzzz_0, g_0_yyyyyzz_0_yyyzzzz_0, g_0_yyyyyzz_0_yyzzzzz_0, g_0_yyyyyzz_0_yzzzzzz_0, g_0_yyyyyzz_0_zzzzzzz_0, g_0_yyyyzz_0_xxxxxxx_0, g_0_yyyyzz_0_xxxxxxx_1, g_0_yyyyzz_0_xxxxxxz_0, g_0_yyyyzz_0_xxxxxxz_1, g_0_yyyyzz_0_xxxxxyz_0, g_0_yyyyzz_0_xxxxxyz_1, g_0_yyyyzz_0_xxxxxz_1, g_0_yyyyzz_0_xxxxxzz_0, g_0_yyyyzz_0_xxxxxzz_1, g_0_yyyyzz_0_xxxxyyz_0, g_0_yyyyzz_0_xxxxyyz_1, g_0_yyyyzz_0_xxxxyz_1, g_0_yyyyzz_0_xxxxyzz_0, g_0_yyyyzz_0_xxxxyzz_1, g_0_yyyyzz_0_xxxxzz_1, g_0_yyyyzz_0_xxxxzzz_0, g_0_yyyyzz_0_xxxxzzz_1, g_0_yyyyzz_0_xxxyyyz_0, g_0_yyyyzz_0_xxxyyyz_1, g_0_yyyyzz_0_xxxyyz_1, g_0_yyyyzz_0_xxxyyzz_0, g_0_yyyyzz_0_xxxyyzz_1, g_0_yyyyzz_0_xxxyzz_1, g_0_yyyyzz_0_xxxyzzz_0, g_0_yyyyzz_0_xxxyzzz_1, g_0_yyyyzz_0_xxxzzz_1, g_0_yyyyzz_0_xxxzzzz_0, g_0_yyyyzz_0_xxxzzzz_1, g_0_yyyyzz_0_xxyyyyz_0, g_0_yyyyzz_0_xxyyyyz_1, g_0_yyyyzz_0_xxyyyz_1, g_0_yyyyzz_0_xxyyyzz_0, g_0_yyyyzz_0_xxyyyzz_1, g_0_yyyyzz_0_xxyyzz_1, g_0_yyyyzz_0_xxyyzzz_0, g_0_yyyyzz_0_xxyyzzz_1, g_0_yyyyzz_0_xxyzzz_1, g_0_yyyyzz_0_xxyzzzz_0, g_0_yyyyzz_0_xxyzzzz_1, g_0_yyyyzz_0_xxzzzz_1, g_0_yyyyzz_0_xxzzzzz_0, g_0_yyyyzz_0_xxzzzzz_1, g_0_yyyyzz_0_xyyyyyz_0, g_0_yyyyzz_0_xyyyyyz_1, g_0_yyyyzz_0_xyyyyz_1, g_0_yyyyzz_0_xyyyyzz_0, g_0_yyyyzz_0_xyyyyzz_1, g_0_yyyyzz_0_xyyyzz_1, g_0_yyyyzz_0_xyyyzzz_0, g_0_yyyyzz_0_xyyyzzz_1, g_0_yyyyzz_0_xyyzzz_1, g_0_yyyyzz_0_xyyzzzz_0, g_0_yyyyzz_0_xyyzzzz_1, g_0_yyyyzz_0_xyzzzz_1, g_0_yyyyzz_0_xyzzzzz_0, g_0_yyyyzz_0_xyzzzzz_1, g_0_yyyyzz_0_xzzzzz_1, g_0_yyyyzz_0_xzzzzzz_0, g_0_yyyyzz_0_xzzzzzz_1, g_0_yyyyzz_0_yyyyyyz_0, g_0_yyyyzz_0_yyyyyyz_1, g_0_yyyyzz_0_yyyyyz_1, g_0_yyyyzz_0_yyyyyzz_0, g_0_yyyyzz_0_yyyyyzz_1, g_0_yyyyzz_0_yyyyzz_1, g_0_yyyyzz_0_yyyyzzz_0, g_0_yyyyzz_0_yyyyzzz_1, g_0_yyyyzz_0_yyyzzz_1, g_0_yyyyzz_0_yyyzzzz_0, g_0_yyyyzz_0_yyyzzzz_1, g_0_yyyyzz_0_yyzzzz_1, g_0_yyyyzz_0_yyzzzzz_0, g_0_yyyyzz_0_yyzzzzz_1, g_0_yyyyzz_0_yzzzzz_1, g_0_yyyyzz_0_yzzzzzz_0, g_0_yyyyzz_0_yzzzzzz_1, g_0_yyyyzz_0_zzzzzz_1, g_0_yyyyzz_0_zzzzzzz_0, g_0_yyyyzz_0_zzzzzzz_1, g_0_yyyzz_0_xxxxxxx_0, g_0_yyyzz_0_xxxxxxx_1, g_0_yyyzz_0_xxxxxxz_0, g_0_yyyzz_0_xxxxxxz_1, g_0_yyyzz_0_xxxxxyz_0, g_0_yyyzz_0_xxxxxyz_1, g_0_yyyzz_0_xxxxxzz_0, g_0_yyyzz_0_xxxxxzz_1, g_0_yyyzz_0_xxxxyyz_0, g_0_yyyzz_0_xxxxyyz_1, g_0_yyyzz_0_xxxxyzz_0, g_0_yyyzz_0_xxxxyzz_1, g_0_yyyzz_0_xxxxzzz_0, g_0_yyyzz_0_xxxxzzz_1, g_0_yyyzz_0_xxxyyyz_0, g_0_yyyzz_0_xxxyyyz_1, g_0_yyyzz_0_xxxyyzz_0, g_0_yyyzz_0_xxxyyzz_1, g_0_yyyzz_0_xxxyzzz_0, g_0_yyyzz_0_xxxyzzz_1, g_0_yyyzz_0_xxxzzzz_0, g_0_yyyzz_0_xxxzzzz_1, g_0_yyyzz_0_xxyyyyz_0, g_0_yyyzz_0_xxyyyyz_1, g_0_yyyzz_0_xxyyyzz_0, g_0_yyyzz_0_xxyyyzz_1, g_0_yyyzz_0_xxyyzzz_0, g_0_yyyzz_0_xxyyzzz_1, g_0_yyyzz_0_xxyzzzz_0, g_0_yyyzz_0_xxyzzzz_1, g_0_yyyzz_0_xxzzzzz_0, g_0_yyyzz_0_xxzzzzz_1, g_0_yyyzz_0_xyyyyyz_0, g_0_yyyzz_0_xyyyyyz_1, g_0_yyyzz_0_xyyyyzz_0, g_0_yyyzz_0_xyyyyzz_1, g_0_yyyzz_0_xyyyzzz_0, g_0_yyyzz_0_xyyyzzz_1, g_0_yyyzz_0_xyyzzzz_0, g_0_yyyzz_0_xyyzzzz_1, g_0_yyyzz_0_xyzzzzz_0, g_0_yyyzz_0_xyzzzzz_1, g_0_yyyzz_0_xzzzzzz_0, g_0_yyyzz_0_xzzzzzz_1, g_0_yyyzz_0_yyyyyyz_0, g_0_yyyzz_0_yyyyyyz_1, g_0_yyyzz_0_yyyyyzz_0, g_0_yyyzz_0_yyyyyzz_1, g_0_yyyzz_0_yyyyzzz_0, g_0_yyyzz_0_yyyyzzz_1, g_0_yyyzz_0_yyyzzzz_0, g_0_yyyzz_0_yyyzzzz_1, g_0_yyyzz_0_yyzzzzz_0, g_0_yyyzz_0_yyzzzzz_1, g_0_yyyzz_0_yzzzzzz_0, g_0_yyyzz_0_yzzzzzz_1, g_0_yyyzz_0_zzzzzzz_0, g_0_yyyzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzz_0_xxxxxxx_0[i] = 4.0 * g_0_yyyzz_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxxxx_0[i] * pb_y + g_0_yyyyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxxxy_0[i] = g_0_yyyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxxxxy_0[i] * pb_z + g_0_yyyyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxxxxz_0[i] = 4.0 * g_0_yyyzz_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxxxz_0[i] * pb_y + g_0_yyyyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxxyy_0[i] = g_0_yyyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxxxyy_0[i] * pb_z + g_0_yyyyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxxxyz_0[i] = 4.0 * g_0_yyyzz_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxyz_0[i] * pb_y + g_0_yyyyzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxxzz_0[i] = 4.0 * g_0_yyyzz_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxxzz_0[i] * pb_y + g_0_yyyyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxyyy_0[i] = g_0_yyyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxxyyy_0[i] * pb_z + g_0_yyyyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyzz_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxyyz_0[i] * pb_y + g_0_yyyyzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyzz_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxyzz_0[i] * pb_y + g_0_yyyyzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyzz_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxzzz_0[i] * pb_y + g_0_yyyyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxyyyy_0[i] = g_0_yyyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxyyyy_0[i] * pb_z + g_0_yyyyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxyyyz_0[i] = 4.0 * g_0_yyyzz_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyyyz_0[i] * pb_y + g_0_yyyyzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxyyzz_0[i] = 4.0 * g_0_yyyzz_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyyzz_0[i] * pb_y + g_0_yyyyzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxyzzz_0[i] = 4.0 * g_0_yyyzz_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyzzz_0[i] * pb_y + g_0_yyyyzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxzzzz_0[i] = 4.0 * g_0_yyyzz_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxzzzz_0[i] * pb_y + g_0_yyyyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyyyyy_0[i] = g_0_yyyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxyyyyy_0[i] * pb_z + g_0_yyyyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxyyyyz_0[i] = 4.0 * g_0_yyyzz_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyyyz_0[i] * pb_y + g_0_yyyyzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyyyzz_0[i] = 4.0 * g_0_yyyzz_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyyzz_0[i] * pb_y + g_0_yyyyzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyyzzz_0[i] = 4.0 * g_0_yyyzz_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyzzz_0[i] * pb_y + g_0_yyyyzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyzzzz_0[i] = 4.0 * g_0_yyyzz_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyzzzz_0[i] * pb_y + g_0_yyyyzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxzzzzz_0[i] = 4.0 * g_0_yyyzz_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxzzzzz_0[i] * pb_y + g_0_yyyyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyyyyy_0[i] = g_0_yyyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xyyyyyy_0[i] * pb_z + g_0_yyyyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xyyyyyz_0[i] = 4.0 * g_0_yyyzz_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyyyz_0[i] * pb_y + g_0_yyyyzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyyyzz_0[i] = 4.0 * g_0_yyyzz_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyyzz_0[i] * pb_y + g_0_yyyyzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyyzzz_0[i] = 4.0 * g_0_yyyzz_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyzzz_0[i] * pb_y + g_0_yyyyzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyzzzz_0[i] = 4.0 * g_0_yyyzz_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyzzzz_0[i] * pb_y + g_0_yyyyzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyzzzzz_0[i] = 4.0 * g_0_yyyzz_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyzzzzz_0[i] * pb_y + g_0_yyyyzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xzzzzzz_0[i] = 4.0 * g_0_yyyzz_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_xzzzzzz_0[i] * pb_y + g_0_yyyyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyyyyy_0[i] = g_0_yyyyy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_yyyyyyy_0[i] * pb_z + g_0_yyyyyz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_yyyyyyz_0[i] = 4.0 * g_0_yyyzz_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyyyyz_0[i] * pb_y + g_0_yyyyzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyyyzz_0[i] = 4.0 * g_0_yyyzz_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyyyzz_0[i] * pb_y + g_0_yyyyzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyyzzz_0[i] = 4.0 * g_0_yyyzz_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyyzzz_0[i] * pb_y + g_0_yyyyzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyzzzz_0[i] = 4.0 * g_0_yyyzz_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyzzzz_0[i] * pb_y + g_0_yyyyzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyzzzzz_0[i] = 4.0 * g_0_yyyzz_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyzzzzz_0[i] * pb_y + g_0_yyyyzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yzzzzzz_0[i] = 4.0 * g_0_yyyzz_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yzzzzzz_0[i] * pb_y + g_0_yyyyzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_zzzzzzz_0[i] = 4.0 * g_0_yyyzz_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyyzz_0_zzzzzzz_0[i] * pb_y + g_0_yyyyzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1116-1152 components of targeted buffer : SKSK

    auto g_0_yyyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 1116);

    auto g_0_yyyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1117);

    auto g_0_yyyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 1118);

    auto g_0_yyyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 1119);

    auto g_0_yyyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 1120);

    auto g_0_yyyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 1121);

    auto g_0_yyyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 1122);

    auto g_0_yyyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 1123);

    auto g_0_yyyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 1124);

    auto g_0_yyyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 1125);

    auto g_0_yyyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1126);

    auto g_0_yyyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1127);

    auto g_0_yyyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1128);

    auto g_0_yyyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1129);

    auto g_0_yyyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1130);

    auto g_0_yyyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1131);

    auto g_0_yyyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1132);

    auto g_0_yyyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1133);

    auto g_0_yyyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1134);

    auto g_0_yyyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1135);

    auto g_0_yyyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1136);

    auto g_0_yyyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1137);

    auto g_0_yyyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1138);

    auto g_0_yyyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1139);

    auto g_0_yyyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1140);

    auto g_0_yyyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1141);

    auto g_0_yyyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1142);

    auto g_0_yyyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1143);

    auto g_0_yyyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1144);

    auto g_0_yyyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1145);

    auto g_0_yyyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1146);

    auto g_0_yyyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1147);

    auto g_0_yyyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1148);

    auto g_0_yyyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1149);

    auto g_0_yyyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1150);

    auto g_0_yyyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1151);

    #pragma omp simd aligned(g_0_yyyyz_0_xxxxxxy_0, g_0_yyyyz_0_xxxxxxy_1, g_0_yyyyz_0_xxxxxyy_0, g_0_yyyyz_0_xxxxxyy_1, g_0_yyyyz_0_xxxxyyy_0, g_0_yyyyz_0_xxxxyyy_1, g_0_yyyyz_0_xxxyyyy_0, g_0_yyyyz_0_xxxyyyy_1, g_0_yyyyz_0_xxyyyyy_0, g_0_yyyyz_0_xxyyyyy_1, g_0_yyyyz_0_xyyyyyy_0, g_0_yyyyz_0_xyyyyyy_1, g_0_yyyyz_0_yyyyyyy_0, g_0_yyyyz_0_yyyyyyy_1, g_0_yyyyzz_0_xxxxxxy_0, g_0_yyyyzz_0_xxxxxxy_1, g_0_yyyyzz_0_xxxxxyy_0, g_0_yyyyzz_0_xxxxxyy_1, g_0_yyyyzz_0_xxxxyyy_0, g_0_yyyyzz_0_xxxxyyy_1, g_0_yyyyzz_0_xxxyyyy_0, g_0_yyyyzz_0_xxxyyyy_1, g_0_yyyyzz_0_xxyyyyy_0, g_0_yyyyzz_0_xxyyyyy_1, g_0_yyyyzz_0_xyyyyyy_0, g_0_yyyyzz_0_xyyyyyy_1, g_0_yyyyzz_0_yyyyyyy_0, g_0_yyyyzz_0_yyyyyyy_1, g_0_yyyyzzz_0_xxxxxxx_0, g_0_yyyyzzz_0_xxxxxxy_0, g_0_yyyyzzz_0_xxxxxxz_0, g_0_yyyyzzz_0_xxxxxyy_0, g_0_yyyyzzz_0_xxxxxyz_0, g_0_yyyyzzz_0_xxxxxzz_0, g_0_yyyyzzz_0_xxxxyyy_0, g_0_yyyyzzz_0_xxxxyyz_0, g_0_yyyyzzz_0_xxxxyzz_0, g_0_yyyyzzz_0_xxxxzzz_0, g_0_yyyyzzz_0_xxxyyyy_0, g_0_yyyyzzz_0_xxxyyyz_0, g_0_yyyyzzz_0_xxxyyzz_0, g_0_yyyyzzz_0_xxxyzzz_0, g_0_yyyyzzz_0_xxxzzzz_0, g_0_yyyyzzz_0_xxyyyyy_0, g_0_yyyyzzz_0_xxyyyyz_0, g_0_yyyyzzz_0_xxyyyzz_0, g_0_yyyyzzz_0_xxyyzzz_0, g_0_yyyyzzz_0_xxyzzzz_0, g_0_yyyyzzz_0_xxzzzzz_0, g_0_yyyyzzz_0_xyyyyyy_0, g_0_yyyyzzz_0_xyyyyyz_0, g_0_yyyyzzz_0_xyyyyzz_0, g_0_yyyyzzz_0_xyyyzzz_0, g_0_yyyyzzz_0_xyyzzzz_0, g_0_yyyyzzz_0_xyzzzzz_0, g_0_yyyyzzz_0_xzzzzzz_0, g_0_yyyyzzz_0_yyyyyyy_0, g_0_yyyyzzz_0_yyyyyyz_0, g_0_yyyyzzz_0_yyyyyzz_0, g_0_yyyyzzz_0_yyyyzzz_0, g_0_yyyyzzz_0_yyyzzzz_0, g_0_yyyyzzz_0_yyzzzzz_0, g_0_yyyyzzz_0_yzzzzzz_0, g_0_yyyyzzz_0_zzzzzzz_0, g_0_yyyzzz_0_xxxxxxx_0, g_0_yyyzzz_0_xxxxxxx_1, g_0_yyyzzz_0_xxxxxxz_0, g_0_yyyzzz_0_xxxxxxz_1, g_0_yyyzzz_0_xxxxxyz_0, g_0_yyyzzz_0_xxxxxyz_1, g_0_yyyzzz_0_xxxxxz_1, g_0_yyyzzz_0_xxxxxzz_0, g_0_yyyzzz_0_xxxxxzz_1, g_0_yyyzzz_0_xxxxyyz_0, g_0_yyyzzz_0_xxxxyyz_1, g_0_yyyzzz_0_xxxxyz_1, g_0_yyyzzz_0_xxxxyzz_0, g_0_yyyzzz_0_xxxxyzz_1, g_0_yyyzzz_0_xxxxzz_1, g_0_yyyzzz_0_xxxxzzz_0, g_0_yyyzzz_0_xxxxzzz_1, g_0_yyyzzz_0_xxxyyyz_0, g_0_yyyzzz_0_xxxyyyz_1, g_0_yyyzzz_0_xxxyyz_1, g_0_yyyzzz_0_xxxyyzz_0, g_0_yyyzzz_0_xxxyyzz_1, g_0_yyyzzz_0_xxxyzz_1, g_0_yyyzzz_0_xxxyzzz_0, g_0_yyyzzz_0_xxxyzzz_1, g_0_yyyzzz_0_xxxzzz_1, g_0_yyyzzz_0_xxxzzzz_0, g_0_yyyzzz_0_xxxzzzz_1, g_0_yyyzzz_0_xxyyyyz_0, g_0_yyyzzz_0_xxyyyyz_1, g_0_yyyzzz_0_xxyyyz_1, g_0_yyyzzz_0_xxyyyzz_0, g_0_yyyzzz_0_xxyyyzz_1, g_0_yyyzzz_0_xxyyzz_1, g_0_yyyzzz_0_xxyyzzz_0, g_0_yyyzzz_0_xxyyzzz_1, g_0_yyyzzz_0_xxyzzz_1, g_0_yyyzzz_0_xxyzzzz_0, g_0_yyyzzz_0_xxyzzzz_1, g_0_yyyzzz_0_xxzzzz_1, g_0_yyyzzz_0_xxzzzzz_0, g_0_yyyzzz_0_xxzzzzz_1, g_0_yyyzzz_0_xyyyyyz_0, g_0_yyyzzz_0_xyyyyyz_1, g_0_yyyzzz_0_xyyyyz_1, g_0_yyyzzz_0_xyyyyzz_0, g_0_yyyzzz_0_xyyyyzz_1, g_0_yyyzzz_0_xyyyzz_1, g_0_yyyzzz_0_xyyyzzz_0, g_0_yyyzzz_0_xyyyzzz_1, g_0_yyyzzz_0_xyyzzz_1, g_0_yyyzzz_0_xyyzzzz_0, g_0_yyyzzz_0_xyyzzzz_1, g_0_yyyzzz_0_xyzzzz_1, g_0_yyyzzz_0_xyzzzzz_0, g_0_yyyzzz_0_xyzzzzz_1, g_0_yyyzzz_0_xzzzzz_1, g_0_yyyzzz_0_xzzzzzz_0, g_0_yyyzzz_0_xzzzzzz_1, g_0_yyyzzz_0_yyyyyyz_0, g_0_yyyzzz_0_yyyyyyz_1, g_0_yyyzzz_0_yyyyyz_1, g_0_yyyzzz_0_yyyyyzz_0, g_0_yyyzzz_0_yyyyyzz_1, g_0_yyyzzz_0_yyyyzz_1, g_0_yyyzzz_0_yyyyzzz_0, g_0_yyyzzz_0_yyyyzzz_1, g_0_yyyzzz_0_yyyzzz_1, g_0_yyyzzz_0_yyyzzzz_0, g_0_yyyzzz_0_yyyzzzz_1, g_0_yyyzzz_0_yyzzzz_1, g_0_yyyzzz_0_yyzzzzz_0, g_0_yyyzzz_0_yyzzzzz_1, g_0_yyyzzz_0_yzzzzz_1, g_0_yyyzzz_0_yzzzzzz_0, g_0_yyyzzz_0_yzzzzzz_1, g_0_yyyzzz_0_zzzzzz_1, g_0_yyyzzz_0_zzzzzzz_0, g_0_yyyzzz_0_zzzzzzz_1, g_0_yyzzz_0_xxxxxxx_0, g_0_yyzzz_0_xxxxxxx_1, g_0_yyzzz_0_xxxxxxz_0, g_0_yyzzz_0_xxxxxxz_1, g_0_yyzzz_0_xxxxxyz_0, g_0_yyzzz_0_xxxxxyz_1, g_0_yyzzz_0_xxxxxzz_0, g_0_yyzzz_0_xxxxxzz_1, g_0_yyzzz_0_xxxxyyz_0, g_0_yyzzz_0_xxxxyyz_1, g_0_yyzzz_0_xxxxyzz_0, g_0_yyzzz_0_xxxxyzz_1, g_0_yyzzz_0_xxxxzzz_0, g_0_yyzzz_0_xxxxzzz_1, g_0_yyzzz_0_xxxyyyz_0, g_0_yyzzz_0_xxxyyyz_1, g_0_yyzzz_0_xxxyyzz_0, g_0_yyzzz_0_xxxyyzz_1, g_0_yyzzz_0_xxxyzzz_0, g_0_yyzzz_0_xxxyzzz_1, g_0_yyzzz_0_xxxzzzz_0, g_0_yyzzz_0_xxxzzzz_1, g_0_yyzzz_0_xxyyyyz_0, g_0_yyzzz_0_xxyyyyz_1, g_0_yyzzz_0_xxyyyzz_0, g_0_yyzzz_0_xxyyyzz_1, g_0_yyzzz_0_xxyyzzz_0, g_0_yyzzz_0_xxyyzzz_1, g_0_yyzzz_0_xxyzzzz_0, g_0_yyzzz_0_xxyzzzz_1, g_0_yyzzz_0_xxzzzzz_0, g_0_yyzzz_0_xxzzzzz_1, g_0_yyzzz_0_xyyyyyz_0, g_0_yyzzz_0_xyyyyyz_1, g_0_yyzzz_0_xyyyyzz_0, g_0_yyzzz_0_xyyyyzz_1, g_0_yyzzz_0_xyyyzzz_0, g_0_yyzzz_0_xyyyzzz_1, g_0_yyzzz_0_xyyzzzz_0, g_0_yyzzz_0_xyyzzzz_1, g_0_yyzzz_0_xyzzzzz_0, g_0_yyzzz_0_xyzzzzz_1, g_0_yyzzz_0_xzzzzzz_0, g_0_yyzzz_0_xzzzzzz_1, g_0_yyzzz_0_yyyyyyz_0, g_0_yyzzz_0_yyyyyyz_1, g_0_yyzzz_0_yyyyyzz_0, g_0_yyzzz_0_yyyyyzz_1, g_0_yyzzz_0_yyyyzzz_0, g_0_yyzzz_0_yyyyzzz_1, g_0_yyzzz_0_yyyzzzz_0, g_0_yyzzz_0_yyyzzzz_1, g_0_yyzzz_0_yyzzzzz_0, g_0_yyzzz_0_yyzzzzz_1, g_0_yyzzz_0_yzzzzzz_0, g_0_yyzzz_0_yzzzzzz_1, g_0_yyzzz_0_zzzzzzz_0, g_0_yyzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzz_0_xxxxxxx_0[i] = 3.0 * g_0_yyzzz_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxxxx_0[i] * pb_y + g_0_yyyzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxxxy_0[i] = 2.0 * g_0_yyyyz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxxxy_0[i] * pb_z + g_0_yyyyzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxxxxz_0[i] = 3.0 * g_0_yyzzz_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxxxz_0[i] * pb_y + g_0_yyyzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_yyyyz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxxyy_0[i] * pb_z + g_0_yyyyzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxxxyz_0[i] = 3.0 * g_0_yyzzz_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxyz_0[i] * pb_y + g_0_yyyzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxxzz_0[i] = 3.0 * g_0_yyzzz_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxxzz_0[i] * pb_y + g_0_yyyzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxyyy_0[i] = 2.0 * g_0_yyyyz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxxyyy_0[i] * pb_z + g_0_yyyyzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxxyyz_0[i] = 3.0 * g_0_yyzzz_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxyyz_0[i] * pb_y + g_0_yyyzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxyzz_0[i] = 3.0 * g_0_yyzzz_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxyzz_0[i] * pb_y + g_0_yyyzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxzzz_0[i] = 3.0 * g_0_yyzzz_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxzzz_0[i] * pb_y + g_0_yyyzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxyyyy_0[i] = 2.0 * g_0_yyyyz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxxyyyy_0[i] * pb_z + g_0_yyyyzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyzzz_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyyyz_0[i] * pb_y + g_0_yyyzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyzzz_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyyzz_0[i] * pb_y + g_0_yyyzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyzzz_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyzzz_0[i] * pb_y + g_0_yyyzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyzzz_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxzzzz_0[i] * pb_y + g_0_yyyzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyyz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xxyyyyy_0[i] * pb_z + g_0_yyyyzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxyyyyz_0[i] = 3.0 * g_0_yyzzz_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyyyz_0[i] * pb_y + g_0_yyyzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyyyzz_0[i] = 3.0 * g_0_yyzzz_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyyzz_0[i] * pb_y + g_0_yyyzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyyzzz_0[i] = 3.0 * g_0_yyzzz_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyzzz_0[i] * pb_y + g_0_yyyzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyzzzz_0[i] = 3.0 * g_0_yyzzz_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyzzzz_0[i] * pb_y + g_0_yyyzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxzzzzz_0[i] = 3.0 * g_0_yyzzz_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxzzzzz_0[i] * pb_y + g_0_yyyzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyyyyy_0[i] = 2.0 * g_0_yyyyz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_xyyyyyy_0[i] * pb_z + g_0_yyyyzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xyyyyyz_0[i] = 3.0 * g_0_yyzzz_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyyyz_0[i] * pb_y + g_0_yyyzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyyyzz_0[i] = 3.0 * g_0_yyzzz_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyyzz_0[i] * pb_y + g_0_yyyzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyyzzz_0[i] = 3.0 * g_0_yyzzz_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyzzz_0[i] * pb_y + g_0_yyyzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyzzzz_0[i] = 3.0 * g_0_yyzzz_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyzzzz_0[i] * pb_y + g_0_yyyzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyzzzzz_0[i] = 3.0 * g_0_yyzzz_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyzzzzz_0[i] * pb_y + g_0_yyyzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xzzzzzz_0[i] = 3.0 * g_0_yyzzz_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_xzzzzzz_0[i] * pb_y + g_0_yyyzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyyyyy_0[i] = 2.0 * g_0_yyyyz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyyzz_0_yyyyyyy_0[i] * pb_z + g_0_yyyyzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_yyyyyyz_0[i] = 3.0 * g_0_yyzzz_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyyyyz_0[i] * pb_y + g_0_yyyzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyyyzz_0[i] = 3.0 * g_0_yyzzz_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyyyzz_0[i] * pb_y + g_0_yyyzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyyzzz_0[i] = 3.0 * g_0_yyzzz_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyyzzz_0[i] * pb_y + g_0_yyyzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyzzzz_0[i] = 3.0 * g_0_yyzzz_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyzzzz_0[i] * pb_y + g_0_yyyzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyzzzzz_0[i] = 3.0 * g_0_yyzzz_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyzzzzz_0[i] * pb_y + g_0_yyyzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yzzzzzz_0[i] = 3.0 * g_0_yyzzz_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yzzzzzz_0[i] * pb_y + g_0_yyyzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_zzzzzzz_0[i] = 3.0 * g_0_yyzzz_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyzzz_0_zzzzzzz_0[i] * pb_y + g_0_yyyzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1152-1188 components of targeted buffer : SKSK

    auto g_0_yyyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 1152);

    auto g_0_yyyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1153);

    auto g_0_yyyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 1154);

    auto g_0_yyyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 1155);

    auto g_0_yyyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 1156);

    auto g_0_yyyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 1157);

    auto g_0_yyyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 1158);

    auto g_0_yyyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 1159);

    auto g_0_yyyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 1160);

    auto g_0_yyyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 1161);

    auto g_0_yyyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1162);

    auto g_0_yyyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1163);

    auto g_0_yyyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1164);

    auto g_0_yyyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1165);

    auto g_0_yyyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1166);

    auto g_0_yyyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1167);

    auto g_0_yyyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1168);

    auto g_0_yyyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1169);

    auto g_0_yyyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1170);

    auto g_0_yyyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1171);

    auto g_0_yyyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1172);

    auto g_0_yyyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1173);

    auto g_0_yyyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1174);

    auto g_0_yyyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1175);

    auto g_0_yyyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1176);

    auto g_0_yyyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1177);

    auto g_0_yyyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1178);

    auto g_0_yyyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1179);

    auto g_0_yyyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1180);

    auto g_0_yyyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1181);

    auto g_0_yyyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1182);

    auto g_0_yyyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1183);

    auto g_0_yyyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1184);

    auto g_0_yyyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1185);

    auto g_0_yyyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1186);

    auto g_0_yyyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1187);

    #pragma omp simd aligned(g_0_yyyzz_0_xxxxxxy_0, g_0_yyyzz_0_xxxxxxy_1, g_0_yyyzz_0_xxxxxyy_0, g_0_yyyzz_0_xxxxxyy_1, g_0_yyyzz_0_xxxxyyy_0, g_0_yyyzz_0_xxxxyyy_1, g_0_yyyzz_0_xxxyyyy_0, g_0_yyyzz_0_xxxyyyy_1, g_0_yyyzz_0_xxyyyyy_0, g_0_yyyzz_0_xxyyyyy_1, g_0_yyyzz_0_xyyyyyy_0, g_0_yyyzz_0_xyyyyyy_1, g_0_yyyzz_0_yyyyyyy_0, g_0_yyyzz_0_yyyyyyy_1, g_0_yyyzzz_0_xxxxxxy_0, g_0_yyyzzz_0_xxxxxxy_1, g_0_yyyzzz_0_xxxxxyy_0, g_0_yyyzzz_0_xxxxxyy_1, g_0_yyyzzz_0_xxxxyyy_0, g_0_yyyzzz_0_xxxxyyy_1, g_0_yyyzzz_0_xxxyyyy_0, g_0_yyyzzz_0_xxxyyyy_1, g_0_yyyzzz_0_xxyyyyy_0, g_0_yyyzzz_0_xxyyyyy_1, g_0_yyyzzz_0_xyyyyyy_0, g_0_yyyzzz_0_xyyyyyy_1, g_0_yyyzzz_0_yyyyyyy_0, g_0_yyyzzz_0_yyyyyyy_1, g_0_yyyzzzz_0_xxxxxxx_0, g_0_yyyzzzz_0_xxxxxxy_0, g_0_yyyzzzz_0_xxxxxxz_0, g_0_yyyzzzz_0_xxxxxyy_0, g_0_yyyzzzz_0_xxxxxyz_0, g_0_yyyzzzz_0_xxxxxzz_0, g_0_yyyzzzz_0_xxxxyyy_0, g_0_yyyzzzz_0_xxxxyyz_0, g_0_yyyzzzz_0_xxxxyzz_0, g_0_yyyzzzz_0_xxxxzzz_0, g_0_yyyzzzz_0_xxxyyyy_0, g_0_yyyzzzz_0_xxxyyyz_0, g_0_yyyzzzz_0_xxxyyzz_0, g_0_yyyzzzz_0_xxxyzzz_0, g_0_yyyzzzz_0_xxxzzzz_0, g_0_yyyzzzz_0_xxyyyyy_0, g_0_yyyzzzz_0_xxyyyyz_0, g_0_yyyzzzz_0_xxyyyzz_0, g_0_yyyzzzz_0_xxyyzzz_0, g_0_yyyzzzz_0_xxyzzzz_0, g_0_yyyzzzz_0_xxzzzzz_0, g_0_yyyzzzz_0_xyyyyyy_0, g_0_yyyzzzz_0_xyyyyyz_0, g_0_yyyzzzz_0_xyyyyzz_0, g_0_yyyzzzz_0_xyyyzzz_0, g_0_yyyzzzz_0_xyyzzzz_0, g_0_yyyzzzz_0_xyzzzzz_0, g_0_yyyzzzz_0_xzzzzzz_0, g_0_yyyzzzz_0_yyyyyyy_0, g_0_yyyzzzz_0_yyyyyyz_0, g_0_yyyzzzz_0_yyyyyzz_0, g_0_yyyzzzz_0_yyyyzzz_0, g_0_yyyzzzz_0_yyyzzzz_0, g_0_yyyzzzz_0_yyzzzzz_0, g_0_yyyzzzz_0_yzzzzzz_0, g_0_yyyzzzz_0_zzzzzzz_0, g_0_yyzzzz_0_xxxxxxx_0, g_0_yyzzzz_0_xxxxxxx_1, g_0_yyzzzz_0_xxxxxxz_0, g_0_yyzzzz_0_xxxxxxz_1, g_0_yyzzzz_0_xxxxxyz_0, g_0_yyzzzz_0_xxxxxyz_1, g_0_yyzzzz_0_xxxxxz_1, g_0_yyzzzz_0_xxxxxzz_0, g_0_yyzzzz_0_xxxxxzz_1, g_0_yyzzzz_0_xxxxyyz_0, g_0_yyzzzz_0_xxxxyyz_1, g_0_yyzzzz_0_xxxxyz_1, g_0_yyzzzz_0_xxxxyzz_0, g_0_yyzzzz_0_xxxxyzz_1, g_0_yyzzzz_0_xxxxzz_1, g_0_yyzzzz_0_xxxxzzz_0, g_0_yyzzzz_0_xxxxzzz_1, g_0_yyzzzz_0_xxxyyyz_0, g_0_yyzzzz_0_xxxyyyz_1, g_0_yyzzzz_0_xxxyyz_1, g_0_yyzzzz_0_xxxyyzz_0, g_0_yyzzzz_0_xxxyyzz_1, g_0_yyzzzz_0_xxxyzz_1, g_0_yyzzzz_0_xxxyzzz_0, g_0_yyzzzz_0_xxxyzzz_1, g_0_yyzzzz_0_xxxzzz_1, g_0_yyzzzz_0_xxxzzzz_0, g_0_yyzzzz_0_xxxzzzz_1, g_0_yyzzzz_0_xxyyyyz_0, g_0_yyzzzz_0_xxyyyyz_1, g_0_yyzzzz_0_xxyyyz_1, g_0_yyzzzz_0_xxyyyzz_0, g_0_yyzzzz_0_xxyyyzz_1, g_0_yyzzzz_0_xxyyzz_1, g_0_yyzzzz_0_xxyyzzz_0, g_0_yyzzzz_0_xxyyzzz_1, g_0_yyzzzz_0_xxyzzz_1, g_0_yyzzzz_0_xxyzzzz_0, g_0_yyzzzz_0_xxyzzzz_1, g_0_yyzzzz_0_xxzzzz_1, g_0_yyzzzz_0_xxzzzzz_0, g_0_yyzzzz_0_xxzzzzz_1, g_0_yyzzzz_0_xyyyyyz_0, g_0_yyzzzz_0_xyyyyyz_1, g_0_yyzzzz_0_xyyyyz_1, g_0_yyzzzz_0_xyyyyzz_0, g_0_yyzzzz_0_xyyyyzz_1, g_0_yyzzzz_0_xyyyzz_1, g_0_yyzzzz_0_xyyyzzz_0, g_0_yyzzzz_0_xyyyzzz_1, g_0_yyzzzz_0_xyyzzz_1, g_0_yyzzzz_0_xyyzzzz_0, g_0_yyzzzz_0_xyyzzzz_1, g_0_yyzzzz_0_xyzzzz_1, g_0_yyzzzz_0_xyzzzzz_0, g_0_yyzzzz_0_xyzzzzz_1, g_0_yyzzzz_0_xzzzzz_1, g_0_yyzzzz_0_xzzzzzz_0, g_0_yyzzzz_0_xzzzzzz_1, g_0_yyzzzz_0_yyyyyyz_0, g_0_yyzzzz_0_yyyyyyz_1, g_0_yyzzzz_0_yyyyyz_1, g_0_yyzzzz_0_yyyyyzz_0, g_0_yyzzzz_0_yyyyyzz_1, g_0_yyzzzz_0_yyyyzz_1, g_0_yyzzzz_0_yyyyzzz_0, g_0_yyzzzz_0_yyyyzzz_1, g_0_yyzzzz_0_yyyzzz_1, g_0_yyzzzz_0_yyyzzzz_0, g_0_yyzzzz_0_yyyzzzz_1, g_0_yyzzzz_0_yyzzzz_1, g_0_yyzzzz_0_yyzzzzz_0, g_0_yyzzzz_0_yyzzzzz_1, g_0_yyzzzz_0_yzzzzz_1, g_0_yyzzzz_0_yzzzzzz_0, g_0_yyzzzz_0_yzzzzzz_1, g_0_yyzzzz_0_zzzzzz_1, g_0_yyzzzz_0_zzzzzzz_0, g_0_yyzzzz_0_zzzzzzz_1, g_0_yzzzz_0_xxxxxxx_0, g_0_yzzzz_0_xxxxxxx_1, g_0_yzzzz_0_xxxxxxz_0, g_0_yzzzz_0_xxxxxxz_1, g_0_yzzzz_0_xxxxxyz_0, g_0_yzzzz_0_xxxxxyz_1, g_0_yzzzz_0_xxxxxzz_0, g_0_yzzzz_0_xxxxxzz_1, g_0_yzzzz_0_xxxxyyz_0, g_0_yzzzz_0_xxxxyyz_1, g_0_yzzzz_0_xxxxyzz_0, g_0_yzzzz_0_xxxxyzz_1, g_0_yzzzz_0_xxxxzzz_0, g_0_yzzzz_0_xxxxzzz_1, g_0_yzzzz_0_xxxyyyz_0, g_0_yzzzz_0_xxxyyyz_1, g_0_yzzzz_0_xxxyyzz_0, g_0_yzzzz_0_xxxyyzz_1, g_0_yzzzz_0_xxxyzzz_0, g_0_yzzzz_0_xxxyzzz_1, g_0_yzzzz_0_xxxzzzz_0, g_0_yzzzz_0_xxxzzzz_1, g_0_yzzzz_0_xxyyyyz_0, g_0_yzzzz_0_xxyyyyz_1, g_0_yzzzz_0_xxyyyzz_0, g_0_yzzzz_0_xxyyyzz_1, g_0_yzzzz_0_xxyyzzz_0, g_0_yzzzz_0_xxyyzzz_1, g_0_yzzzz_0_xxyzzzz_0, g_0_yzzzz_0_xxyzzzz_1, g_0_yzzzz_0_xxzzzzz_0, g_0_yzzzz_0_xxzzzzz_1, g_0_yzzzz_0_xyyyyyz_0, g_0_yzzzz_0_xyyyyyz_1, g_0_yzzzz_0_xyyyyzz_0, g_0_yzzzz_0_xyyyyzz_1, g_0_yzzzz_0_xyyyzzz_0, g_0_yzzzz_0_xyyyzzz_1, g_0_yzzzz_0_xyyzzzz_0, g_0_yzzzz_0_xyyzzzz_1, g_0_yzzzz_0_xyzzzzz_0, g_0_yzzzz_0_xyzzzzz_1, g_0_yzzzz_0_xzzzzzz_0, g_0_yzzzz_0_xzzzzzz_1, g_0_yzzzz_0_yyyyyyz_0, g_0_yzzzz_0_yyyyyyz_1, g_0_yzzzz_0_yyyyyzz_0, g_0_yzzzz_0_yyyyyzz_1, g_0_yzzzz_0_yyyyzzz_0, g_0_yzzzz_0_yyyyzzz_1, g_0_yzzzz_0_yyyzzzz_0, g_0_yzzzz_0_yyyzzzz_1, g_0_yzzzz_0_yyzzzzz_0, g_0_yzzzz_0_yyzzzzz_1, g_0_yzzzz_0_yzzzzzz_0, g_0_yzzzz_0_yzzzzzz_1, g_0_yzzzz_0_zzzzzzz_0, g_0_yzzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzz_0_xxxxxxx_0[i] = 2.0 * g_0_yzzzz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxxxx_0[i] * pb_y + g_0_yyzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxxxy_0[i] = 3.0 * g_0_yyyzz_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxxxy_0[i] * pb_z + g_0_yyyzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxxxxz_0[i] = 2.0 * g_0_yzzzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxxxz_0[i] * pb_y + g_0_yyzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxxyy_0[i] = 3.0 * g_0_yyyzz_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxxyy_0[i] * pb_z + g_0_yyyzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxxxyz_0[i] = 2.0 * g_0_yzzzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxyz_0[i] * pb_y + g_0_yyzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxxzz_0[i] = 2.0 * g_0_yzzzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxxzz_0[i] * pb_y + g_0_yyzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_yyyzz_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxxyyy_0[i] * pb_z + g_0_yyyzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_yzzzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxyyz_0[i] * pb_y + g_0_yyzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxyzz_0[i] = 2.0 * g_0_yzzzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxyzz_0[i] * pb_y + g_0_yyzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxzzz_0[i] = 2.0 * g_0_yzzzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxzzz_0[i] * pb_y + g_0_yyzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyyzz_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxxyyyy_0[i] * pb_z + g_0_yyyzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxyyyz_0[i] = 2.0 * g_0_yzzzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyyyz_0[i] * pb_y + g_0_yyzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_yzzzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyyzz_0[i] * pb_y + g_0_yyzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxyzzz_0[i] = 2.0 * g_0_yzzzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyzzz_0[i] * pb_y + g_0_yyzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxzzzz_0[i] = 2.0 * g_0_yzzzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxzzzz_0[i] * pb_y + g_0_yyzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyyyyy_0[i] = 3.0 * g_0_yyyzz_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xxyyyyy_0[i] * pb_z + g_0_yyyzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzzzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyyyz_0[i] * pb_y + g_0_yyzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzzzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyyzz_0[i] * pb_y + g_0_yyzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzzzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyzzz_0[i] * pb_y + g_0_yyzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzzzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyzzzz_0[i] * pb_y + g_0_yyzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yzzzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxzzzzz_0[i] * pb_y + g_0_yyzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyyyyy_0[i] = 3.0 * g_0_yyyzz_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_xyyyyyy_0[i] * pb_z + g_0_yyyzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xyyyyyz_0[i] = 2.0 * g_0_yzzzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyyyz_0[i] * pb_y + g_0_yyzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyyyzz_0[i] = 2.0 * g_0_yzzzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyyzz_0[i] * pb_y + g_0_yyzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyyzzz_0[i] = 2.0 * g_0_yzzzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyzzz_0[i] * pb_y + g_0_yyzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_yzzzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyzzzz_0[i] * pb_y + g_0_yyzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyzzzzz_0[i] = 2.0 * g_0_yzzzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyzzzzz_0[i] * pb_y + g_0_yyzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xzzzzzz_0[i] = 2.0 * g_0_yzzzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_xzzzzzz_0[i] * pb_y + g_0_yyzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyyyyy_0[i] = 3.0 * g_0_yyyzz_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyzzz_0_yyyyyyy_0[i] * pb_z + g_0_yyyzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_yyyyyyz_0[i] = 2.0 * g_0_yzzzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyyyyz_0[i] * pb_y + g_0_yyzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyyyzz_0[i] = 2.0 * g_0_yzzzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyyyzz_0[i] * pb_y + g_0_yyzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyyzzz_0[i] = 2.0 * g_0_yzzzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyyzzz_0[i] * pb_y + g_0_yyzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyzzzz_0[i] = 2.0 * g_0_yzzzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyzzzz_0[i] * pb_y + g_0_yyzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_yzzzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyzzzzz_0[i] * pb_y + g_0_yyzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yzzzzzz_0[i] = 2.0 * g_0_yzzzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yzzzzzz_0[i] * pb_y + g_0_yyzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_zzzzzzz_0[i] = 2.0 * g_0_yzzzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyzzzz_0_zzzzzzz_0[i] * pb_y + g_0_yyzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1188-1224 components of targeted buffer : SKSK

    auto g_0_yyzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 1188);

    auto g_0_yyzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1189);

    auto g_0_yyzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 1190);

    auto g_0_yyzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 1191);

    auto g_0_yyzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 1192);

    auto g_0_yyzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 1193);

    auto g_0_yyzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 1194);

    auto g_0_yyzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 1195);

    auto g_0_yyzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 1196);

    auto g_0_yyzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 1197);

    auto g_0_yyzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1198);

    auto g_0_yyzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1199);

    auto g_0_yyzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1200);

    auto g_0_yyzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1201);

    auto g_0_yyzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1202);

    auto g_0_yyzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1203);

    auto g_0_yyzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1204);

    auto g_0_yyzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1205);

    auto g_0_yyzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1206);

    auto g_0_yyzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1207);

    auto g_0_yyzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1208);

    auto g_0_yyzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1209);

    auto g_0_yyzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1210);

    auto g_0_yyzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1211);

    auto g_0_yyzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1212);

    auto g_0_yyzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1213);

    auto g_0_yyzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1214);

    auto g_0_yyzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1215);

    auto g_0_yyzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1216);

    auto g_0_yyzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1217);

    auto g_0_yyzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1218);

    auto g_0_yyzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1219);

    auto g_0_yyzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1220);

    auto g_0_yyzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1221);

    auto g_0_yyzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1222);

    auto g_0_yyzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1223);

    #pragma omp simd aligned(g_0_yyzzz_0_xxxxxxy_0, g_0_yyzzz_0_xxxxxxy_1, g_0_yyzzz_0_xxxxxyy_0, g_0_yyzzz_0_xxxxxyy_1, g_0_yyzzz_0_xxxxyyy_0, g_0_yyzzz_0_xxxxyyy_1, g_0_yyzzz_0_xxxyyyy_0, g_0_yyzzz_0_xxxyyyy_1, g_0_yyzzz_0_xxyyyyy_0, g_0_yyzzz_0_xxyyyyy_1, g_0_yyzzz_0_xyyyyyy_0, g_0_yyzzz_0_xyyyyyy_1, g_0_yyzzz_0_yyyyyyy_0, g_0_yyzzz_0_yyyyyyy_1, g_0_yyzzzz_0_xxxxxxy_0, g_0_yyzzzz_0_xxxxxxy_1, g_0_yyzzzz_0_xxxxxyy_0, g_0_yyzzzz_0_xxxxxyy_1, g_0_yyzzzz_0_xxxxyyy_0, g_0_yyzzzz_0_xxxxyyy_1, g_0_yyzzzz_0_xxxyyyy_0, g_0_yyzzzz_0_xxxyyyy_1, g_0_yyzzzz_0_xxyyyyy_0, g_0_yyzzzz_0_xxyyyyy_1, g_0_yyzzzz_0_xyyyyyy_0, g_0_yyzzzz_0_xyyyyyy_1, g_0_yyzzzz_0_yyyyyyy_0, g_0_yyzzzz_0_yyyyyyy_1, g_0_yyzzzzz_0_xxxxxxx_0, g_0_yyzzzzz_0_xxxxxxy_0, g_0_yyzzzzz_0_xxxxxxz_0, g_0_yyzzzzz_0_xxxxxyy_0, g_0_yyzzzzz_0_xxxxxyz_0, g_0_yyzzzzz_0_xxxxxzz_0, g_0_yyzzzzz_0_xxxxyyy_0, g_0_yyzzzzz_0_xxxxyyz_0, g_0_yyzzzzz_0_xxxxyzz_0, g_0_yyzzzzz_0_xxxxzzz_0, g_0_yyzzzzz_0_xxxyyyy_0, g_0_yyzzzzz_0_xxxyyyz_0, g_0_yyzzzzz_0_xxxyyzz_0, g_0_yyzzzzz_0_xxxyzzz_0, g_0_yyzzzzz_0_xxxzzzz_0, g_0_yyzzzzz_0_xxyyyyy_0, g_0_yyzzzzz_0_xxyyyyz_0, g_0_yyzzzzz_0_xxyyyzz_0, g_0_yyzzzzz_0_xxyyzzz_0, g_0_yyzzzzz_0_xxyzzzz_0, g_0_yyzzzzz_0_xxzzzzz_0, g_0_yyzzzzz_0_xyyyyyy_0, g_0_yyzzzzz_0_xyyyyyz_0, g_0_yyzzzzz_0_xyyyyzz_0, g_0_yyzzzzz_0_xyyyzzz_0, g_0_yyzzzzz_0_xyyzzzz_0, g_0_yyzzzzz_0_xyzzzzz_0, g_0_yyzzzzz_0_xzzzzzz_0, g_0_yyzzzzz_0_yyyyyyy_0, g_0_yyzzzzz_0_yyyyyyz_0, g_0_yyzzzzz_0_yyyyyzz_0, g_0_yyzzzzz_0_yyyyzzz_0, g_0_yyzzzzz_0_yyyzzzz_0, g_0_yyzzzzz_0_yyzzzzz_0, g_0_yyzzzzz_0_yzzzzzz_0, g_0_yyzzzzz_0_zzzzzzz_0, g_0_yzzzzz_0_xxxxxxx_0, g_0_yzzzzz_0_xxxxxxx_1, g_0_yzzzzz_0_xxxxxxz_0, g_0_yzzzzz_0_xxxxxxz_1, g_0_yzzzzz_0_xxxxxyz_0, g_0_yzzzzz_0_xxxxxyz_1, g_0_yzzzzz_0_xxxxxz_1, g_0_yzzzzz_0_xxxxxzz_0, g_0_yzzzzz_0_xxxxxzz_1, g_0_yzzzzz_0_xxxxyyz_0, g_0_yzzzzz_0_xxxxyyz_1, g_0_yzzzzz_0_xxxxyz_1, g_0_yzzzzz_0_xxxxyzz_0, g_0_yzzzzz_0_xxxxyzz_1, g_0_yzzzzz_0_xxxxzz_1, g_0_yzzzzz_0_xxxxzzz_0, g_0_yzzzzz_0_xxxxzzz_1, g_0_yzzzzz_0_xxxyyyz_0, g_0_yzzzzz_0_xxxyyyz_1, g_0_yzzzzz_0_xxxyyz_1, g_0_yzzzzz_0_xxxyyzz_0, g_0_yzzzzz_0_xxxyyzz_1, g_0_yzzzzz_0_xxxyzz_1, g_0_yzzzzz_0_xxxyzzz_0, g_0_yzzzzz_0_xxxyzzz_1, g_0_yzzzzz_0_xxxzzz_1, g_0_yzzzzz_0_xxxzzzz_0, g_0_yzzzzz_0_xxxzzzz_1, g_0_yzzzzz_0_xxyyyyz_0, g_0_yzzzzz_0_xxyyyyz_1, g_0_yzzzzz_0_xxyyyz_1, g_0_yzzzzz_0_xxyyyzz_0, g_0_yzzzzz_0_xxyyyzz_1, g_0_yzzzzz_0_xxyyzz_1, g_0_yzzzzz_0_xxyyzzz_0, g_0_yzzzzz_0_xxyyzzz_1, g_0_yzzzzz_0_xxyzzz_1, g_0_yzzzzz_0_xxyzzzz_0, g_0_yzzzzz_0_xxyzzzz_1, g_0_yzzzzz_0_xxzzzz_1, g_0_yzzzzz_0_xxzzzzz_0, g_0_yzzzzz_0_xxzzzzz_1, g_0_yzzzzz_0_xyyyyyz_0, g_0_yzzzzz_0_xyyyyyz_1, g_0_yzzzzz_0_xyyyyz_1, g_0_yzzzzz_0_xyyyyzz_0, g_0_yzzzzz_0_xyyyyzz_1, g_0_yzzzzz_0_xyyyzz_1, g_0_yzzzzz_0_xyyyzzz_0, g_0_yzzzzz_0_xyyyzzz_1, g_0_yzzzzz_0_xyyzzz_1, g_0_yzzzzz_0_xyyzzzz_0, g_0_yzzzzz_0_xyyzzzz_1, g_0_yzzzzz_0_xyzzzz_1, g_0_yzzzzz_0_xyzzzzz_0, g_0_yzzzzz_0_xyzzzzz_1, g_0_yzzzzz_0_xzzzzz_1, g_0_yzzzzz_0_xzzzzzz_0, g_0_yzzzzz_0_xzzzzzz_1, g_0_yzzzzz_0_yyyyyyz_0, g_0_yzzzzz_0_yyyyyyz_1, g_0_yzzzzz_0_yyyyyz_1, g_0_yzzzzz_0_yyyyyzz_0, g_0_yzzzzz_0_yyyyyzz_1, g_0_yzzzzz_0_yyyyzz_1, g_0_yzzzzz_0_yyyyzzz_0, g_0_yzzzzz_0_yyyyzzz_1, g_0_yzzzzz_0_yyyzzz_1, g_0_yzzzzz_0_yyyzzzz_0, g_0_yzzzzz_0_yyyzzzz_1, g_0_yzzzzz_0_yyzzzz_1, g_0_yzzzzz_0_yyzzzzz_0, g_0_yzzzzz_0_yyzzzzz_1, g_0_yzzzzz_0_yzzzzz_1, g_0_yzzzzz_0_yzzzzzz_0, g_0_yzzzzz_0_yzzzzzz_1, g_0_yzzzzz_0_zzzzzz_1, g_0_yzzzzz_0_zzzzzzz_0, g_0_yzzzzz_0_zzzzzzz_1, g_0_zzzzz_0_xxxxxxx_0, g_0_zzzzz_0_xxxxxxx_1, g_0_zzzzz_0_xxxxxxz_0, g_0_zzzzz_0_xxxxxxz_1, g_0_zzzzz_0_xxxxxyz_0, g_0_zzzzz_0_xxxxxyz_1, g_0_zzzzz_0_xxxxxzz_0, g_0_zzzzz_0_xxxxxzz_1, g_0_zzzzz_0_xxxxyyz_0, g_0_zzzzz_0_xxxxyyz_1, g_0_zzzzz_0_xxxxyzz_0, g_0_zzzzz_0_xxxxyzz_1, g_0_zzzzz_0_xxxxzzz_0, g_0_zzzzz_0_xxxxzzz_1, g_0_zzzzz_0_xxxyyyz_0, g_0_zzzzz_0_xxxyyyz_1, g_0_zzzzz_0_xxxyyzz_0, g_0_zzzzz_0_xxxyyzz_1, g_0_zzzzz_0_xxxyzzz_0, g_0_zzzzz_0_xxxyzzz_1, g_0_zzzzz_0_xxxzzzz_0, g_0_zzzzz_0_xxxzzzz_1, g_0_zzzzz_0_xxyyyyz_0, g_0_zzzzz_0_xxyyyyz_1, g_0_zzzzz_0_xxyyyzz_0, g_0_zzzzz_0_xxyyyzz_1, g_0_zzzzz_0_xxyyzzz_0, g_0_zzzzz_0_xxyyzzz_1, g_0_zzzzz_0_xxyzzzz_0, g_0_zzzzz_0_xxyzzzz_1, g_0_zzzzz_0_xxzzzzz_0, g_0_zzzzz_0_xxzzzzz_1, g_0_zzzzz_0_xyyyyyz_0, g_0_zzzzz_0_xyyyyyz_1, g_0_zzzzz_0_xyyyyzz_0, g_0_zzzzz_0_xyyyyzz_1, g_0_zzzzz_0_xyyyzzz_0, g_0_zzzzz_0_xyyyzzz_1, g_0_zzzzz_0_xyyzzzz_0, g_0_zzzzz_0_xyyzzzz_1, g_0_zzzzz_0_xyzzzzz_0, g_0_zzzzz_0_xyzzzzz_1, g_0_zzzzz_0_xzzzzzz_0, g_0_zzzzz_0_xzzzzzz_1, g_0_zzzzz_0_yyyyyyz_0, g_0_zzzzz_0_yyyyyyz_1, g_0_zzzzz_0_yyyyyzz_0, g_0_zzzzz_0_yyyyyzz_1, g_0_zzzzz_0_yyyyzzz_0, g_0_zzzzz_0_yyyyzzz_1, g_0_zzzzz_0_yyyzzzz_0, g_0_zzzzz_0_yyyzzzz_1, g_0_zzzzz_0_yyzzzzz_0, g_0_zzzzz_0_yyzzzzz_1, g_0_zzzzz_0_yzzzzzz_0, g_0_zzzzz_0_yzzzzzz_1, g_0_zzzzz_0_zzzzzzz_0, g_0_zzzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzz_0_xxxxxxx_0[i] = g_0_zzzzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxxxx_0[i] * pb_y + g_0_yzzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxxxy_0[i] = 4.0 * g_0_yyzzz_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxxxy_0[i] * pb_z + g_0_yyzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxxxxz_0[i] = g_0_zzzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxxxz_0[i] * pb_y + g_0_yzzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxxyy_0[i] = 4.0 * g_0_yyzzz_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxxyy_0[i] * pb_z + g_0_yyzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxxxyz_0[i] = g_0_zzzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxxyz_0[i] * pb_y + g_0_yzzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxxzz_0[i] = g_0_zzzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxxzz_0[i] * pb_y + g_0_yzzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyzzz_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxxyyy_0[i] * pb_z + g_0_yyzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxxyyz_0[i] = g_0_zzzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxyyz_0[i] * pb_y + g_0_yzzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxyzz_0[i] = g_0_zzzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxyzz_0[i] * pb_y + g_0_yzzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxzzz_0[i] = g_0_zzzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxzzz_0[i] * pb_y + g_0_yzzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_yyzzz_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxxyyyy_0[i] * pb_z + g_0_yyzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxyyyz_0[i] = g_0_zzzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyyyz_0[i] * pb_y + g_0_yzzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxyyzz_0[i] = g_0_zzzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyyzz_0[i] * pb_y + g_0_yzzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxyzzz_0[i] = g_0_zzzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyzzz_0[i] * pb_y + g_0_yzzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxzzzz_0[i] = g_0_zzzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxzzzz_0[i] * pb_y + g_0_yzzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyyyyy_0[i] = 4.0 * g_0_yyzzz_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xxyyyyy_0[i] * pb_z + g_0_yyzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxyyyyz_0[i] = g_0_zzzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyyyz_0[i] * pb_y + g_0_yzzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyyyzz_0[i] = g_0_zzzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyyzz_0[i] * pb_y + g_0_yzzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyyzzz_0[i] = g_0_zzzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyzzz_0[i] * pb_y + g_0_yzzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyzzzz_0[i] = g_0_zzzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyzzzz_0[i] * pb_y + g_0_yzzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxzzzzz_0[i] = g_0_zzzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxzzzzz_0[i] * pb_y + g_0_yzzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyyyyy_0[i] = 4.0 * g_0_yyzzz_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_xyyyyyy_0[i] * pb_z + g_0_yyzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xyyyyyz_0[i] = g_0_zzzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyyyz_0[i] * pb_y + g_0_yzzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyyyzz_0[i] = g_0_zzzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyyzz_0[i] * pb_y + g_0_yzzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyyzzz_0[i] = g_0_zzzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyzzz_0[i] * pb_y + g_0_yzzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyzzzz_0[i] = g_0_zzzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyzzzz_0[i] * pb_y + g_0_yzzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyzzzzz_0[i] = g_0_zzzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyzzzzz_0[i] * pb_y + g_0_yzzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xzzzzzz_0[i] = g_0_zzzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzzzzzz_0[i] * pb_y + g_0_yzzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyyyyy_0[i] = 4.0 * g_0_yyzzz_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyzzzz_0_yyyyyyy_0[i] * pb_z + g_0_yyzzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_yyyyyyz_0[i] = g_0_zzzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyyyyz_0[i] * pb_y + g_0_yzzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyyyzz_0[i] = g_0_zzzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyyyzz_0[i] * pb_y + g_0_yzzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyyzzz_0[i] = g_0_zzzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyyzzz_0[i] * pb_y + g_0_yzzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyzzzz_0[i] = g_0_zzzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyzzzz_0[i] * pb_y + g_0_yzzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyzzzzz_0[i] = g_0_zzzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyzzzzz_0[i] * pb_y + g_0_yzzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yzzzzzz_0[i] = g_0_zzzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yzzzzzz_0[i] * pb_y + g_0_yzzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_zzzzzzz_0[i] = g_0_zzzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzzzzzz_0[i] * pb_y + g_0_yzzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1224-1260 components of targeted buffer : SKSK

    auto g_0_yzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 1224);

    auto g_0_yzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1225);

    auto g_0_yzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 1226);

    auto g_0_yzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 1227);

    auto g_0_yzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 1228);

    auto g_0_yzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 1229);

    auto g_0_yzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 1230);

    auto g_0_yzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 1231);

    auto g_0_yzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 1232);

    auto g_0_yzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 1233);

    auto g_0_yzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1234);

    auto g_0_yzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1235);

    auto g_0_yzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1236);

    auto g_0_yzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1237);

    auto g_0_yzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1238);

    auto g_0_yzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1239);

    auto g_0_yzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1240);

    auto g_0_yzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1241);

    auto g_0_yzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1242);

    auto g_0_yzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1243);

    auto g_0_yzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1244);

    auto g_0_yzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1245);

    auto g_0_yzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1246);

    auto g_0_yzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1247);

    auto g_0_yzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1248);

    auto g_0_yzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1249);

    auto g_0_yzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1250);

    auto g_0_yzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1251);

    auto g_0_yzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1252);

    auto g_0_yzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1253);

    auto g_0_yzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1254);

    auto g_0_yzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1255);

    auto g_0_yzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1256);

    auto g_0_yzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1257);

    auto g_0_yzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1258);

    auto g_0_yzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1259);

    #pragma omp simd aligned(g_0_yzzzzzz_0_xxxxxxx_0, g_0_yzzzzzz_0_xxxxxxy_0, g_0_yzzzzzz_0_xxxxxxz_0, g_0_yzzzzzz_0_xxxxxyy_0, g_0_yzzzzzz_0_xxxxxyz_0, g_0_yzzzzzz_0_xxxxxzz_0, g_0_yzzzzzz_0_xxxxyyy_0, g_0_yzzzzzz_0_xxxxyyz_0, g_0_yzzzzzz_0_xxxxyzz_0, g_0_yzzzzzz_0_xxxxzzz_0, g_0_yzzzzzz_0_xxxyyyy_0, g_0_yzzzzzz_0_xxxyyyz_0, g_0_yzzzzzz_0_xxxyyzz_0, g_0_yzzzzzz_0_xxxyzzz_0, g_0_yzzzzzz_0_xxxzzzz_0, g_0_yzzzzzz_0_xxyyyyy_0, g_0_yzzzzzz_0_xxyyyyz_0, g_0_yzzzzzz_0_xxyyyzz_0, g_0_yzzzzzz_0_xxyyzzz_0, g_0_yzzzzzz_0_xxyzzzz_0, g_0_yzzzzzz_0_xxzzzzz_0, g_0_yzzzzzz_0_xyyyyyy_0, g_0_yzzzzzz_0_xyyyyyz_0, g_0_yzzzzzz_0_xyyyyzz_0, g_0_yzzzzzz_0_xyyyzzz_0, g_0_yzzzzzz_0_xyyzzzz_0, g_0_yzzzzzz_0_xyzzzzz_0, g_0_yzzzzzz_0_xzzzzzz_0, g_0_yzzzzzz_0_yyyyyyy_0, g_0_yzzzzzz_0_yyyyyyz_0, g_0_yzzzzzz_0_yyyyyzz_0, g_0_yzzzzzz_0_yyyyzzz_0, g_0_yzzzzzz_0_yyyzzzz_0, g_0_yzzzzzz_0_yyzzzzz_0, g_0_yzzzzzz_0_yzzzzzz_0, g_0_yzzzzzz_0_zzzzzzz_0, g_0_zzzzzz_0_xxxxxx_1, g_0_zzzzzz_0_xxxxxxx_0, g_0_zzzzzz_0_xxxxxxx_1, g_0_zzzzzz_0_xxxxxxy_0, g_0_zzzzzz_0_xxxxxxy_1, g_0_zzzzzz_0_xxxxxxz_0, g_0_zzzzzz_0_xxxxxxz_1, g_0_zzzzzz_0_xxxxxy_1, g_0_zzzzzz_0_xxxxxyy_0, g_0_zzzzzz_0_xxxxxyy_1, g_0_zzzzzz_0_xxxxxyz_0, g_0_zzzzzz_0_xxxxxyz_1, g_0_zzzzzz_0_xxxxxz_1, g_0_zzzzzz_0_xxxxxzz_0, g_0_zzzzzz_0_xxxxxzz_1, g_0_zzzzzz_0_xxxxyy_1, g_0_zzzzzz_0_xxxxyyy_0, g_0_zzzzzz_0_xxxxyyy_1, g_0_zzzzzz_0_xxxxyyz_0, g_0_zzzzzz_0_xxxxyyz_1, g_0_zzzzzz_0_xxxxyz_1, g_0_zzzzzz_0_xxxxyzz_0, g_0_zzzzzz_0_xxxxyzz_1, g_0_zzzzzz_0_xxxxzz_1, g_0_zzzzzz_0_xxxxzzz_0, g_0_zzzzzz_0_xxxxzzz_1, g_0_zzzzzz_0_xxxyyy_1, g_0_zzzzzz_0_xxxyyyy_0, g_0_zzzzzz_0_xxxyyyy_1, g_0_zzzzzz_0_xxxyyyz_0, g_0_zzzzzz_0_xxxyyyz_1, g_0_zzzzzz_0_xxxyyz_1, g_0_zzzzzz_0_xxxyyzz_0, g_0_zzzzzz_0_xxxyyzz_1, g_0_zzzzzz_0_xxxyzz_1, g_0_zzzzzz_0_xxxyzzz_0, g_0_zzzzzz_0_xxxyzzz_1, g_0_zzzzzz_0_xxxzzz_1, g_0_zzzzzz_0_xxxzzzz_0, g_0_zzzzzz_0_xxxzzzz_1, g_0_zzzzzz_0_xxyyyy_1, g_0_zzzzzz_0_xxyyyyy_0, g_0_zzzzzz_0_xxyyyyy_1, g_0_zzzzzz_0_xxyyyyz_0, g_0_zzzzzz_0_xxyyyyz_1, g_0_zzzzzz_0_xxyyyz_1, g_0_zzzzzz_0_xxyyyzz_0, g_0_zzzzzz_0_xxyyyzz_1, g_0_zzzzzz_0_xxyyzz_1, g_0_zzzzzz_0_xxyyzzz_0, g_0_zzzzzz_0_xxyyzzz_1, g_0_zzzzzz_0_xxyzzz_1, g_0_zzzzzz_0_xxyzzzz_0, g_0_zzzzzz_0_xxyzzzz_1, g_0_zzzzzz_0_xxzzzz_1, g_0_zzzzzz_0_xxzzzzz_0, g_0_zzzzzz_0_xxzzzzz_1, g_0_zzzzzz_0_xyyyyy_1, g_0_zzzzzz_0_xyyyyyy_0, g_0_zzzzzz_0_xyyyyyy_1, g_0_zzzzzz_0_xyyyyyz_0, g_0_zzzzzz_0_xyyyyyz_1, g_0_zzzzzz_0_xyyyyz_1, g_0_zzzzzz_0_xyyyyzz_0, g_0_zzzzzz_0_xyyyyzz_1, g_0_zzzzzz_0_xyyyzz_1, g_0_zzzzzz_0_xyyyzzz_0, g_0_zzzzzz_0_xyyyzzz_1, g_0_zzzzzz_0_xyyzzz_1, g_0_zzzzzz_0_xyyzzzz_0, g_0_zzzzzz_0_xyyzzzz_1, g_0_zzzzzz_0_xyzzzz_1, g_0_zzzzzz_0_xyzzzzz_0, g_0_zzzzzz_0_xyzzzzz_1, g_0_zzzzzz_0_xzzzzz_1, g_0_zzzzzz_0_xzzzzzz_0, g_0_zzzzzz_0_xzzzzzz_1, g_0_zzzzzz_0_yyyyyy_1, g_0_zzzzzz_0_yyyyyyy_0, g_0_zzzzzz_0_yyyyyyy_1, g_0_zzzzzz_0_yyyyyyz_0, g_0_zzzzzz_0_yyyyyyz_1, g_0_zzzzzz_0_yyyyyz_1, g_0_zzzzzz_0_yyyyyzz_0, g_0_zzzzzz_0_yyyyyzz_1, g_0_zzzzzz_0_yyyyzz_1, g_0_zzzzzz_0_yyyyzzz_0, g_0_zzzzzz_0_yyyyzzz_1, g_0_zzzzzz_0_yyyzzz_1, g_0_zzzzzz_0_yyyzzzz_0, g_0_zzzzzz_0_yyyzzzz_1, g_0_zzzzzz_0_yyzzzz_1, g_0_zzzzzz_0_yyzzzzz_0, g_0_zzzzzz_0_yyzzzzz_1, g_0_zzzzzz_0_yzzzzz_1, g_0_zzzzzz_0_yzzzzzz_0, g_0_zzzzzz_0_yzzzzzz_1, g_0_zzzzzz_0_zzzzzz_1, g_0_zzzzzz_0_zzzzzzz_0, g_0_zzzzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzz_0_xxxxxxx_0[i] = g_0_zzzzzz_0_xxxxxxx_0[i] * pb_y + g_0_zzzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxxxy_0[i] = g_0_zzzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxxy_0[i] * pb_y + g_0_zzzzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxxxz_0[i] = g_0_zzzzzz_0_xxxxxxz_0[i] * pb_y + g_0_zzzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxxyy_0[i] = 2.0 * g_0_zzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxyy_0[i] * pb_y + g_0_zzzzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxxyz_0[i] = g_0_zzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxyz_0[i] * pb_y + g_0_zzzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxxzz_0[i] = g_0_zzzzzz_0_xxxxxzz_0[i] * pb_y + g_0_zzzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_zzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyyy_0[i] * pb_y + g_0_zzzzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_zzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyyz_0[i] * pb_y + g_0_zzzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxyzz_0[i] = g_0_zzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyzz_0[i] * pb_y + g_0_zzzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxzzz_0[i] = g_0_zzzzzz_0_xxxxzzz_0[i] * pb_y + g_0_zzzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_zzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyyy_0[i] * pb_y + g_0_zzzzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyyz_0[i] * pb_y + g_0_zzzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_zzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyzz_0[i] * pb_y + g_0_zzzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyzzz_0[i] = g_0_zzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyzzz_0[i] * pb_y + g_0_zzzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxzzzz_0[i] = g_0_zzzzzz_0_xxxzzzz_0[i] * pb_y + g_0_zzzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_zzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyyy_0[i] * pb_y + g_0_zzzzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyyyz_0[i] = 4.0 * g_0_zzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyyz_0[i] * pb_y + g_0_zzzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_zzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyzz_0[i] * pb_y + g_0_zzzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyzzz_0[i] * pb_y + g_0_zzzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyzzzz_0[i] = g_0_zzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzzzz_0[i] * pb_y + g_0_zzzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxzzzzz_0[i] = g_0_zzzzzz_0_xxzzzzz_0[i] * pb_y + g_0_zzzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyyyy_0[i] = 6.0 * g_0_zzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyyy_0[i] * pb_y + g_0_zzzzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyyyz_0[i] = 5.0 * g_0_zzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyyz_0[i] * pb_y + g_0_zzzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyyzz_0[i] = 4.0 * g_0_zzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyzz_0[i] * pb_y + g_0_zzzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_zzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyzzz_0[i] * pb_y + g_0_zzzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_zzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzzzz_0[i] * pb_y + g_0_zzzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyzzzzz_0[i] = g_0_zzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzzzz_0[i] * pb_y + g_0_zzzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xzzzzzz_0[i] = g_0_zzzzzz_0_xzzzzzz_0[i] * pb_y + g_0_zzzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyyyy_0[i] = 7.0 * g_0_zzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyyyy_0[i] * pb_y + g_0_zzzzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyyyz_0[i] = 6.0 * g_0_zzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyyyz_0[i] * pb_y + g_0_zzzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyyzz_0[i] = 5.0 * g_0_zzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyyzz_0[i] * pb_y + g_0_zzzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyzzz_0[i] = 4.0 * g_0_zzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyzzz_0[i] * pb_y + g_0_zzzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_zzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyzzzz_0[i] * pb_y + g_0_zzzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_zzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyzzzzz_0[i] * pb_y + g_0_zzzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yzzzzzz_0[i] = g_0_zzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzzzzzz_0[i] * pb_y + g_0_zzzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_zzzzzzz_0[i] = g_0_zzzzzz_0_zzzzzzz_0[i] * pb_y + g_0_zzzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1260-1296 components of targeted buffer : SKSK

    auto g_0_zzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 1260);

    auto g_0_zzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 1261);

    auto g_0_zzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 1262);

    auto g_0_zzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 1263);

    auto g_0_zzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 1264);

    auto g_0_zzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 1265);

    auto g_0_zzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 1266);

    auto g_0_zzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 1267);

    auto g_0_zzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 1268);

    auto g_0_zzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 1269);

    auto g_0_zzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1270);

    auto g_0_zzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1271);

    auto g_0_zzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1272);

    auto g_0_zzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1273);

    auto g_0_zzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1274);

    auto g_0_zzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1275);

    auto g_0_zzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1276);

    auto g_0_zzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1277);

    auto g_0_zzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1278);

    auto g_0_zzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1279);

    auto g_0_zzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1280);

    auto g_0_zzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1281);

    auto g_0_zzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1282);

    auto g_0_zzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1283);

    auto g_0_zzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1284);

    auto g_0_zzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1285);

    auto g_0_zzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1286);

    auto g_0_zzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1287);

    auto g_0_zzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 1288);

    auto g_0_zzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 1289);

    auto g_0_zzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 1290);

    auto g_0_zzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 1291);

    auto g_0_zzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1292);

    auto g_0_zzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1293);

    auto g_0_zzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1294);

    auto g_0_zzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 1295);

    #pragma omp simd aligned(g_0_zzzzz_0_xxxxxxx_0, g_0_zzzzz_0_xxxxxxx_1, g_0_zzzzz_0_xxxxxxy_0, g_0_zzzzz_0_xxxxxxy_1, g_0_zzzzz_0_xxxxxxz_0, g_0_zzzzz_0_xxxxxxz_1, g_0_zzzzz_0_xxxxxyy_0, g_0_zzzzz_0_xxxxxyy_1, g_0_zzzzz_0_xxxxxyz_0, g_0_zzzzz_0_xxxxxyz_1, g_0_zzzzz_0_xxxxxzz_0, g_0_zzzzz_0_xxxxxzz_1, g_0_zzzzz_0_xxxxyyy_0, g_0_zzzzz_0_xxxxyyy_1, g_0_zzzzz_0_xxxxyyz_0, g_0_zzzzz_0_xxxxyyz_1, g_0_zzzzz_0_xxxxyzz_0, g_0_zzzzz_0_xxxxyzz_1, g_0_zzzzz_0_xxxxzzz_0, g_0_zzzzz_0_xxxxzzz_1, g_0_zzzzz_0_xxxyyyy_0, g_0_zzzzz_0_xxxyyyy_1, g_0_zzzzz_0_xxxyyyz_0, g_0_zzzzz_0_xxxyyyz_1, g_0_zzzzz_0_xxxyyzz_0, g_0_zzzzz_0_xxxyyzz_1, g_0_zzzzz_0_xxxyzzz_0, g_0_zzzzz_0_xxxyzzz_1, g_0_zzzzz_0_xxxzzzz_0, g_0_zzzzz_0_xxxzzzz_1, g_0_zzzzz_0_xxyyyyy_0, g_0_zzzzz_0_xxyyyyy_1, g_0_zzzzz_0_xxyyyyz_0, g_0_zzzzz_0_xxyyyyz_1, g_0_zzzzz_0_xxyyyzz_0, g_0_zzzzz_0_xxyyyzz_1, g_0_zzzzz_0_xxyyzzz_0, g_0_zzzzz_0_xxyyzzz_1, g_0_zzzzz_0_xxyzzzz_0, g_0_zzzzz_0_xxyzzzz_1, g_0_zzzzz_0_xxzzzzz_0, g_0_zzzzz_0_xxzzzzz_1, g_0_zzzzz_0_xyyyyyy_0, g_0_zzzzz_0_xyyyyyy_1, g_0_zzzzz_0_xyyyyyz_0, g_0_zzzzz_0_xyyyyyz_1, g_0_zzzzz_0_xyyyyzz_0, g_0_zzzzz_0_xyyyyzz_1, g_0_zzzzz_0_xyyyzzz_0, g_0_zzzzz_0_xyyyzzz_1, g_0_zzzzz_0_xyyzzzz_0, g_0_zzzzz_0_xyyzzzz_1, g_0_zzzzz_0_xyzzzzz_0, g_0_zzzzz_0_xyzzzzz_1, g_0_zzzzz_0_xzzzzzz_0, g_0_zzzzz_0_xzzzzzz_1, g_0_zzzzz_0_yyyyyyy_0, g_0_zzzzz_0_yyyyyyy_1, g_0_zzzzz_0_yyyyyyz_0, g_0_zzzzz_0_yyyyyyz_1, g_0_zzzzz_0_yyyyyzz_0, g_0_zzzzz_0_yyyyyzz_1, g_0_zzzzz_0_yyyyzzz_0, g_0_zzzzz_0_yyyyzzz_1, g_0_zzzzz_0_yyyzzzz_0, g_0_zzzzz_0_yyyzzzz_1, g_0_zzzzz_0_yyzzzzz_0, g_0_zzzzz_0_yyzzzzz_1, g_0_zzzzz_0_yzzzzzz_0, g_0_zzzzz_0_yzzzzzz_1, g_0_zzzzz_0_zzzzzzz_0, g_0_zzzzz_0_zzzzzzz_1, g_0_zzzzzz_0_xxxxxx_1, g_0_zzzzzz_0_xxxxxxx_0, g_0_zzzzzz_0_xxxxxxx_1, g_0_zzzzzz_0_xxxxxxy_0, g_0_zzzzzz_0_xxxxxxy_1, g_0_zzzzzz_0_xxxxxxz_0, g_0_zzzzzz_0_xxxxxxz_1, g_0_zzzzzz_0_xxxxxy_1, g_0_zzzzzz_0_xxxxxyy_0, g_0_zzzzzz_0_xxxxxyy_1, g_0_zzzzzz_0_xxxxxyz_0, g_0_zzzzzz_0_xxxxxyz_1, g_0_zzzzzz_0_xxxxxz_1, g_0_zzzzzz_0_xxxxxzz_0, g_0_zzzzzz_0_xxxxxzz_1, g_0_zzzzzz_0_xxxxyy_1, g_0_zzzzzz_0_xxxxyyy_0, g_0_zzzzzz_0_xxxxyyy_1, g_0_zzzzzz_0_xxxxyyz_0, g_0_zzzzzz_0_xxxxyyz_1, g_0_zzzzzz_0_xxxxyz_1, g_0_zzzzzz_0_xxxxyzz_0, g_0_zzzzzz_0_xxxxyzz_1, g_0_zzzzzz_0_xxxxzz_1, g_0_zzzzzz_0_xxxxzzz_0, g_0_zzzzzz_0_xxxxzzz_1, g_0_zzzzzz_0_xxxyyy_1, g_0_zzzzzz_0_xxxyyyy_0, g_0_zzzzzz_0_xxxyyyy_1, g_0_zzzzzz_0_xxxyyyz_0, g_0_zzzzzz_0_xxxyyyz_1, g_0_zzzzzz_0_xxxyyz_1, g_0_zzzzzz_0_xxxyyzz_0, g_0_zzzzzz_0_xxxyyzz_1, g_0_zzzzzz_0_xxxyzz_1, g_0_zzzzzz_0_xxxyzzz_0, g_0_zzzzzz_0_xxxyzzz_1, g_0_zzzzzz_0_xxxzzz_1, g_0_zzzzzz_0_xxxzzzz_0, g_0_zzzzzz_0_xxxzzzz_1, g_0_zzzzzz_0_xxyyyy_1, g_0_zzzzzz_0_xxyyyyy_0, g_0_zzzzzz_0_xxyyyyy_1, g_0_zzzzzz_0_xxyyyyz_0, g_0_zzzzzz_0_xxyyyyz_1, g_0_zzzzzz_0_xxyyyz_1, g_0_zzzzzz_0_xxyyyzz_0, g_0_zzzzzz_0_xxyyyzz_1, g_0_zzzzzz_0_xxyyzz_1, g_0_zzzzzz_0_xxyyzzz_0, g_0_zzzzzz_0_xxyyzzz_1, g_0_zzzzzz_0_xxyzzz_1, g_0_zzzzzz_0_xxyzzzz_0, g_0_zzzzzz_0_xxyzzzz_1, g_0_zzzzzz_0_xxzzzz_1, g_0_zzzzzz_0_xxzzzzz_0, g_0_zzzzzz_0_xxzzzzz_1, g_0_zzzzzz_0_xyyyyy_1, g_0_zzzzzz_0_xyyyyyy_0, g_0_zzzzzz_0_xyyyyyy_1, g_0_zzzzzz_0_xyyyyyz_0, g_0_zzzzzz_0_xyyyyyz_1, g_0_zzzzzz_0_xyyyyz_1, g_0_zzzzzz_0_xyyyyzz_0, g_0_zzzzzz_0_xyyyyzz_1, g_0_zzzzzz_0_xyyyzz_1, g_0_zzzzzz_0_xyyyzzz_0, g_0_zzzzzz_0_xyyyzzz_1, g_0_zzzzzz_0_xyyzzz_1, g_0_zzzzzz_0_xyyzzzz_0, g_0_zzzzzz_0_xyyzzzz_1, g_0_zzzzzz_0_xyzzzz_1, g_0_zzzzzz_0_xyzzzzz_0, g_0_zzzzzz_0_xyzzzzz_1, g_0_zzzzzz_0_xzzzzz_1, g_0_zzzzzz_0_xzzzzzz_0, g_0_zzzzzz_0_xzzzzzz_1, g_0_zzzzzz_0_yyyyyy_1, g_0_zzzzzz_0_yyyyyyy_0, g_0_zzzzzz_0_yyyyyyy_1, g_0_zzzzzz_0_yyyyyyz_0, g_0_zzzzzz_0_yyyyyyz_1, g_0_zzzzzz_0_yyyyyz_1, g_0_zzzzzz_0_yyyyyzz_0, g_0_zzzzzz_0_yyyyyzz_1, g_0_zzzzzz_0_yyyyzz_1, g_0_zzzzzz_0_yyyyzzz_0, g_0_zzzzzz_0_yyyyzzz_1, g_0_zzzzzz_0_yyyzzz_1, g_0_zzzzzz_0_yyyzzzz_0, g_0_zzzzzz_0_yyyzzzz_1, g_0_zzzzzz_0_yyzzzz_1, g_0_zzzzzz_0_yyzzzzz_0, g_0_zzzzzz_0_yyzzzzz_1, g_0_zzzzzz_0_yzzzzz_1, g_0_zzzzzz_0_yzzzzzz_0, g_0_zzzzzz_0_yzzzzzz_1, g_0_zzzzzz_0_zzzzzz_1, g_0_zzzzzz_0_zzzzzzz_0, g_0_zzzzzz_0_zzzzzzz_1, g_0_zzzzzzz_0_xxxxxxx_0, g_0_zzzzzzz_0_xxxxxxy_0, g_0_zzzzzzz_0_xxxxxxz_0, g_0_zzzzzzz_0_xxxxxyy_0, g_0_zzzzzzz_0_xxxxxyz_0, g_0_zzzzzzz_0_xxxxxzz_0, g_0_zzzzzzz_0_xxxxyyy_0, g_0_zzzzzzz_0_xxxxyyz_0, g_0_zzzzzzz_0_xxxxyzz_0, g_0_zzzzzzz_0_xxxxzzz_0, g_0_zzzzzzz_0_xxxyyyy_0, g_0_zzzzzzz_0_xxxyyyz_0, g_0_zzzzzzz_0_xxxyyzz_0, g_0_zzzzzzz_0_xxxyzzz_0, g_0_zzzzzzz_0_xxxzzzz_0, g_0_zzzzzzz_0_xxyyyyy_0, g_0_zzzzzzz_0_xxyyyyz_0, g_0_zzzzzzz_0_xxyyyzz_0, g_0_zzzzzzz_0_xxyyzzz_0, g_0_zzzzzzz_0_xxyzzzz_0, g_0_zzzzzzz_0_xxzzzzz_0, g_0_zzzzzzz_0_xyyyyyy_0, g_0_zzzzzzz_0_xyyyyyz_0, g_0_zzzzzzz_0_xyyyyzz_0, g_0_zzzzzzz_0_xyyyzzz_0, g_0_zzzzzzz_0_xyyzzzz_0, g_0_zzzzzzz_0_xyzzzzz_0, g_0_zzzzzzz_0_xzzzzzz_0, g_0_zzzzzzz_0_yyyyyyy_0, g_0_zzzzzzz_0_yyyyyyz_0, g_0_zzzzzzz_0_yyyyyzz_0, g_0_zzzzzzz_0_yyyyzzz_0, g_0_zzzzzzz_0_yyyzzzz_0, g_0_zzzzzzz_0_yyzzzzz_0, g_0_zzzzzzz_0_yzzzzzz_0, g_0_zzzzzzz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzz_0_xxxxxxx_0[i] = 6.0 * g_0_zzzzz_0_xxxxxxx_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxxxx_0[i] * pb_z + g_0_zzzzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_zzzzz_0_xxxxxxy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxxxy_0[i] * pb_z + g_0_zzzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxxxz_0[i] = 6.0 * g_0_zzzzz_0_xxxxxxz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxxz_0[i] * pb_z + g_0_zzzzzz_0_xxxxxxz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxxyy_0[i] = 6.0 * g_0_zzzzz_0_xxxxxyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxxyy_0[i] * pb_z + g_0_zzzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxxyz_0[i] = 6.0 * g_0_zzzzz_0_xxxxxyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxyz_0[i] * pb_z + g_0_zzzzzz_0_xxxxxyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxxzz_0[i] = 6.0 * g_0_zzzzz_0_xxxxxzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxzz_0[i] * pb_z + g_0_zzzzzz_0_xxxxxzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxyyy_0[i] = 6.0 * g_0_zzzzz_0_xxxxyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxyyy_0[i] * pb_z + g_0_zzzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxyyz_0[i] = 6.0 * g_0_zzzzz_0_xxxxyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyyz_0[i] * pb_z + g_0_zzzzzz_0_xxxxyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxyzz_0[i] = 6.0 * g_0_zzzzz_0_xxxxyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyzz_0[i] * pb_z + g_0_zzzzzz_0_xxxxyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxzzz_0[i] = 6.0 * g_0_zzzzz_0_xxxxzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxzzz_0[i] * pb_z + g_0_zzzzzz_0_xxxxzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyyyy_0[i] = 6.0 * g_0_zzzzz_0_xxxyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxyyyy_0[i] * pb_z + g_0_zzzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyyyz_0[i] = 6.0 * g_0_zzzzz_0_xxxyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyyz_0[i] * pb_z + g_0_zzzzzz_0_xxxyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyyzz_0[i] = 6.0 * g_0_zzzzz_0_xxxyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyzz_0[i] * pb_z + g_0_zzzzzz_0_xxxyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyzzz_0[i] = 6.0 * g_0_zzzzz_0_xxxyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyzzz_0[i] * pb_z + g_0_zzzzzz_0_xxxyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxzzzz_0[i] = 6.0 * g_0_zzzzz_0_xxxzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxzzzz_0[i] * pb_z + g_0_zzzzzz_0_xxxzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyyyy_0[i] = 6.0 * g_0_zzzzz_0_xxyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxyyyyy_0[i] * pb_z + g_0_zzzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyyyz_0[i] = 6.0 * g_0_zzzzz_0_xxyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyyz_0[i] * pb_z + g_0_zzzzzz_0_xxyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyyzz_0[i] = 6.0 * g_0_zzzzz_0_xxyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyzz_0[i] * pb_z + g_0_zzzzzz_0_xxyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyzzz_0[i] = 6.0 * g_0_zzzzz_0_xxyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyzzz_0[i] * pb_z + g_0_zzzzzz_0_xxyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyzzzz_0[i] = 6.0 * g_0_zzzzz_0_xxyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzzzz_0[i] * pb_z + g_0_zzzzzz_0_xxyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxzzzzz_0[i] = 6.0 * g_0_zzzzz_0_xxzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxzzzzz_0[i] * pb_z + g_0_zzzzzz_0_xxzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyyyy_0[i] = 6.0 * g_0_zzzzz_0_xyyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_xyyyyyy_0[i] * pb_z + g_0_zzzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyyyz_0[i] = 6.0 * g_0_zzzzz_0_xyyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyyz_0[i] * pb_z + g_0_zzzzzz_0_xyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyyzz_0[i] = 6.0 * g_0_zzzzz_0_xyyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyzz_0[i] * pb_z + g_0_zzzzzz_0_xyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyzzz_0[i] = 6.0 * g_0_zzzzz_0_xyyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyzzz_0[i] * pb_z + g_0_zzzzzz_0_xyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyzzzz_0[i] = 6.0 * g_0_zzzzz_0_xyyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzzzz_0[i] * pb_z + g_0_zzzzzz_0_xyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyzzzzz_0[i] = 6.0 * g_0_zzzzz_0_xyzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzzzz_0[i] * pb_z + g_0_zzzzzz_0_xyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xzzzzzz_0[i] = 6.0 * g_0_zzzzz_0_xzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzzzzzz_0[i] * pb_z + g_0_zzzzzz_0_xzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyyyy_0[i] = 6.0 * g_0_zzzzz_0_yyyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_zzzzzz_0_yyyyyyy_0[i] * pb_z + g_0_zzzzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyyyz_0[i] = 6.0 * g_0_zzzzz_0_yyyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_zzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyyyz_0[i] * pb_z + g_0_zzzzzz_0_yyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyyzz_0[i] = 6.0 * g_0_zzzzz_0_yyyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyyzz_0[i] * pb_z + g_0_zzzzzz_0_yyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyzzz_0[i] = 6.0 * g_0_zzzzz_0_yyyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyzzz_0[i] * pb_z + g_0_zzzzzz_0_yyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyzzzz_0[i] = 6.0 * g_0_zzzzz_0_yyyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyzzzz_0[i] * pb_z + g_0_zzzzzz_0_yyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyzzzzz_0[i] = 6.0 * g_0_zzzzz_0_yyzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyzzzzz_0[i] * pb_z + g_0_zzzzzz_0_yyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yzzzzzz_0[i] = 6.0 * g_0_zzzzz_0_yzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzzzzzz_0[i] * pb_z + g_0_zzzzzz_0_yzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_zzzzzzz_0[i] = 6.0 * g_0_zzzzz_0_zzzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_zzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_zzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_zzzzzzz_0[i] * pb_z + g_0_zzzzzz_0_zzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

