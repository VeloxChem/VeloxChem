#include "ElectronRepulsionPrimRecSLSK.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_slsk(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_slsk,
                                  size_t idx_eri_0_sisk,
                                  size_t idx_eri_1_sisk,
                                  size_t idx_eri_1_sksi,
                                  size_t idx_eri_0_sksk,
                                  size_t idx_eri_1_sksk,
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

    auto g_0_xxxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 38);

    auto g_0_xxxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 41);

    auto g_0_xxxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 45);

    auto g_0_xxxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 50);

    auto g_0_xxxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 56);

    auto g_0_xxxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 63);

    auto g_0_xxxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sisk + 72);

    auto g_0_xxxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sisk + 73);

    auto g_0_xxxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 75);

    auto g_0_xxxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 78);

    auto g_0_xxxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 82);

    auto g_0_xxxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 87);

    auto g_0_xxxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 93);

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

    auto g_0_yyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sisk + 795);

    auto g_0_yyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sisk + 798);

    auto g_0_yyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sisk + 802);

    auto g_0_yyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 807);

    auto g_0_yyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 813);

    auto g_0_yyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sisk + 820);

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

    auto g_0_yzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sisk + 938);

    auto g_0_yzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sisk + 940);

    auto g_0_yzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sisk + 941);

    auto g_0_yzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sisk + 943);

    auto g_0_yzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sisk + 944);

    auto g_0_yzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sisk + 945);

    auto g_0_yzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sisk + 947);

    auto g_0_yzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sisk + 948);

    auto g_0_yzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sisk + 949);

    auto g_0_yzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sisk + 950);

    auto g_0_yzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 952);

    auto g_0_yzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 953);

    auto g_0_yzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 954);

    auto g_0_yzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 955);

    auto g_0_yzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 956);

    auto g_0_yzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sisk + 958);

    auto g_0_yzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sisk + 959);

    auto g_0_yzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sisk + 960);

    auto g_0_yzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sisk + 961);

    auto g_0_yzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 962);

    auto g_0_yzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sisk + 963);

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

    auto g_0_xxxxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 38);

    auto g_0_xxxxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 41);

    auto g_0_xxxxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 45);

    auto g_0_xxxxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 50);

    auto g_0_xxxxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 56);

    auto g_0_xxxxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 63);

    auto g_0_xxxxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sisk + 72);

    auto g_0_xxxxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sisk + 73);

    auto g_0_xxxxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 75);

    auto g_0_xxxxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 78);

    auto g_0_xxxxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 82);

    auto g_0_xxxxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 87);

    auto g_0_xxxxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 93);

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

    auto g_0_yyyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sisk + 795);

    auto g_0_yyyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sisk + 798);

    auto g_0_yyyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sisk + 802);

    auto g_0_yyyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 807);

    auto g_0_yyyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 813);

    auto g_0_yyyyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sisk + 820);

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

    auto g_0_yzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sisk + 938);

    auto g_0_yzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sisk + 940);

    auto g_0_yzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sisk + 941);

    auto g_0_yzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sisk + 943);

    auto g_0_yzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sisk + 944);

    auto g_0_yzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sisk + 945);

    auto g_0_yzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sisk + 947);

    auto g_0_yzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sisk + 948);

    auto g_0_yzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sisk + 949);

    auto g_0_yzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sisk + 950);

    auto g_0_yzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 952);

    auto g_0_yzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 953);

    auto g_0_yzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 954);

    auto g_0_yzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 955);

    auto g_0_yzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 956);

    auto g_0_yzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sisk + 958);

    auto g_0_yzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sisk + 959);

    auto g_0_yzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sisk + 960);

    auto g_0_yzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sisk + 961);

    auto g_0_yzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 962);

    auto g_0_yzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sisk + 963);

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

    auto g_0_xxxxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 58);

    auto g_0_xxxxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 60);

    auto g_0_xxxxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 61);

    auto g_0_xxxxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 63);

    auto g_0_xxxxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 64);

    auto g_0_xxxxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 65);

    auto g_0_xxxxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 67);

    auto g_0_xxxxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 68);

    auto g_0_xxxxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 69);

    auto g_0_xxxxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 70);

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

    auto g_0_xxxyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 340);

    auto g_0_xxxyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 343);

    auto g_0_xxxyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 344);

    auto g_0_xxxyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 347);

    auto g_0_xxxyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 348);

    auto g_0_xxxyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 349);

    auto g_0_xxxyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 352);

    auto g_0_xxxyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 353);

    auto g_0_xxxyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 354);

    auto g_0_xxxyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 355);

    auto g_0_xxxyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 358);

    auto g_0_xxxyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 359);

    auto g_0_xxxyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 360);

    auto g_0_xxxyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 361);

    auto g_0_xxxyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 362);

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

    auto g_0_xxyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 480);

    auto g_0_xxyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 483);

    auto g_0_xxyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 484);

    auto g_0_xxyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 487);

    auto g_0_xxyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 488);

    auto g_0_xxyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 489);

    auto g_0_xxyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 492);

    auto g_0_xxyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 493);

    auto g_0_xxyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 494);

    auto g_0_xxyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 495);

    auto g_0_xxyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 498);

    auto g_0_xxyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 499);

    auto g_0_xxyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 500);

    auto g_0_xxyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 501);

    auto g_0_xxyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 502);

    auto g_0_xxyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 508);

    auto g_0_xxyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 511);

    auto g_0_xxyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 512);

    auto g_0_xxyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 515);

    auto g_0_xxyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 516);

    auto g_0_xxyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 517);

    auto g_0_xxyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 520);

    auto g_0_xxyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 521);

    auto g_0_xxyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 522);

    auto g_0_xxyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 523);

    auto g_0_xxyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 526);

    auto g_0_xxyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 527);

    auto g_0_xxyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 528);

    auto g_0_xxyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 529);

    auto g_0_xxyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 530);

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

    auto g_0_xyyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 666);

    auto g_0_xyyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 667);

    auto g_0_xyyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 668);

    auto g_0_xyyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 669);

    auto g_0_xyyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 670);

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

    auto g_0_xyyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 694);

    auto g_0_xyyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 695);

    auto g_0_xyyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 696);

    auto g_0_xyyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 697);

    auto g_0_xyyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 698);

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

    auto g_0_xyyzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 722);

    auto g_0_xyyzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 723);

    auto g_0_xyyzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 724);

    auto g_0_xyyzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 725);

    auto g_0_xyyzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 726);

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

    auto g_0_yyyyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 814);

    auto g_0_yyyyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 816);

    auto g_0_yyyyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 817);

    auto g_0_yyyyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 819);

    auto g_0_yyyyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 820);

    auto g_0_yyyyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 821);

    auto g_0_yyyyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 823);

    auto g_0_yyyyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 824);

    auto g_0_yyyyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 825);

    auto g_0_yyyyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 826);

    auto g_0_yyyyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 828);

    auto g_0_yyyyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 829);

    auto g_0_yyyyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 830);

    auto g_0_yyyyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 831);

    auto g_0_yyyyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 832);

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

    /// Set up components of auxilary buffer : SKSK

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

    auto g_0_xxxxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 36);

    auto g_0_xxxxxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 37);

    auto g_0_xxxxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 38);

    auto g_0_xxxxxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 39);

    auto g_0_xxxxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 41);

    auto g_0_xxxxxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 42);

    auto g_0_xxxxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 45);

    auto g_0_xxxxxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 46);

    auto g_0_xxxxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 50);

    auto g_0_xxxxxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 51);

    auto g_0_xxxxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 56);

    auto g_0_xxxxxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 57);

    auto g_0_xxxxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 63);

    auto g_0_xxxxxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 64);

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

    auto g_0_xxxxxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 101);

    auto g_0_xxxxxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 102);

    auto g_0_xxxxxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 103);

    auto g_0_xxxxxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 104);

    auto g_0_xxxxxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 105);

    auto g_0_xxxxxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 106);

    auto g_0_xxxxxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 107);

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

    auto g_0_xxxxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 253);

    auto g_0_xxxxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 255);

    auto g_0_xxxxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 258);

    auto g_0_xxxxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 262);

    auto g_0_xxxxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 267);

    auto g_0_xxxxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 273);

    auto g_0_xxxxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 288);

    auto g_0_xxxxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 290);

    auto g_0_xxxxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 293);

    auto g_0_xxxxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 297);

    auto g_0_xxxxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 302);

    auto g_0_xxxxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 308);

    auto g_0_xxxxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 315);

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

    auto g_0_xxxyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 397);

    auto g_0_xxxyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 399);

    auto g_0_xxxyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 402);

    auto g_0_xxxyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 406);

    auto g_0_xxxyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 411);

    auto g_0_xxxyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 417);

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

    auto g_0_xxxyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 468);

    auto g_0_xxxyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 470);

    auto g_0_xxxyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 473);

    auto g_0_xxxyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 477);

    auto g_0_xxxyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 482);

    auto g_0_xxxyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 488);

    auto g_0_xxxyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 495);

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

    auto g_0_xxyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 577);

    auto g_0_xxyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 579);

    auto g_0_xxyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 582);

    auto g_0_xxyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 586);

    auto g_0_xxyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 591);

    auto g_0_xxyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 597);

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

    auto g_0_xxyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 684);

    auto g_0_xxyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 686);

    auto g_0_xxyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 689);

    auto g_0_xxyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 693);

    auto g_0_xxyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 698);

    auto g_0_xxyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 704);

    auto g_0_xxyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 711);

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

    auto g_0_xyyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 756);

    auto g_0_xyyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sksk + 757);

    auto g_0_xyyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sksk + 759);

    auto g_0_xyyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 760);

    auto g_0_xyyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sksk + 762);

    auto g_0_xyyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 763);

    auto g_0_xyyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 764);

    auto g_0_xyyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sksk + 766);

    auto g_0_xyyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 767);

    auto g_0_xyyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 768);

    auto g_0_xyyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 769);

    auto g_0_xyyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 771);

    auto g_0_xyyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 772);

    auto g_0_xyyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 773);

    auto g_0_xyyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 774);

    auto g_0_xyyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 775);

    auto g_0_xyyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 777);

    auto g_0_xyyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 778);

    auto g_0_xyyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 779);

    auto g_0_xyyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 780);

    auto g_0_xyyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 781);

    auto g_0_xyyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 782);

    auto g_0_xyyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 784);

    auto g_0_xyyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 785);

    auto g_0_xyyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 786);

    auto g_0_xyyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 787);

    auto g_0_xyyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 788);

    auto g_0_xyyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 789);

    auto g_0_xyyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 790);

    auto g_0_xyyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 791);

    auto g_0_xyyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 832);

    auto g_0_xyyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 835);

    auto g_0_xyyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 836);

    auto g_0_xyyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 839);

    auto g_0_xyyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 840);

    auto g_0_xyyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 841);

    auto g_0_xyyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 844);

    auto g_0_xyyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 845);

    auto g_0_xyyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 846);

    auto g_0_xyyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 847);

    auto g_0_xyyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 850);

    auto g_0_xyyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 851);

    auto g_0_xyyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 852);

    auto g_0_xyyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 853);

    auto g_0_xyyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 854);

    auto g_0_xyyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 856);

    auto g_0_xyyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 857);

    auto g_0_xyyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 858);

    auto g_0_xyyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 859);

    auto g_0_xyyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 860);

    auto g_0_xyyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 861);

    auto g_0_xyyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 862);

    auto g_0_xyyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 863);

    auto g_0_xyyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 868);

    auto g_0_xyyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 871);

    auto g_0_xyyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 872);

    auto g_0_xyyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 875);

    auto g_0_xyyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 876);

    auto g_0_xyyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 877);

    auto g_0_xyyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 880);

    auto g_0_xyyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 881);

    auto g_0_xyyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 882);

    auto g_0_xyyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 883);

    auto g_0_xyyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 886);

    auto g_0_xyyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 887);

    auto g_0_xyyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 888);

    auto g_0_xyyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 889);

    auto g_0_xyyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 890);

    auto g_0_xyyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 892);

    auto g_0_xyyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 893);

    auto g_0_xyyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 894);

    auto g_0_xyyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 895);

    auto g_0_xyyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 896);

    auto g_0_xyyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 897);

    auto g_0_xyyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 898);

    auto g_0_xyyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 899);

    auto g_0_xyyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 904);

    auto g_0_xyyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 907);

    auto g_0_xyyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 908);

    auto g_0_xyyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 911);

    auto g_0_xyyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 912);

    auto g_0_xyyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 913);

    auto g_0_xyyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 916);

    auto g_0_xyyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 917);

    auto g_0_xyyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 918);

    auto g_0_xyyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 919);

    auto g_0_xyyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 922);

    auto g_0_xyyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 923);

    auto g_0_xyyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 924);

    auto g_0_xyyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 925);

    auto g_0_xyyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 926);

    auto g_0_xyyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sksk + 928);

    auto g_0_xyyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 929);

    auto g_0_xyyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 930);

    auto g_0_xyyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 931);

    auto g_0_xyyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 932);

    auto g_0_xyyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 933);

    auto g_0_xyyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 934);

    auto g_0_xyyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 935);

    auto g_0_xzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sksk + 972);

    auto g_0_xzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sksk + 974);

    auto g_0_xzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sksk + 976);

    auto g_0_xzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sksk + 977);

    auto g_0_xzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sksk + 979);

    auto g_0_xzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sksk + 980);

    auto g_0_xzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sksk + 981);

    auto g_0_xzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sksk + 983);

    auto g_0_xzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sksk + 984);

    auto g_0_xzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sksk + 985);

    auto g_0_xzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sksk + 986);

    auto g_0_xzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sksk + 988);

    auto g_0_xzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sksk + 989);

    auto g_0_xzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sksk + 990);

    auto g_0_xzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sksk + 991);

    auto g_0_xzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sksk + 992);

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

    /// Set up components of auxilary buffer : SKSK

    auto g_0_xxxxxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk);

    auto g_0_xxxxxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1);

    auto g_0_xxxxxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 2);

    auto g_0_xxxxxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 3);

    auto g_0_xxxxxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 4);

    auto g_0_xxxxxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 5);

    auto g_0_xxxxxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 6);

    auto g_0_xxxxxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 7);

    auto g_0_xxxxxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 8);

    auto g_0_xxxxxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 9);

    auto g_0_xxxxxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 10);

    auto g_0_xxxxxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 11);

    auto g_0_xxxxxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 12);

    auto g_0_xxxxxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 13);

    auto g_0_xxxxxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 14);

    auto g_0_xxxxxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 15);

    auto g_0_xxxxxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 16);

    auto g_0_xxxxxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 17);

    auto g_0_xxxxxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 18);

    auto g_0_xxxxxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 19);

    auto g_0_xxxxxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 20);

    auto g_0_xxxxxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 21);

    auto g_0_xxxxxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 22);

    auto g_0_xxxxxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 23);

    auto g_0_xxxxxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 24);

    auto g_0_xxxxxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 25);

    auto g_0_xxxxxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 26);

    auto g_0_xxxxxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 27);

    auto g_0_xxxxxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 28);

    auto g_0_xxxxxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 29);

    auto g_0_xxxxxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 30);

    auto g_0_xxxxxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 31);

    auto g_0_xxxxxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 32);

    auto g_0_xxxxxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 33);

    auto g_0_xxxxxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 34);

    auto g_0_xxxxxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 35);

    auto g_0_xxxxxxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 36);

    auto g_0_xxxxxxy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 37);

    auto g_0_xxxxxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 38);

    auto g_0_xxxxxxy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 39);

    auto g_0_xxxxxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 41);

    auto g_0_xxxxxxy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 42);

    auto g_0_xxxxxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 45);

    auto g_0_xxxxxxy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 46);

    auto g_0_xxxxxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 50);

    auto g_0_xxxxxxy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 51);

    auto g_0_xxxxxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 56);

    auto g_0_xxxxxxy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 57);

    auto g_0_xxxxxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 63);

    auto g_0_xxxxxxy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 64);

    auto g_0_xxxxxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 72);

    auto g_0_xxxxxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 73);

    auto g_0_xxxxxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 74);

    auto g_0_xxxxxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 75);

    auto g_0_xxxxxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 76);

    auto g_0_xxxxxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 77);

    auto g_0_xxxxxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 78);

    auto g_0_xxxxxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 79);

    auto g_0_xxxxxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 80);

    auto g_0_xxxxxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 81);

    auto g_0_xxxxxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 82);

    auto g_0_xxxxxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 83);

    auto g_0_xxxxxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 84);

    auto g_0_xxxxxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 85);

    auto g_0_xxxxxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 86);

    auto g_0_xxxxxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 87);

    auto g_0_xxxxxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 88);

    auto g_0_xxxxxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 89);

    auto g_0_xxxxxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 90);

    auto g_0_xxxxxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 91);

    auto g_0_xxxxxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 92);

    auto g_0_xxxxxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 93);

    auto g_0_xxxxxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 94);

    auto g_0_xxxxxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 95);

    auto g_0_xxxxxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 96);

    auto g_0_xxxxxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 97);

    auto g_0_xxxxxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 98);

    auto g_0_xxxxxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 99);

    auto g_0_xxxxxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 101);

    auto g_0_xxxxxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 102);

    auto g_0_xxxxxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 103);

    auto g_0_xxxxxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 104);

    auto g_0_xxxxxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 105);

    auto g_0_xxxxxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 106);

    auto g_0_xxxxxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 107);

    auto g_0_xxxxxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 108);

    auto g_0_xxxxxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 109);

    auto g_0_xxxxxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 110);

    auto g_0_xxxxxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 111);

    auto g_0_xxxxxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 112);

    auto g_0_xxxxxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 113);

    auto g_0_xxxxxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 114);

    auto g_0_xxxxxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 115);

    auto g_0_xxxxxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 116);

    auto g_0_xxxxxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 117);

    auto g_0_xxxxxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 118);

    auto g_0_xxxxxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 119);

    auto g_0_xxxxxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 120);

    auto g_0_xxxxxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 121);

    auto g_0_xxxxxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 122);

    auto g_0_xxxxxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 123);

    auto g_0_xxxxxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 124);

    auto g_0_xxxxxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 125);

    auto g_0_xxxxxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 126);

    auto g_0_xxxxxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 127);

    auto g_0_xxxxxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 128);

    auto g_0_xxxxxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 129);

    auto g_0_xxxxxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 130);

    auto g_0_xxxxxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 131);

    auto g_0_xxxxxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 132);

    auto g_0_xxxxxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 133);

    auto g_0_xxxxxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 134);

    auto g_0_xxxxxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 135);

    auto g_0_xxxxxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 136);

    auto g_0_xxxxxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 137);

    auto g_0_xxxxxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 138);

    auto g_0_xxxxxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 139);

    auto g_0_xxxxxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 140);

    auto g_0_xxxxxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 141);

    auto g_0_xxxxxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 142);

    auto g_0_xxxxxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 143);

    auto g_0_xxxxxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 180);

    auto g_0_xxxxxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 181);

    auto g_0_xxxxxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 182);

    auto g_0_xxxxxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 183);

    auto g_0_xxxxxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 184);

    auto g_0_xxxxxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 185);

    auto g_0_xxxxxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 186);

    auto g_0_xxxxxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 187);

    auto g_0_xxxxxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 188);

    auto g_0_xxxxxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 189);

    auto g_0_xxxxxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 190);

    auto g_0_xxxxxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 191);

    auto g_0_xxxxxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 192);

    auto g_0_xxxxxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 193);

    auto g_0_xxxxxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 194);

    auto g_0_xxxxxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 195);

    auto g_0_xxxxxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 196);

    auto g_0_xxxxxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 197);

    auto g_0_xxxxxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 198);

    auto g_0_xxxxxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 199);

    auto g_0_xxxxxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 200);

    auto g_0_xxxxxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 201);

    auto g_0_xxxxxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 202);

    auto g_0_xxxxxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 203);

    auto g_0_xxxxxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 204);

    auto g_0_xxxxxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 205);

    auto g_0_xxxxxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 206);

    auto g_0_xxxxxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 207);

    auto g_0_xxxxxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 208);

    auto g_0_xxxxxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 209);

    auto g_0_xxxxxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 210);

    auto g_0_xxxxxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 211);

    auto g_0_xxxxxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 212);

    auto g_0_xxxxxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 213);

    auto g_0_xxxxxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 214);

    auto g_0_xxxxxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 215);

    auto g_0_xxxxyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 216);

    auto g_0_xxxxyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 217);

    auto g_0_xxxxyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 218);

    auto g_0_xxxxyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 219);

    auto g_0_xxxxyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 220);

    auto g_0_xxxxyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 221);

    auto g_0_xxxxyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 222);

    auto g_0_xxxxyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 223);

    auto g_0_xxxxyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 224);

    auto g_0_xxxxyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 225);

    auto g_0_xxxxyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 226);

    auto g_0_xxxxyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 227);

    auto g_0_xxxxyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 228);

    auto g_0_xxxxyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 229);

    auto g_0_xxxxyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 230);

    auto g_0_xxxxyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 231);

    auto g_0_xxxxyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 232);

    auto g_0_xxxxyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 233);

    auto g_0_xxxxyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 234);

    auto g_0_xxxxyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 235);

    auto g_0_xxxxyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 236);

    auto g_0_xxxxyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 237);

    auto g_0_xxxxyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 238);

    auto g_0_xxxxyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 239);

    auto g_0_xxxxyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 240);

    auto g_0_xxxxyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 241);

    auto g_0_xxxxyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 242);

    auto g_0_xxxxyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 243);

    auto g_0_xxxxyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 244);

    auto g_0_xxxxyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 245);

    auto g_0_xxxxyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 246);

    auto g_0_xxxxyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 247);

    auto g_0_xxxxyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 248);

    auto g_0_xxxxyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 249);

    auto g_0_xxxxyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 250);

    auto g_0_xxxxyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 251);

    auto g_0_xxxxyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 253);

    auto g_0_xxxxyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 255);

    auto g_0_xxxxyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 258);

    auto g_0_xxxxyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 262);

    auto g_0_xxxxyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 267);

    auto g_0_xxxxyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 273);

    auto g_0_xxxxyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 288);

    auto g_0_xxxxyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 290);

    auto g_0_xxxxyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 293);

    auto g_0_xxxxyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 297);

    auto g_0_xxxxyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 302);

    auto g_0_xxxxyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 308);

    auto g_0_xxxxyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 315);

    auto g_0_xxxxzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 324);

    auto g_0_xxxxzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 325);

    auto g_0_xxxxzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 326);

    auto g_0_xxxxzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 327);

    auto g_0_xxxxzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 328);

    auto g_0_xxxxzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 329);

    auto g_0_xxxxzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 330);

    auto g_0_xxxxzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 331);

    auto g_0_xxxxzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 332);

    auto g_0_xxxxzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 333);

    auto g_0_xxxxzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 334);

    auto g_0_xxxxzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 335);

    auto g_0_xxxxzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 336);

    auto g_0_xxxxzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 337);

    auto g_0_xxxxzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 338);

    auto g_0_xxxxzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 339);

    auto g_0_xxxxzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 340);

    auto g_0_xxxxzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 341);

    auto g_0_xxxxzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 342);

    auto g_0_xxxxzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 343);

    auto g_0_xxxxzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 344);

    auto g_0_xxxxzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 345);

    auto g_0_xxxxzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 346);

    auto g_0_xxxxzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 347);

    auto g_0_xxxxzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 348);

    auto g_0_xxxxzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 349);

    auto g_0_xxxxzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 350);

    auto g_0_xxxxzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 351);

    auto g_0_xxxxzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 352);

    auto g_0_xxxxzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 353);

    auto g_0_xxxxzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 354);

    auto g_0_xxxxzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 355);

    auto g_0_xxxxzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 356);

    auto g_0_xxxxzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 357);

    auto g_0_xxxxzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 358);

    auto g_0_xxxxzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 359);

    auto g_0_xxxyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 360);

    auto g_0_xxxyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 361);

    auto g_0_xxxyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 362);

    auto g_0_xxxyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 363);

    auto g_0_xxxyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 364);

    auto g_0_xxxyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 365);

    auto g_0_xxxyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 366);

    auto g_0_xxxyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 367);

    auto g_0_xxxyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 368);

    auto g_0_xxxyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 369);

    auto g_0_xxxyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 370);

    auto g_0_xxxyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 371);

    auto g_0_xxxyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 372);

    auto g_0_xxxyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 373);

    auto g_0_xxxyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 374);

    auto g_0_xxxyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 375);

    auto g_0_xxxyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 376);

    auto g_0_xxxyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 377);

    auto g_0_xxxyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 378);

    auto g_0_xxxyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 379);

    auto g_0_xxxyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 380);

    auto g_0_xxxyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 381);

    auto g_0_xxxyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 382);

    auto g_0_xxxyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 383);

    auto g_0_xxxyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 384);

    auto g_0_xxxyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 385);

    auto g_0_xxxyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 386);

    auto g_0_xxxyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 387);

    auto g_0_xxxyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 388);

    auto g_0_xxxyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 389);

    auto g_0_xxxyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 390);

    auto g_0_xxxyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 391);

    auto g_0_xxxyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 392);

    auto g_0_xxxyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 393);

    auto g_0_xxxyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 394);

    auto g_0_xxxyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 395);

    auto g_0_xxxyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 397);

    auto g_0_xxxyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 399);

    auto g_0_xxxyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 402);

    auto g_0_xxxyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 406);

    auto g_0_xxxyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 411);

    auto g_0_xxxyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 417);

    auto g_0_xxxyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 432);

    auto g_0_xxxyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 433);

    auto g_0_xxxyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 434);

    auto g_0_xxxyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 435);

    auto g_0_xxxyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 436);

    auto g_0_xxxyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 437);

    auto g_0_xxxyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 438);

    auto g_0_xxxyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 439);

    auto g_0_xxxyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 440);

    auto g_0_xxxyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 441);

    auto g_0_xxxyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 442);

    auto g_0_xxxyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 443);

    auto g_0_xxxyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 444);

    auto g_0_xxxyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 445);

    auto g_0_xxxyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 446);

    auto g_0_xxxyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 447);

    auto g_0_xxxyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 448);

    auto g_0_xxxyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 449);

    auto g_0_xxxyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 450);

    auto g_0_xxxyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 451);

    auto g_0_xxxyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 452);

    auto g_0_xxxyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 453);

    auto g_0_xxxyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 454);

    auto g_0_xxxyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 455);

    auto g_0_xxxyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 456);

    auto g_0_xxxyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 457);

    auto g_0_xxxyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 458);

    auto g_0_xxxyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 459);

    auto g_0_xxxyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 460);

    auto g_0_xxxyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 461);

    auto g_0_xxxyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 462);

    auto g_0_xxxyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 463);

    auto g_0_xxxyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 464);

    auto g_0_xxxyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 465);

    auto g_0_xxxyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 466);

    auto g_0_xxxyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 467);

    auto g_0_xxxyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 468);

    auto g_0_xxxyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 470);

    auto g_0_xxxyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 473);

    auto g_0_xxxyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 477);

    auto g_0_xxxyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 482);

    auto g_0_xxxyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 488);

    auto g_0_xxxyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 495);

    auto g_0_xxxzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 504);

    auto g_0_xxxzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 505);

    auto g_0_xxxzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 506);

    auto g_0_xxxzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 507);

    auto g_0_xxxzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 508);

    auto g_0_xxxzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 509);

    auto g_0_xxxzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 510);

    auto g_0_xxxzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 511);

    auto g_0_xxxzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 512);

    auto g_0_xxxzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 513);

    auto g_0_xxxzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 514);

    auto g_0_xxxzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 515);

    auto g_0_xxxzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 516);

    auto g_0_xxxzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 517);

    auto g_0_xxxzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 518);

    auto g_0_xxxzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 519);

    auto g_0_xxxzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 520);

    auto g_0_xxxzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 521);

    auto g_0_xxxzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 522);

    auto g_0_xxxzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 523);

    auto g_0_xxxzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 524);

    auto g_0_xxxzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 525);

    auto g_0_xxxzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 526);

    auto g_0_xxxzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 527);

    auto g_0_xxxzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 528);

    auto g_0_xxxzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 529);

    auto g_0_xxxzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 530);

    auto g_0_xxxzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 531);

    auto g_0_xxxzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 532);

    auto g_0_xxxzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 533);

    auto g_0_xxxzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 534);

    auto g_0_xxxzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 535);

    auto g_0_xxxzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 536);

    auto g_0_xxxzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 537);

    auto g_0_xxxzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 538);

    auto g_0_xxxzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 539);

    auto g_0_xxyyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 540);

    auto g_0_xxyyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 541);

    auto g_0_xxyyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 542);

    auto g_0_xxyyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 543);

    auto g_0_xxyyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 544);

    auto g_0_xxyyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 545);

    auto g_0_xxyyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 546);

    auto g_0_xxyyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 547);

    auto g_0_xxyyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 548);

    auto g_0_xxyyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 549);

    auto g_0_xxyyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 550);

    auto g_0_xxyyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 551);

    auto g_0_xxyyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 552);

    auto g_0_xxyyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 553);

    auto g_0_xxyyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 554);

    auto g_0_xxyyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 555);

    auto g_0_xxyyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 556);

    auto g_0_xxyyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 557);

    auto g_0_xxyyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 558);

    auto g_0_xxyyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 559);

    auto g_0_xxyyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 560);

    auto g_0_xxyyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 561);

    auto g_0_xxyyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 562);

    auto g_0_xxyyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 563);

    auto g_0_xxyyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 564);

    auto g_0_xxyyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 565);

    auto g_0_xxyyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 566);

    auto g_0_xxyyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 567);

    auto g_0_xxyyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 568);

    auto g_0_xxyyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 569);

    auto g_0_xxyyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 570);

    auto g_0_xxyyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 571);

    auto g_0_xxyyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 572);

    auto g_0_xxyyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 573);

    auto g_0_xxyyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 574);

    auto g_0_xxyyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 575);

    auto g_0_xxyyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 577);

    auto g_0_xxyyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 579);

    auto g_0_xxyyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 582);

    auto g_0_xxyyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 586);

    auto g_0_xxyyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 591);

    auto g_0_xxyyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 597);

    auto g_0_xxyyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 612);

    auto g_0_xxyyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 613);

    auto g_0_xxyyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 614);

    auto g_0_xxyyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 615);

    auto g_0_xxyyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 616);

    auto g_0_xxyyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 617);

    auto g_0_xxyyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 618);

    auto g_0_xxyyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 619);

    auto g_0_xxyyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 620);

    auto g_0_xxyyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 621);

    auto g_0_xxyyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 622);

    auto g_0_xxyyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 623);

    auto g_0_xxyyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 624);

    auto g_0_xxyyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 625);

    auto g_0_xxyyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 626);

    auto g_0_xxyyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 627);

    auto g_0_xxyyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 628);

    auto g_0_xxyyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 629);

    auto g_0_xxyyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 630);

    auto g_0_xxyyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 631);

    auto g_0_xxyyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 632);

    auto g_0_xxyyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 633);

    auto g_0_xxyyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 634);

    auto g_0_xxyyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 635);

    auto g_0_xxyyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 636);

    auto g_0_xxyyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 637);

    auto g_0_xxyyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 638);

    auto g_0_xxyyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 639);

    auto g_0_xxyyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 640);

    auto g_0_xxyyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 641);

    auto g_0_xxyyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 642);

    auto g_0_xxyyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 643);

    auto g_0_xxyyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 644);

    auto g_0_xxyyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 645);

    auto g_0_xxyyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 646);

    auto g_0_xxyyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 647);

    auto g_0_xxyyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 648);

    auto g_0_xxyyzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 649);

    auto g_0_xxyyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 650);

    auto g_0_xxyyzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 651);

    auto g_0_xxyyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 652);

    auto g_0_xxyyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 653);

    auto g_0_xxyyzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 654);

    auto g_0_xxyyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 655);

    auto g_0_xxyyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 656);

    auto g_0_xxyyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 657);

    auto g_0_xxyyzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 658);

    auto g_0_xxyyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 659);

    auto g_0_xxyyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 660);

    auto g_0_xxyyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 661);

    auto g_0_xxyyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 662);

    auto g_0_xxyyzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 663);

    auto g_0_xxyyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 664);

    auto g_0_xxyyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 665);

    auto g_0_xxyyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 666);

    auto g_0_xxyyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 667);

    auto g_0_xxyyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 668);

    auto g_0_xxyyzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 669);

    auto g_0_xxyyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 670);

    auto g_0_xxyyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 671);

    auto g_0_xxyyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 672);

    auto g_0_xxyyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 673);

    auto g_0_xxyyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 674);

    auto g_0_xxyyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 675);

    auto g_0_xxyyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 676);

    auto g_0_xxyyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 677);

    auto g_0_xxyyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 678);

    auto g_0_xxyyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 679);

    auto g_0_xxyyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 680);

    auto g_0_xxyyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 681);

    auto g_0_xxyyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 682);

    auto g_0_xxyyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 683);

    auto g_0_xxyzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 684);

    auto g_0_xxyzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 686);

    auto g_0_xxyzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 689);

    auto g_0_xxyzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 693);

    auto g_0_xxyzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 698);

    auto g_0_xxyzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 704);

    auto g_0_xxyzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 711);

    auto g_0_xxzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 720);

    auto g_0_xxzzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 721);

    auto g_0_xxzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 722);

    auto g_0_xxzzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 723);

    auto g_0_xxzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 724);

    auto g_0_xxzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 725);

    auto g_0_xxzzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 726);

    auto g_0_xxzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 727);

    auto g_0_xxzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 728);

    auto g_0_xxzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 729);

    auto g_0_xxzzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 730);

    auto g_0_xxzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 731);

    auto g_0_xxzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 732);

    auto g_0_xxzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 733);

    auto g_0_xxzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 734);

    auto g_0_xxzzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 735);

    auto g_0_xxzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 736);

    auto g_0_xxzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 737);

    auto g_0_xxzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 738);

    auto g_0_xxzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 739);

    auto g_0_xxzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 740);

    auto g_0_xxzzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 741);

    auto g_0_xxzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 742);

    auto g_0_xxzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 743);

    auto g_0_xxzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 744);

    auto g_0_xxzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 745);

    auto g_0_xxzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 746);

    auto g_0_xxzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 747);

    auto g_0_xxzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 748);

    auto g_0_xxzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 749);

    auto g_0_xxzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 750);

    auto g_0_xxzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 751);

    auto g_0_xxzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 752);

    auto g_0_xxzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 753);

    auto g_0_xxzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 754);

    auto g_0_xxzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 755);

    auto g_0_xyyyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 756);

    auto g_0_xyyyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 757);

    auto g_0_xyyyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 759);

    auto g_0_xyyyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 760);

    auto g_0_xyyyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 762);

    auto g_0_xyyyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 763);

    auto g_0_xyyyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 764);

    auto g_0_xyyyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 766);

    auto g_0_xyyyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 767);

    auto g_0_xyyyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 768);

    auto g_0_xyyyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 769);

    auto g_0_xyyyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 771);

    auto g_0_xyyyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 772);

    auto g_0_xyyyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 773);

    auto g_0_xyyyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 774);

    auto g_0_xyyyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 775);

    auto g_0_xyyyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 777);

    auto g_0_xyyyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 778);

    auto g_0_xyyyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 779);

    auto g_0_xyyyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 780);

    auto g_0_xyyyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 781);

    auto g_0_xyyyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 782);

    auto g_0_xyyyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 784);

    auto g_0_xyyyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 785);

    auto g_0_xyyyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 786);

    auto g_0_xyyyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 787);

    auto g_0_xyyyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 788);

    auto g_0_xyyyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 789);

    auto g_0_xyyyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 790);

    auto g_0_xyyyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 791);

    auto g_0_xyyyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 832);

    auto g_0_xyyyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 835);

    auto g_0_xyyyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 836);

    auto g_0_xyyyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 839);

    auto g_0_xyyyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 840);

    auto g_0_xyyyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 841);

    auto g_0_xyyyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 844);

    auto g_0_xyyyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 845);

    auto g_0_xyyyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 846);

    auto g_0_xyyyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 847);

    auto g_0_xyyyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 850);

    auto g_0_xyyyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 851);

    auto g_0_xyyyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 852);

    auto g_0_xyyyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 853);

    auto g_0_xyyyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 854);

    auto g_0_xyyyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 856);

    auto g_0_xyyyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 857);

    auto g_0_xyyyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 858);

    auto g_0_xyyyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 859);

    auto g_0_xyyyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 860);

    auto g_0_xyyyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 861);

    auto g_0_xyyyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 862);

    auto g_0_xyyyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 863);

    auto g_0_xyyyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 868);

    auto g_0_xyyyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 871);

    auto g_0_xyyyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 872);

    auto g_0_xyyyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 875);

    auto g_0_xyyyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 876);

    auto g_0_xyyyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 877);

    auto g_0_xyyyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 880);

    auto g_0_xyyyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 881);

    auto g_0_xyyyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 882);

    auto g_0_xyyyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 883);

    auto g_0_xyyyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 886);

    auto g_0_xyyyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 887);

    auto g_0_xyyyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 888);

    auto g_0_xyyyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 889);

    auto g_0_xyyyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 890);

    auto g_0_xyyyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 892);

    auto g_0_xyyyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 893);

    auto g_0_xyyyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 894);

    auto g_0_xyyyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 895);

    auto g_0_xyyyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 896);

    auto g_0_xyyyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 897);

    auto g_0_xyyyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 898);

    auto g_0_xyyyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 899);

    auto g_0_xyyzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 904);

    auto g_0_xyyzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 907);

    auto g_0_xyyzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 908);

    auto g_0_xyyzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 911);

    auto g_0_xyyzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 912);

    auto g_0_xyyzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 913);

    auto g_0_xyyzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 916);

    auto g_0_xyyzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 917);

    auto g_0_xyyzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 918);

    auto g_0_xyyzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 919);

    auto g_0_xyyzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 922);

    auto g_0_xyyzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 923);

    auto g_0_xyyzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 924);

    auto g_0_xyyzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 925);

    auto g_0_xyyzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 926);

    auto g_0_xyyzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 928);

    auto g_0_xyyzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 929);

    auto g_0_xyyzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 930);

    auto g_0_xyyzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 931);

    auto g_0_xyyzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 932);

    auto g_0_xyyzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 933);

    auto g_0_xyyzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 934);

    auto g_0_xyyzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 935);

    auto g_0_xzzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 972);

    auto g_0_xzzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 974);

    auto g_0_xzzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 976);

    auto g_0_xzzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 977);

    auto g_0_xzzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 979);

    auto g_0_xzzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 980);

    auto g_0_xzzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 981);

    auto g_0_xzzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 983);

    auto g_0_xzzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 984);

    auto g_0_xzzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 985);

    auto g_0_xzzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 986);

    auto g_0_xzzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 988);

    auto g_0_xzzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 989);

    auto g_0_xzzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 990);

    auto g_0_xzzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 991);

    auto g_0_xzzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 992);

    auto g_0_xzzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 994);

    auto g_0_xzzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 995);

    auto g_0_xzzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 996);

    auto g_0_xzzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 997);

    auto g_0_xzzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 998);

    auto g_0_xzzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 999);

    auto g_0_xzzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1000);

    auto g_0_xzzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1001);

    auto g_0_xzzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1002);

    auto g_0_xzzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1003);

    auto g_0_xzzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1004);

    auto g_0_xzzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1005);

    auto g_0_xzzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1006);

    auto g_0_xzzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1007);

    auto g_0_yyyyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 1008);

    auto g_0_yyyyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1009);

    auto g_0_yyyyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 1010);

    auto g_0_yyyyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 1011);

    auto g_0_yyyyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 1012);

    auto g_0_yyyyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 1013);

    auto g_0_yyyyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 1014);

    auto g_0_yyyyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 1015);

    auto g_0_yyyyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 1016);

    auto g_0_yyyyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 1017);

    auto g_0_yyyyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1018);

    auto g_0_yyyyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1019);

    auto g_0_yyyyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1020);

    auto g_0_yyyyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1021);

    auto g_0_yyyyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1022);

    auto g_0_yyyyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1023);

    auto g_0_yyyyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1024);

    auto g_0_yyyyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1025);

    auto g_0_yyyyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1026);

    auto g_0_yyyyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1027);

    auto g_0_yyyyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1028);

    auto g_0_yyyyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1029);

    auto g_0_yyyyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1030);

    auto g_0_yyyyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1031);

    auto g_0_yyyyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1032);

    auto g_0_yyyyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1033);

    auto g_0_yyyyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1034);

    auto g_0_yyyyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1035);

    auto g_0_yyyyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1036);

    auto g_0_yyyyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1037);

    auto g_0_yyyyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1038);

    auto g_0_yyyyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1039);

    auto g_0_yyyyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1040);

    auto g_0_yyyyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1041);

    auto g_0_yyyyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1042);

    auto g_0_yyyyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1043);

    auto g_0_yyyyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1045);

    auto g_0_yyyyyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 1046);

    auto g_0_yyyyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 1047);

    auto g_0_yyyyyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 1048);

    auto g_0_yyyyyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 1049);

    auto g_0_yyyyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 1050);

    auto g_0_yyyyyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 1051);

    auto g_0_yyyyyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 1052);

    auto g_0_yyyyyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 1053);

    auto g_0_yyyyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1054);

    auto g_0_yyyyyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1055);

    auto g_0_yyyyyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1056);

    auto g_0_yyyyyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1057);

    auto g_0_yyyyyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1058);

    auto g_0_yyyyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1059);

    auto g_0_yyyyyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1060);

    auto g_0_yyyyyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1061);

    auto g_0_yyyyyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1062);

    auto g_0_yyyyyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1063);

    auto g_0_yyyyyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1064);

    auto g_0_yyyyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1065);

    auto g_0_yyyyyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1066);

    auto g_0_yyyyyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1067);

    auto g_0_yyyyyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1068);

    auto g_0_yyyyyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1069);

    auto g_0_yyyyyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1070);

    auto g_0_yyyyyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1071);

    auto g_0_yyyyyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1072);

    auto g_0_yyyyyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1073);

    auto g_0_yyyyyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1074);

    auto g_0_yyyyyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1075);

    auto g_0_yyyyyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1076);

    auto g_0_yyyyyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1077);

    auto g_0_yyyyyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1078);

    auto g_0_yyyyyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1079);

    auto g_0_yyyyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 1080);

    auto g_0_yyyyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1081);

    auto g_0_yyyyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 1082);

    auto g_0_yyyyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 1083);

    auto g_0_yyyyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 1084);

    auto g_0_yyyyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 1085);

    auto g_0_yyyyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 1086);

    auto g_0_yyyyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 1087);

    auto g_0_yyyyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 1088);

    auto g_0_yyyyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 1089);

    auto g_0_yyyyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1090);

    auto g_0_yyyyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1091);

    auto g_0_yyyyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1092);

    auto g_0_yyyyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1093);

    auto g_0_yyyyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1094);

    auto g_0_yyyyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1095);

    auto g_0_yyyyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1096);

    auto g_0_yyyyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1097);

    auto g_0_yyyyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1098);

    auto g_0_yyyyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1099);

    auto g_0_yyyyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1100);

    auto g_0_yyyyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1101);

    auto g_0_yyyyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1102);

    auto g_0_yyyyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1103);

    auto g_0_yyyyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1104);

    auto g_0_yyyyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1105);

    auto g_0_yyyyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1106);

    auto g_0_yyyyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1107);

    auto g_0_yyyyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1108);

    auto g_0_yyyyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1109);

    auto g_0_yyyyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1110);

    auto g_0_yyyyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1111);

    auto g_0_yyyyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1112);

    auto g_0_yyyyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1113);

    auto g_0_yyyyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1114);

    auto g_0_yyyyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1115);

    auto g_0_yyyyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 1116);

    auto g_0_yyyyzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1117);

    auto g_0_yyyyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 1118);

    auto g_0_yyyyzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 1119);

    auto g_0_yyyyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 1120);

    auto g_0_yyyyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 1121);

    auto g_0_yyyyzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 1122);

    auto g_0_yyyyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 1123);

    auto g_0_yyyyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 1124);

    auto g_0_yyyyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 1125);

    auto g_0_yyyyzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1126);

    auto g_0_yyyyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1127);

    auto g_0_yyyyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1128);

    auto g_0_yyyyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1129);

    auto g_0_yyyyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1130);

    auto g_0_yyyyzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1131);

    auto g_0_yyyyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1132);

    auto g_0_yyyyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1133);

    auto g_0_yyyyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1134);

    auto g_0_yyyyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1135);

    auto g_0_yyyyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1136);

    auto g_0_yyyyzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1137);

    auto g_0_yyyyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1138);

    auto g_0_yyyyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1139);

    auto g_0_yyyyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1140);

    auto g_0_yyyyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1141);

    auto g_0_yyyyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1142);

    auto g_0_yyyyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1143);

    auto g_0_yyyyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1144);

    auto g_0_yyyyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1145);

    auto g_0_yyyyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1146);

    auto g_0_yyyyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1147);

    auto g_0_yyyyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1148);

    auto g_0_yyyyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1149);

    auto g_0_yyyyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1150);

    auto g_0_yyyyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1151);

    auto g_0_yyyzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 1152);

    auto g_0_yyyzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1153);

    auto g_0_yyyzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 1154);

    auto g_0_yyyzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 1155);

    auto g_0_yyyzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 1156);

    auto g_0_yyyzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 1157);

    auto g_0_yyyzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 1158);

    auto g_0_yyyzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 1159);

    auto g_0_yyyzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 1160);

    auto g_0_yyyzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 1161);

    auto g_0_yyyzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1162);

    auto g_0_yyyzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1163);

    auto g_0_yyyzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1164);

    auto g_0_yyyzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1165);

    auto g_0_yyyzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1166);

    auto g_0_yyyzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1167);

    auto g_0_yyyzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1168);

    auto g_0_yyyzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1169);

    auto g_0_yyyzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1170);

    auto g_0_yyyzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1171);

    auto g_0_yyyzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1172);

    auto g_0_yyyzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1173);

    auto g_0_yyyzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1174);

    auto g_0_yyyzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1175);

    auto g_0_yyyzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1176);

    auto g_0_yyyzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1177);

    auto g_0_yyyzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1178);

    auto g_0_yyyzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1179);

    auto g_0_yyyzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1180);

    auto g_0_yyyzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1181);

    auto g_0_yyyzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1182);

    auto g_0_yyyzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1183);

    auto g_0_yyyzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1184);

    auto g_0_yyyzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1185);

    auto g_0_yyyzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1186);

    auto g_0_yyyzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1187);

    auto g_0_yyzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 1188);

    auto g_0_yyzzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1189);

    auto g_0_yyzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 1190);

    auto g_0_yyzzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 1191);

    auto g_0_yyzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 1192);

    auto g_0_yyzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 1193);

    auto g_0_yyzzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 1194);

    auto g_0_yyzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 1195);

    auto g_0_yyzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 1196);

    auto g_0_yyzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 1197);

    auto g_0_yyzzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1198);

    auto g_0_yyzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1199);

    auto g_0_yyzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1200);

    auto g_0_yyzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1201);

    auto g_0_yyzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1202);

    auto g_0_yyzzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1203);

    auto g_0_yyzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1204);

    auto g_0_yyzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1205);

    auto g_0_yyzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1206);

    auto g_0_yyzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1207);

    auto g_0_yyzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1208);

    auto g_0_yyzzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1209);

    auto g_0_yyzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1210);

    auto g_0_yyzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1211);

    auto g_0_yyzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1212);

    auto g_0_yyzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1213);

    auto g_0_yyzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1214);

    auto g_0_yyzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1215);

    auto g_0_yyzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1216);

    auto g_0_yyzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1217);

    auto g_0_yyzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1218);

    auto g_0_yyzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1219);

    auto g_0_yyzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1220);

    auto g_0_yyzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1221);

    auto g_0_yyzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1222);

    auto g_0_yyzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1223);

    auto g_0_yzzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 1224);

    auto g_0_yzzzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1225);

    auto g_0_yzzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 1226);

    auto g_0_yzzzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 1227);

    auto g_0_yzzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 1228);

    auto g_0_yzzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 1229);

    auto g_0_yzzzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 1230);

    auto g_0_yzzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 1231);

    auto g_0_yzzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 1232);

    auto g_0_yzzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 1233);

    auto g_0_yzzzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1234);

    auto g_0_yzzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1235);

    auto g_0_yzzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1236);

    auto g_0_yzzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1237);

    auto g_0_yzzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1238);

    auto g_0_yzzzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1239);

    auto g_0_yzzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1240);

    auto g_0_yzzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1241);

    auto g_0_yzzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1242);

    auto g_0_yzzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1243);

    auto g_0_yzzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1244);

    auto g_0_yzzzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1245);

    auto g_0_yzzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1246);

    auto g_0_yzzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1247);

    auto g_0_yzzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1248);

    auto g_0_yzzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1249);

    auto g_0_yzzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1250);

    auto g_0_yzzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1251);

    auto g_0_yzzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1252);

    auto g_0_yzzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1253);

    auto g_0_yzzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1254);

    auto g_0_yzzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1255);

    auto g_0_yzzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1256);

    auto g_0_yzzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1257);

    auto g_0_yzzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1258);

    auto g_0_yzzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1259);

    auto g_0_zzzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sksk + 1260);

    auto g_0_zzzzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sksk + 1261);

    auto g_0_zzzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sksk + 1262);

    auto g_0_zzzzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sksk + 1263);

    auto g_0_zzzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sksk + 1264);

    auto g_0_zzzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sksk + 1265);

    auto g_0_zzzzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sksk + 1266);

    auto g_0_zzzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sksk + 1267);

    auto g_0_zzzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sksk + 1268);

    auto g_0_zzzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sksk + 1269);

    auto g_0_zzzzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1270);

    auto g_0_zzzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1271);

    auto g_0_zzzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1272);

    auto g_0_zzzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1273);

    auto g_0_zzzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1274);

    auto g_0_zzzzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1275);

    auto g_0_zzzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1276);

    auto g_0_zzzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1277);

    auto g_0_zzzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1278);

    auto g_0_zzzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1279);

    auto g_0_zzzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1280);

    auto g_0_zzzzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1281);

    auto g_0_zzzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1282);

    auto g_0_zzzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1283);

    auto g_0_zzzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1284);

    auto g_0_zzzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1285);

    auto g_0_zzzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1286);

    auto g_0_zzzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1287);

    auto g_0_zzzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sksk + 1288);

    auto g_0_zzzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sksk + 1289);

    auto g_0_zzzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sksk + 1290);

    auto g_0_zzzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sksk + 1291);

    auto g_0_zzzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1292);

    auto g_0_zzzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1293);

    auto g_0_zzzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1294);

    auto g_0_zzzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sksk + 1295);

    /// Set up 0-36 components of targeted buffer : SLSK

    auto g_0_xxxxxxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk);

    auto g_0_xxxxxxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1);

    auto g_0_xxxxxxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 2);

    auto g_0_xxxxxxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 3);

    auto g_0_xxxxxxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 4);

    auto g_0_xxxxxxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 5);

    auto g_0_xxxxxxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 6);

    auto g_0_xxxxxxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 7);

    auto g_0_xxxxxxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 8);

    auto g_0_xxxxxxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 9);

    auto g_0_xxxxxxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 10);

    auto g_0_xxxxxxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 11);

    auto g_0_xxxxxxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 12);

    auto g_0_xxxxxxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 13);

    auto g_0_xxxxxxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 14);

    auto g_0_xxxxxxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 15);

    auto g_0_xxxxxxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 16);

    auto g_0_xxxxxxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 17);

    auto g_0_xxxxxxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 18);

    auto g_0_xxxxxxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 19);

    auto g_0_xxxxxxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 20);

    auto g_0_xxxxxxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 21);

    auto g_0_xxxxxxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 22);

    auto g_0_xxxxxxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 23);

    auto g_0_xxxxxxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 24);

    auto g_0_xxxxxxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 25);

    auto g_0_xxxxxxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 26);

    auto g_0_xxxxxxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 27);

    auto g_0_xxxxxxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 28);

    auto g_0_xxxxxxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 29);

    auto g_0_xxxxxxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 30);

    auto g_0_xxxxxxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 31);

    auto g_0_xxxxxxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 32);

    auto g_0_xxxxxxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 33);

    auto g_0_xxxxxxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 34);

    auto g_0_xxxxxxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 35);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxxxxx_0, g_0_xxxxxx_0_xxxxxxx_1, g_0_xxxxxx_0_xxxxxxy_0, g_0_xxxxxx_0_xxxxxxy_1, g_0_xxxxxx_0_xxxxxxz_0, g_0_xxxxxx_0_xxxxxxz_1, g_0_xxxxxx_0_xxxxxyy_0, g_0_xxxxxx_0_xxxxxyy_1, g_0_xxxxxx_0_xxxxxyz_0, g_0_xxxxxx_0_xxxxxyz_1, g_0_xxxxxx_0_xxxxxzz_0, g_0_xxxxxx_0_xxxxxzz_1, g_0_xxxxxx_0_xxxxyyy_0, g_0_xxxxxx_0_xxxxyyy_1, g_0_xxxxxx_0_xxxxyyz_0, g_0_xxxxxx_0_xxxxyyz_1, g_0_xxxxxx_0_xxxxyzz_0, g_0_xxxxxx_0_xxxxyzz_1, g_0_xxxxxx_0_xxxxzzz_0, g_0_xxxxxx_0_xxxxzzz_1, g_0_xxxxxx_0_xxxyyyy_0, g_0_xxxxxx_0_xxxyyyy_1, g_0_xxxxxx_0_xxxyyyz_0, g_0_xxxxxx_0_xxxyyyz_1, g_0_xxxxxx_0_xxxyyzz_0, g_0_xxxxxx_0_xxxyyzz_1, g_0_xxxxxx_0_xxxyzzz_0, g_0_xxxxxx_0_xxxyzzz_1, g_0_xxxxxx_0_xxxzzzz_0, g_0_xxxxxx_0_xxxzzzz_1, g_0_xxxxxx_0_xxyyyyy_0, g_0_xxxxxx_0_xxyyyyy_1, g_0_xxxxxx_0_xxyyyyz_0, g_0_xxxxxx_0_xxyyyyz_1, g_0_xxxxxx_0_xxyyyzz_0, g_0_xxxxxx_0_xxyyyzz_1, g_0_xxxxxx_0_xxyyzzz_0, g_0_xxxxxx_0_xxyyzzz_1, g_0_xxxxxx_0_xxyzzzz_0, g_0_xxxxxx_0_xxyzzzz_1, g_0_xxxxxx_0_xxzzzzz_0, g_0_xxxxxx_0_xxzzzzz_1, g_0_xxxxxx_0_xyyyyyy_0, g_0_xxxxxx_0_xyyyyyy_1, g_0_xxxxxx_0_xyyyyyz_0, g_0_xxxxxx_0_xyyyyyz_1, g_0_xxxxxx_0_xyyyyzz_0, g_0_xxxxxx_0_xyyyyzz_1, g_0_xxxxxx_0_xyyyzzz_0, g_0_xxxxxx_0_xyyyzzz_1, g_0_xxxxxx_0_xyyzzzz_0, g_0_xxxxxx_0_xyyzzzz_1, g_0_xxxxxx_0_xyzzzzz_0, g_0_xxxxxx_0_xyzzzzz_1, g_0_xxxxxx_0_xzzzzzz_0, g_0_xxxxxx_0_xzzzzzz_1, g_0_xxxxxx_0_yyyyyyy_0, g_0_xxxxxx_0_yyyyyyy_1, g_0_xxxxxx_0_yyyyyyz_0, g_0_xxxxxx_0_yyyyyyz_1, g_0_xxxxxx_0_yyyyyzz_0, g_0_xxxxxx_0_yyyyyzz_1, g_0_xxxxxx_0_yyyyzzz_0, g_0_xxxxxx_0_yyyyzzz_1, g_0_xxxxxx_0_yyyzzzz_0, g_0_xxxxxx_0_yyyzzzz_1, g_0_xxxxxx_0_yyzzzzz_0, g_0_xxxxxx_0_yyzzzzz_1, g_0_xxxxxx_0_yzzzzzz_0, g_0_xxxxxx_0_yzzzzzz_1, g_0_xxxxxx_0_zzzzzzz_0, g_0_xxxxxx_0_zzzzzzz_1, g_0_xxxxxxx_0_xxxxxx_1, g_0_xxxxxxx_0_xxxxxxx_0, g_0_xxxxxxx_0_xxxxxxx_1, g_0_xxxxxxx_0_xxxxxxy_0, g_0_xxxxxxx_0_xxxxxxy_1, g_0_xxxxxxx_0_xxxxxxz_0, g_0_xxxxxxx_0_xxxxxxz_1, g_0_xxxxxxx_0_xxxxxy_1, g_0_xxxxxxx_0_xxxxxyy_0, g_0_xxxxxxx_0_xxxxxyy_1, g_0_xxxxxxx_0_xxxxxyz_0, g_0_xxxxxxx_0_xxxxxyz_1, g_0_xxxxxxx_0_xxxxxz_1, g_0_xxxxxxx_0_xxxxxzz_0, g_0_xxxxxxx_0_xxxxxzz_1, g_0_xxxxxxx_0_xxxxyy_1, g_0_xxxxxxx_0_xxxxyyy_0, g_0_xxxxxxx_0_xxxxyyy_1, g_0_xxxxxxx_0_xxxxyyz_0, g_0_xxxxxxx_0_xxxxyyz_1, g_0_xxxxxxx_0_xxxxyz_1, g_0_xxxxxxx_0_xxxxyzz_0, g_0_xxxxxxx_0_xxxxyzz_1, g_0_xxxxxxx_0_xxxxzz_1, g_0_xxxxxxx_0_xxxxzzz_0, g_0_xxxxxxx_0_xxxxzzz_1, g_0_xxxxxxx_0_xxxyyy_1, g_0_xxxxxxx_0_xxxyyyy_0, g_0_xxxxxxx_0_xxxyyyy_1, g_0_xxxxxxx_0_xxxyyyz_0, g_0_xxxxxxx_0_xxxyyyz_1, g_0_xxxxxxx_0_xxxyyz_1, g_0_xxxxxxx_0_xxxyyzz_0, g_0_xxxxxxx_0_xxxyyzz_1, g_0_xxxxxxx_0_xxxyzz_1, g_0_xxxxxxx_0_xxxyzzz_0, g_0_xxxxxxx_0_xxxyzzz_1, g_0_xxxxxxx_0_xxxzzz_1, g_0_xxxxxxx_0_xxxzzzz_0, g_0_xxxxxxx_0_xxxzzzz_1, g_0_xxxxxxx_0_xxyyyy_1, g_0_xxxxxxx_0_xxyyyyy_0, g_0_xxxxxxx_0_xxyyyyy_1, g_0_xxxxxxx_0_xxyyyyz_0, g_0_xxxxxxx_0_xxyyyyz_1, g_0_xxxxxxx_0_xxyyyz_1, g_0_xxxxxxx_0_xxyyyzz_0, g_0_xxxxxxx_0_xxyyyzz_1, g_0_xxxxxxx_0_xxyyzz_1, g_0_xxxxxxx_0_xxyyzzz_0, g_0_xxxxxxx_0_xxyyzzz_1, g_0_xxxxxxx_0_xxyzzz_1, g_0_xxxxxxx_0_xxyzzzz_0, g_0_xxxxxxx_0_xxyzzzz_1, g_0_xxxxxxx_0_xxzzzz_1, g_0_xxxxxxx_0_xxzzzzz_0, g_0_xxxxxxx_0_xxzzzzz_1, g_0_xxxxxxx_0_xyyyyy_1, g_0_xxxxxxx_0_xyyyyyy_0, g_0_xxxxxxx_0_xyyyyyy_1, g_0_xxxxxxx_0_xyyyyyz_0, g_0_xxxxxxx_0_xyyyyyz_1, g_0_xxxxxxx_0_xyyyyz_1, g_0_xxxxxxx_0_xyyyyzz_0, g_0_xxxxxxx_0_xyyyyzz_1, g_0_xxxxxxx_0_xyyyzz_1, g_0_xxxxxxx_0_xyyyzzz_0, g_0_xxxxxxx_0_xyyyzzz_1, g_0_xxxxxxx_0_xyyzzz_1, g_0_xxxxxxx_0_xyyzzzz_0, g_0_xxxxxxx_0_xyyzzzz_1, g_0_xxxxxxx_0_xyzzzz_1, g_0_xxxxxxx_0_xyzzzzz_0, g_0_xxxxxxx_0_xyzzzzz_1, g_0_xxxxxxx_0_xzzzzz_1, g_0_xxxxxxx_0_xzzzzzz_0, g_0_xxxxxxx_0_xzzzzzz_1, g_0_xxxxxxx_0_yyyyyy_1, g_0_xxxxxxx_0_yyyyyyy_0, g_0_xxxxxxx_0_yyyyyyy_1, g_0_xxxxxxx_0_yyyyyyz_0, g_0_xxxxxxx_0_yyyyyyz_1, g_0_xxxxxxx_0_yyyyyz_1, g_0_xxxxxxx_0_yyyyyzz_0, g_0_xxxxxxx_0_yyyyyzz_1, g_0_xxxxxxx_0_yyyyzz_1, g_0_xxxxxxx_0_yyyyzzz_0, g_0_xxxxxxx_0_yyyyzzz_1, g_0_xxxxxxx_0_yyyzzz_1, g_0_xxxxxxx_0_yyyzzzz_0, g_0_xxxxxxx_0_yyyzzzz_1, g_0_xxxxxxx_0_yyzzzz_1, g_0_xxxxxxx_0_yyzzzzz_0, g_0_xxxxxxx_0_yyzzzzz_1, g_0_xxxxxxx_0_yzzzzz_1, g_0_xxxxxxx_0_yzzzzzz_0, g_0_xxxxxxx_0_yzzzzzz_1, g_0_xxxxxxx_0_zzzzzz_1, g_0_xxxxxxx_0_zzzzzzz_0, g_0_xxxxxxx_0_zzzzzzz_1, g_0_xxxxxxxx_0_xxxxxxx_0, g_0_xxxxxxxx_0_xxxxxxy_0, g_0_xxxxxxxx_0_xxxxxxz_0, g_0_xxxxxxxx_0_xxxxxyy_0, g_0_xxxxxxxx_0_xxxxxyz_0, g_0_xxxxxxxx_0_xxxxxzz_0, g_0_xxxxxxxx_0_xxxxyyy_0, g_0_xxxxxxxx_0_xxxxyyz_0, g_0_xxxxxxxx_0_xxxxyzz_0, g_0_xxxxxxxx_0_xxxxzzz_0, g_0_xxxxxxxx_0_xxxyyyy_0, g_0_xxxxxxxx_0_xxxyyyz_0, g_0_xxxxxxxx_0_xxxyyzz_0, g_0_xxxxxxxx_0_xxxyzzz_0, g_0_xxxxxxxx_0_xxxzzzz_0, g_0_xxxxxxxx_0_xxyyyyy_0, g_0_xxxxxxxx_0_xxyyyyz_0, g_0_xxxxxxxx_0_xxyyyzz_0, g_0_xxxxxxxx_0_xxyyzzz_0, g_0_xxxxxxxx_0_xxyzzzz_0, g_0_xxxxxxxx_0_xxzzzzz_0, g_0_xxxxxxxx_0_xyyyyyy_0, g_0_xxxxxxxx_0_xyyyyyz_0, g_0_xxxxxxxx_0_xyyyyzz_0, g_0_xxxxxxxx_0_xyyyzzz_0, g_0_xxxxxxxx_0_xyyzzzz_0, g_0_xxxxxxxx_0_xyzzzzz_0, g_0_xxxxxxxx_0_xzzzzzz_0, g_0_xxxxxxxx_0_yyyyyyy_0, g_0_xxxxxxxx_0_yyyyyyz_0, g_0_xxxxxxxx_0_yyyyyzz_0, g_0_xxxxxxxx_0_yyyyzzz_0, g_0_xxxxxxxx_0_yyyzzzz_0, g_0_xxxxxxxx_0_yyzzzzz_0, g_0_xxxxxxxx_0_yzzzzzz_0, g_0_xxxxxxxx_0_zzzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxxx_0_xxxxxxx_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxxx_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxxx_1[i] * fti_ab_0 + 7.0 * g_0_xxxxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxxx_0[i] * pb_x + g_0_xxxxxxx_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxxxy_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxxy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxxy_0[i] * pb_x + g_0_xxxxxxx_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxxxz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxxz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxxz_0[i] * pb_x + g_0_xxxxxxx_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxxyy_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxyy_0[i] * pb_x + g_0_xxxxxxx_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxxyz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxyz_0[i] * pb_x + g_0_xxxxxxx_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxxzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxzz_0[i] * pb_x + g_0_xxxxxxx_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxyyy_0[i] = 7.0 * g_0_xxxxxx_0_xxxxyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyyy_0[i] * pb_x + g_0_xxxxxxx_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxyyz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyyz_0[i] * pb_x + g_0_xxxxxxx_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxyzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyzz_0[i] * pb_x + g_0_xxxxxxx_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxzzz_0[i] * pb_x + g_0_xxxxxxx_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyyyy_0[i] = 7.0 * g_0_xxxxxx_0_xxxyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyyy_0[i] * pb_x + g_0_xxxxxxx_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyyyz_0[i] = 7.0 * g_0_xxxxxx_0_xxxyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyyz_0[i] * pb_x + g_0_xxxxxxx_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyyzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyzz_0[i] * pb_x + g_0_xxxxxxx_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyzzz_0[i] * pb_x + g_0_xxxxxxx_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyyyy_0[i] = 7.0 * g_0_xxxxxx_0_xxyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyyy_0[i] * pb_x + g_0_xxxxxxx_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyyyz_0[i] = 7.0 * g_0_xxxxxx_0_xxyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyyz_0[i] * pb_x + g_0_xxxxxxx_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyyzz_0[i] = 7.0 * g_0_xxxxxx_0_xxyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyzz_0[i] * pb_x + g_0_xxxxxxx_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyzzz_0[i] * pb_x + g_0_xxxxxxx_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyyyy_0[i] = 7.0 * g_0_xxxxxx_0_xyyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyyy_0[i] * pb_x + g_0_xxxxxxx_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyyyz_0[i] = 7.0 * g_0_xxxxxx_0_xyyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyyz_0[i] * pb_x + g_0_xxxxxxx_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyyzz_0[i] = 7.0 * g_0_xxxxxx_0_xyyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyzz_0[i] * pb_x + g_0_xxxxxxx_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyzzz_0[i] = 7.0 * g_0_xxxxxx_0_xyyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyzzz_0[i] * pb_x + g_0_xxxxxxx_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xyyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xyzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xzzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyyyy_0[i] = 7.0 * g_0_xxxxxx_0_yyyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyyyy_0[i] * pb_x + g_0_xxxxxxx_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyyyz_0[i] = 7.0 * g_0_xxxxxx_0_yyyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyyyz_0[i] * pb_x + g_0_xxxxxxx_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyyzz_0[i] = 7.0 * g_0_xxxxxx_0_yyyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyyzz_0[i] * pb_x + g_0_xxxxxxx_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyzzz_0[i] = 7.0 * g_0_xxxxxx_0_yyyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyzzz_0[i] * pb_x + g_0_xxxxxxx_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyzzzz_0[i] = 7.0 * g_0_xxxxxx_0_yyyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyzzzz_0[i] * pb_x + g_0_xxxxxxx_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_yyzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yzzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_yzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yzzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_zzzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_zzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zzzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 36-72 components of targeted buffer : SLSK

    auto g_0_xxxxxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 36);

    auto g_0_xxxxxxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 37);

    auto g_0_xxxxxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 38);

    auto g_0_xxxxxxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 39);

    auto g_0_xxxxxxxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 40);

    auto g_0_xxxxxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 41);

    auto g_0_xxxxxxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 42);

    auto g_0_xxxxxxxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 43);

    auto g_0_xxxxxxxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 44);

    auto g_0_xxxxxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 45);

    auto g_0_xxxxxxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 46);

    auto g_0_xxxxxxxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 47);

    auto g_0_xxxxxxxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 48);

    auto g_0_xxxxxxxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 49);

    auto g_0_xxxxxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 50);

    auto g_0_xxxxxxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 51);

    auto g_0_xxxxxxxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 52);

    auto g_0_xxxxxxxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 53);

    auto g_0_xxxxxxxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 54);

    auto g_0_xxxxxxxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 55);

    auto g_0_xxxxxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 56);

    auto g_0_xxxxxxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 57);

    auto g_0_xxxxxxxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 58);

    auto g_0_xxxxxxxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 59);

    auto g_0_xxxxxxxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 60);

    auto g_0_xxxxxxxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 61);

    auto g_0_xxxxxxxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 62);

    auto g_0_xxxxxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 63);

    auto g_0_xxxxxxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 64);

    auto g_0_xxxxxxxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 65);

    auto g_0_xxxxxxxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 66);

    auto g_0_xxxxxxxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 67);

    auto g_0_xxxxxxxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 68);

    auto g_0_xxxxxxxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 69);

    auto g_0_xxxxxxxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 70);

    auto g_0_xxxxxxxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 71);

    #pragma omp simd aligned(g_0_xxxxxxx_0_xxxxxx_1, g_0_xxxxxxx_0_xxxxxxx_0, g_0_xxxxxxx_0_xxxxxxx_1, g_0_xxxxxxx_0_xxxxxxy_0, g_0_xxxxxxx_0_xxxxxxy_1, g_0_xxxxxxx_0_xxxxxxz_0, g_0_xxxxxxx_0_xxxxxxz_1, g_0_xxxxxxx_0_xxxxxy_1, g_0_xxxxxxx_0_xxxxxyy_0, g_0_xxxxxxx_0_xxxxxyy_1, g_0_xxxxxxx_0_xxxxxyz_0, g_0_xxxxxxx_0_xxxxxyz_1, g_0_xxxxxxx_0_xxxxxz_1, g_0_xxxxxxx_0_xxxxxzz_0, g_0_xxxxxxx_0_xxxxxzz_1, g_0_xxxxxxx_0_xxxxyy_1, g_0_xxxxxxx_0_xxxxyyy_0, g_0_xxxxxxx_0_xxxxyyy_1, g_0_xxxxxxx_0_xxxxyyz_0, g_0_xxxxxxx_0_xxxxyyz_1, g_0_xxxxxxx_0_xxxxyz_1, g_0_xxxxxxx_0_xxxxyzz_0, g_0_xxxxxxx_0_xxxxyzz_1, g_0_xxxxxxx_0_xxxxzz_1, g_0_xxxxxxx_0_xxxxzzz_0, g_0_xxxxxxx_0_xxxxzzz_1, g_0_xxxxxxx_0_xxxyyy_1, g_0_xxxxxxx_0_xxxyyyy_0, g_0_xxxxxxx_0_xxxyyyy_1, g_0_xxxxxxx_0_xxxyyyz_0, g_0_xxxxxxx_0_xxxyyyz_1, g_0_xxxxxxx_0_xxxyyz_1, g_0_xxxxxxx_0_xxxyyzz_0, g_0_xxxxxxx_0_xxxyyzz_1, g_0_xxxxxxx_0_xxxyzz_1, g_0_xxxxxxx_0_xxxyzzz_0, g_0_xxxxxxx_0_xxxyzzz_1, g_0_xxxxxxx_0_xxxzzz_1, g_0_xxxxxxx_0_xxxzzzz_0, g_0_xxxxxxx_0_xxxzzzz_1, g_0_xxxxxxx_0_xxyyyy_1, g_0_xxxxxxx_0_xxyyyyy_0, g_0_xxxxxxx_0_xxyyyyy_1, g_0_xxxxxxx_0_xxyyyyz_0, g_0_xxxxxxx_0_xxyyyyz_1, g_0_xxxxxxx_0_xxyyyz_1, g_0_xxxxxxx_0_xxyyyzz_0, g_0_xxxxxxx_0_xxyyyzz_1, g_0_xxxxxxx_0_xxyyzz_1, g_0_xxxxxxx_0_xxyyzzz_0, g_0_xxxxxxx_0_xxyyzzz_1, g_0_xxxxxxx_0_xxyzzz_1, g_0_xxxxxxx_0_xxyzzzz_0, g_0_xxxxxxx_0_xxyzzzz_1, g_0_xxxxxxx_0_xxzzzz_1, g_0_xxxxxxx_0_xxzzzzz_0, g_0_xxxxxxx_0_xxzzzzz_1, g_0_xxxxxxx_0_xyyyyy_1, g_0_xxxxxxx_0_xyyyyyy_0, g_0_xxxxxxx_0_xyyyyyy_1, g_0_xxxxxxx_0_xyyyyyz_0, g_0_xxxxxxx_0_xyyyyyz_1, g_0_xxxxxxx_0_xyyyyz_1, g_0_xxxxxxx_0_xyyyyzz_0, g_0_xxxxxxx_0_xyyyyzz_1, g_0_xxxxxxx_0_xyyyzz_1, g_0_xxxxxxx_0_xyyyzzz_0, g_0_xxxxxxx_0_xyyyzzz_1, g_0_xxxxxxx_0_xyyzzz_1, g_0_xxxxxxx_0_xyyzzzz_0, g_0_xxxxxxx_0_xyyzzzz_1, g_0_xxxxxxx_0_xyzzzz_1, g_0_xxxxxxx_0_xyzzzzz_0, g_0_xxxxxxx_0_xyzzzzz_1, g_0_xxxxxxx_0_xzzzzz_1, g_0_xxxxxxx_0_xzzzzzz_0, g_0_xxxxxxx_0_xzzzzzz_1, g_0_xxxxxxx_0_yyyyyy_1, g_0_xxxxxxx_0_yyyyyyy_0, g_0_xxxxxxx_0_yyyyyyy_1, g_0_xxxxxxx_0_yyyyyyz_0, g_0_xxxxxxx_0_yyyyyyz_1, g_0_xxxxxxx_0_yyyyyz_1, g_0_xxxxxxx_0_yyyyyzz_0, g_0_xxxxxxx_0_yyyyyzz_1, g_0_xxxxxxx_0_yyyyzz_1, g_0_xxxxxxx_0_yyyyzzz_0, g_0_xxxxxxx_0_yyyyzzz_1, g_0_xxxxxxx_0_yyyzzz_1, g_0_xxxxxxx_0_yyyzzzz_0, g_0_xxxxxxx_0_yyyzzzz_1, g_0_xxxxxxx_0_yyzzzz_1, g_0_xxxxxxx_0_yyzzzzz_0, g_0_xxxxxxx_0_yyzzzzz_1, g_0_xxxxxxx_0_yzzzzz_1, g_0_xxxxxxx_0_yzzzzzz_0, g_0_xxxxxxx_0_yzzzzzz_1, g_0_xxxxxxx_0_zzzzzz_1, g_0_xxxxxxx_0_zzzzzzz_0, g_0_xxxxxxx_0_zzzzzzz_1, g_0_xxxxxxxy_0_xxxxxxx_0, g_0_xxxxxxxy_0_xxxxxxy_0, g_0_xxxxxxxy_0_xxxxxxz_0, g_0_xxxxxxxy_0_xxxxxyy_0, g_0_xxxxxxxy_0_xxxxxyz_0, g_0_xxxxxxxy_0_xxxxxzz_0, g_0_xxxxxxxy_0_xxxxyyy_0, g_0_xxxxxxxy_0_xxxxyyz_0, g_0_xxxxxxxy_0_xxxxyzz_0, g_0_xxxxxxxy_0_xxxxzzz_0, g_0_xxxxxxxy_0_xxxyyyy_0, g_0_xxxxxxxy_0_xxxyyyz_0, g_0_xxxxxxxy_0_xxxyyzz_0, g_0_xxxxxxxy_0_xxxyzzz_0, g_0_xxxxxxxy_0_xxxzzzz_0, g_0_xxxxxxxy_0_xxyyyyy_0, g_0_xxxxxxxy_0_xxyyyyz_0, g_0_xxxxxxxy_0_xxyyyzz_0, g_0_xxxxxxxy_0_xxyyzzz_0, g_0_xxxxxxxy_0_xxyzzzz_0, g_0_xxxxxxxy_0_xxzzzzz_0, g_0_xxxxxxxy_0_xyyyyyy_0, g_0_xxxxxxxy_0_xyyyyyz_0, g_0_xxxxxxxy_0_xyyyyzz_0, g_0_xxxxxxxy_0_xyyyzzz_0, g_0_xxxxxxxy_0_xyyzzzz_0, g_0_xxxxxxxy_0_xyzzzzz_0, g_0_xxxxxxxy_0_xzzzzzz_0, g_0_xxxxxxxy_0_yyyyyyy_0, g_0_xxxxxxxy_0_yyyyyyz_0, g_0_xxxxxxxy_0_yyyyyzz_0, g_0_xxxxxxxy_0_yyyyzzz_0, g_0_xxxxxxxy_0_yyyzzzz_0, g_0_xxxxxxxy_0_yyzzzzz_0, g_0_xxxxxxxy_0_yzzzzzz_0, g_0_xxxxxxxy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxy_0_xxxxxxx_0[i] = g_0_xxxxxxx_0_xxxxxxx_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxxxy_0[i] = g_0_xxxxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxxy_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxxxz_0[i] = g_0_xxxxxxx_0_xxxxxxz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxxyy_0[i] = 2.0 * g_0_xxxxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxyy_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxxyz_0[i] = g_0_xxxxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxyz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxxzz_0[i] = g_0_xxxxxxx_0_xxxxxzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxyyy_0[i] = 3.0 * g_0_xxxxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyyy_0[i] * pb_y + g_0_xxxxxxx_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyyz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxyzz_0[i] = g_0_xxxxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxzzz_0[i] = g_0_xxxxxxx_0_xxxxzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyyyy_0[i] = 4.0 * g_0_xxxxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyyy_0[i] * pb_y + g_0_xxxxxxx_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyyz_0[i] * pb_y + g_0_xxxxxxx_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyzzz_0[i] = g_0_xxxxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxzzzz_0[i] = g_0_xxxxxxx_0_xxxzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyyyy_0[i] = 5.0 * g_0_xxxxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyyy_0[i] * pb_y + g_0_xxxxxxx_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyyz_0[i] * pb_y + g_0_xxxxxxx_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyzz_0[i] * pb_y + g_0_xxxxxxx_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyzzzz_0[i] = g_0_xxxxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxzzzzz_0[i] = g_0_xxxxxxx_0_xxzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyyyy_0[i] = 6.0 * g_0_xxxxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyyy_0[i] * pb_y + g_0_xxxxxxx_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyyz_0[i] * pb_y + g_0_xxxxxxx_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyzz_0[i] * pb_y + g_0_xxxxxxx_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyzzz_0[i] * pb_y + g_0_xxxxxxx_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyzzzzz_0[i] = g_0_xxxxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xzzzzzz_0[i] = g_0_xxxxxxx_0_xzzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyyyy_0[i] = 7.0 * g_0_xxxxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyyyy_0[i] * pb_y + g_0_xxxxxxx_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyyyz_0[i] * pb_y + g_0_xxxxxxx_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyyzz_0[i] * pb_y + g_0_xxxxxxx_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyzzz_0[i] * pb_y + g_0_xxxxxxx_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyzzzz_0[i] * pb_y + g_0_xxxxxxx_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yzzzzzz_0[i] = g_0_xxxxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_zzzzzzz_0[i] = g_0_xxxxxxx_0_zzzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 72-108 components of targeted buffer : SLSK

    auto g_0_xxxxxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 72);

    auto g_0_xxxxxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 73);

    auto g_0_xxxxxxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 74);

    auto g_0_xxxxxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 75);

    auto g_0_xxxxxxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 76);

    auto g_0_xxxxxxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 77);

    auto g_0_xxxxxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 78);

    auto g_0_xxxxxxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 79);

    auto g_0_xxxxxxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 80);

    auto g_0_xxxxxxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 81);

    auto g_0_xxxxxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 82);

    auto g_0_xxxxxxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 83);

    auto g_0_xxxxxxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 84);

    auto g_0_xxxxxxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 85);

    auto g_0_xxxxxxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 86);

    auto g_0_xxxxxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 87);

    auto g_0_xxxxxxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 88);

    auto g_0_xxxxxxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 89);

    auto g_0_xxxxxxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 90);

    auto g_0_xxxxxxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 91);

    auto g_0_xxxxxxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 92);

    auto g_0_xxxxxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 93);

    auto g_0_xxxxxxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 94);

    auto g_0_xxxxxxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 95);

    auto g_0_xxxxxxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 96);

    auto g_0_xxxxxxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 97);

    auto g_0_xxxxxxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 98);

    auto g_0_xxxxxxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 99);

    auto g_0_xxxxxxxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 100);

    auto g_0_xxxxxxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 101);

    auto g_0_xxxxxxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 102);

    auto g_0_xxxxxxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 103);

    auto g_0_xxxxxxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 104);

    auto g_0_xxxxxxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 105);

    auto g_0_xxxxxxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 106);

    auto g_0_xxxxxxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 107);

    #pragma omp simd aligned(g_0_xxxxxxx_0_xxxxxx_1, g_0_xxxxxxx_0_xxxxxxx_0, g_0_xxxxxxx_0_xxxxxxx_1, g_0_xxxxxxx_0_xxxxxxy_0, g_0_xxxxxxx_0_xxxxxxy_1, g_0_xxxxxxx_0_xxxxxxz_0, g_0_xxxxxxx_0_xxxxxxz_1, g_0_xxxxxxx_0_xxxxxy_1, g_0_xxxxxxx_0_xxxxxyy_0, g_0_xxxxxxx_0_xxxxxyy_1, g_0_xxxxxxx_0_xxxxxyz_0, g_0_xxxxxxx_0_xxxxxyz_1, g_0_xxxxxxx_0_xxxxxz_1, g_0_xxxxxxx_0_xxxxxzz_0, g_0_xxxxxxx_0_xxxxxzz_1, g_0_xxxxxxx_0_xxxxyy_1, g_0_xxxxxxx_0_xxxxyyy_0, g_0_xxxxxxx_0_xxxxyyy_1, g_0_xxxxxxx_0_xxxxyyz_0, g_0_xxxxxxx_0_xxxxyyz_1, g_0_xxxxxxx_0_xxxxyz_1, g_0_xxxxxxx_0_xxxxyzz_0, g_0_xxxxxxx_0_xxxxyzz_1, g_0_xxxxxxx_0_xxxxzz_1, g_0_xxxxxxx_0_xxxxzzz_0, g_0_xxxxxxx_0_xxxxzzz_1, g_0_xxxxxxx_0_xxxyyy_1, g_0_xxxxxxx_0_xxxyyyy_0, g_0_xxxxxxx_0_xxxyyyy_1, g_0_xxxxxxx_0_xxxyyyz_0, g_0_xxxxxxx_0_xxxyyyz_1, g_0_xxxxxxx_0_xxxyyz_1, g_0_xxxxxxx_0_xxxyyzz_0, g_0_xxxxxxx_0_xxxyyzz_1, g_0_xxxxxxx_0_xxxyzz_1, g_0_xxxxxxx_0_xxxyzzz_0, g_0_xxxxxxx_0_xxxyzzz_1, g_0_xxxxxxx_0_xxxzzz_1, g_0_xxxxxxx_0_xxxzzzz_0, g_0_xxxxxxx_0_xxxzzzz_1, g_0_xxxxxxx_0_xxyyyy_1, g_0_xxxxxxx_0_xxyyyyy_0, g_0_xxxxxxx_0_xxyyyyy_1, g_0_xxxxxxx_0_xxyyyyz_0, g_0_xxxxxxx_0_xxyyyyz_1, g_0_xxxxxxx_0_xxyyyz_1, g_0_xxxxxxx_0_xxyyyzz_0, g_0_xxxxxxx_0_xxyyyzz_1, g_0_xxxxxxx_0_xxyyzz_1, g_0_xxxxxxx_0_xxyyzzz_0, g_0_xxxxxxx_0_xxyyzzz_1, g_0_xxxxxxx_0_xxyzzz_1, g_0_xxxxxxx_0_xxyzzzz_0, g_0_xxxxxxx_0_xxyzzzz_1, g_0_xxxxxxx_0_xxzzzz_1, g_0_xxxxxxx_0_xxzzzzz_0, g_0_xxxxxxx_0_xxzzzzz_1, g_0_xxxxxxx_0_xyyyyy_1, g_0_xxxxxxx_0_xyyyyyy_0, g_0_xxxxxxx_0_xyyyyyy_1, g_0_xxxxxxx_0_xyyyyyz_0, g_0_xxxxxxx_0_xyyyyyz_1, g_0_xxxxxxx_0_xyyyyz_1, g_0_xxxxxxx_0_xyyyyzz_0, g_0_xxxxxxx_0_xyyyyzz_1, g_0_xxxxxxx_0_xyyyzz_1, g_0_xxxxxxx_0_xyyyzzz_0, g_0_xxxxxxx_0_xyyyzzz_1, g_0_xxxxxxx_0_xyyzzz_1, g_0_xxxxxxx_0_xyyzzzz_0, g_0_xxxxxxx_0_xyyzzzz_1, g_0_xxxxxxx_0_xyzzzz_1, g_0_xxxxxxx_0_xyzzzzz_0, g_0_xxxxxxx_0_xyzzzzz_1, g_0_xxxxxxx_0_xzzzzz_1, g_0_xxxxxxx_0_xzzzzzz_0, g_0_xxxxxxx_0_xzzzzzz_1, g_0_xxxxxxx_0_yyyyyy_1, g_0_xxxxxxx_0_yyyyyyy_0, g_0_xxxxxxx_0_yyyyyyy_1, g_0_xxxxxxx_0_yyyyyyz_0, g_0_xxxxxxx_0_yyyyyyz_1, g_0_xxxxxxx_0_yyyyyz_1, g_0_xxxxxxx_0_yyyyyzz_0, g_0_xxxxxxx_0_yyyyyzz_1, g_0_xxxxxxx_0_yyyyzz_1, g_0_xxxxxxx_0_yyyyzzz_0, g_0_xxxxxxx_0_yyyyzzz_1, g_0_xxxxxxx_0_yyyzzz_1, g_0_xxxxxxx_0_yyyzzzz_0, g_0_xxxxxxx_0_yyyzzzz_1, g_0_xxxxxxx_0_yyzzzz_1, g_0_xxxxxxx_0_yyzzzzz_0, g_0_xxxxxxx_0_yyzzzzz_1, g_0_xxxxxxx_0_yzzzzz_1, g_0_xxxxxxx_0_yzzzzzz_0, g_0_xxxxxxx_0_yzzzzzz_1, g_0_xxxxxxx_0_zzzzzz_1, g_0_xxxxxxx_0_zzzzzzz_0, g_0_xxxxxxx_0_zzzzzzz_1, g_0_xxxxxxxz_0_xxxxxxx_0, g_0_xxxxxxxz_0_xxxxxxy_0, g_0_xxxxxxxz_0_xxxxxxz_0, g_0_xxxxxxxz_0_xxxxxyy_0, g_0_xxxxxxxz_0_xxxxxyz_0, g_0_xxxxxxxz_0_xxxxxzz_0, g_0_xxxxxxxz_0_xxxxyyy_0, g_0_xxxxxxxz_0_xxxxyyz_0, g_0_xxxxxxxz_0_xxxxyzz_0, g_0_xxxxxxxz_0_xxxxzzz_0, g_0_xxxxxxxz_0_xxxyyyy_0, g_0_xxxxxxxz_0_xxxyyyz_0, g_0_xxxxxxxz_0_xxxyyzz_0, g_0_xxxxxxxz_0_xxxyzzz_0, g_0_xxxxxxxz_0_xxxzzzz_0, g_0_xxxxxxxz_0_xxyyyyy_0, g_0_xxxxxxxz_0_xxyyyyz_0, g_0_xxxxxxxz_0_xxyyyzz_0, g_0_xxxxxxxz_0_xxyyzzz_0, g_0_xxxxxxxz_0_xxyzzzz_0, g_0_xxxxxxxz_0_xxzzzzz_0, g_0_xxxxxxxz_0_xyyyyyy_0, g_0_xxxxxxxz_0_xyyyyyz_0, g_0_xxxxxxxz_0_xyyyyzz_0, g_0_xxxxxxxz_0_xyyyzzz_0, g_0_xxxxxxxz_0_xyyzzzz_0, g_0_xxxxxxxz_0_xyzzzzz_0, g_0_xxxxxxxz_0_xzzzzzz_0, g_0_xxxxxxxz_0_yyyyyyy_0, g_0_xxxxxxxz_0_yyyyyyz_0, g_0_xxxxxxxz_0_yyyyyzz_0, g_0_xxxxxxxz_0_yyyyzzz_0, g_0_xxxxxxxz_0_yyyzzzz_0, g_0_xxxxxxxz_0_yyzzzzz_0, g_0_xxxxxxxz_0_yzzzzzz_0, g_0_xxxxxxxz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxz_0_xxxxxxx_0[i] = g_0_xxxxxxx_0_xxxxxxx_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxxxy_0[i] = g_0_xxxxxxx_0_xxxxxxy_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxxxz_0[i] = g_0_xxxxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxxz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxxyy_0[i] = g_0_xxxxxxx_0_xxxxxyy_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxxyz_0[i] = g_0_xxxxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxyz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxyyy_0[i] = g_0_xxxxxxx_0_xxxxyyy_0[i] * pb_z + g_0_xxxxxxx_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxyyz_0[i] = g_0_xxxxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyyz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyyyy_0[i] = g_0_xxxxxxx_0_xxxyyyy_0[i] * pb_z + g_0_xxxxxxx_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyyyz_0[i] = g_0_xxxxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyyz_0[i] * pb_z + g_0_xxxxxxx_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyyyy_0[i] = g_0_xxxxxxx_0_xxyyyyy_0[i] * pb_z + g_0_xxxxxxx_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyyyz_0[i] = g_0_xxxxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyyz_0[i] * pb_z + g_0_xxxxxxx_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyzz_0[i] * pb_z + g_0_xxxxxxx_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyyyy_0[i] = g_0_xxxxxxx_0_xyyyyyy_0[i] * pb_z + g_0_xxxxxxx_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyyyz_0[i] = g_0_xxxxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyyz_0[i] * pb_z + g_0_xxxxxxx_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyzz_0[i] * pb_z + g_0_xxxxxxx_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyzzz_0[i] * pb_z + g_0_xxxxxxx_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyyyy_0[i] = g_0_xxxxxxx_0_yyyyyyy_0[i] * pb_z + g_0_xxxxxxx_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyyyz_0[i] = g_0_xxxxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyyyz_0[i] * pb_z + g_0_xxxxxxx_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyyzz_0[i] * pb_z + g_0_xxxxxxx_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyzzz_0[i] * pb_z + g_0_xxxxxxx_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyzzzz_0[i] * pb_z + g_0_xxxxxxx_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_zzzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 108-144 components of targeted buffer : SLSK

    auto g_0_xxxxxxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 108);

    auto g_0_xxxxxxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 109);

    auto g_0_xxxxxxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 110);

    auto g_0_xxxxxxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 111);

    auto g_0_xxxxxxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 112);

    auto g_0_xxxxxxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 113);

    auto g_0_xxxxxxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 114);

    auto g_0_xxxxxxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 115);

    auto g_0_xxxxxxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 116);

    auto g_0_xxxxxxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 117);

    auto g_0_xxxxxxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 118);

    auto g_0_xxxxxxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 119);

    auto g_0_xxxxxxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 120);

    auto g_0_xxxxxxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 121);

    auto g_0_xxxxxxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 122);

    auto g_0_xxxxxxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 123);

    auto g_0_xxxxxxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 124);

    auto g_0_xxxxxxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 125);

    auto g_0_xxxxxxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 126);

    auto g_0_xxxxxxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 127);

    auto g_0_xxxxxxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 128);

    auto g_0_xxxxxxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 129);

    auto g_0_xxxxxxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 130);

    auto g_0_xxxxxxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 131);

    auto g_0_xxxxxxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 132);

    auto g_0_xxxxxxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 133);

    auto g_0_xxxxxxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 134);

    auto g_0_xxxxxxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 135);

    auto g_0_xxxxxxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 136);

    auto g_0_xxxxxxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 137);

    auto g_0_xxxxxxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 138);

    auto g_0_xxxxxxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 139);

    auto g_0_xxxxxxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 140);

    auto g_0_xxxxxxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 141);

    auto g_0_xxxxxxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 142);

    auto g_0_xxxxxxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 143);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxxxxx_0, g_0_xxxxxx_0_xxxxxxx_1, g_0_xxxxxx_0_xxxxxxz_0, g_0_xxxxxx_0_xxxxxxz_1, g_0_xxxxxx_0_xxxxxzz_0, g_0_xxxxxx_0_xxxxxzz_1, g_0_xxxxxx_0_xxxxzzz_0, g_0_xxxxxx_0_xxxxzzz_1, g_0_xxxxxx_0_xxxzzzz_0, g_0_xxxxxx_0_xxxzzzz_1, g_0_xxxxxx_0_xxzzzzz_0, g_0_xxxxxx_0_xxzzzzz_1, g_0_xxxxxx_0_xzzzzzz_0, g_0_xxxxxx_0_xzzzzzz_1, g_0_xxxxxxy_0_xxxxxxx_0, g_0_xxxxxxy_0_xxxxxxx_1, g_0_xxxxxxy_0_xxxxxxz_0, g_0_xxxxxxy_0_xxxxxxz_1, g_0_xxxxxxy_0_xxxxxzz_0, g_0_xxxxxxy_0_xxxxxzz_1, g_0_xxxxxxy_0_xxxxzzz_0, g_0_xxxxxxy_0_xxxxzzz_1, g_0_xxxxxxy_0_xxxzzzz_0, g_0_xxxxxxy_0_xxxzzzz_1, g_0_xxxxxxy_0_xxzzzzz_0, g_0_xxxxxxy_0_xxzzzzz_1, g_0_xxxxxxy_0_xzzzzzz_0, g_0_xxxxxxy_0_xzzzzzz_1, g_0_xxxxxxyy_0_xxxxxxx_0, g_0_xxxxxxyy_0_xxxxxxy_0, g_0_xxxxxxyy_0_xxxxxxz_0, g_0_xxxxxxyy_0_xxxxxyy_0, g_0_xxxxxxyy_0_xxxxxyz_0, g_0_xxxxxxyy_0_xxxxxzz_0, g_0_xxxxxxyy_0_xxxxyyy_0, g_0_xxxxxxyy_0_xxxxyyz_0, g_0_xxxxxxyy_0_xxxxyzz_0, g_0_xxxxxxyy_0_xxxxzzz_0, g_0_xxxxxxyy_0_xxxyyyy_0, g_0_xxxxxxyy_0_xxxyyyz_0, g_0_xxxxxxyy_0_xxxyyzz_0, g_0_xxxxxxyy_0_xxxyzzz_0, g_0_xxxxxxyy_0_xxxzzzz_0, g_0_xxxxxxyy_0_xxyyyyy_0, g_0_xxxxxxyy_0_xxyyyyz_0, g_0_xxxxxxyy_0_xxyyyzz_0, g_0_xxxxxxyy_0_xxyyzzz_0, g_0_xxxxxxyy_0_xxyzzzz_0, g_0_xxxxxxyy_0_xxzzzzz_0, g_0_xxxxxxyy_0_xyyyyyy_0, g_0_xxxxxxyy_0_xyyyyyz_0, g_0_xxxxxxyy_0_xyyyyzz_0, g_0_xxxxxxyy_0_xyyyzzz_0, g_0_xxxxxxyy_0_xyyzzzz_0, g_0_xxxxxxyy_0_xyzzzzz_0, g_0_xxxxxxyy_0_xzzzzzz_0, g_0_xxxxxxyy_0_yyyyyyy_0, g_0_xxxxxxyy_0_yyyyyyz_0, g_0_xxxxxxyy_0_yyyyyzz_0, g_0_xxxxxxyy_0_yyyyzzz_0, g_0_xxxxxxyy_0_yyyzzzz_0, g_0_xxxxxxyy_0_yyzzzzz_0, g_0_xxxxxxyy_0_yzzzzzz_0, g_0_xxxxxxyy_0_zzzzzzz_0, g_0_xxxxxyy_0_xxxxxxy_0, g_0_xxxxxyy_0_xxxxxxy_1, g_0_xxxxxyy_0_xxxxxy_1, g_0_xxxxxyy_0_xxxxxyy_0, g_0_xxxxxyy_0_xxxxxyy_1, g_0_xxxxxyy_0_xxxxxyz_0, g_0_xxxxxyy_0_xxxxxyz_1, g_0_xxxxxyy_0_xxxxyy_1, g_0_xxxxxyy_0_xxxxyyy_0, g_0_xxxxxyy_0_xxxxyyy_1, g_0_xxxxxyy_0_xxxxyyz_0, g_0_xxxxxyy_0_xxxxyyz_1, g_0_xxxxxyy_0_xxxxyz_1, g_0_xxxxxyy_0_xxxxyzz_0, g_0_xxxxxyy_0_xxxxyzz_1, g_0_xxxxxyy_0_xxxyyy_1, g_0_xxxxxyy_0_xxxyyyy_0, g_0_xxxxxyy_0_xxxyyyy_1, g_0_xxxxxyy_0_xxxyyyz_0, g_0_xxxxxyy_0_xxxyyyz_1, g_0_xxxxxyy_0_xxxyyz_1, g_0_xxxxxyy_0_xxxyyzz_0, g_0_xxxxxyy_0_xxxyyzz_1, g_0_xxxxxyy_0_xxxyzz_1, g_0_xxxxxyy_0_xxxyzzz_0, g_0_xxxxxyy_0_xxxyzzz_1, g_0_xxxxxyy_0_xxyyyy_1, g_0_xxxxxyy_0_xxyyyyy_0, g_0_xxxxxyy_0_xxyyyyy_1, g_0_xxxxxyy_0_xxyyyyz_0, g_0_xxxxxyy_0_xxyyyyz_1, g_0_xxxxxyy_0_xxyyyz_1, g_0_xxxxxyy_0_xxyyyzz_0, g_0_xxxxxyy_0_xxyyyzz_1, g_0_xxxxxyy_0_xxyyzz_1, g_0_xxxxxyy_0_xxyyzzz_0, g_0_xxxxxyy_0_xxyyzzz_1, g_0_xxxxxyy_0_xxyzzz_1, g_0_xxxxxyy_0_xxyzzzz_0, g_0_xxxxxyy_0_xxyzzzz_1, g_0_xxxxxyy_0_xyyyyy_1, g_0_xxxxxyy_0_xyyyyyy_0, g_0_xxxxxyy_0_xyyyyyy_1, g_0_xxxxxyy_0_xyyyyyz_0, g_0_xxxxxyy_0_xyyyyyz_1, g_0_xxxxxyy_0_xyyyyz_1, g_0_xxxxxyy_0_xyyyyzz_0, g_0_xxxxxyy_0_xyyyyzz_1, g_0_xxxxxyy_0_xyyyzz_1, g_0_xxxxxyy_0_xyyyzzz_0, g_0_xxxxxyy_0_xyyyzzz_1, g_0_xxxxxyy_0_xyyzzz_1, g_0_xxxxxyy_0_xyyzzzz_0, g_0_xxxxxyy_0_xyyzzzz_1, g_0_xxxxxyy_0_xyzzzz_1, g_0_xxxxxyy_0_xyzzzzz_0, g_0_xxxxxyy_0_xyzzzzz_1, g_0_xxxxxyy_0_yyyyyy_1, g_0_xxxxxyy_0_yyyyyyy_0, g_0_xxxxxyy_0_yyyyyyy_1, g_0_xxxxxyy_0_yyyyyyz_0, g_0_xxxxxyy_0_yyyyyyz_1, g_0_xxxxxyy_0_yyyyyz_1, g_0_xxxxxyy_0_yyyyyzz_0, g_0_xxxxxyy_0_yyyyyzz_1, g_0_xxxxxyy_0_yyyyzz_1, g_0_xxxxxyy_0_yyyyzzz_0, g_0_xxxxxyy_0_yyyyzzz_1, g_0_xxxxxyy_0_yyyzzz_1, g_0_xxxxxyy_0_yyyzzzz_0, g_0_xxxxxyy_0_yyyzzzz_1, g_0_xxxxxyy_0_yyzzzz_1, g_0_xxxxxyy_0_yyzzzzz_0, g_0_xxxxxyy_0_yyzzzzz_1, g_0_xxxxxyy_0_yzzzzz_1, g_0_xxxxxyy_0_yzzzzzz_0, g_0_xxxxxyy_0_yzzzzzz_1, g_0_xxxxxyy_0_zzzzzzz_0, g_0_xxxxxyy_0_zzzzzzz_1, g_0_xxxxyy_0_xxxxxxy_0, g_0_xxxxyy_0_xxxxxxy_1, g_0_xxxxyy_0_xxxxxyy_0, g_0_xxxxyy_0_xxxxxyy_1, g_0_xxxxyy_0_xxxxxyz_0, g_0_xxxxyy_0_xxxxxyz_1, g_0_xxxxyy_0_xxxxyyy_0, g_0_xxxxyy_0_xxxxyyy_1, g_0_xxxxyy_0_xxxxyyz_0, g_0_xxxxyy_0_xxxxyyz_1, g_0_xxxxyy_0_xxxxyzz_0, g_0_xxxxyy_0_xxxxyzz_1, g_0_xxxxyy_0_xxxyyyy_0, g_0_xxxxyy_0_xxxyyyy_1, g_0_xxxxyy_0_xxxyyyz_0, g_0_xxxxyy_0_xxxyyyz_1, g_0_xxxxyy_0_xxxyyzz_0, g_0_xxxxyy_0_xxxyyzz_1, g_0_xxxxyy_0_xxxyzzz_0, g_0_xxxxyy_0_xxxyzzz_1, g_0_xxxxyy_0_xxyyyyy_0, g_0_xxxxyy_0_xxyyyyy_1, g_0_xxxxyy_0_xxyyyyz_0, g_0_xxxxyy_0_xxyyyyz_1, g_0_xxxxyy_0_xxyyyzz_0, g_0_xxxxyy_0_xxyyyzz_1, g_0_xxxxyy_0_xxyyzzz_0, g_0_xxxxyy_0_xxyyzzz_1, g_0_xxxxyy_0_xxyzzzz_0, g_0_xxxxyy_0_xxyzzzz_1, g_0_xxxxyy_0_xyyyyyy_0, g_0_xxxxyy_0_xyyyyyy_1, g_0_xxxxyy_0_xyyyyyz_0, g_0_xxxxyy_0_xyyyyyz_1, g_0_xxxxyy_0_xyyyyzz_0, g_0_xxxxyy_0_xyyyyzz_1, g_0_xxxxyy_0_xyyyzzz_0, g_0_xxxxyy_0_xyyyzzz_1, g_0_xxxxyy_0_xyyzzzz_0, g_0_xxxxyy_0_xyyzzzz_1, g_0_xxxxyy_0_xyzzzzz_0, g_0_xxxxyy_0_xyzzzzz_1, g_0_xxxxyy_0_yyyyyyy_0, g_0_xxxxyy_0_yyyyyyy_1, g_0_xxxxyy_0_yyyyyyz_0, g_0_xxxxyy_0_yyyyyyz_1, g_0_xxxxyy_0_yyyyyzz_0, g_0_xxxxyy_0_yyyyyzz_1, g_0_xxxxyy_0_yyyyzzz_0, g_0_xxxxyy_0_yyyyzzz_1, g_0_xxxxyy_0_yyyzzzz_0, g_0_xxxxyy_0_yyyzzzz_1, g_0_xxxxyy_0_yyzzzzz_0, g_0_xxxxyy_0_yyzzzzz_1, g_0_xxxxyy_0_yzzzzzz_0, g_0_xxxxyy_0_yzzzzzz_1, g_0_xxxxyy_0_zzzzzzz_0, g_0_xxxxyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxyy_0_xxxxxxx_0[i] = g_0_xxxxxx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxxxx_0[i] * pb_y + g_0_xxxxxxy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxxxxy_0[i] = 5.0 * g_0_xxxxyy_0_xxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxxxyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxxxy_0[i] * pb_x + g_0_xxxxxyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxxxz_0[i] = g_0_xxxxxx_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxxxz_0[i] * pb_y + g_0_xxxxxxy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxxxyy_0[i] = 5.0 * g_0_xxxxyy_0_xxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxxyy_0[i] * pb_x + g_0_xxxxxyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxxyz_0[i] = 5.0 * g_0_xxxxyy_0_xxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxxyz_0[i] * pb_x + g_0_xxxxxyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxxzz_0[i] = g_0_xxxxxx_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxxzz_0[i] * pb_y + g_0_xxxxxxy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxxyyy_0[i] = 5.0 * g_0_xxxxyy_0_xxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxyyy_0[i] * pb_x + g_0_xxxxxyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxyyz_0[i] = 5.0 * g_0_xxxxyy_0_xxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxyyz_0[i] * pb_x + g_0_xxxxxyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxyzz_0[i] = 5.0 * g_0_xxxxyy_0_xxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxyzz_0[i] * pb_x + g_0_xxxxxyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxzzz_0[i] = g_0_xxxxxx_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxzzz_0[i] * pb_y + g_0_xxxxxxy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxyyyy_0[i] = 5.0 * g_0_xxxxyy_0_xxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyyyy_0[i] * pb_x + g_0_xxxxxyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxyyyz_0[i] = 5.0 * g_0_xxxxyy_0_xxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyyyz_0[i] * pb_x + g_0_xxxxxyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxyyzz_0[i] = 5.0 * g_0_xxxxyy_0_xxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyyzz_0[i] * pb_x + g_0_xxxxxyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxyzzz_0[i] = 5.0 * g_0_xxxxyy_0_xxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyzzz_0[i] * pb_x + g_0_xxxxxyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxzzzz_0[i] = g_0_xxxxxx_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxzzzz_0[i] * pb_y + g_0_xxxxxxy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxyyyyy_0[i] = 5.0 * g_0_xxxxyy_0_xxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyyyy_0[i] * pb_x + g_0_xxxxxyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyyyyz_0[i] = 5.0 * g_0_xxxxyy_0_xxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyyyz_0[i] * pb_x + g_0_xxxxxyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyyyzz_0[i] = 5.0 * g_0_xxxxyy_0_xxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyyzz_0[i] * pb_x + g_0_xxxxxyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyyzzz_0[i] = 5.0 * g_0_xxxxyy_0_xxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyzzz_0[i] * pb_x + g_0_xxxxxyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyzzzz_0[i] = 5.0 * g_0_xxxxyy_0_xxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyzzzz_0[i] * pb_x + g_0_xxxxxyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxzzzzz_0[i] = g_0_xxxxxx_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxzzzzz_0[i] * pb_y + g_0_xxxxxxy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xyyyyyy_0[i] = 5.0 * g_0_xxxxyy_0_xyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyyyy_0[i] * pb_x + g_0_xxxxxyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxyy_0_xyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyyyz_0[i] * pb_x + g_0_xxxxxyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyyyzz_0[i] = 5.0 * g_0_xxxxyy_0_xyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyyzz_0[i] * pb_x + g_0_xxxxxyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyyzzz_0[i] = 5.0 * g_0_xxxxyy_0_xyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyzzz_0[i] * pb_x + g_0_xxxxxyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyzzzz_0[i] = 5.0 * g_0_xxxxyy_0_xyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyzzzz_0[i] * pb_x + g_0_xxxxxyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_xyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyzzzzz_0[i] * pb_x + g_0_xxxxxyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xzzzzzz_0[i] = g_0_xxxxxx_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xzzzzzz_0[i] * pb_y + g_0_xxxxxxy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_yyyyyyy_0[i] = 5.0 * g_0_xxxxyy_0_yyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyyyy_0[i] * pb_x + g_0_xxxxxyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyyyyz_0[i] = 5.0 * g_0_xxxxyy_0_yyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyyyz_0[i] * pb_x + g_0_xxxxxyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxyy_0_yyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyyzz_0[i] * pb_x + g_0_xxxxxyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyyzzz_0[i] = 5.0 * g_0_xxxxyy_0_yyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyzzz_0[i] * pb_x + g_0_xxxxxyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyzzzz_0[i] = 5.0 * g_0_xxxxyy_0_yyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyzzzz_0[i] * pb_x + g_0_xxxxxyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_yyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyzzzzz_0[i] * pb_x + g_0_xxxxxyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yzzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_yzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yzzzzzz_0[i] * pb_x + g_0_xxxxxyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_zzzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_zzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_zzzzzzz_0[i] * pb_x + g_0_xxxxxyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 144-180 components of targeted buffer : SLSK

    auto g_0_xxxxxxyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 144);

    auto g_0_xxxxxxyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 145);

    auto g_0_xxxxxxyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 146);

    auto g_0_xxxxxxyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 147);

    auto g_0_xxxxxxyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 148);

    auto g_0_xxxxxxyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 149);

    auto g_0_xxxxxxyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 150);

    auto g_0_xxxxxxyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 151);

    auto g_0_xxxxxxyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 152);

    auto g_0_xxxxxxyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 153);

    auto g_0_xxxxxxyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 154);

    auto g_0_xxxxxxyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 155);

    auto g_0_xxxxxxyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 156);

    auto g_0_xxxxxxyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 157);

    auto g_0_xxxxxxyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 158);

    auto g_0_xxxxxxyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 159);

    auto g_0_xxxxxxyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 160);

    auto g_0_xxxxxxyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 161);

    auto g_0_xxxxxxyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 162);

    auto g_0_xxxxxxyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 163);

    auto g_0_xxxxxxyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 164);

    auto g_0_xxxxxxyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 165);

    auto g_0_xxxxxxyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 166);

    auto g_0_xxxxxxyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 167);

    auto g_0_xxxxxxyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 168);

    auto g_0_xxxxxxyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 169);

    auto g_0_xxxxxxyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 170);

    auto g_0_xxxxxxyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 171);

    auto g_0_xxxxxxyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 172);

    auto g_0_xxxxxxyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 173);

    auto g_0_xxxxxxyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 174);

    auto g_0_xxxxxxyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 175);

    auto g_0_xxxxxxyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 176);

    auto g_0_xxxxxxyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 177);

    auto g_0_xxxxxxyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 178);

    auto g_0_xxxxxxyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 179);

    #pragma omp simd aligned(g_0_xxxxxxy_0_xxxxxxy_0, g_0_xxxxxxy_0_xxxxxxy_1, g_0_xxxxxxy_0_xxxxxyy_0, g_0_xxxxxxy_0_xxxxxyy_1, g_0_xxxxxxy_0_xxxxyyy_0, g_0_xxxxxxy_0_xxxxyyy_1, g_0_xxxxxxy_0_xxxyyyy_0, g_0_xxxxxxy_0_xxxyyyy_1, g_0_xxxxxxy_0_xxyyyyy_0, g_0_xxxxxxy_0_xxyyyyy_1, g_0_xxxxxxy_0_xyyyyyy_0, g_0_xxxxxxy_0_xyyyyyy_1, g_0_xxxxxxy_0_yyyyyyy_0, g_0_xxxxxxy_0_yyyyyyy_1, g_0_xxxxxxyz_0_xxxxxxx_0, g_0_xxxxxxyz_0_xxxxxxy_0, g_0_xxxxxxyz_0_xxxxxxz_0, g_0_xxxxxxyz_0_xxxxxyy_0, g_0_xxxxxxyz_0_xxxxxyz_0, g_0_xxxxxxyz_0_xxxxxzz_0, g_0_xxxxxxyz_0_xxxxyyy_0, g_0_xxxxxxyz_0_xxxxyyz_0, g_0_xxxxxxyz_0_xxxxyzz_0, g_0_xxxxxxyz_0_xxxxzzz_0, g_0_xxxxxxyz_0_xxxyyyy_0, g_0_xxxxxxyz_0_xxxyyyz_0, g_0_xxxxxxyz_0_xxxyyzz_0, g_0_xxxxxxyz_0_xxxyzzz_0, g_0_xxxxxxyz_0_xxxzzzz_0, g_0_xxxxxxyz_0_xxyyyyy_0, g_0_xxxxxxyz_0_xxyyyyz_0, g_0_xxxxxxyz_0_xxyyyzz_0, g_0_xxxxxxyz_0_xxyyzzz_0, g_0_xxxxxxyz_0_xxyzzzz_0, g_0_xxxxxxyz_0_xxzzzzz_0, g_0_xxxxxxyz_0_xyyyyyy_0, g_0_xxxxxxyz_0_xyyyyyz_0, g_0_xxxxxxyz_0_xyyyyzz_0, g_0_xxxxxxyz_0_xyyyzzz_0, g_0_xxxxxxyz_0_xyyzzzz_0, g_0_xxxxxxyz_0_xyzzzzz_0, g_0_xxxxxxyz_0_xzzzzzz_0, g_0_xxxxxxyz_0_yyyyyyy_0, g_0_xxxxxxyz_0_yyyyyyz_0, g_0_xxxxxxyz_0_yyyyyzz_0, g_0_xxxxxxyz_0_yyyyzzz_0, g_0_xxxxxxyz_0_yyyzzzz_0, g_0_xxxxxxyz_0_yyzzzzz_0, g_0_xxxxxxyz_0_yzzzzzz_0, g_0_xxxxxxyz_0_zzzzzzz_0, g_0_xxxxxxz_0_xxxxxxx_0, g_0_xxxxxxz_0_xxxxxxx_1, g_0_xxxxxxz_0_xxxxxxz_0, g_0_xxxxxxz_0_xxxxxxz_1, g_0_xxxxxxz_0_xxxxxyz_0, g_0_xxxxxxz_0_xxxxxyz_1, g_0_xxxxxxz_0_xxxxxz_1, g_0_xxxxxxz_0_xxxxxzz_0, g_0_xxxxxxz_0_xxxxxzz_1, g_0_xxxxxxz_0_xxxxyyz_0, g_0_xxxxxxz_0_xxxxyyz_1, g_0_xxxxxxz_0_xxxxyz_1, g_0_xxxxxxz_0_xxxxyzz_0, g_0_xxxxxxz_0_xxxxyzz_1, g_0_xxxxxxz_0_xxxxzz_1, g_0_xxxxxxz_0_xxxxzzz_0, g_0_xxxxxxz_0_xxxxzzz_1, g_0_xxxxxxz_0_xxxyyyz_0, g_0_xxxxxxz_0_xxxyyyz_1, g_0_xxxxxxz_0_xxxyyz_1, g_0_xxxxxxz_0_xxxyyzz_0, g_0_xxxxxxz_0_xxxyyzz_1, g_0_xxxxxxz_0_xxxyzz_1, g_0_xxxxxxz_0_xxxyzzz_0, g_0_xxxxxxz_0_xxxyzzz_1, g_0_xxxxxxz_0_xxxzzz_1, g_0_xxxxxxz_0_xxxzzzz_0, g_0_xxxxxxz_0_xxxzzzz_1, g_0_xxxxxxz_0_xxyyyyz_0, g_0_xxxxxxz_0_xxyyyyz_1, g_0_xxxxxxz_0_xxyyyz_1, g_0_xxxxxxz_0_xxyyyzz_0, g_0_xxxxxxz_0_xxyyyzz_1, g_0_xxxxxxz_0_xxyyzz_1, g_0_xxxxxxz_0_xxyyzzz_0, g_0_xxxxxxz_0_xxyyzzz_1, g_0_xxxxxxz_0_xxyzzz_1, g_0_xxxxxxz_0_xxyzzzz_0, g_0_xxxxxxz_0_xxyzzzz_1, g_0_xxxxxxz_0_xxzzzz_1, g_0_xxxxxxz_0_xxzzzzz_0, g_0_xxxxxxz_0_xxzzzzz_1, g_0_xxxxxxz_0_xyyyyyz_0, g_0_xxxxxxz_0_xyyyyyz_1, g_0_xxxxxxz_0_xyyyyz_1, g_0_xxxxxxz_0_xyyyyzz_0, g_0_xxxxxxz_0_xyyyyzz_1, g_0_xxxxxxz_0_xyyyzz_1, g_0_xxxxxxz_0_xyyyzzz_0, g_0_xxxxxxz_0_xyyyzzz_1, g_0_xxxxxxz_0_xyyzzz_1, g_0_xxxxxxz_0_xyyzzzz_0, g_0_xxxxxxz_0_xyyzzzz_1, g_0_xxxxxxz_0_xyzzzz_1, g_0_xxxxxxz_0_xyzzzzz_0, g_0_xxxxxxz_0_xyzzzzz_1, g_0_xxxxxxz_0_xzzzzz_1, g_0_xxxxxxz_0_xzzzzzz_0, g_0_xxxxxxz_0_xzzzzzz_1, g_0_xxxxxxz_0_yyyyyyz_0, g_0_xxxxxxz_0_yyyyyyz_1, g_0_xxxxxxz_0_yyyyyz_1, g_0_xxxxxxz_0_yyyyyzz_0, g_0_xxxxxxz_0_yyyyyzz_1, g_0_xxxxxxz_0_yyyyzz_1, g_0_xxxxxxz_0_yyyyzzz_0, g_0_xxxxxxz_0_yyyyzzz_1, g_0_xxxxxxz_0_yyyzzz_1, g_0_xxxxxxz_0_yyyzzzz_0, g_0_xxxxxxz_0_yyyzzzz_1, g_0_xxxxxxz_0_yyzzzz_1, g_0_xxxxxxz_0_yyzzzzz_0, g_0_xxxxxxz_0_yyzzzzz_1, g_0_xxxxxxz_0_yzzzzz_1, g_0_xxxxxxz_0_yzzzzzz_0, g_0_xxxxxxz_0_yzzzzzz_1, g_0_xxxxxxz_0_zzzzzz_1, g_0_xxxxxxz_0_zzzzzzz_0, g_0_xxxxxxz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxyz_0_xxxxxxx_0[i] = g_0_xxxxxxz_0_xxxxxxx_0[i] * pb_y + g_0_xxxxxxz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxxxy_0[i] = g_0_xxxxxxy_0_xxxxxxy_0[i] * pb_z + g_0_xxxxxxy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxxxxz_0[i] = g_0_xxxxxxz_0_xxxxxxz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxxyy_0[i] = g_0_xxxxxxy_0_xxxxxyy_0[i] * pb_z + g_0_xxxxxxy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxxxyz_0[i] = g_0_xxxxxxz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxxxyz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxxzz_0[i] = g_0_xxxxxxz_0_xxxxxzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxyyy_0[i] = g_0_xxxxxxy_0_xxxxyyy_0[i] * pb_z + g_0_xxxxxxy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxxxz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxxyyz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxyzz_0[i] = g_0_xxxxxxz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxxyzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxzzz_0[i] = g_0_xxxxxxz_0_xxxxzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxyyyy_0[i] = g_0_xxxxxxy_0_xxxyyyy_0[i] * pb_z + g_0_xxxxxxy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxxxz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxyyyz_0[i] * pb_y + g_0_xxxxxxz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxxxz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxyyzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxyzzz_0[i] = g_0_xxxxxxz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxyzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxzzzz_0[i] = g_0_xxxxxxz_0_xxxzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyyyyy_0[i] = g_0_xxxxxxy_0_xxyyyyy_0[i] * pb_z + g_0_xxxxxxy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxxxz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyyyyz_0[i] * pb_y + g_0_xxxxxxz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxxxz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyyyzz_0[i] * pb_y + g_0_xxxxxxz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxxxz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyyzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyzzzz_0[i] = g_0_xxxxxxz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxzzzzz_0[i] = g_0_xxxxxxz_0_xxzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyyyyy_0[i] = g_0_xxxxxxy_0_xyyyyyy_0[i] * pb_z + g_0_xxxxxxy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxxxz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyyyyz_0[i] * pb_y + g_0_xxxxxxz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxxxz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyyyzz_0[i] * pb_y + g_0_xxxxxxz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxxxz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyyzzz_0[i] * pb_y + g_0_xxxxxxz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxxxz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyzzzzz_0[i] = g_0_xxxxxxz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xzzzzzz_0[i] = g_0_xxxxxxz_0_xzzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyyyyy_0[i] = g_0_xxxxxxy_0_yyyyyyy_0[i] * pb_z + g_0_xxxxxxy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxxxz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyyyyz_0[i] * pb_y + g_0_xxxxxxz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxxxz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyyyzz_0[i] * pb_y + g_0_xxxxxxz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxxxz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyyzzz_0[i] * pb_y + g_0_xxxxxxz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxxxz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyzzzz_0[i] * pb_y + g_0_xxxxxxz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxxxz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yzzzzzz_0[i] = g_0_xxxxxxz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yzzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_zzzzzzz_0[i] = g_0_xxxxxxz_0_zzzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 180-216 components of targeted buffer : SLSK

    auto g_0_xxxxxxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 180);

    auto g_0_xxxxxxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 181);

    auto g_0_xxxxxxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 182);

    auto g_0_xxxxxxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 183);

    auto g_0_xxxxxxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 184);

    auto g_0_xxxxxxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 185);

    auto g_0_xxxxxxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 186);

    auto g_0_xxxxxxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 187);

    auto g_0_xxxxxxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 188);

    auto g_0_xxxxxxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 189);

    auto g_0_xxxxxxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 190);

    auto g_0_xxxxxxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 191);

    auto g_0_xxxxxxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 192);

    auto g_0_xxxxxxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 193);

    auto g_0_xxxxxxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 194);

    auto g_0_xxxxxxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 195);

    auto g_0_xxxxxxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 196);

    auto g_0_xxxxxxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 197);

    auto g_0_xxxxxxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 198);

    auto g_0_xxxxxxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 199);

    auto g_0_xxxxxxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 200);

    auto g_0_xxxxxxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 201);

    auto g_0_xxxxxxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 202);

    auto g_0_xxxxxxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 203);

    auto g_0_xxxxxxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 204);

    auto g_0_xxxxxxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 205);

    auto g_0_xxxxxxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 206);

    auto g_0_xxxxxxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 207);

    auto g_0_xxxxxxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 208);

    auto g_0_xxxxxxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 209);

    auto g_0_xxxxxxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 210);

    auto g_0_xxxxxxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 211);

    auto g_0_xxxxxxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 212);

    auto g_0_xxxxxxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 213);

    auto g_0_xxxxxxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 214);

    auto g_0_xxxxxxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 215);

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxxxxx_0, g_0_xxxxxx_0_xxxxxxx_1, g_0_xxxxxx_0_xxxxxxy_0, g_0_xxxxxx_0_xxxxxxy_1, g_0_xxxxxx_0_xxxxxyy_0, g_0_xxxxxx_0_xxxxxyy_1, g_0_xxxxxx_0_xxxxyyy_0, g_0_xxxxxx_0_xxxxyyy_1, g_0_xxxxxx_0_xxxyyyy_0, g_0_xxxxxx_0_xxxyyyy_1, g_0_xxxxxx_0_xxyyyyy_0, g_0_xxxxxx_0_xxyyyyy_1, g_0_xxxxxx_0_xyyyyyy_0, g_0_xxxxxx_0_xyyyyyy_1, g_0_xxxxxxz_0_xxxxxxx_0, g_0_xxxxxxz_0_xxxxxxx_1, g_0_xxxxxxz_0_xxxxxxy_0, g_0_xxxxxxz_0_xxxxxxy_1, g_0_xxxxxxz_0_xxxxxyy_0, g_0_xxxxxxz_0_xxxxxyy_1, g_0_xxxxxxz_0_xxxxyyy_0, g_0_xxxxxxz_0_xxxxyyy_1, g_0_xxxxxxz_0_xxxyyyy_0, g_0_xxxxxxz_0_xxxyyyy_1, g_0_xxxxxxz_0_xxyyyyy_0, g_0_xxxxxxz_0_xxyyyyy_1, g_0_xxxxxxz_0_xyyyyyy_0, g_0_xxxxxxz_0_xyyyyyy_1, g_0_xxxxxxzz_0_xxxxxxx_0, g_0_xxxxxxzz_0_xxxxxxy_0, g_0_xxxxxxzz_0_xxxxxxz_0, g_0_xxxxxxzz_0_xxxxxyy_0, g_0_xxxxxxzz_0_xxxxxyz_0, g_0_xxxxxxzz_0_xxxxxzz_0, g_0_xxxxxxzz_0_xxxxyyy_0, g_0_xxxxxxzz_0_xxxxyyz_0, g_0_xxxxxxzz_0_xxxxyzz_0, g_0_xxxxxxzz_0_xxxxzzz_0, g_0_xxxxxxzz_0_xxxyyyy_0, g_0_xxxxxxzz_0_xxxyyyz_0, g_0_xxxxxxzz_0_xxxyyzz_0, g_0_xxxxxxzz_0_xxxyzzz_0, g_0_xxxxxxzz_0_xxxzzzz_0, g_0_xxxxxxzz_0_xxyyyyy_0, g_0_xxxxxxzz_0_xxyyyyz_0, g_0_xxxxxxzz_0_xxyyyzz_0, g_0_xxxxxxzz_0_xxyyzzz_0, g_0_xxxxxxzz_0_xxyzzzz_0, g_0_xxxxxxzz_0_xxzzzzz_0, g_0_xxxxxxzz_0_xyyyyyy_0, g_0_xxxxxxzz_0_xyyyyyz_0, g_0_xxxxxxzz_0_xyyyyzz_0, g_0_xxxxxxzz_0_xyyyzzz_0, g_0_xxxxxxzz_0_xyyzzzz_0, g_0_xxxxxxzz_0_xyzzzzz_0, g_0_xxxxxxzz_0_xzzzzzz_0, g_0_xxxxxxzz_0_yyyyyyy_0, g_0_xxxxxxzz_0_yyyyyyz_0, g_0_xxxxxxzz_0_yyyyyzz_0, g_0_xxxxxxzz_0_yyyyzzz_0, g_0_xxxxxxzz_0_yyyzzzz_0, g_0_xxxxxxzz_0_yyzzzzz_0, g_0_xxxxxxzz_0_yzzzzzz_0, g_0_xxxxxxzz_0_zzzzzzz_0, g_0_xxxxxzz_0_xxxxxxz_0, g_0_xxxxxzz_0_xxxxxxz_1, g_0_xxxxxzz_0_xxxxxyz_0, g_0_xxxxxzz_0_xxxxxyz_1, g_0_xxxxxzz_0_xxxxxz_1, g_0_xxxxxzz_0_xxxxxzz_0, g_0_xxxxxzz_0_xxxxxzz_1, g_0_xxxxxzz_0_xxxxyyz_0, g_0_xxxxxzz_0_xxxxyyz_1, g_0_xxxxxzz_0_xxxxyz_1, g_0_xxxxxzz_0_xxxxyzz_0, g_0_xxxxxzz_0_xxxxyzz_1, g_0_xxxxxzz_0_xxxxzz_1, g_0_xxxxxzz_0_xxxxzzz_0, g_0_xxxxxzz_0_xxxxzzz_1, g_0_xxxxxzz_0_xxxyyyz_0, g_0_xxxxxzz_0_xxxyyyz_1, g_0_xxxxxzz_0_xxxyyz_1, g_0_xxxxxzz_0_xxxyyzz_0, g_0_xxxxxzz_0_xxxyyzz_1, g_0_xxxxxzz_0_xxxyzz_1, g_0_xxxxxzz_0_xxxyzzz_0, g_0_xxxxxzz_0_xxxyzzz_1, g_0_xxxxxzz_0_xxxzzz_1, g_0_xxxxxzz_0_xxxzzzz_0, g_0_xxxxxzz_0_xxxzzzz_1, g_0_xxxxxzz_0_xxyyyyz_0, g_0_xxxxxzz_0_xxyyyyz_1, g_0_xxxxxzz_0_xxyyyz_1, g_0_xxxxxzz_0_xxyyyzz_0, g_0_xxxxxzz_0_xxyyyzz_1, g_0_xxxxxzz_0_xxyyzz_1, g_0_xxxxxzz_0_xxyyzzz_0, g_0_xxxxxzz_0_xxyyzzz_1, g_0_xxxxxzz_0_xxyzzz_1, g_0_xxxxxzz_0_xxyzzzz_0, g_0_xxxxxzz_0_xxyzzzz_1, g_0_xxxxxzz_0_xxzzzz_1, g_0_xxxxxzz_0_xxzzzzz_0, g_0_xxxxxzz_0_xxzzzzz_1, g_0_xxxxxzz_0_xyyyyyz_0, g_0_xxxxxzz_0_xyyyyyz_1, g_0_xxxxxzz_0_xyyyyz_1, g_0_xxxxxzz_0_xyyyyzz_0, g_0_xxxxxzz_0_xyyyyzz_1, g_0_xxxxxzz_0_xyyyzz_1, g_0_xxxxxzz_0_xyyyzzz_0, g_0_xxxxxzz_0_xyyyzzz_1, g_0_xxxxxzz_0_xyyzzz_1, g_0_xxxxxzz_0_xyyzzzz_0, g_0_xxxxxzz_0_xyyzzzz_1, g_0_xxxxxzz_0_xyzzzz_1, g_0_xxxxxzz_0_xyzzzzz_0, g_0_xxxxxzz_0_xyzzzzz_1, g_0_xxxxxzz_0_xzzzzz_1, g_0_xxxxxzz_0_xzzzzzz_0, g_0_xxxxxzz_0_xzzzzzz_1, g_0_xxxxxzz_0_yyyyyyy_0, g_0_xxxxxzz_0_yyyyyyy_1, g_0_xxxxxzz_0_yyyyyyz_0, g_0_xxxxxzz_0_yyyyyyz_1, g_0_xxxxxzz_0_yyyyyz_1, g_0_xxxxxzz_0_yyyyyzz_0, g_0_xxxxxzz_0_yyyyyzz_1, g_0_xxxxxzz_0_yyyyzz_1, g_0_xxxxxzz_0_yyyyzzz_0, g_0_xxxxxzz_0_yyyyzzz_1, g_0_xxxxxzz_0_yyyzzz_1, g_0_xxxxxzz_0_yyyzzzz_0, g_0_xxxxxzz_0_yyyzzzz_1, g_0_xxxxxzz_0_yyzzzz_1, g_0_xxxxxzz_0_yyzzzzz_0, g_0_xxxxxzz_0_yyzzzzz_1, g_0_xxxxxzz_0_yzzzzz_1, g_0_xxxxxzz_0_yzzzzzz_0, g_0_xxxxxzz_0_yzzzzzz_1, g_0_xxxxxzz_0_zzzzzz_1, g_0_xxxxxzz_0_zzzzzzz_0, g_0_xxxxxzz_0_zzzzzzz_1, g_0_xxxxzz_0_xxxxxxz_0, g_0_xxxxzz_0_xxxxxxz_1, g_0_xxxxzz_0_xxxxxyz_0, g_0_xxxxzz_0_xxxxxyz_1, g_0_xxxxzz_0_xxxxxzz_0, g_0_xxxxzz_0_xxxxxzz_1, g_0_xxxxzz_0_xxxxyyz_0, g_0_xxxxzz_0_xxxxyyz_1, g_0_xxxxzz_0_xxxxyzz_0, g_0_xxxxzz_0_xxxxyzz_1, g_0_xxxxzz_0_xxxxzzz_0, g_0_xxxxzz_0_xxxxzzz_1, g_0_xxxxzz_0_xxxyyyz_0, g_0_xxxxzz_0_xxxyyyz_1, g_0_xxxxzz_0_xxxyyzz_0, g_0_xxxxzz_0_xxxyyzz_1, g_0_xxxxzz_0_xxxyzzz_0, g_0_xxxxzz_0_xxxyzzz_1, g_0_xxxxzz_0_xxxzzzz_0, g_0_xxxxzz_0_xxxzzzz_1, g_0_xxxxzz_0_xxyyyyz_0, g_0_xxxxzz_0_xxyyyyz_1, g_0_xxxxzz_0_xxyyyzz_0, g_0_xxxxzz_0_xxyyyzz_1, g_0_xxxxzz_0_xxyyzzz_0, g_0_xxxxzz_0_xxyyzzz_1, g_0_xxxxzz_0_xxyzzzz_0, g_0_xxxxzz_0_xxyzzzz_1, g_0_xxxxzz_0_xxzzzzz_0, g_0_xxxxzz_0_xxzzzzz_1, g_0_xxxxzz_0_xyyyyyz_0, g_0_xxxxzz_0_xyyyyyz_1, g_0_xxxxzz_0_xyyyyzz_0, g_0_xxxxzz_0_xyyyyzz_1, g_0_xxxxzz_0_xyyyzzz_0, g_0_xxxxzz_0_xyyyzzz_1, g_0_xxxxzz_0_xyyzzzz_0, g_0_xxxxzz_0_xyyzzzz_1, g_0_xxxxzz_0_xyzzzzz_0, g_0_xxxxzz_0_xyzzzzz_1, g_0_xxxxzz_0_xzzzzzz_0, g_0_xxxxzz_0_xzzzzzz_1, g_0_xxxxzz_0_yyyyyyy_0, g_0_xxxxzz_0_yyyyyyy_1, g_0_xxxxzz_0_yyyyyyz_0, g_0_xxxxzz_0_yyyyyyz_1, g_0_xxxxzz_0_yyyyyzz_0, g_0_xxxxzz_0_yyyyyzz_1, g_0_xxxxzz_0_yyyyzzz_0, g_0_xxxxzz_0_yyyyzzz_1, g_0_xxxxzz_0_yyyzzzz_0, g_0_xxxxzz_0_yyyzzzz_1, g_0_xxxxzz_0_yyzzzzz_0, g_0_xxxxzz_0_yyzzzzz_1, g_0_xxxxzz_0_yzzzzzz_0, g_0_xxxxzz_0_yzzzzzz_1, g_0_xxxxzz_0_zzzzzzz_0, g_0_xxxxzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxzz_0_xxxxxxx_0[i] = g_0_xxxxxx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxxxx_0[i] * pb_z + g_0_xxxxxxz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxxxy_0[i] = g_0_xxxxxx_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxxxy_0[i] * pb_z + g_0_xxxxxxz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxxxz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxxxzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxxxz_0[i] * pb_x + g_0_xxxxxzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxxxyy_0[i] = g_0_xxxxxx_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxxyy_0[i] * pb_z + g_0_xxxxxxz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxxyz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxxyz_0[i] * pb_x + g_0_xxxxxzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxxxzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxxzz_0[i] * pb_x + g_0_xxxxxzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxxyyy_0[i] = g_0_xxxxxx_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxyyy_0[i] * pb_z + g_0_xxxxxxz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxyyz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxyyz_0[i] * pb_x + g_0_xxxxxzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxxyzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxyzz_0[i] * pb_x + g_0_xxxxxzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxxzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxzzz_0[i] * pb_x + g_0_xxxxxzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxyyyy_0[i] = g_0_xxxxxx_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxyyyy_0[i] * pb_z + g_0_xxxxxxz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxyyyz_0[i] = 5.0 * g_0_xxxxzz_0_xxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyyyz_0[i] * pb_x + g_0_xxxxxzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxyyzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyyzz_0[i] * pb_x + g_0_xxxxxzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxyzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyzzz_0[i] * pb_x + g_0_xxxxxzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyyyyy_0[i] = g_0_xxxxxx_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxyyyyy_0[i] * pb_z + g_0_xxxxxxz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxyyyyz_0[i] = 5.0 * g_0_xxxxzz_0_xxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyyyz_0[i] * pb_x + g_0_xxxxxzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyyyzz_0[i] = 5.0 * g_0_xxxxzz_0_xxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyyzz_0[i] * pb_x + g_0_xxxxxzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyyzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyzzz_0[i] * pb_x + g_0_xxxxxzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyyyyy_0[i] = g_0_xxxxxx_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xyyyyyy_0[i] * pb_z + g_0_xxxxxxz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxzz_0_xyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyyyz_0[i] * pb_x + g_0_xxxxxzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyyyzz_0[i] = 5.0 * g_0_xxxxzz_0_xyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyyzz_0[i] * pb_x + g_0_xxxxxzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyyzzz_0[i] = 5.0 * g_0_xxxxzz_0_xyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyzzz_0[i] * pb_x + g_0_xxxxxzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xzzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xzzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyyyy_0[i] = 5.0 * g_0_xxxxzz_0_yyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyyyyy_0[i] * pb_x + g_0_xxxxxzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyyyz_0[i] = 5.0 * g_0_xxxxzz_0_yyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyyyyz_0[i] * pb_x + g_0_xxxxxzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxzz_0_yyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyyyzz_0[i] * pb_x + g_0_xxxxxzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyzzz_0[i] = 5.0 * g_0_xxxxzz_0_yyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyyzzz_0[i] * pb_x + g_0_xxxxxzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyzzzz_0[i] = 5.0 * g_0_xxxxzz_0_yyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyzzzz_0[i] * pb_x + g_0_xxxxxzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_yyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yzzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_yzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yzzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_zzzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_zzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zzzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 216-252 components of targeted buffer : SLSK

    auto g_0_xxxxxyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 216);

    auto g_0_xxxxxyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 217);

    auto g_0_xxxxxyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 218);

    auto g_0_xxxxxyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 219);

    auto g_0_xxxxxyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 220);

    auto g_0_xxxxxyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 221);

    auto g_0_xxxxxyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 222);

    auto g_0_xxxxxyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 223);

    auto g_0_xxxxxyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 224);

    auto g_0_xxxxxyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 225);

    auto g_0_xxxxxyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 226);

    auto g_0_xxxxxyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 227);

    auto g_0_xxxxxyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 228);

    auto g_0_xxxxxyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 229);

    auto g_0_xxxxxyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 230);

    auto g_0_xxxxxyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 231);

    auto g_0_xxxxxyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 232);

    auto g_0_xxxxxyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 233);

    auto g_0_xxxxxyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 234);

    auto g_0_xxxxxyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 235);

    auto g_0_xxxxxyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 236);

    auto g_0_xxxxxyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 237);

    auto g_0_xxxxxyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 238);

    auto g_0_xxxxxyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 239);

    auto g_0_xxxxxyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 240);

    auto g_0_xxxxxyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 241);

    auto g_0_xxxxxyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 242);

    auto g_0_xxxxxyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 243);

    auto g_0_xxxxxyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 244);

    auto g_0_xxxxxyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 245);

    auto g_0_xxxxxyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 246);

    auto g_0_xxxxxyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 247);

    auto g_0_xxxxxyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 248);

    auto g_0_xxxxxyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 249);

    auto g_0_xxxxxyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 250);

    auto g_0_xxxxxyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 251);

    #pragma omp simd aligned(g_0_xxxxxy_0_xxxxxxx_0, g_0_xxxxxy_0_xxxxxxx_1, g_0_xxxxxy_0_xxxxxxz_0, g_0_xxxxxy_0_xxxxxxz_1, g_0_xxxxxy_0_xxxxxzz_0, g_0_xxxxxy_0_xxxxxzz_1, g_0_xxxxxy_0_xxxxzzz_0, g_0_xxxxxy_0_xxxxzzz_1, g_0_xxxxxy_0_xxxzzzz_0, g_0_xxxxxy_0_xxxzzzz_1, g_0_xxxxxy_0_xxzzzzz_0, g_0_xxxxxy_0_xxzzzzz_1, g_0_xxxxxy_0_xzzzzzz_0, g_0_xxxxxy_0_xzzzzzz_1, g_0_xxxxxyy_0_xxxxxxx_0, g_0_xxxxxyy_0_xxxxxxx_1, g_0_xxxxxyy_0_xxxxxxz_0, g_0_xxxxxyy_0_xxxxxxz_1, g_0_xxxxxyy_0_xxxxxzz_0, g_0_xxxxxyy_0_xxxxxzz_1, g_0_xxxxxyy_0_xxxxzzz_0, g_0_xxxxxyy_0_xxxxzzz_1, g_0_xxxxxyy_0_xxxzzzz_0, g_0_xxxxxyy_0_xxxzzzz_1, g_0_xxxxxyy_0_xxzzzzz_0, g_0_xxxxxyy_0_xxzzzzz_1, g_0_xxxxxyy_0_xzzzzzz_0, g_0_xxxxxyy_0_xzzzzzz_1, g_0_xxxxxyyy_0_xxxxxxx_0, g_0_xxxxxyyy_0_xxxxxxy_0, g_0_xxxxxyyy_0_xxxxxxz_0, g_0_xxxxxyyy_0_xxxxxyy_0, g_0_xxxxxyyy_0_xxxxxyz_0, g_0_xxxxxyyy_0_xxxxxzz_0, g_0_xxxxxyyy_0_xxxxyyy_0, g_0_xxxxxyyy_0_xxxxyyz_0, g_0_xxxxxyyy_0_xxxxyzz_0, g_0_xxxxxyyy_0_xxxxzzz_0, g_0_xxxxxyyy_0_xxxyyyy_0, g_0_xxxxxyyy_0_xxxyyyz_0, g_0_xxxxxyyy_0_xxxyyzz_0, g_0_xxxxxyyy_0_xxxyzzz_0, g_0_xxxxxyyy_0_xxxzzzz_0, g_0_xxxxxyyy_0_xxyyyyy_0, g_0_xxxxxyyy_0_xxyyyyz_0, g_0_xxxxxyyy_0_xxyyyzz_0, g_0_xxxxxyyy_0_xxyyzzz_0, g_0_xxxxxyyy_0_xxyzzzz_0, g_0_xxxxxyyy_0_xxzzzzz_0, g_0_xxxxxyyy_0_xyyyyyy_0, g_0_xxxxxyyy_0_xyyyyyz_0, g_0_xxxxxyyy_0_xyyyyzz_0, g_0_xxxxxyyy_0_xyyyzzz_0, g_0_xxxxxyyy_0_xyyzzzz_0, g_0_xxxxxyyy_0_xyzzzzz_0, g_0_xxxxxyyy_0_xzzzzzz_0, g_0_xxxxxyyy_0_yyyyyyy_0, g_0_xxxxxyyy_0_yyyyyyz_0, g_0_xxxxxyyy_0_yyyyyzz_0, g_0_xxxxxyyy_0_yyyyzzz_0, g_0_xxxxxyyy_0_yyyzzzz_0, g_0_xxxxxyyy_0_yyzzzzz_0, g_0_xxxxxyyy_0_yzzzzzz_0, g_0_xxxxxyyy_0_zzzzzzz_0, g_0_xxxxyyy_0_xxxxxxy_0, g_0_xxxxyyy_0_xxxxxxy_1, g_0_xxxxyyy_0_xxxxxy_1, g_0_xxxxyyy_0_xxxxxyy_0, g_0_xxxxyyy_0_xxxxxyy_1, g_0_xxxxyyy_0_xxxxxyz_0, g_0_xxxxyyy_0_xxxxxyz_1, g_0_xxxxyyy_0_xxxxyy_1, g_0_xxxxyyy_0_xxxxyyy_0, g_0_xxxxyyy_0_xxxxyyy_1, g_0_xxxxyyy_0_xxxxyyz_0, g_0_xxxxyyy_0_xxxxyyz_1, g_0_xxxxyyy_0_xxxxyz_1, g_0_xxxxyyy_0_xxxxyzz_0, g_0_xxxxyyy_0_xxxxyzz_1, g_0_xxxxyyy_0_xxxyyy_1, g_0_xxxxyyy_0_xxxyyyy_0, g_0_xxxxyyy_0_xxxyyyy_1, g_0_xxxxyyy_0_xxxyyyz_0, g_0_xxxxyyy_0_xxxyyyz_1, g_0_xxxxyyy_0_xxxyyz_1, g_0_xxxxyyy_0_xxxyyzz_0, g_0_xxxxyyy_0_xxxyyzz_1, g_0_xxxxyyy_0_xxxyzz_1, g_0_xxxxyyy_0_xxxyzzz_0, g_0_xxxxyyy_0_xxxyzzz_1, g_0_xxxxyyy_0_xxyyyy_1, g_0_xxxxyyy_0_xxyyyyy_0, g_0_xxxxyyy_0_xxyyyyy_1, g_0_xxxxyyy_0_xxyyyyz_0, g_0_xxxxyyy_0_xxyyyyz_1, g_0_xxxxyyy_0_xxyyyz_1, g_0_xxxxyyy_0_xxyyyzz_0, g_0_xxxxyyy_0_xxyyyzz_1, g_0_xxxxyyy_0_xxyyzz_1, g_0_xxxxyyy_0_xxyyzzz_0, g_0_xxxxyyy_0_xxyyzzz_1, g_0_xxxxyyy_0_xxyzzz_1, g_0_xxxxyyy_0_xxyzzzz_0, g_0_xxxxyyy_0_xxyzzzz_1, g_0_xxxxyyy_0_xyyyyy_1, g_0_xxxxyyy_0_xyyyyyy_0, g_0_xxxxyyy_0_xyyyyyy_1, g_0_xxxxyyy_0_xyyyyyz_0, g_0_xxxxyyy_0_xyyyyyz_1, g_0_xxxxyyy_0_xyyyyz_1, g_0_xxxxyyy_0_xyyyyzz_0, g_0_xxxxyyy_0_xyyyyzz_1, g_0_xxxxyyy_0_xyyyzz_1, g_0_xxxxyyy_0_xyyyzzz_0, g_0_xxxxyyy_0_xyyyzzz_1, g_0_xxxxyyy_0_xyyzzz_1, g_0_xxxxyyy_0_xyyzzzz_0, g_0_xxxxyyy_0_xyyzzzz_1, g_0_xxxxyyy_0_xyzzzz_1, g_0_xxxxyyy_0_xyzzzzz_0, g_0_xxxxyyy_0_xyzzzzz_1, g_0_xxxxyyy_0_yyyyyy_1, g_0_xxxxyyy_0_yyyyyyy_0, g_0_xxxxyyy_0_yyyyyyy_1, g_0_xxxxyyy_0_yyyyyyz_0, g_0_xxxxyyy_0_yyyyyyz_1, g_0_xxxxyyy_0_yyyyyz_1, g_0_xxxxyyy_0_yyyyyzz_0, g_0_xxxxyyy_0_yyyyyzz_1, g_0_xxxxyyy_0_yyyyzz_1, g_0_xxxxyyy_0_yyyyzzz_0, g_0_xxxxyyy_0_yyyyzzz_1, g_0_xxxxyyy_0_yyyzzz_1, g_0_xxxxyyy_0_yyyzzzz_0, g_0_xxxxyyy_0_yyyzzzz_1, g_0_xxxxyyy_0_yyzzzz_1, g_0_xxxxyyy_0_yyzzzzz_0, g_0_xxxxyyy_0_yyzzzzz_1, g_0_xxxxyyy_0_yzzzzz_1, g_0_xxxxyyy_0_yzzzzzz_0, g_0_xxxxyyy_0_yzzzzzz_1, g_0_xxxxyyy_0_zzzzzzz_0, g_0_xxxxyyy_0_zzzzzzz_1, g_0_xxxyyy_0_xxxxxxy_0, g_0_xxxyyy_0_xxxxxxy_1, g_0_xxxyyy_0_xxxxxyy_0, g_0_xxxyyy_0_xxxxxyy_1, g_0_xxxyyy_0_xxxxxyz_0, g_0_xxxyyy_0_xxxxxyz_1, g_0_xxxyyy_0_xxxxyyy_0, g_0_xxxyyy_0_xxxxyyy_1, g_0_xxxyyy_0_xxxxyyz_0, g_0_xxxyyy_0_xxxxyyz_1, g_0_xxxyyy_0_xxxxyzz_0, g_0_xxxyyy_0_xxxxyzz_1, g_0_xxxyyy_0_xxxyyyy_0, g_0_xxxyyy_0_xxxyyyy_1, g_0_xxxyyy_0_xxxyyyz_0, g_0_xxxyyy_0_xxxyyyz_1, g_0_xxxyyy_0_xxxyyzz_0, g_0_xxxyyy_0_xxxyyzz_1, g_0_xxxyyy_0_xxxyzzz_0, g_0_xxxyyy_0_xxxyzzz_1, g_0_xxxyyy_0_xxyyyyy_0, g_0_xxxyyy_0_xxyyyyy_1, g_0_xxxyyy_0_xxyyyyz_0, g_0_xxxyyy_0_xxyyyyz_1, g_0_xxxyyy_0_xxyyyzz_0, g_0_xxxyyy_0_xxyyyzz_1, g_0_xxxyyy_0_xxyyzzz_0, g_0_xxxyyy_0_xxyyzzz_1, g_0_xxxyyy_0_xxyzzzz_0, g_0_xxxyyy_0_xxyzzzz_1, g_0_xxxyyy_0_xyyyyyy_0, g_0_xxxyyy_0_xyyyyyy_1, g_0_xxxyyy_0_xyyyyyz_0, g_0_xxxyyy_0_xyyyyyz_1, g_0_xxxyyy_0_xyyyyzz_0, g_0_xxxyyy_0_xyyyyzz_1, g_0_xxxyyy_0_xyyyzzz_0, g_0_xxxyyy_0_xyyyzzz_1, g_0_xxxyyy_0_xyyzzzz_0, g_0_xxxyyy_0_xyyzzzz_1, g_0_xxxyyy_0_xyzzzzz_0, g_0_xxxyyy_0_xyzzzzz_1, g_0_xxxyyy_0_yyyyyyy_0, g_0_xxxyyy_0_yyyyyyy_1, g_0_xxxyyy_0_yyyyyyz_0, g_0_xxxyyy_0_yyyyyyz_1, g_0_xxxyyy_0_yyyyyzz_0, g_0_xxxyyy_0_yyyyyzz_1, g_0_xxxyyy_0_yyyyzzz_0, g_0_xxxyyy_0_yyyyzzz_1, g_0_xxxyyy_0_yyyzzzz_0, g_0_xxxyyy_0_yyyzzzz_1, g_0_xxxyyy_0_yyzzzzz_0, g_0_xxxyyy_0_yyzzzzz_1, g_0_xxxyyy_0_yzzzzzz_0, g_0_xxxyyy_0_yzzzzzz_1, g_0_xxxyyy_0_zzzzzzz_0, g_0_xxxyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyyy_0_xxxxxxx_0[i] = 2.0 * g_0_xxxxxy_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxxxxx_0[i] * pb_y + g_0_xxxxxyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxxxxy_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxxyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxxxy_0[i] * pb_x + g_0_xxxxyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxxxz_0[i] = 2.0 * g_0_xxxxxy_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxxxxz_0[i] * pb_y + g_0_xxxxxyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxxxyy_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxxyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxxyy_0[i] * pb_x + g_0_xxxxyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxxyz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxxyz_0[i] * pb_x + g_0_xxxxyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxxzz_0[i] = 2.0 * g_0_xxxxxy_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxxxzz_0[i] * pb_y + g_0_xxxxxyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxxyyy_0[i] = 4.0 * g_0_xxxyyy_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxyyy_0[i] * pb_x + g_0_xxxxyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxyyz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxyyz_0[i] * pb_x + g_0_xxxxyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxyzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxyzz_0[i] * pb_x + g_0_xxxxyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxzzz_0[i] = 2.0 * g_0_xxxxxy_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxxzzz_0[i] * pb_y + g_0_xxxxxyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxyyyy_0[i] = 4.0 * g_0_xxxyyy_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyyyy_0[i] * pb_x + g_0_xxxxyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxyyyz_0[i] = 4.0 * g_0_xxxyyy_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyyyz_0[i] * pb_x + g_0_xxxxyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxyyzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyyzz_0[i] * pb_x + g_0_xxxxyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxyzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyzzz_0[i] * pb_x + g_0_xxxxyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxzzzz_0[i] = 2.0 * g_0_xxxxxy_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxzzzz_0[i] * pb_y + g_0_xxxxxyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxyyyyy_0[i] = 4.0 * g_0_xxxyyy_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyyyy_0[i] * pb_x + g_0_xxxxyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyyyyz_0[i] = 4.0 * g_0_xxxyyy_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyyyz_0[i] * pb_x + g_0_xxxxyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyyyzz_0[i] = 4.0 * g_0_xxxyyy_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyyzz_0[i] * pb_x + g_0_xxxxyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyyzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyzzz_0[i] * pb_x + g_0_xxxxyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyzzzz_0[i] * pb_x + g_0_xxxxyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxzzzzz_0[i] = 2.0 * g_0_xxxxxy_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxzzzzz_0[i] * pb_y + g_0_xxxxxyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xyyyyyy_0[i] = 4.0 * g_0_xxxyyy_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyyyy_0[i] * pb_x + g_0_xxxxyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyyyyz_0[i] = 4.0 * g_0_xxxyyy_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyyyz_0[i] * pb_x + g_0_xxxxyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyyyzz_0[i] = 4.0 * g_0_xxxyyy_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyyzz_0[i] * pb_x + g_0_xxxxyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyyzzz_0[i] = 4.0 * g_0_xxxyyy_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyzzz_0[i] * pb_x + g_0_xxxxyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyzzzz_0[i] * pb_x + g_0_xxxxyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyzzzzz_0[i] * pb_x + g_0_xxxxyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xzzzzzz_0[i] = 2.0 * g_0_xxxxxy_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xzzzzzz_0[i] * pb_y + g_0_xxxxxyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_yyyyyyy_0[i] = 4.0 * g_0_xxxyyy_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyyyy_0[i] * pb_x + g_0_xxxxyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyyyyz_0[i] = 4.0 * g_0_xxxyyy_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyyyz_0[i] * pb_x + g_0_xxxxyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyyyzz_0[i] = 4.0 * g_0_xxxyyy_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyyzz_0[i] * pb_x + g_0_xxxxyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyyzzz_0[i] = 4.0 * g_0_xxxyyy_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyzzz_0[i] * pb_x + g_0_xxxxyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyzzzz_0[i] * pb_x + g_0_xxxxyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyzzzzz_0[i] * pb_x + g_0_xxxxyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yzzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yzzzzzz_0[i] * pb_x + g_0_xxxxyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_zzzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_zzzzzzz_0[i] * pb_x + g_0_xxxxyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 252-288 components of targeted buffer : SLSK

    auto g_0_xxxxxyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 252);

    auto g_0_xxxxxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 253);

    auto g_0_xxxxxyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 254);

    auto g_0_xxxxxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 255);

    auto g_0_xxxxxyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 256);

    auto g_0_xxxxxyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 257);

    auto g_0_xxxxxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 258);

    auto g_0_xxxxxyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 259);

    auto g_0_xxxxxyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 260);

    auto g_0_xxxxxyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 261);

    auto g_0_xxxxxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 262);

    auto g_0_xxxxxyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 263);

    auto g_0_xxxxxyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 264);

    auto g_0_xxxxxyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 265);

    auto g_0_xxxxxyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 266);

    auto g_0_xxxxxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 267);

    auto g_0_xxxxxyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 268);

    auto g_0_xxxxxyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 269);

    auto g_0_xxxxxyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 270);

    auto g_0_xxxxxyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 271);

    auto g_0_xxxxxyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 272);

    auto g_0_xxxxxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 273);

    auto g_0_xxxxxyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 274);

    auto g_0_xxxxxyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 275);

    auto g_0_xxxxxyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 276);

    auto g_0_xxxxxyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 277);

    auto g_0_xxxxxyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 278);

    auto g_0_xxxxxyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 279);

    auto g_0_xxxxxyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 280);

    auto g_0_xxxxxyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 281);

    auto g_0_xxxxxyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 282);

    auto g_0_xxxxxyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 283);

    auto g_0_xxxxxyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 284);

    auto g_0_xxxxxyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 285);

    auto g_0_xxxxxyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 286);

    auto g_0_xxxxxyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 287);

    #pragma omp simd aligned(g_0_xxxxxyy_0_xxxxxx_1, g_0_xxxxxyy_0_xxxxxxx_0, g_0_xxxxxyy_0_xxxxxxx_1, g_0_xxxxxyy_0_xxxxxxy_0, g_0_xxxxxyy_0_xxxxxxy_1, g_0_xxxxxyy_0_xxxxxxz_0, g_0_xxxxxyy_0_xxxxxxz_1, g_0_xxxxxyy_0_xxxxxy_1, g_0_xxxxxyy_0_xxxxxyy_0, g_0_xxxxxyy_0_xxxxxyy_1, g_0_xxxxxyy_0_xxxxxyz_0, g_0_xxxxxyy_0_xxxxxyz_1, g_0_xxxxxyy_0_xxxxxz_1, g_0_xxxxxyy_0_xxxxxzz_0, g_0_xxxxxyy_0_xxxxxzz_1, g_0_xxxxxyy_0_xxxxyy_1, g_0_xxxxxyy_0_xxxxyyy_0, g_0_xxxxxyy_0_xxxxyyy_1, g_0_xxxxxyy_0_xxxxyyz_0, g_0_xxxxxyy_0_xxxxyyz_1, g_0_xxxxxyy_0_xxxxyz_1, g_0_xxxxxyy_0_xxxxyzz_0, g_0_xxxxxyy_0_xxxxyzz_1, g_0_xxxxxyy_0_xxxxzz_1, g_0_xxxxxyy_0_xxxxzzz_0, g_0_xxxxxyy_0_xxxxzzz_1, g_0_xxxxxyy_0_xxxyyy_1, g_0_xxxxxyy_0_xxxyyyy_0, g_0_xxxxxyy_0_xxxyyyy_1, g_0_xxxxxyy_0_xxxyyyz_0, g_0_xxxxxyy_0_xxxyyyz_1, g_0_xxxxxyy_0_xxxyyz_1, g_0_xxxxxyy_0_xxxyyzz_0, g_0_xxxxxyy_0_xxxyyzz_1, g_0_xxxxxyy_0_xxxyzz_1, g_0_xxxxxyy_0_xxxyzzz_0, g_0_xxxxxyy_0_xxxyzzz_1, g_0_xxxxxyy_0_xxxzzz_1, g_0_xxxxxyy_0_xxxzzzz_0, g_0_xxxxxyy_0_xxxzzzz_1, g_0_xxxxxyy_0_xxyyyy_1, g_0_xxxxxyy_0_xxyyyyy_0, g_0_xxxxxyy_0_xxyyyyy_1, g_0_xxxxxyy_0_xxyyyyz_0, g_0_xxxxxyy_0_xxyyyyz_1, g_0_xxxxxyy_0_xxyyyz_1, g_0_xxxxxyy_0_xxyyyzz_0, g_0_xxxxxyy_0_xxyyyzz_1, g_0_xxxxxyy_0_xxyyzz_1, g_0_xxxxxyy_0_xxyyzzz_0, g_0_xxxxxyy_0_xxyyzzz_1, g_0_xxxxxyy_0_xxyzzz_1, g_0_xxxxxyy_0_xxyzzzz_0, g_0_xxxxxyy_0_xxyzzzz_1, g_0_xxxxxyy_0_xxzzzz_1, g_0_xxxxxyy_0_xxzzzzz_0, g_0_xxxxxyy_0_xxzzzzz_1, g_0_xxxxxyy_0_xyyyyy_1, g_0_xxxxxyy_0_xyyyyyy_0, g_0_xxxxxyy_0_xyyyyyy_1, g_0_xxxxxyy_0_xyyyyyz_0, g_0_xxxxxyy_0_xyyyyyz_1, g_0_xxxxxyy_0_xyyyyz_1, g_0_xxxxxyy_0_xyyyyzz_0, g_0_xxxxxyy_0_xyyyyzz_1, g_0_xxxxxyy_0_xyyyzz_1, g_0_xxxxxyy_0_xyyyzzz_0, g_0_xxxxxyy_0_xyyyzzz_1, g_0_xxxxxyy_0_xyyzzz_1, g_0_xxxxxyy_0_xyyzzzz_0, g_0_xxxxxyy_0_xyyzzzz_1, g_0_xxxxxyy_0_xyzzzz_1, g_0_xxxxxyy_0_xyzzzzz_0, g_0_xxxxxyy_0_xyzzzzz_1, g_0_xxxxxyy_0_xzzzzz_1, g_0_xxxxxyy_0_xzzzzzz_0, g_0_xxxxxyy_0_xzzzzzz_1, g_0_xxxxxyy_0_yyyyyy_1, g_0_xxxxxyy_0_yyyyyyy_0, g_0_xxxxxyy_0_yyyyyyy_1, g_0_xxxxxyy_0_yyyyyyz_0, g_0_xxxxxyy_0_yyyyyyz_1, g_0_xxxxxyy_0_yyyyyz_1, g_0_xxxxxyy_0_yyyyyzz_0, g_0_xxxxxyy_0_yyyyyzz_1, g_0_xxxxxyy_0_yyyyzz_1, g_0_xxxxxyy_0_yyyyzzz_0, g_0_xxxxxyy_0_yyyyzzz_1, g_0_xxxxxyy_0_yyyzzz_1, g_0_xxxxxyy_0_yyyzzzz_0, g_0_xxxxxyy_0_yyyzzzz_1, g_0_xxxxxyy_0_yyzzzz_1, g_0_xxxxxyy_0_yyzzzzz_0, g_0_xxxxxyy_0_yyzzzzz_1, g_0_xxxxxyy_0_yzzzzz_1, g_0_xxxxxyy_0_yzzzzzz_0, g_0_xxxxxyy_0_yzzzzzz_1, g_0_xxxxxyy_0_zzzzzz_1, g_0_xxxxxyy_0_zzzzzzz_0, g_0_xxxxxyy_0_zzzzzzz_1, g_0_xxxxxyyz_0_xxxxxxx_0, g_0_xxxxxyyz_0_xxxxxxy_0, g_0_xxxxxyyz_0_xxxxxxz_0, g_0_xxxxxyyz_0_xxxxxyy_0, g_0_xxxxxyyz_0_xxxxxyz_0, g_0_xxxxxyyz_0_xxxxxzz_0, g_0_xxxxxyyz_0_xxxxyyy_0, g_0_xxxxxyyz_0_xxxxyyz_0, g_0_xxxxxyyz_0_xxxxyzz_0, g_0_xxxxxyyz_0_xxxxzzz_0, g_0_xxxxxyyz_0_xxxyyyy_0, g_0_xxxxxyyz_0_xxxyyyz_0, g_0_xxxxxyyz_0_xxxyyzz_0, g_0_xxxxxyyz_0_xxxyzzz_0, g_0_xxxxxyyz_0_xxxzzzz_0, g_0_xxxxxyyz_0_xxyyyyy_0, g_0_xxxxxyyz_0_xxyyyyz_0, g_0_xxxxxyyz_0_xxyyyzz_0, g_0_xxxxxyyz_0_xxyyzzz_0, g_0_xxxxxyyz_0_xxyzzzz_0, g_0_xxxxxyyz_0_xxzzzzz_0, g_0_xxxxxyyz_0_xyyyyyy_0, g_0_xxxxxyyz_0_xyyyyyz_0, g_0_xxxxxyyz_0_xyyyyzz_0, g_0_xxxxxyyz_0_xyyyzzz_0, g_0_xxxxxyyz_0_xyyzzzz_0, g_0_xxxxxyyz_0_xyzzzzz_0, g_0_xxxxxyyz_0_xzzzzzz_0, g_0_xxxxxyyz_0_yyyyyyy_0, g_0_xxxxxyyz_0_yyyyyyz_0, g_0_xxxxxyyz_0_yyyyyzz_0, g_0_xxxxxyyz_0_yyyyzzz_0, g_0_xxxxxyyz_0_yyyzzzz_0, g_0_xxxxxyyz_0_yyzzzzz_0, g_0_xxxxxyyz_0_yzzzzzz_0, g_0_xxxxxyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyyz_0_xxxxxxx_0[i] = g_0_xxxxxyy_0_xxxxxxx_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxxxy_0[i] = g_0_xxxxxyy_0_xxxxxxy_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxxxz_0[i] = g_0_xxxxxyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxxxz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxxyy_0[i] = g_0_xxxxxyy_0_xxxxxyy_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxxyz_0[i] = g_0_xxxxxyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxxyz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxxxyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxxzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxyyy_0[i] = g_0_xxxxxyy_0_xxxxyyy_0[i] * pb_z + g_0_xxxxxyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxyyz_0[i] = g_0_xxxxxyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxyyz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxxxyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxyzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxxxyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyyyy_0[i] = g_0_xxxxxyy_0_xxxyyyy_0[i] * pb_z + g_0_xxxxxyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyyyz_0[i] = g_0_xxxxxyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyyyz_0[i] * pb_z + g_0_xxxxxyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxxyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyyzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxxxyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxxxyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyyyy_0[i] = g_0_xxxxxyy_0_xxyyyyy_0[i] * pb_z + g_0_xxxxxyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyyyz_0[i] = g_0_xxxxxyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyyyz_0[i] * pb_z + g_0_xxxxxyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxxxyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyyzz_0[i] * pb_z + g_0_xxxxxyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxxxyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxxxyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxxxyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyyyy_0[i] = g_0_xxxxxyy_0_xyyyyyy_0[i] * pb_z + g_0_xxxxxyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyyyz_0[i] = g_0_xxxxxyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyyyz_0[i] * pb_z + g_0_xxxxxyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxxxyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyyzz_0[i] * pb_z + g_0_xxxxxyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxxyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyzzz_0[i] * pb_z + g_0_xxxxxyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxxxyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxxxyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxxxyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xzzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyyyy_0[i] = g_0_xxxxxyy_0_yyyyyyy_0[i] * pb_z + g_0_xxxxxyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyyyz_0[i] = g_0_xxxxxyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyyyyz_0[i] * pb_z + g_0_xxxxxyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxxxyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyyyzz_0[i] * pb_z + g_0_xxxxxyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxxxyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyyzzz_0[i] * pb_z + g_0_xxxxxyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxxxyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyzzzz_0[i] * pb_z + g_0_xxxxxyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxxxyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxxxyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yzzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxxxyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_zzzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 288-324 components of targeted buffer : SLSK

    auto g_0_xxxxxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 288);

    auto g_0_xxxxxyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 289);

    auto g_0_xxxxxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 290);

    auto g_0_xxxxxyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 291);

    auto g_0_xxxxxyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 292);

    auto g_0_xxxxxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 293);

    auto g_0_xxxxxyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 294);

    auto g_0_xxxxxyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 295);

    auto g_0_xxxxxyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 296);

    auto g_0_xxxxxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 297);

    auto g_0_xxxxxyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 298);

    auto g_0_xxxxxyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 299);

    auto g_0_xxxxxyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 300);

    auto g_0_xxxxxyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 301);

    auto g_0_xxxxxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 302);

    auto g_0_xxxxxyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 303);

    auto g_0_xxxxxyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 304);

    auto g_0_xxxxxyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 305);

    auto g_0_xxxxxyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 306);

    auto g_0_xxxxxyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 307);

    auto g_0_xxxxxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 308);

    auto g_0_xxxxxyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 309);

    auto g_0_xxxxxyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 310);

    auto g_0_xxxxxyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 311);

    auto g_0_xxxxxyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 312);

    auto g_0_xxxxxyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 313);

    auto g_0_xxxxxyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 314);

    auto g_0_xxxxxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 315);

    auto g_0_xxxxxyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 316);

    auto g_0_xxxxxyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 317);

    auto g_0_xxxxxyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 318);

    auto g_0_xxxxxyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 319);

    auto g_0_xxxxxyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 320);

    auto g_0_xxxxxyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 321);

    auto g_0_xxxxxyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 322);

    auto g_0_xxxxxyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 323);

    #pragma omp simd aligned(g_0_xxxxxyzz_0_xxxxxxx_0, g_0_xxxxxyzz_0_xxxxxxy_0, g_0_xxxxxyzz_0_xxxxxxz_0, g_0_xxxxxyzz_0_xxxxxyy_0, g_0_xxxxxyzz_0_xxxxxyz_0, g_0_xxxxxyzz_0_xxxxxzz_0, g_0_xxxxxyzz_0_xxxxyyy_0, g_0_xxxxxyzz_0_xxxxyyz_0, g_0_xxxxxyzz_0_xxxxyzz_0, g_0_xxxxxyzz_0_xxxxzzz_0, g_0_xxxxxyzz_0_xxxyyyy_0, g_0_xxxxxyzz_0_xxxyyyz_0, g_0_xxxxxyzz_0_xxxyyzz_0, g_0_xxxxxyzz_0_xxxyzzz_0, g_0_xxxxxyzz_0_xxxzzzz_0, g_0_xxxxxyzz_0_xxyyyyy_0, g_0_xxxxxyzz_0_xxyyyyz_0, g_0_xxxxxyzz_0_xxyyyzz_0, g_0_xxxxxyzz_0_xxyyzzz_0, g_0_xxxxxyzz_0_xxyzzzz_0, g_0_xxxxxyzz_0_xxzzzzz_0, g_0_xxxxxyzz_0_xyyyyyy_0, g_0_xxxxxyzz_0_xyyyyyz_0, g_0_xxxxxyzz_0_xyyyyzz_0, g_0_xxxxxyzz_0_xyyyzzz_0, g_0_xxxxxyzz_0_xyyzzzz_0, g_0_xxxxxyzz_0_xyzzzzz_0, g_0_xxxxxyzz_0_xzzzzzz_0, g_0_xxxxxyzz_0_yyyyyyy_0, g_0_xxxxxyzz_0_yyyyyyz_0, g_0_xxxxxyzz_0_yyyyyzz_0, g_0_xxxxxyzz_0_yyyyzzz_0, g_0_xxxxxyzz_0_yyyzzzz_0, g_0_xxxxxyzz_0_yyzzzzz_0, g_0_xxxxxyzz_0_yzzzzzz_0, g_0_xxxxxyzz_0_zzzzzzz_0, g_0_xxxxxzz_0_xxxxxx_1, g_0_xxxxxzz_0_xxxxxxx_0, g_0_xxxxxzz_0_xxxxxxx_1, g_0_xxxxxzz_0_xxxxxxy_0, g_0_xxxxxzz_0_xxxxxxy_1, g_0_xxxxxzz_0_xxxxxxz_0, g_0_xxxxxzz_0_xxxxxxz_1, g_0_xxxxxzz_0_xxxxxy_1, g_0_xxxxxzz_0_xxxxxyy_0, g_0_xxxxxzz_0_xxxxxyy_1, g_0_xxxxxzz_0_xxxxxyz_0, g_0_xxxxxzz_0_xxxxxyz_1, g_0_xxxxxzz_0_xxxxxz_1, g_0_xxxxxzz_0_xxxxxzz_0, g_0_xxxxxzz_0_xxxxxzz_1, g_0_xxxxxzz_0_xxxxyy_1, g_0_xxxxxzz_0_xxxxyyy_0, g_0_xxxxxzz_0_xxxxyyy_1, g_0_xxxxxzz_0_xxxxyyz_0, g_0_xxxxxzz_0_xxxxyyz_1, g_0_xxxxxzz_0_xxxxyz_1, g_0_xxxxxzz_0_xxxxyzz_0, g_0_xxxxxzz_0_xxxxyzz_1, g_0_xxxxxzz_0_xxxxzz_1, g_0_xxxxxzz_0_xxxxzzz_0, g_0_xxxxxzz_0_xxxxzzz_1, g_0_xxxxxzz_0_xxxyyy_1, g_0_xxxxxzz_0_xxxyyyy_0, g_0_xxxxxzz_0_xxxyyyy_1, g_0_xxxxxzz_0_xxxyyyz_0, g_0_xxxxxzz_0_xxxyyyz_1, g_0_xxxxxzz_0_xxxyyz_1, g_0_xxxxxzz_0_xxxyyzz_0, g_0_xxxxxzz_0_xxxyyzz_1, g_0_xxxxxzz_0_xxxyzz_1, g_0_xxxxxzz_0_xxxyzzz_0, g_0_xxxxxzz_0_xxxyzzz_1, g_0_xxxxxzz_0_xxxzzz_1, g_0_xxxxxzz_0_xxxzzzz_0, g_0_xxxxxzz_0_xxxzzzz_1, g_0_xxxxxzz_0_xxyyyy_1, g_0_xxxxxzz_0_xxyyyyy_0, g_0_xxxxxzz_0_xxyyyyy_1, g_0_xxxxxzz_0_xxyyyyz_0, g_0_xxxxxzz_0_xxyyyyz_1, g_0_xxxxxzz_0_xxyyyz_1, g_0_xxxxxzz_0_xxyyyzz_0, g_0_xxxxxzz_0_xxyyyzz_1, g_0_xxxxxzz_0_xxyyzz_1, g_0_xxxxxzz_0_xxyyzzz_0, g_0_xxxxxzz_0_xxyyzzz_1, g_0_xxxxxzz_0_xxyzzz_1, g_0_xxxxxzz_0_xxyzzzz_0, g_0_xxxxxzz_0_xxyzzzz_1, g_0_xxxxxzz_0_xxzzzz_1, g_0_xxxxxzz_0_xxzzzzz_0, g_0_xxxxxzz_0_xxzzzzz_1, g_0_xxxxxzz_0_xyyyyy_1, g_0_xxxxxzz_0_xyyyyyy_0, g_0_xxxxxzz_0_xyyyyyy_1, g_0_xxxxxzz_0_xyyyyyz_0, g_0_xxxxxzz_0_xyyyyyz_1, g_0_xxxxxzz_0_xyyyyz_1, g_0_xxxxxzz_0_xyyyyzz_0, g_0_xxxxxzz_0_xyyyyzz_1, g_0_xxxxxzz_0_xyyyzz_1, g_0_xxxxxzz_0_xyyyzzz_0, g_0_xxxxxzz_0_xyyyzzz_1, g_0_xxxxxzz_0_xyyzzz_1, g_0_xxxxxzz_0_xyyzzzz_0, g_0_xxxxxzz_0_xyyzzzz_1, g_0_xxxxxzz_0_xyzzzz_1, g_0_xxxxxzz_0_xyzzzzz_0, g_0_xxxxxzz_0_xyzzzzz_1, g_0_xxxxxzz_0_xzzzzz_1, g_0_xxxxxzz_0_xzzzzzz_0, g_0_xxxxxzz_0_xzzzzzz_1, g_0_xxxxxzz_0_yyyyyy_1, g_0_xxxxxzz_0_yyyyyyy_0, g_0_xxxxxzz_0_yyyyyyy_1, g_0_xxxxxzz_0_yyyyyyz_0, g_0_xxxxxzz_0_yyyyyyz_1, g_0_xxxxxzz_0_yyyyyz_1, g_0_xxxxxzz_0_yyyyyzz_0, g_0_xxxxxzz_0_yyyyyzz_1, g_0_xxxxxzz_0_yyyyzz_1, g_0_xxxxxzz_0_yyyyzzz_0, g_0_xxxxxzz_0_yyyyzzz_1, g_0_xxxxxzz_0_yyyzzz_1, g_0_xxxxxzz_0_yyyzzzz_0, g_0_xxxxxzz_0_yyyzzzz_1, g_0_xxxxxzz_0_yyzzzz_1, g_0_xxxxxzz_0_yyzzzzz_0, g_0_xxxxxzz_0_yyzzzzz_1, g_0_xxxxxzz_0_yzzzzz_1, g_0_xxxxxzz_0_yzzzzzz_0, g_0_xxxxxzz_0_yzzzzzz_1, g_0_xxxxxzz_0_zzzzzz_1, g_0_xxxxxzz_0_zzzzzzz_0, g_0_xxxxxzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyzz_0_xxxxxxx_0[i] = g_0_xxxxxzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxxxy_0[i] = g_0_xxxxxzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxxxy_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxxxz_0[i] = g_0_xxxxxzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxxxzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxxyy_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxxyz_0[i] = g_0_xxxxxzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxxyz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxxzz_0[i] = g_0_xxxxxzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxxxxzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxyyy_0[i] * pb_y + g_0_xxxxxzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxxzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxyyz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxyzz_0[i] = g_0_xxxxxzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxyzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxzzz_0[i] = g_0_xxxxxzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxxxxzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyyyy_0[i] * pb_y + g_0_xxxxxzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxxzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyyyz_0[i] * pb_y + g_0_xxxxxzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxxzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyyzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyzzz_0[i] = g_0_xxxxxzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxzzzz_0[i] = g_0_xxxxxzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxxxxzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyyyy_0[i] * pb_y + g_0_xxxxxzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxxzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyyyz_0[i] * pb_y + g_0_xxxxxzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxxzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyyzz_0[i] * pb_y + g_0_xxxxxzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxxzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyzzzz_0[i] = g_0_xxxxxzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxzzzzz_0[i] = g_0_xxxxxzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxxxxzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyyyy_0[i] * pb_y + g_0_xxxxxzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxxzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyyyz_0[i] * pb_y + g_0_xxxxxzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxxzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyyzz_0[i] * pb_y + g_0_xxxxxzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxxzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyzzz_0[i] * pb_y + g_0_xxxxxzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxxzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyzzzzz_0[i] = g_0_xxxxxzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xzzzzzz_0[i] = g_0_xxxxxzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxxxxzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyyyy_0[i] * pb_y + g_0_xxxxxzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxxzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyyyz_0[i] * pb_y + g_0_xxxxxzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxxzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyyzz_0[i] * pb_y + g_0_xxxxxzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxxzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyzzz_0[i] * pb_y + g_0_xxxxxzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxxzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyzzzz_0[i] * pb_y + g_0_xxxxxzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxxzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yzzzzzz_0[i] = g_0_xxxxxzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yzzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_zzzzzzz_0[i] = g_0_xxxxxzz_0_zzzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 324-360 components of targeted buffer : SLSK

    auto g_0_xxxxxzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 324);

    auto g_0_xxxxxzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 325);

    auto g_0_xxxxxzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 326);

    auto g_0_xxxxxzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 327);

    auto g_0_xxxxxzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 328);

    auto g_0_xxxxxzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 329);

    auto g_0_xxxxxzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 330);

    auto g_0_xxxxxzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 331);

    auto g_0_xxxxxzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 332);

    auto g_0_xxxxxzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 333);

    auto g_0_xxxxxzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 334);

    auto g_0_xxxxxzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 335);

    auto g_0_xxxxxzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 336);

    auto g_0_xxxxxzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 337);

    auto g_0_xxxxxzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 338);

    auto g_0_xxxxxzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 339);

    auto g_0_xxxxxzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 340);

    auto g_0_xxxxxzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 341);

    auto g_0_xxxxxzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 342);

    auto g_0_xxxxxzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 343);

    auto g_0_xxxxxzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 344);

    auto g_0_xxxxxzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 345);

    auto g_0_xxxxxzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 346);

    auto g_0_xxxxxzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 347);

    auto g_0_xxxxxzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 348);

    auto g_0_xxxxxzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 349);

    auto g_0_xxxxxzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 350);

    auto g_0_xxxxxzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 351);

    auto g_0_xxxxxzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 352);

    auto g_0_xxxxxzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 353);

    auto g_0_xxxxxzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 354);

    auto g_0_xxxxxzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 355);

    auto g_0_xxxxxzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 356);

    auto g_0_xxxxxzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 357);

    auto g_0_xxxxxzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 358);

    auto g_0_xxxxxzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 359);

    #pragma omp simd aligned(g_0_xxxxxz_0_xxxxxxx_0, g_0_xxxxxz_0_xxxxxxx_1, g_0_xxxxxz_0_xxxxxxy_0, g_0_xxxxxz_0_xxxxxxy_1, g_0_xxxxxz_0_xxxxxyy_0, g_0_xxxxxz_0_xxxxxyy_1, g_0_xxxxxz_0_xxxxyyy_0, g_0_xxxxxz_0_xxxxyyy_1, g_0_xxxxxz_0_xxxyyyy_0, g_0_xxxxxz_0_xxxyyyy_1, g_0_xxxxxz_0_xxyyyyy_0, g_0_xxxxxz_0_xxyyyyy_1, g_0_xxxxxz_0_xyyyyyy_0, g_0_xxxxxz_0_xyyyyyy_1, g_0_xxxxxzz_0_xxxxxxx_0, g_0_xxxxxzz_0_xxxxxxx_1, g_0_xxxxxzz_0_xxxxxxy_0, g_0_xxxxxzz_0_xxxxxxy_1, g_0_xxxxxzz_0_xxxxxyy_0, g_0_xxxxxzz_0_xxxxxyy_1, g_0_xxxxxzz_0_xxxxyyy_0, g_0_xxxxxzz_0_xxxxyyy_1, g_0_xxxxxzz_0_xxxyyyy_0, g_0_xxxxxzz_0_xxxyyyy_1, g_0_xxxxxzz_0_xxyyyyy_0, g_0_xxxxxzz_0_xxyyyyy_1, g_0_xxxxxzz_0_xyyyyyy_0, g_0_xxxxxzz_0_xyyyyyy_1, g_0_xxxxxzzz_0_xxxxxxx_0, g_0_xxxxxzzz_0_xxxxxxy_0, g_0_xxxxxzzz_0_xxxxxxz_0, g_0_xxxxxzzz_0_xxxxxyy_0, g_0_xxxxxzzz_0_xxxxxyz_0, g_0_xxxxxzzz_0_xxxxxzz_0, g_0_xxxxxzzz_0_xxxxyyy_0, g_0_xxxxxzzz_0_xxxxyyz_0, g_0_xxxxxzzz_0_xxxxyzz_0, g_0_xxxxxzzz_0_xxxxzzz_0, g_0_xxxxxzzz_0_xxxyyyy_0, g_0_xxxxxzzz_0_xxxyyyz_0, g_0_xxxxxzzz_0_xxxyyzz_0, g_0_xxxxxzzz_0_xxxyzzz_0, g_0_xxxxxzzz_0_xxxzzzz_0, g_0_xxxxxzzz_0_xxyyyyy_0, g_0_xxxxxzzz_0_xxyyyyz_0, g_0_xxxxxzzz_0_xxyyyzz_0, g_0_xxxxxzzz_0_xxyyzzz_0, g_0_xxxxxzzz_0_xxyzzzz_0, g_0_xxxxxzzz_0_xxzzzzz_0, g_0_xxxxxzzz_0_xyyyyyy_0, g_0_xxxxxzzz_0_xyyyyyz_0, g_0_xxxxxzzz_0_xyyyyzz_0, g_0_xxxxxzzz_0_xyyyzzz_0, g_0_xxxxxzzz_0_xyyzzzz_0, g_0_xxxxxzzz_0_xyzzzzz_0, g_0_xxxxxzzz_0_xzzzzzz_0, g_0_xxxxxzzz_0_yyyyyyy_0, g_0_xxxxxzzz_0_yyyyyyz_0, g_0_xxxxxzzz_0_yyyyyzz_0, g_0_xxxxxzzz_0_yyyyzzz_0, g_0_xxxxxzzz_0_yyyzzzz_0, g_0_xxxxxzzz_0_yyzzzzz_0, g_0_xxxxxzzz_0_yzzzzzz_0, g_0_xxxxxzzz_0_zzzzzzz_0, g_0_xxxxzzz_0_xxxxxxz_0, g_0_xxxxzzz_0_xxxxxxz_1, g_0_xxxxzzz_0_xxxxxyz_0, g_0_xxxxzzz_0_xxxxxyz_1, g_0_xxxxzzz_0_xxxxxz_1, g_0_xxxxzzz_0_xxxxxzz_0, g_0_xxxxzzz_0_xxxxxzz_1, g_0_xxxxzzz_0_xxxxyyz_0, g_0_xxxxzzz_0_xxxxyyz_1, g_0_xxxxzzz_0_xxxxyz_1, g_0_xxxxzzz_0_xxxxyzz_0, g_0_xxxxzzz_0_xxxxyzz_1, g_0_xxxxzzz_0_xxxxzz_1, g_0_xxxxzzz_0_xxxxzzz_0, g_0_xxxxzzz_0_xxxxzzz_1, g_0_xxxxzzz_0_xxxyyyz_0, g_0_xxxxzzz_0_xxxyyyz_1, g_0_xxxxzzz_0_xxxyyz_1, g_0_xxxxzzz_0_xxxyyzz_0, g_0_xxxxzzz_0_xxxyyzz_1, g_0_xxxxzzz_0_xxxyzz_1, g_0_xxxxzzz_0_xxxyzzz_0, g_0_xxxxzzz_0_xxxyzzz_1, g_0_xxxxzzz_0_xxxzzz_1, g_0_xxxxzzz_0_xxxzzzz_0, g_0_xxxxzzz_0_xxxzzzz_1, g_0_xxxxzzz_0_xxyyyyz_0, g_0_xxxxzzz_0_xxyyyyz_1, g_0_xxxxzzz_0_xxyyyz_1, g_0_xxxxzzz_0_xxyyyzz_0, g_0_xxxxzzz_0_xxyyyzz_1, g_0_xxxxzzz_0_xxyyzz_1, g_0_xxxxzzz_0_xxyyzzz_0, g_0_xxxxzzz_0_xxyyzzz_1, g_0_xxxxzzz_0_xxyzzz_1, g_0_xxxxzzz_0_xxyzzzz_0, g_0_xxxxzzz_0_xxyzzzz_1, g_0_xxxxzzz_0_xxzzzz_1, g_0_xxxxzzz_0_xxzzzzz_0, g_0_xxxxzzz_0_xxzzzzz_1, g_0_xxxxzzz_0_xyyyyyz_0, g_0_xxxxzzz_0_xyyyyyz_1, g_0_xxxxzzz_0_xyyyyz_1, g_0_xxxxzzz_0_xyyyyzz_0, g_0_xxxxzzz_0_xyyyyzz_1, g_0_xxxxzzz_0_xyyyzz_1, g_0_xxxxzzz_0_xyyyzzz_0, g_0_xxxxzzz_0_xyyyzzz_1, g_0_xxxxzzz_0_xyyzzz_1, g_0_xxxxzzz_0_xyyzzzz_0, g_0_xxxxzzz_0_xyyzzzz_1, g_0_xxxxzzz_0_xyzzzz_1, g_0_xxxxzzz_0_xyzzzzz_0, g_0_xxxxzzz_0_xyzzzzz_1, g_0_xxxxzzz_0_xzzzzz_1, g_0_xxxxzzz_0_xzzzzzz_0, g_0_xxxxzzz_0_xzzzzzz_1, g_0_xxxxzzz_0_yyyyyyy_0, g_0_xxxxzzz_0_yyyyyyy_1, g_0_xxxxzzz_0_yyyyyyz_0, g_0_xxxxzzz_0_yyyyyyz_1, g_0_xxxxzzz_0_yyyyyz_1, g_0_xxxxzzz_0_yyyyyzz_0, g_0_xxxxzzz_0_yyyyyzz_1, g_0_xxxxzzz_0_yyyyzz_1, g_0_xxxxzzz_0_yyyyzzz_0, g_0_xxxxzzz_0_yyyyzzz_1, g_0_xxxxzzz_0_yyyzzz_1, g_0_xxxxzzz_0_yyyzzzz_0, g_0_xxxxzzz_0_yyyzzzz_1, g_0_xxxxzzz_0_yyzzzz_1, g_0_xxxxzzz_0_yyzzzzz_0, g_0_xxxxzzz_0_yyzzzzz_1, g_0_xxxxzzz_0_yzzzzz_1, g_0_xxxxzzz_0_yzzzzzz_0, g_0_xxxxzzz_0_yzzzzzz_1, g_0_xxxxzzz_0_zzzzzz_1, g_0_xxxxzzz_0_zzzzzzz_0, g_0_xxxxzzz_0_zzzzzzz_1, g_0_xxxzzz_0_xxxxxxz_0, g_0_xxxzzz_0_xxxxxxz_1, g_0_xxxzzz_0_xxxxxyz_0, g_0_xxxzzz_0_xxxxxyz_1, g_0_xxxzzz_0_xxxxxzz_0, g_0_xxxzzz_0_xxxxxzz_1, g_0_xxxzzz_0_xxxxyyz_0, g_0_xxxzzz_0_xxxxyyz_1, g_0_xxxzzz_0_xxxxyzz_0, g_0_xxxzzz_0_xxxxyzz_1, g_0_xxxzzz_0_xxxxzzz_0, g_0_xxxzzz_0_xxxxzzz_1, g_0_xxxzzz_0_xxxyyyz_0, g_0_xxxzzz_0_xxxyyyz_1, g_0_xxxzzz_0_xxxyyzz_0, g_0_xxxzzz_0_xxxyyzz_1, g_0_xxxzzz_0_xxxyzzz_0, g_0_xxxzzz_0_xxxyzzz_1, g_0_xxxzzz_0_xxxzzzz_0, g_0_xxxzzz_0_xxxzzzz_1, g_0_xxxzzz_0_xxyyyyz_0, g_0_xxxzzz_0_xxyyyyz_1, g_0_xxxzzz_0_xxyyyzz_0, g_0_xxxzzz_0_xxyyyzz_1, g_0_xxxzzz_0_xxyyzzz_0, g_0_xxxzzz_0_xxyyzzz_1, g_0_xxxzzz_0_xxyzzzz_0, g_0_xxxzzz_0_xxyzzzz_1, g_0_xxxzzz_0_xxzzzzz_0, g_0_xxxzzz_0_xxzzzzz_1, g_0_xxxzzz_0_xyyyyyz_0, g_0_xxxzzz_0_xyyyyyz_1, g_0_xxxzzz_0_xyyyyzz_0, g_0_xxxzzz_0_xyyyyzz_1, g_0_xxxzzz_0_xyyyzzz_0, g_0_xxxzzz_0_xyyyzzz_1, g_0_xxxzzz_0_xyyzzzz_0, g_0_xxxzzz_0_xyyzzzz_1, g_0_xxxzzz_0_xyzzzzz_0, g_0_xxxzzz_0_xyzzzzz_1, g_0_xxxzzz_0_xzzzzzz_0, g_0_xxxzzz_0_xzzzzzz_1, g_0_xxxzzz_0_yyyyyyy_0, g_0_xxxzzz_0_yyyyyyy_1, g_0_xxxzzz_0_yyyyyyz_0, g_0_xxxzzz_0_yyyyyyz_1, g_0_xxxzzz_0_yyyyyzz_0, g_0_xxxzzz_0_yyyyyzz_1, g_0_xxxzzz_0_yyyyzzz_0, g_0_xxxzzz_0_yyyyzzz_1, g_0_xxxzzz_0_yyyzzzz_0, g_0_xxxzzz_0_yyyzzzz_1, g_0_xxxzzz_0_yyzzzzz_0, g_0_xxxzzz_0_yyzzzzz_1, g_0_xxxzzz_0_yzzzzzz_0, g_0_xxxzzz_0_yzzzzzz_1, g_0_xxxzzz_0_zzzzzzz_0, g_0_xxxzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzzz_0_xxxxxxx_0[i] = 2.0 * g_0_xxxxxz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxxxxx_0[i] * pb_z + g_0_xxxxxzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxxxy_0[i] = 2.0 * g_0_xxxxxz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxxxxy_0[i] * pb_z + g_0_xxxxxzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxxxz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxxzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxxxz_0[i] * pb_x + g_0_xxxxzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxxxz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxxxyy_0[i] * pb_z + g_0_xxxxxzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxxyz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxxyz_0[i] * pb_x + g_0_xxxxzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxxxzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxxzz_0[i] * pb_x + g_0_xxxxzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxxyyy_0[i] = 2.0 * g_0_xxxxxz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxxyyy_0[i] * pb_z + g_0_xxxxxzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxyyz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxyyz_0[i] * pb_x + g_0_xxxxzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxxyzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxyzz_0[i] * pb_x + g_0_xxxxzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxxzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxzzz_0[i] * pb_x + g_0_xxxxzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxyyyy_0[i] = 2.0 * g_0_xxxxxz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxyyyy_0[i] * pb_z + g_0_xxxxxzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxyyyz_0[i] = 4.0 * g_0_xxxzzz_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyyyz_0[i] * pb_x + g_0_xxxxzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxyyzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyyzz_0[i] * pb_x + g_0_xxxxzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxyzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyzzz_0[i] * pb_x + g_0_xxxxzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyyyyy_0[i] = 2.0 * g_0_xxxxxz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxyyyyy_0[i] * pb_z + g_0_xxxxxzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxzzz_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyyyz_0[i] * pb_x + g_0_xxxxzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyyyzz_0[i] = 4.0 * g_0_xxxzzz_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyyzz_0[i] * pb_x + g_0_xxxxzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyyzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyzzz_0[i] * pb_x + g_0_xxxxzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyyyyy_0[i] = 2.0 * g_0_xxxxxz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xyyyyyy_0[i] * pb_z + g_0_xxxxxzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xyyyyyz_0[i] = 4.0 * g_0_xxxzzz_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyyyz_0[i] * pb_x + g_0_xxxxzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxzzz_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyyzz_0[i] * pb_x + g_0_xxxxzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyyzzz_0[i] = 4.0 * g_0_xxxzzz_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyzzz_0[i] * pb_x + g_0_xxxxzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xzzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xzzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyyyy_0[i] = 4.0 * g_0_xxxzzz_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyyyyy_0[i] * pb_x + g_0_xxxxzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyyyz_0[i] = 4.0 * g_0_xxxzzz_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyyyyz_0[i] * pb_x + g_0_xxxxzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyyzz_0[i] = 4.0 * g_0_xxxzzz_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyyyzz_0[i] * pb_x + g_0_xxxxzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxzzz_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyyzzz_0[i] * pb_x + g_0_xxxxzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxzzz_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyzzzz_0[i] * pb_x + g_0_xxxxzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yzzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yzzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_zzzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zzzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 360-396 components of targeted buffer : SLSK

    auto g_0_xxxxyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 360);

    auto g_0_xxxxyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 361);

    auto g_0_xxxxyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 362);

    auto g_0_xxxxyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 363);

    auto g_0_xxxxyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 364);

    auto g_0_xxxxyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 365);

    auto g_0_xxxxyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 366);

    auto g_0_xxxxyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 367);

    auto g_0_xxxxyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 368);

    auto g_0_xxxxyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 369);

    auto g_0_xxxxyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 370);

    auto g_0_xxxxyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 371);

    auto g_0_xxxxyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 372);

    auto g_0_xxxxyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 373);

    auto g_0_xxxxyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 374);

    auto g_0_xxxxyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 375);

    auto g_0_xxxxyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 376);

    auto g_0_xxxxyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 377);

    auto g_0_xxxxyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 378);

    auto g_0_xxxxyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 379);

    auto g_0_xxxxyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 380);

    auto g_0_xxxxyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 381);

    auto g_0_xxxxyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 382);

    auto g_0_xxxxyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 383);

    auto g_0_xxxxyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 384);

    auto g_0_xxxxyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 385);

    auto g_0_xxxxyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 386);

    auto g_0_xxxxyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 387);

    auto g_0_xxxxyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 388);

    auto g_0_xxxxyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 389);

    auto g_0_xxxxyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 390);

    auto g_0_xxxxyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 391);

    auto g_0_xxxxyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 392);

    auto g_0_xxxxyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 393);

    auto g_0_xxxxyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 394);

    auto g_0_xxxxyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 395);

    #pragma omp simd aligned(g_0_xxxxyy_0_xxxxxxx_0, g_0_xxxxyy_0_xxxxxxx_1, g_0_xxxxyy_0_xxxxxxz_0, g_0_xxxxyy_0_xxxxxxz_1, g_0_xxxxyy_0_xxxxxzz_0, g_0_xxxxyy_0_xxxxxzz_1, g_0_xxxxyy_0_xxxxzzz_0, g_0_xxxxyy_0_xxxxzzz_1, g_0_xxxxyy_0_xxxzzzz_0, g_0_xxxxyy_0_xxxzzzz_1, g_0_xxxxyy_0_xxzzzzz_0, g_0_xxxxyy_0_xxzzzzz_1, g_0_xxxxyy_0_xzzzzzz_0, g_0_xxxxyy_0_xzzzzzz_1, g_0_xxxxyyy_0_xxxxxxx_0, g_0_xxxxyyy_0_xxxxxxx_1, g_0_xxxxyyy_0_xxxxxxz_0, g_0_xxxxyyy_0_xxxxxxz_1, g_0_xxxxyyy_0_xxxxxzz_0, g_0_xxxxyyy_0_xxxxxzz_1, g_0_xxxxyyy_0_xxxxzzz_0, g_0_xxxxyyy_0_xxxxzzz_1, g_0_xxxxyyy_0_xxxzzzz_0, g_0_xxxxyyy_0_xxxzzzz_1, g_0_xxxxyyy_0_xxzzzzz_0, g_0_xxxxyyy_0_xxzzzzz_1, g_0_xxxxyyy_0_xzzzzzz_0, g_0_xxxxyyy_0_xzzzzzz_1, g_0_xxxxyyyy_0_xxxxxxx_0, g_0_xxxxyyyy_0_xxxxxxy_0, g_0_xxxxyyyy_0_xxxxxxz_0, g_0_xxxxyyyy_0_xxxxxyy_0, g_0_xxxxyyyy_0_xxxxxyz_0, g_0_xxxxyyyy_0_xxxxxzz_0, g_0_xxxxyyyy_0_xxxxyyy_0, g_0_xxxxyyyy_0_xxxxyyz_0, g_0_xxxxyyyy_0_xxxxyzz_0, g_0_xxxxyyyy_0_xxxxzzz_0, g_0_xxxxyyyy_0_xxxyyyy_0, g_0_xxxxyyyy_0_xxxyyyz_0, g_0_xxxxyyyy_0_xxxyyzz_0, g_0_xxxxyyyy_0_xxxyzzz_0, g_0_xxxxyyyy_0_xxxzzzz_0, g_0_xxxxyyyy_0_xxyyyyy_0, g_0_xxxxyyyy_0_xxyyyyz_0, g_0_xxxxyyyy_0_xxyyyzz_0, g_0_xxxxyyyy_0_xxyyzzz_0, g_0_xxxxyyyy_0_xxyzzzz_0, g_0_xxxxyyyy_0_xxzzzzz_0, g_0_xxxxyyyy_0_xyyyyyy_0, g_0_xxxxyyyy_0_xyyyyyz_0, g_0_xxxxyyyy_0_xyyyyzz_0, g_0_xxxxyyyy_0_xyyyzzz_0, g_0_xxxxyyyy_0_xyyzzzz_0, g_0_xxxxyyyy_0_xyzzzzz_0, g_0_xxxxyyyy_0_xzzzzzz_0, g_0_xxxxyyyy_0_yyyyyyy_0, g_0_xxxxyyyy_0_yyyyyyz_0, g_0_xxxxyyyy_0_yyyyyzz_0, g_0_xxxxyyyy_0_yyyyzzz_0, g_0_xxxxyyyy_0_yyyzzzz_0, g_0_xxxxyyyy_0_yyzzzzz_0, g_0_xxxxyyyy_0_yzzzzzz_0, g_0_xxxxyyyy_0_zzzzzzz_0, g_0_xxxyyyy_0_xxxxxxy_0, g_0_xxxyyyy_0_xxxxxxy_1, g_0_xxxyyyy_0_xxxxxy_1, g_0_xxxyyyy_0_xxxxxyy_0, g_0_xxxyyyy_0_xxxxxyy_1, g_0_xxxyyyy_0_xxxxxyz_0, g_0_xxxyyyy_0_xxxxxyz_1, g_0_xxxyyyy_0_xxxxyy_1, g_0_xxxyyyy_0_xxxxyyy_0, g_0_xxxyyyy_0_xxxxyyy_1, g_0_xxxyyyy_0_xxxxyyz_0, g_0_xxxyyyy_0_xxxxyyz_1, g_0_xxxyyyy_0_xxxxyz_1, g_0_xxxyyyy_0_xxxxyzz_0, g_0_xxxyyyy_0_xxxxyzz_1, g_0_xxxyyyy_0_xxxyyy_1, g_0_xxxyyyy_0_xxxyyyy_0, g_0_xxxyyyy_0_xxxyyyy_1, g_0_xxxyyyy_0_xxxyyyz_0, g_0_xxxyyyy_0_xxxyyyz_1, g_0_xxxyyyy_0_xxxyyz_1, g_0_xxxyyyy_0_xxxyyzz_0, g_0_xxxyyyy_0_xxxyyzz_1, g_0_xxxyyyy_0_xxxyzz_1, g_0_xxxyyyy_0_xxxyzzz_0, g_0_xxxyyyy_0_xxxyzzz_1, g_0_xxxyyyy_0_xxyyyy_1, g_0_xxxyyyy_0_xxyyyyy_0, g_0_xxxyyyy_0_xxyyyyy_1, g_0_xxxyyyy_0_xxyyyyz_0, g_0_xxxyyyy_0_xxyyyyz_1, g_0_xxxyyyy_0_xxyyyz_1, g_0_xxxyyyy_0_xxyyyzz_0, g_0_xxxyyyy_0_xxyyyzz_1, g_0_xxxyyyy_0_xxyyzz_1, g_0_xxxyyyy_0_xxyyzzz_0, g_0_xxxyyyy_0_xxyyzzz_1, g_0_xxxyyyy_0_xxyzzz_1, g_0_xxxyyyy_0_xxyzzzz_0, g_0_xxxyyyy_0_xxyzzzz_1, g_0_xxxyyyy_0_xyyyyy_1, g_0_xxxyyyy_0_xyyyyyy_0, g_0_xxxyyyy_0_xyyyyyy_1, g_0_xxxyyyy_0_xyyyyyz_0, g_0_xxxyyyy_0_xyyyyyz_1, g_0_xxxyyyy_0_xyyyyz_1, g_0_xxxyyyy_0_xyyyyzz_0, g_0_xxxyyyy_0_xyyyyzz_1, g_0_xxxyyyy_0_xyyyzz_1, g_0_xxxyyyy_0_xyyyzzz_0, g_0_xxxyyyy_0_xyyyzzz_1, g_0_xxxyyyy_0_xyyzzz_1, g_0_xxxyyyy_0_xyyzzzz_0, g_0_xxxyyyy_0_xyyzzzz_1, g_0_xxxyyyy_0_xyzzzz_1, g_0_xxxyyyy_0_xyzzzzz_0, g_0_xxxyyyy_0_xyzzzzz_1, g_0_xxxyyyy_0_yyyyyy_1, g_0_xxxyyyy_0_yyyyyyy_0, g_0_xxxyyyy_0_yyyyyyy_1, g_0_xxxyyyy_0_yyyyyyz_0, g_0_xxxyyyy_0_yyyyyyz_1, g_0_xxxyyyy_0_yyyyyz_1, g_0_xxxyyyy_0_yyyyyzz_0, g_0_xxxyyyy_0_yyyyyzz_1, g_0_xxxyyyy_0_yyyyzz_1, g_0_xxxyyyy_0_yyyyzzz_0, g_0_xxxyyyy_0_yyyyzzz_1, g_0_xxxyyyy_0_yyyzzz_1, g_0_xxxyyyy_0_yyyzzzz_0, g_0_xxxyyyy_0_yyyzzzz_1, g_0_xxxyyyy_0_yyzzzz_1, g_0_xxxyyyy_0_yyzzzzz_0, g_0_xxxyyyy_0_yyzzzzz_1, g_0_xxxyyyy_0_yzzzzz_1, g_0_xxxyyyy_0_yzzzzzz_0, g_0_xxxyyyy_0_yzzzzzz_1, g_0_xxxyyyy_0_zzzzzzz_0, g_0_xxxyyyy_0_zzzzzzz_1, g_0_xxyyyy_0_xxxxxxy_0, g_0_xxyyyy_0_xxxxxxy_1, g_0_xxyyyy_0_xxxxxyy_0, g_0_xxyyyy_0_xxxxxyy_1, g_0_xxyyyy_0_xxxxxyz_0, g_0_xxyyyy_0_xxxxxyz_1, g_0_xxyyyy_0_xxxxyyy_0, g_0_xxyyyy_0_xxxxyyy_1, g_0_xxyyyy_0_xxxxyyz_0, g_0_xxyyyy_0_xxxxyyz_1, g_0_xxyyyy_0_xxxxyzz_0, g_0_xxyyyy_0_xxxxyzz_1, g_0_xxyyyy_0_xxxyyyy_0, g_0_xxyyyy_0_xxxyyyy_1, g_0_xxyyyy_0_xxxyyyz_0, g_0_xxyyyy_0_xxxyyyz_1, g_0_xxyyyy_0_xxxyyzz_0, g_0_xxyyyy_0_xxxyyzz_1, g_0_xxyyyy_0_xxxyzzz_0, g_0_xxyyyy_0_xxxyzzz_1, g_0_xxyyyy_0_xxyyyyy_0, g_0_xxyyyy_0_xxyyyyy_1, g_0_xxyyyy_0_xxyyyyz_0, g_0_xxyyyy_0_xxyyyyz_1, g_0_xxyyyy_0_xxyyyzz_0, g_0_xxyyyy_0_xxyyyzz_1, g_0_xxyyyy_0_xxyyzzz_0, g_0_xxyyyy_0_xxyyzzz_1, g_0_xxyyyy_0_xxyzzzz_0, g_0_xxyyyy_0_xxyzzzz_1, g_0_xxyyyy_0_xyyyyyy_0, g_0_xxyyyy_0_xyyyyyy_1, g_0_xxyyyy_0_xyyyyyz_0, g_0_xxyyyy_0_xyyyyyz_1, g_0_xxyyyy_0_xyyyyzz_0, g_0_xxyyyy_0_xyyyyzz_1, g_0_xxyyyy_0_xyyyzzz_0, g_0_xxyyyy_0_xyyyzzz_1, g_0_xxyyyy_0_xyyzzzz_0, g_0_xxyyyy_0_xyyzzzz_1, g_0_xxyyyy_0_xyzzzzz_0, g_0_xxyyyy_0_xyzzzzz_1, g_0_xxyyyy_0_yyyyyyy_0, g_0_xxyyyy_0_yyyyyyy_1, g_0_xxyyyy_0_yyyyyyz_0, g_0_xxyyyy_0_yyyyyyz_1, g_0_xxyyyy_0_yyyyyzz_0, g_0_xxyyyy_0_yyyyyzz_1, g_0_xxyyyy_0_yyyyzzz_0, g_0_xxyyyy_0_yyyyzzz_1, g_0_xxyyyy_0_yyyzzzz_0, g_0_xxyyyy_0_yyyzzzz_1, g_0_xxyyyy_0_yyzzzzz_0, g_0_xxyyyy_0_yyzzzzz_1, g_0_xxyyyy_0_yzzzzzz_0, g_0_xxyyyy_0_yzzzzzz_1, g_0_xxyyyy_0_zzzzzzz_0, g_0_xxyyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyyy_0_xxxxxxx_0[i] = 3.0 * g_0_xxxxyy_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxxxxx_0[i] * pb_y + g_0_xxxxyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxxxxy_0[i] = 3.0 * g_0_xxyyyy_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxxxy_0[i] * pb_x + g_0_xxxyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxxxz_0[i] = 3.0 * g_0_xxxxyy_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxxxxz_0[i] * pb_y + g_0_xxxxyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxxxyy_0[i] = 3.0 * g_0_xxyyyy_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxxyy_0[i] * pb_x + g_0_xxxyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxxyz_0[i] = 3.0 * g_0_xxyyyy_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxxyz_0[i] * pb_x + g_0_xxxyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxxzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxxxzz_0[i] * pb_y + g_0_xxxxyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxxyyy_0[i] = 3.0 * g_0_xxyyyy_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxyyy_0[i] * pb_x + g_0_xxxyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxyyz_0[i] = 3.0 * g_0_xxyyyy_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxyyz_0[i] * pb_x + g_0_xxxyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxyzz_0[i] = 3.0 * g_0_xxyyyy_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxyzz_0[i] * pb_x + g_0_xxxyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxxzzz_0[i] * pb_y + g_0_xxxxyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxyyyy_0[i] = 3.0 * g_0_xxyyyy_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyyyy_0[i] * pb_x + g_0_xxxyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxyyyz_0[i] = 3.0 * g_0_xxyyyy_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyyyz_0[i] * pb_x + g_0_xxxyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxyyzz_0[i] = 3.0 * g_0_xxyyyy_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyyzz_0[i] * pb_x + g_0_xxxyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyzzz_0[i] * pb_x + g_0_xxxyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxzzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxzzzz_0[i] * pb_y + g_0_xxxxyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxyyyyy_0[i] = 3.0 * g_0_xxyyyy_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyyyy_0[i] * pb_x + g_0_xxxyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyyyyz_0[i] = 3.0 * g_0_xxyyyy_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyyyz_0[i] * pb_x + g_0_xxxyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyyyzz_0[i] = 3.0 * g_0_xxyyyy_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyyzz_0[i] * pb_x + g_0_xxxyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyzzz_0[i] * pb_x + g_0_xxxyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyzzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyzzzz_0[i] * pb_x + g_0_xxxyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxzzzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxzzzzz_0[i] * pb_y + g_0_xxxxyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xyyyyyy_0[i] = 3.0 * g_0_xxyyyy_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyyyy_0[i] * pb_x + g_0_xxxyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyyyyz_0[i] = 3.0 * g_0_xxyyyy_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyyyz_0[i] * pb_x + g_0_xxxyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyyyzz_0[i] = 3.0 * g_0_xxyyyy_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyyzz_0[i] * pb_x + g_0_xxxyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyzzz_0[i] * pb_x + g_0_xxxyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyzzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyzzzz_0[i] * pb_x + g_0_xxxyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyzzzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyzzzzz_0[i] * pb_x + g_0_xxxyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xzzzzzz_0[i] = 3.0 * g_0_xxxxyy_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xzzzzzz_0[i] * pb_y + g_0_xxxxyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_yyyyyyy_0[i] = 3.0 * g_0_xxyyyy_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyyyy_0[i] * pb_x + g_0_xxxyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyyyyz_0[i] = 3.0 * g_0_xxyyyy_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyyyz_0[i] * pb_x + g_0_xxxyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyyyzz_0[i] = 3.0 * g_0_xxyyyy_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyyzz_0[i] * pb_x + g_0_xxxyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyzzz_0[i] * pb_x + g_0_xxxyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyzzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyzzzz_0[i] * pb_x + g_0_xxxyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyzzzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyzzzzz_0[i] * pb_x + g_0_xxxyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yzzzzzz_0[i] = 3.0 * g_0_xxyyyy_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yzzzzzz_0[i] * pb_x + g_0_xxxyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_zzzzzzz_0[i] = 3.0 * g_0_xxyyyy_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_zzzzzzz_0[i] * pb_x + g_0_xxxyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 396-432 components of targeted buffer : SLSK

    auto g_0_xxxxyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 396);

    auto g_0_xxxxyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 397);

    auto g_0_xxxxyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 398);

    auto g_0_xxxxyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 399);

    auto g_0_xxxxyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 400);

    auto g_0_xxxxyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 401);

    auto g_0_xxxxyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 402);

    auto g_0_xxxxyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 403);

    auto g_0_xxxxyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 404);

    auto g_0_xxxxyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 405);

    auto g_0_xxxxyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 406);

    auto g_0_xxxxyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 407);

    auto g_0_xxxxyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 408);

    auto g_0_xxxxyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 409);

    auto g_0_xxxxyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 410);

    auto g_0_xxxxyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 411);

    auto g_0_xxxxyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 412);

    auto g_0_xxxxyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 413);

    auto g_0_xxxxyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 414);

    auto g_0_xxxxyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 415);

    auto g_0_xxxxyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 416);

    auto g_0_xxxxyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 417);

    auto g_0_xxxxyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 418);

    auto g_0_xxxxyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 419);

    auto g_0_xxxxyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 420);

    auto g_0_xxxxyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 421);

    auto g_0_xxxxyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 422);

    auto g_0_xxxxyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 423);

    auto g_0_xxxxyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 424);

    auto g_0_xxxxyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 425);

    auto g_0_xxxxyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 426);

    auto g_0_xxxxyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 427);

    auto g_0_xxxxyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 428);

    auto g_0_xxxxyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 429);

    auto g_0_xxxxyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 430);

    auto g_0_xxxxyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 431);

    #pragma omp simd aligned(g_0_xxxxyyy_0_xxxxxx_1, g_0_xxxxyyy_0_xxxxxxx_0, g_0_xxxxyyy_0_xxxxxxx_1, g_0_xxxxyyy_0_xxxxxxy_0, g_0_xxxxyyy_0_xxxxxxy_1, g_0_xxxxyyy_0_xxxxxxz_0, g_0_xxxxyyy_0_xxxxxxz_1, g_0_xxxxyyy_0_xxxxxy_1, g_0_xxxxyyy_0_xxxxxyy_0, g_0_xxxxyyy_0_xxxxxyy_1, g_0_xxxxyyy_0_xxxxxyz_0, g_0_xxxxyyy_0_xxxxxyz_1, g_0_xxxxyyy_0_xxxxxz_1, g_0_xxxxyyy_0_xxxxxzz_0, g_0_xxxxyyy_0_xxxxxzz_1, g_0_xxxxyyy_0_xxxxyy_1, g_0_xxxxyyy_0_xxxxyyy_0, g_0_xxxxyyy_0_xxxxyyy_1, g_0_xxxxyyy_0_xxxxyyz_0, g_0_xxxxyyy_0_xxxxyyz_1, g_0_xxxxyyy_0_xxxxyz_1, g_0_xxxxyyy_0_xxxxyzz_0, g_0_xxxxyyy_0_xxxxyzz_1, g_0_xxxxyyy_0_xxxxzz_1, g_0_xxxxyyy_0_xxxxzzz_0, g_0_xxxxyyy_0_xxxxzzz_1, g_0_xxxxyyy_0_xxxyyy_1, g_0_xxxxyyy_0_xxxyyyy_0, g_0_xxxxyyy_0_xxxyyyy_1, g_0_xxxxyyy_0_xxxyyyz_0, g_0_xxxxyyy_0_xxxyyyz_1, g_0_xxxxyyy_0_xxxyyz_1, g_0_xxxxyyy_0_xxxyyzz_0, g_0_xxxxyyy_0_xxxyyzz_1, g_0_xxxxyyy_0_xxxyzz_1, g_0_xxxxyyy_0_xxxyzzz_0, g_0_xxxxyyy_0_xxxyzzz_1, g_0_xxxxyyy_0_xxxzzz_1, g_0_xxxxyyy_0_xxxzzzz_0, g_0_xxxxyyy_0_xxxzzzz_1, g_0_xxxxyyy_0_xxyyyy_1, g_0_xxxxyyy_0_xxyyyyy_0, g_0_xxxxyyy_0_xxyyyyy_1, g_0_xxxxyyy_0_xxyyyyz_0, g_0_xxxxyyy_0_xxyyyyz_1, g_0_xxxxyyy_0_xxyyyz_1, g_0_xxxxyyy_0_xxyyyzz_0, g_0_xxxxyyy_0_xxyyyzz_1, g_0_xxxxyyy_0_xxyyzz_1, g_0_xxxxyyy_0_xxyyzzz_0, g_0_xxxxyyy_0_xxyyzzz_1, g_0_xxxxyyy_0_xxyzzz_1, g_0_xxxxyyy_0_xxyzzzz_0, g_0_xxxxyyy_0_xxyzzzz_1, g_0_xxxxyyy_0_xxzzzz_1, g_0_xxxxyyy_0_xxzzzzz_0, g_0_xxxxyyy_0_xxzzzzz_1, g_0_xxxxyyy_0_xyyyyy_1, g_0_xxxxyyy_0_xyyyyyy_0, g_0_xxxxyyy_0_xyyyyyy_1, g_0_xxxxyyy_0_xyyyyyz_0, g_0_xxxxyyy_0_xyyyyyz_1, g_0_xxxxyyy_0_xyyyyz_1, g_0_xxxxyyy_0_xyyyyzz_0, g_0_xxxxyyy_0_xyyyyzz_1, g_0_xxxxyyy_0_xyyyzz_1, g_0_xxxxyyy_0_xyyyzzz_0, g_0_xxxxyyy_0_xyyyzzz_1, g_0_xxxxyyy_0_xyyzzz_1, g_0_xxxxyyy_0_xyyzzzz_0, g_0_xxxxyyy_0_xyyzzzz_1, g_0_xxxxyyy_0_xyzzzz_1, g_0_xxxxyyy_0_xyzzzzz_0, g_0_xxxxyyy_0_xyzzzzz_1, g_0_xxxxyyy_0_xzzzzz_1, g_0_xxxxyyy_0_xzzzzzz_0, g_0_xxxxyyy_0_xzzzzzz_1, g_0_xxxxyyy_0_yyyyyy_1, g_0_xxxxyyy_0_yyyyyyy_0, g_0_xxxxyyy_0_yyyyyyy_1, g_0_xxxxyyy_0_yyyyyyz_0, g_0_xxxxyyy_0_yyyyyyz_1, g_0_xxxxyyy_0_yyyyyz_1, g_0_xxxxyyy_0_yyyyyzz_0, g_0_xxxxyyy_0_yyyyyzz_1, g_0_xxxxyyy_0_yyyyzz_1, g_0_xxxxyyy_0_yyyyzzz_0, g_0_xxxxyyy_0_yyyyzzz_1, g_0_xxxxyyy_0_yyyzzz_1, g_0_xxxxyyy_0_yyyzzzz_0, g_0_xxxxyyy_0_yyyzzzz_1, g_0_xxxxyyy_0_yyzzzz_1, g_0_xxxxyyy_0_yyzzzzz_0, g_0_xxxxyyy_0_yyzzzzz_1, g_0_xxxxyyy_0_yzzzzz_1, g_0_xxxxyyy_0_yzzzzzz_0, g_0_xxxxyyy_0_yzzzzzz_1, g_0_xxxxyyy_0_zzzzzz_1, g_0_xxxxyyy_0_zzzzzzz_0, g_0_xxxxyyy_0_zzzzzzz_1, g_0_xxxxyyyz_0_xxxxxxx_0, g_0_xxxxyyyz_0_xxxxxxy_0, g_0_xxxxyyyz_0_xxxxxxz_0, g_0_xxxxyyyz_0_xxxxxyy_0, g_0_xxxxyyyz_0_xxxxxyz_0, g_0_xxxxyyyz_0_xxxxxzz_0, g_0_xxxxyyyz_0_xxxxyyy_0, g_0_xxxxyyyz_0_xxxxyyz_0, g_0_xxxxyyyz_0_xxxxyzz_0, g_0_xxxxyyyz_0_xxxxzzz_0, g_0_xxxxyyyz_0_xxxyyyy_0, g_0_xxxxyyyz_0_xxxyyyz_0, g_0_xxxxyyyz_0_xxxyyzz_0, g_0_xxxxyyyz_0_xxxyzzz_0, g_0_xxxxyyyz_0_xxxzzzz_0, g_0_xxxxyyyz_0_xxyyyyy_0, g_0_xxxxyyyz_0_xxyyyyz_0, g_0_xxxxyyyz_0_xxyyyzz_0, g_0_xxxxyyyz_0_xxyyzzz_0, g_0_xxxxyyyz_0_xxyzzzz_0, g_0_xxxxyyyz_0_xxzzzzz_0, g_0_xxxxyyyz_0_xyyyyyy_0, g_0_xxxxyyyz_0_xyyyyyz_0, g_0_xxxxyyyz_0_xyyyyzz_0, g_0_xxxxyyyz_0_xyyyzzz_0, g_0_xxxxyyyz_0_xyyzzzz_0, g_0_xxxxyyyz_0_xyzzzzz_0, g_0_xxxxyyyz_0_xzzzzzz_0, g_0_xxxxyyyz_0_yyyyyyy_0, g_0_xxxxyyyz_0_yyyyyyz_0, g_0_xxxxyyyz_0_yyyyyzz_0, g_0_xxxxyyyz_0_yyyyzzz_0, g_0_xxxxyyyz_0_yyyzzzz_0, g_0_xxxxyyyz_0_yyzzzzz_0, g_0_xxxxyyyz_0_yzzzzzz_0, g_0_xxxxyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyyz_0_xxxxxxx_0[i] = g_0_xxxxyyy_0_xxxxxxx_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxxxy_0[i] = g_0_xxxxyyy_0_xxxxxxy_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxxxz_0[i] = g_0_xxxxyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxxxz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxxyy_0[i] = g_0_xxxxyyy_0_xxxxxyy_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxxyz_0[i] = g_0_xxxxyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxxyz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxxyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxxzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxyyy_0[i] = g_0_xxxxyyy_0_xxxxyyy_0[i] * pb_z + g_0_xxxxyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxyyz_0[i] = g_0_xxxxyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxyyz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxxyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxyzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxxyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyyyy_0[i] = g_0_xxxxyyy_0_xxxyyyy_0[i] * pb_z + g_0_xxxxyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyyyz_0[i] = g_0_xxxxyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyyyz_0[i] * pb_z + g_0_xxxxyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyyzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxxyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxxyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyyyy_0[i] = g_0_xxxxyyy_0_xxyyyyy_0[i] * pb_z + g_0_xxxxyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyyyz_0[i] = g_0_xxxxyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyyyz_0[i] * pb_z + g_0_xxxxyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxxyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyyzz_0[i] * pb_z + g_0_xxxxyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxxyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxxyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxxyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyyyy_0[i] = g_0_xxxxyyy_0_xyyyyyy_0[i] * pb_z + g_0_xxxxyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyyyz_0[i] = g_0_xxxxyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyyyz_0[i] * pb_z + g_0_xxxxyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxxyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyyzz_0[i] * pb_z + g_0_xxxxyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyzzz_0[i] * pb_z + g_0_xxxxyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxxyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxxyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxxyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xzzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyyyy_0[i] = g_0_xxxxyyy_0_yyyyyyy_0[i] * pb_z + g_0_xxxxyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyyyz_0[i] = g_0_xxxxyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyyyyz_0[i] * pb_z + g_0_xxxxyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxxyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyyyzz_0[i] * pb_z + g_0_xxxxyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxxyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyyzzz_0[i] * pb_z + g_0_xxxxyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxxyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyzzzz_0[i] * pb_z + g_0_xxxxyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxxyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxxyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yzzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxxyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_zzzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 432-468 components of targeted buffer : SLSK

    auto g_0_xxxxyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 432);

    auto g_0_xxxxyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 433);

    auto g_0_xxxxyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 434);

    auto g_0_xxxxyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 435);

    auto g_0_xxxxyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 436);

    auto g_0_xxxxyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 437);

    auto g_0_xxxxyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 438);

    auto g_0_xxxxyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 439);

    auto g_0_xxxxyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 440);

    auto g_0_xxxxyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 441);

    auto g_0_xxxxyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 442);

    auto g_0_xxxxyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 443);

    auto g_0_xxxxyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 444);

    auto g_0_xxxxyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 445);

    auto g_0_xxxxyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 446);

    auto g_0_xxxxyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 447);

    auto g_0_xxxxyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 448);

    auto g_0_xxxxyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 449);

    auto g_0_xxxxyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 450);

    auto g_0_xxxxyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 451);

    auto g_0_xxxxyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 452);

    auto g_0_xxxxyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 453);

    auto g_0_xxxxyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 454);

    auto g_0_xxxxyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 455);

    auto g_0_xxxxyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 456);

    auto g_0_xxxxyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 457);

    auto g_0_xxxxyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 458);

    auto g_0_xxxxyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 459);

    auto g_0_xxxxyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 460);

    auto g_0_xxxxyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 461);

    auto g_0_xxxxyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 462);

    auto g_0_xxxxyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 463);

    auto g_0_xxxxyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 464);

    auto g_0_xxxxyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 465);

    auto g_0_xxxxyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 466);

    auto g_0_xxxxyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 467);

    #pragma omp simd aligned(g_0_xxxxyy_0_xxxxxxy_0, g_0_xxxxyy_0_xxxxxxy_1, g_0_xxxxyy_0_xxxxxyy_0, g_0_xxxxyy_0_xxxxxyy_1, g_0_xxxxyy_0_xxxxyyy_0, g_0_xxxxyy_0_xxxxyyy_1, g_0_xxxxyy_0_xxxyyyy_0, g_0_xxxxyy_0_xxxyyyy_1, g_0_xxxxyy_0_xxyyyyy_0, g_0_xxxxyy_0_xxyyyyy_1, g_0_xxxxyy_0_xyyyyyy_0, g_0_xxxxyy_0_xyyyyyy_1, g_0_xxxxyyz_0_xxxxxxy_0, g_0_xxxxyyz_0_xxxxxxy_1, g_0_xxxxyyz_0_xxxxxyy_0, g_0_xxxxyyz_0_xxxxxyy_1, g_0_xxxxyyz_0_xxxxyyy_0, g_0_xxxxyyz_0_xxxxyyy_1, g_0_xxxxyyz_0_xxxyyyy_0, g_0_xxxxyyz_0_xxxyyyy_1, g_0_xxxxyyz_0_xxyyyyy_0, g_0_xxxxyyz_0_xxyyyyy_1, g_0_xxxxyyz_0_xyyyyyy_0, g_0_xxxxyyz_0_xyyyyyy_1, g_0_xxxxyyzz_0_xxxxxxx_0, g_0_xxxxyyzz_0_xxxxxxy_0, g_0_xxxxyyzz_0_xxxxxxz_0, g_0_xxxxyyzz_0_xxxxxyy_0, g_0_xxxxyyzz_0_xxxxxyz_0, g_0_xxxxyyzz_0_xxxxxzz_0, g_0_xxxxyyzz_0_xxxxyyy_0, g_0_xxxxyyzz_0_xxxxyyz_0, g_0_xxxxyyzz_0_xxxxyzz_0, g_0_xxxxyyzz_0_xxxxzzz_0, g_0_xxxxyyzz_0_xxxyyyy_0, g_0_xxxxyyzz_0_xxxyyyz_0, g_0_xxxxyyzz_0_xxxyyzz_0, g_0_xxxxyyzz_0_xxxyzzz_0, g_0_xxxxyyzz_0_xxxzzzz_0, g_0_xxxxyyzz_0_xxyyyyy_0, g_0_xxxxyyzz_0_xxyyyyz_0, g_0_xxxxyyzz_0_xxyyyzz_0, g_0_xxxxyyzz_0_xxyyzzz_0, g_0_xxxxyyzz_0_xxyzzzz_0, g_0_xxxxyyzz_0_xxzzzzz_0, g_0_xxxxyyzz_0_xyyyyyy_0, g_0_xxxxyyzz_0_xyyyyyz_0, g_0_xxxxyyzz_0_xyyyyzz_0, g_0_xxxxyyzz_0_xyyyzzz_0, g_0_xxxxyyzz_0_xyyzzzz_0, g_0_xxxxyyzz_0_xyzzzzz_0, g_0_xxxxyyzz_0_xzzzzzz_0, g_0_xxxxyyzz_0_yyyyyyy_0, g_0_xxxxyyzz_0_yyyyyyz_0, g_0_xxxxyyzz_0_yyyyyzz_0, g_0_xxxxyyzz_0_yyyyzzz_0, g_0_xxxxyyzz_0_yyyzzzz_0, g_0_xxxxyyzz_0_yyzzzzz_0, g_0_xxxxyyzz_0_yzzzzzz_0, g_0_xxxxyyzz_0_zzzzzzz_0, g_0_xxxxyzz_0_xxxxxxx_0, g_0_xxxxyzz_0_xxxxxxx_1, g_0_xxxxyzz_0_xxxxxxz_0, g_0_xxxxyzz_0_xxxxxxz_1, g_0_xxxxyzz_0_xxxxxzz_0, g_0_xxxxyzz_0_xxxxxzz_1, g_0_xxxxyzz_0_xxxxzzz_0, g_0_xxxxyzz_0_xxxxzzz_1, g_0_xxxxyzz_0_xxxzzzz_0, g_0_xxxxyzz_0_xxxzzzz_1, g_0_xxxxyzz_0_xxzzzzz_0, g_0_xxxxyzz_0_xxzzzzz_1, g_0_xxxxyzz_0_xzzzzzz_0, g_0_xxxxyzz_0_xzzzzzz_1, g_0_xxxxzz_0_xxxxxxx_0, g_0_xxxxzz_0_xxxxxxx_1, g_0_xxxxzz_0_xxxxxxz_0, g_0_xxxxzz_0_xxxxxxz_1, g_0_xxxxzz_0_xxxxxzz_0, g_0_xxxxzz_0_xxxxxzz_1, g_0_xxxxzz_0_xxxxzzz_0, g_0_xxxxzz_0_xxxxzzz_1, g_0_xxxxzz_0_xxxzzzz_0, g_0_xxxxzz_0_xxxzzzz_1, g_0_xxxxzz_0_xxzzzzz_0, g_0_xxxxzz_0_xxzzzzz_1, g_0_xxxxzz_0_xzzzzzz_0, g_0_xxxxzz_0_xzzzzzz_1, g_0_xxxyyzz_0_xxxxxyz_0, g_0_xxxyyzz_0_xxxxxyz_1, g_0_xxxyyzz_0_xxxxyyz_0, g_0_xxxyyzz_0_xxxxyyz_1, g_0_xxxyyzz_0_xxxxyz_1, g_0_xxxyyzz_0_xxxxyzz_0, g_0_xxxyyzz_0_xxxxyzz_1, g_0_xxxyyzz_0_xxxyyyz_0, g_0_xxxyyzz_0_xxxyyyz_1, g_0_xxxyyzz_0_xxxyyz_1, g_0_xxxyyzz_0_xxxyyzz_0, g_0_xxxyyzz_0_xxxyyzz_1, g_0_xxxyyzz_0_xxxyzz_1, g_0_xxxyyzz_0_xxxyzzz_0, g_0_xxxyyzz_0_xxxyzzz_1, g_0_xxxyyzz_0_xxyyyyz_0, g_0_xxxyyzz_0_xxyyyyz_1, g_0_xxxyyzz_0_xxyyyz_1, g_0_xxxyyzz_0_xxyyyzz_0, g_0_xxxyyzz_0_xxyyyzz_1, g_0_xxxyyzz_0_xxyyzz_1, g_0_xxxyyzz_0_xxyyzzz_0, g_0_xxxyyzz_0_xxyyzzz_1, g_0_xxxyyzz_0_xxyzzz_1, g_0_xxxyyzz_0_xxyzzzz_0, g_0_xxxyyzz_0_xxyzzzz_1, g_0_xxxyyzz_0_xyyyyyz_0, g_0_xxxyyzz_0_xyyyyyz_1, g_0_xxxyyzz_0_xyyyyz_1, g_0_xxxyyzz_0_xyyyyzz_0, g_0_xxxyyzz_0_xyyyyzz_1, g_0_xxxyyzz_0_xyyyzz_1, g_0_xxxyyzz_0_xyyyzzz_0, g_0_xxxyyzz_0_xyyyzzz_1, g_0_xxxyyzz_0_xyyzzz_1, g_0_xxxyyzz_0_xyyzzzz_0, g_0_xxxyyzz_0_xyyzzzz_1, g_0_xxxyyzz_0_xyzzzz_1, g_0_xxxyyzz_0_xyzzzzz_0, g_0_xxxyyzz_0_xyzzzzz_1, g_0_xxxyyzz_0_yyyyyyy_0, g_0_xxxyyzz_0_yyyyyyy_1, g_0_xxxyyzz_0_yyyyyyz_0, g_0_xxxyyzz_0_yyyyyyz_1, g_0_xxxyyzz_0_yyyyyz_1, g_0_xxxyyzz_0_yyyyyzz_0, g_0_xxxyyzz_0_yyyyyzz_1, g_0_xxxyyzz_0_yyyyzz_1, g_0_xxxyyzz_0_yyyyzzz_0, g_0_xxxyyzz_0_yyyyzzz_1, g_0_xxxyyzz_0_yyyzzz_1, g_0_xxxyyzz_0_yyyzzzz_0, g_0_xxxyyzz_0_yyyzzzz_1, g_0_xxxyyzz_0_yyzzzz_1, g_0_xxxyyzz_0_yyzzzzz_0, g_0_xxxyyzz_0_yyzzzzz_1, g_0_xxxyyzz_0_yzzzzz_1, g_0_xxxyyzz_0_yzzzzzz_0, g_0_xxxyyzz_0_yzzzzzz_1, g_0_xxxyyzz_0_zzzzzzz_0, g_0_xxxyyzz_0_zzzzzzz_1, g_0_xxyyzz_0_xxxxxyz_0, g_0_xxyyzz_0_xxxxxyz_1, g_0_xxyyzz_0_xxxxyyz_0, g_0_xxyyzz_0_xxxxyyz_1, g_0_xxyyzz_0_xxxxyzz_0, g_0_xxyyzz_0_xxxxyzz_1, g_0_xxyyzz_0_xxxyyyz_0, g_0_xxyyzz_0_xxxyyyz_1, g_0_xxyyzz_0_xxxyyzz_0, g_0_xxyyzz_0_xxxyyzz_1, g_0_xxyyzz_0_xxxyzzz_0, g_0_xxyyzz_0_xxxyzzz_1, g_0_xxyyzz_0_xxyyyyz_0, g_0_xxyyzz_0_xxyyyyz_1, g_0_xxyyzz_0_xxyyyzz_0, g_0_xxyyzz_0_xxyyyzz_1, g_0_xxyyzz_0_xxyyzzz_0, g_0_xxyyzz_0_xxyyzzz_1, g_0_xxyyzz_0_xxyzzzz_0, g_0_xxyyzz_0_xxyzzzz_1, g_0_xxyyzz_0_xyyyyyz_0, g_0_xxyyzz_0_xyyyyyz_1, g_0_xxyyzz_0_xyyyyzz_0, g_0_xxyyzz_0_xyyyyzz_1, g_0_xxyyzz_0_xyyyzzz_0, g_0_xxyyzz_0_xyyyzzz_1, g_0_xxyyzz_0_xyyzzzz_0, g_0_xxyyzz_0_xyyzzzz_1, g_0_xxyyzz_0_xyzzzzz_0, g_0_xxyyzz_0_xyzzzzz_1, g_0_xxyyzz_0_yyyyyyy_0, g_0_xxyyzz_0_yyyyyyy_1, g_0_xxyyzz_0_yyyyyyz_0, g_0_xxyyzz_0_yyyyyyz_1, g_0_xxyyzz_0_yyyyyzz_0, g_0_xxyyzz_0_yyyyyzz_1, g_0_xxyyzz_0_yyyyzzz_0, g_0_xxyyzz_0_yyyyzzz_1, g_0_xxyyzz_0_yyyzzzz_0, g_0_xxyyzz_0_yyyzzzz_1, g_0_xxyyzz_0_yyzzzzz_0, g_0_xxyyzz_0_yyzzzzz_1, g_0_xxyyzz_0_yzzzzzz_0, g_0_xxyyzz_0_yzzzzzz_1, g_0_xxyyzz_0_zzzzzzz_0, g_0_xxyyzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyzz_0_xxxxxxx_0[i] = g_0_xxxxzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxxyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxxxxy_0[i] = g_0_xxxxyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxxxxy_0[i] * pb_z + g_0_xxxxyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxxxxz_0[i] = g_0_xxxxzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxxyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxxxyy_0[i] = g_0_xxxxyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxxxyy_0[i] * pb_z + g_0_xxxxyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxxxyz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxxxyz_0[i] * pb_x + g_0_xxxyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxxxzz_0[i] = g_0_xxxxzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxxyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxxyyy_0[i] = g_0_xxxxyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxxyyy_0[i] * pb_z + g_0_xxxxyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxxyyz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxxyyz_0[i] * pb_x + g_0_xxxyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxxyzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxxyzz_0[i] * pb_x + g_0_xxxyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxxzzz_0[i] = g_0_xxxxzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxxyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxyyyy_0[i] = g_0_xxxxyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxyyyy_0[i] * pb_z + g_0_xxxxyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxyyzz_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxyyyz_0[i] * pb_x + g_0_xxxyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxyyzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxyyzz_0[i] * pb_x + g_0_xxxyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxyzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxyzzz_0[i] * pb_x + g_0_xxxyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxzzzz_0[i] = g_0_xxxxzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxxyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxyyyyy_0[i] = g_0_xxxxyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxyyyyy_0[i] * pb_z + g_0_xxxxyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxyyyyz_0[i] = 3.0 * g_0_xxyyzz_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyyyyz_0[i] * pb_x + g_0_xxxyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxyyzz_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyyyzz_0[i] * pb_x + g_0_xxxyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxyyzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyyzzz_0[i] * pb_x + g_0_xxxyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxyzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyzzzz_0[i] * pb_x + g_0_xxxyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxzzzzz_0[i] = g_0_xxxxzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxxyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xyyyyyy_0[i] = g_0_xxxxyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xyyyyyy_0[i] * pb_z + g_0_xxxxyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xyyyyyz_0[i] = 3.0 * g_0_xxyyzz_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyyyyz_0[i] * pb_x + g_0_xxxyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyyyyzz_0[i] = 3.0 * g_0_xxyyzz_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyyyzz_0[i] * pb_x + g_0_xxxyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxyyzz_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyyzzz_0[i] * pb_x + g_0_xxxyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyyzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyzzzz_0[i] * pb_x + g_0_xxxyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyzzzzz_0[i] * pb_x + g_0_xxxyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xzzzzzz_0[i] = g_0_xxxxzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxxyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_yyyyyyy_0[i] = 3.0 * g_0_xxyyzz_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyyyyy_0[i] * pb_x + g_0_xxxyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyyyyz_0[i] = 3.0 * g_0_xxyyzz_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyyyyz_0[i] * pb_x + g_0_xxxyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyyyzz_0[i] = 3.0 * g_0_xxyyzz_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyyyzz_0[i] * pb_x + g_0_xxxyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyyzzz_0[i] = 3.0 * g_0_xxyyzz_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyyzzz_0[i] * pb_x + g_0_xxxyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxyyzz_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyzzzz_0[i] * pb_x + g_0_xxxyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyzzzzz_0[i] * pb_x + g_0_xxxyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yzzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yzzzzzz_0[i] * pb_x + g_0_xxxyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_zzzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_zzzzzzz_0[i] * pb_x + g_0_xxxyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 468-504 components of targeted buffer : SLSK

    auto g_0_xxxxyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 468);

    auto g_0_xxxxyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 469);

    auto g_0_xxxxyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 470);

    auto g_0_xxxxyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 471);

    auto g_0_xxxxyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 472);

    auto g_0_xxxxyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 473);

    auto g_0_xxxxyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 474);

    auto g_0_xxxxyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 475);

    auto g_0_xxxxyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 476);

    auto g_0_xxxxyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 477);

    auto g_0_xxxxyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 478);

    auto g_0_xxxxyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 479);

    auto g_0_xxxxyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 480);

    auto g_0_xxxxyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 481);

    auto g_0_xxxxyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 482);

    auto g_0_xxxxyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 483);

    auto g_0_xxxxyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 484);

    auto g_0_xxxxyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 485);

    auto g_0_xxxxyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 486);

    auto g_0_xxxxyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 487);

    auto g_0_xxxxyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 488);

    auto g_0_xxxxyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 489);

    auto g_0_xxxxyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 490);

    auto g_0_xxxxyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 491);

    auto g_0_xxxxyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 492);

    auto g_0_xxxxyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 493);

    auto g_0_xxxxyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 494);

    auto g_0_xxxxyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 495);

    auto g_0_xxxxyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 496);

    auto g_0_xxxxyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 497);

    auto g_0_xxxxyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 498);

    auto g_0_xxxxyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 499);

    auto g_0_xxxxyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 500);

    auto g_0_xxxxyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 501);

    auto g_0_xxxxyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 502);

    auto g_0_xxxxyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 503);

    #pragma omp simd aligned(g_0_xxxxyzzz_0_xxxxxxx_0, g_0_xxxxyzzz_0_xxxxxxy_0, g_0_xxxxyzzz_0_xxxxxxz_0, g_0_xxxxyzzz_0_xxxxxyy_0, g_0_xxxxyzzz_0_xxxxxyz_0, g_0_xxxxyzzz_0_xxxxxzz_0, g_0_xxxxyzzz_0_xxxxyyy_0, g_0_xxxxyzzz_0_xxxxyyz_0, g_0_xxxxyzzz_0_xxxxyzz_0, g_0_xxxxyzzz_0_xxxxzzz_0, g_0_xxxxyzzz_0_xxxyyyy_0, g_0_xxxxyzzz_0_xxxyyyz_0, g_0_xxxxyzzz_0_xxxyyzz_0, g_0_xxxxyzzz_0_xxxyzzz_0, g_0_xxxxyzzz_0_xxxzzzz_0, g_0_xxxxyzzz_0_xxyyyyy_0, g_0_xxxxyzzz_0_xxyyyyz_0, g_0_xxxxyzzz_0_xxyyyzz_0, g_0_xxxxyzzz_0_xxyyzzz_0, g_0_xxxxyzzz_0_xxyzzzz_0, g_0_xxxxyzzz_0_xxzzzzz_0, g_0_xxxxyzzz_0_xyyyyyy_0, g_0_xxxxyzzz_0_xyyyyyz_0, g_0_xxxxyzzz_0_xyyyyzz_0, g_0_xxxxyzzz_0_xyyyzzz_0, g_0_xxxxyzzz_0_xyyzzzz_0, g_0_xxxxyzzz_0_xyzzzzz_0, g_0_xxxxyzzz_0_xzzzzzz_0, g_0_xxxxyzzz_0_yyyyyyy_0, g_0_xxxxyzzz_0_yyyyyyz_0, g_0_xxxxyzzz_0_yyyyyzz_0, g_0_xxxxyzzz_0_yyyyzzz_0, g_0_xxxxyzzz_0_yyyzzzz_0, g_0_xxxxyzzz_0_yyzzzzz_0, g_0_xxxxyzzz_0_yzzzzzz_0, g_0_xxxxyzzz_0_zzzzzzz_0, g_0_xxxxzzz_0_xxxxxx_1, g_0_xxxxzzz_0_xxxxxxx_0, g_0_xxxxzzz_0_xxxxxxx_1, g_0_xxxxzzz_0_xxxxxxy_0, g_0_xxxxzzz_0_xxxxxxy_1, g_0_xxxxzzz_0_xxxxxxz_0, g_0_xxxxzzz_0_xxxxxxz_1, g_0_xxxxzzz_0_xxxxxy_1, g_0_xxxxzzz_0_xxxxxyy_0, g_0_xxxxzzz_0_xxxxxyy_1, g_0_xxxxzzz_0_xxxxxyz_0, g_0_xxxxzzz_0_xxxxxyz_1, g_0_xxxxzzz_0_xxxxxz_1, g_0_xxxxzzz_0_xxxxxzz_0, g_0_xxxxzzz_0_xxxxxzz_1, g_0_xxxxzzz_0_xxxxyy_1, g_0_xxxxzzz_0_xxxxyyy_0, g_0_xxxxzzz_0_xxxxyyy_1, g_0_xxxxzzz_0_xxxxyyz_0, g_0_xxxxzzz_0_xxxxyyz_1, g_0_xxxxzzz_0_xxxxyz_1, g_0_xxxxzzz_0_xxxxyzz_0, g_0_xxxxzzz_0_xxxxyzz_1, g_0_xxxxzzz_0_xxxxzz_1, g_0_xxxxzzz_0_xxxxzzz_0, g_0_xxxxzzz_0_xxxxzzz_1, g_0_xxxxzzz_0_xxxyyy_1, g_0_xxxxzzz_0_xxxyyyy_0, g_0_xxxxzzz_0_xxxyyyy_1, g_0_xxxxzzz_0_xxxyyyz_0, g_0_xxxxzzz_0_xxxyyyz_1, g_0_xxxxzzz_0_xxxyyz_1, g_0_xxxxzzz_0_xxxyyzz_0, g_0_xxxxzzz_0_xxxyyzz_1, g_0_xxxxzzz_0_xxxyzz_1, g_0_xxxxzzz_0_xxxyzzz_0, g_0_xxxxzzz_0_xxxyzzz_1, g_0_xxxxzzz_0_xxxzzz_1, g_0_xxxxzzz_0_xxxzzzz_0, g_0_xxxxzzz_0_xxxzzzz_1, g_0_xxxxzzz_0_xxyyyy_1, g_0_xxxxzzz_0_xxyyyyy_0, g_0_xxxxzzz_0_xxyyyyy_1, g_0_xxxxzzz_0_xxyyyyz_0, g_0_xxxxzzz_0_xxyyyyz_1, g_0_xxxxzzz_0_xxyyyz_1, g_0_xxxxzzz_0_xxyyyzz_0, g_0_xxxxzzz_0_xxyyyzz_1, g_0_xxxxzzz_0_xxyyzz_1, g_0_xxxxzzz_0_xxyyzzz_0, g_0_xxxxzzz_0_xxyyzzz_1, g_0_xxxxzzz_0_xxyzzz_1, g_0_xxxxzzz_0_xxyzzzz_0, g_0_xxxxzzz_0_xxyzzzz_1, g_0_xxxxzzz_0_xxzzzz_1, g_0_xxxxzzz_0_xxzzzzz_0, g_0_xxxxzzz_0_xxzzzzz_1, g_0_xxxxzzz_0_xyyyyy_1, g_0_xxxxzzz_0_xyyyyyy_0, g_0_xxxxzzz_0_xyyyyyy_1, g_0_xxxxzzz_0_xyyyyyz_0, g_0_xxxxzzz_0_xyyyyyz_1, g_0_xxxxzzz_0_xyyyyz_1, g_0_xxxxzzz_0_xyyyyzz_0, g_0_xxxxzzz_0_xyyyyzz_1, g_0_xxxxzzz_0_xyyyzz_1, g_0_xxxxzzz_0_xyyyzzz_0, g_0_xxxxzzz_0_xyyyzzz_1, g_0_xxxxzzz_0_xyyzzz_1, g_0_xxxxzzz_0_xyyzzzz_0, g_0_xxxxzzz_0_xyyzzzz_1, g_0_xxxxzzz_0_xyzzzz_1, g_0_xxxxzzz_0_xyzzzzz_0, g_0_xxxxzzz_0_xyzzzzz_1, g_0_xxxxzzz_0_xzzzzz_1, g_0_xxxxzzz_0_xzzzzzz_0, g_0_xxxxzzz_0_xzzzzzz_1, g_0_xxxxzzz_0_yyyyyy_1, g_0_xxxxzzz_0_yyyyyyy_0, g_0_xxxxzzz_0_yyyyyyy_1, g_0_xxxxzzz_0_yyyyyyz_0, g_0_xxxxzzz_0_yyyyyyz_1, g_0_xxxxzzz_0_yyyyyz_1, g_0_xxxxzzz_0_yyyyyzz_0, g_0_xxxxzzz_0_yyyyyzz_1, g_0_xxxxzzz_0_yyyyzz_1, g_0_xxxxzzz_0_yyyyzzz_0, g_0_xxxxzzz_0_yyyyzzz_1, g_0_xxxxzzz_0_yyyzzz_1, g_0_xxxxzzz_0_yyyzzzz_0, g_0_xxxxzzz_0_yyyzzzz_1, g_0_xxxxzzz_0_yyzzzz_1, g_0_xxxxzzz_0_yyzzzzz_0, g_0_xxxxzzz_0_yyzzzzz_1, g_0_xxxxzzz_0_yzzzzz_1, g_0_xxxxzzz_0_yzzzzzz_0, g_0_xxxxzzz_0_yzzzzzz_1, g_0_xxxxzzz_0_zzzzzz_1, g_0_xxxxzzz_0_zzzzzzz_0, g_0_xxxxzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzzz_0_xxxxxxx_0[i] = g_0_xxxxzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxxxy_0[i] = g_0_xxxxzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxxxy_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxxxz_0[i] = g_0_xxxxzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxxzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxxyy_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxxyz_0[i] = g_0_xxxxzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxxyz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxxzz_0[i] = g_0_xxxxzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxxxzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxyyy_0[i] * pb_y + g_0_xxxxzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxyyz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxyzz_0[i] = g_0_xxxxzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxyzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxzzz_0[i] = g_0_xxxxzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxxxzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyyyy_0[i] * pb_y + g_0_xxxxzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyyyz_0[i] * pb_y + g_0_xxxxzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyyzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyzzz_0[i] = g_0_xxxxzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxzzzz_0[i] = g_0_xxxxzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxxxzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyyyy_0[i] * pb_y + g_0_xxxxzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyyyz_0[i] * pb_y + g_0_xxxxzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyyzz_0[i] * pb_y + g_0_xxxxzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyzzzz_0[i] = g_0_xxxxzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxzzzzz_0[i] = g_0_xxxxzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxxxzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyyyy_0[i] * pb_y + g_0_xxxxzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyyyz_0[i] * pb_y + g_0_xxxxzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyyzz_0[i] * pb_y + g_0_xxxxzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyzzz_0[i] * pb_y + g_0_xxxxzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyzzzzz_0[i] = g_0_xxxxzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xzzzzzz_0[i] = g_0_xxxxzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxxxzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyyyy_0[i] * pb_y + g_0_xxxxzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyyyz_0[i] * pb_y + g_0_xxxxzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyyzz_0[i] * pb_y + g_0_xxxxzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyzzz_0[i] * pb_y + g_0_xxxxzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyzzzz_0[i] * pb_y + g_0_xxxxzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yzzzzzz_0[i] = g_0_xxxxzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yzzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_zzzzzzz_0[i] = g_0_xxxxzzz_0_zzzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 504-540 components of targeted buffer : SLSK

    auto g_0_xxxxzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 504);

    auto g_0_xxxxzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 505);

    auto g_0_xxxxzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 506);

    auto g_0_xxxxzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 507);

    auto g_0_xxxxzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 508);

    auto g_0_xxxxzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 509);

    auto g_0_xxxxzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 510);

    auto g_0_xxxxzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 511);

    auto g_0_xxxxzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 512);

    auto g_0_xxxxzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 513);

    auto g_0_xxxxzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 514);

    auto g_0_xxxxzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 515);

    auto g_0_xxxxzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 516);

    auto g_0_xxxxzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 517);

    auto g_0_xxxxzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 518);

    auto g_0_xxxxzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 519);

    auto g_0_xxxxzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 520);

    auto g_0_xxxxzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 521);

    auto g_0_xxxxzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 522);

    auto g_0_xxxxzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 523);

    auto g_0_xxxxzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 524);

    auto g_0_xxxxzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 525);

    auto g_0_xxxxzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 526);

    auto g_0_xxxxzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 527);

    auto g_0_xxxxzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 528);

    auto g_0_xxxxzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 529);

    auto g_0_xxxxzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 530);

    auto g_0_xxxxzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 531);

    auto g_0_xxxxzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 532);

    auto g_0_xxxxzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 533);

    auto g_0_xxxxzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 534);

    auto g_0_xxxxzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 535);

    auto g_0_xxxxzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 536);

    auto g_0_xxxxzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 537);

    auto g_0_xxxxzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 538);

    auto g_0_xxxxzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 539);

    #pragma omp simd aligned(g_0_xxxxzz_0_xxxxxxx_0, g_0_xxxxzz_0_xxxxxxx_1, g_0_xxxxzz_0_xxxxxxy_0, g_0_xxxxzz_0_xxxxxxy_1, g_0_xxxxzz_0_xxxxxyy_0, g_0_xxxxzz_0_xxxxxyy_1, g_0_xxxxzz_0_xxxxyyy_0, g_0_xxxxzz_0_xxxxyyy_1, g_0_xxxxzz_0_xxxyyyy_0, g_0_xxxxzz_0_xxxyyyy_1, g_0_xxxxzz_0_xxyyyyy_0, g_0_xxxxzz_0_xxyyyyy_1, g_0_xxxxzz_0_xyyyyyy_0, g_0_xxxxzz_0_xyyyyyy_1, g_0_xxxxzzz_0_xxxxxxx_0, g_0_xxxxzzz_0_xxxxxxx_1, g_0_xxxxzzz_0_xxxxxxy_0, g_0_xxxxzzz_0_xxxxxxy_1, g_0_xxxxzzz_0_xxxxxyy_0, g_0_xxxxzzz_0_xxxxxyy_1, g_0_xxxxzzz_0_xxxxyyy_0, g_0_xxxxzzz_0_xxxxyyy_1, g_0_xxxxzzz_0_xxxyyyy_0, g_0_xxxxzzz_0_xxxyyyy_1, g_0_xxxxzzz_0_xxyyyyy_0, g_0_xxxxzzz_0_xxyyyyy_1, g_0_xxxxzzz_0_xyyyyyy_0, g_0_xxxxzzz_0_xyyyyyy_1, g_0_xxxxzzzz_0_xxxxxxx_0, g_0_xxxxzzzz_0_xxxxxxy_0, g_0_xxxxzzzz_0_xxxxxxz_0, g_0_xxxxzzzz_0_xxxxxyy_0, g_0_xxxxzzzz_0_xxxxxyz_0, g_0_xxxxzzzz_0_xxxxxzz_0, g_0_xxxxzzzz_0_xxxxyyy_0, g_0_xxxxzzzz_0_xxxxyyz_0, g_0_xxxxzzzz_0_xxxxyzz_0, g_0_xxxxzzzz_0_xxxxzzz_0, g_0_xxxxzzzz_0_xxxyyyy_0, g_0_xxxxzzzz_0_xxxyyyz_0, g_0_xxxxzzzz_0_xxxyyzz_0, g_0_xxxxzzzz_0_xxxyzzz_0, g_0_xxxxzzzz_0_xxxzzzz_0, g_0_xxxxzzzz_0_xxyyyyy_0, g_0_xxxxzzzz_0_xxyyyyz_0, g_0_xxxxzzzz_0_xxyyyzz_0, g_0_xxxxzzzz_0_xxyyzzz_0, g_0_xxxxzzzz_0_xxyzzzz_0, g_0_xxxxzzzz_0_xxzzzzz_0, g_0_xxxxzzzz_0_xyyyyyy_0, g_0_xxxxzzzz_0_xyyyyyz_0, g_0_xxxxzzzz_0_xyyyyzz_0, g_0_xxxxzzzz_0_xyyyzzz_0, g_0_xxxxzzzz_0_xyyzzzz_0, g_0_xxxxzzzz_0_xyzzzzz_0, g_0_xxxxzzzz_0_xzzzzzz_0, g_0_xxxxzzzz_0_yyyyyyy_0, g_0_xxxxzzzz_0_yyyyyyz_0, g_0_xxxxzzzz_0_yyyyyzz_0, g_0_xxxxzzzz_0_yyyyzzz_0, g_0_xxxxzzzz_0_yyyzzzz_0, g_0_xxxxzzzz_0_yyzzzzz_0, g_0_xxxxzzzz_0_yzzzzzz_0, g_0_xxxxzzzz_0_zzzzzzz_0, g_0_xxxzzzz_0_xxxxxxz_0, g_0_xxxzzzz_0_xxxxxxz_1, g_0_xxxzzzz_0_xxxxxyz_0, g_0_xxxzzzz_0_xxxxxyz_1, g_0_xxxzzzz_0_xxxxxz_1, g_0_xxxzzzz_0_xxxxxzz_0, g_0_xxxzzzz_0_xxxxxzz_1, g_0_xxxzzzz_0_xxxxyyz_0, g_0_xxxzzzz_0_xxxxyyz_1, g_0_xxxzzzz_0_xxxxyz_1, g_0_xxxzzzz_0_xxxxyzz_0, g_0_xxxzzzz_0_xxxxyzz_1, g_0_xxxzzzz_0_xxxxzz_1, g_0_xxxzzzz_0_xxxxzzz_0, g_0_xxxzzzz_0_xxxxzzz_1, g_0_xxxzzzz_0_xxxyyyz_0, g_0_xxxzzzz_0_xxxyyyz_1, g_0_xxxzzzz_0_xxxyyz_1, g_0_xxxzzzz_0_xxxyyzz_0, g_0_xxxzzzz_0_xxxyyzz_1, g_0_xxxzzzz_0_xxxyzz_1, g_0_xxxzzzz_0_xxxyzzz_0, g_0_xxxzzzz_0_xxxyzzz_1, g_0_xxxzzzz_0_xxxzzz_1, g_0_xxxzzzz_0_xxxzzzz_0, g_0_xxxzzzz_0_xxxzzzz_1, g_0_xxxzzzz_0_xxyyyyz_0, g_0_xxxzzzz_0_xxyyyyz_1, g_0_xxxzzzz_0_xxyyyz_1, g_0_xxxzzzz_0_xxyyyzz_0, g_0_xxxzzzz_0_xxyyyzz_1, g_0_xxxzzzz_0_xxyyzz_1, g_0_xxxzzzz_0_xxyyzzz_0, g_0_xxxzzzz_0_xxyyzzz_1, g_0_xxxzzzz_0_xxyzzz_1, g_0_xxxzzzz_0_xxyzzzz_0, g_0_xxxzzzz_0_xxyzzzz_1, g_0_xxxzzzz_0_xxzzzz_1, g_0_xxxzzzz_0_xxzzzzz_0, g_0_xxxzzzz_0_xxzzzzz_1, g_0_xxxzzzz_0_xyyyyyz_0, g_0_xxxzzzz_0_xyyyyyz_1, g_0_xxxzzzz_0_xyyyyz_1, g_0_xxxzzzz_0_xyyyyzz_0, g_0_xxxzzzz_0_xyyyyzz_1, g_0_xxxzzzz_0_xyyyzz_1, g_0_xxxzzzz_0_xyyyzzz_0, g_0_xxxzzzz_0_xyyyzzz_1, g_0_xxxzzzz_0_xyyzzz_1, g_0_xxxzzzz_0_xyyzzzz_0, g_0_xxxzzzz_0_xyyzzzz_1, g_0_xxxzzzz_0_xyzzzz_1, g_0_xxxzzzz_0_xyzzzzz_0, g_0_xxxzzzz_0_xyzzzzz_1, g_0_xxxzzzz_0_xzzzzz_1, g_0_xxxzzzz_0_xzzzzzz_0, g_0_xxxzzzz_0_xzzzzzz_1, g_0_xxxzzzz_0_yyyyyyy_0, g_0_xxxzzzz_0_yyyyyyy_1, g_0_xxxzzzz_0_yyyyyyz_0, g_0_xxxzzzz_0_yyyyyyz_1, g_0_xxxzzzz_0_yyyyyz_1, g_0_xxxzzzz_0_yyyyyzz_0, g_0_xxxzzzz_0_yyyyyzz_1, g_0_xxxzzzz_0_yyyyzz_1, g_0_xxxzzzz_0_yyyyzzz_0, g_0_xxxzzzz_0_yyyyzzz_1, g_0_xxxzzzz_0_yyyzzz_1, g_0_xxxzzzz_0_yyyzzzz_0, g_0_xxxzzzz_0_yyyzzzz_1, g_0_xxxzzzz_0_yyzzzz_1, g_0_xxxzzzz_0_yyzzzzz_0, g_0_xxxzzzz_0_yyzzzzz_1, g_0_xxxzzzz_0_yzzzzz_1, g_0_xxxzzzz_0_yzzzzzz_0, g_0_xxxzzzz_0_yzzzzzz_1, g_0_xxxzzzz_0_zzzzzz_1, g_0_xxxzzzz_0_zzzzzzz_0, g_0_xxxzzzz_0_zzzzzzz_1, g_0_xxzzzz_0_xxxxxxz_0, g_0_xxzzzz_0_xxxxxxz_1, g_0_xxzzzz_0_xxxxxyz_0, g_0_xxzzzz_0_xxxxxyz_1, g_0_xxzzzz_0_xxxxxzz_0, g_0_xxzzzz_0_xxxxxzz_1, g_0_xxzzzz_0_xxxxyyz_0, g_0_xxzzzz_0_xxxxyyz_1, g_0_xxzzzz_0_xxxxyzz_0, g_0_xxzzzz_0_xxxxyzz_1, g_0_xxzzzz_0_xxxxzzz_0, g_0_xxzzzz_0_xxxxzzz_1, g_0_xxzzzz_0_xxxyyyz_0, g_0_xxzzzz_0_xxxyyyz_1, g_0_xxzzzz_0_xxxyyzz_0, g_0_xxzzzz_0_xxxyyzz_1, g_0_xxzzzz_0_xxxyzzz_0, g_0_xxzzzz_0_xxxyzzz_1, g_0_xxzzzz_0_xxxzzzz_0, g_0_xxzzzz_0_xxxzzzz_1, g_0_xxzzzz_0_xxyyyyz_0, g_0_xxzzzz_0_xxyyyyz_1, g_0_xxzzzz_0_xxyyyzz_0, g_0_xxzzzz_0_xxyyyzz_1, g_0_xxzzzz_0_xxyyzzz_0, g_0_xxzzzz_0_xxyyzzz_1, g_0_xxzzzz_0_xxyzzzz_0, g_0_xxzzzz_0_xxyzzzz_1, g_0_xxzzzz_0_xxzzzzz_0, g_0_xxzzzz_0_xxzzzzz_1, g_0_xxzzzz_0_xyyyyyz_0, g_0_xxzzzz_0_xyyyyyz_1, g_0_xxzzzz_0_xyyyyzz_0, g_0_xxzzzz_0_xyyyyzz_1, g_0_xxzzzz_0_xyyyzzz_0, g_0_xxzzzz_0_xyyyzzz_1, g_0_xxzzzz_0_xyyzzzz_0, g_0_xxzzzz_0_xyyzzzz_1, g_0_xxzzzz_0_xyzzzzz_0, g_0_xxzzzz_0_xyzzzzz_1, g_0_xxzzzz_0_xzzzzzz_0, g_0_xxzzzz_0_xzzzzzz_1, g_0_xxzzzz_0_yyyyyyy_0, g_0_xxzzzz_0_yyyyyyy_1, g_0_xxzzzz_0_yyyyyyz_0, g_0_xxzzzz_0_yyyyyyz_1, g_0_xxzzzz_0_yyyyyzz_0, g_0_xxzzzz_0_yyyyyzz_1, g_0_xxzzzz_0_yyyyzzz_0, g_0_xxzzzz_0_yyyyzzz_1, g_0_xxzzzz_0_yyyzzzz_0, g_0_xxzzzz_0_yyyzzzz_1, g_0_xxzzzz_0_yyzzzzz_0, g_0_xxzzzz_0_yyzzzzz_1, g_0_xxzzzz_0_yzzzzzz_0, g_0_xxzzzz_0_yzzzzzz_1, g_0_xxzzzz_0_zzzzzzz_0, g_0_xxzzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzzz_0_xxxxxxx_0[i] = 3.0 * g_0_xxxxzz_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxxxxx_0[i] * pb_z + g_0_xxxxzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxxxy_0[i] = 3.0 * g_0_xxxxzz_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxxxxy_0[i] * pb_z + g_0_xxxxzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxxxz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxxxz_0[i] * pb_x + g_0_xxxzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxxxyy_0[i] = 3.0 * g_0_xxxxzz_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxxxyy_0[i] * pb_z + g_0_xxxxzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxxyz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxxyz_0[i] * pb_x + g_0_xxxzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxxxzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxxzz_0[i] * pb_x + g_0_xxxzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxxyyy_0[i] * pb_z + g_0_xxxxzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxyyz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxyyz_0[i] * pb_x + g_0_xxxzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxxyzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxyzz_0[i] * pb_x + g_0_xxxzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxxzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxzzz_0[i] * pb_x + g_0_xxxzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxyyyy_0[i] * pb_z + g_0_xxxxzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyyyz_0[i] * pb_x + g_0_xxxzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyyzz_0[i] * pb_x + g_0_xxxzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyzzz_0[i] * pb_x + g_0_xxxzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyyyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxyyyyy_0[i] * pb_z + g_0_xxxxzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxyyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyyyz_0[i] * pb_x + g_0_xxxzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxzzzz_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyyzz_0[i] * pb_x + g_0_xxxzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyyzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyzzz_0[i] * pb_x + g_0_xxxzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyyyyy_0[i] = 3.0 * g_0_xxxxzz_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xyyyyyy_0[i] * pb_z + g_0_xxxxzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xyyyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyyyz_0[i] * pb_x + g_0_xxxzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyyyzz_0[i] = 3.0 * g_0_xxzzzz_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyyzz_0[i] * pb_x + g_0_xxxzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxzzzz_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyzzz_0[i] * pb_x + g_0_xxxzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xzzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xzzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyyyy_0[i] = 3.0 * g_0_xxzzzz_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyyyyy_0[i] * pb_x + g_0_xxxzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyyyz_0[i] = 3.0 * g_0_xxzzzz_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyyyyz_0[i] * pb_x + g_0_xxxzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyyzz_0[i] = 3.0 * g_0_xxzzzz_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyyyzz_0[i] * pb_x + g_0_xxxzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyzzz_0[i] = 3.0 * g_0_xxzzzz_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyyzzz_0[i] * pb_x + g_0_xxxzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxzzzz_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyzzzz_0[i] * pb_x + g_0_xxxzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yzzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yzzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_zzzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zzzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 540-576 components of targeted buffer : SLSK

    auto g_0_xxxyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 540);

    auto g_0_xxxyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 541);

    auto g_0_xxxyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 542);

    auto g_0_xxxyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 543);

    auto g_0_xxxyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 544);

    auto g_0_xxxyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 545);

    auto g_0_xxxyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 546);

    auto g_0_xxxyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 547);

    auto g_0_xxxyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 548);

    auto g_0_xxxyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 549);

    auto g_0_xxxyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 550);

    auto g_0_xxxyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 551);

    auto g_0_xxxyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 552);

    auto g_0_xxxyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 553);

    auto g_0_xxxyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 554);

    auto g_0_xxxyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 555);

    auto g_0_xxxyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 556);

    auto g_0_xxxyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 557);

    auto g_0_xxxyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 558);

    auto g_0_xxxyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 559);

    auto g_0_xxxyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 560);

    auto g_0_xxxyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 561);

    auto g_0_xxxyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 562);

    auto g_0_xxxyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 563);

    auto g_0_xxxyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 564);

    auto g_0_xxxyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 565);

    auto g_0_xxxyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 566);

    auto g_0_xxxyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 567);

    auto g_0_xxxyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 568);

    auto g_0_xxxyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 569);

    auto g_0_xxxyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 570);

    auto g_0_xxxyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 571);

    auto g_0_xxxyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 572);

    auto g_0_xxxyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 573);

    auto g_0_xxxyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 574);

    auto g_0_xxxyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 575);

    #pragma omp simd aligned(g_0_xxxyyy_0_xxxxxxx_0, g_0_xxxyyy_0_xxxxxxx_1, g_0_xxxyyy_0_xxxxxxz_0, g_0_xxxyyy_0_xxxxxxz_1, g_0_xxxyyy_0_xxxxxzz_0, g_0_xxxyyy_0_xxxxxzz_1, g_0_xxxyyy_0_xxxxzzz_0, g_0_xxxyyy_0_xxxxzzz_1, g_0_xxxyyy_0_xxxzzzz_0, g_0_xxxyyy_0_xxxzzzz_1, g_0_xxxyyy_0_xxzzzzz_0, g_0_xxxyyy_0_xxzzzzz_1, g_0_xxxyyy_0_xzzzzzz_0, g_0_xxxyyy_0_xzzzzzz_1, g_0_xxxyyyy_0_xxxxxxx_0, g_0_xxxyyyy_0_xxxxxxx_1, g_0_xxxyyyy_0_xxxxxxz_0, g_0_xxxyyyy_0_xxxxxxz_1, g_0_xxxyyyy_0_xxxxxzz_0, g_0_xxxyyyy_0_xxxxxzz_1, g_0_xxxyyyy_0_xxxxzzz_0, g_0_xxxyyyy_0_xxxxzzz_1, g_0_xxxyyyy_0_xxxzzzz_0, g_0_xxxyyyy_0_xxxzzzz_1, g_0_xxxyyyy_0_xxzzzzz_0, g_0_xxxyyyy_0_xxzzzzz_1, g_0_xxxyyyy_0_xzzzzzz_0, g_0_xxxyyyy_0_xzzzzzz_1, g_0_xxxyyyyy_0_xxxxxxx_0, g_0_xxxyyyyy_0_xxxxxxy_0, g_0_xxxyyyyy_0_xxxxxxz_0, g_0_xxxyyyyy_0_xxxxxyy_0, g_0_xxxyyyyy_0_xxxxxyz_0, g_0_xxxyyyyy_0_xxxxxzz_0, g_0_xxxyyyyy_0_xxxxyyy_0, g_0_xxxyyyyy_0_xxxxyyz_0, g_0_xxxyyyyy_0_xxxxyzz_0, g_0_xxxyyyyy_0_xxxxzzz_0, g_0_xxxyyyyy_0_xxxyyyy_0, g_0_xxxyyyyy_0_xxxyyyz_0, g_0_xxxyyyyy_0_xxxyyzz_0, g_0_xxxyyyyy_0_xxxyzzz_0, g_0_xxxyyyyy_0_xxxzzzz_0, g_0_xxxyyyyy_0_xxyyyyy_0, g_0_xxxyyyyy_0_xxyyyyz_0, g_0_xxxyyyyy_0_xxyyyzz_0, g_0_xxxyyyyy_0_xxyyzzz_0, g_0_xxxyyyyy_0_xxyzzzz_0, g_0_xxxyyyyy_0_xxzzzzz_0, g_0_xxxyyyyy_0_xyyyyyy_0, g_0_xxxyyyyy_0_xyyyyyz_0, g_0_xxxyyyyy_0_xyyyyzz_0, g_0_xxxyyyyy_0_xyyyzzz_0, g_0_xxxyyyyy_0_xyyzzzz_0, g_0_xxxyyyyy_0_xyzzzzz_0, g_0_xxxyyyyy_0_xzzzzzz_0, g_0_xxxyyyyy_0_yyyyyyy_0, g_0_xxxyyyyy_0_yyyyyyz_0, g_0_xxxyyyyy_0_yyyyyzz_0, g_0_xxxyyyyy_0_yyyyzzz_0, g_0_xxxyyyyy_0_yyyzzzz_0, g_0_xxxyyyyy_0_yyzzzzz_0, g_0_xxxyyyyy_0_yzzzzzz_0, g_0_xxxyyyyy_0_zzzzzzz_0, g_0_xxyyyyy_0_xxxxxxy_0, g_0_xxyyyyy_0_xxxxxxy_1, g_0_xxyyyyy_0_xxxxxy_1, g_0_xxyyyyy_0_xxxxxyy_0, g_0_xxyyyyy_0_xxxxxyy_1, g_0_xxyyyyy_0_xxxxxyz_0, g_0_xxyyyyy_0_xxxxxyz_1, g_0_xxyyyyy_0_xxxxyy_1, g_0_xxyyyyy_0_xxxxyyy_0, g_0_xxyyyyy_0_xxxxyyy_1, g_0_xxyyyyy_0_xxxxyyz_0, g_0_xxyyyyy_0_xxxxyyz_1, g_0_xxyyyyy_0_xxxxyz_1, g_0_xxyyyyy_0_xxxxyzz_0, g_0_xxyyyyy_0_xxxxyzz_1, g_0_xxyyyyy_0_xxxyyy_1, g_0_xxyyyyy_0_xxxyyyy_0, g_0_xxyyyyy_0_xxxyyyy_1, g_0_xxyyyyy_0_xxxyyyz_0, g_0_xxyyyyy_0_xxxyyyz_1, g_0_xxyyyyy_0_xxxyyz_1, g_0_xxyyyyy_0_xxxyyzz_0, g_0_xxyyyyy_0_xxxyyzz_1, g_0_xxyyyyy_0_xxxyzz_1, g_0_xxyyyyy_0_xxxyzzz_0, g_0_xxyyyyy_0_xxxyzzz_1, g_0_xxyyyyy_0_xxyyyy_1, g_0_xxyyyyy_0_xxyyyyy_0, g_0_xxyyyyy_0_xxyyyyy_1, g_0_xxyyyyy_0_xxyyyyz_0, g_0_xxyyyyy_0_xxyyyyz_1, g_0_xxyyyyy_0_xxyyyz_1, g_0_xxyyyyy_0_xxyyyzz_0, g_0_xxyyyyy_0_xxyyyzz_1, g_0_xxyyyyy_0_xxyyzz_1, g_0_xxyyyyy_0_xxyyzzz_0, g_0_xxyyyyy_0_xxyyzzz_1, g_0_xxyyyyy_0_xxyzzz_1, g_0_xxyyyyy_0_xxyzzzz_0, g_0_xxyyyyy_0_xxyzzzz_1, g_0_xxyyyyy_0_xyyyyy_1, g_0_xxyyyyy_0_xyyyyyy_0, g_0_xxyyyyy_0_xyyyyyy_1, g_0_xxyyyyy_0_xyyyyyz_0, g_0_xxyyyyy_0_xyyyyyz_1, g_0_xxyyyyy_0_xyyyyz_1, g_0_xxyyyyy_0_xyyyyzz_0, g_0_xxyyyyy_0_xyyyyzz_1, g_0_xxyyyyy_0_xyyyzz_1, g_0_xxyyyyy_0_xyyyzzz_0, g_0_xxyyyyy_0_xyyyzzz_1, g_0_xxyyyyy_0_xyyzzz_1, g_0_xxyyyyy_0_xyyzzzz_0, g_0_xxyyyyy_0_xyyzzzz_1, g_0_xxyyyyy_0_xyzzzz_1, g_0_xxyyyyy_0_xyzzzzz_0, g_0_xxyyyyy_0_xyzzzzz_1, g_0_xxyyyyy_0_yyyyyy_1, g_0_xxyyyyy_0_yyyyyyy_0, g_0_xxyyyyy_0_yyyyyyy_1, g_0_xxyyyyy_0_yyyyyyz_0, g_0_xxyyyyy_0_yyyyyyz_1, g_0_xxyyyyy_0_yyyyyz_1, g_0_xxyyyyy_0_yyyyyzz_0, g_0_xxyyyyy_0_yyyyyzz_1, g_0_xxyyyyy_0_yyyyzz_1, g_0_xxyyyyy_0_yyyyzzz_0, g_0_xxyyyyy_0_yyyyzzz_1, g_0_xxyyyyy_0_yyyzzz_1, g_0_xxyyyyy_0_yyyzzzz_0, g_0_xxyyyyy_0_yyyzzzz_1, g_0_xxyyyyy_0_yyzzzz_1, g_0_xxyyyyy_0_yyzzzzz_0, g_0_xxyyyyy_0_yyzzzzz_1, g_0_xxyyyyy_0_yzzzzz_1, g_0_xxyyyyy_0_yzzzzzz_0, g_0_xxyyyyy_0_yzzzzzz_1, g_0_xxyyyyy_0_zzzzzzz_0, g_0_xxyyyyy_0_zzzzzzz_1, g_0_xyyyyy_0_xxxxxxy_0, g_0_xyyyyy_0_xxxxxxy_1, g_0_xyyyyy_0_xxxxxyy_0, g_0_xyyyyy_0_xxxxxyy_1, g_0_xyyyyy_0_xxxxxyz_0, g_0_xyyyyy_0_xxxxxyz_1, g_0_xyyyyy_0_xxxxyyy_0, g_0_xyyyyy_0_xxxxyyy_1, g_0_xyyyyy_0_xxxxyyz_0, g_0_xyyyyy_0_xxxxyyz_1, g_0_xyyyyy_0_xxxxyzz_0, g_0_xyyyyy_0_xxxxyzz_1, g_0_xyyyyy_0_xxxyyyy_0, g_0_xyyyyy_0_xxxyyyy_1, g_0_xyyyyy_0_xxxyyyz_0, g_0_xyyyyy_0_xxxyyyz_1, g_0_xyyyyy_0_xxxyyzz_0, g_0_xyyyyy_0_xxxyyzz_1, g_0_xyyyyy_0_xxxyzzz_0, g_0_xyyyyy_0_xxxyzzz_1, g_0_xyyyyy_0_xxyyyyy_0, g_0_xyyyyy_0_xxyyyyy_1, g_0_xyyyyy_0_xxyyyyz_0, g_0_xyyyyy_0_xxyyyyz_1, g_0_xyyyyy_0_xxyyyzz_0, g_0_xyyyyy_0_xxyyyzz_1, g_0_xyyyyy_0_xxyyzzz_0, g_0_xyyyyy_0_xxyyzzz_1, g_0_xyyyyy_0_xxyzzzz_0, g_0_xyyyyy_0_xxyzzzz_1, g_0_xyyyyy_0_xyyyyyy_0, g_0_xyyyyy_0_xyyyyyy_1, g_0_xyyyyy_0_xyyyyyz_0, g_0_xyyyyy_0_xyyyyyz_1, g_0_xyyyyy_0_xyyyyzz_0, g_0_xyyyyy_0_xyyyyzz_1, g_0_xyyyyy_0_xyyyzzz_0, g_0_xyyyyy_0_xyyyzzz_1, g_0_xyyyyy_0_xyyzzzz_0, g_0_xyyyyy_0_xyyzzzz_1, g_0_xyyyyy_0_xyzzzzz_0, g_0_xyyyyy_0_xyzzzzz_1, g_0_xyyyyy_0_yyyyyyy_0, g_0_xyyyyy_0_yyyyyyy_1, g_0_xyyyyy_0_yyyyyyz_0, g_0_xyyyyy_0_yyyyyyz_1, g_0_xyyyyy_0_yyyyyzz_0, g_0_xyyyyy_0_yyyyyzz_1, g_0_xyyyyy_0_yyyyzzz_0, g_0_xyyyyy_0_yyyyzzz_1, g_0_xyyyyy_0_yyyzzzz_0, g_0_xyyyyy_0_yyyzzzz_1, g_0_xyyyyy_0_yyzzzzz_0, g_0_xyyyyy_0_yyzzzzz_1, g_0_xyyyyy_0_yzzzzzz_0, g_0_xyyyyy_0_yzzzzzz_1, g_0_xyyyyy_0_zzzzzzz_0, g_0_xyyyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyyy_0_xxxxxxx_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxxxxx_0[i] * pb_y + g_0_xxxyyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxxxxy_0[i] = 2.0 * g_0_xyyyyy_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxxxy_0[i] * pb_x + g_0_xxyyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxxxz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxxxxz_0[i] * pb_y + g_0_xxxyyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxxxyy_0[i] = 2.0 * g_0_xyyyyy_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxxyy_0[i] * pb_x + g_0_xxyyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxxyz_0[i] = 2.0 * g_0_xyyyyy_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxxyz_0[i] * pb_x + g_0_xxyyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxxzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxxxzz_0[i] * pb_y + g_0_xxxyyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxxyyy_0[i] = 2.0 * g_0_xyyyyy_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxyyy_0[i] * pb_x + g_0_xxyyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxyyz_0[i] = 2.0 * g_0_xyyyyy_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxyyz_0[i] * pb_x + g_0_xxyyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxyzz_0[i] = 2.0 * g_0_xyyyyy_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxyzz_0[i] * pb_x + g_0_xxyyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxxzzz_0[i] * pb_y + g_0_xxxyyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxyyyy_0[i] = 2.0 * g_0_xyyyyy_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyyyy_0[i] * pb_x + g_0_xxyyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxyyyz_0[i] = 2.0 * g_0_xyyyyy_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyyyz_0[i] * pb_x + g_0_xxyyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxyyzz_0[i] = 2.0 * g_0_xyyyyy_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyyzz_0[i] * pb_x + g_0_xxyyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxyzzz_0[i] = 2.0 * g_0_xyyyyy_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyzzz_0[i] * pb_x + g_0_xxyyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxzzzz_0[i] * pb_y + g_0_xxxyyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxyyyyy_0[i] = 2.0 * g_0_xyyyyy_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyyyy_0[i] * pb_x + g_0_xxyyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyyyyz_0[i] = 2.0 * g_0_xyyyyy_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyyyz_0[i] * pb_x + g_0_xxyyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyyyzz_0[i] = 2.0 * g_0_xyyyyy_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyyzz_0[i] * pb_x + g_0_xxyyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyyzzz_0[i] = 2.0 * g_0_xyyyyy_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyzzz_0[i] * pb_x + g_0_xxyyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyzzzz_0[i] = 2.0 * g_0_xyyyyy_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyzzzz_0[i] * pb_x + g_0_xxyyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxzzzzz_0[i] * pb_y + g_0_xxxyyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xyyyyyy_0[i] = 2.0 * g_0_xyyyyy_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyyyy_0[i] * pb_x + g_0_xxyyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyyyyz_0[i] = 2.0 * g_0_xyyyyy_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyyyz_0[i] * pb_x + g_0_xxyyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyyyzz_0[i] = 2.0 * g_0_xyyyyy_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyyzz_0[i] * pb_x + g_0_xxyyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyyzzz_0[i] = 2.0 * g_0_xyyyyy_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyzzz_0[i] * pb_x + g_0_xxyyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyzzzz_0[i] = 2.0 * g_0_xyyyyy_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyzzzz_0[i] * pb_x + g_0_xxyyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyzzzzz_0[i] = 2.0 * g_0_xyyyyy_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyzzzzz_0[i] * pb_x + g_0_xxyyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xzzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xzzzzzz_0[i] * pb_y + g_0_xxxyyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_yyyyyyy_0[i] = 2.0 * g_0_xyyyyy_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyyyy_0[i] * pb_x + g_0_xxyyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyyyyz_0[i] = 2.0 * g_0_xyyyyy_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyyyz_0[i] * pb_x + g_0_xxyyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyyyzz_0[i] = 2.0 * g_0_xyyyyy_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyyzz_0[i] * pb_x + g_0_xxyyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyyzzz_0[i] = 2.0 * g_0_xyyyyy_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyzzz_0[i] * pb_x + g_0_xxyyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyzzzz_0[i] = 2.0 * g_0_xyyyyy_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyzzzz_0[i] * pb_x + g_0_xxyyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyzzzzz_0[i] = 2.0 * g_0_xyyyyy_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyzzzzz_0[i] * pb_x + g_0_xxyyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yzzzzzz_0[i] = 2.0 * g_0_xyyyyy_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yzzzzzz_0[i] * pb_x + g_0_xxyyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_zzzzzzz_0[i] = 2.0 * g_0_xyyyyy_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_zzzzzzz_0[i] * pb_x + g_0_xxyyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 576-612 components of targeted buffer : SLSK

    auto g_0_xxxyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 576);

    auto g_0_xxxyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 577);

    auto g_0_xxxyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 578);

    auto g_0_xxxyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 579);

    auto g_0_xxxyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 580);

    auto g_0_xxxyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 581);

    auto g_0_xxxyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 582);

    auto g_0_xxxyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 583);

    auto g_0_xxxyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 584);

    auto g_0_xxxyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 585);

    auto g_0_xxxyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 586);

    auto g_0_xxxyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 587);

    auto g_0_xxxyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 588);

    auto g_0_xxxyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 589);

    auto g_0_xxxyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 590);

    auto g_0_xxxyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 591);

    auto g_0_xxxyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 592);

    auto g_0_xxxyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 593);

    auto g_0_xxxyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 594);

    auto g_0_xxxyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 595);

    auto g_0_xxxyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 596);

    auto g_0_xxxyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 597);

    auto g_0_xxxyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 598);

    auto g_0_xxxyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 599);

    auto g_0_xxxyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 600);

    auto g_0_xxxyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 601);

    auto g_0_xxxyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 602);

    auto g_0_xxxyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 603);

    auto g_0_xxxyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 604);

    auto g_0_xxxyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 605);

    auto g_0_xxxyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 606);

    auto g_0_xxxyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 607);

    auto g_0_xxxyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 608);

    auto g_0_xxxyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 609);

    auto g_0_xxxyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 610);

    auto g_0_xxxyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 611);

    #pragma omp simd aligned(g_0_xxxyyyy_0_xxxxxx_1, g_0_xxxyyyy_0_xxxxxxx_0, g_0_xxxyyyy_0_xxxxxxx_1, g_0_xxxyyyy_0_xxxxxxy_0, g_0_xxxyyyy_0_xxxxxxy_1, g_0_xxxyyyy_0_xxxxxxz_0, g_0_xxxyyyy_0_xxxxxxz_1, g_0_xxxyyyy_0_xxxxxy_1, g_0_xxxyyyy_0_xxxxxyy_0, g_0_xxxyyyy_0_xxxxxyy_1, g_0_xxxyyyy_0_xxxxxyz_0, g_0_xxxyyyy_0_xxxxxyz_1, g_0_xxxyyyy_0_xxxxxz_1, g_0_xxxyyyy_0_xxxxxzz_0, g_0_xxxyyyy_0_xxxxxzz_1, g_0_xxxyyyy_0_xxxxyy_1, g_0_xxxyyyy_0_xxxxyyy_0, g_0_xxxyyyy_0_xxxxyyy_1, g_0_xxxyyyy_0_xxxxyyz_0, g_0_xxxyyyy_0_xxxxyyz_1, g_0_xxxyyyy_0_xxxxyz_1, g_0_xxxyyyy_0_xxxxyzz_0, g_0_xxxyyyy_0_xxxxyzz_1, g_0_xxxyyyy_0_xxxxzz_1, g_0_xxxyyyy_0_xxxxzzz_0, g_0_xxxyyyy_0_xxxxzzz_1, g_0_xxxyyyy_0_xxxyyy_1, g_0_xxxyyyy_0_xxxyyyy_0, g_0_xxxyyyy_0_xxxyyyy_1, g_0_xxxyyyy_0_xxxyyyz_0, g_0_xxxyyyy_0_xxxyyyz_1, g_0_xxxyyyy_0_xxxyyz_1, g_0_xxxyyyy_0_xxxyyzz_0, g_0_xxxyyyy_0_xxxyyzz_1, g_0_xxxyyyy_0_xxxyzz_1, g_0_xxxyyyy_0_xxxyzzz_0, g_0_xxxyyyy_0_xxxyzzz_1, g_0_xxxyyyy_0_xxxzzz_1, g_0_xxxyyyy_0_xxxzzzz_0, g_0_xxxyyyy_0_xxxzzzz_1, g_0_xxxyyyy_0_xxyyyy_1, g_0_xxxyyyy_0_xxyyyyy_0, g_0_xxxyyyy_0_xxyyyyy_1, g_0_xxxyyyy_0_xxyyyyz_0, g_0_xxxyyyy_0_xxyyyyz_1, g_0_xxxyyyy_0_xxyyyz_1, g_0_xxxyyyy_0_xxyyyzz_0, g_0_xxxyyyy_0_xxyyyzz_1, g_0_xxxyyyy_0_xxyyzz_1, g_0_xxxyyyy_0_xxyyzzz_0, g_0_xxxyyyy_0_xxyyzzz_1, g_0_xxxyyyy_0_xxyzzz_1, g_0_xxxyyyy_0_xxyzzzz_0, g_0_xxxyyyy_0_xxyzzzz_1, g_0_xxxyyyy_0_xxzzzz_1, g_0_xxxyyyy_0_xxzzzzz_0, g_0_xxxyyyy_0_xxzzzzz_1, g_0_xxxyyyy_0_xyyyyy_1, g_0_xxxyyyy_0_xyyyyyy_0, g_0_xxxyyyy_0_xyyyyyy_1, g_0_xxxyyyy_0_xyyyyyz_0, g_0_xxxyyyy_0_xyyyyyz_1, g_0_xxxyyyy_0_xyyyyz_1, g_0_xxxyyyy_0_xyyyyzz_0, g_0_xxxyyyy_0_xyyyyzz_1, g_0_xxxyyyy_0_xyyyzz_1, g_0_xxxyyyy_0_xyyyzzz_0, g_0_xxxyyyy_0_xyyyzzz_1, g_0_xxxyyyy_0_xyyzzz_1, g_0_xxxyyyy_0_xyyzzzz_0, g_0_xxxyyyy_0_xyyzzzz_1, g_0_xxxyyyy_0_xyzzzz_1, g_0_xxxyyyy_0_xyzzzzz_0, g_0_xxxyyyy_0_xyzzzzz_1, g_0_xxxyyyy_0_xzzzzz_1, g_0_xxxyyyy_0_xzzzzzz_0, g_0_xxxyyyy_0_xzzzzzz_1, g_0_xxxyyyy_0_yyyyyy_1, g_0_xxxyyyy_0_yyyyyyy_0, g_0_xxxyyyy_0_yyyyyyy_1, g_0_xxxyyyy_0_yyyyyyz_0, g_0_xxxyyyy_0_yyyyyyz_1, g_0_xxxyyyy_0_yyyyyz_1, g_0_xxxyyyy_0_yyyyyzz_0, g_0_xxxyyyy_0_yyyyyzz_1, g_0_xxxyyyy_0_yyyyzz_1, g_0_xxxyyyy_0_yyyyzzz_0, g_0_xxxyyyy_0_yyyyzzz_1, g_0_xxxyyyy_0_yyyzzz_1, g_0_xxxyyyy_0_yyyzzzz_0, g_0_xxxyyyy_0_yyyzzzz_1, g_0_xxxyyyy_0_yyzzzz_1, g_0_xxxyyyy_0_yyzzzzz_0, g_0_xxxyyyy_0_yyzzzzz_1, g_0_xxxyyyy_0_yzzzzz_1, g_0_xxxyyyy_0_yzzzzzz_0, g_0_xxxyyyy_0_yzzzzzz_1, g_0_xxxyyyy_0_zzzzzz_1, g_0_xxxyyyy_0_zzzzzzz_0, g_0_xxxyyyy_0_zzzzzzz_1, g_0_xxxyyyyz_0_xxxxxxx_0, g_0_xxxyyyyz_0_xxxxxxy_0, g_0_xxxyyyyz_0_xxxxxxz_0, g_0_xxxyyyyz_0_xxxxxyy_0, g_0_xxxyyyyz_0_xxxxxyz_0, g_0_xxxyyyyz_0_xxxxxzz_0, g_0_xxxyyyyz_0_xxxxyyy_0, g_0_xxxyyyyz_0_xxxxyyz_0, g_0_xxxyyyyz_0_xxxxyzz_0, g_0_xxxyyyyz_0_xxxxzzz_0, g_0_xxxyyyyz_0_xxxyyyy_0, g_0_xxxyyyyz_0_xxxyyyz_0, g_0_xxxyyyyz_0_xxxyyzz_0, g_0_xxxyyyyz_0_xxxyzzz_0, g_0_xxxyyyyz_0_xxxzzzz_0, g_0_xxxyyyyz_0_xxyyyyy_0, g_0_xxxyyyyz_0_xxyyyyz_0, g_0_xxxyyyyz_0_xxyyyzz_0, g_0_xxxyyyyz_0_xxyyzzz_0, g_0_xxxyyyyz_0_xxyzzzz_0, g_0_xxxyyyyz_0_xxzzzzz_0, g_0_xxxyyyyz_0_xyyyyyy_0, g_0_xxxyyyyz_0_xyyyyyz_0, g_0_xxxyyyyz_0_xyyyyzz_0, g_0_xxxyyyyz_0_xyyyzzz_0, g_0_xxxyyyyz_0_xyyzzzz_0, g_0_xxxyyyyz_0_xyzzzzz_0, g_0_xxxyyyyz_0_xzzzzzz_0, g_0_xxxyyyyz_0_yyyyyyy_0, g_0_xxxyyyyz_0_yyyyyyz_0, g_0_xxxyyyyz_0_yyyyyzz_0, g_0_xxxyyyyz_0_yyyyzzz_0, g_0_xxxyyyyz_0_yyyzzzz_0, g_0_xxxyyyyz_0_yyzzzzz_0, g_0_xxxyyyyz_0_yzzzzzz_0, g_0_xxxyyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyyz_0_xxxxxxx_0[i] = g_0_xxxyyyy_0_xxxxxxx_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxxxy_0[i] = g_0_xxxyyyy_0_xxxxxxy_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxxxz_0[i] = g_0_xxxyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxxxz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxxyy_0[i] = g_0_xxxyyyy_0_xxxxxyy_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxxyz_0[i] = g_0_xxxyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxxyz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxxzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxyyy_0[i] = g_0_xxxyyyy_0_xxxxyyy_0[i] * pb_z + g_0_xxxyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxyyz_0[i] = g_0_xxxyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxyyz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxyzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyyyy_0[i] = g_0_xxxyyyy_0_xxxyyyy_0[i] * pb_z + g_0_xxxyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyyyz_0[i] = g_0_xxxyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyyyz_0[i] * pb_z + g_0_xxxyyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyyzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyyyy_0[i] = g_0_xxxyyyy_0_xxyyyyy_0[i] * pb_z + g_0_xxxyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyyyz_0[i] = g_0_xxxyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyyyz_0[i] * pb_z + g_0_xxxyyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyyzz_0[i] * pb_z + g_0_xxxyyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyyyy_0[i] = g_0_xxxyyyy_0_xyyyyyy_0[i] * pb_z + g_0_xxxyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyyyz_0[i] = g_0_xxxyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyyyz_0[i] * pb_z + g_0_xxxyyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyyzz_0[i] * pb_z + g_0_xxxyyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyzzz_0[i] * pb_z + g_0_xxxyyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xzzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyyyy_0[i] = g_0_xxxyyyy_0_yyyyyyy_0[i] * pb_z + g_0_xxxyyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyyyz_0[i] = g_0_xxxyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyyyyz_0[i] * pb_z + g_0_xxxyyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyyyzz_0[i] * pb_z + g_0_xxxyyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyyzzz_0[i] * pb_z + g_0_xxxyyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyzzzz_0[i] * pb_z + g_0_xxxyyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yzzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_zzzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 612-648 components of targeted buffer : SLSK

    auto g_0_xxxyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 612);

    auto g_0_xxxyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 613);

    auto g_0_xxxyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 614);

    auto g_0_xxxyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 615);

    auto g_0_xxxyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 616);

    auto g_0_xxxyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 617);

    auto g_0_xxxyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 618);

    auto g_0_xxxyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 619);

    auto g_0_xxxyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 620);

    auto g_0_xxxyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 621);

    auto g_0_xxxyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 622);

    auto g_0_xxxyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 623);

    auto g_0_xxxyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 624);

    auto g_0_xxxyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 625);

    auto g_0_xxxyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 626);

    auto g_0_xxxyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 627);

    auto g_0_xxxyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 628);

    auto g_0_xxxyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 629);

    auto g_0_xxxyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 630);

    auto g_0_xxxyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 631);

    auto g_0_xxxyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 632);

    auto g_0_xxxyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 633);

    auto g_0_xxxyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 634);

    auto g_0_xxxyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 635);

    auto g_0_xxxyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 636);

    auto g_0_xxxyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 637);

    auto g_0_xxxyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 638);

    auto g_0_xxxyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 639);

    auto g_0_xxxyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 640);

    auto g_0_xxxyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 641);

    auto g_0_xxxyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 642);

    auto g_0_xxxyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 643);

    auto g_0_xxxyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 644);

    auto g_0_xxxyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 645);

    auto g_0_xxxyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 646);

    auto g_0_xxxyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 647);

    #pragma omp simd aligned(g_0_xxxyyy_0_xxxxxxy_0, g_0_xxxyyy_0_xxxxxxy_1, g_0_xxxyyy_0_xxxxxyy_0, g_0_xxxyyy_0_xxxxxyy_1, g_0_xxxyyy_0_xxxxyyy_0, g_0_xxxyyy_0_xxxxyyy_1, g_0_xxxyyy_0_xxxyyyy_0, g_0_xxxyyy_0_xxxyyyy_1, g_0_xxxyyy_0_xxyyyyy_0, g_0_xxxyyy_0_xxyyyyy_1, g_0_xxxyyy_0_xyyyyyy_0, g_0_xxxyyy_0_xyyyyyy_1, g_0_xxxyyyz_0_xxxxxxy_0, g_0_xxxyyyz_0_xxxxxxy_1, g_0_xxxyyyz_0_xxxxxyy_0, g_0_xxxyyyz_0_xxxxxyy_1, g_0_xxxyyyz_0_xxxxyyy_0, g_0_xxxyyyz_0_xxxxyyy_1, g_0_xxxyyyz_0_xxxyyyy_0, g_0_xxxyyyz_0_xxxyyyy_1, g_0_xxxyyyz_0_xxyyyyy_0, g_0_xxxyyyz_0_xxyyyyy_1, g_0_xxxyyyz_0_xyyyyyy_0, g_0_xxxyyyz_0_xyyyyyy_1, g_0_xxxyyyzz_0_xxxxxxx_0, g_0_xxxyyyzz_0_xxxxxxy_0, g_0_xxxyyyzz_0_xxxxxxz_0, g_0_xxxyyyzz_0_xxxxxyy_0, g_0_xxxyyyzz_0_xxxxxyz_0, g_0_xxxyyyzz_0_xxxxxzz_0, g_0_xxxyyyzz_0_xxxxyyy_0, g_0_xxxyyyzz_0_xxxxyyz_0, g_0_xxxyyyzz_0_xxxxyzz_0, g_0_xxxyyyzz_0_xxxxzzz_0, g_0_xxxyyyzz_0_xxxyyyy_0, g_0_xxxyyyzz_0_xxxyyyz_0, g_0_xxxyyyzz_0_xxxyyzz_0, g_0_xxxyyyzz_0_xxxyzzz_0, g_0_xxxyyyzz_0_xxxzzzz_0, g_0_xxxyyyzz_0_xxyyyyy_0, g_0_xxxyyyzz_0_xxyyyyz_0, g_0_xxxyyyzz_0_xxyyyzz_0, g_0_xxxyyyzz_0_xxyyzzz_0, g_0_xxxyyyzz_0_xxyzzzz_0, g_0_xxxyyyzz_0_xxzzzzz_0, g_0_xxxyyyzz_0_xyyyyyy_0, g_0_xxxyyyzz_0_xyyyyyz_0, g_0_xxxyyyzz_0_xyyyyzz_0, g_0_xxxyyyzz_0_xyyyzzz_0, g_0_xxxyyyzz_0_xyyzzzz_0, g_0_xxxyyyzz_0_xyzzzzz_0, g_0_xxxyyyzz_0_xzzzzzz_0, g_0_xxxyyyzz_0_yyyyyyy_0, g_0_xxxyyyzz_0_yyyyyyz_0, g_0_xxxyyyzz_0_yyyyyzz_0, g_0_xxxyyyzz_0_yyyyzzz_0, g_0_xxxyyyzz_0_yyyzzzz_0, g_0_xxxyyyzz_0_yyzzzzz_0, g_0_xxxyyyzz_0_yzzzzzz_0, g_0_xxxyyyzz_0_zzzzzzz_0, g_0_xxxyyzz_0_xxxxxxx_0, g_0_xxxyyzz_0_xxxxxxx_1, g_0_xxxyyzz_0_xxxxxxz_0, g_0_xxxyyzz_0_xxxxxxz_1, g_0_xxxyyzz_0_xxxxxzz_0, g_0_xxxyyzz_0_xxxxxzz_1, g_0_xxxyyzz_0_xxxxzzz_0, g_0_xxxyyzz_0_xxxxzzz_1, g_0_xxxyyzz_0_xxxzzzz_0, g_0_xxxyyzz_0_xxxzzzz_1, g_0_xxxyyzz_0_xxzzzzz_0, g_0_xxxyyzz_0_xxzzzzz_1, g_0_xxxyyzz_0_xzzzzzz_0, g_0_xxxyyzz_0_xzzzzzz_1, g_0_xxxyzz_0_xxxxxxx_0, g_0_xxxyzz_0_xxxxxxx_1, g_0_xxxyzz_0_xxxxxxz_0, g_0_xxxyzz_0_xxxxxxz_1, g_0_xxxyzz_0_xxxxxzz_0, g_0_xxxyzz_0_xxxxxzz_1, g_0_xxxyzz_0_xxxxzzz_0, g_0_xxxyzz_0_xxxxzzz_1, g_0_xxxyzz_0_xxxzzzz_0, g_0_xxxyzz_0_xxxzzzz_1, g_0_xxxyzz_0_xxzzzzz_0, g_0_xxxyzz_0_xxzzzzz_1, g_0_xxxyzz_0_xzzzzzz_0, g_0_xxxyzz_0_xzzzzzz_1, g_0_xxyyyzz_0_xxxxxyz_0, g_0_xxyyyzz_0_xxxxxyz_1, g_0_xxyyyzz_0_xxxxyyz_0, g_0_xxyyyzz_0_xxxxyyz_1, g_0_xxyyyzz_0_xxxxyz_1, g_0_xxyyyzz_0_xxxxyzz_0, g_0_xxyyyzz_0_xxxxyzz_1, g_0_xxyyyzz_0_xxxyyyz_0, g_0_xxyyyzz_0_xxxyyyz_1, g_0_xxyyyzz_0_xxxyyz_1, g_0_xxyyyzz_0_xxxyyzz_0, g_0_xxyyyzz_0_xxxyyzz_1, g_0_xxyyyzz_0_xxxyzz_1, g_0_xxyyyzz_0_xxxyzzz_0, g_0_xxyyyzz_0_xxxyzzz_1, g_0_xxyyyzz_0_xxyyyyz_0, g_0_xxyyyzz_0_xxyyyyz_1, g_0_xxyyyzz_0_xxyyyz_1, g_0_xxyyyzz_0_xxyyyzz_0, g_0_xxyyyzz_0_xxyyyzz_1, g_0_xxyyyzz_0_xxyyzz_1, g_0_xxyyyzz_0_xxyyzzz_0, g_0_xxyyyzz_0_xxyyzzz_1, g_0_xxyyyzz_0_xxyzzz_1, g_0_xxyyyzz_0_xxyzzzz_0, g_0_xxyyyzz_0_xxyzzzz_1, g_0_xxyyyzz_0_xyyyyyz_0, g_0_xxyyyzz_0_xyyyyyz_1, g_0_xxyyyzz_0_xyyyyz_1, g_0_xxyyyzz_0_xyyyyzz_0, g_0_xxyyyzz_0_xyyyyzz_1, g_0_xxyyyzz_0_xyyyzz_1, g_0_xxyyyzz_0_xyyyzzz_0, g_0_xxyyyzz_0_xyyyzzz_1, g_0_xxyyyzz_0_xyyzzz_1, g_0_xxyyyzz_0_xyyzzzz_0, g_0_xxyyyzz_0_xyyzzzz_1, g_0_xxyyyzz_0_xyzzzz_1, g_0_xxyyyzz_0_xyzzzzz_0, g_0_xxyyyzz_0_xyzzzzz_1, g_0_xxyyyzz_0_yyyyyyy_0, g_0_xxyyyzz_0_yyyyyyy_1, g_0_xxyyyzz_0_yyyyyyz_0, g_0_xxyyyzz_0_yyyyyyz_1, g_0_xxyyyzz_0_yyyyyz_1, g_0_xxyyyzz_0_yyyyyzz_0, g_0_xxyyyzz_0_yyyyyzz_1, g_0_xxyyyzz_0_yyyyzz_1, g_0_xxyyyzz_0_yyyyzzz_0, g_0_xxyyyzz_0_yyyyzzz_1, g_0_xxyyyzz_0_yyyzzz_1, g_0_xxyyyzz_0_yyyzzzz_0, g_0_xxyyyzz_0_yyyzzzz_1, g_0_xxyyyzz_0_yyzzzz_1, g_0_xxyyyzz_0_yyzzzzz_0, g_0_xxyyyzz_0_yyzzzzz_1, g_0_xxyyyzz_0_yzzzzz_1, g_0_xxyyyzz_0_yzzzzzz_0, g_0_xxyyyzz_0_yzzzzzz_1, g_0_xxyyyzz_0_zzzzzzz_0, g_0_xxyyyzz_0_zzzzzzz_1, g_0_xyyyzz_0_xxxxxyz_0, g_0_xyyyzz_0_xxxxxyz_1, g_0_xyyyzz_0_xxxxyyz_0, g_0_xyyyzz_0_xxxxyyz_1, g_0_xyyyzz_0_xxxxyzz_0, g_0_xyyyzz_0_xxxxyzz_1, g_0_xyyyzz_0_xxxyyyz_0, g_0_xyyyzz_0_xxxyyyz_1, g_0_xyyyzz_0_xxxyyzz_0, g_0_xyyyzz_0_xxxyyzz_1, g_0_xyyyzz_0_xxxyzzz_0, g_0_xyyyzz_0_xxxyzzz_1, g_0_xyyyzz_0_xxyyyyz_0, g_0_xyyyzz_0_xxyyyyz_1, g_0_xyyyzz_0_xxyyyzz_0, g_0_xyyyzz_0_xxyyyzz_1, g_0_xyyyzz_0_xxyyzzz_0, g_0_xyyyzz_0_xxyyzzz_1, g_0_xyyyzz_0_xxyzzzz_0, g_0_xyyyzz_0_xxyzzzz_1, g_0_xyyyzz_0_xyyyyyz_0, g_0_xyyyzz_0_xyyyyyz_1, g_0_xyyyzz_0_xyyyyzz_0, g_0_xyyyzz_0_xyyyyzz_1, g_0_xyyyzz_0_xyyyzzz_0, g_0_xyyyzz_0_xyyyzzz_1, g_0_xyyyzz_0_xyyzzzz_0, g_0_xyyyzz_0_xyyzzzz_1, g_0_xyyyzz_0_xyzzzzz_0, g_0_xyyyzz_0_xyzzzzz_1, g_0_xyyyzz_0_yyyyyyy_0, g_0_xyyyzz_0_yyyyyyy_1, g_0_xyyyzz_0_yyyyyyz_0, g_0_xyyyzz_0_yyyyyyz_1, g_0_xyyyzz_0_yyyyyzz_0, g_0_xyyyzz_0_yyyyyzz_1, g_0_xyyyzz_0_yyyyzzz_0, g_0_xyyyzz_0_yyyyzzz_1, g_0_xyyyzz_0_yyyzzzz_0, g_0_xyyyzz_0_yyyzzzz_1, g_0_xyyyzz_0_yyzzzzz_0, g_0_xyyyzz_0_yyzzzzz_1, g_0_xyyyzz_0_yzzzzzz_0, g_0_xyyyzz_0_yzzzzzz_1, g_0_xyyyzz_0_zzzzzzz_0, g_0_xyyyzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyzz_0_xxxxxxx_0[i] = 2.0 * g_0_xxxyzz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxyyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxxxxy_0[i] = g_0_xxxyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxxxxy_0[i] * pb_z + g_0_xxxyyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxxxxz_0[i] = 2.0 * g_0_xxxyzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxyyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxxxyy_0[i] = g_0_xxxyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxxxyy_0[i] * pb_z + g_0_xxxyyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxxxyz_0[i] = 2.0 * g_0_xyyyzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxyyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxxxyz_0[i] * pb_x + g_0_xxyyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxyzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxyyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxxyyy_0[i] = g_0_xxxyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxxyyy_0[i] * pb_z + g_0_xxxyyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxxyyz_0[i] = 2.0 * g_0_xyyyzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxxyyz_0[i] * pb_x + g_0_xxyyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxxyzz_0[i] = 2.0 * g_0_xyyyzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxxyzz_0[i] * pb_x + g_0_xxyyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxxzzz_0[i] = 2.0 * g_0_xxxyzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxyyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxyyyy_0[i] = g_0_xxxyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxyyyy_0[i] * pb_z + g_0_xxxyyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxyyyz_0[i] = 2.0 * g_0_xyyyzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxyyyz_0[i] * pb_x + g_0_xxyyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxyyzz_0[i] = 2.0 * g_0_xyyyzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxyyzz_0[i] * pb_x + g_0_xxyyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxyzzz_0[i] = 2.0 * g_0_xyyyzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxyzzz_0[i] * pb_x + g_0_xxyyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxzzzz_0[i] = 2.0 * g_0_xxxyzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxyyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxyyyyy_0[i] = g_0_xxxyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxyyyyy_0[i] * pb_z + g_0_xxxyyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxyyyyz_0[i] = 2.0 * g_0_xyyyzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyyyyz_0[i] * pb_x + g_0_xxyyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxyyyzz_0[i] = 2.0 * g_0_xyyyzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyyyzz_0[i] * pb_x + g_0_xxyyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxyyzzz_0[i] = 2.0 * g_0_xyyyzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyyzzz_0[i] * pb_x + g_0_xxyyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxyzzzz_0[i] = 2.0 * g_0_xyyyzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyzzzz_0[i] * pb_x + g_0_xxyyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxzzzzz_0[i] = 2.0 * g_0_xxxyzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxyyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xyyyyyy_0[i] = g_0_xxxyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xyyyyyy_0[i] * pb_z + g_0_xxxyyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xyyyyyz_0[i] = 2.0 * g_0_xyyyzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyyyyz_0[i] * pb_x + g_0_xxyyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyyyyzz_0[i] = 2.0 * g_0_xyyyzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyyyzz_0[i] * pb_x + g_0_xxyyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyyyzzz_0[i] = 2.0 * g_0_xyyyzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyyzzz_0[i] * pb_x + g_0_xxyyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyyzzzz_0[i] = 2.0 * g_0_xyyyzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyzzzz_0[i] * pb_x + g_0_xxyyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyzzzzz_0[i] = 2.0 * g_0_xyyyzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyzzzzz_0[i] * pb_x + g_0_xxyyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xzzzzzz_0[i] = 2.0 * g_0_xxxyzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxyyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_yyyyyyy_0[i] = 2.0 * g_0_xyyyzz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyyyyy_0[i] * pb_x + g_0_xxyyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyyyyz_0[i] = 2.0 * g_0_xyyyzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyyyyz_0[i] * pb_x + g_0_xxyyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyyyzz_0[i] = 2.0 * g_0_xyyyzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyyyzz_0[i] * pb_x + g_0_xxyyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyyzzz_0[i] = 2.0 * g_0_xyyyzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyyzzz_0[i] * pb_x + g_0_xxyyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyzzzz_0[i] = 2.0 * g_0_xyyyzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyzzzz_0[i] * pb_x + g_0_xxyyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyzzzzz_0[i] = 2.0 * g_0_xyyyzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyzzzzz_0[i] * pb_x + g_0_xxyyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yzzzzzz_0[i] = 2.0 * g_0_xyyyzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yzzzzzz_0[i] * pb_x + g_0_xxyyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_zzzzzzz_0[i] = 2.0 * g_0_xyyyzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_zzzzzzz_0[i] * pb_x + g_0_xxyyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 648-684 components of targeted buffer : SLSK

    auto g_0_xxxyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 648);

    auto g_0_xxxyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 649);

    auto g_0_xxxyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 650);

    auto g_0_xxxyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 651);

    auto g_0_xxxyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 652);

    auto g_0_xxxyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 653);

    auto g_0_xxxyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 654);

    auto g_0_xxxyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 655);

    auto g_0_xxxyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 656);

    auto g_0_xxxyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 657);

    auto g_0_xxxyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 658);

    auto g_0_xxxyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 659);

    auto g_0_xxxyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 660);

    auto g_0_xxxyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 661);

    auto g_0_xxxyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 662);

    auto g_0_xxxyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 663);

    auto g_0_xxxyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 664);

    auto g_0_xxxyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 665);

    auto g_0_xxxyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 666);

    auto g_0_xxxyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 667);

    auto g_0_xxxyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 668);

    auto g_0_xxxyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 669);

    auto g_0_xxxyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 670);

    auto g_0_xxxyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 671);

    auto g_0_xxxyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 672);

    auto g_0_xxxyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 673);

    auto g_0_xxxyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 674);

    auto g_0_xxxyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 675);

    auto g_0_xxxyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 676);

    auto g_0_xxxyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 677);

    auto g_0_xxxyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 678);

    auto g_0_xxxyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 679);

    auto g_0_xxxyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 680);

    auto g_0_xxxyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 681);

    auto g_0_xxxyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 682);

    auto g_0_xxxyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 683);

    #pragma omp simd aligned(g_0_xxxyyz_0_xxxxxxy_0, g_0_xxxyyz_0_xxxxxxy_1, g_0_xxxyyz_0_xxxxxyy_0, g_0_xxxyyz_0_xxxxxyy_1, g_0_xxxyyz_0_xxxxyyy_0, g_0_xxxyyz_0_xxxxyyy_1, g_0_xxxyyz_0_xxxyyyy_0, g_0_xxxyyz_0_xxxyyyy_1, g_0_xxxyyz_0_xxyyyyy_0, g_0_xxxyyz_0_xxyyyyy_1, g_0_xxxyyz_0_xyyyyyy_0, g_0_xxxyyz_0_xyyyyyy_1, g_0_xxxyyzz_0_xxxxxxy_0, g_0_xxxyyzz_0_xxxxxxy_1, g_0_xxxyyzz_0_xxxxxyy_0, g_0_xxxyyzz_0_xxxxxyy_1, g_0_xxxyyzz_0_xxxxyyy_0, g_0_xxxyyzz_0_xxxxyyy_1, g_0_xxxyyzz_0_xxxyyyy_0, g_0_xxxyyzz_0_xxxyyyy_1, g_0_xxxyyzz_0_xxyyyyy_0, g_0_xxxyyzz_0_xxyyyyy_1, g_0_xxxyyzz_0_xyyyyyy_0, g_0_xxxyyzz_0_xyyyyyy_1, g_0_xxxyyzzz_0_xxxxxxx_0, g_0_xxxyyzzz_0_xxxxxxy_0, g_0_xxxyyzzz_0_xxxxxxz_0, g_0_xxxyyzzz_0_xxxxxyy_0, g_0_xxxyyzzz_0_xxxxxyz_0, g_0_xxxyyzzz_0_xxxxxzz_0, g_0_xxxyyzzz_0_xxxxyyy_0, g_0_xxxyyzzz_0_xxxxyyz_0, g_0_xxxyyzzz_0_xxxxyzz_0, g_0_xxxyyzzz_0_xxxxzzz_0, g_0_xxxyyzzz_0_xxxyyyy_0, g_0_xxxyyzzz_0_xxxyyyz_0, g_0_xxxyyzzz_0_xxxyyzz_0, g_0_xxxyyzzz_0_xxxyzzz_0, g_0_xxxyyzzz_0_xxxzzzz_0, g_0_xxxyyzzz_0_xxyyyyy_0, g_0_xxxyyzzz_0_xxyyyyz_0, g_0_xxxyyzzz_0_xxyyyzz_0, g_0_xxxyyzzz_0_xxyyzzz_0, g_0_xxxyyzzz_0_xxyzzzz_0, g_0_xxxyyzzz_0_xxzzzzz_0, g_0_xxxyyzzz_0_xyyyyyy_0, g_0_xxxyyzzz_0_xyyyyyz_0, g_0_xxxyyzzz_0_xyyyyzz_0, g_0_xxxyyzzz_0_xyyyzzz_0, g_0_xxxyyzzz_0_xyyzzzz_0, g_0_xxxyyzzz_0_xyzzzzz_0, g_0_xxxyyzzz_0_xzzzzzz_0, g_0_xxxyyzzz_0_yyyyyyy_0, g_0_xxxyyzzz_0_yyyyyyz_0, g_0_xxxyyzzz_0_yyyyyzz_0, g_0_xxxyyzzz_0_yyyyzzz_0, g_0_xxxyyzzz_0_yyyzzzz_0, g_0_xxxyyzzz_0_yyzzzzz_0, g_0_xxxyyzzz_0_yzzzzzz_0, g_0_xxxyyzzz_0_zzzzzzz_0, g_0_xxxyzzz_0_xxxxxxx_0, g_0_xxxyzzz_0_xxxxxxx_1, g_0_xxxyzzz_0_xxxxxxz_0, g_0_xxxyzzz_0_xxxxxxz_1, g_0_xxxyzzz_0_xxxxxzz_0, g_0_xxxyzzz_0_xxxxxzz_1, g_0_xxxyzzz_0_xxxxzzz_0, g_0_xxxyzzz_0_xxxxzzz_1, g_0_xxxyzzz_0_xxxzzzz_0, g_0_xxxyzzz_0_xxxzzzz_1, g_0_xxxyzzz_0_xxzzzzz_0, g_0_xxxyzzz_0_xxzzzzz_1, g_0_xxxyzzz_0_xzzzzzz_0, g_0_xxxyzzz_0_xzzzzzz_1, g_0_xxxzzz_0_xxxxxxx_0, g_0_xxxzzz_0_xxxxxxx_1, g_0_xxxzzz_0_xxxxxxz_0, g_0_xxxzzz_0_xxxxxxz_1, g_0_xxxzzz_0_xxxxxzz_0, g_0_xxxzzz_0_xxxxxzz_1, g_0_xxxzzz_0_xxxxzzz_0, g_0_xxxzzz_0_xxxxzzz_1, g_0_xxxzzz_0_xxxzzzz_0, g_0_xxxzzz_0_xxxzzzz_1, g_0_xxxzzz_0_xxzzzzz_0, g_0_xxxzzz_0_xxzzzzz_1, g_0_xxxzzz_0_xzzzzzz_0, g_0_xxxzzz_0_xzzzzzz_1, g_0_xxyyzzz_0_xxxxxyz_0, g_0_xxyyzzz_0_xxxxxyz_1, g_0_xxyyzzz_0_xxxxyyz_0, g_0_xxyyzzz_0_xxxxyyz_1, g_0_xxyyzzz_0_xxxxyz_1, g_0_xxyyzzz_0_xxxxyzz_0, g_0_xxyyzzz_0_xxxxyzz_1, g_0_xxyyzzz_0_xxxyyyz_0, g_0_xxyyzzz_0_xxxyyyz_1, g_0_xxyyzzz_0_xxxyyz_1, g_0_xxyyzzz_0_xxxyyzz_0, g_0_xxyyzzz_0_xxxyyzz_1, g_0_xxyyzzz_0_xxxyzz_1, g_0_xxyyzzz_0_xxxyzzz_0, g_0_xxyyzzz_0_xxxyzzz_1, g_0_xxyyzzz_0_xxyyyyz_0, g_0_xxyyzzz_0_xxyyyyz_1, g_0_xxyyzzz_0_xxyyyz_1, g_0_xxyyzzz_0_xxyyyzz_0, g_0_xxyyzzz_0_xxyyyzz_1, g_0_xxyyzzz_0_xxyyzz_1, g_0_xxyyzzz_0_xxyyzzz_0, g_0_xxyyzzz_0_xxyyzzz_1, g_0_xxyyzzz_0_xxyzzz_1, g_0_xxyyzzz_0_xxyzzzz_0, g_0_xxyyzzz_0_xxyzzzz_1, g_0_xxyyzzz_0_xyyyyyz_0, g_0_xxyyzzz_0_xyyyyyz_1, g_0_xxyyzzz_0_xyyyyz_1, g_0_xxyyzzz_0_xyyyyzz_0, g_0_xxyyzzz_0_xyyyyzz_1, g_0_xxyyzzz_0_xyyyzz_1, g_0_xxyyzzz_0_xyyyzzz_0, g_0_xxyyzzz_0_xyyyzzz_1, g_0_xxyyzzz_0_xyyzzz_1, g_0_xxyyzzz_0_xyyzzzz_0, g_0_xxyyzzz_0_xyyzzzz_1, g_0_xxyyzzz_0_xyzzzz_1, g_0_xxyyzzz_0_xyzzzzz_0, g_0_xxyyzzz_0_xyzzzzz_1, g_0_xxyyzzz_0_yyyyyyy_0, g_0_xxyyzzz_0_yyyyyyy_1, g_0_xxyyzzz_0_yyyyyyz_0, g_0_xxyyzzz_0_yyyyyyz_1, g_0_xxyyzzz_0_yyyyyz_1, g_0_xxyyzzz_0_yyyyyzz_0, g_0_xxyyzzz_0_yyyyyzz_1, g_0_xxyyzzz_0_yyyyzz_1, g_0_xxyyzzz_0_yyyyzzz_0, g_0_xxyyzzz_0_yyyyzzz_1, g_0_xxyyzzz_0_yyyzzz_1, g_0_xxyyzzz_0_yyyzzzz_0, g_0_xxyyzzz_0_yyyzzzz_1, g_0_xxyyzzz_0_yyzzzz_1, g_0_xxyyzzz_0_yyzzzzz_0, g_0_xxyyzzz_0_yyzzzzz_1, g_0_xxyyzzz_0_yzzzzz_1, g_0_xxyyzzz_0_yzzzzzz_0, g_0_xxyyzzz_0_yzzzzzz_1, g_0_xxyyzzz_0_zzzzzzz_0, g_0_xxyyzzz_0_zzzzzzz_1, g_0_xyyzzz_0_xxxxxyz_0, g_0_xyyzzz_0_xxxxxyz_1, g_0_xyyzzz_0_xxxxyyz_0, g_0_xyyzzz_0_xxxxyyz_1, g_0_xyyzzz_0_xxxxyzz_0, g_0_xyyzzz_0_xxxxyzz_1, g_0_xyyzzz_0_xxxyyyz_0, g_0_xyyzzz_0_xxxyyyz_1, g_0_xyyzzz_0_xxxyyzz_0, g_0_xyyzzz_0_xxxyyzz_1, g_0_xyyzzz_0_xxxyzzz_0, g_0_xyyzzz_0_xxxyzzz_1, g_0_xyyzzz_0_xxyyyyz_0, g_0_xyyzzz_0_xxyyyyz_1, g_0_xyyzzz_0_xxyyyzz_0, g_0_xyyzzz_0_xxyyyzz_1, g_0_xyyzzz_0_xxyyzzz_0, g_0_xyyzzz_0_xxyyzzz_1, g_0_xyyzzz_0_xxyzzzz_0, g_0_xyyzzz_0_xxyzzzz_1, g_0_xyyzzz_0_xyyyyyz_0, g_0_xyyzzz_0_xyyyyyz_1, g_0_xyyzzz_0_xyyyyzz_0, g_0_xyyzzz_0_xyyyyzz_1, g_0_xyyzzz_0_xyyyzzz_0, g_0_xyyzzz_0_xyyyzzz_1, g_0_xyyzzz_0_xyyzzzz_0, g_0_xyyzzz_0_xyyzzzz_1, g_0_xyyzzz_0_xyzzzzz_0, g_0_xyyzzz_0_xyzzzzz_1, g_0_xyyzzz_0_yyyyyyy_0, g_0_xyyzzz_0_yyyyyyy_1, g_0_xyyzzz_0_yyyyyyz_0, g_0_xyyzzz_0_yyyyyyz_1, g_0_xyyzzz_0_yyyyyzz_0, g_0_xyyzzz_0_yyyyyzz_1, g_0_xyyzzz_0_yyyyzzz_0, g_0_xyyzzz_0_yyyyzzz_1, g_0_xyyzzz_0_yyyzzzz_0, g_0_xyyzzz_0_yyyzzzz_1, g_0_xyyzzz_0_yyzzzzz_0, g_0_xyyzzz_0_yyzzzzz_1, g_0_xyyzzz_0_yzzzzzz_0, g_0_xyyzzz_0_yzzzzzz_1, g_0_xyyzzz_0_zzzzzzz_0, g_0_xyyzzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzzz_0_xxxxxxx_0[i] = g_0_xxxzzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxyzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxxxxy_0[i] = 2.0 * g_0_xxxyyz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxxxy_0[i] * pb_z + g_0_xxxyyzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxxxxz_0[i] = g_0_xxxzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxyzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxyyz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxxyy_0[i] * pb_z + g_0_xxxyyzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxxxyz_0[i] = 2.0 * g_0_xyyzzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxyyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxxxyz_0[i] * pb_x + g_0_xxyyzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxxxzz_0[i] = g_0_xxxzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxyzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxxyyy_0[i] = 2.0 * g_0_xxxyyz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxyyy_0[i] * pb_z + g_0_xxxyyzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xyyzzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxxyyz_0[i] * pb_x + g_0_xxyyzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxxyzz_0[i] = 2.0 * g_0_xyyzzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxxyzz_0[i] * pb_x + g_0_xxyyzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxxzzz_0[i] = g_0_xxxzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxyzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxyyyy_0[i] = 2.0 * g_0_xxxyyz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxyyyy_0[i] * pb_z + g_0_xxxyyzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxyyyz_0[i] = 2.0 * g_0_xyyzzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxyyyz_0[i] * pb_x + g_0_xxyyzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xyyzzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxyyzz_0[i] * pb_x + g_0_xxyyzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxyzzz_0[i] = 2.0 * g_0_xyyzzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxyzzz_0[i] * pb_x + g_0_xxyyzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxzzzz_0[i] = g_0_xxxzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxyzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_xxxyyz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxyyyyy_0[i] * pb_z + g_0_xxxyyzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxyyyyz_0[i] = 2.0 * g_0_xyyzzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyyyyz_0[i] * pb_x + g_0_xxyyzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxyyyzz_0[i] = 2.0 * g_0_xyyzzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyyyzz_0[i] * pb_x + g_0_xxyyzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xyyzzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyyzzz_0[i] * pb_x + g_0_xxyyzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxyzzzz_0[i] = 2.0 * g_0_xyyzzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyzzzz_0[i] * pb_x + g_0_xxyyzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxzzzzz_0[i] = g_0_xxxzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxyzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xyyyyyy_0[i] = 2.0 * g_0_xxxyyz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xyyyyyy_0[i] * pb_z + g_0_xxxyyzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xyyyyyz_0[i] = 2.0 * g_0_xyyzzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyyyyz_0[i] * pb_x + g_0_xxyyzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyyyyzz_0[i] = 2.0 * g_0_xyyzzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyyyzz_0[i] * pb_x + g_0_xxyyzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyyyzzz_0[i] = 2.0 * g_0_xyyzzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyyzzz_0[i] * pb_x + g_0_xxyyzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xyyzzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyzzzz_0[i] * pb_x + g_0_xxyyzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyzzzzz_0[i] = 2.0 * g_0_xyyzzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyzzzzz_0[i] * pb_x + g_0_xxyyzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xzzzzzz_0[i] = g_0_xxxzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxyzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_yyyyyyy_0[i] = 2.0 * g_0_xyyzzz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyyyyy_0[i] * pb_x + g_0_xxyyzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyyyyz_0[i] = 2.0 * g_0_xyyzzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyyyyz_0[i] * pb_x + g_0_xxyyzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyyyzz_0[i] = 2.0 * g_0_xyyzzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyyyzz_0[i] * pb_x + g_0_xxyyzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyyzzz_0[i] = 2.0 * g_0_xyyzzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyyzzz_0[i] * pb_x + g_0_xxyyzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyzzzz_0[i] = 2.0 * g_0_xyyzzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyzzzz_0[i] * pb_x + g_0_xxyyzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xyyzzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyzzzzz_0[i] * pb_x + g_0_xxyyzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yzzzzzz_0[i] = 2.0 * g_0_xyyzzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yzzzzzz_0[i] * pb_x + g_0_xxyyzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_zzzzzzz_0[i] = 2.0 * g_0_xyyzzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_zzzzzzz_0[i] * pb_x + g_0_xxyyzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 684-720 components of targeted buffer : SLSK

    auto g_0_xxxyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 684);

    auto g_0_xxxyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 685);

    auto g_0_xxxyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 686);

    auto g_0_xxxyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 687);

    auto g_0_xxxyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 688);

    auto g_0_xxxyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 689);

    auto g_0_xxxyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 690);

    auto g_0_xxxyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 691);

    auto g_0_xxxyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 692);

    auto g_0_xxxyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 693);

    auto g_0_xxxyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 694);

    auto g_0_xxxyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 695);

    auto g_0_xxxyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 696);

    auto g_0_xxxyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 697);

    auto g_0_xxxyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 698);

    auto g_0_xxxyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 699);

    auto g_0_xxxyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 700);

    auto g_0_xxxyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 701);

    auto g_0_xxxyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 702);

    auto g_0_xxxyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 703);

    auto g_0_xxxyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 704);

    auto g_0_xxxyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 705);

    auto g_0_xxxyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 706);

    auto g_0_xxxyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 707);

    auto g_0_xxxyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 708);

    auto g_0_xxxyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 709);

    auto g_0_xxxyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 710);

    auto g_0_xxxyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 711);

    auto g_0_xxxyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 712);

    auto g_0_xxxyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 713);

    auto g_0_xxxyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 714);

    auto g_0_xxxyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 715);

    auto g_0_xxxyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 716);

    auto g_0_xxxyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 717);

    auto g_0_xxxyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 718);

    auto g_0_xxxyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 719);

    #pragma omp simd aligned(g_0_xxxyzzzz_0_xxxxxxx_0, g_0_xxxyzzzz_0_xxxxxxy_0, g_0_xxxyzzzz_0_xxxxxxz_0, g_0_xxxyzzzz_0_xxxxxyy_0, g_0_xxxyzzzz_0_xxxxxyz_0, g_0_xxxyzzzz_0_xxxxxzz_0, g_0_xxxyzzzz_0_xxxxyyy_0, g_0_xxxyzzzz_0_xxxxyyz_0, g_0_xxxyzzzz_0_xxxxyzz_0, g_0_xxxyzzzz_0_xxxxzzz_0, g_0_xxxyzzzz_0_xxxyyyy_0, g_0_xxxyzzzz_0_xxxyyyz_0, g_0_xxxyzzzz_0_xxxyyzz_0, g_0_xxxyzzzz_0_xxxyzzz_0, g_0_xxxyzzzz_0_xxxzzzz_0, g_0_xxxyzzzz_0_xxyyyyy_0, g_0_xxxyzzzz_0_xxyyyyz_0, g_0_xxxyzzzz_0_xxyyyzz_0, g_0_xxxyzzzz_0_xxyyzzz_0, g_0_xxxyzzzz_0_xxyzzzz_0, g_0_xxxyzzzz_0_xxzzzzz_0, g_0_xxxyzzzz_0_xyyyyyy_0, g_0_xxxyzzzz_0_xyyyyyz_0, g_0_xxxyzzzz_0_xyyyyzz_0, g_0_xxxyzzzz_0_xyyyzzz_0, g_0_xxxyzzzz_0_xyyzzzz_0, g_0_xxxyzzzz_0_xyzzzzz_0, g_0_xxxyzzzz_0_xzzzzzz_0, g_0_xxxyzzzz_0_yyyyyyy_0, g_0_xxxyzzzz_0_yyyyyyz_0, g_0_xxxyzzzz_0_yyyyyzz_0, g_0_xxxyzzzz_0_yyyyzzz_0, g_0_xxxyzzzz_0_yyyzzzz_0, g_0_xxxyzzzz_0_yyzzzzz_0, g_0_xxxyzzzz_0_yzzzzzz_0, g_0_xxxyzzzz_0_zzzzzzz_0, g_0_xxxzzzz_0_xxxxxx_1, g_0_xxxzzzz_0_xxxxxxx_0, g_0_xxxzzzz_0_xxxxxxx_1, g_0_xxxzzzz_0_xxxxxxy_0, g_0_xxxzzzz_0_xxxxxxy_1, g_0_xxxzzzz_0_xxxxxxz_0, g_0_xxxzzzz_0_xxxxxxz_1, g_0_xxxzzzz_0_xxxxxy_1, g_0_xxxzzzz_0_xxxxxyy_0, g_0_xxxzzzz_0_xxxxxyy_1, g_0_xxxzzzz_0_xxxxxyz_0, g_0_xxxzzzz_0_xxxxxyz_1, g_0_xxxzzzz_0_xxxxxz_1, g_0_xxxzzzz_0_xxxxxzz_0, g_0_xxxzzzz_0_xxxxxzz_1, g_0_xxxzzzz_0_xxxxyy_1, g_0_xxxzzzz_0_xxxxyyy_0, g_0_xxxzzzz_0_xxxxyyy_1, g_0_xxxzzzz_0_xxxxyyz_0, g_0_xxxzzzz_0_xxxxyyz_1, g_0_xxxzzzz_0_xxxxyz_1, g_0_xxxzzzz_0_xxxxyzz_0, g_0_xxxzzzz_0_xxxxyzz_1, g_0_xxxzzzz_0_xxxxzz_1, g_0_xxxzzzz_0_xxxxzzz_0, g_0_xxxzzzz_0_xxxxzzz_1, g_0_xxxzzzz_0_xxxyyy_1, g_0_xxxzzzz_0_xxxyyyy_0, g_0_xxxzzzz_0_xxxyyyy_1, g_0_xxxzzzz_0_xxxyyyz_0, g_0_xxxzzzz_0_xxxyyyz_1, g_0_xxxzzzz_0_xxxyyz_1, g_0_xxxzzzz_0_xxxyyzz_0, g_0_xxxzzzz_0_xxxyyzz_1, g_0_xxxzzzz_0_xxxyzz_1, g_0_xxxzzzz_0_xxxyzzz_0, g_0_xxxzzzz_0_xxxyzzz_1, g_0_xxxzzzz_0_xxxzzz_1, g_0_xxxzzzz_0_xxxzzzz_0, g_0_xxxzzzz_0_xxxzzzz_1, g_0_xxxzzzz_0_xxyyyy_1, g_0_xxxzzzz_0_xxyyyyy_0, g_0_xxxzzzz_0_xxyyyyy_1, g_0_xxxzzzz_0_xxyyyyz_0, g_0_xxxzzzz_0_xxyyyyz_1, g_0_xxxzzzz_0_xxyyyz_1, g_0_xxxzzzz_0_xxyyyzz_0, g_0_xxxzzzz_0_xxyyyzz_1, g_0_xxxzzzz_0_xxyyzz_1, g_0_xxxzzzz_0_xxyyzzz_0, g_0_xxxzzzz_0_xxyyzzz_1, g_0_xxxzzzz_0_xxyzzz_1, g_0_xxxzzzz_0_xxyzzzz_0, g_0_xxxzzzz_0_xxyzzzz_1, g_0_xxxzzzz_0_xxzzzz_1, g_0_xxxzzzz_0_xxzzzzz_0, g_0_xxxzzzz_0_xxzzzzz_1, g_0_xxxzzzz_0_xyyyyy_1, g_0_xxxzzzz_0_xyyyyyy_0, g_0_xxxzzzz_0_xyyyyyy_1, g_0_xxxzzzz_0_xyyyyyz_0, g_0_xxxzzzz_0_xyyyyyz_1, g_0_xxxzzzz_0_xyyyyz_1, g_0_xxxzzzz_0_xyyyyzz_0, g_0_xxxzzzz_0_xyyyyzz_1, g_0_xxxzzzz_0_xyyyzz_1, g_0_xxxzzzz_0_xyyyzzz_0, g_0_xxxzzzz_0_xyyyzzz_1, g_0_xxxzzzz_0_xyyzzz_1, g_0_xxxzzzz_0_xyyzzzz_0, g_0_xxxzzzz_0_xyyzzzz_1, g_0_xxxzzzz_0_xyzzzz_1, g_0_xxxzzzz_0_xyzzzzz_0, g_0_xxxzzzz_0_xyzzzzz_1, g_0_xxxzzzz_0_xzzzzz_1, g_0_xxxzzzz_0_xzzzzzz_0, g_0_xxxzzzz_0_xzzzzzz_1, g_0_xxxzzzz_0_yyyyyy_1, g_0_xxxzzzz_0_yyyyyyy_0, g_0_xxxzzzz_0_yyyyyyy_1, g_0_xxxzzzz_0_yyyyyyz_0, g_0_xxxzzzz_0_yyyyyyz_1, g_0_xxxzzzz_0_yyyyyz_1, g_0_xxxzzzz_0_yyyyyzz_0, g_0_xxxzzzz_0_yyyyyzz_1, g_0_xxxzzzz_0_yyyyzz_1, g_0_xxxzzzz_0_yyyyzzz_0, g_0_xxxzzzz_0_yyyyzzz_1, g_0_xxxzzzz_0_yyyzzz_1, g_0_xxxzzzz_0_yyyzzzz_0, g_0_xxxzzzz_0_yyyzzzz_1, g_0_xxxzzzz_0_yyzzzz_1, g_0_xxxzzzz_0_yyzzzzz_0, g_0_xxxzzzz_0_yyzzzzz_1, g_0_xxxzzzz_0_yzzzzz_1, g_0_xxxzzzz_0_yzzzzzz_0, g_0_xxxzzzz_0_yzzzzzz_1, g_0_xxxzzzz_0_zzzzzz_1, g_0_xxxzzzz_0_zzzzzzz_0, g_0_xxxzzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzzz_0_xxxxxxx_0[i] = g_0_xxxzzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxxxy_0[i] = g_0_xxxzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxxxy_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxxxz_0[i] = g_0_xxxzzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxxyy_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxxyz_0[i] = g_0_xxxzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxxyz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxxzz_0[i] = g_0_xxxzzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxxzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxyyy_0[i] * pb_y + g_0_xxxzzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxyyz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxyzz_0[i] = g_0_xxxzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxyzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxzzz_0[i] = g_0_xxxzzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxxzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyyyy_0[i] * pb_y + g_0_xxxzzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyyyz_0[i] * pb_y + g_0_xxxzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyyzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyzzz_0[i] = g_0_xxxzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxzzzz_0[i] = g_0_xxxzzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxxzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyyyy_0[i] * pb_y + g_0_xxxzzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyyyz_0[i] * pb_y + g_0_xxxzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyyzz_0[i] * pb_y + g_0_xxxzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyzzzz_0[i] = g_0_xxxzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxzzzzz_0[i] = g_0_xxxzzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxxzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyyyy_0[i] * pb_y + g_0_xxxzzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyyyz_0[i] * pb_y + g_0_xxxzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyyzz_0[i] * pb_y + g_0_xxxzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyzzz_0[i] * pb_y + g_0_xxxzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyzzzzz_0[i] = g_0_xxxzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xzzzzzz_0[i] = g_0_xxxzzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxxzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyyyy_0[i] * pb_y + g_0_xxxzzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyyyz_0[i] * pb_y + g_0_xxxzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyyzz_0[i] * pb_y + g_0_xxxzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyzzz_0[i] * pb_y + g_0_xxxzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyzzzz_0[i] * pb_y + g_0_xxxzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yzzzzzz_0[i] = g_0_xxxzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yzzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_zzzzzzz_0[i] = g_0_xxxzzzz_0_zzzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 720-756 components of targeted buffer : SLSK

    auto g_0_xxxzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 720);

    auto g_0_xxxzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 721);

    auto g_0_xxxzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 722);

    auto g_0_xxxzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 723);

    auto g_0_xxxzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 724);

    auto g_0_xxxzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 725);

    auto g_0_xxxzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 726);

    auto g_0_xxxzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 727);

    auto g_0_xxxzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 728);

    auto g_0_xxxzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 729);

    auto g_0_xxxzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 730);

    auto g_0_xxxzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 731);

    auto g_0_xxxzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 732);

    auto g_0_xxxzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 733);

    auto g_0_xxxzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 734);

    auto g_0_xxxzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 735);

    auto g_0_xxxzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 736);

    auto g_0_xxxzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 737);

    auto g_0_xxxzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 738);

    auto g_0_xxxzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 739);

    auto g_0_xxxzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 740);

    auto g_0_xxxzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 741);

    auto g_0_xxxzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 742);

    auto g_0_xxxzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 743);

    auto g_0_xxxzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 744);

    auto g_0_xxxzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 745);

    auto g_0_xxxzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 746);

    auto g_0_xxxzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 747);

    auto g_0_xxxzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 748);

    auto g_0_xxxzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 749);

    auto g_0_xxxzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 750);

    auto g_0_xxxzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 751);

    auto g_0_xxxzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 752);

    auto g_0_xxxzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 753);

    auto g_0_xxxzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 754);

    auto g_0_xxxzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 755);

    #pragma omp simd aligned(g_0_xxxzzz_0_xxxxxxx_0, g_0_xxxzzz_0_xxxxxxx_1, g_0_xxxzzz_0_xxxxxxy_0, g_0_xxxzzz_0_xxxxxxy_1, g_0_xxxzzz_0_xxxxxyy_0, g_0_xxxzzz_0_xxxxxyy_1, g_0_xxxzzz_0_xxxxyyy_0, g_0_xxxzzz_0_xxxxyyy_1, g_0_xxxzzz_0_xxxyyyy_0, g_0_xxxzzz_0_xxxyyyy_1, g_0_xxxzzz_0_xxyyyyy_0, g_0_xxxzzz_0_xxyyyyy_1, g_0_xxxzzz_0_xyyyyyy_0, g_0_xxxzzz_0_xyyyyyy_1, g_0_xxxzzzz_0_xxxxxxx_0, g_0_xxxzzzz_0_xxxxxxx_1, g_0_xxxzzzz_0_xxxxxxy_0, g_0_xxxzzzz_0_xxxxxxy_1, g_0_xxxzzzz_0_xxxxxyy_0, g_0_xxxzzzz_0_xxxxxyy_1, g_0_xxxzzzz_0_xxxxyyy_0, g_0_xxxzzzz_0_xxxxyyy_1, g_0_xxxzzzz_0_xxxyyyy_0, g_0_xxxzzzz_0_xxxyyyy_1, g_0_xxxzzzz_0_xxyyyyy_0, g_0_xxxzzzz_0_xxyyyyy_1, g_0_xxxzzzz_0_xyyyyyy_0, g_0_xxxzzzz_0_xyyyyyy_1, g_0_xxxzzzzz_0_xxxxxxx_0, g_0_xxxzzzzz_0_xxxxxxy_0, g_0_xxxzzzzz_0_xxxxxxz_0, g_0_xxxzzzzz_0_xxxxxyy_0, g_0_xxxzzzzz_0_xxxxxyz_0, g_0_xxxzzzzz_0_xxxxxzz_0, g_0_xxxzzzzz_0_xxxxyyy_0, g_0_xxxzzzzz_0_xxxxyyz_0, g_0_xxxzzzzz_0_xxxxyzz_0, g_0_xxxzzzzz_0_xxxxzzz_0, g_0_xxxzzzzz_0_xxxyyyy_0, g_0_xxxzzzzz_0_xxxyyyz_0, g_0_xxxzzzzz_0_xxxyyzz_0, g_0_xxxzzzzz_0_xxxyzzz_0, g_0_xxxzzzzz_0_xxxzzzz_0, g_0_xxxzzzzz_0_xxyyyyy_0, g_0_xxxzzzzz_0_xxyyyyz_0, g_0_xxxzzzzz_0_xxyyyzz_0, g_0_xxxzzzzz_0_xxyyzzz_0, g_0_xxxzzzzz_0_xxyzzzz_0, g_0_xxxzzzzz_0_xxzzzzz_0, g_0_xxxzzzzz_0_xyyyyyy_0, g_0_xxxzzzzz_0_xyyyyyz_0, g_0_xxxzzzzz_0_xyyyyzz_0, g_0_xxxzzzzz_0_xyyyzzz_0, g_0_xxxzzzzz_0_xyyzzzz_0, g_0_xxxzzzzz_0_xyzzzzz_0, g_0_xxxzzzzz_0_xzzzzzz_0, g_0_xxxzzzzz_0_yyyyyyy_0, g_0_xxxzzzzz_0_yyyyyyz_0, g_0_xxxzzzzz_0_yyyyyzz_0, g_0_xxxzzzzz_0_yyyyzzz_0, g_0_xxxzzzzz_0_yyyzzzz_0, g_0_xxxzzzzz_0_yyzzzzz_0, g_0_xxxzzzzz_0_yzzzzzz_0, g_0_xxxzzzzz_0_zzzzzzz_0, g_0_xxzzzzz_0_xxxxxxz_0, g_0_xxzzzzz_0_xxxxxxz_1, g_0_xxzzzzz_0_xxxxxyz_0, g_0_xxzzzzz_0_xxxxxyz_1, g_0_xxzzzzz_0_xxxxxz_1, g_0_xxzzzzz_0_xxxxxzz_0, g_0_xxzzzzz_0_xxxxxzz_1, g_0_xxzzzzz_0_xxxxyyz_0, g_0_xxzzzzz_0_xxxxyyz_1, g_0_xxzzzzz_0_xxxxyz_1, g_0_xxzzzzz_0_xxxxyzz_0, g_0_xxzzzzz_0_xxxxyzz_1, g_0_xxzzzzz_0_xxxxzz_1, g_0_xxzzzzz_0_xxxxzzz_0, g_0_xxzzzzz_0_xxxxzzz_1, g_0_xxzzzzz_0_xxxyyyz_0, g_0_xxzzzzz_0_xxxyyyz_1, g_0_xxzzzzz_0_xxxyyz_1, g_0_xxzzzzz_0_xxxyyzz_0, g_0_xxzzzzz_0_xxxyyzz_1, g_0_xxzzzzz_0_xxxyzz_1, g_0_xxzzzzz_0_xxxyzzz_0, g_0_xxzzzzz_0_xxxyzzz_1, g_0_xxzzzzz_0_xxxzzz_1, g_0_xxzzzzz_0_xxxzzzz_0, g_0_xxzzzzz_0_xxxzzzz_1, g_0_xxzzzzz_0_xxyyyyz_0, g_0_xxzzzzz_0_xxyyyyz_1, g_0_xxzzzzz_0_xxyyyz_1, g_0_xxzzzzz_0_xxyyyzz_0, g_0_xxzzzzz_0_xxyyyzz_1, g_0_xxzzzzz_0_xxyyzz_1, g_0_xxzzzzz_0_xxyyzzz_0, g_0_xxzzzzz_0_xxyyzzz_1, g_0_xxzzzzz_0_xxyzzz_1, g_0_xxzzzzz_0_xxyzzzz_0, g_0_xxzzzzz_0_xxyzzzz_1, g_0_xxzzzzz_0_xxzzzz_1, g_0_xxzzzzz_0_xxzzzzz_0, g_0_xxzzzzz_0_xxzzzzz_1, g_0_xxzzzzz_0_xyyyyyz_0, g_0_xxzzzzz_0_xyyyyyz_1, g_0_xxzzzzz_0_xyyyyz_1, g_0_xxzzzzz_0_xyyyyzz_0, g_0_xxzzzzz_0_xyyyyzz_1, g_0_xxzzzzz_0_xyyyzz_1, g_0_xxzzzzz_0_xyyyzzz_0, g_0_xxzzzzz_0_xyyyzzz_1, g_0_xxzzzzz_0_xyyzzz_1, g_0_xxzzzzz_0_xyyzzzz_0, g_0_xxzzzzz_0_xyyzzzz_1, g_0_xxzzzzz_0_xyzzzz_1, g_0_xxzzzzz_0_xyzzzzz_0, g_0_xxzzzzz_0_xyzzzzz_1, g_0_xxzzzzz_0_xzzzzz_1, g_0_xxzzzzz_0_xzzzzzz_0, g_0_xxzzzzz_0_xzzzzzz_1, g_0_xxzzzzz_0_yyyyyyy_0, g_0_xxzzzzz_0_yyyyyyy_1, g_0_xxzzzzz_0_yyyyyyz_0, g_0_xxzzzzz_0_yyyyyyz_1, g_0_xxzzzzz_0_yyyyyz_1, g_0_xxzzzzz_0_yyyyyzz_0, g_0_xxzzzzz_0_yyyyyzz_1, g_0_xxzzzzz_0_yyyyzz_1, g_0_xxzzzzz_0_yyyyzzz_0, g_0_xxzzzzz_0_yyyyzzz_1, g_0_xxzzzzz_0_yyyzzz_1, g_0_xxzzzzz_0_yyyzzzz_0, g_0_xxzzzzz_0_yyyzzzz_1, g_0_xxzzzzz_0_yyzzzz_1, g_0_xxzzzzz_0_yyzzzzz_0, g_0_xxzzzzz_0_yyzzzzz_1, g_0_xxzzzzz_0_yzzzzz_1, g_0_xxzzzzz_0_yzzzzzz_0, g_0_xxzzzzz_0_yzzzzzz_1, g_0_xxzzzzz_0_zzzzzz_1, g_0_xxzzzzz_0_zzzzzzz_0, g_0_xxzzzzz_0_zzzzzzz_1, g_0_xzzzzz_0_xxxxxxz_0, g_0_xzzzzz_0_xxxxxxz_1, g_0_xzzzzz_0_xxxxxyz_0, g_0_xzzzzz_0_xxxxxyz_1, g_0_xzzzzz_0_xxxxxzz_0, g_0_xzzzzz_0_xxxxxzz_1, g_0_xzzzzz_0_xxxxyyz_0, g_0_xzzzzz_0_xxxxyyz_1, g_0_xzzzzz_0_xxxxyzz_0, g_0_xzzzzz_0_xxxxyzz_1, g_0_xzzzzz_0_xxxxzzz_0, g_0_xzzzzz_0_xxxxzzz_1, g_0_xzzzzz_0_xxxyyyz_0, g_0_xzzzzz_0_xxxyyyz_1, g_0_xzzzzz_0_xxxyyzz_0, g_0_xzzzzz_0_xxxyyzz_1, g_0_xzzzzz_0_xxxyzzz_0, g_0_xzzzzz_0_xxxyzzz_1, g_0_xzzzzz_0_xxxzzzz_0, g_0_xzzzzz_0_xxxzzzz_1, g_0_xzzzzz_0_xxyyyyz_0, g_0_xzzzzz_0_xxyyyyz_1, g_0_xzzzzz_0_xxyyyzz_0, g_0_xzzzzz_0_xxyyyzz_1, g_0_xzzzzz_0_xxyyzzz_0, g_0_xzzzzz_0_xxyyzzz_1, g_0_xzzzzz_0_xxyzzzz_0, g_0_xzzzzz_0_xxyzzzz_1, g_0_xzzzzz_0_xxzzzzz_0, g_0_xzzzzz_0_xxzzzzz_1, g_0_xzzzzz_0_xyyyyyz_0, g_0_xzzzzz_0_xyyyyyz_1, g_0_xzzzzz_0_xyyyyzz_0, g_0_xzzzzz_0_xyyyyzz_1, g_0_xzzzzz_0_xyyyzzz_0, g_0_xzzzzz_0_xyyyzzz_1, g_0_xzzzzz_0_xyyzzzz_0, g_0_xzzzzz_0_xyyzzzz_1, g_0_xzzzzz_0_xyzzzzz_0, g_0_xzzzzz_0_xyzzzzz_1, g_0_xzzzzz_0_xzzzzzz_0, g_0_xzzzzz_0_xzzzzzz_1, g_0_xzzzzz_0_yyyyyyy_0, g_0_xzzzzz_0_yyyyyyy_1, g_0_xzzzzz_0_yyyyyyz_0, g_0_xzzzzz_0_yyyyyyz_1, g_0_xzzzzz_0_yyyyyzz_0, g_0_xzzzzz_0_yyyyyzz_1, g_0_xzzzzz_0_yyyyzzz_0, g_0_xzzzzz_0_yyyyzzz_1, g_0_xzzzzz_0_yyyzzzz_0, g_0_xzzzzz_0_yyyzzzz_1, g_0_xzzzzz_0_yyzzzzz_0, g_0_xzzzzz_0_yyzzzzz_1, g_0_xzzzzz_0_yzzzzzz_0, g_0_xzzzzz_0_yzzzzzz_1, g_0_xzzzzz_0_zzzzzzz_0, g_0_xzzzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzzz_0_xxxxxxx_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxxxxx_0[i] * pb_z + g_0_xxxzzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxxxy_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxxxxy_0[i] * pb_z + g_0_xxxzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxxxz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxxxz_0[i] * pb_x + g_0_xxzzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxxxyy_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxxxyy_0[i] * pb_z + g_0_xxxzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxxyz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxxyz_0[i] * pb_x + g_0_xxzzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxxxzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxxzz_0[i] * pb_x + g_0_xxzzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_xxxzzz_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxxyyy_0[i] * pb_z + g_0_xxxzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxyyz_0[i] * pb_x + g_0_xxzzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxxyzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxyzz_0[i] * pb_x + g_0_xxzzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxxzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxzzz_0[i] * pb_x + g_0_xxzzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxyyyy_0[i] * pb_z + g_0_xxxzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxyyyz_0[i] = 2.0 * g_0_xzzzzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyyyz_0[i] * pb_x + g_0_xxzzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyyzz_0[i] * pb_x + g_0_xxzzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxyzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyzzz_0[i] * pb_x + g_0_xxzzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxyyyyy_0[i] * pb_z + g_0_xxxzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_xzzzzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyyyz_0[i] * pb_x + g_0_xxzzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_xzzzzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyyzz_0[i] * pb_x + g_0_xxzzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyzzz_0[i] * pb_x + g_0_xxzzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xyyyyyy_0[i] * pb_z + g_0_xxxzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xyyyyyz_0[i] = 2.0 * g_0_xzzzzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyyyz_0[i] * pb_x + g_0_xxzzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyyyzz_0[i] = 2.0 * g_0_xzzzzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyyzz_0[i] * pb_x + g_0_xxzzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyyzzz_0[i] = 2.0 * g_0_xzzzzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyzzz_0[i] * pb_x + g_0_xxzzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xzzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xzzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyyyy_0[i] = 2.0 * g_0_xzzzzz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyyyyy_0[i] * pb_x + g_0_xxzzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyyyz_0[i] = 2.0 * g_0_xzzzzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyyyyz_0[i] * pb_x + g_0_xxzzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyyzz_0[i] = 2.0 * g_0_xzzzzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyyyzz_0[i] * pb_x + g_0_xxzzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyzzz_0[i] = 2.0 * g_0_xzzzzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyyzzz_0[i] * pb_x + g_0_xxzzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyzzzz_0[i] = 2.0 * g_0_xzzzzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyzzzz_0[i] * pb_x + g_0_xxzzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yzzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yzzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_zzzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zzzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 756-792 components of targeted buffer : SLSK

    auto g_0_xxyyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 756);

    auto g_0_xxyyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 757);

    auto g_0_xxyyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 758);

    auto g_0_xxyyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 759);

    auto g_0_xxyyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 760);

    auto g_0_xxyyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 761);

    auto g_0_xxyyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 762);

    auto g_0_xxyyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 763);

    auto g_0_xxyyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 764);

    auto g_0_xxyyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 765);

    auto g_0_xxyyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 766);

    auto g_0_xxyyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 767);

    auto g_0_xxyyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 768);

    auto g_0_xxyyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 769);

    auto g_0_xxyyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 770);

    auto g_0_xxyyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 771);

    auto g_0_xxyyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 772);

    auto g_0_xxyyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 773);

    auto g_0_xxyyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 774);

    auto g_0_xxyyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 775);

    auto g_0_xxyyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 776);

    auto g_0_xxyyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 777);

    auto g_0_xxyyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 778);

    auto g_0_xxyyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 779);

    auto g_0_xxyyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 780);

    auto g_0_xxyyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 781);

    auto g_0_xxyyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 782);

    auto g_0_xxyyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 783);

    auto g_0_xxyyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 784);

    auto g_0_xxyyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 785);

    auto g_0_xxyyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 786);

    auto g_0_xxyyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 787);

    auto g_0_xxyyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 788);

    auto g_0_xxyyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 789);

    auto g_0_xxyyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 790);

    auto g_0_xxyyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 791);

    #pragma omp simd aligned(g_0_xxyyyy_0_xxxxxxx_0, g_0_xxyyyy_0_xxxxxxx_1, g_0_xxyyyy_0_xxxxxxz_0, g_0_xxyyyy_0_xxxxxxz_1, g_0_xxyyyy_0_xxxxxzz_0, g_0_xxyyyy_0_xxxxxzz_1, g_0_xxyyyy_0_xxxxzzz_0, g_0_xxyyyy_0_xxxxzzz_1, g_0_xxyyyy_0_xxxzzzz_0, g_0_xxyyyy_0_xxxzzzz_1, g_0_xxyyyy_0_xxzzzzz_0, g_0_xxyyyy_0_xxzzzzz_1, g_0_xxyyyy_0_xzzzzzz_0, g_0_xxyyyy_0_xzzzzzz_1, g_0_xxyyyyy_0_xxxxxxx_0, g_0_xxyyyyy_0_xxxxxxx_1, g_0_xxyyyyy_0_xxxxxxz_0, g_0_xxyyyyy_0_xxxxxxz_1, g_0_xxyyyyy_0_xxxxxzz_0, g_0_xxyyyyy_0_xxxxxzz_1, g_0_xxyyyyy_0_xxxxzzz_0, g_0_xxyyyyy_0_xxxxzzz_1, g_0_xxyyyyy_0_xxxzzzz_0, g_0_xxyyyyy_0_xxxzzzz_1, g_0_xxyyyyy_0_xxzzzzz_0, g_0_xxyyyyy_0_xxzzzzz_1, g_0_xxyyyyy_0_xzzzzzz_0, g_0_xxyyyyy_0_xzzzzzz_1, g_0_xxyyyyyy_0_xxxxxxx_0, g_0_xxyyyyyy_0_xxxxxxy_0, g_0_xxyyyyyy_0_xxxxxxz_0, g_0_xxyyyyyy_0_xxxxxyy_0, g_0_xxyyyyyy_0_xxxxxyz_0, g_0_xxyyyyyy_0_xxxxxzz_0, g_0_xxyyyyyy_0_xxxxyyy_0, g_0_xxyyyyyy_0_xxxxyyz_0, g_0_xxyyyyyy_0_xxxxyzz_0, g_0_xxyyyyyy_0_xxxxzzz_0, g_0_xxyyyyyy_0_xxxyyyy_0, g_0_xxyyyyyy_0_xxxyyyz_0, g_0_xxyyyyyy_0_xxxyyzz_0, g_0_xxyyyyyy_0_xxxyzzz_0, g_0_xxyyyyyy_0_xxxzzzz_0, g_0_xxyyyyyy_0_xxyyyyy_0, g_0_xxyyyyyy_0_xxyyyyz_0, g_0_xxyyyyyy_0_xxyyyzz_0, g_0_xxyyyyyy_0_xxyyzzz_0, g_0_xxyyyyyy_0_xxyzzzz_0, g_0_xxyyyyyy_0_xxzzzzz_0, g_0_xxyyyyyy_0_xyyyyyy_0, g_0_xxyyyyyy_0_xyyyyyz_0, g_0_xxyyyyyy_0_xyyyyzz_0, g_0_xxyyyyyy_0_xyyyzzz_0, g_0_xxyyyyyy_0_xyyzzzz_0, g_0_xxyyyyyy_0_xyzzzzz_0, g_0_xxyyyyyy_0_xzzzzzz_0, g_0_xxyyyyyy_0_yyyyyyy_0, g_0_xxyyyyyy_0_yyyyyyz_0, g_0_xxyyyyyy_0_yyyyyzz_0, g_0_xxyyyyyy_0_yyyyzzz_0, g_0_xxyyyyyy_0_yyyzzzz_0, g_0_xxyyyyyy_0_yyzzzzz_0, g_0_xxyyyyyy_0_yzzzzzz_0, g_0_xxyyyyyy_0_zzzzzzz_0, g_0_xyyyyyy_0_xxxxxxy_0, g_0_xyyyyyy_0_xxxxxxy_1, g_0_xyyyyyy_0_xxxxxy_1, g_0_xyyyyyy_0_xxxxxyy_0, g_0_xyyyyyy_0_xxxxxyy_1, g_0_xyyyyyy_0_xxxxxyz_0, g_0_xyyyyyy_0_xxxxxyz_1, g_0_xyyyyyy_0_xxxxyy_1, g_0_xyyyyyy_0_xxxxyyy_0, g_0_xyyyyyy_0_xxxxyyy_1, g_0_xyyyyyy_0_xxxxyyz_0, g_0_xyyyyyy_0_xxxxyyz_1, g_0_xyyyyyy_0_xxxxyz_1, g_0_xyyyyyy_0_xxxxyzz_0, g_0_xyyyyyy_0_xxxxyzz_1, g_0_xyyyyyy_0_xxxyyy_1, g_0_xyyyyyy_0_xxxyyyy_0, g_0_xyyyyyy_0_xxxyyyy_1, g_0_xyyyyyy_0_xxxyyyz_0, g_0_xyyyyyy_0_xxxyyyz_1, g_0_xyyyyyy_0_xxxyyz_1, g_0_xyyyyyy_0_xxxyyzz_0, g_0_xyyyyyy_0_xxxyyzz_1, g_0_xyyyyyy_0_xxxyzz_1, g_0_xyyyyyy_0_xxxyzzz_0, g_0_xyyyyyy_0_xxxyzzz_1, g_0_xyyyyyy_0_xxyyyy_1, g_0_xyyyyyy_0_xxyyyyy_0, g_0_xyyyyyy_0_xxyyyyy_1, g_0_xyyyyyy_0_xxyyyyz_0, g_0_xyyyyyy_0_xxyyyyz_1, g_0_xyyyyyy_0_xxyyyz_1, g_0_xyyyyyy_0_xxyyyzz_0, g_0_xyyyyyy_0_xxyyyzz_1, g_0_xyyyyyy_0_xxyyzz_1, g_0_xyyyyyy_0_xxyyzzz_0, g_0_xyyyyyy_0_xxyyzzz_1, g_0_xyyyyyy_0_xxyzzz_1, g_0_xyyyyyy_0_xxyzzzz_0, g_0_xyyyyyy_0_xxyzzzz_1, g_0_xyyyyyy_0_xyyyyy_1, g_0_xyyyyyy_0_xyyyyyy_0, g_0_xyyyyyy_0_xyyyyyy_1, g_0_xyyyyyy_0_xyyyyyz_0, g_0_xyyyyyy_0_xyyyyyz_1, g_0_xyyyyyy_0_xyyyyz_1, g_0_xyyyyyy_0_xyyyyzz_0, g_0_xyyyyyy_0_xyyyyzz_1, g_0_xyyyyyy_0_xyyyzz_1, g_0_xyyyyyy_0_xyyyzzz_0, g_0_xyyyyyy_0_xyyyzzz_1, g_0_xyyyyyy_0_xyyzzz_1, g_0_xyyyyyy_0_xyyzzzz_0, g_0_xyyyyyy_0_xyyzzzz_1, g_0_xyyyyyy_0_xyzzzz_1, g_0_xyyyyyy_0_xyzzzzz_0, g_0_xyyyyyy_0_xyzzzzz_1, g_0_xyyyyyy_0_yyyyyy_1, g_0_xyyyyyy_0_yyyyyyy_0, g_0_xyyyyyy_0_yyyyyyy_1, g_0_xyyyyyy_0_yyyyyyz_0, g_0_xyyyyyy_0_yyyyyyz_1, g_0_xyyyyyy_0_yyyyyz_1, g_0_xyyyyyy_0_yyyyyzz_0, g_0_xyyyyyy_0_yyyyyzz_1, g_0_xyyyyyy_0_yyyyzz_1, g_0_xyyyyyy_0_yyyyzzz_0, g_0_xyyyyyy_0_yyyyzzz_1, g_0_xyyyyyy_0_yyyzzz_1, g_0_xyyyyyy_0_yyyzzzz_0, g_0_xyyyyyy_0_yyyzzzz_1, g_0_xyyyyyy_0_yyzzzz_1, g_0_xyyyyyy_0_yyzzzzz_0, g_0_xyyyyyy_0_yyzzzzz_1, g_0_xyyyyyy_0_yzzzzz_1, g_0_xyyyyyy_0_yzzzzzz_0, g_0_xyyyyyy_0_yzzzzzz_1, g_0_xyyyyyy_0_zzzzzzz_0, g_0_xyyyyyy_0_zzzzzzz_1, g_0_yyyyyy_0_xxxxxxy_0, g_0_yyyyyy_0_xxxxxxy_1, g_0_yyyyyy_0_xxxxxyy_0, g_0_yyyyyy_0_xxxxxyy_1, g_0_yyyyyy_0_xxxxxyz_0, g_0_yyyyyy_0_xxxxxyz_1, g_0_yyyyyy_0_xxxxyyy_0, g_0_yyyyyy_0_xxxxyyy_1, g_0_yyyyyy_0_xxxxyyz_0, g_0_yyyyyy_0_xxxxyyz_1, g_0_yyyyyy_0_xxxxyzz_0, g_0_yyyyyy_0_xxxxyzz_1, g_0_yyyyyy_0_xxxyyyy_0, g_0_yyyyyy_0_xxxyyyy_1, g_0_yyyyyy_0_xxxyyyz_0, g_0_yyyyyy_0_xxxyyyz_1, g_0_yyyyyy_0_xxxyyzz_0, g_0_yyyyyy_0_xxxyyzz_1, g_0_yyyyyy_0_xxxyzzz_0, g_0_yyyyyy_0_xxxyzzz_1, g_0_yyyyyy_0_xxyyyyy_0, g_0_yyyyyy_0_xxyyyyy_1, g_0_yyyyyy_0_xxyyyyz_0, g_0_yyyyyy_0_xxyyyyz_1, g_0_yyyyyy_0_xxyyyzz_0, g_0_yyyyyy_0_xxyyyzz_1, g_0_yyyyyy_0_xxyyzzz_0, g_0_yyyyyy_0_xxyyzzz_1, g_0_yyyyyy_0_xxyzzzz_0, g_0_yyyyyy_0_xxyzzzz_1, g_0_yyyyyy_0_xyyyyyy_0, g_0_yyyyyy_0_xyyyyyy_1, g_0_yyyyyy_0_xyyyyyz_0, g_0_yyyyyy_0_xyyyyyz_1, g_0_yyyyyy_0_xyyyyzz_0, g_0_yyyyyy_0_xyyyyzz_1, g_0_yyyyyy_0_xyyyzzz_0, g_0_yyyyyy_0_xyyyzzz_1, g_0_yyyyyy_0_xyyzzzz_0, g_0_yyyyyy_0_xyyzzzz_1, g_0_yyyyyy_0_xyzzzzz_0, g_0_yyyyyy_0_xyzzzzz_1, g_0_yyyyyy_0_yyyyyyy_0, g_0_yyyyyy_0_yyyyyyy_1, g_0_yyyyyy_0_yyyyyyz_0, g_0_yyyyyy_0_yyyyyyz_1, g_0_yyyyyy_0_yyyyyzz_0, g_0_yyyyyy_0_yyyyyzz_1, g_0_yyyyyy_0_yyyyzzz_0, g_0_yyyyyy_0_yyyyzzz_1, g_0_yyyyyy_0_yyyzzzz_0, g_0_yyyyyy_0_yyyzzzz_1, g_0_yyyyyy_0_yyzzzzz_0, g_0_yyyyyy_0_yyzzzzz_1, g_0_yyyyyy_0_yzzzzzz_0, g_0_yyyyyy_0_yzzzzzz_1, g_0_yyyyyy_0_zzzzzzz_0, g_0_yyyyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyyy_0_xxxxxxx_0[i] = 5.0 * g_0_xxyyyy_0_xxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxxxxx_0[i] * pb_y + g_0_xxyyyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxxxxy_0[i] = g_0_yyyyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xyyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxxxy_0[i] * pb_x + g_0_xyyyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxxxz_0[i] = 5.0 * g_0_xxyyyy_0_xxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxxxxz_0[i] * pb_y + g_0_xxyyyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxxxyy_0[i] = g_0_yyyyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xyyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxxyy_0[i] * pb_x + g_0_xyyyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxxyz_0[i] = g_0_yyyyyy_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxxyz_0[i] * pb_x + g_0_xyyyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxxzz_0[i] = 5.0 * g_0_xxyyyy_0_xxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxxxzz_0[i] * pb_y + g_0_xxyyyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxxyyy_0[i] = g_0_yyyyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxyyy_0[i] * pb_x + g_0_xyyyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxyyz_0[i] = g_0_yyyyyy_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxyyz_0[i] * pb_x + g_0_xyyyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxyzz_0[i] = g_0_yyyyyy_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxyzz_0[i] * pb_x + g_0_xyyyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxzzz_0[i] = 5.0 * g_0_xxyyyy_0_xxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxxzzz_0[i] * pb_y + g_0_xxyyyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxyyyy_0[i] = g_0_yyyyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyyyy_0[i] * pb_x + g_0_xyyyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxyyyz_0[i] = g_0_yyyyyy_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyyyz_0[i] * pb_x + g_0_xyyyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxyyzz_0[i] = g_0_yyyyyy_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyyzz_0[i] * pb_x + g_0_xyyyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxyzzz_0[i] = g_0_yyyyyy_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyzzz_0[i] * pb_x + g_0_xyyyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxzzzz_0[i] * pb_y + g_0_xxyyyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxyyyyy_0[i] = g_0_yyyyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyyyy_0[i] * pb_x + g_0_xyyyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyyyyz_0[i] = g_0_yyyyyy_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyyyz_0[i] * pb_x + g_0_xyyyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyyyzz_0[i] = g_0_yyyyyy_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyyzz_0[i] * pb_x + g_0_xyyyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyyzzz_0[i] = g_0_yyyyyy_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyzzz_0[i] * pb_x + g_0_xyyyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyzzzz_0[i] = g_0_yyyyyy_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyzzzz_0[i] * pb_x + g_0_xyyyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxzzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxzzzzz_0[i] * pb_y + g_0_xxyyyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xyyyyyy_0[i] = g_0_yyyyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyyyyy_0[i] * pb_x + g_0_xyyyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyyyyz_0[i] = g_0_yyyyyy_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyyyyz_0[i] * pb_x + g_0_xyyyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyyyzz_0[i] = g_0_yyyyyy_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyyyzz_0[i] * pb_x + g_0_xyyyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyyzzz_0[i] = g_0_yyyyyy_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyyzzz_0[i] * pb_x + g_0_xyyyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyzzzz_0[i] = g_0_yyyyyy_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyzzzz_0[i] * pb_x + g_0_xyyyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyzzzzz_0[i] = g_0_yyyyyy_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyzzzzz_0[i] * pb_x + g_0_xyyyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xzzzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xzzzzzz_0[i] * pb_y + g_0_xxyyyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_yyyyyyy_0[i] = g_0_yyyyyy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyyyy_0[i] * pb_x + g_0_xyyyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyyyyz_0[i] = g_0_yyyyyy_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyyyz_0[i] * pb_x + g_0_xyyyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyyyzz_0[i] = g_0_yyyyyy_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyyzz_0[i] * pb_x + g_0_xyyyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyyzzz_0[i] = g_0_yyyyyy_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyzzz_0[i] * pb_x + g_0_xyyyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyzzzz_0[i] = g_0_yyyyyy_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyzzzz_0[i] * pb_x + g_0_xyyyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyzzzzz_0[i] = g_0_yyyyyy_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyzzzzz_0[i] * pb_x + g_0_xyyyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yzzzzzz_0[i] = g_0_yyyyyy_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzzzzzz_0[i] * pb_x + g_0_xyyyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_zzzzzzz_0[i] = g_0_yyyyyy_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_zzzzzzz_0[i] * pb_x + g_0_xyyyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 792-828 components of targeted buffer : SLSK

    auto g_0_xxyyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 792);

    auto g_0_xxyyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 793);

    auto g_0_xxyyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 794);

    auto g_0_xxyyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 795);

    auto g_0_xxyyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 796);

    auto g_0_xxyyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 797);

    auto g_0_xxyyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 798);

    auto g_0_xxyyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 799);

    auto g_0_xxyyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 800);

    auto g_0_xxyyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 801);

    auto g_0_xxyyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 802);

    auto g_0_xxyyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 803);

    auto g_0_xxyyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 804);

    auto g_0_xxyyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 805);

    auto g_0_xxyyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 806);

    auto g_0_xxyyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 807);

    auto g_0_xxyyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 808);

    auto g_0_xxyyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 809);

    auto g_0_xxyyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 810);

    auto g_0_xxyyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 811);

    auto g_0_xxyyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 812);

    auto g_0_xxyyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 813);

    auto g_0_xxyyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 814);

    auto g_0_xxyyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 815);

    auto g_0_xxyyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 816);

    auto g_0_xxyyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 817);

    auto g_0_xxyyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 818);

    auto g_0_xxyyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 819);

    auto g_0_xxyyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 820);

    auto g_0_xxyyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 821);

    auto g_0_xxyyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 822);

    auto g_0_xxyyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 823);

    auto g_0_xxyyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 824);

    auto g_0_xxyyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 825);

    auto g_0_xxyyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 826);

    auto g_0_xxyyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 827);

    #pragma omp simd aligned(g_0_xxyyyyy_0_xxxxxx_1, g_0_xxyyyyy_0_xxxxxxx_0, g_0_xxyyyyy_0_xxxxxxx_1, g_0_xxyyyyy_0_xxxxxxy_0, g_0_xxyyyyy_0_xxxxxxy_1, g_0_xxyyyyy_0_xxxxxxz_0, g_0_xxyyyyy_0_xxxxxxz_1, g_0_xxyyyyy_0_xxxxxy_1, g_0_xxyyyyy_0_xxxxxyy_0, g_0_xxyyyyy_0_xxxxxyy_1, g_0_xxyyyyy_0_xxxxxyz_0, g_0_xxyyyyy_0_xxxxxyz_1, g_0_xxyyyyy_0_xxxxxz_1, g_0_xxyyyyy_0_xxxxxzz_0, g_0_xxyyyyy_0_xxxxxzz_1, g_0_xxyyyyy_0_xxxxyy_1, g_0_xxyyyyy_0_xxxxyyy_0, g_0_xxyyyyy_0_xxxxyyy_1, g_0_xxyyyyy_0_xxxxyyz_0, g_0_xxyyyyy_0_xxxxyyz_1, g_0_xxyyyyy_0_xxxxyz_1, g_0_xxyyyyy_0_xxxxyzz_0, g_0_xxyyyyy_0_xxxxyzz_1, g_0_xxyyyyy_0_xxxxzz_1, g_0_xxyyyyy_0_xxxxzzz_0, g_0_xxyyyyy_0_xxxxzzz_1, g_0_xxyyyyy_0_xxxyyy_1, g_0_xxyyyyy_0_xxxyyyy_0, g_0_xxyyyyy_0_xxxyyyy_1, g_0_xxyyyyy_0_xxxyyyz_0, g_0_xxyyyyy_0_xxxyyyz_1, g_0_xxyyyyy_0_xxxyyz_1, g_0_xxyyyyy_0_xxxyyzz_0, g_0_xxyyyyy_0_xxxyyzz_1, g_0_xxyyyyy_0_xxxyzz_1, g_0_xxyyyyy_0_xxxyzzz_0, g_0_xxyyyyy_0_xxxyzzz_1, g_0_xxyyyyy_0_xxxzzz_1, g_0_xxyyyyy_0_xxxzzzz_0, g_0_xxyyyyy_0_xxxzzzz_1, g_0_xxyyyyy_0_xxyyyy_1, g_0_xxyyyyy_0_xxyyyyy_0, g_0_xxyyyyy_0_xxyyyyy_1, g_0_xxyyyyy_0_xxyyyyz_0, g_0_xxyyyyy_0_xxyyyyz_1, g_0_xxyyyyy_0_xxyyyz_1, g_0_xxyyyyy_0_xxyyyzz_0, g_0_xxyyyyy_0_xxyyyzz_1, g_0_xxyyyyy_0_xxyyzz_1, g_0_xxyyyyy_0_xxyyzzz_0, g_0_xxyyyyy_0_xxyyzzz_1, g_0_xxyyyyy_0_xxyzzz_1, g_0_xxyyyyy_0_xxyzzzz_0, g_0_xxyyyyy_0_xxyzzzz_1, g_0_xxyyyyy_0_xxzzzz_1, g_0_xxyyyyy_0_xxzzzzz_0, g_0_xxyyyyy_0_xxzzzzz_1, g_0_xxyyyyy_0_xyyyyy_1, g_0_xxyyyyy_0_xyyyyyy_0, g_0_xxyyyyy_0_xyyyyyy_1, g_0_xxyyyyy_0_xyyyyyz_0, g_0_xxyyyyy_0_xyyyyyz_1, g_0_xxyyyyy_0_xyyyyz_1, g_0_xxyyyyy_0_xyyyyzz_0, g_0_xxyyyyy_0_xyyyyzz_1, g_0_xxyyyyy_0_xyyyzz_1, g_0_xxyyyyy_0_xyyyzzz_0, g_0_xxyyyyy_0_xyyyzzz_1, g_0_xxyyyyy_0_xyyzzz_1, g_0_xxyyyyy_0_xyyzzzz_0, g_0_xxyyyyy_0_xyyzzzz_1, g_0_xxyyyyy_0_xyzzzz_1, g_0_xxyyyyy_0_xyzzzzz_0, g_0_xxyyyyy_0_xyzzzzz_1, g_0_xxyyyyy_0_xzzzzz_1, g_0_xxyyyyy_0_xzzzzzz_0, g_0_xxyyyyy_0_xzzzzzz_1, g_0_xxyyyyy_0_yyyyyy_1, g_0_xxyyyyy_0_yyyyyyy_0, g_0_xxyyyyy_0_yyyyyyy_1, g_0_xxyyyyy_0_yyyyyyz_0, g_0_xxyyyyy_0_yyyyyyz_1, g_0_xxyyyyy_0_yyyyyz_1, g_0_xxyyyyy_0_yyyyyzz_0, g_0_xxyyyyy_0_yyyyyzz_1, g_0_xxyyyyy_0_yyyyzz_1, g_0_xxyyyyy_0_yyyyzzz_0, g_0_xxyyyyy_0_yyyyzzz_1, g_0_xxyyyyy_0_yyyzzz_1, g_0_xxyyyyy_0_yyyzzzz_0, g_0_xxyyyyy_0_yyyzzzz_1, g_0_xxyyyyy_0_yyzzzz_1, g_0_xxyyyyy_0_yyzzzzz_0, g_0_xxyyyyy_0_yyzzzzz_1, g_0_xxyyyyy_0_yzzzzz_1, g_0_xxyyyyy_0_yzzzzzz_0, g_0_xxyyyyy_0_yzzzzzz_1, g_0_xxyyyyy_0_zzzzzz_1, g_0_xxyyyyy_0_zzzzzzz_0, g_0_xxyyyyy_0_zzzzzzz_1, g_0_xxyyyyyz_0_xxxxxxx_0, g_0_xxyyyyyz_0_xxxxxxy_0, g_0_xxyyyyyz_0_xxxxxxz_0, g_0_xxyyyyyz_0_xxxxxyy_0, g_0_xxyyyyyz_0_xxxxxyz_0, g_0_xxyyyyyz_0_xxxxxzz_0, g_0_xxyyyyyz_0_xxxxyyy_0, g_0_xxyyyyyz_0_xxxxyyz_0, g_0_xxyyyyyz_0_xxxxyzz_0, g_0_xxyyyyyz_0_xxxxzzz_0, g_0_xxyyyyyz_0_xxxyyyy_0, g_0_xxyyyyyz_0_xxxyyyz_0, g_0_xxyyyyyz_0_xxxyyzz_0, g_0_xxyyyyyz_0_xxxyzzz_0, g_0_xxyyyyyz_0_xxxzzzz_0, g_0_xxyyyyyz_0_xxyyyyy_0, g_0_xxyyyyyz_0_xxyyyyz_0, g_0_xxyyyyyz_0_xxyyyzz_0, g_0_xxyyyyyz_0_xxyyzzz_0, g_0_xxyyyyyz_0_xxyzzzz_0, g_0_xxyyyyyz_0_xxzzzzz_0, g_0_xxyyyyyz_0_xyyyyyy_0, g_0_xxyyyyyz_0_xyyyyyz_0, g_0_xxyyyyyz_0_xyyyyzz_0, g_0_xxyyyyyz_0_xyyyzzz_0, g_0_xxyyyyyz_0_xyyzzzz_0, g_0_xxyyyyyz_0_xyzzzzz_0, g_0_xxyyyyyz_0_xzzzzzz_0, g_0_xxyyyyyz_0_yyyyyyy_0, g_0_xxyyyyyz_0_yyyyyyz_0, g_0_xxyyyyyz_0_yyyyyzz_0, g_0_xxyyyyyz_0_yyyyzzz_0, g_0_xxyyyyyz_0_yyyzzzz_0, g_0_xxyyyyyz_0_yyzzzzz_0, g_0_xxyyyyyz_0_yzzzzzz_0, g_0_xxyyyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyyz_0_xxxxxxx_0[i] = g_0_xxyyyyy_0_xxxxxxx_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxxxy_0[i] = g_0_xxyyyyy_0_xxxxxxy_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxxxz_0[i] = g_0_xxyyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxxxz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxxyy_0[i] = g_0_xxyyyyy_0_xxxxxyy_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxxyz_0[i] = g_0_xxyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxxyz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxyyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxxzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxyyy_0[i] = g_0_xxyyyyy_0_xxxxyyy_0[i] * pb_z + g_0_xxyyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxyyz_0[i] = g_0_xxyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxyyz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxyzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxyyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyyyy_0[i] = g_0_xxyyyyy_0_xxxyyyy_0[i] * pb_z + g_0_xxyyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyyyz_0[i] = g_0_xxyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyyyz_0[i] * pb_z + g_0_xxyyyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyyzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxyyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyyyy_0[i] = g_0_xxyyyyy_0_xxyyyyy_0[i] * pb_z + g_0_xxyyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyyyz_0[i] = g_0_xxyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyyyz_0[i] * pb_z + g_0_xxyyyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyyzz_0[i] * pb_z + g_0_xxyyyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxyyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyyyy_0[i] = g_0_xxyyyyy_0_xyyyyyy_0[i] * pb_z + g_0_xxyyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyyyz_0[i] = g_0_xxyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyyyz_0[i] * pb_z + g_0_xxyyyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyyzz_0[i] * pb_z + g_0_xxyyyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyzzz_0[i] * pb_z + g_0_xxyyyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxyyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xzzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyyyy_0[i] = g_0_xxyyyyy_0_yyyyyyy_0[i] * pb_z + g_0_xxyyyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyyyz_0[i] = g_0_xxyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyyyyz_0[i] * pb_z + g_0_xxyyyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyyyzz_0[i] * pb_z + g_0_xxyyyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyyzzz_0[i] * pb_z + g_0_xxyyyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyzzzz_0[i] * pb_z + g_0_xxyyyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yzzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxyyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_zzzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 828-864 components of targeted buffer : SLSK

    auto g_0_xxyyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 828);

    auto g_0_xxyyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 829);

    auto g_0_xxyyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 830);

    auto g_0_xxyyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 831);

    auto g_0_xxyyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 832);

    auto g_0_xxyyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 833);

    auto g_0_xxyyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 834);

    auto g_0_xxyyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 835);

    auto g_0_xxyyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 836);

    auto g_0_xxyyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 837);

    auto g_0_xxyyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 838);

    auto g_0_xxyyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 839);

    auto g_0_xxyyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 840);

    auto g_0_xxyyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 841);

    auto g_0_xxyyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 842);

    auto g_0_xxyyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 843);

    auto g_0_xxyyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 844);

    auto g_0_xxyyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 845);

    auto g_0_xxyyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 846);

    auto g_0_xxyyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 847);

    auto g_0_xxyyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 848);

    auto g_0_xxyyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 849);

    auto g_0_xxyyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 850);

    auto g_0_xxyyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 851);

    auto g_0_xxyyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 852);

    auto g_0_xxyyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 853);

    auto g_0_xxyyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 854);

    auto g_0_xxyyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 855);

    auto g_0_xxyyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 856);

    auto g_0_xxyyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 857);

    auto g_0_xxyyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 858);

    auto g_0_xxyyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 859);

    auto g_0_xxyyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 860);

    auto g_0_xxyyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 861);

    auto g_0_xxyyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 862);

    auto g_0_xxyyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 863);

    #pragma omp simd aligned(g_0_xxyyyy_0_xxxxxxy_0, g_0_xxyyyy_0_xxxxxxy_1, g_0_xxyyyy_0_xxxxxyy_0, g_0_xxyyyy_0_xxxxxyy_1, g_0_xxyyyy_0_xxxxyyy_0, g_0_xxyyyy_0_xxxxyyy_1, g_0_xxyyyy_0_xxxyyyy_0, g_0_xxyyyy_0_xxxyyyy_1, g_0_xxyyyy_0_xxyyyyy_0, g_0_xxyyyy_0_xxyyyyy_1, g_0_xxyyyy_0_xyyyyyy_0, g_0_xxyyyy_0_xyyyyyy_1, g_0_xxyyyyz_0_xxxxxxy_0, g_0_xxyyyyz_0_xxxxxxy_1, g_0_xxyyyyz_0_xxxxxyy_0, g_0_xxyyyyz_0_xxxxxyy_1, g_0_xxyyyyz_0_xxxxyyy_0, g_0_xxyyyyz_0_xxxxyyy_1, g_0_xxyyyyz_0_xxxyyyy_0, g_0_xxyyyyz_0_xxxyyyy_1, g_0_xxyyyyz_0_xxyyyyy_0, g_0_xxyyyyz_0_xxyyyyy_1, g_0_xxyyyyz_0_xyyyyyy_0, g_0_xxyyyyz_0_xyyyyyy_1, g_0_xxyyyyzz_0_xxxxxxx_0, g_0_xxyyyyzz_0_xxxxxxy_0, g_0_xxyyyyzz_0_xxxxxxz_0, g_0_xxyyyyzz_0_xxxxxyy_0, g_0_xxyyyyzz_0_xxxxxyz_0, g_0_xxyyyyzz_0_xxxxxzz_0, g_0_xxyyyyzz_0_xxxxyyy_0, g_0_xxyyyyzz_0_xxxxyyz_0, g_0_xxyyyyzz_0_xxxxyzz_0, g_0_xxyyyyzz_0_xxxxzzz_0, g_0_xxyyyyzz_0_xxxyyyy_0, g_0_xxyyyyzz_0_xxxyyyz_0, g_0_xxyyyyzz_0_xxxyyzz_0, g_0_xxyyyyzz_0_xxxyzzz_0, g_0_xxyyyyzz_0_xxxzzzz_0, g_0_xxyyyyzz_0_xxyyyyy_0, g_0_xxyyyyzz_0_xxyyyyz_0, g_0_xxyyyyzz_0_xxyyyzz_0, g_0_xxyyyyzz_0_xxyyzzz_0, g_0_xxyyyyzz_0_xxyzzzz_0, g_0_xxyyyyzz_0_xxzzzzz_0, g_0_xxyyyyzz_0_xyyyyyy_0, g_0_xxyyyyzz_0_xyyyyyz_0, g_0_xxyyyyzz_0_xyyyyzz_0, g_0_xxyyyyzz_0_xyyyzzz_0, g_0_xxyyyyzz_0_xyyzzzz_0, g_0_xxyyyyzz_0_xyzzzzz_0, g_0_xxyyyyzz_0_xzzzzzz_0, g_0_xxyyyyzz_0_yyyyyyy_0, g_0_xxyyyyzz_0_yyyyyyz_0, g_0_xxyyyyzz_0_yyyyyzz_0, g_0_xxyyyyzz_0_yyyyzzz_0, g_0_xxyyyyzz_0_yyyzzzz_0, g_0_xxyyyyzz_0_yyzzzzz_0, g_0_xxyyyyzz_0_yzzzzzz_0, g_0_xxyyyyzz_0_zzzzzzz_0, g_0_xxyyyzz_0_xxxxxxx_0, g_0_xxyyyzz_0_xxxxxxx_1, g_0_xxyyyzz_0_xxxxxxz_0, g_0_xxyyyzz_0_xxxxxxz_1, g_0_xxyyyzz_0_xxxxxzz_0, g_0_xxyyyzz_0_xxxxxzz_1, g_0_xxyyyzz_0_xxxxzzz_0, g_0_xxyyyzz_0_xxxxzzz_1, g_0_xxyyyzz_0_xxxzzzz_0, g_0_xxyyyzz_0_xxxzzzz_1, g_0_xxyyyzz_0_xxzzzzz_0, g_0_xxyyyzz_0_xxzzzzz_1, g_0_xxyyyzz_0_xzzzzzz_0, g_0_xxyyyzz_0_xzzzzzz_1, g_0_xxyyzz_0_xxxxxxx_0, g_0_xxyyzz_0_xxxxxxx_1, g_0_xxyyzz_0_xxxxxxz_0, g_0_xxyyzz_0_xxxxxxz_1, g_0_xxyyzz_0_xxxxxzz_0, g_0_xxyyzz_0_xxxxxzz_1, g_0_xxyyzz_0_xxxxzzz_0, g_0_xxyyzz_0_xxxxzzz_1, g_0_xxyyzz_0_xxxzzzz_0, g_0_xxyyzz_0_xxxzzzz_1, g_0_xxyyzz_0_xxzzzzz_0, g_0_xxyyzz_0_xxzzzzz_1, g_0_xxyyzz_0_xzzzzzz_0, g_0_xxyyzz_0_xzzzzzz_1, g_0_xyyyyzz_0_xxxxxyz_0, g_0_xyyyyzz_0_xxxxxyz_1, g_0_xyyyyzz_0_xxxxyyz_0, g_0_xyyyyzz_0_xxxxyyz_1, g_0_xyyyyzz_0_xxxxyz_1, g_0_xyyyyzz_0_xxxxyzz_0, g_0_xyyyyzz_0_xxxxyzz_1, g_0_xyyyyzz_0_xxxyyyz_0, g_0_xyyyyzz_0_xxxyyyz_1, g_0_xyyyyzz_0_xxxyyz_1, g_0_xyyyyzz_0_xxxyyzz_0, g_0_xyyyyzz_0_xxxyyzz_1, g_0_xyyyyzz_0_xxxyzz_1, g_0_xyyyyzz_0_xxxyzzz_0, g_0_xyyyyzz_0_xxxyzzz_1, g_0_xyyyyzz_0_xxyyyyz_0, g_0_xyyyyzz_0_xxyyyyz_1, g_0_xyyyyzz_0_xxyyyz_1, g_0_xyyyyzz_0_xxyyyzz_0, g_0_xyyyyzz_0_xxyyyzz_1, g_0_xyyyyzz_0_xxyyzz_1, g_0_xyyyyzz_0_xxyyzzz_0, g_0_xyyyyzz_0_xxyyzzz_1, g_0_xyyyyzz_0_xxyzzz_1, g_0_xyyyyzz_0_xxyzzzz_0, g_0_xyyyyzz_0_xxyzzzz_1, g_0_xyyyyzz_0_xyyyyyz_0, g_0_xyyyyzz_0_xyyyyyz_1, g_0_xyyyyzz_0_xyyyyz_1, g_0_xyyyyzz_0_xyyyyzz_0, g_0_xyyyyzz_0_xyyyyzz_1, g_0_xyyyyzz_0_xyyyzz_1, g_0_xyyyyzz_0_xyyyzzz_0, g_0_xyyyyzz_0_xyyyzzz_1, g_0_xyyyyzz_0_xyyzzz_1, g_0_xyyyyzz_0_xyyzzzz_0, g_0_xyyyyzz_0_xyyzzzz_1, g_0_xyyyyzz_0_xyzzzz_1, g_0_xyyyyzz_0_xyzzzzz_0, g_0_xyyyyzz_0_xyzzzzz_1, g_0_xyyyyzz_0_yyyyyyy_0, g_0_xyyyyzz_0_yyyyyyy_1, g_0_xyyyyzz_0_yyyyyyz_0, g_0_xyyyyzz_0_yyyyyyz_1, g_0_xyyyyzz_0_yyyyyz_1, g_0_xyyyyzz_0_yyyyyzz_0, g_0_xyyyyzz_0_yyyyyzz_1, g_0_xyyyyzz_0_yyyyzz_1, g_0_xyyyyzz_0_yyyyzzz_0, g_0_xyyyyzz_0_yyyyzzz_1, g_0_xyyyyzz_0_yyyzzz_1, g_0_xyyyyzz_0_yyyzzzz_0, g_0_xyyyyzz_0_yyyzzzz_1, g_0_xyyyyzz_0_yyzzzz_1, g_0_xyyyyzz_0_yyzzzzz_0, g_0_xyyyyzz_0_yyzzzzz_1, g_0_xyyyyzz_0_yzzzzz_1, g_0_xyyyyzz_0_yzzzzzz_0, g_0_xyyyyzz_0_yzzzzzz_1, g_0_xyyyyzz_0_zzzzzzz_0, g_0_xyyyyzz_0_zzzzzzz_1, g_0_yyyyzz_0_xxxxxyz_0, g_0_yyyyzz_0_xxxxxyz_1, g_0_yyyyzz_0_xxxxyyz_0, g_0_yyyyzz_0_xxxxyyz_1, g_0_yyyyzz_0_xxxxyzz_0, g_0_yyyyzz_0_xxxxyzz_1, g_0_yyyyzz_0_xxxyyyz_0, g_0_yyyyzz_0_xxxyyyz_1, g_0_yyyyzz_0_xxxyyzz_0, g_0_yyyyzz_0_xxxyyzz_1, g_0_yyyyzz_0_xxxyzzz_0, g_0_yyyyzz_0_xxxyzzz_1, g_0_yyyyzz_0_xxyyyyz_0, g_0_yyyyzz_0_xxyyyyz_1, g_0_yyyyzz_0_xxyyyzz_0, g_0_yyyyzz_0_xxyyyzz_1, g_0_yyyyzz_0_xxyyzzz_0, g_0_yyyyzz_0_xxyyzzz_1, g_0_yyyyzz_0_xxyzzzz_0, g_0_yyyyzz_0_xxyzzzz_1, g_0_yyyyzz_0_xyyyyyz_0, g_0_yyyyzz_0_xyyyyyz_1, g_0_yyyyzz_0_xyyyyzz_0, g_0_yyyyzz_0_xyyyyzz_1, g_0_yyyyzz_0_xyyyzzz_0, g_0_yyyyzz_0_xyyyzzz_1, g_0_yyyyzz_0_xyyzzzz_0, g_0_yyyyzz_0_xyyzzzz_1, g_0_yyyyzz_0_xyzzzzz_0, g_0_yyyyzz_0_xyzzzzz_1, g_0_yyyyzz_0_yyyyyyy_0, g_0_yyyyzz_0_yyyyyyy_1, g_0_yyyyzz_0_yyyyyyz_0, g_0_yyyyzz_0_yyyyyyz_1, g_0_yyyyzz_0_yyyyyzz_0, g_0_yyyyzz_0_yyyyyzz_1, g_0_yyyyzz_0_yyyyzzz_0, g_0_yyyyzz_0_yyyyzzz_1, g_0_yyyyzz_0_yyyzzzz_0, g_0_yyyyzz_0_yyyzzzz_1, g_0_yyyyzz_0_yyzzzzz_0, g_0_yyyyzz_0_yyzzzzz_1, g_0_yyyyzz_0_yzzzzzz_0, g_0_yyyyzz_0_yzzzzzz_1, g_0_yyyyzz_0_zzzzzzz_0, g_0_yyyyzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyzz_0_xxxxxxx_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxxxx_0[i] * pb_y + g_0_xxyyyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxxxxy_0[i] = g_0_xxyyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxxxxy_0[i] * pb_z + g_0_xxyyyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxxxxz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxxxz_0[i] * pb_y + g_0_xxyyyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxxxyy_0[i] = g_0_xxyyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxxxyy_0[i] * pb_z + g_0_xxyyyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxxxyz_0[i] = g_0_yyyyzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxxxyz_0[i] * pb_x + g_0_xyyyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxxxzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxxzz_0[i] * pb_y + g_0_xxyyyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxxyyy_0[i] = g_0_xxyyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxxyyy_0[i] * pb_z + g_0_xxyyyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxxyyz_0[i] = g_0_yyyyzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxxyyz_0[i] * pb_x + g_0_xyyyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxxyzz_0[i] = g_0_yyyyzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxxyzz_0[i] * pb_x + g_0_xyyyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxxzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxzzz_0[i] * pb_y + g_0_xxyyyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxyyyy_0[i] = g_0_xxyyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxyyyy_0[i] * pb_z + g_0_xxyyyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxyyyz_0[i] = g_0_yyyyzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxyyyz_0[i] * pb_x + g_0_xyyyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxyyzz_0[i] = g_0_yyyyzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxyyzz_0[i] * pb_x + g_0_xyyyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxyzzz_0[i] = g_0_yyyyzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxyzzz_0[i] * pb_x + g_0_xyyyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxzzzz_0[i] * pb_y + g_0_xxyyyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxyyyyy_0[i] = g_0_xxyyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxyyyyy_0[i] * pb_z + g_0_xxyyyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxyyyyz_0[i] = g_0_yyyyzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyyyyz_0[i] * pb_x + g_0_xyyyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxyyyzz_0[i] = g_0_yyyyzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyyyzz_0[i] * pb_x + g_0_xyyyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxyyzzz_0[i] = g_0_yyyyzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyyzzz_0[i] * pb_x + g_0_xyyyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxyzzzz_0[i] = g_0_yyyyzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyzzzz_0[i] * pb_x + g_0_xyyyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxzzzzz_0[i] * pb_y + g_0_xxyyyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xyyyyyy_0[i] = g_0_xxyyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xyyyyyy_0[i] * pb_z + g_0_xxyyyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xyyyyyz_0[i] = g_0_yyyyzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyyyyyz_0[i] * pb_x + g_0_xyyyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyyyyzz_0[i] = g_0_yyyyzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyyyyzz_0[i] * pb_x + g_0_xyyyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyyyzzz_0[i] = g_0_yyyyzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyyyzzz_0[i] * pb_x + g_0_xyyyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyyzzzz_0[i] = g_0_yyyyzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyyzzzz_0[i] * pb_x + g_0_xyyyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyzzzzz_0[i] = g_0_yyyyzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyzzzzz_0[i] * pb_x + g_0_xyyyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xzzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xzzzzzz_0[i] * pb_y + g_0_xxyyyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_yyyyyyy_0[i] = g_0_yyyyzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyyyy_0[i] * pb_x + g_0_xyyyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyyyyz_0[i] = g_0_yyyyzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyyyz_0[i] * pb_x + g_0_xyyyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyyyzz_0[i] = g_0_yyyyzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyyzz_0[i] * pb_x + g_0_xyyyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyyzzz_0[i] = g_0_yyyyzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyzzz_0[i] * pb_x + g_0_xyyyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyzzzz_0[i] = g_0_yyyyzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyzzzz_0[i] * pb_x + g_0_xyyyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyzzzzz_0[i] = g_0_yyyyzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyzzzzz_0[i] * pb_x + g_0_xyyyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yzzzzzz_0[i] = g_0_yyyyzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzzzzzz_0[i] * pb_x + g_0_xyyyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_zzzzzzz_0[i] = g_0_yyyyzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_zzzzzzz_0[i] * pb_x + g_0_xyyyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 864-900 components of targeted buffer : SLSK

    auto g_0_xxyyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 864);

    auto g_0_xxyyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 865);

    auto g_0_xxyyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 866);

    auto g_0_xxyyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 867);

    auto g_0_xxyyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 868);

    auto g_0_xxyyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 869);

    auto g_0_xxyyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 870);

    auto g_0_xxyyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 871);

    auto g_0_xxyyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 872);

    auto g_0_xxyyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 873);

    auto g_0_xxyyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 874);

    auto g_0_xxyyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 875);

    auto g_0_xxyyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 876);

    auto g_0_xxyyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 877);

    auto g_0_xxyyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 878);

    auto g_0_xxyyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 879);

    auto g_0_xxyyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 880);

    auto g_0_xxyyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 881);

    auto g_0_xxyyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 882);

    auto g_0_xxyyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 883);

    auto g_0_xxyyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 884);

    auto g_0_xxyyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 885);

    auto g_0_xxyyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 886);

    auto g_0_xxyyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 887);

    auto g_0_xxyyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 888);

    auto g_0_xxyyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 889);

    auto g_0_xxyyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 890);

    auto g_0_xxyyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 891);

    auto g_0_xxyyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 892);

    auto g_0_xxyyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 893);

    auto g_0_xxyyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 894);

    auto g_0_xxyyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 895);

    auto g_0_xxyyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 896);

    auto g_0_xxyyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 897);

    auto g_0_xxyyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 898);

    auto g_0_xxyyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 899);

    #pragma omp simd aligned(g_0_xxyyyz_0_xxxxxxy_0, g_0_xxyyyz_0_xxxxxxy_1, g_0_xxyyyz_0_xxxxxyy_0, g_0_xxyyyz_0_xxxxxyy_1, g_0_xxyyyz_0_xxxxyyy_0, g_0_xxyyyz_0_xxxxyyy_1, g_0_xxyyyz_0_xxxyyyy_0, g_0_xxyyyz_0_xxxyyyy_1, g_0_xxyyyz_0_xxyyyyy_0, g_0_xxyyyz_0_xxyyyyy_1, g_0_xxyyyz_0_xyyyyyy_0, g_0_xxyyyz_0_xyyyyyy_1, g_0_xxyyyzz_0_xxxxxxy_0, g_0_xxyyyzz_0_xxxxxxy_1, g_0_xxyyyzz_0_xxxxxyy_0, g_0_xxyyyzz_0_xxxxxyy_1, g_0_xxyyyzz_0_xxxxyyy_0, g_0_xxyyyzz_0_xxxxyyy_1, g_0_xxyyyzz_0_xxxyyyy_0, g_0_xxyyyzz_0_xxxyyyy_1, g_0_xxyyyzz_0_xxyyyyy_0, g_0_xxyyyzz_0_xxyyyyy_1, g_0_xxyyyzz_0_xyyyyyy_0, g_0_xxyyyzz_0_xyyyyyy_1, g_0_xxyyyzzz_0_xxxxxxx_0, g_0_xxyyyzzz_0_xxxxxxy_0, g_0_xxyyyzzz_0_xxxxxxz_0, g_0_xxyyyzzz_0_xxxxxyy_0, g_0_xxyyyzzz_0_xxxxxyz_0, g_0_xxyyyzzz_0_xxxxxzz_0, g_0_xxyyyzzz_0_xxxxyyy_0, g_0_xxyyyzzz_0_xxxxyyz_0, g_0_xxyyyzzz_0_xxxxyzz_0, g_0_xxyyyzzz_0_xxxxzzz_0, g_0_xxyyyzzz_0_xxxyyyy_0, g_0_xxyyyzzz_0_xxxyyyz_0, g_0_xxyyyzzz_0_xxxyyzz_0, g_0_xxyyyzzz_0_xxxyzzz_0, g_0_xxyyyzzz_0_xxxzzzz_0, g_0_xxyyyzzz_0_xxyyyyy_0, g_0_xxyyyzzz_0_xxyyyyz_0, g_0_xxyyyzzz_0_xxyyyzz_0, g_0_xxyyyzzz_0_xxyyzzz_0, g_0_xxyyyzzz_0_xxyzzzz_0, g_0_xxyyyzzz_0_xxzzzzz_0, g_0_xxyyyzzz_0_xyyyyyy_0, g_0_xxyyyzzz_0_xyyyyyz_0, g_0_xxyyyzzz_0_xyyyyzz_0, g_0_xxyyyzzz_0_xyyyzzz_0, g_0_xxyyyzzz_0_xyyzzzz_0, g_0_xxyyyzzz_0_xyzzzzz_0, g_0_xxyyyzzz_0_xzzzzzz_0, g_0_xxyyyzzz_0_yyyyyyy_0, g_0_xxyyyzzz_0_yyyyyyz_0, g_0_xxyyyzzz_0_yyyyyzz_0, g_0_xxyyyzzz_0_yyyyzzz_0, g_0_xxyyyzzz_0_yyyzzzz_0, g_0_xxyyyzzz_0_yyzzzzz_0, g_0_xxyyyzzz_0_yzzzzzz_0, g_0_xxyyyzzz_0_zzzzzzz_0, g_0_xxyyzzz_0_xxxxxxx_0, g_0_xxyyzzz_0_xxxxxxx_1, g_0_xxyyzzz_0_xxxxxxz_0, g_0_xxyyzzz_0_xxxxxxz_1, g_0_xxyyzzz_0_xxxxxzz_0, g_0_xxyyzzz_0_xxxxxzz_1, g_0_xxyyzzz_0_xxxxzzz_0, g_0_xxyyzzz_0_xxxxzzz_1, g_0_xxyyzzz_0_xxxzzzz_0, g_0_xxyyzzz_0_xxxzzzz_1, g_0_xxyyzzz_0_xxzzzzz_0, g_0_xxyyzzz_0_xxzzzzz_1, g_0_xxyyzzz_0_xzzzzzz_0, g_0_xxyyzzz_0_xzzzzzz_1, g_0_xxyzzz_0_xxxxxxx_0, g_0_xxyzzz_0_xxxxxxx_1, g_0_xxyzzz_0_xxxxxxz_0, g_0_xxyzzz_0_xxxxxxz_1, g_0_xxyzzz_0_xxxxxzz_0, g_0_xxyzzz_0_xxxxxzz_1, g_0_xxyzzz_0_xxxxzzz_0, g_0_xxyzzz_0_xxxxzzz_1, g_0_xxyzzz_0_xxxzzzz_0, g_0_xxyzzz_0_xxxzzzz_1, g_0_xxyzzz_0_xxzzzzz_0, g_0_xxyzzz_0_xxzzzzz_1, g_0_xxyzzz_0_xzzzzzz_0, g_0_xxyzzz_0_xzzzzzz_1, g_0_xyyyzzz_0_xxxxxyz_0, g_0_xyyyzzz_0_xxxxxyz_1, g_0_xyyyzzz_0_xxxxyyz_0, g_0_xyyyzzz_0_xxxxyyz_1, g_0_xyyyzzz_0_xxxxyz_1, g_0_xyyyzzz_0_xxxxyzz_0, g_0_xyyyzzz_0_xxxxyzz_1, g_0_xyyyzzz_0_xxxyyyz_0, g_0_xyyyzzz_0_xxxyyyz_1, g_0_xyyyzzz_0_xxxyyz_1, g_0_xyyyzzz_0_xxxyyzz_0, g_0_xyyyzzz_0_xxxyyzz_1, g_0_xyyyzzz_0_xxxyzz_1, g_0_xyyyzzz_0_xxxyzzz_0, g_0_xyyyzzz_0_xxxyzzz_1, g_0_xyyyzzz_0_xxyyyyz_0, g_0_xyyyzzz_0_xxyyyyz_1, g_0_xyyyzzz_0_xxyyyz_1, g_0_xyyyzzz_0_xxyyyzz_0, g_0_xyyyzzz_0_xxyyyzz_1, g_0_xyyyzzz_0_xxyyzz_1, g_0_xyyyzzz_0_xxyyzzz_0, g_0_xyyyzzz_0_xxyyzzz_1, g_0_xyyyzzz_0_xxyzzz_1, g_0_xyyyzzz_0_xxyzzzz_0, g_0_xyyyzzz_0_xxyzzzz_1, g_0_xyyyzzz_0_xyyyyyz_0, g_0_xyyyzzz_0_xyyyyyz_1, g_0_xyyyzzz_0_xyyyyz_1, g_0_xyyyzzz_0_xyyyyzz_0, g_0_xyyyzzz_0_xyyyyzz_1, g_0_xyyyzzz_0_xyyyzz_1, g_0_xyyyzzz_0_xyyyzzz_0, g_0_xyyyzzz_0_xyyyzzz_1, g_0_xyyyzzz_0_xyyzzz_1, g_0_xyyyzzz_0_xyyzzzz_0, g_0_xyyyzzz_0_xyyzzzz_1, g_0_xyyyzzz_0_xyzzzz_1, g_0_xyyyzzz_0_xyzzzzz_0, g_0_xyyyzzz_0_xyzzzzz_1, g_0_xyyyzzz_0_yyyyyyy_0, g_0_xyyyzzz_0_yyyyyyy_1, g_0_xyyyzzz_0_yyyyyyz_0, g_0_xyyyzzz_0_yyyyyyz_1, g_0_xyyyzzz_0_yyyyyz_1, g_0_xyyyzzz_0_yyyyyzz_0, g_0_xyyyzzz_0_yyyyyzz_1, g_0_xyyyzzz_0_yyyyzz_1, g_0_xyyyzzz_0_yyyyzzz_0, g_0_xyyyzzz_0_yyyyzzz_1, g_0_xyyyzzz_0_yyyzzz_1, g_0_xyyyzzz_0_yyyzzzz_0, g_0_xyyyzzz_0_yyyzzzz_1, g_0_xyyyzzz_0_yyzzzz_1, g_0_xyyyzzz_0_yyzzzzz_0, g_0_xyyyzzz_0_yyzzzzz_1, g_0_xyyyzzz_0_yzzzzz_1, g_0_xyyyzzz_0_yzzzzzz_0, g_0_xyyyzzz_0_yzzzzzz_1, g_0_xyyyzzz_0_zzzzzzz_0, g_0_xyyyzzz_0_zzzzzzz_1, g_0_yyyzzz_0_xxxxxyz_0, g_0_yyyzzz_0_xxxxxyz_1, g_0_yyyzzz_0_xxxxyyz_0, g_0_yyyzzz_0_xxxxyyz_1, g_0_yyyzzz_0_xxxxyzz_0, g_0_yyyzzz_0_xxxxyzz_1, g_0_yyyzzz_0_xxxyyyz_0, g_0_yyyzzz_0_xxxyyyz_1, g_0_yyyzzz_0_xxxyyzz_0, g_0_yyyzzz_0_xxxyyzz_1, g_0_yyyzzz_0_xxxyzzz_0, g_0_yyyzzz_0_xxxyzzz_1, g_0_yyyzzz_0_xxyyyyz_0, g_0_yyyzzz_0_xxyyyyz_1, g_0_yyyzzz_0_xxyyyzz_0, g_0_yyyzzz_0_xxyyyzz_1, g_0_yyyzzz_0_xxyyzzz_0, g_0_yyyzzz_0_xxyyzzz_1, g_0_yyyzzz_0_xxyzzzz_0, g_0_yyyzzz_0_xxyzzzz_1, g_0_yyyzzz_0_xyyyyyz_0, g_0_yyyzzz_0_xyyyyyz_1, g_0_yyyzzz_0_xyyyyzz_0, g_0_yyyzzz_0_xyyyyzz_1, g_0_yyyzzz_0_xyyyzzz_0, g_0_yyyzzz_0_xyyyzzz_1, g_0_yyyzzz_0_xyyzzzz_0, g_0_yyyzzz_0_xyyzzzz_1, g_0_yyyzzz_0_xyzzzzz_0, g_0_yyyzzz_0_xyzzzzz_1, g_0_yyyzzz_0_yyyyyyy_0, g_0_yyyzzz_0_yyyyyyy_1, g_0_yyyzzz_0_yyyyyyz_0, g_0_yyyzzz_0_yyyyyyz_1, g_0_yyyzzz_0_yyyyyzz_0, g_0_yyyzzz_0_yyyyyzz_1, g_0_yyyzzz_0_yyyyzzz_0, g_0_yyyzzz_0_yyyyzzz_1, g_0_yyyzzz_0_yyyzzzz_0, g_0_yyyzzz_0_yyyzzzz_1, g_0_yyyzzz_0_yyzzzzz_0, g_0_yyyzzz_0_yyzzzzz_1, g_0_yyyzzz_0_yzzzzzz_0, g_0_yyyzzz_0_yzzzzzz_1, g_0_yyyzzz_0_zzzzzzz_0, g_0_yyyzzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzzz_0_xxxxxxx_0[i] = 2.0 * g_0_xxyzzz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxyyzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxxxxy_0[i] = 2.0 * g_0_xxyyyz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxxxy_0[i] * pb_z + g_0_xxyyyzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxxxxz_0[i] = 2.0 * g_0_xxyzzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxyyzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxyyyz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxxyy_0[i] * pb_z + g_0_xxyyyzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxxxyz_0[i] = g_0_yyyzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxxxyz_0[i] * pb_x + g_0_xyyyzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxxxzz_0[i] = 2.0 * g_0_xxyzzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxyyzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxxyyy_0[i] = 2.0 * g_0_xxyyyz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxyyy_0[i] * pb_z + g_0_xxyyyzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxxyyz_0[i] = g_0_yyyzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxxyyz_0[i] * pb_x + g_0_xyyyzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxxyzz_0[i] = g_0_yyyzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxxyzz_0[i] * pb_x + g_0_xyyyzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxxzzz_0[i] = 2.0 * g_0_xxyzzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxyyzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxyyyy_0[i] = 2.0 * g_0_xxyyyz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxyyyy_0[i] * pb_z + g_0_xxyyyzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxyyyz_0[i] = g_0_yyyzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxyyyz_0[i] * pb_x + g_0_xyyyzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxyyzz_0[i] = g_0_yyyzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxyyzz_0[i] * pb_x + g_0_xyyyzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxyzzz_0[i] = g_0_yyyzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxyzzz_0[i] * pb_x + g_0_xyyyzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxzzzz_0[i] = 2.0 * g_0_xxyzzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxyyzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_xxyyyz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxyyyyy_0[i] * pb_z + g_0_xxyyyzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxyyyyz_0[i] = g_0_yyyzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyyyyz_0[i] * pb_x + g_0_xyyyzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxyyyzz_0[i] = g_0_yyyzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyyyzz_0[i] * pb_x + g_0_xyyyzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxyyzzz_0[i] = g_0_yyyzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyyzzz_0[i] * pb_x + g_0_xyyyzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxyzzzz_0[i] = g_0_yyyzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyzzzz_0[i] * pb_x + g_0_xyyyzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxzzzzz_0[i] = 2.0 * g_0_xxyzzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxyyzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xyyyyyy_0[i] = 2.0 * g_0_xxyyyz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xyyyyyy_0[i] * pb_z + g_0_xxyyyzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xyyyyyz_0[i] = g_0_yyyzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyyyyyz_0[i] * pb_x + g_0_xyyyzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyyyyzz_0[i] = g_0_yyyzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyyyyzz_0[i] * pb_x + g_0_xyyyzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyyyzzz_0[i] = g_0_yyyzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyyyzzz_0[i] * pb_x + g_0_xyyyzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyyzzzz_0[i] = g_0_yyyzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyyzzzz_0[i] * pb_x + g_0_xyyyzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyzzzzz_0[i] = g_0_yyyzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyzzzzz_0[i] * pb_x + g_0_xyyyzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xzzzzzz_0[i] = 2.0 * g_0_xxyzzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxyyzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_yyyyyyy_0[i] = g_0_yyyzzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyyyy_0[i] * pb_x + g_0_xyyyzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyyyyz_0[i] = g_0_yyyzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyyyz_0[i] * pb_x + g_0_xyyyzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyyyzz_0[i] = g_0_yyyzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyyzz_0[i] * pb_x + g_0_xyyyzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyyzzz_0[i] = g_0_yyyzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyzzz_0[i] * pb_x + g_0_xyyyzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyzzzz_0[i] = g_0_yyyzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyzzzz_0[i] * pb_x + g_0_xyyyzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyzzzzz_0[i] = g_0_yyyzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyzzzzz_0[i] * pb_x + g_0_xyyyzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yzzzzzz_0[i] = g_0_yyyzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzzzzzz_0[i] * pb_x + g_0_xyyyzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_zzzzzzz_0[i] = g_0_yyyzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_zzzzzzz_0[i] * pb_x + g_0_xyyyzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 900-936 components of targeted buffer : SLSK

    auto g_0_xxyyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 900);

    auto g_0_xxyyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 901);

    auto g_0_xxyyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 902);

    auto g_0_xxyyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 903);

    auto g_0_xxyyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 904);

    auto g_0_xxyyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 905);

    auto g_0_xxyyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 906);

    auto g_0_xxyyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 907);

    auto g_0_xxyyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 908);

    auto g_0_xxyyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 909);

    auto g_0_xxyyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 910);

    auto g_0_xxyyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 911);

    auto g_0_xxyyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 912);

    auto g_0_xxyyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 913);

    auto g_0_xxyyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 914);

    auto g_0_xxyyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 915);

    auto g_0_xxyyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 916);

    auto g_0_xxyyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 917);

    auto g_0_xxyyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 918);

    auto g_0_xxyyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 919);

    auto g_0_xxyyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 920);

    auto g_0_xxyyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 921);

    auto g_0_xxyyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 922);

    auto g_0_xxyyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 923);

    auto g_0_xxyyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 924);

    auto g_0_xxyyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 925);

    auto g_0_xxyyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 926);

    auto g_0_xxyyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 927);

    auto g_0_xxyyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 928);

    auto g_0_xxyyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 929);

    auto g_0_xxyyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 930);

    auto g_0_xxyyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 931);

    auto g_0_xxyyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 932);

    auto g_0_xxyyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 933);

    auto g_0_xxyyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 934);

    auto g_0_xxyyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 935);

    #pragma omp simd aligned(g_0_xxyyzz_0_xxxxxxy_0, g_0_xxyyzz_0_xxxxxxy_1, g_0_xxyyzz_0_xxxxxyy_0, g_0_xxyyzz_0_xxxxxyy_1, g_0_xxyyzz_0_xxxxyyy_0, g_0_xxyyzz_0_xxxxyyy_1, g_0_xxyyzz_0_xxxyyyy_0, g_0_xxyyzz_0_xxxyyyy_1, g_0_xxyyzz_0_xxyyyyy_0, g_0_xxyyzz_0_xxyyyyy_1, g_0_xxyyzz_0_xyyyyyy_0, g_0_xxyyzz_0_xyyyyyy_1, g_0_xxyyzzz_0_xxxxxxy_0, g_0_xxyyzzz_0_xxxxxxy_1, g_0_xxyyzzz_0_xxxxxyy_0, g_0_xxyyzzz_0_xxxxxyy_1, g_0_xxyyzzz_0_xxxxyyy_0, g_0_xxyyzzz_0_xxxxyyy_1, g_0_xxyyzzz_0_xxxyyyy_0, g_0_xxyyzzz_0_xxxyyyy_1, g_0_xxyyzzz_0_xxyyyyy_0, g_0_xxyyzzz_0_xxyyyyy_1, g_0_xxyyzzz_0_xyyyyyy_0, g_0_xxyyzzz_0_xyyyyyy_1, g_0_xxyyzzzz_0_xxxxxxx_0, g_0_xxyyzzzz_0_xxxxxxy_0, g_0_xxyyzzzz_0_xxxxxxz_0, g_0_xxyyzzzz_0_xxxxxyy_0, g_0_xxyyzzzz_0_xxxxxyz_0, g_0_xxyyzzzz_0_xxxxxzz_0, g_0_xxyyzzzz_0_xxxxyyy_0, g_0_xxyyzzzz_0_xxxxyyz_0, g_0_xxyyzzzz_0_xxxxyzz_0, g_0_xxyyzzzz_0_xxxxzzz_0, g_0_xxyyzzzz_0_xxxyyyy_0, g_0_xxyyzzzz_0_xxxyyyz_0, g_0_xxyyzzzz_0_xxxyyzz_0, g_0_xxyyzzzz_0_xxxyzzz_0, g_0_xxyyzzzz_0_xxxzzzz_0, g_0_xxyyzzzz_0_xxyyyyy_0, g_0_xxyyzzzz_0_xxyyyyz_0, g_0_xxyyzzzz_0_xxyyyzz_0, g_0_xxyyzzzz_0_xxyyzzz_0, g_0_xxyyzzzz_0_xxyzzzz_0, g_0_xxyyzzzz_0_xxzzzzz_0, g_0_xxyyzzzz_0_xyyyyyy_0, g_0_xxyyzzzz_0_xyyyyyz_0, g_0_xxyyzzzz_0_xyyyyzz_0, g_0_xxyyzzzz_0_xyyyzzz_0, g_0_xxyyzzzz_0_xyyzzzz_0, g_0_xxyyzzzz_0_xyzzzzz_0, g_0_xxyyzzzz_0_xzzzzzz_0, g_0_xxyyzzzz_0_yyyyyyy_0, g_0_xxyyzzzz_0_yyyyyyz_0, g_0_xxyyzzzz_0_yyyyyzz_0, g_0_xxyyzzzz_0_yyyyzzz_0, g_0_xxyyzzzz_0_yyyzzzz_0, g_0_xxyyzzzz_0_yyzzzzz_0, g_0_xxyyzzzz_0_yzzzzzz_0, g_0_xxyyzzzz_0_zzzzzzz_0, g_0_xxyzzzz_0_xxxxxxx_0, g_0_xxyzzzz_0_xxxxxxx_1, g_0_xxyzzzz_0_xxxxxxz_0, g_0_xxyzzzz_0_xxxxxxz_1, g_0_xxyzzzz_0_xxxxxzz_0, g_0_xxyzzzz_0_xxxxxzz_1, g_0_xxyzzzz_0_xxxxzzz_0, g_0_xxyzzzz_0_xxxxzzz_1, g_0_xxyzzzz_0_xxxzzzz_0, g_0_xxyzzzz_0_xxxzzzz_1, g_0_xxyzzzz_0_xxzzzzz_0, g_0_xxyzzzz_0_xxzzzzz_1, g_0_xxyzzzz_0_xzzzzzz_0, g_0_xxyzzzz_0_xzzzzzz_1, g_0_xxzzzz_0_xxxxxxx_0, g_0_xxzzzz_0_xxxxxxx_1, g_0_xxzzzz_0_xxxxxxz_0, g_0_xxzzzz_0_xxxxxxz_1, g_0_xxzzzz_0_xxxxxzz_0, g_0_xxzzzz_0_xxxxxzz_1, g_0_xxzzzz_0_xxxxzzz_0, g_0_xxzzzz_0_xxxxzzz_1, g_0_xxzzzz_0_xxxzzzz_0, g_0_xxzzzz_0_xxxzzzz_1, g_0_xxzzzz_0_xxzzzzz_0, g_0_xxzzzz_0_xxzzzzz_1, g_0_xxzzzz_0_xzzzzzz_0, g_0_xxzzzz_0_xzzzzzz_1, g_0_xyyzzzz_0_xxxxxyz_0, g_0_xyyzzzz_0_xxxxxyz_1, g_0_xyyzzzz_0_xxxxyyz_0, g_0_xyyzzzz_0_xxxxyyz_1, g_0_xyyzzzz_0_xxxxyz_1, g_0_xyyzzzz_0_xxxxyzz_0, g_0_xyyzzzz_0_xxxxyzz_1, g_0_xyyzzzz_0_xxxyyyz_0, g_0_xyyzzzz_0_xxxyyyz_1, g_0_xyyzzzz_0_xxxyyz_1, g_0_xyyzzzz_0_xxxyyzz_0, g_0_xyyzzzz_0_xxxyyzz_1, g_0_xyyzzzz_0_xxxyzz_1, g_0_xyyzzzz_0_xxxyzzz_0, g_0_xyyzzzz_0_xxxyzzz_1, g_0_xyyzzzz_0_xxyyyyz_0, g_0_xyyzzzz_0_xxyyyyz_1, g_0_xyyzzzz_0_xxyyyz_1, g_0_xyyzzzz_0_xxyyyzz_0, g_0_xyyzzzz_0_xxyyyzz_1, g_0_xyyzzzz_0_xxyyzz_1, g_0_xyyzzzz_0_xxyyzzz_0, g_0_xyyzzzz_0_xxyyzzz_1, g_0_xyyzzzz_0_xxyzzz_1, g_0_xyyzzzz_0_xxyzzzz_0, g_0_xyyzzzz_0_xxyzzzz_1, g_0_xyyzzzz_0_xyyyyyz_0, g_0_xyyzzzz_0_xyyyyyz_1, g_0_xyyzzzz_0_xyyyyz_1, g_0_xyyzzzz_0_xyyyyzz_0, g_0_xyyzzzz_0_xyyyyzz_1, g_0_xyyzzzz_0_xyyyzz_1, g_0_xyyzzzz_0_xyyyzzz_0, g_0_xyyzzzz_0_xyyyzzz_1, g_0_xyyzzzz_0_xyyzzz_1, g_0_xyyzzzz_0_xyyzzzz_0, g_0_xyyzzzz_0_xyyzzzz_1, g_0_xyyzzzz_0_xyzzzz_1, g_0_xyyzzzz_0_xyzzzzz_0, g_0_xyyzzzz_0_xyzzzzz_1, g_0_xyyzzzz_0_yyyyyyy_0, g_0_xyyzzzz_0_yyyyyyy_1, g_0_xyyzzzz_0_yyyyyyz_0, g_0_xyyzzzz_0_yyyyyyz_1, g_0_xyyzzzz_0_yyyyyz_1, g_0_xyyzzzz_0_yyyyyzz_0, g_0_xyyzzzz_0_yyyyyzz_1, g_0_xyyzzzz_0_yyyyzz_1, g_0_xyyzzzz_0_yyyyzzz_0, g_0_xyyzzzz_0_yyyyzzz_1, g_0_xyyzzzz_0_yyyzzz_1, g_0_xyyzzzz_0_yyyzzzz_0, g_0_xyyzzzz_0_yyyzzzz_1, g_0_xyyzzzz_0_yyzzzz_1, g_0_xyyzzzz_0_yyzzzzz_0, g_0_xyyzzzz_0_yyzzzzz_1, g_0_xyyzzzz_0_yzzzzz_1, g_0_xyyzzzz_0_yzzzzzz_0, g_0_xyyzzzz_0_yzzzzzz_1, g_0_xyyzzzz_0_zzzzzzz_0, g_0_xyyzzzz_0_zzzzzzz_1, g_0_yyzzzz_0_xxxxxyz_0, g_0_yyzzzz_0_xxxxxyz_1, g_0_yyzzzz_0_xxxxyyz_0, g_0_yyzzzz_0_xxxxyyz_1, g_0_yyzzzz_0_xxxxyzz_0, g_0_yyzzzz_0_xxxxyzz_1, g_0_yyzzzz_0_xxxyyyz_0, g_0_yyzzzz_0_xxxyyyz_1, g_0_yyzzzz_0_xxxyyzz_0, g_0_yyzzzz_0_xxxyyzz_1, g_0_yyzzzz_0_xxxyzzz_0, g_0_yyzzzz_0_xxxyzzz_1, g_0_yyzzzz_0_xxyyyyz_0, g_0_yyzzzz_0_xxyyyyz_1, g_0_yyzzzz_0_xxyyyzz_0, g_0_yyzzzz_0_xxyyyzz_1, g_0_yyzzzz_0_xxyyzzz_0, g_0_yyzzzz_0_xxyyzzz_1, g_0_yyzzzz_0_xxyzzzz_0, g_0_yyzzzz_0_xxyzzzz_1, g_0_yyzzzz_0_xyyyyyz_0, g_0_yyzzzz_0_xyyyyyz_1, g_0_yyzzzz_0_xyyyyzz_0, g_0_yyzzzz_0_xyyyyzz_1, g_0_yyzzzz_0_xyyyzzz_0, g_0_yyzzzz_0_xyyyzzz_1, g_0_yyzzzz_0_xyyzzzz_0, g_0_yyzzzz_0_xyyzzzz_1, g_0_yyzzzz_0_xyzzzzz_0, g_0_yyzzzz_0_xyzzzzz_1, g_0_yyzzzz_0_yyyyyyy_0, g_0_yyzzzz_0_yyyyyyy_1, g_0_yyzzzz_0_yyyyyyz_0, g_0_yyzzzz_0_yyyyyyz_1, g_0_yyzzzz_0_yyyyyzz_0, g_0_yyzzzz_0_yyyyyzz_1, g_0_yyzzzz_0_yyyyzzz_0, g_0_yyzzzz_0_yyyyzzz_1, g_0_yyzzzz_0_yyyzzzz_0, g_0_yyzzzz_0_yyyzzzz_1, g_0_yyzzzz_0_yyzzzzz_0, g_0_yyzzzz_0_yyzzzzz_1, g_0_yyzzzz_0_yzzzzzz_0, g_0_yyzzzz_0_yzzzzzz_1, g_0_yyzzzz_0_zzzzzzz_0, g_0_yyzzzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzzz_0_xxxxxxx_0[i] = g_0_xxzzzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxyzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxxxxy_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxxxy_0[i] * pb_z + g_0_xxyyzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxxxxz_0[i] = g_0_xxzzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxyzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxxxyy_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxxyy_0[i] * pb_z + g_0_xxyyzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxxxyz_0[i] = g_0_yyzzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxxxyz_0[i] * pb_x + g_0_xyyzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxxxzz_0[i] = g_0_xxzzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxyzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxyyzz_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxyyy_0[i] * pb_z + g_0_xxyyzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxxyyz_0[i] = g_0_yyzzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxxyyz_0[i] * pb_x + g_0_xyyzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxxyzz_0[i] = g_0_yyzzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxxyzz_0[i] * pb_x + g_0_xyyzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxxzzz_0[i] = g_0_xxzzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxyzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_xxyyzz_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxyyyy_0[i] * pb_z + g_0_xxyyzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxyyyz_0[i] = g_0_yyzzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxyyyz_0[i] * pb_x + g_0_xyyzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxyyzz_0[i] = g_0_yyzzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxyyzz_0[i] * pb_x + g_0_xyyzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxyzzz_0[i] = g_0_yyzzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxyzzz_0[i] * pb_x + g_0_xyyzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxzzzz_0[i] = g_0_xxzzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxyzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxyyyyy_0[i] = 3.0 * g_0_xxyyzz_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxyyyyy_0[i] * pb_z + g_0_xxyyzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxyyyyz_0[i] = g_0_yyzzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyyyyz_0[i] * pb_x + g_0_xyyzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxyyyzz_0[i] = g_0_yyzzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyyyzz_0[i] * pb_x + g_0_xyyzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxyyzzz_0[i] = g_0_yyzzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyyzzz_0[i] * pb_x + g_0_xyyzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxyzzzz_0[i] = g_0_yyzzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyzzzz_0[i] * pb_x + g_0_xyyzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxzzzzz_0[i] = g_0_xxzzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxyzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xyyyyyy_0[i] = 3.0 * g_0_xxyyzz_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xyyyyyy_0[i] * pb_z + g_0_xxyyzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xyyyyyz_0[i] = g_0_yyzzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyyyyyz_0[i] * pb_x + g_0_xyyzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyyyyzz_0[i] = g_0_yyzzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyyyyzz_0[i] * pb_x + g_0_xyyzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyyyzzz_0[i] = g_0_yyzzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyyyzzz_0[i] * pb_x + g_0_xyyzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyyzzzz_0[i] = g_0_yyzzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyyzzzz_0[i] * pb_x + g_0_xyyzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyzzzzz_0[i] = g_0_yyzzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyzzzzz_0[i] * pb_x + g_0_xyyzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xzzzzzz_0[i] = g_0_xxzzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxyzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_yyyyyyy_0[i] = g_0_yyzzzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyyyy_0[i] * pb_x + g_0_xyyzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyyyyz_0[i] = g_0_yyzzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyyyz_0[i] * pb_x + g_0_xyyzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyyyzz_0[i] = g_0_yyzzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyyzz_0[i] * pb_x + g_0_xyyzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyyzzz_0[i] = g_0_yyzzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyzzz_0[i] * pb_x + g_0_xyyzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyzzzz_0[i] = g_0_yyzzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyzzzz_0[i] * pb_x + g_0_xyyzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyzzzzz_0[i] = g_0_yyzzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyzzzzz_0[i] * pb_x + g_0_xyyzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yzzzzzz_0[i] = g_0_yyzzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzzzzzz_0[i] * pb_x + g_0_xyyzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_zzzzzzz_0[i] = g_0_yyzzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_zzzzzzz_0[i] * pb_x + g_0_xyyzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 936-972 components of targeted buffer : SLSK

    auto g_0_xxyzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 936);

    auto g_0_xxyzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 937);

    auto g_0_xxyzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 938);

    auto g_0_xxyzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 939);

    auto g_0_xxyzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 940);

    auto g_0_xxyzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 941);

    auto g_0_xxyzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 942);

    auto g_0_xxyzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 943);

    auto g_0_xxyzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 944);

    auto g_0_xxyzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 945);

    auto g_0_xxyzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 946);

    auto g_0_xxyzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 947);

    auto g_0_xxyzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 948);

    auto g_0_xxyzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 949);

    auto g_0_xxyzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 950);

    auto g_0_xxyzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 951);

    auto g_0_xxyzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 952);

    auto g_0_xxyzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 953);

    auto g_0_xxyzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 954);

    auto g_0_xxyzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 955);

    auto g_0_xxyzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 956);

    auto g_0_xxyzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 957);

    auto g_0_xxyzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 958);

    auto g_0_xxyzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 959);

    auto g_0_xxyzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 960);

    auto g_0_xxyzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 961);

    auto g_0_xxyzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 962);

    auto g_0_xxyzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 963);

    auto g_0_xxyzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 964);

    auto g_0_xxyzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 965);

    auto g_0_xxyzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 966);

    auto g_0_xxyzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 967);

    auto g_0_xxyzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 968);

    auto g_0_xxyzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 969);

    auto g_0_xxyzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 970);

    auto g_0_xxyzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 971);

    #pragma omp simd aligned(g_0_xxyzzzzz_0_xxxxxxx_0, g_0_xxyzzzzz_0_xxxxxxy_0, g_0_xxyzzzzz_0_xxxxxxz_0, g_0_xxyzzzzz_0_xxxxxyy_0, g_0_xxyzzzzz_0_xxxxxyz_0, g_0_xxyzzzzz_0_xxxxxzz_0, g_0_xxyzzzzz_0_xxxxyyy_0, g_0_xxyzzzzz_0_xxxxyyz_0, g_0_xxyzzzzz_0_xxxxyzz_0, g_0_xxyzzzzz_0_xxxxzzz_0, g_0_xxyzzzzz_0_xxxyyyy_0, g_0_xxyzzzzz_0_xxxyyyz_0, g_0_xxyzzzzz_0_xxxyyzz_0, g_0_xxyzzzzz_0_xxxyzzz_0, g_0_xxyzzzzz_0_xxxzzzz_0, g_0_xxyzzzzz_0_xxyyyyy_0, g_0_xxyzzzzz_0_xxyyyyz_0, g_0_xxyzzzzz_0_xxyyyzz_0, g_0_xxyzzzzz_0_xxyyzzz_0, g_0_xxyzzzzz_0_xxyzzzz_0, g_0_xxyzzzzz_0_xxzzzzz_0, g_0_xxyzzzzz_0_xyyyyyy_0, g_0_xxyzzzzz_0_xyyyyyz_0, g_0_xxyzzzzz_0_xyyyyzz_0, g_0_xxyzzzzz_0_xyyyzzz_0, g_0_xxyzzzzz_0_xyyzzzz_0, g_0_xxyzzzzz_0_xyzzzzz_0, g_0_xxyzzzzz_0_xzzzzzz_0, g_0_xxyzzzzz_0_yyyyyyy_0, g_0_xxyzzzzz_0_yyyyyyz_0, g_0_xxyzzzzz_0_yyyyyzz_0, g_0_xxyzzzzz_0_yyyyzzz_0, g_0_xxyzzzzz_0_yyyzzzz_0, g_0_xxyzzzzz_0_yyzzzzz_0, g_0_xxyzzzzz_0_yzzzzzz_0, g_0_xxyzzzzz_0_zzzzzzz_0, g_0_xxzzzzz_0_xxxxxx_1, g_0_xxzzzzz_0_xxxxxxx_0, g_0_xxzzzzz_0_xxxxxxx_1, g_0_xxzzzzz_0_xxxxxxy_0, g_0_xxzzzzz_0_xxxxxxy_1, g_0_xxzzzzz_0_xxxxxxz_0, g_0_xxzzzzz_0_xxxxxxz_1, g_0_xxzzzzz_0_xxxxxy_1, g_0_xxzzzzz_0_xxxxxyy_0, g_0_xxzzzzz_0_xxxxxyy_1, g_0_xxzzzzz_0_xxxxxyz_0, g_0_xxzzzzz_0_xxxxxyz_1, g_0_xxzzzzz_0_xxxxxz_1, g_0_xxzzzzz_0_xxxxxzz_0, g_0_xxzzzzz_0_xxxxxzz_1, g_0_xxzzzzz_0_xxxxyy_1, g_0_xxzzzzz_0_xxxxyyy_0, g_0_xxzzzzz_0_xxxxyyy_1, g_0_xxzzzzz_0_xxxxyyz_0, g_0_xxzzzzz_0_xxxxyyz_1, g_0_xxzzzzz_0_xxxxyz_1, g_0_xxzzzzz_0_xxxxyzz_0, g_0_xxzzzzz_0_xxxxyzz_1, g_0_xxzzzzz_0_xxxxzz_1, g_0_xxzzzzz_0_xxxxzzz_0, g_0_xxzzzzz_0_xxxxzzz_1, g_0_xxzzzzz_0_xxxyyy_1, g_0_xxzzzzz_0_xxxyyyy_0, g_0_xxzzzzz_0_xxxyyyy_1, g_0_xxzzzzz_0_xxxyyyz_0, g_0_xxzzzzz_0_xxxyyyz_1, g_0_xxzzzzz_0_xxxyyz_1, g_0_xxzzzzz_0_xxxyyzz_0, g_0_xxzzzzz_0_xxxyyzz_1, g_0_xxzzzzz_0_xxxyzz_1, g_0_xxzzzzz_0_xxxyzzz_0, g_0_xxzzzzz_0_xxxyzzz_1, g_0_xxzzzzz_0_xxxzzz_1, g_0_xxzzzzz_0_xxxzzzz_0, g_0_xxzzzzz_0_xxxzzzz_1, g_0_xxzzzzz_0_xxyyyy_1, g_0_xxzzzzz_0_xxyyyyy_0, g_0_xxzzzzz_0_xxyyyyy_1, g_0_xxzzzzz_0_xxyyyyz_0, g_0_xxzzzzz_0_xxyyyyz_1, g_0_xxzzzzz_0_xxyyyz_1, g_0_xxzzzzz_0_xxyyyzz_0, g_0_xxzzzzz_0_xxyyyzz_1, g_0_xxzzzzz_0_xxyyzz_1, g_0_xxzzzzz_0_xxyyzzz_0, g_0_xxzzzzz_0_xxyyzzz_1, g_0_xxzzzzz_0_xxyzzz_1, g_0_xxzzzzz_0_xxyzzzz_0, g_0_xxzzzzz_0_xxyzzzz_1, g_0_xxzzzzz_0_xxzzzz_1, g_0_xxzzzzz_0_xxzzzzz_0, g_0_xxzzzzz_0_xxzzzzz_1, g_0_xxzzzzz_0_xyyyyy_1, g_0_xxzzzzz_0_xyyyyyy_0, g_0_xxzzzzz_0_xyyyyyy_1, g_0_xxzzzzz_0_xyyyyyz_0, g_0_xxzzzzz_0_xyyyyyz_1, g_0_xxzzzzz_0_xyyyyz_1, g_0_xxzzzzz_0_xyyyyzz_0, g_0_xxzzzzz_0_xyyyyzz_1, g_0_xxzzzzz_0_xyyyzz_1, g_0_xxzzzzz_0_xyyyzzz_0, g_0_xxzzzzz_0_xyyyzzz_1, g_0_xxzzzzz_0_xyyzzz_1, g_0_xxzzzzz_0_xyyzzzz_0, g_0_xxzzzzz_0_xyyzzzz_1, g_0_xxzzzzz_0_xyzzzz_1, g_0_xxzzzzz_0_xyzzzzz_0, g_0_xxzzzzz_0_xyzzzzz_1, g_0_xxzzzzz_0_xzzzzz_1, g_0_xxzzzzz_0_xzzzzzz_0, g_0_xxzzzzz_0_xzzzzzz_1, g_0_xxzzzzz_0_yyyyyy_1, g_0_xxzzzzz_0_yyyyyyy_0, g_0_xxzzzzz_0_yyyyyyy_1, g_0_xxzzzzz_0_yyyyyyz_0, g_0_xxzzzzz_0_yyyyyyz_1, g_0_xxzzzzz_0_yyyyyz_1, g_0_xxzzzzz_0_yyyyyzz_0, g_0_xxzzzzz_0_yyyyyzz_1, g_0_xxzzzzz_0_yyyyzz_1, g_0_xxzzzzz_0_yyyyzzz_0, g_0_xxzzzzz_0_yyyyzzz_1, g_0_xxzzzzz_0_yyyzzz_1, g_0_xxzzzzz_0_yyyzzzz_0, g_0_xxzzzzz_0_yyyzzzz_1, g_0_xxzzzzz_0_yyzzzz_1, g_0_xxzzzzz_0_yyzzzzz_0, g_0_xxzzzzz_0_yyzzzzz_1, g_0_xxzzzzz_0_yzzzzz_1, g_0_xxzzzzz_0_yzzzzzz_0, g_0_xxzzzzz_0_yzzzzzz_1, g_0_xxzzzzz_0_zzzzzz_1, g_0_xxzzzzz_0_zzzzzzz_0, g_0_xxzzzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzzz_0_xxxxxxx_0[i] = g_0_xxzzzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxxxy_0[i] = g_0_xxzzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxxxy_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxxxz_0[i] = g_0_xxzzzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxxyy_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxxyz_0[i] = g_0_xxzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxxyz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxxzz_0[i] = g_0_xxzzzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxyyy_0[i] * pb_y + g_0_xxzzzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxyyz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxyzz_0[i] = g_0_xxzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxyzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxzzz_0[i] = g_0_xxzzzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyyyy_0[i] * pb_y + g_0_xxzzzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyyyz_0[i] * pb_y + g_0_xxzzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyyzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyzzz_0[i] = g_0_xxzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxzzzz_0[i] = g_0_xxzzzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyyyy_0[i] * pb_y + g_0_xxzzzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyyyz_0[i] * pb_y + g_0_xxzzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyyzz_0[i] * pb_y + g_0_xxzzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyzzzz_0[i] = g_0_xxzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxzzzzz_0[i] = g_0_xxzzzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyyyy_0[i] * pb_y + g_0_xxzzzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyyyz_0[i] * pb_y + g_0_xxzzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyyzz_0[i] * pb_y + g_0_xxzzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyzzz_0[i] * pb_y + g_0_xxzzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyzzzzz_0[i] = g_0_xxzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xzzzzzz_0[i] = g_0_xxzzzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyyyy_0[i] * pb_y + g_0_xxzzzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyyyz_0[i] * pb_y + g_0_xxzzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyyzz_0[i] * pb_y + g_0_xxzzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyzzz_0[i] * pb_y + g_0_xxzzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyzzzz_0[i] * pb_y + g_0_xxzzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yzzzzzz_0[i] = g_0_xxzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yzzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_zzzzzzz_0[i] = g_0_xxzzzzz_0_zzzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 972-1008 components of targeted buffer : SLSK

    auto g_0_xxzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 972);

    auto g_0_xxzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 973);

    auto g_0_xxzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 974);

    auto g_0_xxzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 975);

    auto g_0_xxzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 976);

    auto g_0_xxzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 977);

    auto g_0_xxzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 978);

    auto g_0_xxzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 979);

    auto g_0_xxzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 980);

    auto g_0_xxzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 981);

    auto g_0_xxzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 982);

    auto g_0_xxzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 983);

    auto g_0_xxzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 984);

    auto g_0_xxzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 985);

    auto g_0_xxzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 986);

    auto g_0_xxzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 987);

    auto g_0_xxzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 988);

    auto g_0_xxzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 989);

    auto g_0_xxzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 990);

    auto g_0_xxzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 991);

    auto g_0_xxzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 992);

    auto g_0_xxzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 993);

    auto g_0_xxzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 994);

    auto g_0_xxzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 995);

    auto g_0_xxzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 996);

    auto g_0_xxzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 997);

    auto g_0_xxzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 998);

    auto g_0_xxzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 999);

    auto g_0_xxzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1000);

    auto g_0_xxzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1001);

    auto g_0_xxzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1002);

    auto g_0_xxzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1003);

    auto g_0_xxzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1004);

    auto g_0_xxzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1005);

    auto g_0_xxzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1006);

    auto g_0_xxzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1007);

    #pragma omp simd aligned(g_0_xxzzzz_0_xxxxxxx_0, g_0_xxzzzz_0_xxxxxxx_1, g_0_xxzzzz_0_xxxxxxy_0, g_0_xxzzzz_0_xxxxxxy_1, g_0_xxzzzz_0_xxxxxyy_0, g_0_xxzzzz_0_xxxxxyy_1, g_0_xxzzzz_0_xxxxyyy_0, g_0_xxzzzz_0_xxxxyyy_1, g_0_xxzzzz_0_xxxyyyy_0, g_0_xxzzzz_0_xxxyyyy_1, g_0_xxzzzz_0_xxyyyyy_0, g_0_xxzzzz_0_xxyyyyy_1, g_0_xxzzzz_0_xyyyyyy_0, g_0_xxzzzz_0_xyyyyyy_1, g_0_xxzzzzz_0_xxxxxxx_0, g_0_xxzzzzz_0_xxxxxxx_1, g_0_xxzzzzz_0_xxxxxxy_0, g_0_xxzzzzz_0_xxxxxxy_1, g_0_xxzzzzz_0_xxxxxyy_0, g_0_xxzzzzz_0_xxxxxyy_1, g_0_xxzzzzz_0_xxxxyyy_0, g_0_xxzzzzz_0_xxxxyyy_1, g_0_xxzzzzz_0_xxxyyyy_0, g_0_xxzzzzz_0_xxxyyyy_1, g_0_xxzzzzz_0_xxyyyyy_0, g_0_xxzzzzz_0_xxyyyyy_1, g_0_xxzzzzz_0_xyyyyyy_0, g_0_xxzzzzz_0_xyyyyyy_1, g_0_xxzzzzzz_0_xxxxxxx_0, g_0_xxzzzzzz_0_xxxxxxy_0, g_0_xxzzzzzz_0_xxxxxxz_0, g_0_xxzzzzzz_0_xxxxxyy_0, g_0_xxzzzzzz_0_xxxxxyz_0, g_0_xxzzzzzz_0_xxxxxzz_0, g_0_xxzzzzzz_0_xxxxyyy_0, g_0_xxzzzzzz_0_xxxxyyz_0, g_0_xxzzzzzz_0_xxxxyzz_0, g_0_xxzzzzzz_0_xxxxzzz_0, g_0_xxzzzzzz_0_xxxyyyy_0, g_0_xxzzzzzz_0_xxxyyyz_0, g_0_xxzzzzzz_0_xxxyyzz_0, g_0_xxzzzzzz_0_xxxyzzz_0, g_0_xxzzzzzz_0_xxxzzzz_0, g_0_xxzzzzzz_0_xxyyyyy_0, g_0_xxzzzzzz_0_xxyyyyz_0, g_0_xxzzzzzz_0_xxyyyzz_0, g_0_xxzzzzzz_0_xxyyzzz_0, g_0_xxzzzzzz_0_xxyzzzz_0, g_0_xxzzzzzz_0_xxzzzzz_0, g_0_xxzzzzzz_0_xyyyyyy_0, g_0_xxzzzzzz_0_xyyyyyz_0, g_0_xxzzzzzz_0_xyyyyzz_0, g_0_xxzzzzzz_0_xyyyzzz_0, g_0_xxzzzzzz_0_xyyzzzz_0, g_0_xxzzzzzz_0_xyzzzzz_0, g_0_xxzzzzzz_0_xzzzzzz_0, g_0_xxzzzzzz_0_yyyyyyy_0, g_0_xxzzzzzz_0_yyyyyyz_0, g_0_xxzzzzzz_0_yyyyyzz_0, g_0_xxzzzzzz_0_yyyyzzz_0, g_0_xxzzzzzz_0_yyyzzzz_0, g_0_xxzzzzzz_0_yyzzzzz_0, g_0_xxzzzzzz_0_yzzzzzz_0, g_0_xxzzzzzz_0_zzzzzzz_0, g_0_xzzzzzz_0_xxxxxxz_0, g_0_xzzzzzz_0_xxxxxxz_1, g_0_xzzzzzz_0_xxxxxyz_0, g_0_xzzzzzz_0_xxxxxyz_1, g_0_xzzzzzz_0_xxxxxz_1, g_0_xzzzzzz_0_xxxxxzz_0, g_0_xzzzzzz_0_xxxxxzz_1, g_0_xzzzzzz_0_xxxxyyz_0, g_0_xzzzzzz_0_xxxxyyz_1, g_0_xzzzzzz_0_xxxxyz_1, g_0_xzzzzzz_0_xxxxyzz_0, g_0_xzzzzzz_0_xxxxyzz_1, g_0_xzzzzzz_0_xxxxzz_1, g_0_xzzzzzz_0_xxxxzzz_0, g_0_xzzzzzz_0_xxxxzzz_1, g_0_xzzzzzz_0_xxxyyyz_0, g_0_xzzzzzz_0_xxxyyyz_1, g_0_xzzzzzz_0_xxxyyz_1, g_0_xzzzzzz_0_xxxyyzz_0, g_0_xzzzzzz_0_xxxyyzz_1, g_0_xzzzzzz_0_xxxyzz_1, g_0_xzzzzzz_0_xxxyzzz_0, g_0_xzzzzzz_0_xxxyzzz_1, g_0_xzzzzzz_0_xxxzzz_1, g_0_xzzzzzz_0_xxxzzzz_0, g_0_xzzzzzz_0_xxxzzzz_1, g_0_xzzzzzz_0_xxyyyyz_0, g_0_xzzzzzz_0_xxyyyyz_1, g_0_xzzzzzz_0_xxyyyz_1, g_0_xzzzzzz_0_xxyyyzz_0, g_0_xzzzzzz_0_xxyyyzz_1, g_0_xzzzzzz_0_xxyyzz_1, g_0_xzzzzzz_0_xxyyzzz_0, g_0_xzzzzzz_0_xxyyzzz_1, g_0_xzzzzzz_0_xxyzzz_1, g_0_xzzzzzz_0_xxyzzzz_0, g_0_xzzzzzz_0_xxyzzzz_1, g_0_xzzzzzz_0_xxzzzz_1, g_0_xzzzzzz_0_xxzzzzz_0, g_0_xzzzzzz_0_xxzzzzz_1, g_0_xzzzzzz_0_xyyyyyz_0, g_0_xzzzzzz_0_xyyyyyz_1, g_0_xzzzzzz_0_xyyyyz_1, g_0_xzzzzzz_0_xyyyyzz_0, g_0_xzzzzzz_0_xyyyyzz_1, g_0_xzzzzzz_0_xyyyzz_1, g_0_xzzzzzz_0_xyyyzzz_0, g_0_xzzzzzz_0_xyyyzzz_1, g_0_xzzzzzz_0_xyyzzz_1, g_0_xzzzzzz_0_xyyzzzz_0, g_0_xzzzzzz_0_xyyzzzz_1, g_0_xzzzzzz_0_xyzzzz_1, g_0_xzzzzzz_0_xyzzzzz_0, g_0_xzzzzzz_0_xyzzzzz_1, g_0_xzzzzzz_0_xzzzzz_1, g_0_xzzzzzz_0_xzzzzzz_0, g_0_xzzzzzz_0_xzzzzzz_1, g_0_xzzzzzz_0_yyyyyyy_0, g_0_xzzzzzz_0_yyyyyyy_1, g_0_xzzzzzz_0_yyyyyyz_0, g_0_xzzzzzz_0_yyyyyyz_1, g_0_xzzzzzz_0_yyyyyz_1, g_0_xzzzzzz_0_yyyyyzz_0, g_0_xzzzzzz_0_yyyyyzz_1, g_0_xzzzzzz_0_yyyyzz_1, g_0_xzzzzzz_0_yyyyzzz_0, g_0_xzzzzzz_0_yyyyzzz_1, g_0_xzzzzzz_0_yyyzzz_1, g_0_xzzzzzz_0_yyyzzzz_0, g_0_xzzzzzz_0_yyyzzzz_1, g_0_xzzzzzz_0_yyzzzz_1, g_0_xzzzzzz_0_yyzzzzz_0, g_0_xzzzzzz_0_yyzzzzz_1, g_0_xzzzzzz_0_yzzzzz_1, g_0_xzzzzzz_0_yzzzzzz_0, g_0_xzzzzzz_0_yzzzzzz_1, g_0_xzzzzzz_0_zzzzzz_1, g_0_xzzzzzz_0_zzzzzzz_0, g_0_xzzzzzz_0_zzzzzzz_1, g_0_zzzzzz_0_xxxxxxz_0, g_0_zzzzzz_0_xxxxxxz_1, g_0_zzzzzz_0_xxxxxyz_0, g_0_zzzzzz_0_xxxxxyz_1, g_0_zzzzzz_0_xxxxxzz_0, g_0_zzzzzz_0_xxxxxzz_1, g_0_zzzzzz_0_xxxxyyz_0, g_0_zzzzzz_0_xxxxyyz_1, g_0_zzzzzz_0_xxxxyzz_0, g_0_zzzzzz_0_xxxxyzz_1, g_0_zzzzzz_0_xxxxzzz_0, g_0_zzzzzz_0_xxxxzzz_1, g_0_zzzzzz_0_xxxyyyz_0, g_0_zzzzzz_0_xxxyyyz_1, g_0_zzzzzz_0_xxxyyzz_0, g_0_zzzzzz_0_xxxyyzz_1, g_0_zzzzzz_0_xxxyzzz_0, g_0_zzzzzz_0_xxxyzzz_1, g_0_zzzzzz_0_xxxzzzz_0, g_0_zzzzzz_0_xxxzzzz_1, g_0_zzzzzz_0_xxyyyyz_0, g_0_zzzzzz_0_xxyyyyz_1, g_0_zzzzzz_0_xxyyyzz_0, g_0_zzzzzz_0_xxyyyzz_1, g_0_zzzzzz_0_xxyyzzz_0, g_0_zzzzzz_0_xxyyzzz_1, g_0_zzzzzz_0_xxyzzzz_0, g_0_zzzzzz_0_xxyzzzz_1, g_0_zzzzzz_0_xxzzzzz_0, g_0_zzzzzz_0_xxzzzzz_1, g_0_zzzzzz_0_xyyyyyz_0, g_0_zzzzzz_0_xyyyyyz_1, g_0_zzzzzz_0_xyyyyzz_0, g_0_zzzzzz_0_xyyyyzz_1, g_0_zzzzzz_0_xyyyzzz_0, g_0_zzzzzz_0_xyyyzzz_1, g_0_zzzzzz_0_xyyzzzz_0, g_0_zzzzzz_0_xyyzzzz_1, g_0_zzzzzz_0_xyzzzzz_0, g_0_zzzzzz_0_xyzzzzz_1, g_0_zzzzzz_0_xzzzzzz_0, g_0_zzzzzz_0_xzzzzzz_1, g_0_zzzzzz_0_yyyyyyy_0, g_0_zzzzzz_0_yyyyyyy_1, g_0_zzzzzz_0_yyyyyyz_0, g_0_zzzzzz_0_yyyyyyz_1, g_0_zzzzzz_0_yyyyyzz_0, g_0_zzzzzz_0_yyyyyzz_1, g_0_zzzzzz_0_yyyyzzz_0, g_0_zzzzzz_0_yyyyzzz_1, g_0_zzzzzz_0_yyyzzzz_0, g_0_zzzzzz_0_yyyzzzz_1, g_0_zzzzzz_0_yyzzzzz_0, g_0_zzzzzz_0_yyzzzzz_1, g_0_zzzzzz_0_yzzzzzz_0, g_0_zzzzzz_0_yzzzzzz_1, g_0_zzzzzz_0_zzzzzzz_0, g_0_zzzzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzzz_0_xxxxxxx_0[i] = 5.0 * g_0_xxzzzz_0_xxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxxxxx_0[i] * pb_z + g_0_xxzzzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxxxy_0[i] = 5.0 * g_0_xxzzzz_0_xxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxxxxy_0[i] * pb_z + g_0_xxzzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxxxz_0[i] = g_0_zzzzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xzzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxxxz_0[i] * pb_x + g_0_xzzzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_xxzzzz_0_xxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxxxyy_0[i] * pb_z + g_0_xxzzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxxyz_0[i] = g_0_zzzzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xzzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxxyz_0[i] * pb_x + g_0_xzzzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxxxzz_0[i] = g_0_zzzzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xzzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxxzz_0[i] * pb_x + g_0_xzzzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxxyyy_0[i] = 5.0 * g_0_xxzzzz_0_xxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxxyyy_0[i] * pb_z + g_0_xxzzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxyyz_0[i] = g_0_zzzzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxyyz_0[i] * pb_x + g_0_xzzzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxxyzz_0[i] = g_0_zzzzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxyzz_0[i] * pb_x + g_0_xzzzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxxzzz_0[i] = g_0_zzzzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxzzz_0[i] * pb_x + g_0_xzzzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxyyyy_0[i] = 5.0 * g_0_xxzzzz_0_xxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxyyyy_0[i] * pb_z + g_0_xxzzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxyyyz_0[i] = g_0_zzzzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxyyyz_0[i] * pb_x + g_0_xzzzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxyyzz_0[i] = g_0_zzzzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxyyzz_0[i] * pb_x + g_0_xzzzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxyzzz_0[i] = g_0_zzzzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxyzzz_0[i] * pb_x + g_0_xzzzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxzzzz_0[i] = g_0_zzzzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxzzzz_0_xxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxyyyyy_0[i] * pb_z + g_0_xxzzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxyyyyz_0[i] = g_0_zzzzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyyyyz_0[i] * pb_x + g_0_xzzzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyyyzz_0[i] = g_0_zzzzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyyyzz_0[i] * pb_x + g_0_xzzzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyyzzz_0[i] = g_0_zzzzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyyzzz_0[i] * pb_x + g_0_xzzzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyzzzz_0[i] = g_0_zzzzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxzzzzz_0[i] = g_0_zzzzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxzzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyyyyy_0[i] = 5.0 * g_0_xxzzzz_0_xyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xyyyyyy_0[i] * pb_z + g_0_xxzzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xyyyyyz_0[i] = g_0_zzzzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyyyyyz_0[i] * pb_x + g_0_xzzzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyyyzz_0[i] = g_0_zzzzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyyyyzz_0[i] * pb_x + g_0_xzzzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyyzzz_0[i] = g_0_zzzzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyyyzzz_0[i] * pb_x + g_0_xzzzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyzzzz_0[i] = g_0_zzzzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyyzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyzzzzz_0[i] = g_0_zzzzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyzzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xzzzzzz_0[i] = g_0_zzzzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xzzzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyyyy_0[i] = g_0_zzzzzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyyyy_0[i] * pb_x + g_0_xzzzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyyyz_0[i] = g_0_zzzzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyyyz_0[i] * pb_x + g_0_xzzzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyyzz_0[i] = g_0_zzzzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyyzz_0[i] * pb_x + g_0_xzzzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyzzz_0[i] = g_0_zzzzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyzzz_0[i] * pb_x + g_0_xzzzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyzzzz_0[i] = g_0_zzzzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyzzzz_0[i] * pb_x + g_0_xzzzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyzzzzz_0[i] = g_0_zzzzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyzzzzz_0[i] * pb_x + g_0_xzzzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yzzzzzz_0[i] = g_0_zzzzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzzzzzz_0[i] * pb_x + g_0_xzzzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_zzzzzzz_0[i] = g_0_zzzzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzzzzzz_0[i] * pb_x + g_0_xzzzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1008-1044 components of targeted buffer : SLSK

    auto g_0_xyyyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1008);

    auto g_0_xyyyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1009);

    auto g_0_xyyyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1010);

    auto g_0_xyyyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1011);

    auto g_0_xyyyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1012);

    auto g_0_xyyyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1013);

    auto g_0_xyyyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1014);

    auto g_0_xyyyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1015);

    auto g_0_xyyyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1016);

    auto g_0_xyyyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1017);

    auto g_0_xyyyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1018);

    auto g_0_xyyyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1019);

    auto g_0_xyyyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1020);

    auto g_0_xyyyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1021);

    auto g_0_xyyyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1022);

    auto g_0_xyyyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1023);

    auto g_0_xyyyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1024);

    auto g_0_xyyyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1025);

    auto g_0_xyyyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1026);

    auto g_0_xyyyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1027);

    auto g_0_xyyyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1028);

    auto g_0_xyyyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1029);

    auto g_0_xyyyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1030);

    auto g_0_xyyyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1031);

    auto g_0_xyyyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1032);

    auto g_0_xyyyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1033);

    auto g_0_xyyyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1034);

    auto g_0_xyyyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1035);

    auto g_0_xyyyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1036);

    auto g_0_xyyyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1037);

    auto g_0_xyyyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1038);

    auto g_0_xyyyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1039);

    auto g_0_xyyyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1040);

    auto g_0_xyyyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1041);

    auto g_0_xyyyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1042);

    auto g_0_xyyyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1043);

    #pragma omp simd aligned(g_0_xyyyyyyy_0_xxxxxxx_0, g_0_xyyyyyyy_0_xxxxxxy_0, g_0_xyyyyyyy_0_xxxxxxz_0, g_0_xyyyyyyy_0_xxxxxyy_0, g_0_xyyyyyyy_0_xxxxxyz_0, g_0_xyyyyyyy_0_xxxxxzz_0, g_0_xyyyyyyy_0_xxxxyyy_0, g_0_xyyyyyyy_0_xxxxyyz_0, g_0_xyyyyyyy_0_xxxxyzz_0, g_0_xyyyyyyy_0_xxxxzzz_0, g_0_xyyyyyyy_0_xxxyyyy_0, g_0_xyyyyyyy_0_xxxyyyz_0, g_0_xyyyyyyy_0_xxxyyzz_0, g_0_xyyyyyyy_0_xxxyzzz_0, g_0_xyyyyyyy_0_xxxzzzz_0, g_0_xyyyyyyy_0_xxyyyyy_0, g_0_xyyyyyyy_0_xxyyyyz_0, g_0_xyyyyyyy_0_xxyyyzz_0, g_0_xyyyyyyy_0_xxyyzzz_0, g_0_xyyyyyyy_0_xxyzzzz_0, g_0_xyyyyyyy_0_xxzzzzz_0, g_0_xyyyyyyy_0_xyyyyyy_0, g_0_xyyyyyyy_0_xyyyyyz_0, g_0_xyyyyyyy_0_xyyyyzz_0, g_0_xyyyyyyy_0_xyyyzzz_0, g_0_xyyyyyyy_0_xyyzzzz_0, g_0_xyyyyyyy_0_xyzzzzz_0, g_0_xyyyyyyy_0_xzzzzzz_0, g_0_xyyyyyyy_0_yyyyyyy_0, g_0_xyyyyyyy_0_yyyyyyz_0, g_0_xyyyyyyy_0_yyyyyzz_0, g_0_xyyyyyyy_0_yyyyzzz_0, g_0_xyyyyyyy_0_yyyzzzz_0, g_0_xyyyyyyy_0_yyzzzzz_0, g_0_xyyyyyyy_0_yzzzzzz_0, g_0_xyyyyyyy_0_zzzzzzz_0, g_0_yyyyyyy_0_xxxxxx_1, g_0_yyyyyyy_0_xxxxxxx_0, g_0_yyyyyyy_0_xxxxxxx_1, g_0_yyyyyyy_0_xxxxxxy_0, g_0_yyyyyyy_0_xxxxxxy_1, g_0_yyyyyyy_0_xxxxxxz_0, g_0_yyyyyyy_0_xxxxxxz_1, g_0_yyyyyyy_0_xxxxxy_1, g_0_yyyyyyy_0_xxxxxyy_0, g_0_yyyyyyy_0_xxxxxyy_1, g_0_yyyyyyy_0_xxxxxyz_0, g_0_yyyyyyy_0_xxxxxyz_1, g_0_yyyyyyy_0_xxxxxz_1, g_0_yyyyyyy_0_xxxxxzz_0, g_0_yyyyyyy_0_xxxxxzz_1, g_0_yyyyyyy_0_xxxxyy_1, g_0_yyyyyyy_0_xxxxyyy_0, g_0_yyyyyyy_0_xxxxyyy_1, g_0_yyyyyyy_0_xxxxyyz_0, g_0_yyyyyyy_0_xxxxyyz_1, g_0_yyyyyyy_0_xxxxyz_1, g_0_yyyyyyy_0_xxxxyzz_0, g_0_yyyyyyy_0_xxxxyzz_1, g_0_yyyyyyy_0_xxxxzz_1, g_0_yyyyyyy_0_xxxxzzz_0, g_0_yyyyyyy_0_xxxxzzz_1, g_0_yyyyyyy_0_xxxyyy_1, g_0_yyyyyyy_0_xxxyyyy_0, g_0_yyyyyyy_0_xxxyyyy_1, g_0_yyyyyyy_0_xxxyyyz_0, g_0_yyyyyyy_0_xxxyyyz_1, g_0_yyyyyyy_0_xxxyyz_1, g_0_yyyyyyy_0_xxxyyzz_0, g_0_yyyyyyy_0_xxxyyzz_1, g_0_yyyyyyy_0_xxxyzz_1, g_0_yyyyyyy_0_xxxyzzz_0, g_0_yyyyyyy_0_xxxyzzz_1, g_0_yyyyyyy_0_xxxzzz_1, g_0_yyyyyyy_0_xxxzzzz_0, g_0_yyyyyyy_0_xxxzzzz_1, g_0_yyyyyyy_0_xxyyyy_1, g_0_yyyyyyy_0_xxyyyyy_0, g_0_yyyyyyy_0_xxyyyyy_1, g_0_yyyyyyy_0_xxyyyyz_0, g_0_yyyyyyy_0_xxyyyyz_1, g_0_yyyyyyy_0_xxyyyz_1, g_0_yyyyyyy_0_xxyyyzz_0, g_0_yyyyyyy_0_xxyyyzz_1, g_0_yyyyyyy_0_xxyyzz_1, g_0_yyyyyyy_0_xxyyzzz_0, g_0_yyyyyyy_0_xxyyzzz_1, g_0_yyyyyyy_0_xxyzzz_1, g_0_yyyyyyy_0_xxyzzzz_0, g_0_yyyyyyy_0_xxyzzzz_1, g_0_yyyyyyy_0_xxzzzz_1, g_0_yyyyyyy_0_xxzzzzz_0, g_0_yyyyyyy_0_xxzzzzz_1, g_0_yyyyyyy_0_xyyyyy_1, g_0_yyyyyyy_0_xyyyyyy_0, g_0_yyyyyyy_0_xyyyyyy_1, g_0_yyyyyyy_0_xyyyyyz_0, g_0_yyyyyyy_0_xyyyyyz_1, g_0_yyyyyyy_0_xyyyyz_1, g_0_yyyyyyy_0_xyyyyzz_0, g_0_yyyyyyy_0_xyyyyzz_1, g_0_yyyyyyy_0_xyyyzz_1, g_0_yyyyyyy_0_xyyyzzz_0, g_0_yyyyyyy_0_xyyyzzz_1, g_0_yyyyyyy_0_xyyzzz_1, g_0_yyyyyyy_0_xyyzzzz_0, g_0_yyyyyyy_0_xyyzzzz_1, g_0_yyyyyyy_0_xyzzzz_1, g_0_yyyyyyy_0_xyzzzzz_0, g_0_yyyyyyy_0_xyzzzzz_1, g_0_yyyyyyy_0_xzzzzz_1, g_0_yyyyyyy_0_xzzzzzz_0, g_0_yyyyyyy_0_xzzzzzz_1, g_0_yyyyyyy_0_yyyyyy_1, g_0_yyyyyyy_0_yyyyyyy_0, g_0_yyyyyyy_0_yyyyyyy_1, g_0_yyyyyyy_0_yyyyyyz_0, g_0_yyyyyyy_0_yyyyyyz_1, g_0_yyyyyyy_0_yyyyyz_1, g_0_yyyyyyy_0_yyyyyzz_0, g_0_yyyyyyy_0_yyyyyzz_1, g_0_yyyyyyy_0_yyyyzz_1, g_0_yyyyyyy_0_yyyyzzz_0, g_0_yyyyyyy_0_yyyyzzz_1, g_0_yyyyyyy_0_yyyzzz_1, g_0_yyyyyyy_0_yyyzzzz_0, g_0_yyyyyyy_0_yyyzzzz_1, g_0_yyyyyyy_0_yyzzzz_1, g_0_yyyyyyy_0_yyzzzzz_0, g_0_yyyyyyy_0_yyzzzzz_1, g_0_yyyyyyy_0_yzzzzz_1, g_0_yyyyyyy_0_yzzzzzz_0, g_0_yyyyyyy_0_yzzzzzz_1, g_0_yyyyyyy_0_zzzzzz_1, g_0_yyyyyyy_0_zzzzzzz_0, g_0_yyyyyyy_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyy_0_xxxxxxx_0[i] = 7.0 * g_0_yyyyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxxx_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxxxy_0[i] = 6.0 * g_0_yyyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxxy_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxxz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxxyy_0[i] = 5.0 * g_0_yyyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxyy_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxyz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxyyy_0[i] = 4.0 * g_0_yyyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyyy_0[i] * pb_x + g_0_yyyyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyyz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyyyy_0[i] = 3.0 * g_0_yyyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyyy_0[i] * pb_x + g_0_yyyyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyyz_0[i] * pb_x + g_0_yyyyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyyyy_0[i] = 2.0 * g_0_yyyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyyy_0[i] * pb_x + g_0_yyyyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyyz_0[i] * pb_x + g_0_yyyyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyzz_0[i] * pb_x + g_0_yyyyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyyyy_0[i] = g_0_yyyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyyy_0[i] * pb_x + g_0_yyyyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyyyz_0[i] = g_0_yyyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyyz_0[i] * pb_x + g_0_yyyyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyyzz_0[i] = g_0_yyyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyzz_0[i] * pb_x + g_0_yyyyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyzzz_0[i] = g_0_yyyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyzzz_0[i] * pb_x + g_0_yyyyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyzzzz_0[i] = g_0_yyyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyzzzzz_0[i] = g_0_yyyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xzzzzzz_0[i] = g_0_yyyyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyyyy_0[i] = g_0_yyyyyyy_0_yyyyyyy_0[i] * pb_x + g_0_yyyyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyyyz_0[i] = g_0_yyyyyyy_0_yyyyyyz_0[i] * pb_x + g_0_yyyyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyyzz_0[i] = g_0_yyyyyyy_0_yyyyyzz_0[i] * pb_x + g_0_yyyyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyzzz_0[i] = g_0_yyyyyyy_0_yyyyzzz_0[i] * pb_x + g_0_yyyyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyzzzz_0[i] = g_0_yyyyyyy_0_yyyzzzz_0[i] * pb_x + g_0_yyyyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyzzzzz_0[i] = g_0_yyyyyyy_0_yyzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yzzzzzz_0[i] = g_0_yyyyyyy_0_yzzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_zzzzzzz_0[i] = g_0_yyyyyyy_0_zzzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1044-1080 components of targeted buffer : SLSK

    auto g_0_xyyyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1044);

    auto g_0_xyyyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1045);

    auto g_0_xyyyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1046);

    auto g_0_xyyyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1047);

    auto g_0_xyyyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1048);

    auto g_0_xyyyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1049);

    auto g_0_xyyyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1050);

    auto g_0_xyyyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1051);

    auto g_0_xyyyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1052);

    auto g_0_xyyyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1053);

    auto g_0_xyyyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1054);

    auto g_0_xyyyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1055);

    auto g_0_xyyyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1056);

    auto g_0_xyyyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1057);

    auto g_0_xyyyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1058);

    auto g_0_xyyyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1059);

    auto g_0_xyyyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1060);

    auto g_0_xyyyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1061);

    auto g_0_xyyyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1062);

    auto g_0_xyyyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1063);

    auto g_0_xyyyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1064);

    auto g_0_xyyyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1065);

    auto g_0_xyyyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1066);

    auto g_0_xyyyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1067);

    auto g_0_xyyyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1068);

    auto g_0_xyyyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1069);

    auto g_0_xyyyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1070);

    auto g_0_xyyyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1071);

    auto g_0_xyyyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1072);

    auto g_0_xyyyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1073);

    auto g_0_xyyyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1074);

    auto g_0_xyyyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1075);

    auto g_0_xyyyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1076);

    auto g_0_xyyyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1077);

    auto g_0_xyyyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1078);

    auto g_0_xyyyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1079);

    #pragma omp simd aligned(g_0_xyyyyyy_0_xxxxxxx_0, g_0_xyyyyyy_0_xxxxxxx_1, g_0_xyyyyyy_0_xxxxxxy_0, g_0_xyyyyyy_0_xxxxxxy_1, g_0_xyyyyyy_0_xxxxxyy_0, g_0_xyyyyyy_0_xxxxxyy_1, g_0_xyyyyyy_0_xxxxyyy_0, g_0_xyyyyyy_0_xxxxyyy_1, g_0_xyyyyyy_0_xxxyyyy_0, g_0_xyyyyyy_0_xxxyyyy_1, g_0_xyyyyyy_0_xxyyyyy_0, g_0_xyyyyyy_0_xxyyyyy_1, g_0_xyyyyyy_0_xyyyyyy_0, g_0_xyyyyyy_0_xyyyyyy_1, g_0_xyyyyyyz_0_xxxxxxx_0, g_0_xyyyyyyz_0_xxxxxxy_0, g_0_xyyyyyyz_0_xxxxxxz_0, g_0_xyyyyyyz_0_xxxxxyy_0, g_0_xyyyyyyz_0_xxxxxyz_0, g_0_xyyyyyyz_0_xxxxxzz_0, g_0_xyyyyyyz_0_xxxxyyy_0, g_0_xyyyyyyz_0_xxxxyyz_0, g_0_xyyyyyyz_0_xxxxyzz_0, g_0_xyyyyyyz_0_xxxxzzz_0, g_0_xyyyyyyz_0_xxxyyyy_0, g_0_xyyyyyyz_0_xxxyyyz_0, g_0_xyyyyyyz_0_xxxyyzz_0, g_0_xyyyyyyz_0_xxxyzzz_0, g_0_xyyyyyyz_0_xxxzzzz_0, g_0_xyyyyyyz_0_xxyyyyy_0, g_0_xyyyyyyz_0_xxyyyyz_0, g_0_xyyyyyyz_0_xxyyyzz_0, g_0_xyyyyyyz_0_xxyyzzz_0, g_0_xyyyyyyz_0_xxyzzzz_0, g_0_xyyyyyyz_0_xxzzzzz_0, g_0_xyyyyyyz_0_xyyyyyy_0, g_0_xyyyyyyz_0_xyyyyyz_0, g_0_xyyyyyyz_0_xyyyyzz_0, g_0_xyyyyyyz_0_xyyyzzz_0, g_0_xyyyyyyz_0_xyyzzzz_0, g_0_xyyyyyyz_0_xyzzzzz_0, g_0_xyyyyyyz_0_xzzzzzz_0, g_0_xyyyyyyz_0_yyyyyyy_0, g_0_xyyyyyyz_0_yyyyyyz_0, g_0_xyyyyyyz_0_yyyyyzz_0, g_0_xyyyyyyz_0_yyyyzzz_0, g_0_xyyyyyyz_0_yyyzzzz_0, g_0_xyyyyyyz_0_yyzzzzz_0, g_0_xyyyyyyz_0_yzzzzzz_0, g_0_xyyyyyyz_0_zzzzzzz_0, g_0_yyyyyyz_0_xxxxxxz_0, g_0_yyyyyyz_0_xxxxxxz_1, g_0_yyyyyyz_0_xxxxxyz_0, g_0_yyyyyyz_0_xxxxxyz_1, g_0_yyyyyyz_0_xxxxxz_1, g_0_yyyyyyz_0_xxxxxzz_0, g_0_yyyyyyz_0_xxxxxzz_1, g_0_yyyyyyz_0_xxxxyyz_0, g_0_yyyyyyz_0_xxxxyyz_1, g_0_yyyyyyz_0_xxxxyz_1, g_0_yyyyyyz_0_xxxxyzz_0, g_0_yyyyyyz_0_xxxxyzz_1, g_0_yyyyyyz_0_xxxxzz_1, g_0_yyyyyyz_0_xxxxzzz_0, g_0_yyyyyyz_0_xxxxzzz_1, g_0_yyyyyyz_0_xxxyyyz_0, g_0_yyyyyyz_0_xxxyyyz_1, g_0_yyyyyyz_0_xxxyyz_1, g_0_yyyyyyz_0_xxxyyzz_0, g_0_yyyyyyz_0_xxxyyzz_1, g_0_yyyyyyz_0_xxxyzz_1, g_0_yyyyyyz_0_xxxyzzz_0, g_0_yyyyyyz_0_xxxyzzz_1, g_0_yyyyyyz_0_xxxzzz_1, g_0_yyyyyyz_0_xxxzzzz_0, g_0_yyyyyyz_0_xxxzzzz_1, g_0_yyyyyyz_0_xxyyyyz_0, g_0_yyyyyyz_0_xxyyyyz_1, g_0_yyyyyyz_0_xxyyyz_1, g_0_yyyyyyz_0_xxyyyzz_0, g_0_yyyyyyz_0_xxyyyzz_1, g_0_yyyyyyz_0_xxyyzz_1, g_0_yyyyyyz_0_xxyyzzz_0, g_0_yyyyyyz_0_xxyyzzz_1, g_0_yyyyyyz_0_xxyzzz_1, g_0_yyyyyyz_0_xxyzzzz_0, g_0_yyyyyyz_0_xxyzzzz_1, g_0_yyyyyyz_0_xxzzzz_1, g_0_yyyyyyz_0_xxzzzzz_0, g_0_yyyyyyz_0_xxzzzzz_1, g_0_yyyyyyz_0_xyyyyyz_0, g_0_yyyyyyz_0_xyyyyyz_1, g_0_yyyyyyz_0_xyyyyz_1, g_0_yyyyyyz_0_xyyyyzz_0, g_0_yyyyyyz_0_xyyyyzz_1, g_0_yyyyyyz_0_xyyyzz_1, g_0_yyyyyyz_0_xyyyzzz_0, g_0_yyyyyyz_0_xyyyzzz_1, g_0_yyyyyyz_0_xyyzzz_1, g_0_yyyyyyz_0_xyyzzzz_0, g_0_yyyyyyz_0_xyyzzzz_1, g_0_yyyyyyz_0_xyzzzz_1, g_0_yyyyyyz_0_xyzzzzz_0, g_0_yyyyyyz_0_xyzzzzz_1, g_0_yyyyyyz_0_xzzzzz_1, g_0_yyyyyyz_0_xzzzzzz_0, g_0_yyyyyyz_0_xzzzzzz_1, g_0_yyyyyyz_0_yyyyyyy_0, g_0_yyyyyyz_0_yyyyyyy_1, g_0_yyyyyyz_0_yyyyyyz_0, g_0_yyyyyyz_0_yyyyyyz_1, g_0_yyyyyyz_0_yyyyyz_1, g_0_yyyyyyz_0_yyyyyzz_0, g_0_yyyyyyz_0_yyyyyzz_1, g_0_yyyyyyz_0_yyyyzz_1, g_0_yyyyyyz_0_yyyyzzz_0, g_0_yyyyyyz_0_yyyyzzz_1, g_0_yyyyyyz_0_yyyzzz_1, g_0_yyyyyyz_0_yyyzzzz_0, g_0_yyyyyyz_0_yyyzzzz_1, g_0_yyyyyyz_0_yyzzzz_1, g_0_yyyyyyz_0_yyzzzzz_0, g_0_yyyyyyz_0_yyzzzzz_1, g_0_yyyyyyz_0_yzzzzz_1, g_0_yyyyyyz_0_yzzzzzz_0, g_0_yyyyyyz_0_yzzzzzz_1, g_0_yyyyyyz_0_zzzzzz_1, g_0_yyyyyyz_0_zzzzzzz_0, g_0_yyyyyyz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyz_0_xxxxxxx_0[i] = g_0_xyyyyyy_0_xxxxxxx_0[i] * pb_z + g_0_xyyyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxxxy_0[i] = g_0_xyyyyyy_0_xxxxxxy_0[i] * pb_z + g_0_xyyyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyyyz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxxxz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxxxyy_0[i] = g_0_xyyyyyy_0_xxxxxyy_0[i] * pb_z + g_0_xyyyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyyyz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxxyz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyyyz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxxzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxxyyy_0[i] = g_0_xyyyyyy_0_xxxxyyy_0[i] * pb_z + g_0_xyyyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyyyz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxyyz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyyyz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxyzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyyyz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxyyyy_0[i] = g_0_xyyyyyy_0_xxxyyyy_0[i] * pb_z + g_0_xyyyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyyyz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxyyyz_0[i] * pb_x + g_0_yyyyyyz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyyyz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxyyzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyyyz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxyzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyyyz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyyyyy_0[i] = g_0_xyyyyyy_0_xxyyyyy_0[i] * pb_z + g_0_xyyyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyyyz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyyyyz_0[i] * pb_x + g_0_yyyyyyz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyyyz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyyyzz_0[i] * pb_x + g_0_yyyyyyz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyyyz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyyzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyyyz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyyyz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyyyyy_0[i] = g_0_xyyyyyy_0_xyyyyyy_0[i] * pb_z + g_0_xyyyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xyyyyyz_0[i] = g_0_yyyyyyz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyyyyz_0[i] * pb_x + g_0_yyyyyyz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyyyzz_0[i] = g_0_yyyyyyz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyyyzz_0[i] * pb_x + g_0_yyyyyyz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyyzzz_0[i] = g_0_yyyyyyz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyyzzz_0[i] * pb_x + g_0_yyyyyyz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyzzzz_0[i] = g_0_yyyyyyz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyzzzzz_0[i] = g_0_yyyyyyz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xzzzzzz_0[i] = g_0_yyyyyyz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xzzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyyyy_0[i] = g_0_yyyyyyz_0_yyyyyyy_0[i] * pb_x + g_0_yyyyyyz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyyyz_0[i] = g_0_yyyyyyz_0_yyyyyyz_0[i] * pb_x + g_0_yyyyyyz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyyzz_0[i] = g_0_yyyyyyz_0_yyyyyzz_0[i] * pb_x + g_0_yyyyyyz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyzzz_0[i] = g_0_yyyyyyz_0_yyyyzzz_0[i] * pb_x + g_0_yyyyyyz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyzzzz_0[i] = g_0_yyyyyyz_0_yyyzzzz_0[i] * pb_x + g_0_yyyyyyz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyzzzzz_0[i] = g_0_yyyyyyz_0_yyzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yzzzzzz_0[i] = g_0_yyyyyyz_0_yzzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_zzzzzzz_0[i] = g_0_yyyyyyz_0_zzzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1080-1116 components of targeted buffer : SLSK

    auto g_0_xyyyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1080);

    auto g_0_xyyyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1081);

    auto g_0_xyyyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1082);

    auto g_0_xyyyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1083);

    auto g_0_xyyyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1084);

    auto g_0_xyyyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1085);

    auto g_0_xyyyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1086);

    auto g_0_xyyyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1087);

    auto g_0_xyyyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1088);

    auto g_0_xyyyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1089);

    auto g_0_xyyyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1090);

    auto g_0_xyyyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1091);

    auto g_0_xyyyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1092);

    auto g_0_xyyyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1093);

    auto g_0_xyyyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1094);

    auto g_0_xyyyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1095);

    auto g_0_xyyyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1096);

    auto g_0_xyyyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1097);

    auto g_0_xyyyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1098);

    auto g_0_xyyyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1099);

    auto g_0_xyyyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1100);

    auto g_0_xyyyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1101);

    auto g_0_xyyyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1102);

    auto g_0_xyyyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1103);

    auto g_0_xyyyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1104);

    auto g_0_xyyyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1105);

    auto g_0_xyyyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1106);

    auto g_0_xyyyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1107);

    auto g_0_xyyyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1108);

    auto g_0_xyyyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1109);

    auto g_0_xyyyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1110);

    auto g_0_xyyyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1111);

    auto g_0_xyyyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1112);

    auto g_0_xyyyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1113);

    auto g_0_xyyyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1114);

    auto g_0_xyyyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1115);

    #pragma omp simd aligned(g_0_xyyyyyzz_0_xxxxxxx_0, g_0_xyyyyyzz_0_xxxxxxy_0, g_0_xyyyyyzz_0_xxxxxxz_0, g_0_xyyyyyzz_0_xxxxxyy_0, g_0_xyyyyyzz_0_xxxxxyz_0, g_0_xyyyyyzz_0_xxxxxzz_0, g_0_xyyyyyzz_0_xxxxyyy_0, g_0_xyyyyyzz_0_xxxxyyz_0, g_0_xyyyyyzz_0_xxxxyzz_0, g_0_xyyyyyzz_0_xxxxzzz_0, g_0_xyyyyyzz_0_xxxyyyy_0, g_0_xyyyyyzz_0_xxxyyyz_0, g_0_xyyyyyzz_0_xxxyyzz_0, g_0_xyyyyyzz_0_xxxyzzz_0, g_0_xyyyyyzz_0_xxxzzzz_0, g_0_xyyyyyzz_0_xxyyyyy_0, g_0_xyyyyyzz_0_xxyyyyz_0, g_0_xyyyyyzz_0_xxyyyzz_0, g_0_xyyyyyzz_0_xxyyzzz_0, g_0_xyyyyyzz_0_xxyzzzz_0, g_0_xyyyyyzz_0_xxzzzzz_0, g_0_xyyyyyzz_0_xyyyyyy_0, g_0_xyyyyyzz_0_xyyyyyz_0, g_0_xyyyyyzz_0_xyyyyzz_0, g_0_xyyyyyzz_0_xyyyzzz_0, g_0_xyyyyyzz_0_xyyzzzz_0, g_0_xyyyyyzz_0_xyzzzzz_0, g_0_xyyyyyzz_0_xzzzzzz_0, g_0_xyyyyyzz_0_yyyyyyy_0, g_0_xyyyyyzz_0_yyyyyyz_0, g_0_xyyyyyzz_0_yyyyyzz_0, g_0_xyyyyyzz_0_yyyyzzz_0, g_0_xyyyyyzz_0_yyyzzzz_0, g_0_xyyyyyzz_0_yyzzzzz_0, g_0_xyyyyyzz_0_yzzzzzz_0, g_0_xyyyyyzz_0_zzzzzzz_0, g_0_yyyyyzz_0_xxxxxx_1, g_0_yyyyyzz_0_xxxxxxx_0, g_0_yyyyyzz_0_xxxxxxx_1, g_0_yyyyyzz_0_xxxxxxy_0, g_0_yyyyyzz_0_xxxxxxy_1, g_0_yyyyyzz_0_xxxxxxz_0, g_0_yyyyyzz_0_xxxxxxz_1, g_0_yyyyyzz_0_xxxxxy_1, g_0_yyyyyzz_0_xxxxxyy_0, g_0_yyyyyzz_0_xxxxxyy_1, g_0_yyyyyzz_0_xxxxxyz_0, g_0_yyyyyzz_0_xxxxxyz_1, g_0_yyyyyzz_0_xxxxxz_1, g_0_yyyyyzz_0_xxxxxzz_0, g_0_yyyyyzz_0_xxxxxzz_1, g_0_yyyyyzz_0_xxxxyy_1, g_0_yyyyyzz_0_xxxxyyy_0, g_0_yyyyyzz_0_xxxxyyy_1, g_0_yyyyyzz_0_xxxxyyz_0, g_0_yyyyyzz_0_xxxxyyz_1, g_0_yyyyyzz_0_xxxxyz_1, g_0_yyyyyzz_0_xxxxyzz_0, g_0_yyyyyzz_0_xxxxyzz_1, g_0_yyyyyzz_0_xxxxzz_1, g_0_yyyyyzz_0_xxxxzzz_0, g_0_yyyyyzz_0_xxxxzzz_1, g_0_yyyyyzz_0_xxxyyy_1, g_0_yyyyyzz_0_xxxyyyy_0, g_0_yyyyyzz_0_xxxyyyy_1, g_0_yyyyyzz_0_xxxyyyz_0, g_0_yyyyyzz_0_xxxyyyz_1, g_0_yyyyyzz_0_xxxyyz_1, g_0_yyyyyzz_0_xxxyyzz_0, g_0_yyyyyzz_0_xxxyyzz_1, g_0_yyyyyzz_0_xxxyzz_1, g_0_yyyyyzz_0_xxxyzzz_0, g_0_yyyyyzz_0_xxxyzzz_1, g_0_yyyyyzz_0_xxxzzz_1, g_0_yyyyyzz_0_xxxzzzz_0, g_0_yyyyyzz_0_xxxzzzz_1, g_0_yyyyyzz_0_xxyyyy_1, g_0_yyyyyzz_0_xxyyyyy_0, g_0_yyyyyzz_0_xxyyyyy_1, g_0_yyyyyzz_0_xxyyyyz_0, g_0_yyyyyzz_0_xxyyyyz_1, g_0_yyyyyzz_0_xxyyyz_1, g_0_yyyyyzz_0_xxyyyzz_0, g_0_yyyyyzz_0_xxyyyzz_1, g_0_yyyyyzz_0_xxyyzz_1, g_0_yyyyyzz_0_xxyyzzz_0, g_0_yyyyyzz_0_xxyyzzz_1, g_0_yyyyyzz_0_xxyzzz_1, g_0_yyyyyzz_0_xxyzzzz_0, g_0_yyyyyzz_0_xxyzzzz_1, g_0_yyyyyzz_0_xxzzzz_1, g_0_yyyyyzz_0_xxzzzzz_0, g_0_yyyyyzz_0_xxzzzzz_1, g_0_yyyyyzz_0_xyyyyy_1, g_0_yyyyyzz_0_xyyyyyy_0, g_0_yyyyyzz_0_xyyyyyy_1, g_0_yyyyyzz_0_xyyyyyz_0, g_0_yyyyyzz_0_xyyyyyz_1, g_0_yyyyyzz_0_xyyyyz_1, g_0_yyyyyzz_0_xyyyyzz_0, g_0_yyyyyzz_0_xyyyyzz_1, g_0_yyyyyzz_0_xyyyzz_1, g_0_yyyyyzz_0_xyyyzzz_0, g_0_yyyyyzz_0_xyyyzzz_1, g_0_yyyyyzz_0_xyyzzz_1, g_0_yyyyyzz_0_xyyzzzz_0, g_0_yyyyyzz_0_xyyzzzz_1, g_0_yyyyyzz_0_xyzzzz_1, g_0_yyyyyzz_0_xyzzzzz_0, g_0_yyyyyzz_0_xyzzzzz_1, g_0_yyyyyzz_0_xzzzzz_1, g_0_yyyyyzz_0_xzzzzzz_0, g_0_yyyyyzz_0_xzzzzzz_1, g_0_yyyyyzz_0_yyyyyy_1, g_0_yyyyyzz_0_yyyyyyy_0, g_0_yyyyyzz_0_yyyyyyy_1, g_0_yyyyyzz_0_yyyyyyz_0, g_0_yyyyyzz_0_yyyyyyz_1, g_0_yyyyyzz_0_yyyyyz_1, g_0_yyyyyzz_0_yyyyyzz_0, g_0_yyyyyzz_0_yyyyyzz_1, g_0_yyyyyzz_0_yyyyzz_1, g_0_yyyyyzz_0_yyyyzzz_0, g_0_yyyyyzz_0_yyyyzzz_1, g_0_yyyyyzz_0_yyyzzz_1, g_0_yyyyyzz_0_yyyzzzz_0, g_0_yyyyyzz_0_yyyzzzz_1, g_0_yyyyyzz_0_yyzzzz_1, g_0_yyyyyzz_0_yyzzzzz_0, g_0_yyyyyzz_0_yyzzzzz_1, g_0_yyyyyzz_0_yzzzzz_1, g_0_yyyyyzz_0_yzzzzzz_0, g_0_yyyyyzz_0_yzzzzzz_1, g_0_yyyyyzz_0_zzzzzz_1, g_0_yyyyyzz_0_zzzzzzz_0, g_0_yyyyyzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyyyyzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxxx_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyyyyzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxxy_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyyzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxxz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyyyyzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxyy_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxyz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyyzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyyyyzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxyyy_0[i] * pb_x + g_0_yyyyyzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxyyz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxyzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyyzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyyyyzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyyyy_0[i] * pb_x + g_0_yyyyyzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyyyz_0[i] * pb_x + g_0_yyyyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyyzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyyzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyyyzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyyyy_0[i] * pb_x + g_0_yyyyyzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyyyz_0[i] * pb_x + g_0_yyyyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyyzz_0[i] * pb_x + g_0_yyyyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyyzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyyyy_0[i] = g_0_yyyyyzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyyyy_0[i] * pb_x + g_0_yyyyyzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyyyz_0[i] = g_0_yyyyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyyyz_0[i] * pb_x + g_0_yyyyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyyzz_0[i] = g_0_yyyyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyyzz_0[i] * pb_x + g_0_yyyyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyzzz_0[i] = g_0_yyyyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyzzz_0[i] * pb_x + g_0_yyyyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyzzzz_0[i] = g_0_yyyyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyzzzzz_0[i] = g_0_yyyyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xzzzzzz_0[i] = g_0_yyyyyzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xzzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyyyy_0[i] = g_0_yyyyyzz_0_yyyyyyy_0[i] * pb_x + g_0_yyyyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyyyz_0[i] = g_0_yyyyyzz_0_yyyyyyz_0[i] * pb_x + g_0_yyyyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyyzz_0[i] = g_0_yyyyyzz_0_yyyyyzz_0[i] * pb_x + g_0_yyyyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyzzz_0[i] = g_0_yyyyyzz_0_yyyyzzz_0[i] * pb_x + g_0_yyyyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyzzzz_0[i] = g_0_yyyyyzz_0_yyyzzzz_0[i] * pb_x + g_0_yyyyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyzzzzz_0[i] = g_0_yyyyyzz_0_yyzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yzzzzzz_0[i] = g_0_yyyyyzz_0_yzzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_zzzzzzz_0[i] = g_0_yyyyyzz_0_zzzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1116-1152 components of targeted buffer : SLSK

    auto g_0_xyyyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1116);

    auto g_0_xyyyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1117);

    auto g_0_xyyyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1118);

    auto g_0_xyyyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1119);

    auto g_0_xyyyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1120);

    auto g_0_xyyyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1121);

    auto g_0_xyyyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1122);

    auto g_0_xyyyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1123);

    auto g_0_xyyyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1124);

    auto g_0_xyyyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1125);

    auto g_0_xyyyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1126);

    auto g_0_xyyyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1127);

    auto g_0_xyyyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1128);

    auto g_0_xyyyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1129);

    auto g_0_xyyyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1130);

    auto g_0_xyyyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1131);

    auto g_0_xyyyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1132);

    auto g_0_xyyyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1133);

    auto g_0_xyyyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1134);

    auto g_0_xyyyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1135);

    auto g_0_xyyyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1136);

    auto g_0_xyyyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1137);

    auto g_0_xyyyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1138);

    auto g_0_xyyyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1139);

    auto g_0_xyyyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1140);

    auto g_0_xyyyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1141);

    auto g_0_xyyyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1142);

    auto g_0_xyyyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1143);

    auto g_0_xyyyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1144);

    auto g_0_xyyyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1145);

    auto g_0_xyyyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1146);

    auto g_0_xyyyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1147);

    auto g_0_xyyyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1148);

    auto g_0_xyyyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1149);

    auto g_0_xyyyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1150);

    auto g_0_xyyyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1151);

    #pragma omp simd aligned(g_0_xyyyyzzz_0_xxxxxxx_0, g_0_xyyyyzzz_0_xxxxxxy_0, g_0_xyyyyzzz_0_xxxxxxz_0, g_0_xyyyyzzz_0_xxxxxyy_0, g_0_xyyyyzzz_0_xxxxxyz_0, g_0_xyyyyzzz_0_xxxxxzz_0, g_0_xyyyyzzz_0_xxxxyyy_0, g_0_xyyyyzzz_0_xxxxyyz_0, g_0_xyyyyzzz_0_xxxxyzz_0, g_0_xyyyyzzz_0_xxxxzzz_0, g_0_xyyyyzzz_0_xxxyyyy_0, g_0_xyyyyzzz_0_xxxyyyz_0, g_0_xyyyyzzz_0_xxxyyzz_0, g_0_xyyyyzzz_0_xxxyzzz_0, g_0_xyyyyzzz_0_xxxzzzz_0, g_0_xyyyyzzz_0_xxyyyyy_0, g_0_xyyyyzzz_0_xxyyyyz_0, g_0_xyyyyzzz_0_xxyyyzz_0, g_0_xyyyyzzz_0_xxyyzzz_0, g_0_xyyyyzzz_0_xxyzzzz_0, g_0_xyyyyzzz_0_xxzzzzz_0, g_0_xyyyyzzz_0_xyyyyyy_0, g_0_xyyyyzzz_0_xyyyyyz_0, g_0_xyyyyzzz_0_xyyyyzz_0, g_0_xyyyyzzz_0_xyyyzzz_0, g_0_xyyyyzzz_0_xyyzzzz_0, g_0_xyyyyzzz_0_xyzzzzz_0, g_0_xyyyyzzz_0_xzzzzzz_0, g_0_xyyyyzzz_0_yyyyyyy_0, g_0_xyyyyzzz_0_yyyyyyz_0, g_0_xyyyyzzz_0_yyyyyzz_0, g_0_xyyyyzzz_0_yyyyzzz_0, g_0_xyyyyzzz_0_yyyzzzz_0, g_0_xyyyyzzz_0_yyzzzzz_0, g_0_xyyyyzzz_0_yzzzzzz_0, g_0_xyyyyzzz_0_zzzzzzz_0, g_0_yyyyzzz_0_xxxxxx_1, g_0_yyyyzzz_0_xxxxxxx_0, g_0_yyyyzzz_0_xxxxxxx_1, g_0_yyyyzzz_0_xxxxxxy_0, g_0_yyyyzzz_0_xxxxxxy_1, g_0_yyyyzzz_0_xxxxxxz_0, g_0_yyyyzzz_0_xxxxxxz_1, g_0_yyyyzzz_0_xxxxxy_1, g_0_yyyyzzz_0_xxxxxyy_0, g_0_yyyyzzz_0_xxxxxyy_1, g_0_yyyyzzz_0_xxxxxyz_0, g_0_yyyyzzz_0_xxxxxyz_1, g_0_yyyyzzz_0_xxxxxz_1, g_0_yyyyzzz_0_xxxxxzz_0, g_0_yyyyzzz_0_xxxxxzz_1, g_0_yyyyzzz_0_xxxxyy_1, g_0_yyyyzzz_0_xxxxyyy_0, g_0_yyyyzzz_0_xxxxyyy_1, g_0_yyyyzzz_0_xxxxyyz_0, g_0_yyyyzzz_0_xxxxyyz_1, g_0_yyyyzzz_0_xxxxyz_1, g_0_yyyyzzz_0_xxxxyzz_0, g_0_yyyyzzz_0_xxxxyzz_1, g_0_yyyyzzz_0_xxxxzz_1, g_0_yyyyzzz_0_xxxxzzz_0, g_0_yyyyzzz_0_xxxxzzz_1, g_0_yyyyzzz_0_xxxyyy_1, g_0_yyyyzzz_0_xxxyyyy_0, g_0_yyyyzzz_0_xxxyyyy_1, g_0_yyyyzzz_0_xxxyyyz_0, g_0_yyyyzzz_0_xxxyyyz_1, g_0_yyyyzzz_0_xxxyyz_1, g_0_yyyyzzz_0_xxxyyzz_0, g_0_yyyyzzz_0_xxxyyzz_1, g_0_yyyyzzz_0_xxxyzz_1, g_0_yyyyzzz_0_xxxyzzz_0, g_0_yyyyzzz_0_xxxyzzz_1, g_0_yyyyzzz_0_xxxzzz_1, g_0_yyyyzzz_0_xxxzzzz_0, g_0_yyyyzzz_0_xxxzzzz_1, g_0_yyyyzzz_0_xxyyyy_1, g_0_yyyyzzz_0_xxyyyyy_0, g_0_yyyyzzz_0_xxyyyyy_1, g_0_yyyyzzz_0_xxyyyyz_0, g_0_yyyyzzz_0_xxyyyyz_1, g_0_yyyyzzz_0_xxyyyz_1, g_0_yyyyzzz_0_xxyyyzz_0, g_0_yyyyzzz_0_xxyyyzz_1, g_0_yyyyzzz_0_xxyyzz_1, g_0_yyyyzzz_0_xxyyzzz_0, g_0_yyyyzzz_0_xxyyzzz_1, g_0_yyyyzzz_0_xxyzzz_1, g_0_yyyyzzz_0_xxyzzzz_0, g_0_yyyyzzz_0_xxyzzzz_1, g_0_yyyyzzz_0_xxzzzz_1, g_0_yyyyzzz_0_xxzzzzz_0, g_0_yyyyzzz_0_xxzzzzz_1, g_0_yyyyzzz_0_xyyyyy_1, g_0_yyyyzzz_0_xyyyyyy_0, g_0_yyyyzzz_0_xyyyyyy_1, g_0_yyyyzzz_0_xyyyyyz_0, g_0_yyyyzzz_0_xyyyyyz_1, g_0_yyyyzzz_0_xyyyyz_1, g_0_yyyyzzz_0_xyyyyzz_0, g_0_yyyyzzz_0_xyyyyzz_1, g_0_yyyyzzz_0_xyyyzz_1, g_0_yyyyzzz_0_xyyyzzz_0, g_0_yyyyzzz_0_xyyyzzz_1, g_0_yyyyzzz_0_xyyzzz_1, g_0_yyyyzzz_0_xyyzzzz_0, g_0_yyyyzzz_0_xyyzzzz_1, g_0_yyyyzzz_0_xyzzzz_1, g_0_yyyyzzz_0_xyzzzzz_0, g_0_yyyyzzz_0_xyzzzzz_1, g_0_yyyyzzz_0_xzzzzz_1, g_0_yyyyzzz_0_xzzzzzz_0, g_0_yyyyzzz_0_xzzzzzz_1, g_0_yyyyzzz_0_yyyyyy_1, g_0_yyyyzzz_0_yyyyyyy_0, g_0_yyyyzzz_0_yyyyyyy_1, g_0_yyyyzzz_0_yyyyyyz_0, g_0_yyyyzzz_0_yyyyyyz_1, g_0_yyyyzzz_0_yyyyyz_1, g_0_yyyyzzz_0_yyyyyzz_0, g_0_yyyyzzz_0_yyyyyzz_1, g_0_yyyyzzz_0_yyyyzz_1, g_0_yyyyzzz_0_yyyyzzz_0, g_0_yyyyzzz_0_yyyyzzz_1, g_0_yyyyzzz_0_yyyzzz_1, g_0_yyyyzzz_0_yyyzzzz_0, g_0_yyyyzzz_0_yyyzzzz_1, g_0_yyyyzzz_0_yyzzzz_1, g_0_yyyyzzz_0_yyzzzzz_0, g_0_yyyyzzz_0_yyzzzzz_1, g_0_yyyyzzz_0_yzzzzz_1, g_0_yyyyzzz_0_yzzzzzz_0, g_0_yyyyzzz_0_yzzzzzz_1, g_0_yyyyzzz_0_zzzzzz_1, g_0_yyyyzzz_0_zzzzzzz_0, g_0_yyyyzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyyyzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxxx_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyyyzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxxy_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxxz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyyyzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxyy_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxyz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyyyzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxyyy_0[i] * pb_x + g_0_yyyyzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxyyz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxyzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyyyzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyyyy_0[i] * pb_x + g_0_yyyyzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyyyz_0[i] * pb_x + g_0_yyyyzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyyzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyyzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyyyy_0[i] * pb_x + g_0_yyyyzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyyyz_0[i] * pb_x + g_0_yyyyzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyyzz_0[i] * pb_x + g_0_yyyyzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyyyy_0[i] = g_0_yyyyzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyyyy_0[i] * pb_x + g_0_yyyyzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyyyz_0[i] = g_0_yyyyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyyyz_0[i] * pb_x + g_0_yyyyzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyyzz_0[i] = g_0_yyyyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyyzz_0[i] * pb_x + g_0_yyyyzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyzzz_0[i] = g_0_yyyyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyzzz_0[i] * pb_x + g_0_yyyyzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyzzzz_0[i] = g_0_yyyyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyzzzzz_0[i] = g_0_yyyyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xzzzzzz_0[i] = g_0_yyyyzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xzzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyyyy_0[i] = g_0_yyyyzzz_0_yyyyyyy_0[i] * pb_x + g_0_yyyyzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyyyz_0[i] = g_0_yyyyzzz_0_yyyyyyz_0[i] * pb_x + g_0_yyyyzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyyzz_0[i] = g_0_yyyyzzz_0_yyyyyzz_0[i] * pb_x + g_0_yyyyzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyzzz_0[i] = g_0_yyyyzzz_0_yyyyzzz_0[i] * pb_x + g_0_yyyyzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyzzzz_0[i] = g_0_yyyyzzz_0_yyyzzzz_0[i] * pb_x + g_0_yyyyzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyzzzzz_0[i] = g_0_yyyyzzz_0_yyzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yzzzzzz_0[i] = g_0_yyyyzzz_0_yzzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_zzzzzzz_0[i] = g_0_yyyyzzz_0_zzzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1152-1188 components of targeted buffer : SLSK

    auto g_0_xyyyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1152);

    auto g_0_xyyyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1153);

    auto g_0_xyyyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1154);

    auto g_0_xyyyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1155);

    auto g_0_xyyyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1156);

    auto g_0_xyyyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1157);

    auto g_0_xyyyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1158);

    auto g_0_xyyyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1159);

    auto g_0_xyyyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1160);

    auto g_0_xyyyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1161);

    auto g_0_xyyyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1162);

    auto g_0_xyyyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1163);

    auto g_0_xyyyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1164);

    auto g_0_xyyyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1165);

    auto g_0_xyyyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1166);

    auto g_0_xyyyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1167);

    auto g_0_xyyyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1168);

    auto g_0_xyyyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1169);

    auto g_0_xyyyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1170);

    auto g_0_xyyyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1171);

    auto g_0_xyyyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1172);

    auto g_0_xyyyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1173);

    auto g_0_xyyyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1174);

    auto g_0_xyyyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1175);

    auto g_0_xyyyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1176);

    auto g_0_xyyyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1177);

    auto g_0_xyyyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1178);

    auto g_0_xyyyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1179);

    auto g_0_xyyyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1180);

    auto g_0_xyyyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1181);

    auto g_0_xyyyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1182);

    auto g_0_xyyyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1183);

    auto g_0_xyyyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1184);

    auto g_0_xyyyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1185);

    auto g_0_xyyyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1186);

    auto g_0_xyyyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1187);

    #pragma omp simd aligned(g_0_xyyyzzzz_0_xxxxxxx_0, g_0_xyyyzzzz_0_xxxxxxy_0, g_0_xyyyzzzz_0_xxxxxxz_0, g_0_xyyyzzzz_0_xxxxxyy_0, g_0_xyyyzzzz_0_xxxxxyz_0, g_0_xyyyzzzz_0_xxxxxzz_0, g_0_xyyyzzzz_0_xxxxyyy_0, g_0_xyyyzzzz_0_xxxxyyz_0, g_0_xyyyzzzz_0_xxxxyzz_0, g_0_xyyyzzzz_0_xxxxzzz_0, g_0_xyyyzzzz_0_xxxyyyy_0, g_0_xyyyzzzz_0_xxxyyyz_0, g_0_xyyyzzzz_0_xxxyyzz_0, g_0_xyyyzzzz_0_xxxyzzz_0, g_0_xyyyzzzz_0_xxxzzzz_0, g_0_xyyyzzzz_0_xxyyyyy_0, g_0_xyyyzzzz_0_xxyyyyz_0, g_0_xyyyzzzz_0_xxyyyzz_0, g_0_xyyyzzzz_0_xxyyzzz_0, g_0_xyyyzzzz_0_xxyzzzz_0, g_0_xyyyzzzz_0_xxzzzzz_0, g_0_xyyyzzzz_0_xyyyyyy_0, g_0_xyyyzzzz_0_xyyyyyz_0, g_0_xyyyzzzz_0_xyyyyzz_0, g_0_xyyyzzzz_0_xyyyzzz_0, g_0_xyyyzzzz_0_xyyzzzz_0, g_0_xyyyzzzz_0_xyzzzzz_0, g_0_xyyyzzzz_0_xzzzzzz_0, g_0_xyyyzzzz_0_yyyyyyy_0, g_0_xyyyzzzz_0_yyyyyyz_0, g_0_xyyyzzzz_0_yyyyyzz_0, g_0_xyyyzzzz_0_yyyyzzz_0, g_0_xyyyzzzz_0_yyyzzzz_0, g_0_xyyyzzzz_0_yyzzzzz_0, g_0_xyyyzzzz_0_yzzzzzz_0, g_0_xyyyzzzz_0_zzzzzzz_0, g_0_yyyzzzz_0_xxxxxx_1, g_0_yyyzzzz_0_xxxxxxx_0, g_0_yyyzzzz_0_xxxxxxx_1, g_0_yyyzzzz_0_xxxxxxy_0, g_0_yyyzzzz_0_xxxxxxy_1, g_0_yyyzzzz_0_xxxxxxz_0, g_0_yyyzzzz_0_xxxxxxz_1, g_0_yyyzzzz_0_xxxxxy_1, g_0_yyyzzzz_0_xxxxxyy_0, g_0_yyyzzzz_0_xxxxxyy_1, g_0_yyyzzzz_0_xxxxxyz_0, g_0_yyyzzzz_0_xxxxxyz_1, g_0_yyyzzzz_0_xxxxxz_1, g_0_yyyzzzz_0_xxxxxzz_0, g_0_yyyzzzz_0_xxxxxzz_1, g_0_yyyzzzz_0_xxxxyy_1, g_0_yyyzzzz_0_xxxxyyy_0, g_0_yyyzzzz_0_xxxxyyy_1, g_0_yyyzzzz_0_xxxxyyz_0, g_0_yyyzzzz_0_xxxxyyz_1, g_0_yyyzzzz_0_xxxxyz_1, g_0_yyyzzzz_0_xxxxyzz_0, g_0_yyyzzzz_0_xxxxyzz_1, g_0_yyyzzzz_0_xxxxzz_1, g_0_yyyzzzz_0_xxxxzzz_0, g_0_yyyzzzz_0_xxxxzzz_1, g_0_yyyzzzz_0_xxxyyy_1, g_0_yyyzzzz_0_xxxyyyy_0, g_0_yyyzzzz_0_xxxyyyy_1, g_0_yyyzzzz_0_xxxyyyz_0, g_0_yyyzzzz_0_xxxyyyz_1, g_0_yyyzzzz_0_xxxyyz_1, g_0_yyyzzzz_0_xxxyyzz_0, g_0_yyyzzzz_0_xxxyyzz_1, g_0_yyyzzzz_0_xxxyzz_1, g_0_yyyzzzz_0_xxxyzzz_0, g_0_yyyzzzz_0_xxxyzzz_1, g_0_yyyzzzz_0_xxxzzz_1, g_0_yyyzzzz_0_xxxzzzz_0, g_0_yyyzzzz_0_xxxzzzz_1, g_0_yyyzzzz_0_xxyyyy_1, g_0_yyyzzzz_0_xxyyyyy_0, g_0_yyyzzzz_0_xxyyyyy_1, g_0_yyyzzzz_0_xxyyyyz_0, g_0_yyyzzzz_0_xxyyyyz_1, g_0_yyyzzzz_0_xxyyyz_1, g_0_yyyzzzz_0_xxyyyzz_0, g_0_yyyzzzz_0_xxyyyzz_1, g_0_yyyzzzz_0_xxyyzz_1, g_0_yyyzzzz_0_xxyyzzz_0, g_0_yyyzzzz_0_xxyyzzz_1, g_0_yyyzzzz_0_xxyzzz_1, g_0_yyyzzzz_0_xxyzzzz_0, g_0_yyyzzzz_0_xxyzzzz_1, g_0_yyyzzzz_0_xxzzzz_1, g_0_yyyzzzz_0_xxzzzzz_0, g_0_yyyzzzz_0_xxzzzzz_1, g_0_yyyzzzz_0_xyyyyy_1, g_0_yyyzzzz_0_xyyyyyy_0, g_0_yyyzzzz_0_xyyyyyy_1, g_0_yyyzzzz_0_xyyyyyz_0, g_0_yyyzzzz_0_xyyyyyz_1, g_0_yyyzzzz_0_xyyyyz_1, g_0_yyyzzzz_0_xyyyyzz_0, g_0_yyyzzzz_0_xyyyyzz_1, g_0_yyyzzzz_0_xyyyzz_1, g_0_yyyzzzz_0_xyyyzzz_0, g_0_yyyzzzz_0_xyyyzzz_1, g_0_yyyzzzz_0_xyyzzz_1, g_0_yyyzzzz_0_xyyzzzz_0, g_0_yyyzzzz_0_xyyzzzz_1, g_0_yyyzzzz_0_xyzzzz_1, g_0_yyyzzzz_0_xyzzzzz_0, g_0_yyyzzzz_0_xyzzzzz_1, g_0_yyyzzzz_0_xzzzzz_1, g_0_yyyzzzz_0_xzzzzzz_0, g_0_yyyzzzz_0_xzzzzzz_1, g_0_yyyzzzz_0_yyyyyy_1, g_0_yyyzzzz_0_yyyyyyy_0, g_0_yyyzzzz_0_yyyyyyy_1, g_0_yyyzzzz_0_yyyyyyz_0, g_0_yyyzzzz_0_yyyyyyz_1, g_0_yyyzzzz_0_yyyyyz_1, g_0_yyyzzzz_0_yyyyyzz_0, g_0_yyyzzzz_0_yyyyyzz_1, g_0_yyyzzzz_0_yyyyzz_1, g_0_yyyzzzz_0_yyyyzzz_0, g_0_yyyzzzz_0_yyyyzzz_1, g_0_yyyzzzz_0_yyyzzz_1, g_0_yyyzzzz_0_yyyzzzz_0, g_0_yyyzzzz_0_yyyzzzz_1, g_0_yyyzzzz_0_yyzzzz_1, g_0_yyyzzzz_0_yyzzzzz_0, g_0_yyyzzzz_0_yyzzzzz_1, g_0_yyyzzzz_0_yzzzzz_1, g_0_yyyzzzz_0_yzzzzzz_0, g_0_yyyzzzz_0_yzzzzzz_1, g_0_yyyzzzz_0_zzzzzz_1, g_0_yyyzzzz_0_zzzzzzz_0, g_0_yyyzzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyyzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxxx_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyyzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxxy_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxxz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyyzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxyy_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxyz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyyzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxyyy_0[i] * pb_x + g_0_yyyzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxyyz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxyzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyyzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyyyy_0[i] * pb_x + g_0_yyyzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyyyz_0[i] * pb_x + g_0_yyyzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyyzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyyyy_0[i] * pb_x + g_0_yyyzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyyyz_0[i] * pb_x + g_0_yyyzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyyzz_0[i] * pb_x + g_0_yyyzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyyyy_0[i] = g_0_yyyzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyyyy_0[i] * pb_x + g_0_yyyzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyyyz_0[i] = g_0_yyyzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyyyz_0[i] * pb_x + g_0_yyyzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyyzz_0[i] = g_0_yyyzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyyzz_0[i] * pb_x + g_0_yyyzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyzzz_0[i] = g_0_yyyzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyzzz_0[i] * pb_x + g_0_yyyzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyzzzz_0[i] = g_0_yyyzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyzzzzz_0[i] = g_0_yyyzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xzzzzzz_0[i] = g_0_yyyzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xzzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyyyy_0[i] = g_0_yyyzzzz_0_yyyyyyy_0[i] * pb_x + g_0_yyyzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyyyz_0[i] = g_0_yyyzzzz_0_yyyyyyz_0[i] * pb_x + g_0_yyyzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyyzz_0[i] = g_0_yyyzzzz_0_yyyyyzz_0[i] * pb_x + g_0_yyyzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyzzz_0[i] = g_0_yyyzzzz_0_yyyyzzz_0[i] * pb_x + g_0_yyyzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyzzzz_0[i] = g_0_yyyzzzz_0_yyyzzzz_0[i] * pb_x + g_0_yyyzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyzzzzz_0[i] = g_0_yyyzzzz_0_yyzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yzzzzzz_0[i] = g_0_yyyzzzz_0_yzzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_zzzzzzz_0[i] = g_0_yyyzzzz_0_zzzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1188-1224 components of targeted buffer : SLSK

    auto g_0_xyyzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1188);

    auto g_0_xyyzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1189);

    auto g_0_xyyzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1190);

    auto g_0_xyyzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1191);

    auto g_0_xyyzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1192);

    auto g_0_xyyzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1193);

    auto g_0_xyyzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1194);

    auto g_0_xyyzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1195);

    auto g_0_xyyzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1196);

    auto g_0_xyyzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1197);

    auto g_0_xyyzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1198);

    auto g_0_xyyzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1199);

    auto g_0_xyyzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1200);

    auto g_0_xyyzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1201);

    auto g_0_xyyzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1202);

    auto g_0_xyyzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1203);

    auto g_0_xyyzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1204);

    auto g_0_xyyzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1205);

    auto g_0_xyyzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1206);

    auto g_0_xyyzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1207);

    auto g_0_xyyzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1208);

    auto g_0_xyyzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1209);

    auto g_0_xyyzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1210);

    auto g_0_xyyzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1211);

    auto g_0_xyyzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1212);

    auto g_0_xyyzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1213);

    auto g_0_xyyzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1214);

    auto g_0_xyyzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1215);

    auto g_0_xyyzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1216);

    auto g_0_xyyzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1217);

    auto g_0_xyyzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1218);

    auto g_0_xyyzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1219);

    auto g_0_xyyzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1220);

    auto g_0_xyyzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1221);

    auto g_0_xyyzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1222);

    auto g_0_xyyzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1223);

    #pragma omp simd aligned(g_0_xyyzzzzz_0_xxxxxxx_0, g_0_xyyzzzzz_0_xxxxxxy_0, g_0_xyyzzzzz_0_xxxxxxz_0, g_0_xyyzzzzz_0_xxxxxyy_0, g_0_xyyzzzzz_0_xxxxxyz_0, g_0_xyyzzzzz_0_xxxxxzz_0, g_0_xyyzzzzz_0_xxxxyyy_0, g_0_xyyzzzzz_0_xxxxyyz_0, g_0_xyyzzzzz_0_xxxxyzz_0, g_0_xyyzzzzz_0_xxxxzzz_0, g_0_xyyzzzzz_0_xxxyyyy_0, g_0_xyyzzzzz_0_xxxyyyz_0, g_0_xyyzzzzz_0_xxxyyzz_0, g_0_xyyzzzzz_0_xxxyzzz_0, g_0_xyyzzzzz_0_xxxzzzz_0, g_0_xyyzzzzz_0_xxyyyyy_0, g_0_xyyzzzzz_0_xxyyyyz_0, g_0_xyyzzzzz_0_xxyyyzz_0, g_0_xyyzzzzz_0_xxyyzzz_0, g_0_xyyzzzzz_0_xxyzzzz_0, g_0_xyyzzzzz_0_xxzzzzz_0, g_0_xyyzzzzz_0_xyyyyyy_0, g_0_xyyzzzzz_0_xyyyyyz_0, g_0_xyyzzzzz_0_xyyyyzz_0, g_0_xyyzzzzz_0_xyyyzzz_0, g_0_xyyzzzzz_0_xyyzzzz_0, g_0_xyyzzzzz_0_xyzzzzz_0, g_0_xyyzzzzz_0_xzzzzzz_0, g_0_xyyzzzzz_0_yyyyyyy_0, g_0_xyyzzzzz_0_yyyyyyz_0, g_0_xyyzzzzz_0_yyyyyzz_0, g_0_xyyzzzzz_0_yyyyzzz_0, g_0_xyyzzzzz_0_yyyzzzz_0, g_0_xyyzzzzz_0_yyzzzzz_0, g_0_xyyzzzzz_0_yzzzzzz_0, g_0_xyyzzzzz_0_zzzzzzz_0, g_0_yyzzzzz_0_xxxxxx_1, g_0_yyzzzzz_0_xxxxxxx_0, g_0_yyzzzzz_0_xxxxxxx_1, g_0_yyzzzzz_0_xxxxxxy_0, g_0_yyzzzzz_0_xxxxxxy_1, g_0_yyzzzzz_0_xxxxxxz_0, g_0_yyzzzzz_0_xxxxxxz_1, g_0_yyzzzzz_0_xxxxxy_1, g_0_yyzzzzz_0_xxxxxyy_0, g_0_yyzzzzz_0_xxxxxyy_1, g_0_yyzzzzz_0_xxxxxyz_0, g_0_yyzzzzz_0_xxxxxyz_1, g_0_yyzzzzz_0_xxxxxz_1, g_0_yyzzzzz_0_xxxxxzz_0, g_0_yyzzzzz_0_xxxxxzz_1, g_0_yyzzzzz_0_xxxxyy_1, g_0_yyzzzzz_0_xxxxyyy_0, g_0_yyzzzzz_0_xxxxyyy_1, g_0_yyzzzzz_0_xxxxyyz_0, g_0_yyzzzzz_0_xxxxyyz_1, g_0_yyzzzzz_0_xxxxyz_1, g_0_yyzzzzz_0_xxxxyzz_0, g_0_yyzzzzz_0_xxxxyzz_1, g_0_yyzzzzz_0_xxxxzz_1, g_0_yyzzzzz_0_xxxxzzz_0, g_0_yyzzzzz_0_xxxxzzz_1, g_0_yyzzzzz_0_xxxyyy_1, g_0_yyzzzzz_0_xxxyyyy_0, g_0_yyzzzzz_0_xxxyyyy_1, g_0_yyzzzzz_0_xxxyyyz_0, g_0_yyzzzzz_0_xxxyyyz_1, g_0_yyzzzzz_0_xxxyyz_1, g_0_yyzzzzz_0_xxxyyzz_0, g_0_yyzzzzz_0_xxxyyzz_1, g_0_yyzzzzz_0_xxxyzz_1, g_0_yyzzzzz_0_xxxyzzz_0, g_0_yyzzzzz_0_xxxyzzz_1, g_0_yyzzzzz_0_xxxzzz_1, g_0_yyzzzzz_0_xxxzzzz_0, g_0_yyzzzzz_0_xxxzzzz_1, g_0_yyzzzzz_0_xxyyyy_1, g_0_yyzzzzz_0_xxyyyyy_0, g_0_yyzzzzz_0_xxyyyyy_1, g_0_yyzzzzz_0_xxyyyyz_0, g_0_yyzzzzz_0_xxyyyyz_1, g_0_yyzzzzz_0_xxyyyz_1, g_0_yyzzzzz_0_xxyyyzz_0, g_0_yyzzzzz_0_xxyyyzz_1, g_0_yyzzzzz_0_xxyyzz_1, g_0_yyzzzzz_0_xxyyzzz_0, g_0_yyzzzzz_0_xxyyzzz_1, g_0_yyzzzzz_0_xxyzzz_1, g_0_yyzzzzz_0_xxyzzzz_0, g_0_yyzzzzz_0_xxyzzzz_1, g_0_yyzzzzz_0_xxzzzz_1, g_0_yyzzzzz_0_xxzzzzz_0, g_0_yyzzzzz_0_xxzzzzz_1, g_0_yyzzzzz_0_xyyyyy_1, g_0_yyzzzzz_0_xyyyyyy_0, g_0_yyzzzzz_0_xyyyyyy_1, g_0_yyzzzzz_0_xyyyyyz_0, g_0_yyzzzzz_0_xyyyyyz_1, g_0_yyzzzzz_0_xyyyyz_1, g_0_yyzzzzz_0_xyyyyzz_0, g_0_yyzzzzz_0_xyyyyzz_1, g_0_yyzzzzz_0_xyyyzz_1, g_0_yyzzzzz_0_xyyyzzz_0, g_0_yyzzzzz_0_xyyyzzz_1, g_0_yyzzzzz_0_xyyzzz_1, g_0_yyzzzzz_0_xyyzzzz_0, g_0_yyzzzzz_0_xyyzzzz_1, g_0_yyzzzzz_0_xyzzzz_1, g_0_yyzzzzz_0_xyzzzzz_0, g_0_yyzzzzz_0_xyzzzzz_1, g_0_yyzzzzz_0_xzzzzz_1, g_0_yyzzzzz_0_xzzzzzz_0, g_0_yyzzzzz_0_xzzzzzz_1, g_0_yyzzzzz_0_yyyyyy_1, g_0_yyzzzzz_0_yyyyyyy_0, g_0_yyzzzzz_0_yyyyyyy_1, g_0_yyzzzzz_0_yyyyyyz_0, g_0_yyzzzzz_0_yyyyyyz_1, g_0_yyzzzzz_0_yyyyyz_1, g_0_yyzzzzz_0_yyyyyzz_0, g_0_yyzzzzz_0_yyyyyzz_1, g_0_yyzzzzz_0_yyyyzz_1, g_0_yyzzzzz_0_yyyyzzz_0, g_0_yyzzzzz_0_yyyyzzz_1, g_0_yyzzzzz_0_yyyzzz_1, g_0_yyzzzzz_0_yyyzzzz_0, g_0_yyzzzzz_0_yyyzzzz_1, g_0_yyzzzzz_0_yyzzzz_1, g_0_yyzzzzz_0_yyzzzzz_0, g_0_yyzzzzz_0_yyzzzzz_1, g_0_yyzzzzz_0_yzzzzz_1, g_0_yyzzzzz_0_yzzzzzz_0, g_0_yyzzzzz_0_yzzzzzz_1, g_0_yyzzzzz_0_zzzzzz_1, g_0_yyzzzzz_0_zzzzzzz_0, g_0_yyzzzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyzzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxxx_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxxy_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxxz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxyy_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxyz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxyyy_0[i] * pb_x + g_0_yyzzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxyyz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxyzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyyyy_0[i] * pb_x + g_0_yyzzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyyyz_0[i] * pb_x + g_0_yyzzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyyzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyyyy_0[i] * pb_x + g_0_yyzzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyyyz_0[i] * pb_x + g_0_yyzzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyyzz_0[i] * pb_x + g_0_yyzzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyyyy_0[i] = g_0_yyzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyyyy_0[i] * pb_x + g_0_yyzzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyyyz_0[i] = g_0_yyzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyyyz_0[i] * pb_x + g_0_yyzzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyyzz_0[i] = g_0_yyzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyyzz_0[i] * pb_x + g_0_yyzzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyzzz_0[i] = g_0_yyzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyzzz_0[i] * pb_x + g_0_yyzzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyzzzz_0[i] = g_0_yyzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyzzzzz_0[i] = g_0_yyzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xzzzzzz_0[i] = g_0_yyzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xzzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyyyy_0[i] = g_0_yyzzzzz_0_yyyyyyy_0[i] * pb_x + g_0_yyzzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyyyz_0[i] = g_0_yyzzzzz_0_yyyyyyz_0[i] * pb_x + g_0_yyzzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyyzz_0[i] = g_0_yyzzzzz_0_yyyyyzz_0[i] * pb_x + g_0_yyzzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyzzz_0[i] = g_0_yyzzzzz_0_yyyyzzz_0[i] * pb_x + g_0_yyzzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyzzzz_0[i] = g_0_yyzzzzz_0_yyyzzzz_0[i] * pb_x + g_0_yyzzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyzzzzz_0[i] = g_0_yyzzzzz_0_yyzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yzzzzzz_0[i] = g_0_yyzzzzz_0_yzzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_zzzzzzz_0[i] = g_0_yyzzzzz_0_zzzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1224-1260 components of targeted buffer : SLSK

    auto g_0_xyzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1224);

    auto g_0_xyzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1225);

    auto g_0_xyzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1226);

    auto g_0_xyzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1227);

    auto g_0_xyzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1228);

    auto g_0_xyzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1229);

    auto g_0_xyzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1230);

    auto g_0_xyzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1231);

    auto g_0_xyzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1232);

    auto g_0_xyzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1233);

    auto g_0_xyzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1234);

    auto g_0_xyzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1235);

    auto g_0_xyzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1236);

    auto g_0_xyzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1237);

    auto g_0_xyzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1238);

    auto g_0_xyzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1239);

    auto g_0_xyzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1240);

    auto g_0_xyzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1241);

    auto g_0_xyzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1242);

    auto g_0_xyzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1243);

    auto g_0_xyzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1244);

    auto g_0_xyzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1245);

    auto g_0_xyzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1246);

    auto g_0_xyzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1247);

    auto g_0_xyzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1248);

    auto g_0_xyzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1249);

    auto g_0_xyzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1250);

    auto g_0_xyzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1251);

    auto g_0_xyzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1252);

    auto g_0_xyzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1253);

    auto g_0_xyzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1254);

    auto g_0_xyzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1255);

    auto g_0_xyzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1256);

    auto g_0_xyzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1257);

    auto g_0_xyzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1258);

    auto g_0_xyzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1259);

    #pragma omp simd aligned(g_0_xyzzzzzz_0_xxxxxxx_0, g_0_xyzzzzzz_0_xxxxxxy_0, g_0_xyzzzzzz_0_xxxxxxz_0, g_0_xyzzzzzz_0_xxxxxyy_0, g_0_xyzzzzzz_0_xxxxxyz_0, g_0_xyzzzzzz_0_xxxxxzz_0, g_0_xyzzzzzz_0_xxxxyyy_0, g_0_xyzzzzzz_0_xxxxyyz_0, g_0_xyzzzzzz_0_xxxxyzz_0, g_0_xyzzzzzz_0_xxxxzzz_0, g_0_xyzzzzzz_0_xxxyyyy_0, g_0_xyzzzzzz_0_xxxyyyz_0, g_0_xyzzzzzz_0_xxxyyzz_0, g_0_xyzzzzzz_0_xxxyzzz_0, g_0_xyzzzzzz_0_xxxzzzz_0, g_0_xyzzzzzz_0_xxyyyyy_0, g_0_xyzzzzzz_0_xxyyyyz_0, g_0_xyzzzzzz_0_xxyyyzz_0, g_0_xyzzzzzz_0_xxyyzzz_0, g_0_xyzzzzzz_0_xxyzzzz_0, g_0_xyzzzzzz_0_xxzzzzz_0, g_0_xyzzzzzz_0_xyyyyyy_0, g_0_xyzzzzzz_0_xyyyyyz_0, g_0_xyzzzzzz_0_xyyyyzz_0, g_0_xyzzzzzz_0_xyyyzzz_0, g_0_xyzzzzzz_0_xyyzzzz_0, g_0_xyzzzzzz_0_xyzzzzz_0, g_0_xyzzzzzz_0_xzzzzzz_0, g_0_xyzzzzzz_0_yyyyyyy_0, g_0_xyzzzzzz_0_yyyyyyz_0, g_0_xyzzzzzz_0_yyyyyzz_0, g_0_xyzzzzzz_0_yyyyzzz_0, g_0_xyzzzzzz_0_yyyzzzz_0, g_0_xyzzzzzz_0_yyzzzzz_0, g_0_xyzzzzzz_0_yzzzzzz_0, g_0_xyzzzzzz_0_zzzzzzz_0, g_0_xzzzzzz_0_xxxxxxx_0, g_0_xzzzzzz_0_xxxxxxx_1, g_0_xzzzzzz_0_xxxxxxz_0, g_0_xzzzzzz_0_xxxxxxz_1, g_0_xzzzzzz_0_xxxxxzz_0, g_0_xzzzzzz_0_xxxxxzz_1, g_0_xzzzzzz_0_xxxxzzz_0, g_0_xzzzzzz_0_xxxxzzz_1, g_0_xzzzzzz_0_xxxzzzz_0, g_0_xzzzzzz_0_xxxzzzz_1, g_0_xzzzzzz_0_xxzzzzz_0, g_0_xzzzzzz_0_xxzzzzz_1, g_0_xzzzzzz_0_xzzzzzz_0, g_0_xzzzzzz_0_xzzzzzz_1, g_0_yzzzzzz_0_xxxxxxy_0, g_0_yzzzzzz_0_xxxxxxy_1, g_0_yzzzzzz_0_xxxxxy_1, g_0_yzzzzzz_0_xxxxxyy_0, g_0_yzzzzzz_0_xxxxxyy_1, g_0_yzzzzzz_0_xxxxxyz_0, g_0_yzzzzzz_0_xxxxxyz_1, g_0_yzzzzzz_0_xxxxyy_1, g_0_yzzzzzz_0_xxxxyyy_0, g_0_yzzzzzz_0_xxxxyyy_1, g_0_yzzzzzz_0_xxxxyyz_0, g_0_yzzzzzz_0_xxxxyyz_1, g_0_yzzzzzz_0_xxxxyz_1, g_0_yzzzzzz_0_xxxxyzz_0, g_0_yzzzzzz_0_xxxxyzz_1, g_0_yzzzzzz_0_xxxyyy_1, g_0_yzzzzzz_0_xxxyyyy_0, g_0_yzzzzzz_0_xxxyyyy_1, g_0_yzzzzzz_0_xxxyyyz_0, g_0_yzzzzzz_0_xxxyyyz_1, g_0_yzzzzzz_0_xxxyyz_1, g_0_yzzzzzz_0_xxxyyzz_0, g_0_yzzzzzz_0_xxxyyzz_1, g_0_yzzzzzz_0_xxxyzz_1, g_0_yzzzzzz_0_xxxyzzz_0, g_0_yzzzzzz_0_xxxyzzz_1, g_0_yzzzzzz_0_xxyyyy_1, g_0_yzzzzzz_0_xxyyyyy_0, g_0_yzzzzzz_0_xxyyyyy_1, g_0_yzzzzzz_0_xxyyyyz_0, g_0_yzzzzzz_0_xxyyyyz_1, g_0_yzzzzzz_0_xxyyyz_1, g_0_yzzzzzz_0_xxyyyzz_0, g_0_yzzzzzz_0_xxyyyzz_1, g_0_yzzzzzz_0_xxyyzz_1, g_0_yzzzzzz_0_xxyyzzz_0, g_0_yzzzzzz_0_xxyyzzz_1, g_0_yzzzzzz_0_xxyzzz_1, g_0_yzzzzzz_0_xxyzzzz_0, g_0_yzzzzzz_0_xxyzzzz_1, g_0_yzzzzzz_0_xyyyyy_1, g_0_yzzzzzz_0_xyyyyyy_0, g_0_yzzzzzz_0_xyyyyyy_1, g_0_yzzzzzz_0_xyyyyyz_0, g_0_yzzzzzz_0_xyyyyyz_1, g_0_yzzzzzz_0_xyyyyz_1, g_0_yzzzzzz_0_xyyyyzz_0, g_0_yzzzzzz_0_xyyyyzz_1, g_0_yzzzzzz_0_xyyyzz_1, g_0_yzzzzzz_0_xyyyzzz_0, g_0_yzzzzzz_0_xyyyzzz_1, g_0_yzzzzzz_0_xyyzzz_1, g_0_yzzzzzz_0_xyyzzzz_0, g_0_yzzzzzz_0_xyyzzzz_1, g_0_yzzzzzz_0_xyzzzz_1, g_0_yzzzzzz_0_xyzzzzz_0, g_0_yzzzzzz_0_xyzzzzz_1, g_0_yzzzzzz_0_yyyyyy_1, g_0_yzzzzzz_0_yyyyyyy_0, g_0_yzzzzzz_0_yyyyyyy_1, g_0_yzzzzzz_0_yyyyyyz_0, g_0_yzzzzzz_0_yyyyyyz_1, g_0_yzzzzzz_0_yyyyyz_1, g_0_yzzzzzz_0_yyyyyzz_0, g_0_yzzzzzz_0_yyyyyzz_1, g_0_yzzzzzz_0_yyyyzz_1, g_0_yzzzzzz_0_yyyyzzz_0, g_0_yzzzzzz_0_yyyyzzz_1, g_0_yzzzzzz_0_yyyzzz_1, g_0_yzzzzzz_0_yyyzzzz_0, g_0_yzzzzzz_0_yyyzzzz_1, g_0_yzzzzzz_0_yyzzzz_1, g_0_yzzzzzz_0_yyzzzzz_0, g_0_yzzzzzz_0_yyzzzzz_1, g_0_yzzzzzz_0_yzzzzz_1, g_0_yzzzzzz_0_yzzzzzz_0, g_0_yzzzzzz_0_yzzzzzz_1, g_0_yzzzzzz_0_zzzzzzz_0, g_0_yzzzzzz_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzzz_0_xxxxxxx_0[i] = g_0_xzzzzzz_0_xxxxxxx_0[i] * pb_y + g_0_xzzzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yzzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxxxy_0[i] * pb_x + g_0_yzzzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxxxz_0[i] = g_0_xzzzzzz_0_xxxxxxz_0[i] * pb_y + g_0_xzzzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yzzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxxyy_0[i] * pb_x + g_0_yzzzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yzzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxxyz_0[i] * pb_x + g_0_yzzzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxxzz_0[i] = g_0_xzzzzzz_0_xxxxxzz_0[i] * pb_y + g_0_xzzzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yzzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxyyy_0[i] * pb_x + g_0_yzzzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yzzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxyyz_0[i] * pb_x + g_0_yzzzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yzzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxyzz_0[i] * pb_x + g_0_yzzzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxzzz_0[i] = g_0_xzzzzzz_0_xxxxzzz_0[i] * pb_y + g_0_xzzzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yzzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyyyy_0[i] * pb_x + g_0_yzzzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yzzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyyyz_0[i] * pb_x + g_0_yzzzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yzzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyyzz_0[i] * pb_x + g_0_yzzzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yzzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyzzz_0[i] * pb_x + g_0_yzzzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxzzzz_0[i] = g_0_xzzzzzz_0_xxxzzzz_0[i] * pb_y + g_0_xzzzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yzzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyyyy_0[i] * pb_x + g_0_yzzzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyyyz_0[i] * pb_x + g_0_yzzzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyyzz_0[i] * pb_x + g_0_yzzzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyzzz_0[i] * pb_x + g_0_yzzzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyzzzz_0[i] * pb_x + g_0_yzzzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxzzzzz_0[i] = g_0_xzzzzzz_0_xxzzzzz_0[i] * pb_y + g_0_xzzzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xyyyyyy_0[i] = g_0_yzzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyyyy_0[i] * pb_x + g_0_yzzzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyyyyz_0[i] = g_0_yzzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyyyz_0[i] * pb_x + g_0_yzzzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyyyzz_0[i] = g_0_yzzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyyzz_0[i] * pb_x + g_0_yzzzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyyzzz_0[i] = g_0_yzzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyzzz_0[i] * pb_x + g_0_yzzzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyzzzz_0[i] = g_0_yzzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyzzzz_0[i] * pb_x + g_0_yzzzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyzzzzz_0[i] = g_0_yzzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyzzzzz_0[i] * pb_x + g_0_yzzzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xzzzzzz_0[i] = g_0_xzzzzzz_0_xzzzzzz_0[i] * pb_y + g_0_xzzzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_yyyyyyy_0[i] = g_0_yzzzzzz_0_yyyyyyy_0[i] * pb_x + g_0_yzzzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyyyyz_0[i] = g_0_yzzzzzz_0_yyyyyyz_0[i] * pb_x + g_0_yzzzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyyyzz_0[i] = g_0_yzzzzzz_0_yyyyyzz_0[i] * pb_x + g_0_yzzzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyyzzz_0[i] = g_0_yzzzzzz_0_yyyyzzz_0[i] * pb_x + g_0_yzzzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyzzzz_0[i] = g_0_yzzzzzz_0_yyyzzzz_0[i] * pb_x + g_0_yzzzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyzzzzz_0[i] = g_0_yzzzzzz_0_yyzzzzz_0[i] * pb_x + g_0_yzzzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yzzzzzz_0[i] = g_0_yzzzzzz_0_yzzzzzz_0[i] * pb_x + g_0_yzzzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_zzzzzzz_0[i] = g_0_yzzzzzz_0_zzzzzzz_0[i] * pb_x + g_0_yzzzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1260-1296 components of targeted buffer : SLSK

    auto g_0_xzzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1260);

    auto g_0_xzzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1261);

    auto g_0_xzzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1262);

    auto g_0_xzzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1263);

    auto g_0_xzzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1264);

    auto g_0_xzzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1265);

    auto g_0_xzzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1266);

    auto g_0_xzzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1267);

    auto g_0_xzzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1268);

    auto g_0_xzzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1269);

    auto g_0_xzzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1270);

    auto g_0_xzzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1271);

    auto g_0_xzzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1272);

    auto g_0_xzzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1273);

    auto g_0_xzzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1274);

    auto g_0_xzzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1275);

    auto g_0_xzzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1276);

    auto g_0_xzzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1277);

    auto g_0_xzzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1278);

    auto g_0_xzzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1279);

    auto g_0_xzzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1280);

    auto g_0_xzzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1281);

    auto g_0_xzzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1282);

    auto g_0_xzzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1283);

    auto g_0_xzzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1284);

    auto g_0_xzzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1285);

    auto g_0_xzzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1286);

    auto g_0_xzzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1287);

    auto g_0_xzzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1288);

    auto g_0_xzzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1289);

    auto g_0_xzzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1290);

    auto g_0_xzzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1291);

    auto g_0_xzzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1292);

    auto g_0_xzzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1293);

    auto g_0_xzzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1294);

    auto g_0_xzzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1295);

    #pragma omp simd aligned(g_0_xzzzzzzz_0_xxxxxxx_0, g_0_xzzzzzzz_0_xxxxxxy_0, g_0_xzzzzzzz_0_xxxxxxz_0, g_0_xzzzzzzz_0_xxxxxyy_0, g_0_xzzzzzzz_0_xxxxxyz_0, g_0_xzzzzzzz_0_xxxxxzz_0, g_0_xzzzzzzz_0_xxxxyyy_0, g_0_xzzzzzzz_0_xxxxyyz_0, g_0_xzzzzzzz_0_xxxxyzz_0, g_0_xzzzzzzz_0_xxxxzzz_0, g_0_xzzzzzzz_0_xxxyyyy_0, g_0_xzzzzzzz_0_xxxyyyz_0, g_0_xzzzzzzz_0_xxxyyzz_0, g_0_xzzzzzzz_0_xxxyzzz_0, g_0_xzzzzzzz_0_xxxzzzz_0, g_0_xzzzzzzz_0_xxyyyyy_0, g_0_xzzzzzzz_0_xxyyyyz_0, g_0_xzzzzzzz_0_xxyyyzz_0, g_0_xzzzzzzz_0_xxyyzzz_0, g_0_xzzzzzzz_0_xxyzzzz_0, g_0_xzzzzzzz_0_xxzzzzz_0, g_0_xzzzzzzz_0_xyyyyyy_0, g_0_xzzzzzzz_0_xyyyyyz_0, g_0_xzzzzzzz_0_xyyyyzz_0, g_0_xzzzzzzz_0_xyyyzzz_0, g_0_xzzzzzzz_0_xyyzzzz_0, g_0_xzzzzzzz_0_xyzzzzz_0, g_0_xzzzzzzz_0_xzzzzzz_0, g_0_xzzzzzzz_0_yyyyyyy_0, g_0_xzzzzzzz_0_yyyyyyz_0, g_0_xzzzzzzz_0_yyyyyzz_0, g_0_xzzzzzzz_0_yyyyzzz_0, g_0_xzzzzzzz_0_yyyzzzz_0, g_0_xzzzzzzz_0_yyzzzzz_0, g_0_xzzzzzzz_0_yzzzzzz_0, g_0_xzzzzzzz_0_zzzzzzz_0, g_0_zzzzzzz_0_xxxxxx_1, g_0_zzzzzzz_0_xxxxxxx_0, g_0_zzzzzzz_0_xxxxxxx_1, g_0_zzzzzzz_0_xxxxxxy_0, g_0_zzzzzzz_0_xxxxxxy_1, g_0_zzzzzzz_0_xxxxxxz_0, g_0_zzzzzzz_0_xxxxxxz_1, g_0_zzzzzzz_0_xxxxxy_1, g_0_zzzzzzz_0_xxxxxyy_0, g_0_zzzzzzz_0_xxxxxyy_1, g_0_zzzzzzz_0_xxxxxyz_0, g_0_zzzzzzz_0_xxxxxyz_1, g_0_zzzzzzz_0_xxxxxz_1, g_0_zzzzzzz_0_xxxxxzz_0, g_0_zzzzzzz_0_xxxxxzz_1, g_0_zzzzzzz_0_xxxxyy_1, g_0_zzzzzzz_0_xxxxyyy_0, g_0_zzzzzzz_0_xxxxyyy_1, g_0_zzzzzzz_0_xxxxyyz_0, g_0_zzzzzzz_0_xxxxyyz_1, g_0_zzzzzzz_0_xxxxyz_1, g_0_zzzzzzz_0_xxxxyzz_0, g_0_zzzzzzz_0_xxxxyzz_1, g_0_zzzzzzz_0_xxxxzz_1, g_0_zzzzzzz_0_xxxxzzz_0, g_0_zzzzzzz_0_xxxxzzz_1, g_0_zzzzzzz_0_xxxyyy_1, g_0_zzzzzzz_0_xxxyyyy_0, g_0_zzzzzzz_0_xxxyyyy_1, g_0_zzzzzzz_0_xxxyyyz_0, g_0_zzzzzzz_0_xxxyyyz_1, g_0_zzzzzzz_0_xxxyyz_1, g_0_zzzzzzz_0_xxxyyzz_0, g_0_zzzzzzz_0_xxxyyzz_1, g_0_zzzzzzz_0_xxxyzz_1, g_0_zzzzzzz_0_xxxyzzz_0, g_0_zzzzzzz_0_xxxyzzz_1, g_0_zzzzzzz_0_xxxzzz_1, g_0_zzzzzzz_0_xxxzzzz_0, g_0_zzzzzzz_0_xxxzzzz_1, g_0_zzzzzzz_0_xxyyyy_1, g_0_zzzzzzz_0_xxyyyyy_0, g_0_zzzzzzz_0_xxyyyyy_1, g_0_zzzzzzz_0_xxyyyyz_0, g_0_zzzzzzz_0_xxyyyyz_1, g_0_zzzzzzz_0_xxyyyz_1, g_0_zzzzzzz_0_xxyyyzz_0, g_0_zzzzzzz_0_xxyyyzz_1, g_0_zzzzzzz_0_xxyyzz_1, g_0_zzzzzzz_0_xxyyzzz_0, g_0_zzzzzzz_0_xxyyzzz_1, g_0_zzzzzzz_0_xxyzzz_1, g_0_zzzzzzz_0_xxyzzzz_0, g_0_zzzzzzz_0_xxyzzzz_1, g_0_zzzzzzz_0_xxzzzz_1, g_0_zzzzzzz_0_xxzzzzz_0, g_0_zzzzzzz_0_xxzzzzz_1, g_0_zzzzzzz_0_xyyyyy_1, g_0_zzzzzzz_0_xyyyyyy_0, g_0_zzzzzzz_0_xyyyyyy_1, g_0_zzzzzzz_0_xyyyyyz_0, g_0_zzzzzzz_0_xyyyyyz_1, g_0_zzzzzzz_0_xyyyyz_1, g_0_zzzzzzz_0_xyyyyzz_0, g_0_zzzzzzz_0_xyyyyzz_1, g_0_zzzzzzz_0_xyyyzz_1, g_0_zzzzzzz_0_xyyyzzz_0, g_0_zzzzzzz_0_xyyyzzz_1, g_0_zzzzzzz_0_xyyzzz_1, g_0_zzzzzzz_0_xyyzzzz_0, g_0_zzzzzzz_0_xyyzzzz_1, g_0_zzzzzzz_0_xyzzzz_1, g_0_zzzzzzz_0_xyzzzzz_0, g_0_zzzzzzz_0_xyzzzzz_1, g_0_zzzzzzz_0_xzzzzz_1, g_0_zzzzzzz_0_xzzzzzz_0, g_0_zzzzzzz_0_xzzzzzz_1, g_0_zzzzzzz_0_yyyyyy_1, g_0_zzzzzzz_0_yyyyyyy_0, g_0_zzzzzzz_0_yyyyyyy_1, g_0_zzzzzzz_0_yyyyyyz_0, g_0_zzzzzzz_0_yyyyyyz_1, g_0_zzzzzzz_0_yyyyyz_1, g_0_zzzzzzz_0_yyyyyzz_0, g_0_zzzzzzz_0_yyyyyzz_1, g_0_zzzzzzz_0_yyyyzz_1, g_0_zzzzzzz_0_yyyyzzz_0, g_0_zzzzzzz_0_yyyyzzz_1, g_0_zzzzzzz_0_yyyzzz_1, g_0_zzzzzzz_0_yyyzzzz_0, g_0_zzzzzzz_0_yyyzzzz_1, g_0_zzzzzzz_0_yyzzzz_1, g_0_zzzzzzz_0_yyzzzzz_0, g_0_zzzzzzz_0_yyzzzzz_1, g_0_zzzzzzz_0_yzzzzz_1, g_0_zzzzzzz_0_yzzzzzz_0, g_0_zzzzzzz_0_yzzzzzz_1, g_0_zzzzzzz_0_zzzzzz_1, g_0_zzzzzzz_0_zzzzzzz_0, g_0_zzzzzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzzz_0_xxxxxxx_0[i] = 7.0 * g_0_zzzzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxxx_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_zzzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxxy_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxxxz_0[i] = 6.0 * g_0_zzzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxxz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_zzzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxyy_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_zzzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxyz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxxzz_0[i] = 5.0 * g_0_zzzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_zzzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyyy_0[i] * pb_x + g_0_zzzzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_zzzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyyz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_zzzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxzzz_0[i] = 4.0 * g_0_zzzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_zzzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyyy_0[i] * pb_x + g_0_zzzzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyyz_0[i] * pb_x + g_0_zzzzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_zzzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_zzzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_zzzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_zzzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyyy_0[i] * pb_x + g_0_zzzzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_zzzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyyz_0[i] * pb_x + g_0_zzzzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_zzzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyzz_0[i] * pb_x + g_0_zzzzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_zzzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_zzzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyyyy_0[i] = g_0_zzzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyyy_0[i] * pb_x + g_0_zzzzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyyyz_0[i] = g_0_zzzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyyz_0[i] * pb_x + g_0_zzzzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyyzz_0[i] = g_0_zzzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyzz_0[i] * pb_x + g_0_zzzzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyzzz_0[i] = g_0_zzzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyzzz_0[i] * pb_x + g_0_zzzzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyzzzz_0[i] = g_0_zzzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyzzzzz_0[i] = g_0_zzzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xzzzzzz_0[i] = g_0_zzzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyyyy_0[i] = g_0_zzzzzzz_0_yyyyyyy_0[i] * pb_x + g_0_zzzzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyyyz_0[i] = g_0_zzzzzzz_0_yyyyyyz_0[i] * pb_x + g_0_zzzzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyyzz_0[i] = g_0_zzzzzzz_0_yyyyyzz_0[i] * pb_x + g_0_zzzzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyzzz_0[i] = g_0_zzzzzzz_0_yyyyzzz_0[i] * pb_x + g_0_zzzzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyzzzz_0[i] = g_0_zzzzzzz_0_yyyzzzz_0[i] * pb_x + g_0_zzzzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyzzzzz_0[i] = g_0_zzzzzzz_0_yyzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yzzzzzz_0[i] = g_0_zzzzzzz_0_yzzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_zzzzzzz_0[i] = g_0_zzzzzzz_0_zzzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1296-1332 components of targeted buffer : SLSK

    auto g_0_yyyyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1296);

    auto g_0_yyyyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1297);

    auto g_0_yyyyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1298);

    auto g_0_yyyyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1299);

    auto g_0_yyyyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1300);

    auto g_0_yyyyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1301);

    auto g_0_yyyyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1302);

    auto g_0_yyyyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1303);

    auto g_0_yyyyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1304);

    auto g_0_yyyyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1305);

    auto g_0_yyyyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1306);

    auto g_0_yyyyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1307);

    auto g_0_yyyyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1308);

    auto g_0_yyyyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1309);

    auto g_0_yyyyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1310);

    auto g_0_yyyyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1311);

    auto g_0_yyyyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1312);

    auto g_0_yyyyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1313);

    auto g_0_yyyyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1314);

    auto g_0_yyyyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1315);

    auto g_0_yyyyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1316);

    auto g_0_yyyyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1317);

    auto g_0_yyyyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1318);

    auto g_0_yyyyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1319);

    auto g_0_yyyyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1320);

    auto g_0_yyyyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1321);

    auto g_0_yyyyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1322);

    auto g_0_yyyyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1323);

    auto g_0_yyyyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1324);

    auto g_0_yyyyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1325);

    auto g_0_yyyyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1326);

    auto g_0_yyyyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1327);

    auto g_0_yyyyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1328);

    auto g_0_yyyyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1329);

    auto g_0_yyyyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1330);

    auto g_0_yyyyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1331);

    #pragma omp simd aligned(g_0_yyyyyy_0_xxxxxxx_0, g_0_yyyyyy_0_xxxxxxx_1, g_0_yyyyyy_0_xxxxxxy_0, g_0_yyyyyy_0_xxxxxxy_1, g_0_yyyyyy_0_xxxxxxz_0, g_0_yyyyyy_0_xxxxxxz_1, g_0_yyyyyy_0_xxxxxyy_0, g_0_yyyyyy_0_xxxxxyy_1, g_0_yyyyyy_0_xxxxxyz_0, g_0_yyyyyy_0_xxxxxyz_1, g_0_yyyyyy_0_xxxxxzz_0, g_0_yyyyyy_0_xxxxxzz_1, g_0_yyyyyy_0_xxxxyyy_0, g_0_yyyyyy_0_xxxxyyy_1, g_0_yyyyyy_0_xxxxyyz_0, g_0_yyyyyy_0_xxxxyyz_1, g_0_yyyyyy_0_xxxxyzz_0, g_0_yyyyyy_0_xxxxyzz_1, g_0_yyyyyy_0_xxxxzzz_0, g_0_yyyyyy_0_xxxxzzz_1, g_0_yyyyyy_0_xxxyyyy_0, g_0_yyyyyy_0_xxxyyyy_1, g_0_yyyyyy_0_xxxyyyz_0, g_0_yyyyyy_0_xxxyyyz_1, g_0_yyyyyy_0_xxxyyzz_0, g_0_yyyyyy_0_xxxyyzz_1, g_0_yyyyyy_0_xxxyzzz_0, g_0_yyyyyy_0_xxxyzzz_1, g_0_yyyyyy_0_xxxzzzz_0, g_0_yyyyyy_0_xxxzzzz_1, g_0_yyyyyy_0_xxyyyyy_0, g_0_yyyyyy_0_xxyyyyy_1, g_0_yyyyyy_0_xxyyyyz_0, g_0_yyyyyy_0_xxyyyyz_1, g_0_yyyyyy_0_xxyyyzz_0, g_0_yyyyyy_0_xxyyyzz_1, g_0_yyyyyy_0_xxyyzzz_0, g_0_yyyyyy_0_xxyyzzz_1, g_0_yyyyyy_0_xxyzzzz_0, g_0_yyyyyy_0_xxyzzzz_1, g_0_yyyyyy_0_xxzzzzz_0, g_0_yyyyyy_0_xxzzzzz_1, g_0_yyyyyy_0_xyyyyyy_0, g_0_yyyyyy_0_xyyyyyy_1, g_0_yyyyyy_0_xyyyyyz_0, g_0_yyyyyy_0_xyyyyyz_1, g_0_yyyyyy_0_xyyyyzz_0, g_0_yyyyyy_0_xyyyyzz_1, g_0_yyyyyy_0_xyyyzzz_0, g_0_yyyyyy_0_xyyyzzz_1, g_0_yyyyyy_0_xyyzzzz_0, g_0_yyyyyy_0_xyyzzzz_1, g_0_yyyyyy_0_xyzzzzz_0, g_0_yyyyyy_0_xyzzzzz_1, g_0_yyyyyy_0_xzzzzzz_0, g_0_yyyyyy_0_xzzzzzz_1, g_0_yyyyyy_0_yyyyyyy_0, g_0_yyyyyy_0_yyyyyyy_1, g_0_yyyyyy_0_yyyyyyz_0, g_0_yyyyyy_0_yyyyyyz_1, g_0_yyyyyy_0_yyyyyzz_0, g_0_yyyyyy_0_yyyyyzz_1, g_0_yyyyyy_0_yyyyzzz_0, g_0_yyyyyy_0_yyyyzzz_1, g_0_yyyyyy_0_yyyzzzz_0, g_0_yyyyyy_0_yyyzzzz_1, g_0_yyyyyy_0_yyzzzzz_0, g_0_yyyyyy_0_yyzzzzz_1, g_0_yyyyyy_0_yzzzzzz_0, g_0_yyyyyy_0_yzzzzzz_1, g_0_yyyyyy_0_zzzzzzz_0, g_0_yyyyyy_0_zzzzzzz_1, g_0_yyyyyyy_0_xxxxxx_1, g_0_yyyyyyy_0_xxxxxxx_0, g_0_yyyyyyy_0_xxxxxxx_1, g_0_yyyyyyy_0_xxxxxxy_0, g_0_yyyyyyy_0_xxxxxxy_1, g_0_yyyyyyy_0_xxxxxxz_0, g_0_yyyyyyy_0_xxxxxxz_1, g_0_yyyyyyy_0_xxxxxy_1, g_0_yyyyyyy_0_xxxxxyy_0, g_0_yyyyyyy_0_xxxxxyy_1, g_0_yyyyyyy_0_xxxxxyz_0, g_0_yyyyyyy_0_xxxxxyz_1, g_0_yyyyyyy_0_xxxxxz_1, g_0_yyyyyyy_0_xxxxxzz_0, g_0_yyyyyyy_0_xxxxxzz_1, g_0_yyyyyyy_0_xxxxyy_1, g_0_yyyyyyy_0_xxxxyyy_0, g_0_yyyyyyy_0_xxxxyyy_1, g_0_yyyyyyy_0_xxxxyyz_0, g_0_yyyyyyy_0_xxxxyyz_1, g_0_yyyyyyy_0_xxxxyz_1, g_0_yyyyyyy_0_xxxxyzz_0, g_0_yyyyyyy_0_xxxxyzz_1, g_0_yyyyyyy_0_xxxxzz_1, g_0_yyyyyyy_0_xxxxzzz_0, g_0_yyyyyyy_0_xxxxzzz_1, g_0_yyyyyyy_0_xxxyyy_1, g_0_yyyyyyy_0_xxxyyyy_0, g_0_yyyyyyy_0_xxxyyyy_1, g_0_yyyyyyy_0_xxxyyyz_0, g_0_yyyyyyy_0_xxxyyyz_1, g_0_yyyyyyy_0_xxxyyz_1, g_0_yyyyyyy_0_xxxyyzz_0, g_0_yyyyyyy_0_xxxyyzz_1, g_0_yyyyyyy_0_xxxyzz_1, g_0_yyyyyyy_0_xxxyzzz_0, g_0_yyyyyyy_0_xxxyzzz_1, g_0_yyyyyyy_0_xxxzzz_1, g_0_yyyyyyy_0_xxxzzzz_0, g_0_yyyyyyy_0_xxxzzzz_1, g_0_yyyyyyy_0_xxyyyy_1, g_0_yyyyyyy_0_xxyyyyy_0, g_0_yyyyyyy_0_xxyyyyy_1, g_0_yyyyyyy_0_xxyyyyz_0, g_0_yyyyyyy_0_xxyyyyz_1, g_0_yyyyyyy_0_xxyyyz_1, g_0_yyyyyyy_0_xxyyyzz_0, g_0_yyyyyyy_0_xxyyyzz_1, g_0_yyyyyyy_0_xxyyzz_1, g_0_yyyyyyy_0_xxyyzzz_0, g_0_yyyyyyy_0_xxyyzzz_1, g_0_yyyyyyy_0_xxyzzz_1, g_0_yyyyyyy_0_xxyzzzz_0, g_0_yyyyyyy_0_xxyzzzz_1, g_0_yyyyyyy_0_xxzzzz_1, g_0_yyyyyyy_0_xxzzzzz_0, g_0_yyyyyyy_0_xxzzzzz_1, g_0_yyyyyyy_0_xyyyyy_1, g_0_yyyyyyy_0_xyyyyyy_0, g_0_yyyyyyy_0_xyyyyyy_1, g_0_yyyyyyy_0_xyyyyyz_0, g_0_yyyyyyy_0_xyyyyyz_1, g_0_yyyyyyy_0_xyyyyz_1, g_0_yyyyyyy_0_xyyyyzz_0, g_0_yyyyyyy_0_xyyyyzz_1, g_0_yyyyyyy_0_xyyyzz_1, g_0_yyyyyyy_0_xyyyzzz_0, g_0_yyyyyyy_0_xyyyzzz_1, g_0_yyyyyyy_0_xyyzzz_1, g_0_yyyyyyy_0_xyyzzzz_0, g_0_yyyyyyy_0_xyyzzzz_1, g_0_yyyyyyy_0_xyzzzz_1, g_0_yyyyyyy_0_xyzzzzz_0, g_0_yyyyyyy_0_xyzzzzz_1, g_0_yyyyyyy_0_xzzzzz_1, g_0_yyyyyyy_0_xzzzzzz_0, g_0_yyyyyyy_0_xzzzzzz_1, g_0_yyyyyyy_0_yyyyyy_1, g_0_yyyyyyy_0_yyyyyyy_0, g_0_yyyyyyy_0_yyyyyyy_1, g_0_yyyyyyy_0_yyyyyyz_0, g_0_yyyyyyy_0_yyyyyyz_1, g_0_yyyyyyy_0_yyyyyz_1, g_0_yyyyyyy_0_yyyyyzz_0, g_0_yyyyyyy_0_yyyyyzz_1, g_0_yyyyyyy_0_yyyyzz_1, g_0_yyyyyyy_0_yyyyzzz_0, g_0_yyyyyyy_0_yyyyzzz_1, g_0_yyyyyyy_0_yyyzzz_1, g_0_yyyyyyy_0_yyyzzzz_0, g_0_yyyyyyy_0_yyyzzzz_1, g_0_yyyyyyy_0_yyzzzz_1, g_0_yyyyyyy_0_yyzzzzz_0, g_0_yyyyyyy_0_yyzzzzz_1, g_0_yyyyyyy_0_yzzzzz_1, g_0_yyyyyyy_0_yzzzzzz_0, g_0_yyyyyyy_0_yzzzzzz_1, g_0_yyyyyyy_0_zzzzzz_1, g_0_yyyyyyy_0_zzzzzzz_0, g_0_yyyyyyy_0_zzzzzzz_1, g_0_yyyyyyyy_0_xxxxxxx_0, g_0_yyyyyyyy_0_xxxxxxy_0, g_0_yyyyyyyy_0_xxxxxxz_0, g_0_yyyyyyyy_0_xxxxxyy_0, g_0_yyyyyyyy_0_xxxxxyz_0, g_0_yyyyyyyy_0_xxxxxzz_0, g_0_yyyyyyyy_0_xxxxyyy_0, g_0_yyyyyyyy_0_xxxxyyz_0, g_0_yyyyyyyy_0_xxxxyzz_0, g_0_yyyyyyyy_0_xxxxzzz_0, g_0_yyyyyyyy_0_xxxyyyy_0, g_0_yyyyyyyy_0_xxxyyyz_0, g_0_yyyyyyyy_0_xxxyyzz_0, g_0_yyyyyyyy_0_xxxyzzz_0, g_0_yyyyyyyy_0_xxxzzzz_0, g_0_yyyyyyyy_0_xxyyyyy_0, g_0_yyyyyyyy_0_xxyyyyz_0, g_0_yyyyyyyy_0_xxyyyzz_0, g_0_yyyyyyyy_0_xxyyzzz_0, g_0_yyyyyyyy_0_xxyzzzz_0, g_0_yyyyyyyy_0_xxzzzzz_0, g_0_yyyyyyyy_0_xyyyyyy_0, g_0_yyyyyyyy_0_xyyyyyz_0, g_0_yyyyyyyy_0_xyyyyzz_0, g_0_yyyyyyyy_0_xyyyzzz_0, g_0_yyyyyyyy_0_xyyzzzz_0, g_0_yyyyyyyy_0_xyzzzzz_0, g_0_yyyyyyyy_0_xzzzzzz_0, g_0_yyyyyyyy_0_yyyyyyy_0, g_0_yyyyyyyy_0_yyyyyyz_0, g_0_yyyyyyyy_0_yyyyyzz_0, g_0_yyyyyyyy_0_yyyyzzz_0, g_0_yyyyyyyy_0_yyyzzzz_0, g_0_yyyyyyyy_0_yyzzzzz_0, g_0_yyyyyyyy_0_yzzzzzz_0, g_0_yyyyyyyy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyyy_0_xxxxxxx_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxxx_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxxxx_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxxxy_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxxy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxxy_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxxxz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxxz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxxxz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxxyy_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxyy_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxxyz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxyz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxxzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxxzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxyyy_0[i] = 7.0 * g_0_yyyyyy_0_xxxxyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyyy_0[i] * pb_y + g_0_yyyyyyy_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxyyz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyyz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxyzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyyyy_0[i] = 7.0 * g_0_yyyyyy_0_xxxyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyyy_0[i] * pb_y + g_0_yyyyyyy_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyyyz_0[i] = 7.0 * g_0_yyyyyy_0_xxxyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyyz_0[i] * pb_y + g_0_yyyyyyy_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyyzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyyyy_0[i] = 7.0 * g_0_yyyyyy_0_xxyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyyy_0[i] * pb_y + g_0_yyyyyyy_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyyyz_0[i] = 7.0 * g_0_yyyyyy_0_xxyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyyz_0[i] * pb_y + g_0_yyyyyyy_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyyzz_0[i] = 7.0 * g_0_yyyyyy_0_xxyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyzz_0[i] * pb_y + g_0_yyyyyyy_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyyyy_0[i] = 7.0 * g_0_yyyyyy_0_xyyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_yyyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyyy_0[i] * pb_y + g_0_yyyyyyy_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyyyz_0[i] = 7.0 * g_0_yyyyyy_0_xyyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyyz_0[i] * pb_y + g_0_yyyyyyy_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyyzz_0[i] = 7.0 * g_0_yyyyyy_0_xyyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyzz_0[i] * pb_y + g_0_yyyyyyy_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyzzz_0[i] = 7.0 * g_0_yyyyyy_0_xyyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyzzz_0[i] * pb_y + g_0_yyyyyyy_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xyyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xyzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xzzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xzzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyyyy_0[i] = 7.0 * g_0_yyyyyy_0_yyyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyyyy_1[i] * fti_ab_0 + 7.0 * g_0_yyyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyyyy_0[i] * pb_y + g_0_yyyyyyy_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyyyz_0[i] = 7.0 * g_0_yyyyyy_0_yyyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyyyz_0[i] * pb_y + g_0_yyyyyyy_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyyzz_0[i] = 7.0 * g_0_yyyyyy_0_yyyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyyzz_0[i] * pb_y + g_0_yyyyyyy_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyzzz_0[i] = 7.0 * g_0_yyyyyy_0_yyyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyzzz_0[i] * pb_y + g_0_yyyyyyy_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyzzzz_0[i] = 7.0 * g_0_yyyyyy_0_yyyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyzzzz_0[i] * pb_y + g_0_yyyyyyy_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_yyzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yzzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_yzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_zzzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_zzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zzzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1332-1368 components of targeted buffer : SLSK

    auto g_0_yyyyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1332);

    auto g_0_yyyyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1333);

    auto g_0_yyyyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1334);

    auto g_0_yyyyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1335);

    auto g_0_yyyyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1336);

    auto g_0_yyyyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1337);

    auto g_0_yyyyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1338);

    auto g_0_yyyyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1339);

    auto g_0_yyyyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1340);

    auto g_0_yyyyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1341);

    auto g_0_yyyyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1342);

    auto g_0_yyyyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1343);

    auto g_0_yyyyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1344);

    auto g_0_yyyyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1345);

    auto g_0_yyyyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1346);

    auto g_0_yyyyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1347);

    auto g_0_yyyyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1348);

    auto g_0_yyyyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1349);

    auto g_0_yyyyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1350);

    auto g_0_yyyyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1351);

    auto g_0_yyyyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1352);

    auto g_0_yyyyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1353);

    auto g_0_yyyyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1354);

    auto g_0_yyyyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1355);

    auto g_0_yyyyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1356);

    auto g_0_yyyyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1357);

    auto g_0_yyyyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1358);

    auto g_0_yyyyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1359);

    auto g_0_yyyyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1360);

    auto g_0_yyyyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1361);

    auto g_0_yyyyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1362);

    auto g_0_yyyyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1363);

    auto g_0_yyyyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1364);

    auto g_0_yyyyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1365);

    auto g_0_yyyyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1366);

    auto g_0_yyyyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1367);

    #pragma omp simd aligned(g_0_yyyyyyy_0_xxxxxx_1, g_0_yyyyyyy_0_xxxxxxx_0, g_0_yyyyyyy_0_xxxxxxx_1, g_0_yyyyyyy_0_xxxxxxy_0, g_0_yyyyyyy_0_xxxxxxy_1, g_0_yyyyyyy_0_xxxxxxz_0, g_0_yyyyyyy_0_xxxxxxz_1, g_0_yyyyyyy_0_xxxxxy_1, g_0_yyyyyyy_0_xxxxxyy_0, g_0_yyyyyyy_0_xxxxxyy_1, g_0_yyyyyyy_0_xxxxxyz_0, g_0_yyyyyyy_0_xxxxxyz_1, g_0_yyyyyyy_0_xxxxxz_1, g_0_yyyyyyy_0_xxxxxzz_0, g_0_yyyyyyy_0_xxxxxzz_1, g_0_yyyyyyy_0_xxxxyy_1, g_0_yyyyyyy_0_xxxxyyy_0, g_0_yyyyyyy_0_xxxxyyy_1, g_0_yyyyyyy_0_xxxxyyz_0, g_0_yyyyyyy_0_xxxxyyz_1, g_0_yyyyyyy_0_xxxxyz_1, g_0_yyyyyyy_0_xxxxyzz_0, g_0_yyyyyyy_0_xxxxyzz_1, g_0_yyyyyyy_0_xxxxzz_1, g_0_yyyyyyy_0_xxxxzzz_0, g_0_yyyyyyy_0_xxxxzzz_1, g_0_yyyyyyy_0_xxxyyy_1, g_0_yyyyyyy_0_xxxyyyy_0, g_0_yyyyyyy_0_xxxyyyy_1, g_0_yyyyyyy_0_xxxyyyz_0, g_0_yyyyyyy_0_xxxyyyz_1, g_0_yyyyyyy_0_xxxyyz_1, g_0_yyyyyyy_0_xxxyyzz_0, g_0_yyyyyyy_0_xxxyyzz_1, g_0_yyyyyyy_0_xxxyzz_1, g_0_yyyyyyy_0_xxxyzzz_0, g_0_yyyyyyy_0_xxxyzzz_1, g_0_yyyyyyy_0_xxxzzz_1, g_0_yyyyyyy_0_xxxzzzz_0, g_0_yyyyyyy_0_xxxzzzz_1, g_0_yyyyyyy_0_xxyyyy_1, g_0_yyyyyyy_0_xxyyyyy_0, g_0_yyyyyyy_0_xxyyyyy_1, g_0_yyyyyyy_0_xxyyyyz_0, g_0_yyyyyyy_0_xxyyyyz_1, g_0_yyyyyyy_0_xxyyyz_1, g_0_yyyyyyy_0_xxyyyzz_0, g_0_yyyyyyy_0_xxyyyzz_1, g_0_yyyyyyy_0_xxyyzz_1, g_0_yyyyyyy_0_xxyyzzz_0, g_0_yyyyyyy_0_xxyyzzz_1, g_0_yyyyyyy_0_xxyzzz_1, g_0_yyyyyyy_0_xxyzzzz_0, g_0_yyyyyyy_0_xxyzzzz_1, g_0_yyyyyyy_0_xxzzzz_1, g_0_yyyyyyy_0_xxzzzzz_0, g_0_yyyyyyy_0_xxzzzzz_1, g_0_yyyyyyy_0_xyyyyy_1, g_0_yyyyyyy_0_xyyyyyy_0, g_0_yyyyyyy_0_xyyyyyy_1, g_0_yyyyyyy_0_xyyyyyz_0, g_0_yyyyyyy_0_xyyyyyz_1, g_0_yyyyyyy_0_xyyyyz_1, g_0_yyyyyyy_0_xyyyyzz_0, g_0_yyyyyyy_0_xyyyyzz_1, g_0_yyyyyyy_0_xyyyzz_1, g_0_yyyyyyy_0_xyyyzzz_0, g_0_yyyyyyy_0_xyyyzzz_1, g_0_yyyyyyy_0_xyyzzz_1, g_0_yyyyyyy_0_xyyzzzz_0, g_0_yyyyyyy_0_xyyzzzz_1, g_0_yyyyyyy_0_xyzzzz_1, g_0_yyyyyyy_0_xyzzzzz_0, g_0_yyyyyyy_0_xyzzzzz_1, g_0_yyyyyyy_0_xzzzzz_1, g_0_yyyyyyy_0_xzzzzzz_0, g_0_yyyyyyy_0_xzzzzzz_1, g_0_yyyyyyy_0_yyyyyy_1, g_0_yyyyyyy_0_yyyyyyy_0, g_0_yyyyyyy_0_yyyyyyy_1, g_0_yyyyyyy_0_yyyyyyz_0, g_0_yyyyyyy_0_yyyyyyz_1, g_0_yyyyyyy_0_yyyyyz_1, g_0_yyyyyyy_0_yyyyyzz_0, g_0_yyyyyyy_0_yyyyyzz_1, g_0_yyyyyyy_0_yyyyzz_1, g_0_yyyyyyy_0_yyyyzzz_0, g_0_yyyyyyy_0_yyyyzzz_1, g_0_yyyyyyy_0_yyyzzz_1, g_0_yyyyyyy_0_yyyzzzz_0, g_0_yyyyyyy_0_yyyzzzz_1, g_0_yyyyyyy_0_yyzzzz_1, g_0_yyyyyyy_0_yyzzzzz_0, g_0_yyyyyyy_0_yyzzzzz_1, g_0_yyyyyyy_0_yzzzzz_1, g_0_yyyyyyy_0_yzzzzzz_0, g_0_yyyyyyy_0_yzzzzzz_1, g_0_yyyyyyy_0_zzzzzz_1, g_0_yyyyyyy_0_zzzzzzz_0, g_0_yyyyyyy_0_zzzzzzz_1, g_0_yyyyyyyz_0_xxxxxxx_0, g_0_yyyyyyyz_0_xxxxxxy_0, g_0_yyyyyyyz_0_xxxxxxz_0, g_0_yyyyyyyz_0_xxxxxyy_0, g_0_yyyyyyyz_0_xxxxxyz_0, g_0_yyyyyyyz_0_xxxxxzz_0, g_0_yyyyyyyz_0_xxxxyyy_0, g_0_yyyyyyyz_0_xxxxyyz_0, g_0_yyyyyyyz_0_xxxxyzz_0, g_0_yyyyyyyz_0_xxxxzzz_0, g_0_yyyyyyyz_0_xxxyyyy_0, g_0_yyyyyyyz_0_xxxyyyz_0, g_0_yyyyyyyz_0_xxxyyzz_0, g_0_yyyyyyyz_0_xxxyzzz_0, g_0_yyyyyyyz_0_xxxzzzz_0, g_0_yyyyyyyz_0_xxyyyyy_0, g_0_yyyyyyyz_0_xxyyyyz_0, g_0_yyyyyyyz_0_xxyyyzz_0, g_0_yyyyyyyz_0_xxyyzzz_0, g_0_yyyyyyyz_0_xxyzzzz_0, g_0_yyyyyyyz_0_xxzzzzz_0, g_0_yyyyyyyz_0_xyyyyyy_0, g_0_yyyyyyyz_0_xyyyyyz_0, g_0_yyyyyyyz_0_xyyyyzz_0, g_0_yyyyyyyz_0_xyyyzzz_0, g_0_yyyyyyyz_0_xyyzzzz_0, g_0_yyyyyyyz_0_xyzzzzz_0, g_0_yyyyyyyz_0_xzzzzzz_0, g_0_yyyyyyyz_0_yyyyyyy_0, g_0_yyyyyyyz_0_yyyyyyz_0, g_0_yyyyyyyz_0_yyyyyzz_0, g_0_yyyyyyyz_0_yyyyzzz_0, g_0_yyyyyyyz_0_yyyzzzz_0, g_0_yyyyyyyz_0_yyzzzzz_0, g_0_yyyyyyyz_0_yzzzzzz_0, g_0_yyyyyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyyz_0_xxxxxxx_0[i] = g_0_yyyyyyy_0_xxxxxxx_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxxxy_0[i] = g_0_yyyyyyy_0_xxxxxxy_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxxxz_0[i] = g_0_yyyyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxxz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxxyy_0[i] = g_0_yyyyyyy_0_xxxxxyy_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxxyz_0[i] = g_0_yyyyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxyz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_yyyyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxyyy_0[i] = g_0_yyyyyyy_0_xxxxyyy_0[i] * pb_z + g_0_yyyyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxyyz_0[i] = g_0_yyyyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyyz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyyyy_0[i] = g_0_yyyyyyy_0_xxxyyyy_0[i] * pb_z + g_0_yyyyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyyyz_0[i] = g_0_yyyyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyyz_0[i] * pb_z + g_0_yyyyyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_yyyyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyyyy_0[i] = g_0_yyyyyyy_0_xxyyyyy_0[i] * pb_z + g_0_yyyyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyyyz_0[i] = g_0_yyyyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyyz_0[i] * pb_z + g_0_yyyyyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyzz_0[i] * pb_z + g_0_yyyyyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_yyyyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_yyyyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyyyy_0[i] = g_0_yyyyyyy_0_xyyyyyy_0[i] * pb_z + g_0_yyyyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyyyz_0[i] = g_0_yyyyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyyz_0[i] * pb_z + g_0_yyyyyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyzz_0[i] * pb_z + g_0_yyyyyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyzzz_0[i] * pb_z + g_0_yyyyyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_yyyyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_yyyyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_yyyyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyyyy_0[i] = g_0_yyyyyyy_0_yyyyyyy_0[i] * pb_z + g_0_yyyyyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyyyz_0[i] = g_0_yyyyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyyyz_0[i] * pb_z + g_0_yyyyyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_yyyyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyyzz_0[i] * pb_z + g_0_yyyyyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_yyyyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyzzz_0[i] * pb_z + g_0_yyyyyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_yyyyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyzzzz_0[i] * pb_z + g_0_yyyyyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_yyyyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_yyyyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_yyyyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_zzzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 1368-1404 components of targeted buffer : SLSK

    auto g_0_yyyyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1368);

    auto g_0_yyyyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1369);

    auto g_0_yyyyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1370);

    auto g_0_yyyyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1371);

    auto g_0_yyyyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1372);

    auto g_0_yyyyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1373);

    auto g_0_yyyyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1374);

    auto g_0_yyyyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1375);

    auto g_0_yyyyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1376);

    auto g_0_yyyyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1377);

    auto g_0_yyyyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1378);

    auto g_0_yyyyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1379);

    auto g_0_yyyyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1380);

    auto g_0_yyyyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1381);

    auto g_0_yyyyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1382);

    auto g_0_yyyyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1383);

    auto g_0_yyyyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1384);

    auto g_0_yyyyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1385);

    auto g_0_yyyyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1386);

    auto g_0_yyyyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1387);

    auto g_0_yyyyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1388);

    auto g_0_yyyyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1389);

    auto g_0_yyyyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1390);

    auto g_0_yyyyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1391);

    auto g_0_yyyyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1392);

    auto g_0_yyyyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1393);

    auto g_0_yyyyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1394);

    auto g_0_yyyyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1395);

    auto g_0_yyyyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1396);

    auto g_0_yyyyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1397);

    auto g_0_yyyyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1398);

    auto g_0_yyyyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1399);

    auto g_0_yyyyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1400);

    auto g_0_yyyyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1401);

    auto g_0_yyyyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1402);

    auto g_0_yyyyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1403);

    #pragma omp simd aligned(g_0_yyyyyy_0_xxxxxxy_0, g_0_yyyyyy_0_xxxxxxy_1, g_0_yyyyyy_0_xxxxxyy_0, g_0_yyyyyy_0_xxxxxyy_1, g_0_yyyyyy_0_xxxxyyy_0, g_0_yyyyyy_0_xxxxyyy_1, g_0_yyyyyy_0_xxxyyyy_0, g_0_yyyyyy_0_xxxyyyy_1, g_0_yyyyyy_0_xxyyyyy_0, g_0_yyyyyy_0_xxyyyyy_1, g_0_yyyyyy_0_xyyyyyy_0, g_0_yyyyyy_0_xyyyyyy_1, g_0_yyyyyy_0_yyyyyyy_0, g_0_yyyyyy_0_yyyyyyy_1, g_0_yyyyyyz_0_xxxxxxy_0, g_0_yyyyyyz_0_xxxxxxy_1, g_0_yyyyyyz_0_xxxxxyy_0, g_0_yyyyyyz_0_xxxxxyy_1, g_0_yyyyyyz_0_xxxxyyy_0, g_0_yyyyyyz_0_xxxxyyy_1, g_0_yyyyyyz_0_xxxyyyy_0, g_0_yyyyyyz_0_xxxyyyy_1, g_0_yyyyyyz_0_xxyyyyy_0, g_0_yyyyyyz_0_xxyyyyy_1, g_0_yyyyyyz_0_xyyyyyy_0, g_0_yyyyyyz_0_xyyyyyy_1, g_0_yyyyyyz_0_yyyyyyy_0, g_0_yyyyyyz_0_yyyyyyy_1, g_0_yyyyyyzz_0_xxxxxxx_0, g_0_yyyyyyzz_0_xxxxxxy_0, g_0_yyyyyyzz_0_xxxxxxz_0, g_0_yyyyyyzz_0_xxxxxyy_0, g_0_yyyyyyzz_0_xxxxxyz_0, g_0_yyyyyyzz_0_xxxxxzz_0, g_0_yyyyyyzz_0_xxxxyyy_0, g_0_yyyyyyzz_0_xxxxyyz_0, g_0_yyyyyyzz_0_xxxxyzz_0, g_0_yyyyyyzz_0_xxxxzzz_0, g_0_yyyyyyzz_0_xxxyyyy_0, g_0_yyyyyyzz_0_xxxyyyz_0, g_0_yyyyyyzz_0_xxxyyzz_0, g_0_yyyyyyzz_0_xxxyzzz_0, g_0_yyyyyyzz_0_xxxzzzz_0, g_0_yyyyyyzz_0_xxyyyyy_0, g_0_yyyyyyzz_0_xxyyyyz_0, g_0_yyyyyyzz_0_xxyyyzz_0, g_0_yyyyyyzz_0_xxyyzzz_0, g_0_yyyyyyzz_0_xxyzzzz_0, g_0_yyyyyyzz_0_xxzzzzz_0, g_0_yyyyyyzz_0_xyyyyyy_0, g_0_yyyyyyzz_0_xyyyyyz_0, g_0_yyyyyyzz_0_xyyyyzz_0, g_0_yyyyyyzz_0_xyyyzzz_0, g_0_yyyyyyzz_0_xyyzzzz_0, g_0_yyyyyyzz_0_xyzzzzz_0, g_0_yyyyyyzz_0_xzzzzzz_0, g_0_yyyyyyzz_0_yyyyyyy_0, g_0_yyyyyyzz_0_yyyyyyz_0, g_0_yyyyyyzz_0_yyyyyzz_0, g_0_yyyyyyzz_0_yyyyzzz_0, g_0_yyyyyyzz_0_yyyzzzz_0, g_0_yyyyyyzz_0_yyzzzzz_0, g_0_yyyyyyzz_0_yzzzzzz_0, g_0_yyyyyyzz_0_zzzzzzz_0, g_0_yyyyyzz_0_xxxxxxx_0, g_0_yyyyyzz_0_xxxxxxx_1, g_0_yyyyyzz_0_xxxxxxz_0, g_0_yyyyyzz_0_xxxxxxz_1, g_0_yyyyyzz_0_xxxxxyz_0, g_0_yyyyyzz_0_xxxxxyz_1, g_0_yyyyyzz_0_xxxxxz_1, g_0_yyyyyzz_0_xxxxxzz_0, g_0_yyyyyzz_0_xxxxxzz_1, g_0_yyyyyzz_0_xxxxyyz_0, g_0_yyyyyzz_0_xxxxyyz_1, g_0_yyyyyzz_0_xxxxyz_1, g_0_yyyyyzz_0_xxxxyzz_0, g_0_yyyyyzz_0_xxxxyzz_1, g_0_yyyyyzz_0_xxxxzz_1, g_0_yyyyyzz_0_xxxxzzz_0, g_0_yyyyyzz_0_xxxxzzz_1, g_0_yyyyyzz_0_xxxyyyz_0, g_0_yyyyyzz_0_xxxyyyz_1, g_0_yyyyyzz_0_xxxyyz_1, g_0_yyyyyzz_0_xxxyyzz_0, g_0_yyyyyzz_0_xxxyyzz_1, g_0_yyyyyzz_0_xxxyzz_1, g_0_yyyyyzz_0_xxxyzzz_0, g_0_yyyyyzz_0_xxxyzzz_1, g_0_yyyyyzz_0_xxxzzz_1, g_0_yyyyyzz_0_xxxzzzz_0, g_0_yyyyyzz_0_xxxzzzz_1, g_0_yyyyyzz_0_xxyyyyz_0, g_0_yyyyyzz_0_xxyyyyz_1, g_0_yyyyyzz_0_xxyyyz_1, g_0_yyyyyzz_0_xxyyyzz_0, g_0_yyyyyzz_0_xxyyyzz_1, g_0_yyyyyzz_0_xxyyzz_1, g_0_yyyyyzz_0_xxyyzzz_0, g_0_yyyyyzz_0_xxyyzzz_1, g_0_yyyyyzz_0_xxyzzz_1, g_0_yyyyyzz_0_xxyzzzz_0, g_0_yyyyyzz_0_xxyzzzz_1, g_0_yyyyyzz_0_xxzzzz_1, g_0_yyyyyzz_0_xxzzzzz_0, g_0_yyyyyzz_0_xxzzzzz_1, g_0_yyyyyzz_0_xyyyyyz_0, g_0_yyyyyzz_0_xyyyyyz_1, g_0_yyyyyzz_0_xyyyyz_1, g_0_yyyyyzz_0_xyyyyzz_0, g_0_yyyyyzz_0_xyyyyzz_1, g_0_yyyyyzz_0_xyyyzz_1, g_0_yyyyyzz_0_xyyyzzz_0, g_0_yyyyyzz_0_xyyyzzz_1, g_0_yyyyyzz_0_xyyzzz_1, g_0_yyyyyzz_0_xyyzzzz_0, g_0_yyyyyzz_0_xyyzzzz_1, g_0_yyyyyzz_0_xyzzzz_1, g_0_yyyyyzz_0_xyzzzzz_0, g_0_yyyyyzz_0_xyzzzzz_1, g_0_yyyyyzz_0_xzzzzz_1, g_0_yyyyyzz_0_xzzzzzz_0, g_0_yyyyyzz_0_xzzzzzz_1, g_0_yyyyyzz_0_yyyyyyz_0, g_0_yyyyyzz_0_yyyyyyz_1, g_0_yyyyyzz_0_yyyyyz_1, g_0_yyyyyzz_0_yyyyyzz_0, g_0_yyyyyzz_0_yyyyyzz_1, g_0_yyyyyzz_0_yyyyzz_1, g_0_yyyyyzz_0_yyyyzzz_0, g_0_yyyyyzz_0_yyyyzzz_1, g_0_yyyyyzz_0_yyyzzz_1, g_0_yyyyyzz_0_yyyzzzz_0, g_0_yyyyyzz_0_yyyzzzz_1, g_0_yyyyyzz_0_yyzzzz_1, g_0_yyyyyzz_0_yyzzzzz_0, g_0_yyyyyzz_0_yyzzzzz_1, g_0_yyyyyzz_0_yzzzzz_1, g_0_yyyyyzz_0_yzzzzzz_0, g_0_yyyyyzz_0_yzzzzzz_1, g_0_yyyyyzz_0_zzzzzz_1, g_0_yyyyyzz_0_zzzzzzz_0, g_0_yyyyyzz_0_zzzzzzz_1, g_0_yyyyzz_0_xxxxxxx_0, g_0_yyyyzz_0_xxxxxxx_1, g_0_yyyyzz_0_xxxxxxz_0, g_0_yyyyzz_0_xxxxxxz_1, g_0_yyyyzz_0_xxxxxyz_0, g_0_yyyyzz_0_xxxxxyz_1, g_0_yyyyzz_0_xxxxxzz_0, g_0_yyyyzz_0_xxxxxzz_1, g_0_yyyyzz_0_xxxxyyz_0, g_0_yyyyzz_0_xxxxyyz_1, g_0_yyyyzz_0_xxxxyzz_0, g_0_yyyyzz_0_xxxxyzz_1, g_0_yyyyzz_0_xxxxzzz_0, g_0_yyyyzz_0_xxxxzzz_1, g_0_yyyyzz_0_xxxyyyz_0, g_0_yyyyzz_0_xxxyyyz_1, g_0_yyyyzz_0_xxxyyzz_0, g_0_yyyyzz_0_xxxyyzz_1, g_0_yyyyzz_0_xxxyzzz_0, g_0_yyyyzz_0_xxxyzzz_1, g_0_yyyyzz_0_xxxzzzz_0, g_0_yyyyzz_0_xxxzzzz_1, g_0_yyyyzz_0_xxyyyyz_0, g_0_yyyyzz_0_xxyyyyz_1, g_0_yyyyzz_0_xxyyyzz_0, g_0_yyyyzz_0_xxyyyzz_1, g_0_yyyyzz_0_xxyyzzz_0, g_0_yyyyzz_0_xxyyzzz_1, g_0_yyyyzz_0_xxyzzzz_0, g_0_yyyyzz_0_xxyzzzz_1, g_0_yyyyzz_0_xxzzzzz_0, g_0_yyyyzz_0_xxzzzzz_1, g_0_yyyyzz_0_xyyyyyz_0, g_0_yyyyzz_0_xyyyyyz_1, g_0_yyyyzz_0_xyyyyzz_0, g_0_yyyyzz_0_xyyyyzz_1, g_0_yyyyzz_0_xyyyzzz_0, g_0_yyyyzz_0_xyyyzzz_1, g_0_yyyyzz_0_xyyzzzz_0, g_0_yyyyzz_0_xyyzzzz_1, g_0_yyyyzz_0_xyzzzzz_0, g_0_yyyyzz_0_xyzzzzz_1, g_0_yyyyzz_0_xzzzzzz_0, g_0_yyyyzz_0_xzzzzzz_1, g_0_yyyyzz_0_yyyyyyz_0, g_0_yyyyzz_0_yyyyyyz_1, g_0_yyyyzz_0_yyyyyzz_0, g_0_yyyyzz_0_yyyyyzz_1, g_0_yyyyzz_0_yyyyzzz_0, g_0_yyyyzz_0_yyyyzzz_1, g_0_yyyyzz_0_yyyzzzz_0, g_0_yyyyzz_0_yyyzzzz_1, g_0_yyyyzz_0_yyzzzzz_0, g_0_yyyyzz_0_yyzzzzz_1, g_0_yyyyzz_0_yzzzzzz_0, g_0_yyyyzz_0_yzzzzzz_1, g_0_yyyyzz_0_zzzzzzz_0, g_0_yyyyzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyzz_0_xxxxxxx_0[i] = 5.0 * g_0_yyyyzz_0_xxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxxxx_0[i] * pb_y + g_0_yyyyyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxxxy_0[i] = g_0_yyyyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxxxxy_0[i] * pb_z + g_0_yyyyyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxxxxz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxxxz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxxyy_0[i] = g_0_yyyyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxxxyy_0[i] * pb_z + g_0_yyyyyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxyz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxxzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxyyy_0[i] = g_0_yyyyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxxyyy_0[i] * pb_z + g_0_yyyyyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxxyyz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxyyz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxyzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxyzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxyyyy_0[i] = g_0_yyyyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxyyyy_0[i] * pb_z + g_0_yyyyyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxyyyz_0[i] = 5.0 * g_0_yyyyzz_0_xxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyyyz_0[i] * pb_y + g_0_yyyyyzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxyyzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyyzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxyzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyyyyy_0[i] = g_0_yyyyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxyyyyy_0[i] * pb_z + g_0_yyyyyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxyyyyz_0[i] = 5.0 * g_0_yyyyzz_0_xxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyyyz_0[i] * pb_y + g_0_yyyyyzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyyyzz_0[i] = 5.0 * g_0_yyyyzz_0_xxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyyzz_0[i] * pb_y + g_0_yyyyyzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyyzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyyyyy_0[i] = g_0_yyyyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xyyyyyy_0[i] * pb_z + g_0_yyyyyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xyyyyyz_0[i] = 5.0 * g_0_yyyyzz_0_xyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyyyz_0[i] * pb_y + g_0_yyyyyzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyyyzz_0[i] = 5.0 * g_0_yyyyzz_0_xyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyyzz_0[i] * pb_y + g_0_yyyyyzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyyzzz_0[i] = 5.0 * g_0_yyyyzz_0_xyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyzzz_0[i] * pb_y + g_0_yyyyyzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xzzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xzzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyyyyy_0[i] = g_0_yyyyyy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_yyyyyyy_0[i] * pb_z + g_0_yyyyyyz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_yyyyyyz_0[i] = 5.0 * g_0_yyyyzz_0_yyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyyyyz_0[i] * pb_y + g_0_yyyyyzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyyyzz_0[i] = 5.0 * g_0_yyyyzz_0_yyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyyyzz_0[i] * pb_y + g_0_yyyyyzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyyzzz_0[i] = 5.0 * g_0_yyyyzz_0_yyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyyzzz_0[i] * pb_y + g_0_yyyyyzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyzzzz_0[i] = 5.0 * g_0_yyyyzz_0_yyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyzzzz_0[i] * pb_y + g_0_yyyyyzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_yyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yzzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_yzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yzzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_zzzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_zzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zzzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1404-1440 components of targeted buffer : SLSK

    auto g_0_yyyyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1404);

    auto g_0_yyyyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1405);

    auto g_0_yyyyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1406);

    auto g_0_yyyyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1407);

    auto g_0_yyyyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1408);

    auto g_0_yyyyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1409);

    auto g_0_yyyyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1410);

    auto g_0_yyyyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1411);

    auto g_0_yyyyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1412);

    auto g_0_yyyyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1413);

    auto g_0_yyyyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1414);

    auto g_0_yyyyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1415);

    auto g_0_yyyyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1416);

    auto g_0_yyyyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1417);

    auto g_0_yyyyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1418);

    auto g_0_yyyyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1419);

    auto g_0_yyyyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1420);

    auto g_0_yyyyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1421);

    auto g_0_yyyyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1422);

    auto g_0_yyyyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1423);

    auto g_0_yyyyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1424);

    auto g_0_yyyyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1425);

    auto g_0_yyyyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1426);

    auto g_0_yyyyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1427);

    auto g_0_yyyyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1428);

    auto g_0_yyyyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1429);

    auto g_0_yyyyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1430);

    auto g_0_yyyyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1431);

    auto g_0_yyyyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1432);

    auto g_0_yyyyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1433);

    auto g_0_yyyyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1434);

    auto g_0_yyyyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1435);

    auto g_0_yyyyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1436);

    auto g_0_yyyyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1437);

    auto g_0_yyyyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1438);

    auto g_0_yyyyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1439);

    #pragma omp simd aligned(g_0_yyyyyz_0_xxxxxxy_0, g_0_yyyyyz_0_xxxxxxy_1, g_0_yyyyyz_0_xxxxxyy_0, g_0_yyyyyz_0_xxxxxyy_1, g_0_yyyyyz_0_xxxxyyy_0, g_0_yyyyyz_0_xxxxyyy_1, g_0_yyyyyz_0_xxxyyyy_0, g_0_yyyyyz_0_xxxyyyy_1, g_0_yyyyyz_0_xxyyyyy_0, g_0_yyyyyz_0_xxyyyyy_1, g_0_yyyyyz_0_xyyyyyy_0, g_0_yyyyyz_0_xyyyyyy_1, g_0_yyyyyz_0_yyyyyyy_0, g_0_yyyyyz_0_yyyyyyy_1, g_0_yyyyyzz_0_xxxxxxy_0, g_0_yyyyyzz_0_xxxxxxy_1, g_0_yyyyyzz_0_xxxxxyy_0, g_0_yyyyyzz_0_xxxxxyy_1, g_0_yyyyyzz_0_xxxxyyy_0, g_0_yyyyyzz_0_xxxxyyy_1, g_0_yyyyyzz_0_xxxyyyy_0, g_0_yyyyyzz_0_xxxyyyy_1, g_0_yyyyyzz_0_xxyyyyy_0, g_0_yyyyyzz_0_xxyyyyy_1, g_0_yyyyyzz_0_xyyyyyy_0, g_0_yyyyyzz_0_xyyyyyy_1, g_0_yyyyyzz_0_yyyyyyy_0, g_0_yyyyyzz_0_yyyyyyy_1, g_0_yyyyyzzz_0_xxxxxxx_0, g_0_yyyyyzzz_0_xxxxxxy_0, g_0_yyyyyzzz_0_xxxxxxz_0, g_0_yyyyyzzz_0_xxxxxyy_0, g_0_yyyyyzzz_0_xxxxxyz_0, g_0_yyyyyzzz_0_xxxxxzz_0, g_0_yyyyyzzz_0_xxxxyyy_0, g_0_yyyyyzzz_0_xxxxyyz_0, g_0_yyyyyzzz_0_xxxxyzz_0, g_0_yyyyyzzz_0_xxxxzzz_0, g_0_yyyyyzzz_0_xxxyyyy_0, g_0_yyyyyzzz_0_xxxyyyz_0, g_0_yyyyyzzz_0_xxxyyzz_0, g_0_yyyyyzzz_0_xxxyzzz_0, g_0_yyyyyzzz_0_xxxzzzz_0, g_0_yyyyyzzz_0_xxyyyyy_0, g_0_yyyyyzzz_0_xxyyyyz_0, g_0_yyyyyzzz_0_xxyyyzz_0, g_0_yyyyyzzz_0_xxyyzzz_0, g_0_yyyyyzzz_0_xxyzzzz_0, g_0_yyyyyzzz_0_xxzzzzz_0, g_0_yyyyyzzz_0_xyyyyyy_0, g_0_yyyyyzzz_0_xyyyyyz_0, g_0_yyyyyzzz_0_xyyyyzz_0, g_0_yyyyyzzz_0_xyyyzzz_0, g_0_yyyyyzzz_0_xyyzzzz_0, g_0_yyyyyzzz_0_xyzzzzz_0, g_0_yyyyyzzz_0_xzzzzzz_0, g_0_yyyyyzzz_0_yyyyyyy_0, g_0_yyyyyzzz_0_yyyyyyz_0, g_0_yyyyyzzz_0_yyyyyzz_0, g_0_yyyyyzzz_0_yyyyzzz_0, g_0_yyyyyzzz_0_yyyzzzz_0, g_0_yyyyyzzz_0_yyzzzzz_0, g_0_yyyyyzzz_0_yzzzzzz_0, g_0_yyyyyzzz_0_zzzzzzz_0, g_0_yyyyzzz_0_xxxxxxx_0, g_0_yyyyzzz_0_xxxxxxx_1, g_0_yyyyzzz_0_xxxxxxz_0, g_0_yyyyzzz_0_xxxxxxz_1, g_0_yyyyzzz_0_xxxxxyz_0, g_0_yyyyzzz_0_xxxxxyz_1, g_0_yyyyzzz_0_xxxxxz_1, g_0_yyyyzzz_0_xxxxxzz_0, g_0_yyyyzzz_0_xxxxxzz_1, g_0_yyyyzzz_0_xxxxyyz_0, g_0_yyyyzzz_0_xxxxyyz_1, g_0_yyyyzzz_0_xxxxyz_1, g_0_yyyyzzz_0_xxxxyzz_0, g_0_yyyyzzz_0_xxxxyzz_1, g_0_yyyyzzz_0_xxxxzz_1, g_0_yyyyzzz_0_xxxxzzz_0, g_0_yyyyzzz_0_xxxxzzz_1, g_0_yyyyzzz_0_xxxyyyz_0, g_0_yyyyzzz_0_xxxyyyz_1, g_0_yyyyzzz_0_xxxyyz_1, g_0_yyyyzzz_0_xxxyyzz_0, g_0_yyyyzzz_0_xxxyyzz_1, g_0_yyyyzzz_0_xxxyzz_1, g_0_yyyyzzz_0_xxxyzzz_0, g_0_yyyyzzz_0_xxxyzzz_1, g_0_yyyyzzz_0_xxxzzz_1, g_0_yyyyzzz_0_xxxzzzz_0, g_0_yyyyzzz_0_xxxzzzz_1, g_0_yyyyzzz_0_xxyyyyz_0, g_0_yyyyzzz_0_xxyyyyz_1, g_0_yyyyzzz_0_xxyyyz_1, g_0_yyyyzzz_0_xxyyyzz_0, g_0_yyyyzzz_0_xxyyyzz_1, g_0_yyyyzzz_0_xxyyzz_1, g_0_yyyyzzz_0_xxyyzzz_0, g_0_yyyyzzz_0_xxyyzzz_1, g_0_yyyyzzz_0_xxyzzz_1, g_0_yyyyzzz_0_xxyzzzz_0, g_0_yyyyzzz_0_xxyzzzz_1, g_0_yyyyzzz_0_xxzzzz_1, g_0_yyyyzzz_0_xxzzzzz_0, g_0_yyyyzzz_0_xxzzzzz_1, g_0_yyyyzzz_0_xyyyyyz_0, g_0_yyyyzzz_0_xyyyyyz_1, g_0_yyyyzzz_0_xyyyyz_1, g_0_yyyyzzz_0_xyyyyzz_0, g_0_yyyyzzz_0_xyyyyzz_1, g_0_yyyyzzz_0_xyyyzz_1, g_0_yyyyzzz_0_xyyyzzz_0, g_0_yyyyzzz_0_xyyyzzz_1, g_0_yyyyzzz_0_xyyzzz_1, g_0_yyyyzzz_0_xyyzzzz_0, g_0_yyyyzzz_0_xyyzzzz_1, g_0_yyyyzzz_0_xyzzzz_1, g_0_yyyyzzz_0_xyzzzzz_0, g_0_yyyyzzz_0_xyzzzzz_1, g_0_yyyyzzz_0_xzzzzz_1, g_0_yyyyzzz_0_xzzzzzz_0, g_0_yyyyzzz_0_xzzzzzz_1, g_0_yyyyzzz_0_yyyyyyz_0, g_0_yyyyzzz_0_yyyyyyz_1, g_0_yyyyzzz_0_yyyyyz_1, g_0_yyyyzzz_0_yyyyyzz_0, g_0_yyyyzzz_0_yyyyyzz_1, g_0_yyyyzzz_0_yyyyzz_1, g_0_yyyyzzz_0_yyyyzzz_0, g_0_yyyyzzz_0_yyyyzzz_1, g_0_yyyyzzz_0_yyyzzz_1, g_0_yyyyzzz_0_yyyzzzz_0, g_0_yyyyzzz_0_yyyzzzz_1, g_0_yyyyzzz_0_yyzzzz_1, g_0_yyyyzzz_0_yyzzzzz_0, g_0_yyyyzzz_0_yyzzzzz_1, g_0_yyyyzzz_0_yzzzzz_1, g_0_yyyyzzz_0_yzzzzzz_0, g_0_yyyyzzz_0_yzzzzzz_1, g_0_yyyyzzz_0_zzzzzz_1, g_0_yyyyzzz_0_zzzzzzz_0, g_0_yyyyzzz_0_zzzzzzz_1, g_0_yyyzzz_0_xxxxxxx_0, g_0_yyyzzz_0_xxxxxxx_1, g_0_yyyzzz_0_xxxxxxz_0, g_0_yyyzzz_0_xxxxxxz_1, g_0_yyyzzz_0_xxxxxyz_0, g_0_yyyzzz_0_xxxxxyz_1, g_0_yyyzzz_0_xxxxxzz_0, g_0_yyyzzz_0_xxxxxzz_1, g_0_yyyzzz_0_xxxxyyz_0, g_0_yyyzzz_0_xxxxyyz_1, g_0_yyyzzz_0_xxxxyzz_0, g_0_yyyzzz_0_xxxxyzz_1, g_0_yyyzzz_0_xxxxzzz_0, g_0_yyyzzz_0_xxxxzzz_1, g_0_yyyzzz_0_xxxyyyz_0, g_0_yyyzzz_0_xxxyyyz_1, g_0_yyyzzz_0_xxxyyzz_0, g_0_yyyzzz_0_xxxyyzz_1, g_0_yyyzzz_0_xxxyzzz_0, g_0_yyyzzz_0_xxxyzzz_1, g_0_yyyzzz_0_xxxzzzz_0, g_0_yyyzzz_0_xxxzzzz_1, g_0_yyyzzz_0_xxyyyyz_0, g_0_yyyzzz_0_xxyyyyz_1, g_0_yyyzzz_0_xxyyyzz_0, g_0_yyyzzz_0_xxyyyzz_1, g_0_yyyzzz_0_xxyyzzz_0, g_0_yyyzzz_0_xxyyzzz_1, g_0_yyyzzz_0_xxyzzzz_0, g_0_yyyzzz_0_xxyzzzz_1, g_0_yyyzzz_0_xxzzzzz_0, g_0_yyyzzz_0_xxzzzzz_1, g_0_yyyzzz_0_xyyyyyz_0, g_0_yyyzzz_0_xyyyyyz_1, g_0_yyyzzz_0_xyyyyzz_0, g_0_yyyzzz_0_xyyyyzz_1, g_0_yyyzzz_0_xyyyzzz_0, g_0_yyyzzz_0_xyyyzzz_1, g_0_yyyzzz_0_xyyzzzz_0, g_0_yyyzzz_0_xyyzzzz_1, g_0_yyyzzz_0_xyzzzzz_0, g_0_yyyzzz_0_xyzzzzz_1, g_0_yyyzzz_0_xzzzzzz_0, g_0_yyyzzz_0_xzzzzzz_1, g_0_yyyzzz_0_yyyyyyz_0, g_0_yyyzzz_0_yyyyyyz_1, g_0_yyyzzz_0_yyyyyzz_0, g_0_yyyzzz_0_yyyyyzz_1, g_0_yyyzzz_0_yyyyzzz_0, g_0_yyyzzz_0_yyyyzzz_1, g_0_yyyzzz_0_yyyzzzz_0, g_0_yyyzzz_0_yyyzzzz_1, g_0_yyyzzz_0_yyzzzzz_0, g_0_yyyzzz_0_yyzzzzz_1, g_0_yyyzzz_0_yzzzzzz_0, g_0_yyyzzz_0_yzzzzzz_1, g_0_yyyzzz_0_zzzzzzz_0, g_0_yyyzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzzz_0_xxxxxxx_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxxxx_0[i] * pb_y + g_0_yyyyzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxxxy_0[i] = 2.0 * g_0_yyyyyz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxxxy_0[i] * pb_z + g_0_yyyyyzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxxxxz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxxxz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_yyyyyz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxxyy_0[i] * pb_z + g_0_yyyyyzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxxxyz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxyz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxxzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxxzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxyyy_0[i] = 2.0 * g_0_yyyyyz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxyyy_0[i] * pb_z + g_0_yyyyyzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxyyz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxyzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxyyyy_0[i] = 2.0 * g_0_yyyyyz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxyyyy_0[i] * pb_z + g_0_yyyyyzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxyyyz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyyyz_0[i] * pb_y + g_0_yyyyzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxyyzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyyzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxyzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyyyz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxyyyyy_0[i] * pb_z + g_0_yyyyyzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxyyyyz_0[i] = 4.0 * g_0_yyyzzz_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyyyz_0[i] * pb_y + g_0_yyyyzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyyyzz_0[i] = 4.0 * g_0_yyyzzz_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyyzz_0[i] * pb_y + g_0_yyyyzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyyzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyyyyy_0[i] = 2.0 * g_0_yyyyyz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xyyyyyy_0[i] * pb_z + g_0_yyyyyzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xyyyyyz_0[i] = 4.0 * g_0_yyyzzz_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyyyz_0[i] * pb_y + g_0_yyyyzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyyyzz_0[i] = 4.0 * g_0_yyyzzz_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyyzz_0[i] * pb_y + g_0_yyyyzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyyzzz_0[i] = 4.0 * g_0_yyyzzz_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyzzz_0[i] * pb_y + g_0_yyyyzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xzzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xzzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyyyyy_0[i] = 2.0 * g_0_yyyyyz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_yyyyyyy_0[i] * pb_z + g_0_yyyyyzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_yyyyyyz_0[i] = 4.0 * g_0_yyyzzz_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyyyyz_0[i] * pb_y + g_0_yyyyzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyyyzz_0[i] = 4.0 * g_0_yyyzzz_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyyyzz_0[i] * pb_y + g_0_yyyyzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyyzzz_0[i] = 4.0 * g_0_yyyzzz_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyyzzz_0[i] * pb_y + g_0_yyyyzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyzzzz_0[i] = 4.0 * g_0_yyyzzz_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyzzzz_0[i] * pb_y + g_0_yyyyzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yzzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yzzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_zzzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zzzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1440-1476 components of targeted buffer : SLSK

    auto g_0_yyyyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1440);

    auto g_0_yyyyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1441);

    auto g_0_yyyyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1442);

    auto g_0_yyyyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1443);

    auto g_0_yyyyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1444);

    auto g_0_yyyyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1445);

    auto g_0_yyyyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1446);

    auto g_0_yyyyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1447);

    auto g_0_yyyyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1448);

    auto g_0_yyyyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1449);

    auto g_0_yyyyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1450);

    auto g_0_yyyyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1451);

    auto g_0_yyyyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1452);

    auto g_0_yyyyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1453);

    auto g_0_yyyyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1454);

    auto g_0_yyyyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1455);

    auto g_0_yyyyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1456);

    auto g_0_yyyyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1457);

    auto g_0_yyyyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1458);

    auto g_0_yyyyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1459);

    auto g_0_yyyyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1460);

    auto g_0_yyyyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1461);

    auto g_0_yyyyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1462);

    auto g_0_yyyyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1463);

    auto g_0_yyyyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1464);

    auto g_0_yyyyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1465);

    auto g_0_yyyyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1466);

    auto g_0_yyyyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1467);

    auto g_0_yyyyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1468);

    auto g_0_yyyyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1469);

    auto g_0_yyyyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1470);

    auto g_0_yyyyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1471);

    auto g_0_yyyyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1472);

    auto g_0_yyyyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1473);

    auto g_0_yyyyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1474);

    auto g_0_yyyyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1475);

    #pragma omp simd aligned(g_0_yyyyzz_0_xxxxxxy_0, g_0_yyyyzz_0_xxxxxxy_1, g_0_yyyyzz_0_xxxxxyy_0, g_0_yyyyzz_0_xxxxxyy_1, g_0_yyyyzz_0_xxxxyyy_0, g_0_yyyyzz_0_xxxxyyy_1, g_0_yyyyzz_0_xxxyyyy_0, g_0_yyyyzz_0_xxxyyyy_1, g_0_yyyyzz_0_xxyyyyy_0, g_0_yyyyzz_0_xxyyyyy_1, g_0_yyyyzz_0_xyyyyyy_0, g_0_yyyyzz_0_xyyyyyy_1, g_0_yyyyzz_0_yyyyyyy_0, g_0_yyyyzz_0_yyyyyyy_1, g_0_yyyyzzz_0_xxxxxxy_0, g_0_yyyyzzz_0_xxxxxxy_1, g_0_yyyyzzz_0_xxxxxyy_0, g_0_yyyyzzz_0_xxxxxyy_1, g_0_yyyyzzz_0_xxxxyyy_0, g_0_yyyyzzz_0_xxxxyyy_1, g_0_yyyyzzz_0_xxxyyyy_0, g_0_yyyyzzz_0_xxxyyyy_1, g_0_yyyyzzz_0_xxyyyyy_0, g_0_yyyyzzz_0_xxyyyyy_1, g_0_yyyyzzz_0_xyyyyyy_0, g_0_yyyyzzz_0_xyyyyyy_1, g_0_yyyyzzz_0_yyyyyyy_0, g_0_yyyyzzz_0_yyyyyyy_1, g_0_yyyyzzzz_0_xxxxxxx_0, g_0_yyyyzzzz_0_xxxxxxy_0, g_0_yyyyzzzz_0_xxxxxxz_0, g_0_yyyyzzzz_0_xxxxxyy_0, g_0_yyyyzzzz_0_xxxxxyz_0, g_0_yyyyzzzz_0_xxxxxzz_0, g_0_yyyyzzzz_0_xxxxyyy_0, g_0_yyyyzzzz_0_xxxxyyz_0, g_0_yyyyzzzz_0_xxxxyzz_0, g_0_yyyyzzzz_0_xxxxzzz_0, g_0_yyyyzzzz_0_xxxyyyy_0, g_0_yyyyzzzz_0_xxxyyyz_0, g_0_yyyyzzzz_0_xxxyyzz_0, g_0_yyyyzzzz_0_xxxyzzz_0, g_0_yyyyzzzz_0_xxxzzzz_0, g_0_yyyyzzzz_0_xxyyyyy_0, g_0_yyyyzzzz_0_xxyyyyz_0, g_0_yyyyzzzz_0_xxyyyzz_0, g_0_yyyyzzzz_0_xxyyzzz_0, g_0_yyyyzzzz_0_xxyzzzz_0, g_0_yyyyzzzz_0_xxzzzzz_0, g_0_yyyyzzzz_0_xyyyyyy_0, g_0_yyyyzzzz_0_xyyyyyz_0, g_0_yyyyzzzz_0_xyyyyzz_0, g_0_yyyyzzzz_0_xyyyzzz_0, g_0_yyyyzzzz_0_xyyzzzz_0, g_0_yyyyzzzz_0_xyzzzzz_0, g_0_yyyyzzzz_0_xzzzzzz_0, g_0_yyyyzzzz_0_yyyyyyy_0, g_0_yyyyzzzz_0_yyyyyyz_0, g_0_yyyyzzzz_0_yyyyyzz_0, g_0_yyyyzzzz_0_yyyyzzz_0, g_0_yyyyzzzz_0_yyyzzzz_0, g_0_yyyyzzzz_0_yyzzzzz_0, g_0_yyyyzzzz_0_yzzzzzz_0, g_0_yyyyzzzz_0_zzzzzzz_0, g_0_yyyzzzz_0_xxxxxxx_0, g_0_yyyzzzz_0_xxxxxxx_1, g_0_yyyzzzz_0_xxxxxxz_0, g_0_yyyzzzz_0_xxxxxxz_1, g_0_yyyzzzz_0_xxxxxyz_0, g_0_yyyzzzz_0_xxxxxyz_1, g_0_yyyzzzz_0_xxxxxz_1, g_0_yyyzzzz_0_xxxxxzz_0, g_0_yyyzzzz_0_xxxxxzz_1, g_0_yyyzzzz_0_xxxxyyz_0, g_0_yyyzzzz_0_xxxxyyz_1, g_0_yyyzzzz_0_xxxxyz_1, g_0_yyyzzzz_0_xxxxyzz_0, g_0_yyyzzzz_0_xxxxyzz_1, g_0_yyyzzzz_0_xxxxzz_1, g_0_yyyzzzz_0_xxxxzzz_0, g_0_yyyzzzz_0_xxxxzzz_1, g_0_yyyzzzz_0_xxxyyyz_0, g_0_yyyzzzz_0_xxxyyyz_1, g_0_yyyzzzz_0_xxxyyz_1, g_0_yyyzzzz_0_xxxyyzz_0, g_0_yyyzzzz_0_xxxyyzz_1, g_0_yyyzzzz_0_xxxyzz_1, g_0_yyyzzzz_0_xxxyzzz_0, g_0_yyyzzzz_0_xxxyzzz_1, g_0_yyyzzzz_0_xxxzzz_1, g_0_yyyzzzz_0_xxxzzzz_0, g_0_yyyzzzz_0_xxxzzzz_1, g_0_yyyzzzz_0_xxyyyyz_0, g_0_yyyzzzz_0_xxyyyyz_1, g_0_yyyzzzz_0_xxyyyz_1, g_0_yyyzzzz_0_xxyyyzz_0, g_0_yyyzzzz_0_xxyyyzz_1, g_0_yyyzzzz_0_xxyyzz_1, g_0_yyyzzzz_0_xxyyzzz_0, g_0_yyyzzzz_0_xxyyzzz_1, g_0_yyyzzzz_0_xxyzzz_1, g_0_yyyzzzz_0_xxyzzzz_0, g_0_yyyzzzz_0_xxyzzzz_1, g_0_yyyzzzz_0_xxzzzz_1, g_0_yyyzzzz_0_xxzzzzz_0, g_0_yyyzzzz_0_xxzzzzz_1, g_0_yyyzzzz_0_xyyyyyz_0, g_0_yyyzzzz_0_xyyyyyz_1, g_0_yyyzzzz_0_xyyyyz_1, g_0_yyyzzzz_0_xyyyyzz_0, g_0_yyyzzzz_0_xyyyyzz_1, g_0_yyyzzzz_0_xyyyzz_1, g_0_yyyzzzz_0_xyyyzzz_0, g_0_yyyzzzz_0_xyyyzzz_1, g_0_yyyzzzz_0_xyyzzz_1, g_0_yyyzzzz_0_xyyzzzz_0, g_0_yyyzzzz_0_xyyzzzz_1, g_0_yyyzzzz_0_xyzzzz_1, g_0_yyyzzzz_0_xyzzzzz_0, g_0_yyyzzzz_0_xyzzzzz_1, g_0_yyyzzzz_0_xzzzzz_1, g_0_yyyzzzz_0_xzzzzzz_0, g_0_yyyzzzz_0_xzzzzzz_1, g_0_yyyzzzz_0_yyyyyyz_0, g_0_yyyzzzz_0_yyyyyyz_1, g_0_yyyzzzz_0_yyyyyz_1, g_0_yyyzzzz_0_yyyyyzz_0, g_0_yyyzzzz_0_yyyyyzz_1, g_0_yyyzzzz_0_yyyyzz_1, g_0_yyyzzzz_0_yyyyzzz_0, g_0_yyyzzzz_0_yyyyzzz_1, g_0_yyyzzzz_0_yyyzzz_1, g_0_yyyzzzz_0_yyyzzzz_0, g_0_yyyzzzz_0_yyyzzzz_1, g_0_yyyzzzz_0_yyzzzz_1, g_0_yyyzzzz_0_yyzzzzz_0, g_0_yyyzzzz_0_yyzzzzz_1, g_0_yyyzzzz_0_yzzzzz_1, g_0_yyyzzzz_0_yzzzzzz_0, g_0_yyyzzzz_0_yzzzzzz_1, g_0_yyyzzzz_0_zzzzzz_1, g_0_yyyzzzz_0_zzzzzzz_0, g_0_yyyzzzz_0_zzzzzzz_1, g_0_yyzzzz_0_xxxxxxx_0, g_0_yyzzzz_0_xxxxxxx_1, g_0_yyzzzz_0_xxxxxxz_0, g_0_yyzzzz_0_xxxxxxz_1, g_0_yyzzzz_0_xxxxxyz_0, g_0_yyzzzz_0_xxxxxyz_1, g_0_yyzzzz_0_xxxxxzz_0, g_0_yyzzzz_0_xxxxxzz_1, g_0_yyzzzz_0_xxxxyyz_0, g_0_yyzzzz_0_xxxxyyz_1, g_0_yyzzzz_0_xxxxyzz_0, g_0_yyzzzz_0_xxxxyzz_1, g_0_yyzzzz_0_xxxxzzz_0, g_0_yyzzzz_0_xxxxzzz_1, g_0_yyzzzz_0_xxxyyyz_0, g_0_yyzzzz_0_xxxyyyz_1, g_0_yyzzzz_0_xxxyyzz_0, g_0_yyzzzz_0_xxxyyzz_1, g_0_yyzzzz_0_xxxyzzz_0, g_0_yyzzzz_0_xxxyzzz_1, g_0_yyzzzz_0_xxxzzzz_0, g_0_yyzzzz_0_xxxzzzz_1, g_0_yyzzzz_0_xxyyyyz_0, g_0_yyzzzz_0_xxyyyyz_1, g_0_yyzzzz_0_xxyyyzz_0, g_0_yyzzzz_0_xxyyyzz_1, g_0_yyzzzz_0_xxyyzzz_0, g_0_yyzzzz_0_xxyyzzz_1, g_0_yyzzzz_0_xxyzzzz_0, g_0_yyzzzz_0_xxyzzzz_1, g_0_yyzzzz_0_xxzzzzz_0, g_0_yyzzzz_0_xxzzzzz_1, g_0_yyzzzz_0_xyyyyyz_0, g_0_yyzzzz_0_xyyyyyz_1, g_0_yyzzzz_0_xyyyyzz_0, g_0_yyzzzz_0_xyyyyzz_1, g_0_yyzzzz_0_xyyyzzz_0, g_0_yyzzzz_0_xyyyzzz_1, g_0_yyzzzz_0_xyyzzzz_0, g_0_yyzzzz_0_xyyzzzz_1, g_0_yyzzzz_0_xyzzzzz_0, g_0_yyzzzz_0_xyzzzzz_1, g_0_yyzzzz_0_xzzzzzz_0, g_0_yyzzzz_0_xzzzzzz_1, g_0_yyzzzz_0_yyyyyyz_0, g_0_yyzzzz_0_yyyyyyz_1, g_0_yyzzzz_0_yyyyyzz_0, g_0_yyzzzz_0_yyyyyzz_1, g_0_yyzzzz_0_yyyyzzz_0, g_0_yyzzzz_0_yyyyzzz_1, g_0_yyzzzz_0_yyyzzzz_0, g_0_yyzzzz_0_yyyzzzz_1, g_0_yyzzzz_0_yyzzzzz_0, g_0_yyzzzz_0_yyzzzzz_1, g_0_yyzzzz_0_yzzzzzz_0, g_0_yyzzzz_0_yzzzzzz_1, g_0_yyzzzz_0_zzzzzzz_0, g_0_yyzzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzzz_0_xxxxxxx_0[i] = 3.0 * g_0_yyzzzz_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxxxx_0[i] * pb_y + g_0_yyyzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxxxy_0[i] = 3.0 * g_0_yyyyzz_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxxxy_0[i] * pb_z + g_0_yyyyzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxxxxz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxxxz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxxyy_0[i] = 3.0 * g_0_yyyyzz_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxxyy_0[i] * pb_z + g_0_yyyyzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxxxyz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxyz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxxzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxxzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_yyyyzz_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxyyy_0[i] * pb_z + g_0_yyyyzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxxyyz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxyyz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxyzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxyzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyyyzz_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxyyyy_0[i] * pb_z + g_0_yyyyzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyzzzz_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyyyz_0[i] * pb_y + g_0_yyyzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyyzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyyyyy_0[i] = 3.0 * g_0_yyyyzz_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxyyyyy_0[i] * pb_z + g_0_yyyyzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxyyyyz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyyyz_0[i] * pb_y + g_0_yyyzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyyzz_0[i] * pb_y + g_0_yyyzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyyzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyyyyy_0[i] = 3.0 * g_0_yyyyzz_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xyyyyyy_0[i] * pb_z + g_0_yyyyzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xyyyyyz_0[i] = 3.0 * g_0_yyzzzz_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyyyz_0[i] * pb_y + g_0_yyyzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyyyzz_0[i] = 3.0 * g_0_yyzzzz_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyyzz_0[i] * pb_y + g_0_yyyzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_yyzzzz_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyzzz_0[i] * pb_y + g_0_yyyzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xzzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xzzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyyyyy_0[i] = 3.0 * g_0_yyyyzz_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_yyyyyyy_0[i] * pb_z + g_0_yyyyzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_yyyyyyz_0[i] = 3.0 * g_0_yyzzzz_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyyyyz_0[i] * pb_y + g_0_yyyzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyyyzz_0[i] = 3.0 * g_0_yyzzzz_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyyyzz_0[i] * pb_y + g_0_yyyzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyyzzz_0[i] = 3.0 * g_0_yyzzzz_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyyzzz_0[i] * pb_y + g_0_yyyzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_yyzzzz_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyzzzz_0[i] * pb_y + g_0_yyyzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yzzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yzzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_zzzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zzzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1476-1512 components of targeted buffer : SLSK

    auto g_0_yyyzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1476);

    auto g_0_yyyzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1477);

    auto g_0_yyyzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1478);

    auto g_0_yyyzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1479);

    auto g_0_yyyzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1480);

    auto g_0_yyyzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1481);

    auto g_0_yyyzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1482);

    auto g_0_yyyzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1483);

    auto g_0_yyyzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1484);

    auto g_0_yyyzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1485);

    auto g_0_yyyzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1486);

    auto g_0_yyyzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1487);

    auto g_0_yyyzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1488);

    auto g_0_yyyzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1489);

    auto g_0_yyyzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1490);

    auto g_0_yyyzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1491);

    auto g_0_yyyzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1492);

    auto g_0_yyyzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1493);

    auto g_0_yyyzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1494);

    auto g_0_yyyzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1495);

    auto g_0_yyyzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1496);

    auto g_0_yyyzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1497);

    auto g_0_yyyzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1498);

    auto g_0_yyyzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1499);

    auto g_0_yyyzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1500);

    auto g_0_yyyzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1501);

    auto g_0_yyyzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1502);

    auto g_0_yyyzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1503);

    auto g_0_yyyzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1504);

    auto g_0_yyyzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1505);

    auto g_0_yyyzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1506);

    auto g_0_yyyzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1507);

    auto g_0_yyyzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1508);

    auto g_0_yyyzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1509);

    auto g_0_yyyzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1510);

    auto g_0_yyyzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1511);

    #pragma omp simd aligned(g_0_yyyzzz_0_xxxxxxy_0, g_0_yyyzzz_0_xxxxxxy_1, g_0_yyyzzz_0_xxxxxyy_0, g_0_yyyzzz_0_xxxxxyy_1, g_0_yyyzzz_0_xxxxyyy_0, g_0_yyyzzz_0_xxxxyyy_1, g_0_yyyzzz_0_xxxyyyy_0, g_0_yyyzzz_0_xxxyyyy_1, g_0_yyyzzz_0_xxyyyyy_0, g_0_yyyzzz_0_xxyyyyy_1, g_0_yyyzzz_0_xyyyyyy_0, g_0_yyyzzz_0_xyyyyyy_1, g_0_yyyzzz_0_yyyyyyy_0, g_0_yyyzzz_0_yyyyyyy_1, g_0_yyyzzzz_0_xxxxxxy_0, g_0_yyyzzzz_0_xxxxxxy_1, g_0_yyyzzzz_0_xxxxxyy_0, g_0_yyyzzzz_0_xxxxxyy_1, g_0_yyyzzzz_0_xxxxyyy_0, g_0_yyyzzzz_0_xxxxyyy_1, g_0_yyyzzzz_0_xxxyyyy_0, g_0_yyyzzzz_0_xxxyyyy_1, g_0_yyyzzzz_0_xxyyyyy_0, g_0_yyyzzzz_0_xxyyyyy_1, g_0_yyyzzzz_0_xyyyyyy_0, g_0_yyyzzzz_0_xyyyyyy_1, g_0_yyyzzzz_0_yyyyyyy_0, g_0_yyyzzzz_0_yyyyyyy_1, g_0_yyyzzzzz_0_xxxxxxx_0, g_0_yyyzzzzz_0_xxxxxxy_0, g_0_yyyzzzzz_0_xxxxxxz_0, g_0_yyyzzzzz_0_xxxxxyy_0, g_0_yyyzzzzz_0_xxxxxyz_0, g_0_yyyzzzzz_0_xxxxxzz_0, g_0_yyyzzzzz_0_xxxxyyy_0, g_0_yyyzzzzz_0_xxxxyyz_0, g_0_yyyzzzzz_0_xxxxyzz_0, g_0_yyyzzzzz_0_xxxxzzz_0, g_0_yyyzzzzz_0_xxxyyyy_0, g_0_yyyzzzzz_0_xxxyyyz_0, g_0_yyyzzzzz_0_xxxyyzz_0, g_0_yyyzzzzz_0_xxxyzzz_0, g_0_yyyzzzzz_0_xxxzzzz_0, g_0_yyyzzzzz_0_xxyyyyy_0, g_0_yyyzzzzz_0_xxyyyyz_0, g_0_yyyzzzzz_0_xxyyyzz_0, g_0_yyyzzzzz_0_xxyyzzz_0, g_0_yyyzzzzz_0_xxyzzzz_0, g_0_yyyzzzzz_0_xxzzzzz_0, g_0_yyyzzzzz_0_xyyyyyy_0, g_0_yyyzzzzz_0_xyyyyyz_0, g_0_yyyzzzzz_0_xyyyyzz_0, g_0_yyyzzzzz_0_xyyyzzz_0, g_0_yyyzzzzz_0_xyyzzzz_0, g_0_yyyzzzzz_0_xyzzzzz_0, g_0_yyyzzzzz_0_xzzzzzz_0, g_0_yyyzzzzz_0_yyyyyyy_0, g_0_yyyzzzzz_0_yyyyyyz_0, g_0_yyyzzzzz_0_yyyyyzz_0, g_0_yyyzzzzz_0_yyyyzzz_0, g_0_yyyzzzzz_0_yyyzzzz_0, g_0_yyyzzzzz_0_yyzzzzz_0, g_0_yyyzzzzz_0_yzzzzzz_0, g_0_yyyzzzzz_0_zzzzzzz_0, g_0_yyzzzzz_0_xxxxxxx_0, g_0_yyzzzzz_0_xxxxxxx_1, g_0_yyzzzzz_0_xxxxxxz_0, g_0_yyzzzzz_0_xxxxxxz_1, g_0_yyzzzzz_0_xxxxxyz_0, g_0_yyzzzzz_0_xxxxxyz_1, g_0_yyzzzzz_0_xxxxxz_1, g_0_yyzzzzz_0_xxxxxzz_0, g_0_yyzzzzz_0_xxxxxzz_1, g_0_yyzzzzz_0_xxxxyyz_0, g_0_yyzzzzz_0_xxxxyyz_1, g_0_yyzzzzz_0_xxxxyz_1, g_0_yyzzzzz_0_xxxxyzz_0, g_0_yyzzzzz_0_xxxxyzz_1, g_0_yyzzzzz_0_xxxxzz_1, g_0_yyzzzzz_0_xxxxzzz_0, g_0_yyzzzzz_0_xxxxzzz_1, g_0_yyzzzzz_0_xxxyyyz_0, g_0_yyzzzzz_0_xxxyyyz_1, g_0_yyzzzzz_0_xxxyyz_1, g_0_yyzzzzz_0_xxxyyzz_0, g_0_yyzzzzz_0_xxxyyzz_1, g_0_yyzzzzz_0_xxxyzz_1, g_0_yyzzzzz_0_xxxyzzz_0, g_0_yyzzzzz_0_xxxyzzz_1, g_0_yyzzzzz_0_xxxzzz_1, g_0_yyzzzzz_0_xxxzzzz_0, g_0_yyzzzzz_0_xxxzzzz_1, g_0_yyzzzzz_0_xxyyyyz_0, g_0_yyzzzzz_0_xxyyyyz_1, g_0_yyzzzzz_0_xxyyyz_1, g_0_yyzzzzz_0_xxyyyzz_0, g_0_yyzzzzz_0_xxyyyzz_1, g_0_yyzzzzz_0_xxyyzz_1, g_0_yyzzzzz_0_xxyyzzz_0, g_0_yyzzzzz_0_xxyyzzz_1, g_0_yyzzzzz_0_xxyzzz_1, g_0_yyzzzzz_0_xxyzzzz_0, g_0_yyzzzzz_0_xxyzzzz_1, g_0_yyzzzzz_0_xxzzzz_1, g_0_yyzzzzz_0_xxzzzzz_0, g_0_yyzzzzz_0_xxzzzzz_1, g_0_yyzzzzz_0_xyyyyyz_0, g_0_yyzzzzz_0_xyyyyyz_1, g_0_yyzzzzz_0_xyyyyz_1, g_0_yyzzzzz_0_xyyyyzz_0, g_0_yyzzzzz_0_xyyyyzz_1, g_0_yyzzzzz_0_xyyyzz_1, g_0_yyzzzzz_0_xyyyzzz_0, g_0_yyzzzzz_0_xyyyzzz_1, g_0_yyzzzzz_0_xyyzzz_1, g_0_yyzzzzz_0_xyyzzzz_0, g_0_yyzzzzz_0_xyyzzzz_1, g_0_yyzzzzz_0_xyzzzz_1, g_0_yyzzzzz_0_xyzzzzz_0, g_0_yyzzzzz_0_xyzzzzz_1, g_0_yyzzzzz_0_xzzzzz_1, g_0_yyzzzzz_0_xzzzzzz_0, g_0_yyzzzzz_0_xzzzzzz_1, g_0_yyzzzzz_0_yyyyyyz_0, g_0_yyzzzzz_0_yyyyyyz_1, g_0_yyzzzzz_0_yyyyyz_1, g_0_yyzzzzz_0_yyyyyzz_0, g_0_yyzzzzz_0_yyyyyzz_1, g_0_yyzzzzz_0_yyyyzz_1, g_0_yyzzzzz_0_yyyyzzz_0, g_0_yyzzzzz_0_yyyyzzz_1, g_0_yyzzzzz_0_yyyzzz_1, g_0_yyzzzzz_0_yyyzzzz_0, g_0_yyzzzzz_0_yyyzzzz_1, g_0_yyzzzzz_0_yyzzzz_1, g_0_yyzzzzz_0_yyzzzzz_0, g_0_yyzzzzz_0_yyzzzzz_1, g_0_yyzzzzz_0_yzzzzz_1, g_0_yyzzzzz_0_yzzzzzz_0, g_0_yyzzzzz_0_yzzzzzz_1, g_0_yyzzzzz_0_zzzzzz_1, g_0_yyzzzzz_0_zzzzzzz_0, g_0_yyzzzzz_0_zzzzzzz_1, g_0_yzzzzz_0_xxxxxxx_0, g_0_yzzzzz_0_xxxxxxx_1, g_0_yzzzzz_0_xxxxxxz_0, g_0_yzzzzz_0_xxxxxxz_1, g_0_yzzzzz_0_xxxxxyz_0, g_0_yzzzzz_0_xxxxxyz_1, g_0_yzzzzz_0_xxxxxzz_0, g_0_yzzzzz_0_xxxxxzz_1, g_0_yzzzzz_0_xxxxyyz_0, g_0_yzzzzz_0_xxxxyyz_1, g_0_yzzzzz_0_xxxxyzz_0, g_0_yzzzzz_0_xxxxyzz_1, g_0_yzzzzz_0_xxxxzzz_0, g_0_yzzzzz_0_xxxxzzz_1, g_0_yzzzzz_0_xxxyyyz_0, g_0_yzzzzz_0_xxxyyyz_1, g_0_yzzzzz_0_xxxyyzz_0, g_0_yzzzzz_0_xxxyyzz_1, g_0_yzzzzz_0_xxxyzzz_0, g_0_yzzzzz_0_xxxyzzz_1, g_0_yzzzzz_0_xxxzzzz_0, g_0_yzzzzz_0_xxxzzzz_1, g_0_yzzzzz_0_xxyyyyz_0, g_0_yzzzzz_0_xxyyyyz_1, g_0_yzzzzz_0_xxyyyzz_0, g_0_yzzzzz_0_xxyyyzz_1, g_0_yzzzzz_0_xxyyzzz_0, g_0_yzzzzz_0_xxyyzzz_1, g_0_yzzzzz_0_xxyzzzz_0, g_0_yzzzzz_0_xxyzzzz_1, g_0_yzzzzz_0_xxzzzzz_0, g_0_yzzzzz_0_xxzzzzz_1, g_0_yzzzzz_0_xyyyyyz_0, g_0_yzzzzz_0_xyyyyyz_1, g_0_yzzzzz_0_xyyyyzz_0, g_0_yzzzzz_0_xyyyyzz_1, g_0_yzzzzz_0_xyyyzzz_0, g_0_yzzzzz_0_xyyyzzz_1, g_0_yzzzzz_0_xyyzzzz_0, g_0_yzzzzz_0_xyyzzzz_1, g_0_yzzzzz_0_xyzzzzz_0, g_0_yzzzzz_0_xyzzzzz_1, g_0_yzzzzz_0_xzzzzzz_0, g_0_yzzzzz_0_xzzzzzz_1, g_0_yzzzzz_0_yyyyyyz_0, g_0_yzzzzz_0_yyyyyyz_1, g_0_yzzzzz_0_yyyyyzz_0, g_0_yzzzzz_0_yyyyyzz_1, g_0_yzzzzz_0_yyyyzzz_0, g_0_yzzzzz_0_yyyyzzz_1, g_0_yzzzzz_0_yyyzzzz_0, g_0_yzzzzz_0_yyyzzzz_1, g_0_yzzzzz_0_yyzzzzz_0, g_0_yzzzzz_0_yyzzzzz_1, g_0_yzzzzz_0_yzzzzzz_0, g_0_yzzzzz_0_yzzzzzz_1, g_0_yzzzzz_0_zzzzzzz_0, g_0_yzzzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzzz_0_xxxxxxx_0[i] = 2.0 * g_0_yzzzzz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxxxx_0[i] * pb_y + g_0_yyzzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxxxy_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxxxy_0[i] * pb_z + g_0_yyyzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxxxxz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxxxz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxxyy_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxxyy_0[i] * pb_z + g_0_yyyzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxxxyz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxyz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxxzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxxzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyyzzz_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxyyy_0[i] * pb_z + g_0_yyyzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxyyz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxyzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxyzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_yyyzzz_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxyyyy_0[i] * pb_z + g_0_yyyzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxyyyz_0[i] = 2.0 * g_0_yzzzzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyyyz_0[i] * pb_y + g_0_yyzzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyyzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxyzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyyyyy_0[i] = 4.0 * g_0_yyyzzz_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxyyyyy_0[i] * pb_z + g_0_yyyzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzzzzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyyyz_0[i] * pb_y + g_0_yyzzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzzzzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyyzz_0[i] * pb_y + g_0_yyzzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyyyyy_0[i] = 4.0 * g_0_yyyzzz_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xyyyyyy_0[i] * pb_z + g_0_yyyzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xyyyyyz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyyyz_0[i] * pb_y + g_0_yyzzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyyyzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyyzz_0[i] * pb_y + g_0_yyzzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyyzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyzzz_0[i] * pb_y + g_0_yyzzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xzzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xzzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyyyyy_0[i] = 4.0 * g_0_yyyzzz_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_yyyyyyy_0[i] * pb_z + g_0_yyyzzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_yyyyyyz_0[i] = 2.0 * g_0_yzzzzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyyyyz_0[i] * pb_y + g_0_yyzzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyyyzz_0[i] = 2.0 * g_0_yzzzzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyyyzz_0[i] * pb_y + g_0_yyzzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyyzzz_0[i] = 2.0 * g_0_yzzzzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyyzzz_0[i] * pb_y + g_0_yyzzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyzzzz_0[i] = 2.0 * g_0_yzzzzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyzzzz_0[i] * pb_y + g_0_yyzzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yzzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yzzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_zzzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zzzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1512-1548 components of targeted buffer : SLSK

    auto g_0_yyzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1512);

    auto g_0_yyzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1513);

    auto g_0_yyzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1514);

    auto g_0_yyzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1515);

    auto g_0_yyzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1516);

    auto g_0_yyzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1517);

    auto g_0_yyzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1518);

    auto g_0_yyzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1519);

    auto g_0_yyzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1520);

    auto g_0_yyzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1521);

    auto g_0_yyzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1522);

    auto g_0_yyzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1523);

    auto g_0_yyzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1524);

    auto g_0_yyzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1525);

    auto g_0_yyzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1526);

    auto g_0_yyzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1527);

    auto g_0_yyzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1528);

    auto g_0_yyzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1529);

    auto g_0_yyzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1530);

    auto g_0_yyzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1531);

    auto g_0_yyzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1532);

    auto g_0_yyzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1533);

    auto g_0_yyzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1534);

    auto g_0_yyzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1535);

    auto g_0_yyzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1536);

    auto g_0_yyzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1537);

    auto g_0_yyzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1538);

    auto g_0_yyzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1539);

    auto g_0_yyzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1540);

    auto g_0_yyzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1541);

    auto g_0_yyzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1542);

    auto g_0_yyzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1543);

    auto g_0_yyzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1544);

    auto g_0_yyzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1545);

    auto g_0_yyzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1546);

    auto g_0_yyzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1547);

    #pragma omp simd aligned(g_0_yyzzzz_0_xxxxxxy_0, g_0_yyzzzz_0_xxxxxxy_1, g_0_yyzzzz_0_xxxxxyy_0, g_0_yyzzzz_0_xxxxxyy_1, g_0_yyzzzz_0_xxxxyyy_0, g_0_yyzzzz_0_xxxxyyy_1, g_0_yyzzzz_0_xxxyyyy_0, g_0_yyzzzz_0_xxxyyyy_1, g_0_yyzzzz_0_xxyyyyy_0, g_0_yyzzzz_0_xxyyyyy_1, g_0_yyzzzz_0_xyyyyyy_0, g_0_yyzzzz_0_xyyyyyy_1, g_0_yyzzzz_0_yyyyyyy_0, g_0_yyzzzz_0_yyyyyyy_1, g_0_yyzzzzz_0_xxxxxxy_0, g_0_yyzzzzz_0_xxxxxxy_1, g_0_yyzzzzz_0_xxxxxyy_0, g_0_yyzzzzz_0_xxxxxyy_1, g_0_yyzzzzz_0_xxxxyyy_0, g_0_yyzzzzz_0_xxxxyyy_1, g_0_yyzzzzz_0_xxxyyyy_0, g_0_yyzzzzz_0_xxxyyyy_1, g_0_yyzzzzz_0_xxyyyyy_0, g_0_yyzzzzz_0_xxyyyyy_1, g_0_yyzzzzz_0_xyyyyyy_0, g_0_yyzzzzz_0_xyyyyyy_1, g_0_yyzzzzz_0_yyyyyyy_0, g_0_yyzzzzz_0_yyyyyyy_1, g_0_yyzzzzzz_0_xxxxxxx_0, g_0_yyzzzzzz_0_xxxxxxy_0, g_0_yyzzzzzz_0_xxxxxxz_0, g_0_yyzzzzzz_0_xxxxxyy_0, g_0_yyzzzzzz_0_xxxxxyz_0, g_0_yyzzzzzz_0_xxxxxzz_0, g_0_yyzzzzzz_0_xxxxyyy_0, g_0_yyzzzzzz_0_xxxxyyz_0, g_0_yyzzzzzz_0_xxxxyzz_0, g_0_yyzzzzzz_0_xxxxzzz_0, g_0_yyzzzzzz_0_xxxyyyy_0, g_0_yyzzzzzz_0_xxxyyyz_0, g_0_yyzzzzzz_0_xxxyyzz_0, g_0_yyzzzzzz_0_xxxyzzz_0, g_0_yyzzzzzz_0_xxxzzzz_0, g_0_yyzzzzzz_0_xxyyyyy_0, g_0_yyzzzzzz_0_xxyyyyz_0, g_0_yyzzzzzz_0_xxyyyzz_0, g_0_yyzzzzzz_0_xxyyzzz_0, g_0_yyzzzzzz_0_xxyzzzz_0, g_0_yyzzzzzz_0_xxzzzzz_0, g_0_yyzzzzzz_0_xyyyyyy_0, g_0_yyzzzzzz_0_xyyyyyz_0, g_0_yyzzzzzz_0_xyyyyzz_0, g_0_yyzzzzzz_0_xyyyzzz_0, g_0_yyzzzzzz_0_xyyzzzz_0, g_0_yyzzzzzz_0_xyzzzzz_0, g_0_yyzzzzzz_0_xzzzzzz_0, g_0_yyzzzzzz_0_yyyyyyy_0, g_0_yyzzzzzz_0_yyyyyyz_0, g_0_yyzzzzzz_0_yyyyyzz_0, g_0_yyzzzzzz_0_yyyyzzz_0, g_0_yyzzzzzz_0_yyyzzzz_0, g_0_yyzzzzzz_0_yyzzzzz_0, g_0_yyzzzzzz_0_yzzzzzz_0, g_0_yyzzzzzz_0_zzzzzzz_0, g_0_yzzzzzz_0_xxxxxxx_0, g_0_yzzzzzz_0_xxxxxxx_1, g_0_yzzzzzz_0_xxxxxxz_0, g_0_yzzzzzz_0_xxxxxxz_1, g_0_yzzzzzz_0_xxxxxyz_0, g_0_yzzzzzz_0_xxxxxyz_1, g_0_yzzzzzz_0_xxxxxz_1, g_0_yzzzzzz_0_xxxxxzz_0, g_0_yzzzzzz_0_xxxxxzz_1, g_0_yzzzzzz_0_xxxxyyz_0, g_0_yzzzzzz_0_xxxxyyz_1, g_0_yzzzzzz_0_xxxxyz_1, g_0_yzzzzzz_0_xxxxyzz_0, g_0_yzzzzzz_0_xxxxyzz_1, g_0_yzzzzzz_0_xxxxzz_1, g_0_yzzzzzz_0_xxxxzzz_0, g_0_yzzzzzz_0_xxxxzzz_1, g_0_yzzzzzz_0_xxxyyyz_0, g_0_yzzzzzz_0_xxxyyyz_1, g_0_yzzzzzz_0_xxxyyz_1, g_0_yzzzzzz_0_xxxyyzz_0, g_0_yzzzzzz_0_xxxyyzz_1, g_0_yzzzzzz_0_xxxyzz_1, g_0_yzzzzzz_0_xxxyzzz_0, g_0_yzzzzzz_0_xxxyzzz_1, g_0_yzzzzzz_0_xxxzzz_1, g_0_yzzzzzz_0_xxxzzzz_0, g_0_yzzzzzz_0_xxxzzzz_1, g_0_yzzzzzz_0_xxyyyyz_0, g_0_yzzzzzz_0_xxyyyyz_1, g_0_yzzzzzz_0_xxyyyz_1, g_0_yzzzzzz_0_xxyyyzz_0, g_0_yzzzzzz_0_xxyyyzz_1, g_0_yzzzzzz_0_xxyyzz_1, g_0_yzzzzzz_0_xxyyzzz_0, g_0_yzzzzzz_0_xxyyzzz_1, g_0_yzzzzzz_0_xxyzzz_1, g_0_yzzzzzz_0_xxyzzzz_0, g_0_yzzzzzz_0_xxyzzzz_1, g_0_yzzzzzz_0_xxzzzz_1, g_0_yzzzzzz_0_xxzzzzz_0, g_0_yzzzzzz_0_xxzzzzz_1, g_0_yzzzzzz_0_xyyyyyz_0, g_0_yzzzzzz_0_xyyyyyz_1, g_0_yzzzzzz_0_xyyyyz_1, g_0_yzzzzzz_0_xyyyyzz_0, g_0_yzzzzzz_0_xyyyyzz_1, g_0_yzzzzzz_0_xyyyzz_1, g_0_yzzzzzz_0_xyyyzzz_0, g_0_yzzzzzz_0_xyyyzzz_1, g_0_yzzzzzz_0_xyyzzz_1, g_0_yzzzzzz_0_xyyzzzz_0, g_0_yzzzzzz_0_xyyzzzz_1, g_0_yzzzzzz_0_xyzzzz_1, g_0_yzzzzzz_0_xyzzzzz_0, g_0_yzzzzzz_0_xyzzzzz_1, g_0_yzzzzzz_0_xzzzzz_1, g_0_yzzzzzz_0_xzzzzzz_0, g_0_yzzzzzz_0_xzzzzzz_1, g_0_yzzzzzz_0_yyyyyyz_0, g_0_yzzzzzz_0_yyyyyyz_1, g_0_yzzzzzz_0_yyyyyz_1, g_0_yzzzzzz_0_yyyyyzz_0, g_0_yzzzzzz_0_yyyyyzz_1, g_0_yzzzzzz_0_yyyyzz_1, g_0_yzzzzzz_0_yyyyzzz_0, g_0_yzzzzzz_0_yyyyzzz_1, g_0_yzzzzzz_0_yyyzzz_1, g_0_yzzzzzz_0_yyyzzzz_0, g_0_yzzzzzz_0_yyyzzzz_1, g_0_yzzzzzz_0_yyzzzz_1, g_0_yzzzzzz_0_yyzzzzz_0, g_0_yzzzzzz_0_yyzzzzz_1, g_0_yzzzzzz_0_yzzzzz_1, g_0_yzzzzzz_0_yzzzzzz_0, g_0_yzzzzzz_0_yzzzzzz_1, g_0_yzzzzzz_0_zzzzzz_1, g_0_yzzzzzz_0_zzzzzzz_0, g_0_yzzzzzz_0_zzzzzzz_1, g_0_zzzzzz_0_xxxxxxx_0, g_0_zzzzzz_0_xxxxxxx_1, g_0_zzzzzz_0_xxxxxxz_0, g_0_zzzzzz_0_xxxxxxz_1, g_0_zzzzzz_0_xxxxxyz_0, g_0_zzzzzz_0_xxxxxyz_1, g_0_zzzzzz_0_xxxxxzz_0, g_0_zzzzzz_0_xxxxxzz_1, g_0_zzzzzz_0_xxxxyyz_0, g_0_zzzzzz_0_xxxxyyz_1, g_0_zzzzzz_0_xxxxyzz_0, g_0_zzzzzz_0_xxxxyzz_1, g_0_zzzzzz_0_xxxxzzz_0, g_0_zzzzzz_0_xxxxzzz_1, g_0_zzzzzz_0_xxxyyyz_0, g_0_zzzzzz_0_xxxyyyz_1, g_0_zzzzzz_0_xxxyyzz_0, g_0_zzzzzz_0_xxxyyzz_1, g_0_zzzzzz_0_xxxyzzz_0, g_0_zzzzzz_0_xxxyzzz_1, g_0_zzzzzz_0_xxxzzzz_0, g_0_zzzzzz_0_xxxzzzz_1, g_0_zzzzzz_0_xxyyyyz_0, g_0_zzzzzz_0_xxyyyyz_1, g_0_zzzzzz_0_xxyyyzz_0, g_0_zzzzzz_0_xxyyyzz_1, g_0_zzzzzz_0_xxyyzzz_0, g_0_zzzzzz_0_xxyyzzz_1, g_0_zzzzzz_0_xxyzzzz_0, g_0_zzzzzz_0_xxyzzzz_1, g_0_zzzzzz_0_xxzzzzz_0, g_0_zzzzzz_0_xxzzzzz_1, g_0_zzzzzz_0_xyyyyyz_0, g_0_zzzzzz_0_xyyyyyz_1, g_0_zzzzzz_0_xyyyyzz_0, g_0_zzzzzz_0_xyyyyzz_1, g_0_zzzzzz_0_xyyyzzz_0, g_0_zzzzzz_0_xyyyzzz_1, g_0_zzzzzz_0_xyyzzzz_0, g_0_zzzzzz_0_xyyzzzz_1, g_0_zzzzzz_0_xyzzzzz_0, g_0_zzzzzz_0_xyzzzzz_1, g_0_zzzzzz_0_xzzzzzz_0, g_0_zzzzzz_0_xzzzzzz_1, g_0_zzzzzz_0_yyyyyyz_0, g_0_zzzzzz_0_yyyyyyz_1, g_0_zzzzzz_0_yyyyyzz_0, g_0_zzzzzz_0_yyyyyzz_1, g_0_zzzzzz_0_yyyyzzz_0, g_0_zzzzzz_0_yyyyzzz_1, g_0_zzzzzz_0_yyyzzzz_0, g_0_zzzzzz_0_yyyzzzz_1, g_0_zzzzzz_0_yyzzzzz_0, g_0_zzzzzz_0_yyzzzzz_1, g_0_zzzzzz_0_yzzzzzz_0, g_0_zzzzzz_0_yzzzzzz_1, g_0_zzzzzz_0_zzzzzzz_0, g_0_zzzzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzzz_0_xxxxxxx_0[i] = g_0_zzzzzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxxxx_0[i] * pb_y + g_0_yzzzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxxxy_0[i] = 5.0 * g_0_yyzzzz_0_xxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxxxy_0[i] * pb_z + g_0_yyzzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxxxxz_0[i] = g_0_zzzzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxxxz_0[i] * pb_y + g_0_yzzzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyzzzz_0_xxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxxyy_0[i] * pb_z + g_0_yyzzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxxxyz_0[i] = g_0_zzzzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxxyz_0[i] * pb_y + g_0_yzzzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxxzz_0[i] = g_0_zzzzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxxzz_0[i] * pb_y + g_0_yzzzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxyyy_0[i] = 5.0 * g_0_yyzzzz_0_xxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxyyy_0[i] * pb_z + g_0_yyzzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxxyyz_0[i] = g_0_zzzzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxyyz_0[i] * pb_y + g_0_yzzzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxyzz_0[i] = g_0_zzzzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxyzz_0[i] * pb_y + g_0_yzzzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxzzz_0[i] = g_0_zzzzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxzzz_0[i] * pb_y + g_0_yzzzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxyyyy_0[i] = 5.0 * g_0_yyzzzz_0_xxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxyyyy_0[i] * pb_z + g_0_yyzzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxyyyz_0[i] = g_0_zzzzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyyyz_0[i] * pb_y + g_0_yzzzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxyyzz_0[i] = g_0_zzzzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyyzz_0[i] * pb_y + g_0_yzzzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxyzzz_0[i] = g_0_zzzzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyzzz_0[i] * pb_y + g_0_yzzzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxzzzz_0[i] = g_0_zzzzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxzzzz_0[i] * pb_y + g_0_yzzzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_yyzzzz_0_xxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxyyyyy_0[i] * pb_z + g_0_yyzzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxyyyyz_0[i] = g_0_zzzzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyyyz_0[i] * pb_y + g_0_yzzzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyyyzz_0[i] = g_0_zzzzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyyzz_0[i] * pb_y + g_0_yzzzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyyzzz_0[i] = g_0_zzzzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyzzz_0[i] * pb_y + g_0_yzzzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyzzzz_0[i] = g_0_zzzzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyzzzz_0[i] * pb_y + g_0_yzzzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxzzzzz_0[i] = g_0_zzzzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxzzzzz_0[i] * pb_y + g_0_yzzzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyyyyy_0[i] = 5.0 * g_0_yyzzzz_0_xyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xyyyyyy_0[i] * pb_z + g_0_yyzzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xyyyyyz_0[i] = g_0_zzzzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyyyz_0[i] * pb_y + g_0_yzzzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyyyzz_0[i] = g_0_zzzzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyyzz_0[i] * pb_y + g_0_yzzzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyyzzz_0[i] = g_0_zzzzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyzzz_0[i] * pb_y + g_0_yzzzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyzzzz_0[i] = g_0_zzzzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyzzzz_0[i] * pb_y + g_0_yzzzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyzzzzz_0[i] = g_0_zzzzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyzzzzz_0[i] * pb_y + g_0_yzzzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xzzzzzz_0[i] = g_0_zzzzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzzzzzz_0[i] * pb_y + g_0_yzzzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyyyyy_0[i] = 5.0 * g_0_yyzzzz_0_yyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_yyyyyyy_0[i] * pb_z + g_0_yyzzzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_yyyyyyz_0[i] = g_0_zzzzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yzzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyyyyz_0[i] * pb_y + g_0_yzzzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyyyzz_0[i] = g_0_zzzzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yzzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyyyzz_0[i] * pb_y + g_0_yzzzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyyzzz_0[i] = g_0_zzzzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyyzzz_0[i] * pb_y + g_0_yzzzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyzzzz_0[i] = g_0_zzzzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyzzzz_0[i] * pb_y + g_0_yzzzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyzzzzz_0[i] = g_0_zzzzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyzzzzz_0[i] * pb_y + g_0_yzzzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yzzzzzz_0[i] = g_0_zzzzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yzzzzzz_0[i] * pb_y + g_0_yzzzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_zzzzzzz_0[i] = g_0_zzzzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzzzzzz_0[i] * pb_y + g_0_yzzzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1548-1584 components of targeted buffer : SLSK

    auto g_0_yzzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1548);

    auto g_0_yzzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1549);

    auto g_0_yzzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1550);

    auto g_0_yzzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1551);

    auto g_0_yzzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1552);

    auto g_0_yzzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1553);

    auto g_0_yzzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1554);

    auto g_0_yzzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1555);

    auto g_0_yzzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1556);

    auto g_0_yzzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1557);

    auto g_0_yzzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1558);

    auto g_0_yzzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1559);

    auto g_0_yzzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1560);

    auto g_0_yzzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1561);

    auto g_0_yzzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1562);

    auto g_0_yzzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1563);

    auto g_0_yzzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1564);

    auto g_0_yzzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1565);

    auto g_0_yzzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1566);

    auto g_0_yzzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1567);

    auto g_0_yzzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1568);

    auto g_0_yzzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1569);

    auto g_0_yzzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1570);

    auto g_0_yzzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1571);

    auto g_0_yzzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1572);

    auto g_0_yzzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1573);

    auto g_0_yzzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1574);

    auto g_0_yzzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1575);

    auto g_0_yzzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1576);

    auto g_0_yzzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1577);

    auto g_0_yzzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1578);

    auto g_0_yzzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1579);

    auto g_0_yzzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1580);

    auto g_0_yzzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1581);

    auto g_0_yzzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1582);

    auto g_0_yzzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1583);

    #pragma omp simd aligned(g_0_yzzzzzzz_0_xxxxxxx_0, g_0_yzzzzzzz_0_xxxxxxy_0, g_0_yzzzzzzz_0_xxxxxxz_0, g_0_yzzzzzzz_0_xxxxxyy_0, g_0_yzzzzzzz_0_xxxxxyz_0, g_0_yzzzzzzz_0_xxxxxzz_0, g_0_yzzzzzzz_0_xxxxyyy_0, g_0_yzzzzzzz_0_xxxxyyz_0, g_0_yzzzzzzz_0_xxxxyzz_0, g_0_yzzzzzzz_0_xxxxzzz_0, g_0_yzzzzzzz_0_xxxyyyy_0, g_0_yzzzzzzz_0_xxxyyyz_0, g_0_yzzzzzzz_0_xxxyyzz_0, g_0_yzzzzzzz_0_xxxyzzz_0, g_0_yzzzzzzz_0_xxxzzzz_0, g_0_yzzzzzzz_0_xxyyyyy_0, g_0_yzzzzzzz_0_xxyyyyz_0, g_0_yzzzzzzz_0_xxyyyzz_0, g_0_yzzzzzzz_0_xxyyzzz_0, g_0_yzzzzzzz_0_xxyzzzz_0, g_0_yzzzzzzz_0_xxzzzzz_0, g_0_yzzzzzzz_0_xyyyyyy_0, g_0_yzzzzzzz_0_xyyyyyz_0, g_0_yzzzzzzz_0_xyyyyzz_0, g_0_yzzzzzzz_0_xyyyzzz_0, g_0_yzzzzzzz_0_xyyzzzz_0, g_0_yzzzzzzz_0_xyzzzzz_0, g_0_yzzzzzzz_0_xzzzzzz_0, g_0_yzzzzzzz_0_yyyyyyy_0, g_0_yzzzzzzz_0_yyyyyyz_0, g_0_yzzzzzzz_0_yyyyyzz_0, g_0_yzzzzzzz_0_yyyyzzz_0, g_0_yzzzzzzz_0_yyyzzzz_0, g_0_yzzzzzzz_0_yyzzzzz_0, g_0_yzzzzzzz_0_yzzzzzz_0, g_0_yzzzzzzz_0_zzzzzzz_0, g_0_zzzzzzz_0_xxxxxx_1, g_0_zzzzzzz_0_xxxxxxx_0, g_0_zzzzzzz_0_xxxxxxx_1, g_0_zzzzzzz_0_xxxxxxy_0, g_0_zzzzzzz_0_xxxxxxy_1, g_0_zzzzzzz_0_xxxxxxz_0, g_0_zzzzzzz_0_xxxxxxz_1, g_0_zzzzzzz_0_xxxxxy_1, g_0_zzzzzzz_0_xxxxxyy_0, g_0_zzzzzzz_0_xxxxxyy_1, g_0_zzzzzzz_0_xxxxxyz_0, g_0_zzzzzzz_0_xxxxxyz_1, g_0_zzzzzzz_0_xxxxxz_1, g_0_zzzzzzz_0_xxxxxzz_0, g_0_zzzzzzz_0_xxxxxzz_1, g_0_zzzzzzz_0_xxxxyy_1, g_0_zzzzzzz_0_xxxxyyy_0, g_0_zzzzzzz_0_xxxxyyy_1, g_0_zzzzzzz_0_xxxxyyz_0, g_0_zzzzzzz_0_xxxxyyz_1, g_0_zzzzzzz_0_xxxxyz_1, g_0_zzzzzzz_0_xxxxyzz_0, g_0_zzzzzzz_0_xxxxyzz_1, g_0_zzzzzzz_0_xxxxzz_1, g_0_zzzzzzz_0_xxxxzzz_0, g_0_zzzzzzz_0_xxxxzzz_1, g_0_zzzzzzz_0_xxxyyy_1, g_0_zzzzzzz_0_xxxyyyy_0, g_0_zzzzzzz_0_xxxyyyy_1, g_0_zzzzzzz_0_xxxyyyz_0, g_0_zzzzzzz_0_xxxyyyz_1, g_0_zzzzzzz_0_xxxyyz_1, g_0_zzzzzzz_0_xxxyyzz_0, g_0_zzzzzzz_0_xxxyyzz_1, g_0_zzzzzzz_0_xxxyzz_1, g_0_zzzzzzz_0_xxxyzzz_0, g_0_zzzzzzz_0_xxxyzzz_1, g_0_zzzzzzz_0_xxxzzz_1, g_0_zzzzzzz_0_xxxzzzz_0, g_0_zzzzzzz_0_xxxzzzz_1, g_0_zzzzzzz_0_xxyyyy_1, g_0_zzzzzzz_0_xxyyyyy_0, g_0_zzzzzzz_0_xxyyyyy_1, g_0_zzzzzzz_0_xxyyyyz_0, g_0_zzzzzzz_0_xxyyyyz_1, g_0_zzzzzzz_0_xxyyyz_1, g_0_zzzzzzz_0_xxyyyzz_0, g_0_zzzzzzz_0_xxyyyzz_1, g_0_zzzzzzz_0_xxyyzz_1, g_0_zzzzzzz_0_xxyyzzz_0, g_0_zzzzzzz_0_xxyyzzz_1, g_0_zzzzzzz_0_xxyzzz_1, g_0_zzzzzzz_0_xxyzzzz_0, g_0_zzzzzzz_0_xxyzzzz_1, g_0_zzzzzzz_0_xxzzzz_1, g_0_zzzzzzz_0_xxzzzzz_0, g_0_zzzzzzz_0_xxzzzzz_1, g_0_zzzzzzz_0_xyyyyy_1, g_0_zzzzzzz_0_xyyyyyy_0, g_0_zzzzzzz_0_xyyyyyy_1, g_0_zzzzzzz_0_xyyyyyz_0, g_0_zzzzzzz_0_xyyyyyz_1, g_0_zzzzzzz_0_xyyyyz_1, g_0_zzzzzzz_0_xyyyyzz_0, g_0_zzzzzzz_0_xyyyyzz_1, g_0_zzzzzzz_0_xyyyzz_1, g_0_zzzzzzz_0_xyyyzzz_0, g_0_zzzzzzz_0_xyyyzzz_1, g_0_zzzzzzz_0_xyyzzz_1, g_0_zzzzzzz_0_xyyzzzz_0, g_0_zzzzzzz_0_xyyzzzz_1, g_0_zzzzzzz_0_xyzzzz_1, g_0_zzzzzzz_0_xyzzzzz_0, g_0_zzzzzzz_0_xyzzzzz_1, g_0_zzzzzzz_0_xzzzzz_1, g_0_zzzzzzz_0_xzzzzzz_0, g_0_zzzzzzz_0_xzzzzzz_1, g_0_zzzzzzz_0_yyyyyy_1, g_0_zzzzzzz_0_yyyyyyy_0, g_0_zzzzzzz_0_yyyyyyy_1, g_0_zzzzzzz_0_yyyyyyz_0, g_0_zzzzzzz_0_yyyyyyz_1, g_0_zzzzzzz_0_yyyyyz_1, g_0_zzzzzzz_0_yyyyyzz_0, g_0_zzzzzzz_0_yyyyyzz_1, g_0_zzzzzzz_0_yyyyzz_1, g_0_zzzzzzz_0_yyyyzzz_0, g_0_zzzzzzz_0_yyyyzzz_1, g_0_zzzzzzz_0_yyyzzz_1, g_0_zzzzzzz_0_yyyzzzz_0, g_0_zzzzzzz_0_yyyzzzz_1, g_0_zzzzzzz_0_yyzzzz_1, g_0_zzzzzzz_0_yyzzzzz_0, g_0_zzzzzzz_0_yyzzzzz_1, g_0_zzzzzzz_0_yzzzzz_1, g_0_zzzzzzz_0_yzzzzzz_0, g_0_zzzzzzz_0_yzzzzzz_1, g_0_zzzzzzz_0_zzzzzz_1, g_0_zzzzzzz_0_zzzzzzz_0, g_0_zzzzzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzzz_0_xxxxxxx_0[i] = g_0_zzzzzzz_0_xxxxxxx_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxxxy_0[i] = g_0_zzzzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxxy_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxxxz_0[i] = g_0_zzzzzzz_0_xxxxxxz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxxyy_0[i] = 2.0 * g_0_zzzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxyy_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxxyz_0[i] = g_0_zzzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxyz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxxzz_0[i] = g_0_zzzzzzz_0_xxxxxzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_zzzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyyy_0[i] * pb_y + g_0_zzzzzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_zzzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyyz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxyzz_0[i] = g_0_zzzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxzzz_0[i] = g_0_zzzzzzz_0_xxxxzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_zzzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyyy_0[i] * pb_y + g_0_zzzzzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyyz_0[i] * pb_y + g_0_zzzzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_zzzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyzzz_0[i] = g_0_zzzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxzzzz_0[i] = g_0_zzzzzzz_0_xxxzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_zzzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyyy_0[i] * pb_y + g_0_zzzzzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyyyz_0[i] = 4.0 * g_0_zzzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyyz_0[i] * pb_y + g_0_zzzzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_zzzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyzz_0[i] * pb_y + g_0_zzzzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyzzzz_0[i] = g_0_zzzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxzzzzz_0[i] = g_0_zzzzzzz_0_xxzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyyyy_0[i] = 6.0 * g_0_zzzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyyy_0[i] * pb_y + g_0_zzzzzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyyyz_0[i] = 5.0 * g_0_zzzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyyz_0[i] * pb_y + g_0_zzzzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyyzz_0[i] = 4.0 * g_0_zzzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyzz_0[i] * pb_y + g_0_zzzzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_zzzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyzzz_0[i] * pb_y + g_0_zzzzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_zzzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyzzzzz_0[i] = g_0_zzzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xzzzzzz_0[i] = g_0_zzzzzzz_0_xzzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyyyy_0[i] = 7.0 * g_0_zzzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyyyy_0[i] * pb_y + g_0_zzzzzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyyyz_0[i] = 6.0 * g_0_zzzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyyyz_0[i] * pb_y + g_0_zzzzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyyzz_0[i] = 5.0 * g_0_zzzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyyzz_0[i] * pb_y + g_0_zzzzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyzzz_0[i] = 4.0 * g_0_zzzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyzzz_0[i] * pb_y + g_0_zzzzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_zzzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyzzzz_0[i] * pb_y + g_0_zzzzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_zzzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yzzzzzz_0[i] = g_0_zzzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_zzzzzzz_0[i] = g_0_zzzzzzz_0_zzzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1584-1620 components of targeted buffer : SLSK

    auto g_0_zzzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_slsk + 1584);

    auto g_0_zzzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_slsk + 1585);

    auto g_0_zzzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_slsk + 1586);

    auto g_0_zzzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_slsk + 1587);

    auto g_0_zzzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_slsk + 1588);

    auto g_0_zzzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_slsk + 1589);

    auto g_0_zzzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_slsk + 1590);

    auto g_0_zzzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_slsk + 1591);

    auto g_0_zzzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_slsk + 1592);

    auto g_0_zzzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_slsk + 1593);

    auto g_0_zzzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1594);

    auto g_0_zzzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1595);

    auto g_0_zzzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1596);

    auto g_0_zzzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1597);

    auto g_0_zzzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1598);

    auto g_0_zzzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1599);

    auto g_0_zzzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1600);

    auto g_0_zzzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1601);

    auto g_0_zzzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1602);

    auto g_0_zzzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1603);

    auto g_0_zzzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1604);

    auto g_0_zzzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1605);

    auto g_0_zzzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1606);

    auto g_0_zzzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1607);

    auto g_0_zzzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1608);

    auto g_0_zzzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1609);

    auto g_0_zzzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1610);

    auto g_0_zzzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1611);

    auto g_0_zzzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_slsk + 1612);

    auto g_0_zzzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_slsk + 1613);

    auto g_0_zzzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_slsk + 1614);

    auto g_0_zzzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_slsk + 1615);

    auto g_0_zzzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1616);

    auto g_0_zzzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1617);

    auto g_0_zzzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1618);

    auto g_0_zzzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_slsk + 1619);

    #pragma omp simd aligned(g_0_zzzzzz_0_xxxxxxx_0, g_0_zzzzzz_0_xxxxxxx_1, g_0_zzzzzz_0_xxxxxxy_0, g_0_zzzzzz_0_xxxxxxy_1, g_0_zzzzzz_0_xxxxxxz_0, g_0_zzzzzz_0_xxxxxxz_1, g_0_zzzzzz_0_xxxxxyy_0, g_0_zzzzzz_0_xxxxxyy_1, g_0_zzzzzz_0_xxxxxyz_0, g_0_zzzzzz_0_xxxxxyz_1, g_0_zzzzzz_0_xxxxxzz_0, g_0_zzzzzz_0_xxxxxzz_1, g_0_zzzzzz_0_xxxxyyy_0, g_0_zzzzzz_0_xxxxyyy_1, g_0_zzzzzz_0_xxxxyyz_0, g_0_zzzzzz_0_xxxxyyz_1, g_0_zzzzzz_0_xxxxyzz_0, g_0_zzzzzz_0_xxxxyzz_1, g_0_zzzzzz_0_xxxxzzz_0, g_0_zzzzzz_0_xxxxzzz_1, g_0_zzzzzz_0_xxxyyyy_0, g_0_zzzzzz_0_xxxyyyy_1, g_0_zzzzzz_0_xxxyyyz_0, g_0_zzzzzz_0_xxxyyyz_1, g_0_zzzzzz_0_xxxyyzz_0, g_0_zzzzzz_0_xxxyyzz_1, g_0_zzzzzz_0_xxxyzzz_0, g_0_zzzzzz_0_xxxyzzz_1, g_0_zzzzzz_0_xxxzzzz_0, g_0_zzzzzz_0_xxxzzzz_1, g_0_zzzzzz_0_xxyyyyy_0, g_0_zzzzzz_0_xxyyyyy_1, g_0_zzzzzz_0_xxyyyyz_0, g_0_zzzzzz_0_xxyyyyz_1, g_0_zzzzzz_0_xxyyyzz_0, g_0_zzzzzz_0_xxyyyzz_1, g_0_zzzzzz_0_xxyyzzz_0, g_0_zzzzzz_0_xxyyzzz_1, g_0_zzzzzz_0_xxyzzzz_0, g_0_zzzzzz_0_xxyzzzz_1, g_0_zzzzzz_0_xxzzzzz_0, g_0_zzzzzz_0_xxzzzzz_1, g_0_zzzzzz_0_xyyyyyy_0, g_0_zzzzzz_0_xyyyyyy_1, g_0_zzzzzz_0_xyyyyyz_0, g_0_zzzzzz_0_xyyyyyz_1, g_0_zzzzzz_0_xyyyyzz_0, g_0_zzzzzz_0_xyyyyzz_1, g_0_zzzzzz_0_xyyyzzz_0, g_0_zzzzzz_0_xyyyzzz_1, g_0_zzzzzz_0_xyyzzzz_0, g_0_zzzzzz_0_xyyzzzz_1, g_0_zzzzzz_0_xyzzzzz_0, g_0_zzzzzz_0_xyzzzzz_1, g_0_zzzzzz_0_xzzzzzz_0, g_0_zzzzzz_0_xzzzzzz_1, g_0_zzzzzz_0_yyyyyyy_0, g_0_zzzzzz_0_yyyyyyy_1, g_0_zzzzzz_0_yyyyyyz_0, g_0_zzzzzz_0_yyyyyyz_1, g_0_zzzzzz_0_yyyyyzz_0, g_0_zzzzzz_0_yyyyyzz_1, g_0_zzzzzz_0_yyyyzzz_0, g_0_zzzzzz_0_yyyyzzz_1, g_0_zzzzzz_0_yyyzzzz_0, g_0_zzzzzz_0_yyyzzzz_1, g_0_zzzzzz_0_yyzzzzz_0, g_0_zzzzzz_0_yyzzzzz_1, g_0_zzzzzz_0_yzzzzzz_0, g_0_zzzzzz_0_yzzzzzz_1, g_0_zzzzzz_0_zzzzzzz_0, g_0_zzzzzz_0_zzzzzzz_1, g_0_zzzzzzz_0_xxxxxx_1, g_0_zzzzzzz_0_xxxxxxx_0, g_0_zzzzzzz_0_xxxxxxx_1, g_0_zzzzzzz_0_xxxxxxy_0, g_0_zzzzzzz_0_xxxxxxy_1, g_0_zzzzzzz_0_xxxxxxz_0, g_0_zzzzzzz_0_xxxxxxz_1, g_0_zzzzzzz_0_xxxxxy_1, g_0_zzzzzzz_0_xxxxxyy_0, g_0_zzzzzzz_0_xxxxxyy_1, g_0_zzzzzzz_0_xxxxxyz_0, g_0_zzzzzzz_0_xxxxxyz_1, g_0_zzzzzzz_0_xxxxxz_1, g_0_zzzzzzz_0_xxxxxzz_0, g_0_zzzzzzz_0_xxxxxzz_1, g_0_zzzzzzz_0_xxxxyy_1, g_0_zzzzzzz_0_xxxxyyy_0, g_0_zzzzzzz_0_xxxxyyy_1, g_0_zzzzzzz_0_xxxxyyz_0, g_0_zzzzzzz_0_xxxxyyz_1, g_0_zzzzzzz_0_xxxxyz_1, g_0_zzzzzzz_0_xxxxyzz_0, g_0_zzzzzzz_0_xxxxyzz_1, g_0_zzzzzzz_0_xxxxzz_1, g_0_zzzzzzz_0_xxxxzzz_0, g_0_zzzzzzz_0_xxxxzzz_1, g_0_zzzzzzz_0_xxxyyy_1, g_0_zzzzzzz_0_xxxyyyy_0, g_0_zzzzzzz_0_xxxyyyy_1, g_0_zzzzzzz_0_xxxyyyz_0, g_0_zzzzzzz_0_xxxyyyz_1, g_0_zzzzzzz_0_xxxyyz_1, g_0_zzzzzzz_0_xxxyyzz_0, g_0_zzzzzzz_0_xxxyyzz_1, g_0_zzzzzzz_0_xxxyzz_1, g_0_zzzzzzz_0_xxxyzzz_0, g_0_zzzzzzz_0_xxxyzzz_1, g_0_zzzzzzz_0_xxxzzz_1, g_0_zzzzzzz_0_xxxzzzz_0, g_0_zzzzzzz_0_xxxzzzz_1, g_0_zzzzzzz_0_xxyyyy_1, g_0_zzzzzzz_0_xxyyyyy_0, g_0_zzzzzzz_0_xxyyyyy_1, g_0_zzzzzzz_0_xxyyyyz_0, g_0_zzzzzzz_0_xxyyyyz_1, g_0_zzzzzzz_0_xxyyyz_1, g_0_zzzzzzz_0_xxyyyzz_0, g_0_zzzzzzz_0_xxyyyzz_1, g_0_zzzzzzz_0_xxyyzz_1, g_0_zzzzzzz_0_xxyyzzz_0, g_0_zzzzzzz_0_xxyyzzz_1, g_0_zzzzzzz_0_xxyzzz_1, g_0_zzzzzzz_0_xxyzzzz_0, g_0_zzzzzzz_0_xxyzzzz_1, g_0_zzzzzzz_0_xxzzzz_1, g_0_zzzzzzz_0_xxzzzzz_0, g_0_zzzzzzz_0_xxzzzzz_1, g_0_zzzzzzz_0_xyyyyy_1, g_0_zzzzzzz_0_xyyyyyy_0, g_0_zzzzzzz_0_xyyyyyy_1, g_0_zzzzzzz_0_xyyyyyz_0, g_0_zzzzzzz_0_xyyyyyz_1, g_0_zzzzzzz_0_xyyyyz_1, g_0_zzzzzzz_0_xyyyyzz_0, g_0_zzzzzzz_0_xyyyyzz_1, g_0_zzzzzzz_0_xyyyzz_1, g_0_zzzzzzz_0_xyyyzzz_0, g_0_zzzzzzz_0_xyyyzzz_1, g_0_zzzzzzz_0_xyyzzz_1, g_0_zzzzzzz_0_xyyzzzz_0, g_0_zzzzzzz_0_xyyzzzz_1, g_0_zzzzzzz_0_xyzzzz_1, g_0_zzzzzzz_0_xyzzzzz_0, g_0_zzzzzzz_0_xyzzzzz_1, g_0_zzzzzzz_0_xzzzzz_1, g_0_zzzzzzz_0_xzzzzzz_0, g_0_zzzzzzz_0_xzzzzzz_1, g_0_zzzzzzz_0_yyyyyy_1, g_0_zzzzzzz_0_yyyyyyy_0, g_0_zzzzzzz_0_yyyyyyy_1, g_0_zzzzzzz_0_yyyyyyz_0, g_0_zzzzzzz_0_yyyyyyz_1, g_0_zzzzzzz_0_yyyyyz_1, g_0_zzzzzzz_0_yyyyyzz_0, g_0_zzzzzzz_0_yyyyyzz_1, g_0_zzzzzzz_0_yyyyzz_1, g_0_zzzzzzz_0_yyyyzzz_0, g_0_zzzzzzz_0_yyyyzzz_1, g_0_zzzzzzz_0_yyyzzz_1, g_0_zzzzzzz_0_yyyzzzz_0, g_0_zzzzzzz_0_yyyzzzz_1, g_0_zzzzzzz_0_yyzzzz_1, g_0_zzzzzzz_0_yyzzzzz_0, g_0_zzzzzzz_0_yyzzzzz_1, g_0_zzzzzzz_0_yzzzzz_1, g_0_zzzzzzz_0_yzzzzzz_0, g_0_zzzzzzz_0_yzzzzzz_1, g_0_zzzzzzz_0_zzzzzz_1, g_0_zzzzzzz_0_zzzzzzz_0, g_0_zzzzzzz_0_zzzzzzz_1, g_0_zzzzzzzz_0_xxxxxxx_0, g_0_zzzzzzzz_0_xxxxxxy_0, g_0_zzzzzzzz_0_xxxxxxz_0, g_0_zzzzzzzz_0_xxxxxyy_0, g_0_zzzzzzzz_0_xxxxxyz_0, g_0_zzzzzzzz_0_xxxxxzz_0, g_0_zzzzzzzz_0_xxxxyyy_0, g_0_zzzzzzzz_0_xxxxyyz_0, g_0_zzzzzzzz_0_xxxxyzz_0, g_0_zzzzzzzz_0_xxxxzzz_0, g_0_zzzzzzzz_0_xxxyyyy_0, g_0_zzzzzzzz_0_xxxyyyz_0, g_0_zzzzzzzz_0_xxxyyzz_0, g_0_zzzzzzzz_0_xxxyzzz_0, g_0_zzzzzzzz_0_xxxzzzz_0, g_0_zzzzzzzz_0_xxyyyyy_0, g_0_zzzzzzzz_0_xxyyyyz_0, g_0_zzzzzzzz_0_xxyyyzz_0, g_0_zzzzzzzz_0_xxyyzzz_0, g_0_zzzzzzzz_0_xxyzzzz_0, g_0_zzzzzzzz_0_xxzzzzz_0, g_0_zzzzzzzz_0_xyyyyyy_0, g_0_zzzzzzzz_0_xyyyyyz_0, g_0_zzzzzzzz_0_xyyyyzz_0, g_0_zzzzzzzz_0_xyyyzzz_0, g_0_zzzzzzzz_0_xyyzzzz_0, g_0_zzzzzzzz_0_xyzzzzz_0, g_0_zzzzzzzz_0_xzzzzzz_0, g_0_zzzzzzzz_0_yyyyyyy_0, g_0_zzzzzzzz_0_yyyyyyz_0, g_0_zzzzzzzz_0_yyyyyzz_0, g_0_zzzzzzzz_0_yyyyzzz_0, g_0_zzzzzzzz_0_yyyzzzz_0, g_0_zzzzzzzz_0_yyzzzzz_0, g_0_zzzzzzzz_0_yzzzzzz_0, g_0_zzzzzzzz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzzz_0_xxxxxxx_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxxx_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxxxx_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxxxy_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxxy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxxxy_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxxxz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxxz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxxz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxxz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxxyy_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxxyy_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxxyz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxyz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxxzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxzz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxyyy_0[i] = 7.0 * g_0_zzzzzz_0_xxxxyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxyyy_0[i] * pb_z + g_0_zzzzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxyyz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyyz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxyzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyzz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxzzz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyyyy_0[i] = 7.0 * g_0_zzzzzz_0_xxxyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxyyyy_0[i] * pb_z + g_0_zzzzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyyyz_0[i] = 7.0 * g_0_zzzzzz_0_xxxyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyyz_0[i] * pb_z + g_0_zzzzzzz_0_xxxyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyyzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyzz_0[i] * pb_z + g_0_zzzzzzz_0_xxxyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyzzz_0[i] * pb_z + g_0_zzzzzzz_0_xxxyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxzzzz_0[i] * pb_z + g_0_zzzzzzz_0_xxxzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyyyy_0[i] = 7.0 * g_0_zzzzzz_0_xxyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxyyyyy_0[i] * pb_z + g_0_zzzzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyyyz_0[i] = 7.0 * g_0_zzzzzz_0_xxyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyyz_0[i] * pb_z + g_0_zzzzzzz_0_xxyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyyzz_0[i] = 7.0 * g_0_zzzzzz_0_xxyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyzz_0[i] * pb_z + g_0_zzzzzzz_0_xxyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyzzz_0[i] * pb_z + g_0_zzzzzzz_0_xxyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzzzz_0[i] * pb_z + g_0_zzzzzzz_0_xxyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxzzzzz_0[i] * pb_z + g_0_zzzzzzz_0_xxzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyyyy_0[i] = 7.0 * g_0_zzzzzz_0_xyyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xyyyyyy_0[i] * pb_z + g_0_zzzzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyyyz_0[i] = 7.0 * g_0_zzzzzz_0_xyyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyyz_0[i] * pb_z + g_0_zzzzzzz_0_xyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyyzz_0[i] = 7.0 * g_0_zzzzzz_0_xyyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyzz_0[i] * pb_z + g_0_zzzzzzz_0_xyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyzzz_0[i] = 7.0 * g_0_zzzzzz_0_xyyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyzzz_0[i] * pb_z + g_0_zzzzzzz_0_xyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xyyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzzzz_0[i] * pb_z + g_0_zzzzzzz_0_xyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xyzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzzzz_0[i] * pb_z + g_0_zzzzzzz_0_xyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xzzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzzzzzz_0[i] * pb_z + g_0_zzzzzzz_0_xzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyyyy_0[i] = 7.0 * g_0_zzzzzz_0_yyyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yyyyyyy_0[i] * pb_z + g_0_zzzzzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyyyz_0[i] = 7.0 * g_0_zzzzzz_0_yyyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyyyz_0[i] * pb_z + g_0_zzzzzzz_0_yyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyyzz_0[i] = 7.0 * g_0_zzzzzz_0_yyyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyyzz_0[i] * pb_z + g_0_zzzzzzz_0_yyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyzzz_0[i] = 7.0 * g_0_zzzzzz_0_yyyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyzzz_0[i] * pb_z + g_0_zzzzzzz_0_yyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyzzzz_0[i] = 7.0 * g_0_zzzzzz_0_yyyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyzzzz_0[i] * pb_z + g_0_zzzzzzz_0_yyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_yyzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyzzzzz_0[i] * pb_z + g_0_zzzzzzz_0_yyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yzzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_yzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzzzzzz_0[i] * pb_z + g_0_zzzzzzz_0_yzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_zzzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_zzzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_zzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_zzzzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_zzzzzzz_0[i] * pb_z + g_0_zzzzzzz_0_zzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace
