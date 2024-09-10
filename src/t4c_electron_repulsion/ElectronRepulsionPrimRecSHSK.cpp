#include "ElectronRepulsionPrimRecSHSK.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_shsk(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_shsk,
                                  size_t idx_eri_0_sfsk,
                                  size_t idx_eri_1_sfsk,
                                  size_t idx_eri_1_sgsi,
                                  size_t idx_eri_0_sgsk,
                                  size_t idx_eri_1_sgsk,
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

    /// Set up components of auxilary buffer : SFSK

    auto g_0_xxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk);

    auto g_0_xxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 1);

    auto g_0_xxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 2);

    auto g_0_xxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 3);

    auto g_0_xxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 4);

    auto g_0_xxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 5);

    auto g_0_xxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 6);

    auto g_0_xxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 7);

    auto g_0_xxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 8);

    auto g_0_xxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 9);

    auto g_0_xxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 10);

    auto g_0_xxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 11);

    auto g_0_xxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 12);

    auto g_0_xxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 13);

    auto g_0_xxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 14);

    auto g_0_xxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 15);

    auto g_0_xxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 16);

    auto g_0_xxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 17);

    auto g_0_xxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 18);

    auto g_0_xxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 19);

    auto g_0_xxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 20);

    auto g_0_xxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 21);

    auto g_0_xxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 22);

    auto g_0_xxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 23);

    auto g_0_xxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 24);

    auto g_0_xxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 25);

    auto g_0_xxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 26);

    auto g_0_xxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 27);

    auto g_0_xxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 28);

    auto g_0_xxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 29);

    auto g_0_xxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 30);

    auto g_0_xxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 31);

    auto g_0_xxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 32);

    auto g_0_xxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 33);

    auto g_0_xxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 34);

    auto g_0_xxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 35);

    auto g_0_xxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 36);

    auto g_0_xxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 38);

    auto g_0_xxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 41);

    auto g_0_xxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 45);

    auto g_0_xxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 50);

    auto g_0_xxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 56);

    auto g_0_xxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 63);

    auto g_0_xxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 72);

    auto g_0_xxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 73);

    auto g_0_xxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 75);

    auto g_0_xxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 78);

    auto g_0_xxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 82);

    auto g_0_xxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 87);

    auto g_0_xxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 93);

    auto g_0_xyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 109);

    auto g_0_xyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 111);

    auto g_0_xyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 112);

    auto g_0_xyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 114);

    auto g_0_xyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 115);

    auto g_0_xyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 116);

    auto g_0_xyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 118);

    auto g_0_xyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 119);

    auto g_0_xyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 120);

    auto g_0_xyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 121);

    auto g_0_xyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 123);

    auto g_0_xyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 124);

    auto g_0_xyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 125);

    auto g_0_xyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 126);

    auto g_0_xyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 127);

    auto g_0_xyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 129);

    auto g_0_xyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 130);

    auto g_0_xyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 131);

    auto g_0_xyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 132);

    auto g_0_xyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 133);

    auto g_0_xyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 134);

    auto g_0_xyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 136);

    auto g_0_xyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 137);

    auto g_0_xyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 138);

    auto g_0_xyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 139);

    auto g_0_xyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 140);

    auto g_0_xyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 141);

    auto g_0_xyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 142);

    auto g_0_xyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 143);

    auto g_0_xzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 182);

    auto g_0_xzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 184);

    auto g_0_xzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 185);

    auto g_0_xzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 187);

    auto g_0_xzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 188);

    auto g_0_xzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 189);

    auto g_0_xzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 191);

    auto g_0_xzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 192);

    auto g_0_xzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 193);

    auto g_0_xzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 194);

    auto g_0_xzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 196);

    auto g_0_xzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 197);

    auto g_0_xzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 198);

    auto g_0_xzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 199);

    auto g_0_xzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 200);

    auto g_0_xzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 202);

    auto g_0_xzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 203);

    auto g_0_xzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 204);

    auto g_0_xzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 205);

    auto g_0_xzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 206);

    auto g_0_xzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 207);

    auto g_0_xzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 208);

    auto g_0_xzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 209);

    auto g_0_xzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 210);

    auto g_0_xzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 211);

    auto g_0_xzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 212);

    auto g_0_xzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 213);

    auto g_0_xzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 214);

    auto g_0_xzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 215);

    auto g_0_yyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 216);

    auto g_0_yyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 217);

    auto g_0_yyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 218);

    auto g_0_yyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 219);

    auto g_0_yyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 220);

    auto g_0_yyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 221);

    auto g_0_yyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 222);

    auto g_0_yyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 223);

    auto g_0_yyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 224);

    auto g_0_yyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 225);

    auto g_0_yyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 226);

    auto g_0_yyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 227);

    auto g_0_yyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 228);

    auto g_0_yyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 229);

    auto g_0_yyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 230);

    auto g_0_yyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 231);

    auto g_0_yyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 232);

    auto g_0_yyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 233);

    auto g_0_yyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 234);

    auto g_0_yyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 235);

    auto g_0_yyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 236);

    auto g_0_yyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 237);

    auto g_0_yyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 238);

    auto g_0_yyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 239);

    auto g_0_yyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 240);

    auto g_0_yyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 241);

    auto g_0_yyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 242);

    auto g_0_yyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 243);

    auto g_0_yyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 244);

    auto g_0_yyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 245);

    auto g_0_yyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 246);

    auto g_0_yyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 247);

    auto g_0_yyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 248);

    auto g_0_yyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 249);

    auto g_0_yyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 250);

    auto g_0_yyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 251);

    auto g_0_yyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 253);

    auto g_0_yyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 255);

    auto g_0_yyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 258);

    auto g_0_yyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 262);

    auto g_0_yyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 267);

    auto g_0_yyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 273);

    auto g_0_yyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 280);

    auto g_0_yzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 288);

    auto g_0_yzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 290);

    auto g_0_yzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 292);

    auto g_0_yzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 293);

    auto g_0_yzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 295);

    auto g_0_yzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 296);

    auto g_0_yzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 297);

    auto g_0_yzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 299);

    auto g_0_yzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 300);

    auto g_0_yzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 301);

    auto g_0_yzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 302);

    auto g_0_yzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 304);

    auto g_0_yzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 305);

    auto g_0_yzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 306);

    auto g_0_yzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 307);

    auto g_0_yzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 308);

    auto g_0_yzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 310);

    auto g_0_yzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 311);

    auto g_0_yzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 312);

    auto g_0_yzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 313);

    auto g_0_yzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 314);

    auto g_0_yzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 315);

    auto g_0_yzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 317);

    auto g_0_yzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 318);

    auto g_0_yzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 319);

    auto g_0_yzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 320);

    auto g_0_yzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 321);

    auto g_0_yzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 322);

    auto g_0_yzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 323);

    auto g_0_zzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sfsk + 324);

    auto g_0_zzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sfsk + 325);

    auto g_0_zzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sfsk + 326);

    auto g_0_zzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sfsk + 327);

    auto g_0_zzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sfsk + 328);

    auto g_0_zzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sfsk + 329);

    auto g_0_zzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sfsk + 330);

    auto g_0_zzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sfsk + 331);

    auto g_0_zzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sfsk + 332);

    auto g_0_zzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sfsk + 333);

    auto g_0_zzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 334);

    auto g_0_zzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 335);

    auto g_0_zzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 336);

    auto g_0_zzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 337);

    auto g_0_zzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 338);

    auto g_0_zzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 339);

    auto g_0_zzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 340);

    auto g_0_zzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 341);

    auto g_0_zzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 342);

    auto g_0_zzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 343);

    auto g_0_zzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 344);

    auto g_0_zzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 345);

    auto g_0_zzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 346);

    auto g_0_zzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 347);

    auto g_0_zzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 348);

    auto g_0_zzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 349);

    auto g_0_zzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 350);

    auto g_0_zzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 351);

    auto g_0_zzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sfsk + 352);

    auto g_0_zzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sfsk + 353);

    auto g_0_zzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sfsk + 354);

    auto g_0_zzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sfsk + 355);

    auto g_0_zzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 356);

    auto g_0_zzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 357);

    auto g_0_zzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 358);

    auto g_0_zzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sfsk + 359);

    /// Set up components of auxilary buffer : SFSK

    auto g_0_xxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk);

    auto g_0_xxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 1);

    auto g_0_xxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 2);

    auto g_0_xxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 3);

    auto g_0_xxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 4);

    auto g_0_xxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 5);

    auto g_0_xxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 6);

    auto g_0_xxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 7);

    auto g_0_xxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 8);

    auto g_0_xxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 9);

    auto g_0_xxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 10);

    auto g_0_xxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 11);

    auto g_0_xxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 12);

    auto g_0_xxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 13);

    auto g_0_xxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 14);

    auto g_0_xxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 15);

    auto g_0_xxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 16);

    auto g_0_xxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 17);

    auto g_0_xxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 18);

    auto g_0_xxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 19);

    auto g_0_xxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 20);

    auto g_0_xxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 21);

    auto g_0_xxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 22);

    auto g_0_xxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 23);

    auto g_0_xxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 24);

    auto g_0_xxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 25);

    auto g_0_xxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 26);

    auto g_0_xxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 27);

    auto g_0_xxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 28);

    auto g_0_xxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 29);

    auto g_0_xxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 30);

    auto g_0_xxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 31);

    auto g_0_xxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 32);

    auto g_0_xxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 33);

    auto g_0_xxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 34);

    auto g_0_xxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 35);

    auto g_0_xxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk + 36);

    auto g_0_xxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 38);

    auto g_0_xxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 41);

    auto g_0_xxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 45);

    auto g_0_xxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 50);

    auto g_0_xxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 56);

    auto g_0_xxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 63);

    auto g_0_xxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk + 72);

    auto g_0_xxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 73);

    auto g_0_xxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 75);

    auto g_0_xxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 78);

    auto g_0_xxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 82);

    auto g_0_xxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 87);

    auto g_0_xxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 93);

    auto g_0_xyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 109);

    auto g_0_xyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 111);

    auto g_0_xyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 112);

    auto g_0_xyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 114);

    auto g_0_xyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 115);

    auto g_0_xyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 116);

    auto g_0_xyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 118);

    auto g_0_xyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 119);

    auto g_0_xyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 120);

    auto g_0_xyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 121);

    auto g_0_xyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 123);

    auto g_0_xyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 124);

    auto g_0_xyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 125);

    auto g_0_xyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 126);

    auto g_0_xyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 127);

    auto g_0_xyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 129);

    auto g_0_xyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 130);

    auto g_0_xyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 131);

    auto g_0_xyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 132);

    auto g_0_xyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 133);

    auto g_0_xyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 134);

    auto g_0_xyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 136);

    auto g_0_xyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 137);

    auto g_0_xyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 138);

    auto g_0_xyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 139);

    auto g_0_xyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 140);

    auto g_0_xyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 141);

    auto g_0_xyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 142);

    auto g_0_xyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 143);

    auto g_0_xzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 182);

    auto g_0_xzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 184);

    auto g_0_xzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 185);

    auto g_0_xzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 187);

    auto g_0_xzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 188);

    auto g_0_xzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 189);

    auto g_0_xzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 191);

    auto g_0_xzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 192);

    auto g_0_xzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 193);

    auto g_0_xzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 194);

    auto g_0_xzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 196);

    auto g_0_xzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 197);

    auto g_0_xzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 198);

    auto g_0_xzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 199);

    auto g_0_xzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 200);

    auto g_0_xzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 202);

    auto g_0_xzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 203);

    auto g_0_xzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 204);

    auto g_0_xzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 205);

    auto g_0_xzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 206);

    auto g_0_xzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 207);

    auto g_0_xzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 208);

    auto g_0_xzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 209);

    auto g_0_xzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 210);

    auto g_0_xzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 211);

    auto g_0_xzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 212);

    auto g_0_xzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 213);

    auto g_0_xzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 214);

    auto g_0_xzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 215);

    auto g_0_yyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk + 216);

    auto g_0_yyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 217);

    auto g_0_yyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 218);

    auto g_0_yyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 219);

    auto g_0_yyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 220);

    auto g_0_yyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 221);

    auto g_0_yyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 222);

    auto g_0_yyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 223);

    auto g_0_yyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 224);

    auto g_0_yyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 225);

    auto g_0_yyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 226);

    auto g_0_yyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 227);

    auto g_0_yyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 228);

    auto g_0_yyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 229);

    auto g_0_yyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 230);

    auto g_0_yyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 231);

    auto g_0_yyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 232);

    auto g_0_yyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 233);

    auto g_0_yyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 234);

    auto g_0_yyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 235);

    auto g_0_yyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 236);

    auto g_0_yyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 237);

    auto g_0_yyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 238);

    auto g_0_yyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 239);

    auto g_0_yyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 240);

    auto g_0_yyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 241);

    auto g_0_yyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 242);

    auto g_0_yyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 243);

    auto g_0_yyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 244);

    auto g_0_yyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 245);

    auto g_0_yyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 246);

    auto g_0_yyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 247);

    auto g_0_yyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 248);

    auto g_0_yyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 249);

    auto g_0_yyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 250);

    auto g_0_yyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 251);

    auto g_0_yyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 253);

    auto g_0_yyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 255);

    auto g_0_yyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 258);

    auto g_0_yyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 262);

    auto g_0_yyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 267);

    auto g_0_yyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 273);

    auto g_0_yyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 280);

    auto g_0_yzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk + 288);

    auto g_0_yzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 290);

    auto g_0_yzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 292);

    auto g_0_yzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 293);

    auto g_0_yzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 295);

    auto g_0_yzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 296);

    auto g_0_yzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 297);

    auto g_0_yzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 299);

    auto g_0_yzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 300);

    auto g_0_yzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 301);

    auto g_0_yzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 302);

    auto g_0_yzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 304);

    auto g_0_yzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 305);

    auto g_0_yzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 306);

    auto g_0_yzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 307);

    auto g_0_yzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 308);

    auto g_0_yzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 310);

    auto g_0_yzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 311);

    auto g_0_yzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 312);

    auto g_0_yzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 313);

    auto g_0_yzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 314);

    auto g_0_yzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 315);

    auto g_0_yzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 317);

    auto g_0_yzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 318);

    auto g_0_yzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 319);

    auto g_0_yzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 320);

    auto g_0_yzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 321);

    auto g_0_yzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 322);

    auto g_0_yzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 323);

    auto g_0_zzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sfsk + 324);

    auto g_0_zzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sfsk + 325);

    auto g_0_zzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sfsk + 326);

    auto g_0_zzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sfsk + 327);

    auto g_0_zzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sfsk + 328);

    auto g_0_zzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sfsk + 329);

    auto g_0_zzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sfsk + 330);

    auto g_0_zzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sfsk + 331);

    auto g_0_zzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sfsk + 332);

    auto g_0_zzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sfsk + 333);

    auto g_0_zzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 334);

    auto g_0_zzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 335);

    auto g_0_zzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 336);

    auto g_0_zzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 337);

    auto g_0_zzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 338);

    auto g_0_zzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 339);

    auto g_0_zzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 340);

    auto g_0_zzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 341);

    auto g_0_zzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 342);

    auto g_0_zzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 343);

    auto g_0_zzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 344);

    auto g_0_zzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 345);

    auto g_0_zzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 346);

    auto g_0_zzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 347);

    auto g_0_zzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 348);

    auto g_0_zzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 349);

    auto g_0_zzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 350);

    auto g_0_zzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 351);

    auto g_0_zzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sfsk + 352);

    auto g_0_zzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sfsk + 353);

    auto g_0_zzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sfsk + 354);

    auto g_0_zzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sfsk + 355);

    auto g_0_zzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 356);

    auto g_0_zzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 357);

    auto g_0_zzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 358);

    auto g_0_zzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sfsk + 359);

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

    auto g_0_xxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 58);

    auto g_0_xxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 60);

    auto g_0_xxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 61);

    auto g_0_xxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 63);

    auto g_0_xxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 64);

    auto g_0_xxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 65);

    auto g_0_xxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 67);

    auto g_0_xxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 68);

    auto g_0_xxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 69);

    auto g_0_xxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 70);

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

    auto g_0_yyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 310);

    auto g_0_yyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 312);

    auto g_0_yyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 313);

    auto g_0_yyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 315);

    auto g_0_yyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 316);

    auto g_0_yyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 317);

    auto g_0_yyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 319);

    auto g_0_yyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 320);

    auto g_0_yyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 321);

    auto g_0_yyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 322);

    auto g_0_yyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 324);

    auto g_0_yyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 325);

    auto g_0_yyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 326);

    auto g_0_yyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 327);

    auto g_0_yyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 328);

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

    /// Set up components of auxilary buffer : SGSK

    auto g_0_xxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk);

    auto g_0_xxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 1);

    auto g_0_xxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 2);

    auto g_0_xxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 3);

    auto g_0_xxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 4);

    auto g_0_xxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 5);

    auto g_0_xxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 6);

    auto g_0_xxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 7);

    auto g_0_xxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 8);

    auto g_0_xxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 9);

    auto g_0_xxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 10);

    auto g_0_xxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 11);

    auto g_0_xxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 12);

    auto g_0_xxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 13);

    auto g_0_xxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 14);

    auto g_0_xxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 15);

    auto g_0_xxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 16);

    auto g_0_xxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 17);

    auto g_0_xxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 18);

    auto g_0_xxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 19);

    auto g_0_xxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 20);

    auto g_0_xxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 21);

    auto g_0_xxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 22);

    auto g_0_xxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 23);

    auto g_0_xxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 24);

    auto g_0_xxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 25);

    auto g_0_xxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 26);

    auto g_0_xxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 27);

    auto g_0_xxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 28);

    auto g_0_xxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 29);

    auto g_0_xxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 30);

    auto g_0_xxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 31);

    auto g_0_xxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 32);

    auto g_0_xxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 33);

    auto g_0_xxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 34);

    auto g_0_xxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 35);

    auto g_0_xxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 36);

    auto g_0_xxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 37);

    auto g_0_xxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 38);

    auto g_0_xxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 39);

    auto g_0_xxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 41);

    auto g_0_xxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 42);

    auto g_0_xxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 45);

    auto g_0_xxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 46);

    auto g_0_xxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 50);

    auto g_0_xxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 51);

    auto g_0_xxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 56);

    auto g_0_xxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 57);

    auto g_0_xxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 63);

    auto g_0_xxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 64);

    auto g_0_xxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 72);

    auto g_0_xxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 73);

    auto g_0_xxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 74);

    auto g_0_xxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 75);

    auto g_0_xxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 76);

    auto g_0_xxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 77);

    auto g_0_xxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 78);

    auto g_0_xxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 79);

    auto g_0_xxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 80);

    auto g_0_xxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 81);

    auto g_0_xxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 82);

    auto g_0_xxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 83);

    auto g_0_xxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 84);

    auto g_0_xxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 85);

    auto g_0_xxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 86);

    auto g_0_xxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 87);

    auto g_0_xxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 88);

    auto g_0_xxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 89);

    auto g_0_xxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 90);

    auto g_0_xxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 91);

    auto g_0_xxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 92);

    auto g_0_xxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 93);

    auto g_0_xxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 94);

    auto g_0_xxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 95);

    auto g_0_xxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 96);

    auto g_0_xxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 97);

    auto g_0_xxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 98);

    auto g_0_xxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 99);

    auto g_0_xxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 101);

    auto g_0_xxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 102);

    auto g_0_xxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 103);

    auto g_0_xxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 104);

    auto g_0_xxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 105);

    auto g_0_xxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 106);

    auto g_0_xxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 107);

    auto g_0_xxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 108);

    auto g_0_xxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 109);

    auto g_0_xxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 110);

    auto g_0_xxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 111);

    auto g_0_xxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 112);

    auto g_0_xxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 113);

    auto g_0_xxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 114);

    auto g_0_xxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 115);

    auto g_0_xxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 116);

    auto g_0_xxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 117);

    auto g_0_xxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 118);

    auto g_0_xxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 119);

    auto g_0_xxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 120);

    auto g_0_xxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 121);

    auto g_0_xxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 122);

    auto g_0_xxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 123);

    auto g_0_xxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 124);

    auto g_0_xxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 125);

    auto g_0_xxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 126);

    auto g_0_xxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 127);

    auto g_0_xxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 128);

    auto g_0_xxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 129);

    auto g_0_xxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 130);

    auto g_0_xxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 131);

    auto g_0_xxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 132);

    auto g_0_xxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 133);

    auto g_0_xxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 134);

    auto g_0_xxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 135);

    auto g_0_xxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 136);

    auto g_0_xxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 137);

    auto g_0_xxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 138);

    auto g_0_xxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 139);

    auto g_0_xxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 140);

    auto g_0_xxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 141);

    auto g_0_xxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 142);

    auto g_0_xxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 143);

    auto g_0_xxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 180);

    auto g_0_xxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 181);

    auto g_0_xxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 182);

    auto g_0_xxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 183);

    auto g_0_xxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 184);

    auto g_0_xxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 185);

    auto g_0_xxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 186);

    auto g_0_xxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 187);

    auto g_0_xxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 188);

    auto g_0_xxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 189);

    auto g_0_xxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 190);

    auto g_0_xxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 191);

    auto g_0_xxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 192);

    auto g_0_xxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 193);

    auto g_0_xxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 194);

    auto g_0_xxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 195);

    auto g_0_xxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 196);

    auto g_0_xxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 197);

    auto g_0_xxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 198);

    auto g_0_xxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 199);

    auto g_0_xxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 200);

    auto g_0_xxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 201);

    auto g_0_xxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 202);

    auto g_0_xxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 203);

    auto g_0_xxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 204);

    auto g_0_xxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 205);

    auto g_0_xxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 206);

    auto g_0_xxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 207);

    auto g_0_xxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 208);

    auto g_0_xxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 209);

    auto g_0_xxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 210);

    auto g_0_xxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 211);

    auto g_0_xxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 212);

    auto g_0_xxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 213);

    auto g_0_xxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 214);

    auto g_0_xxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 215);

    auto g_0_xyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 216);

    auto g_0_xyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 217);

    auto g_0_xyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 219);

    auto g_0_xyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 220);

    auto g_0_xyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 222);

    auto g_0_xyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 223);

    auto g_0_xyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 224);

    auto g_0_xyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 226);

    auto g_0_xyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 227);

    auto g_0_xyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 228);

    auto g_0_xyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 229);

    auto g_0_xyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 231);

    auto g_0_xyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 232);

    auto g_0_xyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 233);

    auto g_0_xyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 234);

    auto g_0_xyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 235);

    auto g_0_xyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 237);

    auto g_0_xyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 238);

    auto g_0_xyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 239);

    auto g_0_xyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 240);

    auto g_0_xyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 241);

    auto g_0_xyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 242);

    auto g_0_xyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 244);

    auto g_0_xyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 245);

    auto g_0_xyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 246);

    auto g_0_xyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 247);

    auto g_0_xyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 248);

    auto g_0_xyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 249);

    auto g_0_xyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 250);

    auto g_0_xyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 251);

    auto g_0_xzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 324);

    auto g_0_xzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 326);

    auto g_0_xzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 328);

    auto g_0_xzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 329);

    auto g_0_xzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 331);

    auto g_0_xzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 332);

    auto g_0_xzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 333);

    auto g_0_xzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 335);

    auto g_0_xzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 336);

    auto g_0_xzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 337);

    auto g_0_xzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 338);

    auto g_0_xzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 340);

    auto g_0_xzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 341);

    auto g_0_xzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 342);

    auto g_0_xzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 343);

    auto g_0_xzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 344);

    auto g_0_xzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 346);

    auto g_0_xzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 347);

    auto g_0_xzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 348);

    auto g_0_xzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 349);

    auto g_0_xzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 350);

    auto g_0_xzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 351);

    auto g_0_xzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 352);

    auto g_0_xzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 353);

    auto g_0_xzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 354);

    auto g_0_xzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 355);

    auto g_0_xzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 356);

    auto g_0_xzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 357);

    auto g_0_xzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 358);

    auto g_0_xzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 359);

    auto g_0_yyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 360);

    auto g_0_yyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 361);

    auto g_0_yyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 362);

    auto g_0_yyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 363);

    auto g_0_yyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 364);

    auto g_0_yyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 365);

    auto g_0_yyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 366);

    auto g_0_yyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 367);

    auto g_0_yyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 368);

    auto g_0_yyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 369);

    auto g_0_yyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 370);

    auto g_0_yyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 371);

    auto g_0_yyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 372);

    auto g_0_yyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 373);

    auto g_0_yyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 374);

    auto g_0_yyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 375);

    auto g_0_yyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 376);

    auto g_0_yyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 377);

    auto g_0_yyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 378);

    auto g_0_yyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 379);

    auto g_0_yyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 380);

    auto g_0_yyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 381);

    auto g_0_yyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 382);

    auto g_0_yyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 383);

    auto g_0_yyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 384);

    auto g_0_yyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 385);

    auto g_0_yyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 386);

    auto g_0_yyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 387);

    auto g_0_yyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 388);

    auto g_0_yyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 389);

    auto g_0_yyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 390);

    auto g_0_yyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 391);

    auto g_0_yyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 392);

    auto g_0_yyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 393);

    auto g_0_yyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 394);

    auto g_0_yyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 395);

    auto g_0_yyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 397);

    auto g_0_yyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 398);

    auto g_0_yyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 399);

    auto g_0_yyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 400);

    auto g_0_yyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 401);

    auto g_0_yyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 402);

    auto g_0_yyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 403);

    auto g_0_yyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 404);

    auto g_0_yyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 405);

    auto g_0_yyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 406);

    auto g_0_yyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 407);

    auto g_0_yyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 408);

    auto g_0_yyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 409);

    auto g_0_yyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 410);

    auto g_0_yyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 411);

    auto g_0_yyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 412);

    auto g_0_yyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 413);

    auto g_0_yyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 414);

    auto g_0_yyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 415);

    auto g_0_yyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 416);

    auto g_0_yyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 417);

    auto g_0_yyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 418);

    auto g_0_yyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 419);

    auto g_0_yyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 420);

    auto g_0_yyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 421);

    auto g_0_yyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 422);

    auto g_0_yyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 423);

    auto g_0_yyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 424);

    auto g_0_yyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 425);

    auto g_0_yyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 426);

    auto g_0_yyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 427);

    auto g_0_yyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 428);

    auto g_0_yyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 429);

    auto g_0_yyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 430);

    auto g_0_yyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 431);

    auto g_0_yyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 432);

    auto g_0_yyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 433);

    auto g_0_yyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 434);

    auto g_0_yyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 435);

    auto g_0_yyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 436);

    auto g_0_yyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 437);

    auto g_0_yyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 438);

    auto g_0_yyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 439);

    auto g_0_yyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 440);

    auto g_0_yyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 441);

    auto g_0_yyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 442);

    auto g_0_yyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 443);

    auto g_0_yyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 444);

    auto g_0_yyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 445);

    auto g_0_yyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 446);

    auto g_0_yyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 447);

    auto g_0_yyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 448);

    auto g_0_yyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 449);

    auto g_0_yyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 450);

    auto g_0_yyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 451);

    auto g_0_yyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 452);

    auto g_0_yyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 453);

    auto g_0_yyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 454);

    auto g_0_yyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 455);

    auto g_0_yyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 456);

    auto g_0_yyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 457);

    auto g_0_yyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 458);

    auto g_0_yyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 459);

    auto g_0_yyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 460);

    auto g_0_yyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 461);

    auto g_0_yyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 462);

    auto g_0_yyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 463);

    auto g_0_yyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 464);

    auto g_0_yyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 465);

    auto g_0_yyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 466);

    auto g_0_yyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 467);

    auto g_0_yzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 468);

    auto g_0_yzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 469);

    auto g_0_yzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 470);

    auto g_0_yzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 471);

    auto g_0_yzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 472);

    auto g_0_yzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 473);

    auto g_0_yzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 474);

    auto g_0_yzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 475);

    auto g_0_yzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 476);

    auto g_0_yzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 477);

    auto g_0_yzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 478);

    auto g_0_yzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 479);

    auto g_0_yzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 480);

    auto g_0_yzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 481);

    auto g_0_yzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 482);

    auto g_0_yzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 483);

    auto g_0_yzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 484);

    auto g_0_yzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 485);

    auto g_0_yzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 486);

    auto g_0_yzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 487);

    auto g_0_yzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 488);

    auto g_0_yzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 489);

    auto g_0_yzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 490);

    auto g_0_yzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 491);

    auto g_0_yzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 492);

    auto g_0_yzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 493);

    auto g_0_yzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 494);

    auto g_0_yzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 495);

    auto g_0_yzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 496);

    auto g_0_yzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 497);

    auto g_0_yzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 498);

    auto g_0_yzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 499);

    auto g_0_yzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 500);

    auto g_0_yzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 501);

    auto g_0_yzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 502);

    auto g_0_yzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 503);

    auto g_0_zzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sgsk + 504);

    auto g_0_zzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sgsk + 505);

    auto g_0_zzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sgsk + 506);

    auto g_0_zzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sgsk + 507);

    auto g_0_zzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sgsk + 508);

    auto g_0_zzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sgsk + 509);

    auto g_0_zzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sgsk + 510);

    auto g_0_zzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sgsk + 511);

    auto g_0_zzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sgsk + 512);

    auto g_0_zzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sgsk + 513);

    auto g_0_zzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 514);

    auto g_0_zzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 515);

    auto g_0_zzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 516);

    auto g_0_zzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 517);

    auto g_0_zzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 518);

    auto g_0_zzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 519);

    auto g_0_zzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 520);

    auto g_0_zzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 521);

    auto g_0_zzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 522);

    auto g_0_zzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 523);

    auto g_0_zzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 524);

    auto g_0_zzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 525);

    auto g_0_zzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 526);

    auto g_0_zzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 527);

    auto g_0_zzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 528);

    auto g_0_zzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 529);

    auto g_0_zzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 530);

    auto g_0_zzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 531);

    auto g_0_zzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sgsk + 532);

    auto g_0_zzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sgsk + 533);

    auto g_0_zzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sgsk + 534);

    auto g_0_zzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sgsk + 535);

    auto g_0_zzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 536);

    auto g_0_zzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 537);

    auto g_0_zzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 538);

    auto g_0_zzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sgsk + 539);

    /// Set up components of auxilary buffer : SGSK

    auto g_0_xxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk);

    auto g_0_xxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 1);

    auto g_0_xxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 2);

    auto g_0_xxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 3);

    auto g_0_xxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 4);

    auto g_0_xxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 5);

    auto g_0_xxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 6);

    auto g_0_xxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 7);

    auto g_0_xxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 8);

    auto g_0_xxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 9);

    auto g_0_xxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 10);

    auto g_0_xxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 11);

    auto g_0_xxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 12);

    auto g_0_xxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 13);

    auto g_0_xxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 14);

    auto g_0_xxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 15);

    auto g_0_xxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 16);

    auto g_0_xxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 17);

    auto g_0_xxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 18);

    auto g_0_xxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 19);

    auto g_0_xxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 20);

    auto g_0_xxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 21);

    auto g_0_xxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 22);

    auto g_0_xxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 23);

    auto g_0_xxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 24);

    auto g_0_xxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 25);

    auto g_0_xxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 26);

    auto g_0_xxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 27);

    auto g_0_xxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 28);

    auto g_0_xxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 29);

    auto g_0_xxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 30);

    auto g_0_xxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 31);

    auto g_0_xxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 32);

    auto g_0_xxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 33);

    auto g_0_xxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 34);

    auto g_0_xxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 35);

    auto g_0_xxxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 36);

    auto g_0_xxxy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 37);

    auto g_0_xxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 38);

    auto g_0_xxxy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 39);

    auto g_0_xxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 41);

    auto g_0_xxxy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 42);

    auto g_0_xxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 45);

    auto g_0_xxxy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 46);

    auto g_0_xxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 50);

    auto g_0_xxxy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 51);

    auto g_0_xxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 56);

    auto g_0_xxxy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 57);

    auto g_0_xxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 63);

    auto g_0_xxxy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 64);

    auto g_0_xxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 72);

    auto g_0_xxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 73);

    auto g_0_xxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 74);

    auto g_0_xxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 75);

    auto g_0_xxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 76);

    auto g_0_xxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 77);

    auto g_0_xxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 78);

    auto g_0_xxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 79);

    auto g_0_xxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 80);

    auto g_0_xxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 81);

    auto g_0_xxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 82);

    auto g_0_xxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 83);

    auto g_0_xxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 84);

    auto g_0_xxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 85);

    auto g_0_xxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 86);

    auto g_0_xxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 87);

    auto g_0_xxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 88);

    auto g_0_xxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 89);

    auto g_0_xxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 90);

    auto g_0_xxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 91);

    auto g_0_xxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 92);

    auto g_0_xxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 93);

    auto g_0_xxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 94);

    auto g_0_xxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 95);

    auto g_0_xxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 96);

    auto g_0_xxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 97);

    auto g_0_xxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 98);

    auto g_0_xxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 99);

    auto g_0_xxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 101);

    auto g_0_xxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 102);

    auto g_0_xxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 103);

    auto g_0_xxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 104);

    auto g_0_xxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 105);

    auto g_0_xxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 106);

    auto g_0_xxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 107);

    auto g_0_xxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 108);

    auto g_0_xxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 109);

    auto g_0_xxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 110);

    auto g_0_xxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 111);

    auto g_0_xxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 112);

    auto g_0_xxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 113);

    auto g_0_xxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 114);

    auto g_0_xxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 115);

    auto g_0_xxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 116);

    auto g_0_xxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 117);

    auto g_0_xxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 118);

    auto g_0_xxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 119);

    auto g_0_xxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 120);

    auto g_0_xxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 121);

    auto g_0_xxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 122);

    auto g_0_xxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 123);

    auto g_0_xxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 124);

    auto g_0_xxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 125);

    auto g_0_xxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 126);

    auto g_0_xxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 127);

    auto g_0_xxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 128);

    auto g_0_xxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 129);

    auto g_0_xxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 130);

    auto g_0_xxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 131);

    auto g_0_xxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 132);

    auto g_0_xxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 133);

    auto g_0_xxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 134);

    auto g_0_xxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 135);

    auto g_0_xxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 136);

    auto g_0_xxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 137);

    auto g_0_xxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 138);

    auto g_0_xxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 139);

    auto g_0_xxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 140);

    auto g_0_xxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 141);

    auto g_0_xxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 142);

    auto g_0_xxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 143);

    auto g_0_xxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 180);

    auto g_0_xxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 181);

    auto g_0_xxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 182);

    auto g_0_xxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 183);

    auto g_0_xxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 184);

    auto g_0_xxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 185);

    auto g_0_xxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 186);

    auto g_0_xxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 187);

    auto g_0_xxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 188);

    auto g_0_xxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 189);

    auto g_0_xxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 190);

    auto g_0_xxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 191);

    auto g_0_xxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 192);

    auto g_0_xxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 193);

    auto g_0_xxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 194);

    auto g_0_xxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 195);

    auto g_0_xxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 196);

    auto g_0_xxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 197);

    auto g_0_xxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 198);

    auto g_0_xxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 199);

    auto g_0_xxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 200);

    auto g_0_xxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 201);

    auto g_0_xxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 202);

    auto g_0_xxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 203);

    auto g_0_xxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 204);

    auto g_0_xxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 205);

    auto g_0_xxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 206);

    auto g_0_xxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 207);

    auto g_0_xxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 208);

    auto g_0_xxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 209);

    auto g_0_xxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 210);

    auto g_0_xxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 211);

    auto g_0_xxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 212);

    auto g_0_xxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 213);

    auto g_0_xxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 214);

    auto g_0_xxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 215);

    auto g_0_xyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 216);

    auto g_0_xyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 217);

    auto g_0_xyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 219);

    auto g_0_xyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 220);

    auto g_0_xyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 222);

    auto g_0_xyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 223);

    auto g_0_xyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 224);

    auto g_0_xyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 226);

    auto g_0_xyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 227);

    auto g_0_xyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 228);

    auto g_0_xyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 229);

    auto g_0_xyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 231);

    auto g_0_xyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 232);

    auto g_0_xyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 233);

    auto g_0_xyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 234);

    auto g_0_xyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 235);

    auto g_0_xyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 237);

    auto g_0_xyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 238);

    auto g_0_xyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 239);

    auto g_0_xyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 240);

    auto g_0_xyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 241);

    auto g_0_xyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 242);

    auto g_0_xyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 244);

    auto g_0_xyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 245);

    auto g_0_xyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 246);

    auto g_0_xyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 247);

    auto g_0_xyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 248);

    auto g_0_xyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 249);

    auto g_0_xyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 250);

    auto g_0_xyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 251);

    auto g_0_xzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 324);

    auto g_0_xzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 326);

    auto g_0_xzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 328);

    auto g_0_xzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 329);

    auto g_0_xzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 331);

    auto g_0_xzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 332);

    auto g_0_xzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 333);

    auto g_0_xzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 335);

    auto g_0_xzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 336);

    auto g_0_xzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 337);

    auto g_0_xzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 338);

    auto g_0_xzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 340);

    auto g_0_xzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 341);

    auto g_0_xzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 342);

    auto g_0_xzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 343);

    auto g_0_xzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 344);

    auto g_0_xzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 346);

    auto g_0_xzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 347);

    auto g_0_xzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 348);

    auto g_0_xzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 349);

    auto g_0_xzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 350);

    auto g_0_xzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 351);

    auto g_0_xzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 352);

    auto g_0_xzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 353);

    auto g_0_xzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 354);

    auto g_0_xzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 355);

    auto g_0_xzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 356);

    auto g_0_xzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 357);

    auto g_0_xzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 358);

    auto g_0_xzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 359);

    auto g_0_yyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 360);

    auto g_0_yyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 361);

    auto g_0_yyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 362);

    auto g_0_yyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 363);

    auto g_0_yyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 364);

    auto g_0_yyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 365);

    auto g_0_yyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 366);

    auto g_0_yyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 367);

    auto g_0_yyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 368);

    auto g_0_yyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 369);

    auto g_0_yyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 370);

    auto g_0_yyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 371);

    auto g_0_yyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 372);

    auto g_0_yyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 373);

    auto g_0_yyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 374);

    auto g_0_yyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 375);

    auto g_0_yyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 376);

    auto g_0_yyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 377);

    auto g_0_yyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 378);

    auto g_0_yyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 379);

    auto g_0_yyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 380);

    auto g_0_yyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 381);

    auto g_0_yyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 382);

    auto g_0_yyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 383);

    auto g_0_yyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 384);

    auto g_0_yyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 385);

    auto g_0_yyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 386);

    auto g_0_yyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 387);

    auto g_0_yyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 388);

    auto g_0_yyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 389);

    auto g_0_yyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 390);

    auto g_0_yyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 391);

    auto g_0_yyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 392);

    auto g_0_yyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 393);

    auto g_0_yyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 394);

    auto g_0_yyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 395);

    auto g_0_yyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 397);

    auto g_0_yyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 398);

    auto g_0_yyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 399);

    auto g_0_yyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 400);

    auto g_0_yyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 401);

    auto g_0_yyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 402);

    auto g_0_yyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 403);

    auto g_0_yyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 404);

    auto g_0_yyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 405);

    auto g_0_yyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 406);

    auto g_0_yyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 407);

    auto g_0_yyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 408);

    auto g_0_yyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 409);

    auto g_0_yyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 410);

    auto g_0_yyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 411);

    auto g_0_yyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 412);

    auto g_0_yyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 413);

    auto g_0_yyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 414);

    auto g_0_yyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 415);

    auto g_0_yyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 416);

    auto g_0_yyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 417);

    auto g_0_yyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 418);

    auto g_0_yyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 419);

    auto g_0_yyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 420);

    auto g_0_yyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 421);

    auto g_0_yyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 422);

    auto g_0_yyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 423);

    auto g_0_yyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 424);

    auto g_0_yyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 425);

    auto g_0_yyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 426);

    auto g_0_yyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 427);

    auto g_0_yyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 428);

    auto g_0_yyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 429);

    auto g_0_yyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 430);

    auto g_0_yyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 431);

    auto g_0_yyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 432);

    auto g_0_yyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 433);

    auto g_0_yyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 434);

    auto g_0_yyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 435);

    auto g_0_yyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 436);

    auto g_0_yyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 437);

    auto g_0_yyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 438);

    auto g_0_yyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 439);

    auto g_0_yyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 440);

    auto g_0_yyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 441);

    auto g_0_yyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 442);

    auto g_0_yyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 443);

    auto g_0_yyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 444);

    auto g_0_yyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 445);

    auto g_0_yyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 446);

    auto g_0_yyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 447);

    auto g_0_yyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 448);

    auto g_0_yyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 449);

    auto g_0_yyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 450);

    auto g_0_yyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 451);

    auto g_0_yyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 452);

    auto g_0_yyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 453);

    auto g_0_yyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 454);

    auto g_0_yyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 455);

    auto g_0_yyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 456);

    auto g_0_yyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 457);

    auto g_0_yyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 458);

    auto g_0_yyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 459);

    auto g_0_yyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 460);

    auto g_0_yyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 461);

    auto g_0_yyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 462);

    auto g_0_yyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 463);

    auto g_0_yyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 464);

    auto g_0_yyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 465);

    auto g_0_yyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 466);

    auto g_0_yyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 467);

    auto g_0_yzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 468);

    auto g_0_yzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 469);

    auto g_0_yzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 470);

    auto g_0_yzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 471);

    auto g_0_yzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 472);

    auto g_0_yzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 473);

    auto g_0_yzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 474);

    auto g_0_yzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 475);

    auto g_0_yzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 476);

    auto g_0_yzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 477);

    auto g_0_yzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 478);

    auto g_0_yzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 479);

    auto g_0_yzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 480);

    auto g_0_yzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 481);

    auto g_0_yzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 482);

    auto g_0_yzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 483);

    auto g_0_yzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 484);

    auto g_0_yzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 485);

    auto g_0_yzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 486);

    auto g_0_yzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 487);

    auto g_0_yzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 488);

    auto g_0_yzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 489);

    auto g_0_yzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 490);

    auto g_0_yzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 491);

    auto g_0_yzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 492);

    auto g_0_yzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 493);

    auto g_0_yzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 494);

    auto g_0_yzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 495);

    auto g_0_yzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 496);

    auto g_0_yzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 497);

    auto g_0_yzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 498);

    auto g_0_yzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 499);

    auto g_0_yzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 500);

    auto g_0_yzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 501);

    auto g_0_yzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 502);

    auto g_0_yzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 503);

    auto g_0_zzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sgsk + 504);

    auto g_0_zzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sgsk + 505);

    auto g_0_zzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sgsk + 506);

    auto g_0_zzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sgsk + 507);

    auto g_0_zzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sgsk + 508);

    auto g_0_zzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sgsk + 509);

    auto g_0_zzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sgsk + 510);

    auto g_0_zzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sgsk + 511);

    auto g_0_zzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sgsk + 512);

    auto g_0_zzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sgsk + 513);

    auto g_0_zzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 514);

    auto g_0_zzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 515);

    auto g_0_zzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 516);

    auto g_0_zzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 517);

    auto g_0_zzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 518);

    auto g_0_zzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 519);

    auto g_0_zzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 520);

    auto g_0_zzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 521);

    auto g_0_zzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 522);

    auto g_0_zzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 523);

    auto g_0_zzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 524);

    auto g_0_zzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 525);

    auto g_0_zzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 526);

    auto g_0_zzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 527);

    auto g_0_zzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 528);

    auto g_0_zzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 529);

    auto g_0_zzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 530);

    auto g_0_zzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 531);

    auto g_0_zzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sgsk + 532);

    auto g_0_zzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sgsk + 533);

    auto g_0_zzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sgsk + 534);

    auto g_0_zzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sgsk + 535);

    auto g_0_zzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 536);

    auto g_0_zzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 537);

    auto g_0_zzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 538);

    auto g_0_zzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sgsk + 539);

    /// Set up 0-36 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_xxx_0_xxxxxxx_0, g_0_xxx_0_xxxxxxx_1, g_0_xxx_0_xxxxxxy_0, g_0_xxx_0_xxxxxxy_1, g_0_xxx_0_xxxxxxz_0, g_0_xxx_0_xxxxxxz_1, g_0_xxx_0_xxxxxyy_0, g_0_xxx_0_xxxxxyy_1, g_0_xxx_0_xxxxxyz_0, g_0_xxx_0_xxxxxyz_1, g_0_xxx_0_xxxxxzz_0, g_0_xxx_0_xxxxxzz_1, g_0_xxx_0_xxxxyyy_0, g_0_xxx_0_xxxxyyy_1, g_0_xxx_0_xxxxyyz_0, g_0_xxx_0_xxxxyyz_1, g_0_xxx_0_xxxxyzz_0, g_0_xxx_0_xxxxyzz_1, g_0_xxx_0_xxxxzzz_0, g_0_xxx_0_xxxxzzz_1, g_0_xxx_0_xxxyyyy_0, g_0_xxx_0_xxxyyyy_1, g_0_xxx_0_xxxyyyz_0, g_0_xxx_0_xxxyyyz_1, g_0_xxx_0_xxxyyzz_0, g_0_xxx_0_xxxyyzz_1, g_0_xxx_0_xxxyzzz_0, g_0_xxx_0_xxxyzzz_1, g_0_xxx_0_xxxzzzz_0, g_0_xxx_0_xxxzzzz_1, g_0_xxx_0_xxyyyyy_0, g_0_xxx_0_xxyyyyy_1, g_0_xxx_0_xxyyyyz_0, g_0_xxx_0_xxyyyyz_1, g_0_xxx_0_xxyyyzz_0, g_0_xxx_0_xxyyyzz_1, g_0_xxx_0_xxyyzzz_0, g_0_xxx_0_xxyyzzz_1, g_0_xxx_0_xxyzzzz_0, g_0_xxx_0_xxyzzzz_1, g_0_xxx_0_xxzzzzz_0, g_0_xxx_0_xxzzzzz_1, g_0_xxx_0_xyyyyyy_0, g_0_xxx_0_xyyyyyy_1, g_0_xxx_0_xyyyyyz_0, g_0_xxx_0_xyyyyyz_1, g_0_xxx_0_xyyyyzz_0, g_0_xxx_0_xyyyyzz_1, g_0_xxx_0_xyyyzzz_0, g_0_xxx_0_xyyyzzz_1, g_0_xxx_0_xyyzzzz_0, g_0_xxx_0_xyyzzzz_1, g_0_xxx_0_xyzzzzz_0, g_0_xxx_0_xyzzzzz_1, g_0_xxx_0_xzzzzzz_0, g_0_xxx_0_xzzzzzz_1, g_0_xxx_0_yyyyyyy_0, g_0_xxx_0_yyyyyyy_1, g_0_xxx_0_yyyyyyz_0, g_0_xxx_0_yyyyyyz_1, g_0_xxx_0_yyyyyzz_0, g_0_xxx_0_yyyyyzz_1, g_0_xxx_0_yyyyzzz_0, g_0_xxx_0_yyyyzzz_1, g_0_xxx_0_yyyzzzz_0, g_0_xxx_0_yyyzzzz_1, g_0_xxx_0_yyzzzzz_0, g_0_xxx_0_yyzzzzz_1, g_0_xxx_0_yzzzzzz_0, g_0_xxx_0_yzzzzzz_1, g_0_xxx_0_zzzzzzz_0, g_0_xxx_0_zzzzzzz_1, g_0_xxxx_0_xxxxxx_1, g_0_xxxx_0_xxxxxxx_0, g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxy_0, g_0_xxxx_0_xxxxxxy_1, g_0_xxxx_0_xxxxxxz_0, g_0_xxxx_0_xxxxxxz_1, g_0_xxxx_0_xxxxxy_1, g_0_xxxx_0_xxxxxyy_0, g_0_xxxx_0_xxxxxyy_1, g_0_xxxx_0_xxxxxyz_0, g_0_xxxx_0_xxxxxyz_1, g_0_xxxx_0_xxxxxz_1, g_0_xxxx_0_xxxxxzz_0, g_0_xxxx_0_xxxxxzz_1, g_0_xxxx_0_xxxxyy_1, g_0_xxxx_0_xxxxyyy_0, g_0_xxxx_0_xxxxyyy_1, g_0_xxxx_0_xxxxyyz_0, g_0_xxxx_0_xxxxyyz_1, g_0_xxxx_0_xxxxyz_1, g_0_xxxx_0_xxxxyzz_0, g_0_xxxx_0_xxxxyzz_1, g_0_xxxx_0_xxxxzz_1, g_0_xxxx_0_xxxxzzz_0, g_0_xxxx_0_xxxxzzz_1, g_0_xxxx_0_xxxyyy_1, g_0_xxxx_0_xxxyyyy_0, g_0_xxxx_0_xxxyyyy_1, g_0_xxxx_0_xxxyyyz_0, g_0_xxxx_0_xxxyyyz_1, g_0_xxxx_0_xxxyyz_1, g_0_xxxx_0_xxxyyzz_0, g_0_xxxx_0_xxxyyzz_1, g_0_xxxx_0_xxxyzz_1, g_0_xxxx_0_xxxyzzz_0, g_0_xxxx_0_xxxyzzz_1, g_0_xxxx_0_xxxzzz_1, g_0_xxxx_0_xxxzzzz_0, g_0_xxxx_0_xxxzzzz_1, g_0_xxxx_0_xxyyyy_1, g_0_xxxx_0_xxyyyyy_0, g_0_xxxx_0_xxyyyyy_1, g_0_xxxx_0_xxyyyyz_0, g_0_xxxx_0_xxyyyyz_1, g_0_xxxx_0_xxyyyz_1, g_0_xxxx_0_xxyyyzz_0, g_0_xxxx_0_xxyyyzz_1, g_0_xxxx_0_xxyyzz_1, g_0_xxxx_0_xxyyzzz_0, g_0_xxxx_0_xxyyzzz_1, g_0_xxxx_0_xxyzzz_1, g_0_xxxx_0_xxyzzzz_0, g_0_xxxx_0_xxyzzzz_1, g_0_xxxx_0_xxzzzz_1, g_0_xxxx_0_xxzzzzz_0, g_0_xxxx_0_xxzzzzz_1, g_0_xxxx_0_xyyyyy_1, g_0_xxxx_0_xyyyyyy_0, g_0_xxxx_0_xyyyyyy_1, g_0_xxxx_0_xyyyyyz_0, g_0_xxxx_0_xyyyyyz_1, g_0_xxxx_0_xyyyyz_1, g_0_xxxx_0_xyyyyzz_0, g_0_xxxx_0_xyyyyzz_1, g_0_xxxx_0_xyyyzz_1, g_0_xxxx_0_xyyyzzz_0, g_0_xxxx_0_xyyyzzz_1, g_0_xxxx_0_xyyzzz_1, g_0_xxxx_0_xyyzzzz_0, g_0_xxxx_0_xyyzzzz_1, g_0_xxxx_0_xyzzzz_1, g_0_xxxx_0_xyzzzzz_0, g_0_xxxx_0_xyzzzzz_1, g_0_xxxx_0_xzzzzz_1, g_0_xxxx_0_xzzzzzz_0, g_0_xxxx_0_xzzzzzz_1, g_0_xxxx_0_yyyyyy_1, g_0_xxxx_0_yyyyyyy_0, g_0_xxxx_0_yyyyyyy_1, g_0_xxxx_0_yyyyyyz_0, g_0_xxxx_0_yyyyyyz_1, g_0_xxxx_0_yyyyyz_1, g_0_xxxx_0_yyyyyzz_0, g_0_xxxx_0_yyyyyzz_1, g_0_xxxx_0_yyyyzz_1, g_0_xxxx_0_yyyyzzz_0, g_0_xxxx_0_yyyyzzz_1, g_0_xxxx_0_yyyzzz_1, g_0_xxxx_0_yyyzzzz_0, g_0_xxxx_0_yyyzzzz_1, g_0_xxxx_0_yyzzzz_1, g_0_xxxx_0_yyzzzzz_0, g_0_xxxx_0_yyzzzzz_1, g_0_xxxx_0_yzzzzz_1, g_0_xxxx_0_yzzzzzz_0, g_0_xxxx_0_yzzzzzz_1, g_0_xxxx_0_zzzzzz_1, g_0_xxxx_0_zzzzzzz_0, g_0_xxxx_0_zzzzzzz_1, g_0_xxxxx_0_xxxxxxx_0, g_0_xxxxx_0_xxxxxxy_0, g_0_xxxxx_0_xxxxxxz_0, g_0_xxxxx_0_xxxxxyy_0, g_0_xxxxx_0_xxxxxyz_0, g_0_xxxxx_0_xxxxxzz_0, g_0_xxxxx_0_xxxxyyy_0, g_0_xxxxx_0_xxxxyyz_0, g_0_xxxxx_0_xxxxyzz_0, g_0_xxxxx_0_xxxxzzz_0, g_0_xxxxx_0_xxxyyyy_0, g_0_xxxxx_0_xxxyyyz_0, g_0_xxxxx_0_xxxyyzz_0, g_0_xxxxx_0_xxxyzzz_0, g_0_xxxxx_0_xxxzzzz_0, g_0_xxxxx_0_xxyyyyy_0, g_0_xxxxx_0_xxyyyyz_0, g_0_xxxxx_0_xxyyyzz_0, g_0_xxxxx_0_xxyyzzz_0, g_0_xxxxx_0_xxyzzzz_0, g_0_xxxxx_0_xxzzzzz_0, g_0_xxxxx_0_xyyyyyy_0, g_0_xxxxx_0_xyyyyyz_0, g_0_xxxxx_0_xyyyyzz_0, g_0_xxxxx_0_xyyyzzz_0, g_0_xxxxx_0_xyyzzzz_0, g_0_xxxxx_0_xyzzzzz_0, g_0_xxxxx_0_xzzzzzz_0, g_0_xxxxx_0_yyyyyyy_0, g_0_xxxxx_0_yyyyyyz_0, g_0_xxxxx_0_yyyyyzz_0, g_0_xxxxx_0_yyyyzzz_0, g_0_xxxxx_0_yyyzzzz_0, g_0_xxxxx_0_yyzzzzz_0, g_0_xxxxx_0_yzzzzzz_0, g_0_xxxxx_0_zzzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_xxxxxxx_0[i] = 4.0 * g_0_xxx_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxx_1[i] * fti_ab_0 + 7.0 * g_0_xxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxx_0[i] * pb_x + g_0_xxxx_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxxy_0[i] = 4.0 * g_0_xxx_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxy_0[i] * pb_x + g_0_xxxx_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxxz_0[i] = 4.0 * g_0_xxx_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxz_0[i] * pb_x + g_0_xxxx_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxyy_0[i] = 4.0 * g_0_xxx_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyy_0[i] * pb_x + g_0_xxxx_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxyz_0[i] = 4.0 * g_0_xxx_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyz_0[i] * pb_x + g_0_xxxx_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxzz_0[i] = 4.0 * g_0_xxx_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxzz_0[i] * pb_x + g_0_xxxx_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyyy_0[i] = 4.0 * g_0_xxx_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyy_0[i] * pb_x + g_0_xxxx_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyyz_0[i] = 4.0 * g_0_xxx_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyz_0[i] * pb_x + g_0_xxxx_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyzz_0[i] = 4.0 * g_0_xxx_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyzz_0[i] * pb_x + g_0_xxxx_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxzzz_0[i] = 4.0 * g_0_xxx_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxzzz_0[i] * pb_x + g_0_xxxx_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyyy_0[i] = 4.0 * g_0_xxx_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyy_0[i] * pb_x + g_0_xxxx_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyyz_0[i] = 4.0 * g_0_xxx_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyz_0[i] * pb_x + g_0_xxxx_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyzz_0[i] = 4.0 * g_0_xxx_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyzz_0[i] * pb_x + g_0_xxxx_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyzzz_0[i] = 4.0 * g_0_xxx_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzzz_0[i] * pb_x + g_0_xxxx_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxzzzz_0[i] = 4.0 * g_0_xxx_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxzzzz_0[i] * pb_x + g_0_xxxx_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyyy_0[i] = 4.0 * g_0_xxx_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyy_0[i] * pb_x + g_0_xxxx_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyyz_0[i] = 4.0 * g_0_xxx_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyz_0[i] * pb_x + g_0_xxxx_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyzz_0[i] = 4.0 * g_0_xxx_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyzz_0[i] * pb_x + g_0_xxxx_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyzzz_0[i] = 4.0 * g_0_xxx_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzzz_0[i] * pb_x + g_0_xxxx_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyzzzz_0[i] = 4.0 * g_0_xxx_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzzz_0[i] * pb_x + g_0_xxxx_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxzzzzz_0[i] = 4.0 * g_0_xxx_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzzzzz_0[i] * pb_x + g_0_xxxx_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyyy_0[i] = 4.0 * g_0_xxx_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyy_0[i] * pb_x + g_0_xxxx_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyyz_0[i] = 4.0 * g_0_xxx_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyz_0[i] * pb_x + g_0_xxxx_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyzz_0[i] = 4.0 * g_0_xxx_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyzz_0[i] * pb_x + g_0_xxxx_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyzzz_0[i] = 4.0 * g_0_xxx_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyzzz_0[i] * pb_x + g_0_xxxx_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyzzzz_0[i] = 4.0 * g_0_xxx_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzzzz_0[i] * pb_x + g_0_xxxx_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyzzzzz_0[i] = 4.0 * g_0_xxx_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzzzz_0[i] * pb_x + g_0_xxxx_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xzzzzzz_0[i] = 4.0 * g_0_xxx_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzzzzzz_0[i] * pb_x + g_0_xxxx_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyyy_0[i] = 4.0 * g_0_xxx_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyyy_0[i] * pb_x + g_0_xxxx_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyyz_0[i] = 4.0 * g_0_xxx_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyyz_0[i] * pb_x + g_0_xxxx_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyzz_0[i] = 4.0 * g_0_xxx_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyzz_0[i] * pb_x + g_0_xxxx_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyzzz_0[i] = 4.0 * g_0_xxx_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyzzz_0[i] * pb_x + g_0_xxxx_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyzzzz_0[i] = 4.0 * g_0_xxx_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyzzzz_0[i] * pb_x + g_0_xxxx_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyzzzzz_0[i] = 4.0 * g_0_xxx_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzzzzz_0[i] * pb_x + g_0_xxxx_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yzzzzzz_0[i] = 4.0 * g_0_xxx_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzzzzz_0[i] * pb_x + g_0_xxxx_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_zzzzzzz_0[i] = 4.0 * g_0_xxx_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzzzzz_0[i] * pb_x + g_0_xxxx_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 36-72 components of targeted buffer : SHSK

    auto g_0_xxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 36);

    auto g_0_xxxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 37);

    auto g_0_xxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 38);

    auto g_0_xxxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 39);

    auto g_0_xxxxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 40);

    auto g_0_xxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 41);

    auto g_0_xxxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 42);

    auto g_0_xxxxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 43);

    auto g_0_xxxxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 44);

    auto g_0_xxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 45);

    auto g_0_xxxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 46);

    auto g_0_xxxxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 47);

    auto g_0_xxxxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 48);

    auto g_0_xxxxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 49);

    auto g_0_xxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 50);

    auto g_0_xxxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 51);

    auto g_0_xxxxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 52);

    auto g_0_xxxxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 53);

    auto g_0_xxxxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 54);

    auto g_0_xxxxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 55);

    auto g_0_xxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 56);

    auto g_0_xxxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 57);

    auto g_0_xxxxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 58);

    auto g_0_xxxxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 59);

    auto g_0_xxxxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 60);

    auto g_0_xxxxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 61);

    auto g_0_xxxxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 62);

    auto g_0_xxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 63);

    auto g_0_xxxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 64);

    auto g_0_xxxxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 65);

    auto g_0_xxxxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 66);

    auto g_0_xxxxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 67);

    auto g_0_xxxxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 68);

    auto g_0_xxxxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 69);

    auto g_0_xxxxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 70);

    auto g_0_xxxxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 71);

    #pragma omp simd aligned(g_0_xxxx_0_xxxxxx_1, g_0_xxxx_0_xxxxxxx_0, g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxy_0, g_0_xxxx_0_xxxxxxy_1, g_0_xxxx_0_xxxxxxz_0, g_0_xxxx_0_xxxxxxz_1, g_0_xxxx_0_xxxxxy_1, g_0_xxxx_0_xxxxxyy_0, g_0_xxxx_0_xxxxxyy_1, g_0_xxxx_0_xxxxxyz_0, g_0_xxxx_0_xxxxxyz_1, g_0_xxxx_0_xxxxxz_1, g_0_xxxx_0_xxxxxzz_0, g_0_xxxx_0_xxxxxzz_1, g_0_xxxx_0_xxxxyy_1, g_0_xxxx_0_xxxxyyy_0, g_0_xxxx_0_xxxxyyy_1, g_0_xxxx_0_xxxxyyz_0, g_0_xxxx_0_xxxxyyz_1, g_0_xxxx_0_xxxxyz_1, g_0_xxxx_0_xxxxyzz_0, g_0_xxxx_0_xxxxyzz_1, g_0_xxxx_0_xxxxzz_1, g_0_xxxx_0_xxxxzzz_0, g_0_xxxx_0_xxxxzzz_1, g_0_xxxx_0_xxxyyy_1, g_0_xxxx_0_xxxyyyy_0, g_0_xxxx_0_xxxyyyy_1, g_0_xxxx_0_xxxyyyz_0, g_0_xxxx_0_xxxyyyz_1, g_0_xxxx_0_xxxyyz_1, g_0_xxxx_0_xxxyyzz_0, g_0_xxxx_0_xxxyyzz_1, g_0_xxxx_0_xxxyzz_1, g_0_xxxx_0_xxxyzzz_0, g_0_xxxx_0_xxxyzzz_1, g_0_xxxx_0_xxxzzz_1, g_0_xxxx_0_xxxzzzz_0, g_0_xxxx_0_xxxzzzz_1, g_0_xxxx_0_xxyyyy_1, g_0_xxxx_0_xxyyyyy_0, g_0_xxxx_0_xxyyyyy_1, g_0_xxxx_0_xxyyyyz_0, g_0_xxxx_0_xxyyyyz_1, g_0_xxxx_0_xxyyyz_1, g_0_xxxx_0_xxyyyzz_0, g_0_xxxx_0_xxyyyzz_1, g_0_xxxx_0_xxyyzz_1, g_0_xxxx_0_xxyyzzz_0, g_0_xxxx_0_xxyyzzz_1, g_0_xxxx_0_xxyzzz_1, g_0_xxxx_0_xxyzzzz_0, g_0_xxxx_0_xxyzzzz_1, g_0_xxxx_0_xxzzzz_1, g_0_xxxx_0_xxzzzzz_0, g_0_xxxx_0_xxzzzzz_1, g_0_xxxx_0_xyyyyy_1, g_0_xxxx_0_xyyyyyy_0, g_0_xxxx_0_xyyyyyy_1, g_0_xxxx_0_xyyyyyz_0, g_0_xxxx_0_xyyyyyz_1, g_0_xxxx_0_xyyyyz_1, g_0_xxxx_0_xyyyyzz_0, g_0_xxxx_0_xyyyyzz_1, g_0_xxxx_0_xyyyzz_1, g_0_xxxx_0_xyyyzzz_0, g_0_xxxx_0_xyyyzzz_1, g_0_xxxx_0_xyyzzz_1, g_0_xxxx_0_xyyzzzz_0, g_0_xxxx_0_xyyzzzz_1, g_0_xxxx_0_xyzzzz_1, g_0_xxxx_0_xyzzzzz_0, g_0_xxxx_0_xyzzzzz_1, g_0_xxxx_0_xzzzzz_1, g_0_xxxx_0_xzzzzzz_0, g_0_xxxx_0_xzzzzzz_1, g_0_xxxx_0_yyyyyy_1, g_0_xxxx_0_yyyyyyy_0, g_0_xxxx_0_yyyyyyy_1, g_0_xxxx_0_yyyyyyz_0, g_0_xxxx_0_yyyyyyz_1, g_0_xxxx_0_yyyyyz_1, g_0_xxxx_0_yyyyyzz_0, g_0_xxxx_0_yyyyyzz_1, g_0_xxxx_0_yyyyzz_1, g_0_xxxx_0_yyyyzzz_0, g_0_xxxx_0_yyyyzzz_1, g_0_xxxx_0_yyyzzz_1, g_0_xxxx_0_yyyzzzz_0, g_0_xxxx_0_yyyzzzz_1, g_0_xxxx_0_yyzzzz_1, g_0_xxxx_0_yyzzzzz_0, g_0_xxxx_0_yyzzzzz_1, g_0_xxxx_0_yzzzzz_1, g_0_xxxx_0_yzzzzzz_0, g_0_xxxx_0_yzzzzzz_1, g_0_xxxx_0_zzzzzz_1, g_0_xxxx_0_zzzzzzz_0, g_0_xxxx_0_zzzzzzz_1, g_0_xxxxy_0_xxxxxxx_0, g_0_xxxxy_0_xxxxxxy_0, g_0_xxxxy_0_xxxxxxz_0, g_0_xxxxy_0_xxxxxyy_0, g_0_xxxxy_0_xxxxxyz_0, g_0_xxxxy_0_xxxxxzz_0, g_0_xxxxy_0_xxxxyyy_0, g_0_xxxxy_0_xxxxyyz_0, g_0_xxxxy_0_xxxxyzz_0, g_0_xxxxy_0_xxxxzzz_0, g_0_xxxxy_0_xxxyyyy_0, g_0_xxxxy_0_xxxyyyz_0, g_0_xxxxy_0_xxxyyzz_0, g_0_xxxxy_0_xxxyzzz_0, g_0_xxxxy_0_xxxzzzz_0, g_0_xxxxy_0_xxyyyyy_0, g_0_xxxxy_0_xxyyyyz_0, g_0_xxxxy_0_xxyyyzz_0, g_0_xxxxy_0_xxyyzzz_0, g_0_xxxxy_0_xxyzzzz_0, g_0_xxxxy_0_xxzzzzz_0, g_0_xxxxy_0_xyyyyyy_0, g_0_xxxxy_0_xyyyyyz_0, g_0_xxxxy_0_xyyyyzz_0, g_0_xxxxy_0_xyyyzzz_0, g_0_xxxxy_0_xyyzzzz_0, g_0_xxxxy_0_xyzzzzz_0, g_0_xxxxy_0_xzzzzzz_0, g_0_xxxxy_0_yyyyyyy_0, g_0_xxxxy_0_yyyyyyz_0, g_0_xxxxy_0_yyyyyzz_0, g_0_xxxxy_0_yyyyzzz_0, g_0_xxxxy_0_yyyzzzz_0, g_0_xxxxy_0_yyzzzzz_0, g_0_xxxxy_0_yzzzzzz_0, g_0_xxxxy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxy_0_xxxxxxx_0[i] = g_0_xxxx_0_xxxxxxx_0[i] * pb_y + g_0_xxxx_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxxy_0[i] = g_0_xxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxy_0[i] * pb_y + g_0_xxxx_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxxz_0[i] = g_0_xxxx_0_xxxxxxz_0[i] * pb_y + g_0_xxxx_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxyy_0[i] = 2.0 * g_0_xxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyy_0[i] * pb_y + g_0_xxxx_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxyz_0[i] = g_0_xxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyz_0[i] * pb_y + g_0_xxxx_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxzz_0[i] = g_0_xxxx_0_xxxxxzz_0[i] * pb_y + g_0_xxxx_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyyy_0[i] = 3.0 * g_0_xxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyy_0[i] * pb_y + g_0_xxxx_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyyz_0[i] = 2.0 * g_0_xxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyz_0[i] * pb_y + g_0_xxxx_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyzz_0[i] = g_0_xxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyzz_0[i] * pb_y + g_0_xxxx_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxzzz_0[i] = g_0_xxxx_0_xxxxzzz_0[i] * pb_y + g_0_xxxx_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyyy_0[i] = 4.0 * g_0_xxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyy_0[i] * pb_y + g_0_xxxx_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyyz_0[i] = 3.0 * g_0_xxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyz_0[i] * pb_y + g_0_xxxx_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyzz_0[i] = 2.0 * g_0_xxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyzz_0[i] * pb_y + g_0_xxxx_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyzzz_0[i] = g_0_xxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzzz_0[i] * pb_y + g_0_xxxx_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxzzzz_0[i] = g_0_xxxx_0_xxxzzzz_0[i] * pb_y + g_0_xxxx_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyyy_0[i] = 5.0 * g_0_xxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyy_0[i] * pb_y + g_0_xxxx_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyyz_0[i] = 4.0 * g_0_xxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyz_0[i] * pb_y + g_0_xxxx_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyzz_0[i] = 3.0 * g_0_xxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyzz_0[i] * pb_y + g_0_xxxx_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyzzz_0[i] = 2.0 * g_0_xxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzzz_0[i] * pb_y + g_0_xxxx_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyzzzz_0[i] = g_0_xxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzzz_0[i] * pb_y + g_0_xxxx_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxzzzzz_0[i] = g_0_xxxx_0_xxzzzzz_0[i] * pb_y + g_0_xxxx_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyyy_0[i] = 6.0 * g_0_xxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyy_0[i] * pb_y + g_0_xxxx_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyyz_0[i] = 5.0 * g_0_xxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyz_0[i] * pb_y + g_0_xxxx_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyzz_0[i] = 4.0 * g_0_xxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyzz_0[i] * pb_y + g_0_xxxx_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyzzz_0[i] = 3.0 * g_0_xxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyzzz_0[i] * pb_y + g_0_xxxx_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyzzzz_0[i] = 2.0 * g_0_xxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzzzz_0[i] * pb_y + g_0_xxxx_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyzzzzz_0[i] = g_0_xxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzzzz_0[i] * pb_y + g_0_xxxx_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xzzzzzz_0[i] = g_0_xxxx_0_xzzzzzz_0[i] * pb_y + g_0_xxxx_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyyy_0[i] = 7.0 * g_0_xxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyyy_0[i] * pb_y + g_0_xxxx_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyyz_0[i] = 6.0 * g_0_xxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyyz_0[i] * pb_y + g_0_xxxx_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyzz_0[i] = 5.0 * g_0_xxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyzz_0[i] * pb_y + g_0_xxxx_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyzzz_0[i] = 4.0 * g_0_xxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyzzz_0[i] * pb_y + g_0_xxxx_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyzzzz_0[i] = 3.0 * g_0_xxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyzzzz_0[i] * pb_y + g_0_xxxx_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyzzzzz_0[i] = 2.0 * g_0_xxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzzzzz_0[i] * pb_y + g_0_xxxx_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yzzzzzz_0[i] = g_0_xxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzzzzz_0[i] * pb_y + g_0_xxxx_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_zzzzzzz_0[i] = g_0_xxxx_0_zzzzzzz_0[i] * pb_y + g_0_xxxx_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 72-108 components of targeted buffer : SHSK

    auto g_0_xxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 72);

    auto g_0_xxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 73);

    auto g_0_xxxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 74);

    auto g_0_xxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 75);

    auto g_0_xxxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 76);

    auto g_0_xxxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 77);

    auto g_0_xxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 78);

    auto g_0_xxxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 79);

    auto g_0_xxxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 80);

    auto g_0_xxxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 81);

    auto g_0_xxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 82);

    auto g_0_xxxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 83);

    auto g_0_xxxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 84);

    auto g_0_xxxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 85);

    auto g_0_xxxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 86);

    auto g_0_xxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 87);

    auto g_0_xxxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 88);

    auto g_0_xxxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 89);

    auto g_0_xxxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 90);

    auto g_0_xxxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 91);

    auto g_0_xxxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 92);

    auto g_0_xxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 93);

    auto g_0_xxxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 94);

    auto g_0_xxxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 95);

    auto g_0_xxxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 96);

    auto g_0_xxxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 97);

    auto g_0_xxxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 98);

    auto g_0_xxxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 99);

    auto g_0_xxxxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 100);

    auto g_0_xxxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 101);

    auto g_0_xxxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 102);

    auto g_0_xxxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 103);

    auto g_0_xxxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 104);

    auto g_0_xxxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 105);

    auto g_0_xxxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 106);

    auto g_0_xxxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 107);

    #pragma omp simd aligned(g_0_xxxx_0_xxxxxx_1, g_0_xxxx_0_xxxxxxx_0, g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxy_0, g_0_xxxx_0_xxxxxxy_1, g_0_xxxx_0_xxxxxxz_0, g_0_xxxx_0_xxxxxxz_1, g_0_xxxx_0_xxxxxy_1, g_0_xxxx_0_xxxxxyy_0, g_0_xxxx_0_xxxxxyy_1, g_0_xxxx_0_xxxxxyz_0, g_0_xxxx_0_xxxxxyz_1, g_0_xxxx_0_xxxxxz_1, g_0_xxxx_0_xxxxxzz_0, g_0_xxxx_0_xxxxxzz_1, g_0_xxxx_0_xxxxyy_1, g_0_xxxx_0_xxxxyyy_0, g_0_xxxx_0_xxxxyyy_1, g_0_xxxx_0_xxxxyyz_0, g_0_xxxx_0_xxxxyyz_1, g_0_xxxx_0_xxxxyz_1, g_0_xxxx_0_xxxxyzz_0, g_0_xxxx_0_xxxxyzz_1, g_0_xxxx_0_xxxxzz_1, g_0_xxxx_0_xxxxzzz_0, g_0_xxxx_0_xxxxzzz_1, g_0_xxxx_0_xxxyyy_1, g_0_xxxx_0_xxxyyyy_0, g_0_xxxx_0_xxxyyyy_1, g_0_xxxx_0_xxxyyyz_0, g_0_xxxx_0_xxxyyyz_1, g_0_xxxx_0_xxxyyz_1, g_0_xxxx_0_xxxyyzz_0, g_0_xxxx_0_xxxyyzz_1, g_0_xxxx_0_xxxyzz_1, g_0_xxxx_0_xxxyzzz_0, g_0_xxxx_0_xxxyzzz_1, g_0_xxxx_0_xxxzzz_1, g_0_xxxx_0_xxxzzzz_0, g_0_xxxx_0_xxxzzzz_1, g_0_xxxx_0_xxyyyy_1, g_0_xxxx_0_xxyyyyy_0, g_0_xxxx_0_xxyyyyy_1, g_0_xxxx_0_xxyyyyz_0, g_0_xxxx_0_xxyyyyz_1, g_0_xxxx_0_xxyyyz_1, g_0_xxxx_0_xxyyyzz_0, g_0_xxxx_0_xxyyyzz_1, g_0_xxxx_0_xxyyzz_1, g_0_xxxx_0_xxyyzzz_0, g_0_xxxx_0_xxyyzzz_1, g_0_xxxx_0_xxyzzz_1, g_0_xxxx_0_xxyzzzz_0, g_0_xxxx_0_xxyzzzz_1, g_0_xxxx_0_xxzzzz_1, g_0_xxxx_0_xxzzzzz_0, g_0_xxxx_0_xxzzzzz_1, g_0_xxxx_0_xyyyyy_1, g_0_xxxx_0_xyyyyyy_0, g_0_xxxx_0_xyyyyyy_1, g_0_xxxx_0_xyyyyyz_0, g_0_xxxx_0_xyyyyyz_1, g_0_xxxx_0_xyyyyz_1, g_0_xxxx_0_xyyyyzz_0, g_0_xxxx_0_xyyyyzz_1, g_0_xxxx_0_xyyyzz_1, g_0_xxxx_0_xyyyzzz_0, g_0_xxxx_0_xyyyzzz_1, g_0_xxxx_0_xyyzzz_1, g_0_xxxx_0_xyyzzzz_0, g_0_xxxx_0_xyyzzzz_1, g_0_xxxx_0_xyzzzz_1, g_0_xxxx_0_xyzzzzz_0, g_0_xxxx_0_xyzzzzz_1, g_0_xxxx_0_xzzzzz_1, g_0_xxxx_0_xzzzzzz_0, g_0_xxxx_0_xzzzzzz_1, g_0_xxxx_0_yyyyyy_1, g_0_xxxx_0_yyyyyyy_0, g_0_xxxx_0_yyyyyyy_1, g_0_xxxx_0_yyyyyyz_0, g_0_xxxx_0_yyyyyyz_1, g_0_xxxx_0_yyyyyz_1, g_0_xxxx_0_yyyyyzz_0, g_0_xxxx_0_yyyyyzz_1, g_0_xxxx_0_yyyyzz_1, g_0_xxxx_0_yyyyzzz_0, g_0_xxxx_0_yyyyzzz_1, g_0_xxxx_0_yyyzzz_1, g_0_xxxx_0_yyyzzzz_0, g_0_xxxx_0_yyyzzzz_1, g_0_xxxx_0_yyzzzz_1, g_0_xxxx_0_yyzzzzz_0, g_0_xxxx_0_yyzzzzz_1, g_0_xxxx_0_yzzzzz_1, g_0_xxxx_0_yzzzzzz_0, g_0_xxxx_0_yzzzzzz_1, g_0_xxxx_0_zzzzzz_1, g_0_xxxx_0_zzzzzzz_0, g_0_xxxx_0_zzzzzzz_1, g_0_xxxxz_0_xxxxxxx_0, g_0_xxxxz_0_xxxxxxy_0, g_0_xxxxz_0_xxxxxxz_0, g_0_xxxxz_0_xxxxxyy_0, g_0_xxxxz_0_xxxxxyz_0, g_0_xxxxz_0_xxxxxzz_0, g_0_xxxxz_0_xxxxyyy_0, g_0_xxxxz_0_xxxxyyz_0, g_0_xxxxz_0_xxxxyzz_0, g_0_xxxxz_0_xxxxzzz_0, g_0_xxxxz_0_xxxyyyy_0, g_0_xxxxz_0_xxxyyyz_0, g_0_xxxxz_0_xxxyyzz_0, g_0_xxxxz_0_xxxyzzz_0, g_0_xxxxz_0_xxxzzzz_0, g_0_xxxxz_0_xxyyyyy_0, g_0_xxxxz_0_xxyyyyz_0, g_0_xxxxz_0_xxyyyzz_0, g_0_xxxxz_0_xxyyzzz_0, g_0_xxxxz_0_xxyzzzz_0, g_0_xxxxz_0_xxzzzzz_0, g_0_xxxxz_0_xyyyyyy_0, g_0_xxxxz_0_xyyyyyz_0, g_0_xxxxz_0_xyyyyzz_0, g_0_xxxxz_0_xyyyzzz_0, g_0_xxxxz_0_xyyzzzz_0, g_0_xxxxz_0_xyzzzzz_0, g_0_xxxxz_0_xzzzzzz_0, g_0_xxxxz_0_yyyyyyy_0, g_0_xxxxz_0_yyyyyyz_0, g_0_xxxxz_0_yyyyyzz_0, g_0_xxxxz_0_yyyyzzz_0, g_0_xxxxz_0_yyyzzzz_0, g_0_xxxxz_0_yyzzzzz_0, g_0_xxxxz_0_yzzzzzz_0, g_0_xxxxz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxz_0_xxxxxxx_0[i] = g_0_xxxx_0_xxxxxxx_0[i] * pb_z + g_0_xxxx_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxxy_0[i] = g_0_xxxx_0_xxxxxxy_0[i] * pb_z + g_0_xxxx_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxxz_0[i] = g_0_xxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxz_0[i] * pb_z + g_0_xxxx_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxyy_0[i] = g_0_xxxx_0_xxxxxyy_0[i] * pb_z + g_0_xxxx_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxyz_0[i] = g_0_xxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyz_0[i] * pb_z + g_0_xxxx_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxzz_0[i] * pb_z + g_0_xxxx_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyyy_0[i] = g_0_xxxx_0_xxxxyyy_0[i] * pb_z + g_0_xxxx_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyyz_0[i] = g_0_xxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyz_0[i] * pb_z + g_0_xxxx_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyzz_0[i] * pb_z + g_0_xxxx_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxzzz_0[i] * pb_z + g_0_xxxx_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyyy_0[i] = g_0_xxxx_0_xxxyyyy_0[i] * pb_z + g_0_xxxx_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyyz_0[i] = g_0_xxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyz_0[i] * pb_z + g_0_xxxx_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyzz_0[i] * pb_z + g_0_xxxx_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzzz_0[i] * pb_z + g_0_xxxx_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxzzzz_0[i] * pb_z + g_0_xxxx_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyyy_0[i] = g_0_xxxx_0_xxyyyyy_0[i] * pb_z + g_0_xxxx_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyyz_0[i] = g_0_xxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyz_0[i] * pb_z + g_0_xxxx_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyzz_0[i] * pb_z + g_0_xxxx_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzzz_0[i] * pb_z + g_0_xxxx_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzzz_0[i] * pb_z + g_0_xxxx_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzzzzz_0[i] * pb_z + g_0_xxxx_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyyy_0[i] = g_0_xxxx_0_xyyyyyy_0[i] * pb_z + g_0_xxxx_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyyz_0[i] = g_0_xxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyz_0[i] * pb_z + g_0_xxxx_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyzz_0[i] * pb_z + g_0_xxxx_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyzzz_0[i] * pb_z + g_0_xxxx_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzzzz_0[i] * pb_z + g_0_xxxx_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzzzz_0[i] * pb_z + g_0_xxxx_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzzzzzz_0[i] * pb_z + g_0_xxxx_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyyy_0[i] = g_0_xxxx_0_yyyyyyy_0[i] * pb_z + g_0_xxxx_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyyz_0[i] = g_0_xxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyyz_0[i] * pb_z + g_0_xxxx_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyzz_0[i] * pb_z + g_0_xxxx_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyzzz_0[i] * pb_z + g_0_xxxx_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyzzzz_0[i] * pb_z + g_0_xxxx_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzzzzz_0[i] * pb_z + g_0_xxxx_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzzzzz_0[i] * pb_z + g_0_xxxx_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_zzzzzzz_0[i] * pb_z + g_0_xxxx_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 108-144 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_xxx_0_xxxxxxx_0, g_0_xxx_0_xxxxxxx_1, g_0_xxx_0_xxxxxxz_0, g_0_xxx_0_xxxxxxz_1, g_0_xxx_0_xxxxxzz_0, g_0_xxx_0_xxxxxzz_1, g_0_xxx_0_xxxxzzz_0, g_0_xxx_0_xxxxzzz_1, g_0_xxx_0_xxxzzzz_0, g_0_xxx_0_xxxzzzz_1, g_0_xxx_0_xxzzzzz_0, g_0_xxx_0_xxzzzzz_1, g_0_xxx_0_xzzzzzz_0, g_0_xxx_0_xzzzzzz_1, g_0_xxxy_0_xxxxxxx_0, g_0_xxxy_0_xxxxxxx_1, g_0_xxxy_0_xxxxxxz_0, g_0_xxxy_0_xxxxxxz_1, g_0_xxxy_0_xxxxxzz_0, g_0_xxxy_0_xxxxxzz_1, g_0_xxxy_0_xxxxzzz_0, g_0_xxxy_0_xxxxzzz_1, g_0_xxxy_0_xxxzzzz_0, g_0_xxxy_0_xxxzzzz_1, g_0_xxxy_0_xxzzzzz_0, g_0_xxxy_0_xxzzzzz_1, g_0_xxxy_0_xzzzzzz_0, g_0_xxxy_0_xzzzzzz_1, g_0_xxxyy_0_xxxxxxx_0, g_0_xxxyy_0_xxxxxxy_0, g_0_xxxyy_0_xxxxxxz_0, g_0_xxxyy_0_xxxxxyy_0, g_0_xxxyy_0_xxxxxyz_0, g_0_xxxyy_0_xxxxxzz_0, g_0_xxxyy_0_xxxxyyy_0, g_0_xxxyy_0_xxxxyyz_0, g_0_xxxyy_0_xxxxyzz_0, g_0_xxxyy_0_xxxxzzz_0, g_0_xxxyy_0_xxxyyyy_0, g_0_xxxyy_0_xxxyyyz_0, g_0_xxxyy_0_xxxyyzz_0, g_0_xxxyy_0_xxxyzzz_0, g_0_xxxyy_0_xxxzzzz_0, g_0_xxxyy_0_xxyyyyy_0, g_0_xxxyy_0_xxyyyyz_0, g_0_xxxyy_0_xxyyyzz_0, g_0_xxxyy_0_xxyyzzz_0, g_0_xxxyy_0_xxyzzzz_0, g_0_xxxyy_0_xxzzzzz_0, g_0_xxxyy_0_xyyyyyy_0, g_0_xxxyy_0_xyyyyyz_0, g_0_xxxyy_0_xyyyyzz_0, g_0_xxxyy_0_xyyyzzz_0, g_0_xxxyy_0_xyyzzzz_0, g_0_xxxyy_0_xyzzzzz_0, g_0_xxxyy_0_xzzzzzz_0, g_0_xxxyy_0_yyyyyyy_0, g_0_xxxyy_0_yyyyyyz_0, g_0_xxxyy_0_yyyyyzz_0, g_0_xxxyy_0_yyyyzzz_0, g_0_xxxyy_0_yyyzzzz_0, g_0_xxxyy_0_yyzzzzz_0, g_0_xxxyy_0_yzzzzzz_0, g_0_xxxyy_0_zzzzzzz_0, g_0_xxyy_0_xxxxxxy_0, g_0_xxyy_0_xxxxxxy_1, g_0_xxyy_0_xxxxxy_1, g_0_xxyy_0_xxxxxyy_0, g_0_xxyy_0_xxxxxyy_1, g_0_xxyy_0_xxxxxyz_0, g_0_xxyy_0_xxxxxyz_1, g_0_xxyy_0_xxxxyy_1, g_0_xxyy_0_xxxxyyy_0, g_0_xxyy_0_xxxxyyy_1, g_0_xxyy_0_xxxxyyz_0, g_0_xxyy_0_xxxxyyz_1, g_0_xxyy_0_xxxxyz_1, g_0_xxyy_0_xxxxyzz_0, g_0_xxyy_0_xxxxyzz_1, g_0_xxyy_0_xxxyyy_1, g_0_xxyy_0_xxxyyyy_0, g_0_xxyy_0_xxxyyyy_1, g_0_xxyy_0_xxxyyyz_0, g_0_xxyy_0_xxxyyyz_1, g_0_xxyy_0_xxxyyz_1, g_0_xxyy_0_xxxyyzz_0, g_0_xxyy_0_xxxyyzz_1, g_0_xxyy_0_xxxyzz_1, g_0_xxyy_0_xxxyzzz_0, g_0_xxyy_0_xxxyzzz_1, g_0_xxyy_0_xxyyyy_1, g_0_xxyy_0_xxyyyyy_0, g_0_xxyy_0_xxyyyyy_1, g_0_xxyy_0_xxyyyyz_0, g_0_xxyy_0_xxyyyyz_1, g_0_xxyy_0_xxyyyz_1, g_0_xxyy_0_xxyyyzz_0, g_0_xxyy_0_xxyyyzz_1, g_0_xxyy_0_xxyyzz_1, g_0_xxyy_0_xxyyzzz_0, g_0_xxyy_0_xxyyzzz_1, g_0_xxyy_0_xxyzzz_1, g_0_xxyy_0_xxyzzzz_0, g_0_xxyy_0_xxyzzzz_1, g_0_xxyy_0_xyyyyy_1, g_0_xxyy_0_xyyyyyy_0, g_0_xxyy_0_xyyyyyy_1, g_0_xxyy_0_xyyyyyz_0, g_0_xxyy_0_xyyyyyz_1, g_0_xxyy_0_xyyyyz_1, g_0_xxyy_0_xyyyyzz_0, g_0_xxyy_0_xyyyyzz_1, g_0_xxyy_0_xyyyzz_1, g_0_xxyy_0_xyyyzzz_0, g_0_xxyy_0_xyyyzzz_1, g_0_xxyy_0_xyyzzz_1, g_0_xxyy_0_xyyzzzz_0, g_0_xxyy_0_xyyzzzz_1, g_0_xxyy_0_xyzzzz_1, g_0_xxyy_0_xyzzzzz_0, g_0_xxyy_0_xyzzzzz_1, g_0_xxyy_0_yyyyyy_1, g_0_xxyy_0_yyyyyyy_0, g_0_xxyy_0_yyyyyyy_1, g_0_xxyy_0_yyyyyyz_0, g_0_xxyy_0_yyyyyyz_1, g_0_xxyy_0_yyyyyz_1, g_0_xxyy_0_yyyyyzz_0, g_0_xxyy_0_yyyyyzz_1, g_0_xxyy_0_yyyyzz_1, g_0_xxyy_0_yyyyzzz_0, g_0_xxyy_0_yyyyzzz_1, g_0_xxyy_0_yyyzzz_1, g_0_xxyy_0_yyyzzzz_0, g_0_xxyy_0_yyyzzzz_1, g_0_xxyy_0_yyzzzz_1, g_0_xxyy_0_yyzzzzz_0, g_0_xxyy_0_yyzzzzz_1, g_0_xxyy_0_yzzzzz_1, g_0_xxyy_0_yzzzzzz_0, g_0_xxyy_0_yzzzzzz_1, g_0_xxyy_0_zzzzzzz_0, g_0_xxyy_0_zzzzzzz_1, g_0_xyy_0_xxxxxxy_0, g_0_xyy_0_xxxxxxy_1, g_0_xyy_0_xxxxxyy_0, g_0_xyy_0_xxxxxyy_1, g_0_xyy_0_xxxxxyz_0, g_0_xyy_0_xxxxxyz_1, g_0_xyy_0_xxxxyyy_0, g_0_xyy_0_xxxxyyy_1, g_0_xyy_0_xxxxyyz_0, g_0_xyy_0_xxxxyyz_1, g_0_xyy_0_xxxxyzz_0, g_0_xyy_0_xxxxyzz_1, g_0_xyy_0_xxxyyyy_0, g_0_xyy_0_xxxyyyy_1, g_0_xyy_0_xxxyyyz_0, g_0_xyy_0_xxxyyyz_1, g_0_xyy_0_xxxyyzz_0, g_0_xyy_0_xxxyyzz_1, g_0_xyy_0_xxxyzzz_0, g_0_xyy_0_xxxyzzz_1, g_0_xyy_0_xxyyyyy_0, g_0_xyy_0_xxyyyyy_1, g_0_xyy_0_xxyyyyz_0, g_0_xyy_0_xxyyyyz_1, g_0_xyy_0_xxyyyzz_0, g_0_xyy_0_xxyyyzz_1, g_0_xyy_0_xxyyzzz_0, g_0_xyy_0_xxyyzzz_1, g_0_xyy_0_xxyzzzz_0, g_0_xyy_0_xxyzzzz_1, g_0_xyy_0_xyyyyyy_0, g_0_xyy_0_xyyyyyy_1, g_0_xyy_0_xyyyyyz_0, g_0_xyy_0_xyyyyyz_1, g_0_xyy_0_xyyyyzz_0, g_0_xyy_0_xyyyyzz_1, g_0_xyy_0_xyyyzzz_0, g_0_xyy_0_xyyyzzz_1, g_0_xyy_0_xyyzzzz_0, g_0_xyy_0_xyyzzzz_1, g_0_xyy_0_xyzzzzz_0, g_0_xyy_0_xyzzzzz_1, g_0_xyy_0_yyyyyyy_0, g_0_xyy_0_yyyyyyy_1, g_0_xyy_0_yyyyyyz_0, g_0_xyy_0_yyyyyyz_1, g_0_xyy_0_yyyyyzz_0, g_0_xyy_0_yyyyyzz_1, g_0_xyy_0_yyyyzzz_0, g_0_xyy_0_yyyyzzz_1, g_0_xyy_0_yyyzzzz_0, g_0_xyy_0_yyyzzzz_1, g_0_xyy_0_yyzzzzz_0, g_0_xyy_0_yyzzzzz_1, g_0_xyy_0_yzzzzzz_0, g_0_xyy_0_yzzzzzz_1, g_0_xyy_0_zzzzzzz_0, g_0_xyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyy_0_xxxxxxx_0[i] = g_0_xxx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxxx_0[i] * pb_y + g_0_xxxy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxxxy_0[i] = 2.0 * g_0_xyy_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxxy_0[i] * pb_x + g_0_xxyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxxz_0[i] = g_0_xxx_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxxz_0[i] * pb_y + g_0_xxxy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxxyy_0[i] = 2.0 * g_0_xyy_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxyy_0[i] * pb_x + g_0_xxyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxyz_0[i] = 2.0 * g_0_xyy_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxyz_0[i] * pb_x + g_0_xxyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxzz_0[i] = g_0_xxx_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxzz_0[i] * pb_y + g_0_xxxy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxyyy_0[i] = 2.0 * g_0_xyy_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyyy_0[i] * pb_x + g_0_xxyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxyyz_0[i] = 2.0 * g_0_xyy_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyyz_0[i] * pb_x + g_0_xxyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxyzz_0[i] = 2.0 * g_0_xyy_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyzz_0[i] * pb_x + g_0_xxyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxzzz_0[i] = g_0_xxx_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxzzz_0[i] * pb_y + g_0_xxxy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxyyyy_0[i] = 2.0 * g_0_xyy_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyyy_0[i] * pb_x + g_0_xxyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyyyz_0[i] = 2.0 * g_0_xyy_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyyz_0[i] * pb_x + g_0_xxyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyyzz_0[i] = 2.0 * g_0_xyy_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyzz_0[i] * pb_x + g_0_xxyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyzzz_0[i] = 2.0 * g_0_xyy_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyzzz_0[i] * pb_x + g_0_xxyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxzzzz_0[i] = g_0_xxx_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxzzzz_0[i] * pb_y + g_0_xxxy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxyyyyy_0[i] = 2.0 * g_0_xyy_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyyy_0[i] * pb_x + g_0_xxyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyyyz_0[i] = 2.0 * g_0_xyy_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyyz_0[i] * pb_x + g_0_xxyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyyzz_0[i] = 2.0 * g_0_xyy_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyzz_0[i] * pb_x + g_0_xxyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyzzz_0[i] = 2.0 * g_0_xyy_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyzzz_0[i] * pb_x + g_0_xxyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyzzzz_0[i] = 2.0 * g_0_xyy_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyzzzz_0[i] * pb_x + g_0_xxyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxzzzzz_0[i] = g_0_xxx_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxzzzzz_0[i] * pb_y + g_0_xxxy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xyyyyyy_0[i] = 2.0 * g_0_xyy_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyyy_0[i] * pb_x + g_0_xxyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyyyz_0[i] = 2.0 * g_0_xyy_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyyz_0[i] * pb_x + g_0_xxyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyyzz_0[i] = 2.0 * g_0_xyy_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyzz_0[i] * pb_x + g_0_xxyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyzzz_0[i] = 2.0 * g_0_xyy_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyzzz_0[i] * pb_x + g_0_xxyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyzzzz_0[i] = 2.0 * g_0_xyy_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyzzzz_0[i] * pb_x + g_0_xxyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyzzzzz_0[i] = 2.0 * g_0_xyy_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyzzzzz_0[i] * pb_x + g_0_xxyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xzzzzzz_0[i] = g_0_xxx_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xzzzzzz_0[i] * pb_y + g_0_xxxy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_yyyyyyy_0[i] = 2.0 * g_0_xyy_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyyy_0[i] * pb_x + g_0_xxyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyyyz_0[i] = 2.0 * g_0_xyy_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyyz_0[i] * pb_x + g_0_xxyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyyzz_0[i] = 2.0 * g_0_xyy_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyzz_0[i] * pb_x + g_0_xxyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyzzz_0[i] = 2.0 * g_0_xyy_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyzzz_0[i] * pb_x + g_0_xxyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyzzzz_0[i] = 2.0 * g_0_xyy_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyzzzz_0[i] * pb_x + g_0_xxyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyzzzzz_0[i] = 2.0 * g_0_xyy_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzzzzz_0[i] * pb_x + g_0_xxyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yzzzzzz_0[i] = 2.0 * g_0_xyy_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzzzzz_0[i] * pb_x + g_0_xxyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_zzzzzzz_0[i] = 2.0 * g_0_xyy_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_zzzzzzz_0[i] * pb_x + g_0_xxyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 144-180 components of targeted buffer : SHSK

    auto g_0_xxxyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 144);

    auto g_0_xxxyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 145);

    auto g_0_xxxyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 146);

    auto g_0_xxxyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 147);

    auto g_0_xxxyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 148);

    auto g_0_xxxyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 149);

    auto g_0_xxxyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 150);

    auto g_0_xxxyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 151);

    auto g_0_xxxyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 152);

    auto g_0_xxxyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 153);

    auto g_0_xxxyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 154);

    auto g_0_xxxyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 155);

    auto g_0_xxxyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 156);

    auto g_0_xxxyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 157);

    auto g_0_xxxyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 158);

    auto g_0_xxxyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 159);

    auto g_0_xxxyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 160);

    auto g_0_xxxyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 161);

    auto g_0_xxxyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 162);

    auto g_0_xxxyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 163);

    auto g_0_xxxyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 164);

    auto g_0_xxxyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 165);

    auto g_0_xxxyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 166);

    auto g_0_xxxyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 167);

    auto g_0_xxxyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 168);

    auto g_0_xxxyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 169);

    auto g_0_xxxyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 170);

    auto g_0_xxxyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 171);

    auto g_0_xxxyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 172);

    auto g_0_xxxyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 173);

    auto g_0_xxxyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 174);

    auto g_0_xxxyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 175);

    auto g_0_xxxyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 176);

    auto g_0_xxxyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 177);

    auto g_0_xxxyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 178);

    auto g_0_xxxyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 179);

    #pragma omp simd aligned(g_0_xxxy_0_xxxxxxy_0, g_0_xxxy_0_xxxxxxy_1, g_0_xxxy_0_xxxxxyy_0, g_0_xxxy_0_xxxxxyy_1, g_0_xxxy_0_xxxxyyy_0, g_0_xxxy_0_xxxxyyy_1, g_0_xxxy_0_xxxyyyy_0, g_0_xxxy_0_xxxyyyy_1, g_0_xxxy_0_xxyyyyy_0, g_0_xxxy_0_xxyyyyy_1, g_0_xxxy_0_xyyyyyy_0, g_0_xxxy_0_xyyyyyy_1, g_0_xxxy_0_yyyyyyy_0, g_0_xxxy_0_yyyyyyy_1, g_0_xxxyz_0_xxxxxxx_0, g_0_xxxyz_0_xxxxxxy_0, g_0_xxxyz_0_xxxxxxz_0, g_0_xxxyz_0_xxxxxyy_0, g_0_xxxyz_0_xxxxxyz_0, g_0_xxxyz_0_xxxxxzz_0, g_0_xxxyz_0_xxxxyyy_0, g_0_xxxyz_0_xxxxyyz_0, g_0_xxxyz_0_xxxxyzz_0, g_0_xxxyz_0_xxxxzzz_0, g_0_xxxyz_0_xxxyyyy_0, g_0_xxxyz_0_xxxyyyz_0, g_0_xxxyz_0_xxxyyzz_0, g_0_xxxyz_0_xxxyzzz_0, g_0_xxxyz_0_xxxzzzz_0, g_0_xxxyz_0_xxyyyyy_0, g_0_xxxyz_0_xxyyyyz_0, g_0_xxxyz_0_xxyyyzz_0, g_0_xxxyz_0_xxyyzzz_0, g_0_xxxyz_0_xxyzzzz_0, g_0_xxxyz_0_xxzzzzz_0, g_0_xxxyz_0_xyyyyyy_0, g_0_xxxyz_0_xyyyyyz_0, g_0_xxxyz_0_xyyyyzz_0, g_0_xxxyz_0_xyyyzzz_0, g_0_xxxyz_0_xyyzzzz_0, g_0_xxxyz_0_xyzzzzz_0, g_0_xxxyz_0_xzzzzzz_0, g_0_xxxyz_0_yyyyyyy_0, g_0_xxxyz_0_yyyyyyz_0, g_0_xxxyz_0_yyyyyzz_0, g_0_xxxyz_0_yyyyzzz_0, g_0_xxxyz_0_yyyzzzz_0, g_0_xxxyz_0_yyzzzzz_0, g_0_xxxyz_0_yzzzzzz_0, g_0_xxxyz_0_zzzzzzz_0, g_0_xxxz_0_xxxxxxx_0, g_0_xxxz_0_xxxxxxx_1, g_0_xxxz_0_xxxxxxz_0, g_0_xxxz_0_xxxxxxz_1, g_0_xxxz_0_xxxxxyz_0, g_0_xxxz_0_xxxxxyz_1, g_0_xxxz_0_xxxxxz_1, g_0_xxxz_0_xxxxxzz_0, g_0_xxxz_0_xxxxxzz_1, g_0_xxxz_0_xxxxyyz_0, g_0_xxxz_0_xxxxyyz_1, g_0_xxxz_0_xxxxyz_1, g_0_xxxz_0_xxxxyzz_0, g_0_xxxz_0_xxxxyzz_1, g_0_xxxz_0_xxxxzz_1, g_0_xxxz_0_xxxxzzz_0, g_0_xxxz_0_xxxxzzz_1, g_0_xxxz_0_xxxyyyz_0, g_0_xxxz_0_xxxyyyz_1, g_0_xxxz_0_xxxyyz_1, g_0_xxxz_0_xxxyyzz_0, g_0_xxxz_0_xxxyyzz_1, g_0_xxxz_0_xxxyzz_1, g_0_xxxz_0_xxxyzzz_0, g_0_xxxz_0_xxxyzzz_1, g_0_xxxz_0_xxxzzz_1, g_0_xxxz_0_xxxzzzz_0, g_0_xxxz_0_xxxzzzz_1, g_0_xxxz_0_xxyyyyz_0, g_0_xxxz_0_xxyyyyz_1, g_0_xxxz_0_xxyyyz_1, g_0_xxxz_0_xxyyyzz_0, g_0_xxxz_0_xxyyyzz_1, g_0_xxxz_0_xxyyzz_1, g_0_xxxz_0_xxyyzzz_0, g_0_xxxz_0_xxyyzzz_1, g_0_xxxz_0_xxyzzz_1, g_0_xxxz_0_xxyzzzz_0, g_0_xxxz_0_xxyzzzz_1, g_0_xxxz_0_xxzzzz_1, g_0_xxxz_0_xxzzzzz_0, g_0_xxxz_0_xxzzzzz_1, g_0_xxxz_0_xyyyyyz_0, g_0_xxxz_0_xyyyyyz_1, g_0_xxxz_0_xyyyyz_1, g_0_xxxz_0_xyyyyzz_0, g_0_xxxz_0_xyyyyzz_1, g_0_xxxz_0_xyyyzz_1, g_0_xxxz_0_xyyyzzz_0, g_0_xxxz_0_xyyyzzz_1, g_0_xxxz_0_xyyzzz_1, g_0_xxxz_0_xyyzzzz_0, g_0_xxxz_0_xyyzzzz_1, g_0_xxxz_0_xyzzzz_1, g_0_xxxz_0_xyzzzzz_0, g_0_xxxz_0_xyzzzzz_1, g_0_xxxz_0_xzzzzz_1, g_0_xxxz_0_xzzzzzz_0, g_0_xxxz_0_xzzzzzz_1, g_0_xxxz_0_yyyyyyz_0, g_0_xxxz_0_yyyyyyz_1, g_0_xxxz_0_yyyyyz_1, g_0_xxxz_0_yyyyyzz_0, g_0_xxxz_0_yyyyyzz_1, g_0_xxxz_0_yyyyzz_1, g_0_xxxz_0_yyyyzzz_0, g_0_xxxz_0_yyyyzzz_1, g_0_xxxz_0_yyyzzz_1, g_0_xxxz_0_yyyzzzz_0, g_0_xxxz_0_yyyzzzz_1, g_0_xxxz_0_yyzzzz_1, g_0_xxxz_0_yyzzzzz_0, g_0_xxxz_0_yyzzzzz_1, g_0_xxxz_0_yzzzzz_1, g_0_xxxz_0_yzzzzzz_0, g_0_xxxz_0_yzzzzzz_1, g_0_xxxz_0_zzzzzz_1, g_0_xxxz_0_zzzzzzz_0, g_0_xxxz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyz_0_xxxxxxx_0[i] = g_0_xxxz_0_xxxxxxx_0[i] * pb_y + g_0_xxxz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxxy_0[i] = g_0_xxxy_0_xxxxxxy_0[i] * pb_z + g_0_xxxy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxxxz_0[i] = g_0_xxxz_0_xxxxxxz_0[i] * pb_y + g_0_xxxz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxyy_0[i] = g_0_xxxy_0_xxxxxyy_0[i] * pb_z + g_0_xxxy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxxyz_0[i] = g_0_xxxz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxxyz_0[i] * pb_y + g_0_xxxz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxzz_0[i] = g_0_xxxz_0_xxxxxzz_0[i] * pb_y + g_0_xxxz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxyyy_0[i] = g_0_xxxy_0_xxxxyyy_0[i] * pb_z + g_0_xxxy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxyyz_0[i] * pb_y + g_0_xxxz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxyzz_0[i] = g_0_xxxz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxyzz_0[i] * pb_y + g_0_xxxz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxzzz_0[i] = g_0_xxxz_0_xxxxzzz_0[i] * pb_y + g_0_xxxz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyyyy_0[i] = g_0_xxxy_0_xxxyyyy_0[i] * pb_z + g_0_xxxy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyyyz_0[i] * pb_y + g_0_xxxz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyyzz_0[i] * pb_y + g_0_xxxz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyzzz_0[i] = g_0_xxxz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyzzz_0[i] * pb_y + g_0_xxxz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxzzzz_0[i] = g_0_xxxz_0_xxxzzzz_0[i] * pb_y + g_0_xxxz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyyyy_0[i] = g_0_xxxy_0_xxyyyyy_0[i] * pb_z + g_0_xxxy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyyyz_0[i] * pb_y + g_0_xxxz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyyzz_0[i] * pb_y + g_0_xxxz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyzzz_0[i] * pb_y + g_0_xxxz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyzzzz_0[i] = g_0_xxxz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyzzzz_0[i] * pb_y + g_0_xxxz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxzzzzz_0[i] = g_0_xxxz_0_xxzzzzz_0[i] * pb_y + g_0_xxxz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyyyy_0[i] = g_0_xxxy_0_xyyyyyy_0[i] * pb_z + g_0_xxxy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyyyz_0[i] * pb_y + g_0_xxxz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyyzz_0[i] * pb_y + g_0_xxxz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyzzz_0[i] * pb_y + g_0_xxxz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyzzzz_0[i] * pb_y + g_0_xxxz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyzzzzz_0[i] = g_0_xxxz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyzzzzz_0[i] * pb_y + g_0_xxxz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xzzzzzz_0[i] = g_0_xxxz_0_xzzzzzz_0[i] * pb_y + g_0_xxxz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyyyy_0[i] = g_0_xxxy_0_yyyyyyy_0[i] * pb_z + g_0_xxxy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyyyz_0[i] * pb_y + g_0_xxxz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyyzz_0[i] * pb_y + g_0_xxxz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyzzz_0[i] * pb_y + g_0_xxxz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyzzzz_0[i] * pb_y + g_0_xxxz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyzzzzz_0[i] * pb_y + g_0_xxxz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yzzzzzz_0[i] = g_0_xxxz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yzzzzzz_0[i] * pb_y + g_0_xxxz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_zzzzzzz_0[i] = g_0_xxxz_0_zzzzzzz_0[i] * pb_y + g_0_xxxz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 180-216 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_xxx_0_xxxxxxx_0, g_0_xxx_0_xxxxxxx_1, g_0_xxx_0_xxxxxxy_0, g_0_xxx_0_xxxxxxy_1, g_0_xxx_0_xxxxxyy_0, g_0_xxx_0_xxxxxyy_1, g_0_xxx_0_xxxxyyy_0, g_0_xxx_0_xxxxyyy_1, g_0_xxx_0_xxxyyyy_0, g_0_xxx_0_xxxyyyy_1, g_0_xxx_0_xxyyyyy_0, g_0_xxx_0_xxyyyyy_1, g_0_xxx_0_xyyyyyy_0, g_0_xxx_0_xyyyyyy_1, g_0_xxxz_0_xxxxxxx_0, g_0_xxxz_0_xxxxxxx_1, g_0_xxxz_0_xxxxxxy_0, g_0_xxxz_0_xxxxxxy_1, g_0_xxxz_0_xxxxxyy_0, g_0_xxxz_0_xxxxxyy_1, g_0_xxxz_0_xxxxyyy_0, g_0_xxxz_0_xxxxyyy_1, g_0_xxxz_0_xxxyyyy_0, g_0_xxxz_0_xxxyyyy_1, g_0_xxxz_0_xxyyyyy_0, g_0_xxxz_0_xxyyyyy_1, g_0_xxxz_0_xyyyyyy_0, g_0_xxxz_0_xyyyyyy_1, g_0_xxxzz_0_xxxxxxx_0, g_0_xxxzz_0_xxxxxxy_0, g_0_xxxzz_0_xxxxxxz_0, g_0_xxxzz_0_xxxxxyy_0, g_0_xxxzz_0_xxxxxyz_0, g_0_xxxzz_0_xxxxxzz_0, g_0_xxxzz_0_xxxxyyy_0, g_0_xxxzz_0_xxxxyyz_0, g_0_xxxzz_0_xxxxyzz_0, g_0_xxxzz_0_xxxxzzz_0, g_0_xxxzz_0_xxxyyyy_0, g_0_xxxzz_0_xxxyyyz_0, g_0_xxxzz_0_xxxyyzz_0, g_0_xxxzz_0_xxxyzzz_0, g_0_xxxzz_0_xxxzzzz_0, g_0_xxxzz_0_xxyyyyy_0, g_0_xxxzz_0_xxyyyyz_0, g_0_xxxzz_0_xxyyyzz_0, g_0_xxxzz_0_xxyyzzz_0, g_0_xxxzz_0_xxyzzzz_0, g_0_xxxzz_0_xxzzzzz_0, g_0_xxxzz_0_xyyyyyy_0, g_0_xxxzz_0_xyyyyyz_0, g_0_xxxzz_0_xyyyyzz_0, g_0_xxxzz_0_xyyyzzz_0, g_0_xxxzz_0_xyyzzzz_0, g_0_xxxzz_0_xyzzzzz_0, g_0_xxxzz_0_xzzzzzz_0, g_0_xxxzz_0_yyyyyyy_0, g_0_xxxzz_0_yyyyyyz_0, g_0_xxxzz_0_yyyyyzz_0, g_0_xxxzz_0_yyyyzzz_0, g_0_xxxzz_0_yyyzzzz_0, g_0_xxxzz_0_yyzzzzz_0, g_0_xxxzz_0_yzzzzzz_0, g_0_xxxzz_0_zzzzzzz_0, g_0_xxzz_0_xxxxxxz_0, g_0_xxzz_0_xxxxxxz_1, g_0_xxzz_0_xxxxxyz_0, g_0_xxzz_0_xxxxxyz_1, g_0_xxzz_0_xxxxxz_1, g_0_xxzz_0_xxxxxzz_0, g_0_xxzz_0_xxxxxzz_1, g_0_xxzz_0_xxxxyyz_0, g_0_xxzz_0_xxxxyyz_1, g_0_xxzz_0_xxxxyz_1, g_0_xxzz_0_xxxxyzz_0, g_0_xxzz_0_xxxxyzz_1, g_0_xxzz_0_xxxxzz_1, g_0_xxzz_0_xxxxzzz_0, g_0_xxzz_0_xxxxzzz_1, g_0_xxzz_0_xxxyyyz_0, g_0_xxzz_0_xxxyyyz_1, g_0_xxzz_0_xxxyyz_1, g_0_xxzz_0_xxxyyzz_0, g_0_xxzz_0_xxxyyzz_1, g_0_xxzz_0_xxxyzz_1, g_0_xxzz_0_xxxyzzz_0, g_0_xxzz_0_xxxyzzz_1, g_0_xxzz_0_xxxzzz_1, g_0_xxzz_0_xxxzzzz_0, g_0_xxzz_0_xxxzzzz_1, g_0_xxzz_0_xxyyyyz_0, g_0_xxzz_0_xxyyyyz_1, g_0_xxzz_0_xxyyyz_1, g_0_xxzz_0_xxyyyzz_0, g_0_xxzz_0_xxyyyzz_1, g_0_xxzz_0_xxyyzz_1, g_0_xxzz_0_xxyyzzz_0, g_0_xxzz_0_xxyyzzz_1, g_0_xxzz_0_xxyzzz_1, g_0_xxzz_0_xxyzzzz_0, g_0_xxzz_0_xxyzzzz_1, g_0_xxzz_0_xxzzzz_1, g_0_xxzz_0_xxzzzzz_0, g_0_xxzz_0_xxzzzzz_1, g_0_xxzz_0_xyyyyyz_0, g_0_xxzz_0_xyyyyyz_1, g_0_xxzz_0_xyyyyz_1, g_0_xxzz_0_xyyyyzz_0, g_0_xxzz_0_xyyyyzz_1, g_0_xxzz_0_xyyyzz_1, g_0_xxzz_0_xyyyzzz_0, g_0_xxzz_0_xyyyzzz_1, g_0_xxzz_0_xyyzzz_1, g_0_xxzz_0_xyyzzzz_0, g_0_xxzz_0_xyyzzzz_1, g_0_xxzz_0_xyzzzz_1, g_0_xxzz_0_xyzzzzz_0, g_0_xxzz_0_xyzzzzz_1, g_0_xxzz_0_xzzzzz_1, g_0_xxzz_0_xzzzzzz_0, g_0_xxzz_0_xzzzzzz_1, g_0_xxzz_0_yyyyyyy_0, g_0_xxzz_0_yyyyyyy_1, g_0_xxzz_0_yyyyyyz_0, g_0_xxzz_0_yyyyyyz_1, g_0_xxzz_0_yyyyyz_1, g_0_xxzz_0_yyyyyzz_0, g_0_xxzz_0_yyyyyzz_1, g_0_xxzz_0_yyyyzz_1, g_0_xxzz_0_yyyyzzz_0, g_0_xxzz_0_yyyyzzz_1, g_0_xxzz_0_yyyzzz_1, g_0_xxzz_0_yyyzzzz_0, g_0_xxzz_0_yyyzzzz_1, g_0_xxzz_0_yyzzzz_1, g_0_xxzz_0_yyzzzzz_0, g_0_xxzz_0_yyzzzzz_1, g_0_xxzz_0_yzzzzz_1, g_0_xxzz_0_yzzzzzz_0, g_0_xxzz_0_yzzzzzz_1, g_0_xxzz_0_zzzzzz_1, g_0_xxzz_0_zzzzzzz_0, g_0_xxzz_0_zzzzzzz_1, g_0_xzz_0_xxxxxxz_0, g_0_xzz_0_xxxxxxz_1, g_0_xzz_0_xxxxxyz_0, g_0_xzz_0_xxxxxyz_1, g_0_xzz_0_xxxxxzz_0, g_0_xzz_0_xxxxxzz_1, g_0_xzz_0_xxxxyyz_0, g_0_xzz_0_xxxxyyz_1, g_0_xzz_0_xxxxyzz_0, g_0_xzz_0_xxxxyzz_1, g_0_xzz_0_xxxxzzz_0, g_0_xzz_0_xxxxzzz_1, g_0_xzz_0_xxxyyyz_0, g_0_xzz_0_xxxyyyz_1, g_0_xzz_0_xxxyyzz_0, g_0_xzz_0_xxxyyzz_1, g_0_xzz_0_xxxyzzz_0, g_0_xzz_0_xxxyzzz_1, g_0_xzz_0_xxxzzzz_0, g_0_xzz_0_xxxzzzz_1, g_0_xzz_0_xxyyyyz_0, g_0_xzz_0_xxyyyyz_1, g_0_xzz_0_xxyyyzz_0, g_0_xzz_0_xxyyyzz_1, g_0_xzz_0_xxyyzzz_0, g_0_xzz_0_xxyyzzz_1, g_0_xzz_0_xxyzzzz_0, g_0_xzz_0_xxyzzzz_1, g_0_xzz_0_xxzzzzz_0, g_0_xzz_0_xxzzzzz_1, g_0_xzz_0_xyyyyyz_0, g_0_xzz_0_xyyyyyz_1, g_0_xzz_0_xyyyyzz_0, g_0_xzz_0_xyyyyzz_1, g_0_xzz_0_xyyyzzz_0, g_0_xzz_0_xyyyzzz_1, g_0_xzz_0_xyyzzzz_0, g_0_xzz_0_xyyzzzz_1, g_0_xzz_0_xyzzzzz_0, g_0_xzz_0_xyzzzzz_1, g_0_xzz_0_xzzzzzz_0, g_0_xzz_0_xzzzzzz_1, g_0_xzz_0_yyyyyyy_0, g_0_xzz_0_yyyyyyy_1, g_0_xzz_0_yyyyyyz_0, g_0_xzz_0_yyyyyyz_1, g_0_xzz_0_yyyyyzz_0, g_0_xzz_0_yyyyyzz_1, g_0_xzz_0_yyyyzzz_0, g_0_xzz_0_yyyyzzz_1, g_0_xzz_0_yyyzzzz_0, g_0_xzz_0_yyyzzzz_1, g_0_xzz_0_yyzzzzz_0, g_0_xzz_0_yyzzzzz_1, g_0_xzz_0_yzzzzzz_0, g_0_xzz_0_yzzzzzz_1, g_0_xzz_0_zzzzzzz_0, g_0_xzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzz_0_xxxxxxx_0[i] = g_0_xxx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxxx_0[i] * pb_z + g_0_xxxz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxxy_0[i] = g_0_xxx_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxxy_0[i] * pb_z + g_0_xxxz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxxz_0[i] = 2.0 * g_0_xzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxxz_0[i] * pb_x + g_0_xxzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxxyy_0[i] = g_0_xxx_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxyy_0[i] * pb_z + g_0_xxxz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxyz_0[i] = 2.0 * g_0_xzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxyz_0[i] * pb_x + g_0_xxzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxxzz_0[i] = 2.0 * g_0_xzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxzz_0[i] * pb_x + g_0_xxzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxyyy_0[i] = g_0_xxx_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxyyy_0[i] * pb_z + g_0_xxxz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxyyz_0[i] = 2.0 * g_0_xzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyyz_0[i] * pb_x + g_0_xxzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxyzz_0[i] = 2.0 * g_0_xzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyzz_0[i] * pb_x + g_0_xxzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxzzz_0[i] = 2.0 * g_0_xzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxzzz_0[i] * pb_x + g_0_xxzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyyyy_0[i] = g_0_xxx_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxyyyy_0[i] * pb_z + g_0_xxxz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxyyyz_0[i] = 2.0 * g_0_xzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyyz_0[i] * pb_x + g_0_xxzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyyzz_0[i] = 2.0 * g_0_xzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyzz_0[i] * pb_x + g_0_xxzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyzzz_0[i] = 2.0 * g_0_xzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyzzz_0[i] * pb_x + g_0_xxzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxzzzz_0[i] = 2.0 * g_0_xzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxzzzz_0[i] * pb_x + g_0_xxzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyyyy_0[i] = g_0_xxx_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxyyyyy_0[i] * pb_z + g_0_xxxz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxyyyyz_0[i] = 2.0 * g_0_xzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyyz_0[i] * pb_x + g_0_xxzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyyzz_0[i] = 2.0 * g_0_xzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyzz_0[i] * pb_x + g_0_xxzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyzzz_0[i] = 2.0 * g_0_xzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyzzz_0[i] * pb_x + g_0_xxzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyzzzz_0[i] = 2.0 * g_0_xzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyzzzz_0[i] * pb_x + g_0_xxzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxzzzzz_0[i] = 2.0 * g_0_xzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxzzzzz_0[i] * pb_x + g_0_xxzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyyyy_0[i] = g_0_xxx_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xyyyyyy_0[i] * pb_z + g_0_xxxz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xyyyyyz_0[i] = 2.0 * g_0_xzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyyz_0[i] * pb_x + g_0_xxzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyyzz_0[i] = 2.0 * g_0_xzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyzz_0[i] * pb_x + g_0_xxzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyzzz_0[i] = 2.0 * g_0_xzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyzzz_0[i] * pb_x + g_0_xxzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyzzzz_0[i] = 2.0 * g_0_xzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyzzzz_0[i] * pb_x + g_0_xxzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyzzzzz_0[i] = 2.0 * g_0_xzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyzzzzz_0[i] * pb_x + g_0_xxzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xzzzzzz_0[i] = 2.0 * g_0_xzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xzzzzzz_0[i] * pb_x + g_0_xxzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyyy_0[i] = 2.0 * g_0_xzz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyyy_0[i] * pb_x + g_0_xxzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyyz_0[i] = 2.0 * g_0_xzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyyz_0[i] * pb_x + g_0_xxzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyzz_0[i] = 2.0 * g_0_xzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyzz_0[i] * pb_x + g_0_xxzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyzzz_0[i] = 2.0 * g_0_xzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyzzz_0[i] * pb_x + g_0_xxzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyzzzz_0[i] = 2.0 * g_0_xzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyzzzz_0[i] * pb_x + g_0_xxzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyzzzzz_0[i] = 2.0 * g_0_xzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzzzzz_0[i] * pb_x + g_0_xxzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yzzzzzz_0[i] = 2.0 * g_0_xzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzzzzz_0[i] * pb_x + g_0_xxzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_zzzzzzz_0[i] = 2.0 * g_0_xzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzzzzz_0[i] * pb_x + g_0_xxzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 216-252 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_xxy_0_xxxxxxx_0, g_0_xxy_0_xxxxxxx_1, g_0_xxy_0_xxxxxxz_0, g_0_xxy_0_xxxxxxz_1, g_0_xxy_0_xxxxxzz_0, g_0_xxy_0_xxxxxzz_1, g_0_xxy_0_xxxxzzz_0, g_0_xxy_0_xxxxzzz_1, g_0_xxy_0_xxxzzzz_0, g_0_xxy_0_xxxzzzz_1, g_0_xxy_0_xxzzzzz_0, g_0_xxy_0_xxzzzzz_1, g_0_xxy_0_xzzzzzz_0, g_0_xxy_0_xzzzzzz_1, g_0_xxyy_0_xxxxxxx_0, g_0_xxyy_0_xxxxxxx_1, g_0_xxyy_0_xxxxxxz_0, g_0_xxyy_0_xxxxxxz_1, g_0_xxyy_0_xxxxxzz_0, g_0_xxyy_0_xxxxxzz_1, g_0_xxyy_0_xxxxzzz_0, g_0_xxyy_0_xxxxzzz_1, g_0_xxyy_0_xxxzzzz_0, g_0_xxyy_0_xxxzzzz_1, g_0_xxyy_0_xxzzzzz_0, g_0_xxyy_0_xxzzzzz_1, g_0_xxyy_0_xzzzzzz_0, g_0_xxyy_0_xzzzzzz_1, g_0_xxyyy_0_xxxxxxx_0, g_0_xxyyy_0_xxxxxxy_0, g_0_xxyyy_0_xxxxxxz_0, g_0_xxyyy_0_xxxxxyy_0, g_0_xxyyy_0_xxxxxyz_0, g_0_xxyyy_0_xxxxxzz_0, g_0_xxyyy_0_xxxxyyy_0, g_0_xxyyy_0_xxxxyyz_0, g_0_xxyyy_0_xxxxyzz_0, g_0_xxyyy_0_xxxxzzz_0, g_0_xxyyy_0_xxxyyyy_0, g_0_xxyyy_0_xxxyyyz_0, g_0_xxyyy_0_xxxyyzz_0, g_0_xxyyy_0_xxxyzzz_0, g_0_xxyyy_0_xxxzzzz_0, g_0_xxyyy_0_xxyyyyy_0, g_0_xxyyy_0_xxyyyyz_0, g_0_xxyyy_0_xxyyyzz_0, g_0_xxyyy_0_xxyyzzz_0, g_0_xxyyy_0_xxyzzzz_0, g_0_xxyyy_0_xxzzzzz_0, g_0_xxyyy_0_xyyyyyy_0, g_0_xxyyy_0_xyyyyyz_0, g_0_xxyyy_0_xyyyyzz_0, g_0_xxyyy_0_xyyyzzz_0, g_0_xxyyy_0_xyyzzzz_0, g_0_xxyyy_0_xyzzzzz_0, g_0_xxyyy_0_xzzzzzz_0, g_0_xxyyy_0_yyyyyyy_0, g_0_xxyyy_0_yyyyyyz_0, g_0_xxyyy_0_yyyyyzz_0, g_0_xxyyy_0_yyyyzzz_0, g_0_xxyyy_0_yyyzzzz_0, g_0_xxyyy_0_yyzzzzz_0, g_0_xxyyy_0_yzzzzzz_0, g_0_xxyyy_0_zzzzzzz_0, g_0_xyyy_0_xxxxxxy_0, g_0_xyyy_0_xxxxxxy_1, g_0_xyyy_0_xxxxxy_1, g_0_xyyy_0_xxxxxyy_0, g_0_xyyy_0_xxxxxyy_1, g_0_xyyy_0_xxxxxyz_0, g_0_xyyy_0_xxxxxyz_1, g_0_xyyy_0_xxxxyy_1, g_0_xyyy_0_xxxxyyy_0, g_0_xyyy_0_xxxxyyy_1, g_0_xyyy_0_xxxxyyz_0, g_0_xyyy_0_xxxxyyz_1, g_0_xyyy_0_xxxxyz_1, g_0_xyyy_0_xxxxyzz_0, g_0_xyyy_0_xxxxyzz_1, g_0_xyyy_0_xxxyyy_1, g_0_xyyy_0_xxxyyyy_0, g_0_xyyy_0_xxxyyyy_1, g_0_xyyy_0_xxxyyyz_0, g_0_xyyy_0_xxxyyyz_1, g_0_xyyy_0_xxxyyz_1, g_0_xyyy_0_xxxyyzz_0, g_0_xyyy_0_xxxyyzz_1, g_0_xyyy_0_xxxyzz_1, g_0_xyyy_0_xxxyzzz_0, g_0_xyyy_0_xxxyzzz_1, g_0_xyyy_0_xxyyyy_1, g_0_xyyy_0_xxyyyyy_0, g_0_xyyy_0_xxyyyyy_1, g_0_xyyy_0_xxyyyyz_0, g_0_xyyy_0_xxyyyyz_1, g_0_xyyy_0_xxyyyz_1, g_0_xyyy_0_xxyyyzz_0, g_0_xyyy_0_xxyyyzz_1, g_0_xyyy_0_xxyyzz_1, g_0_xyyy_0_xxyyzzz_0, g_0_xyyy_0_xxyyzzz_1, g_0_xyyy_0_xxyzzz_1, g_0_xyyy_0_xxyzzzz_0, g_0_xyyy_0_xxyzzzz_1, g_0_xyyy_0_xyyyyy_1, g_0_xyyy_0_xyyyyyy_0, g_0_xyyy_0_xyyyyyy_1, g_0_xyyy_0_xyyyyyz_0, g_0_xyyy_0_xyyyyyz_1, g_0_xyyy_0_xyyyyz_1, g_0_xyyy_0_xyyyyzz_0, g_0_xyyy_0_xyyyyzz_1, g_0_xyyy_0_xyyyzz_1, g_0_xyyy_0_xyyyzzz_0, g_0_xyyy_0_xyyyzzz_1, g_0_xyyy_0_xyyzzz_1, g_0_xyyy_0_xyyzzzz_0, g_0_xyyy_0_xyyzzzz_1, g_0_xyyy_0_xyzzzz_1, g_0_xyyy_0_xyzzzzz_0, g_0_xyyy_0_xyzzzzz_1, g_0_xyyy_0_yyyyyy_1, g_0_xyyy_0_yyyyyyy_0, g_0_xyyy_0_yyyyyyy_1, g_0_xyyy_0_yyyyyyz_0, g_0_xyyy_0_yyyyyyz_1, g_0_xyyy_0_yyyyyz_1, g_0_xyyy_0_yyyyyzz_0, g_0_xyyy_0_yyyyyzz_1, g_0_xyyy_0_yyyyzz_1, g_0_xyyy_0_yyyyzzz_0, g_0_xyyy_0_yyyyzzz_1, g_0_xyyy_0_yyyzzz_1, g_0_xyyy_0_yyyzzzz_0, g_0_xyyy_0_yyyzzzz_1, g_0_xyyy_0_yyzzzz_1, g_0_xyyy_0_yyzzzzz_0, g_0_xyyy_0_yyzzzzz_1, g_0_xyyy_0_yzzzzz_1, g_0_xyyy_0_yzzzzzz_0, g_0_xyyy_0_yzzzzzz_1, g_0_xyyy_0_zzzzzzz_0, g_0_xyyy_0_zzzzzzz_1, g_0_yyy_0_xxxxxxy_0, g_0_yyy_0_xxxxxxy_1, g_0_yyy_0_xxxxxyy_0, g_0_yyy_0_xxxxxyy_1, g_0_yyy_0_xxxxxyz_0, g_0_yyy_0_xxxxxyz_1, g_0_yyy_0_xxxxyyy_0, g_0_yyy_0_xxxxyyy_1, g_0_yyy_0_xxxxyyz_0, g_0_yyy_0_xxxxyyz_1, g_0_yyy_0_xxxxyzz_0, g_0_yyy_0_xxxxyzz_1, g_0_yyy_0_xxxyyyy_0, g_0_yyy_0_xxxyyyy_1, g_0_yyy_0_xxxyyyz_0, g_0_yyy_0_xxxyyyz_1, g_0_yyy_0_xxxyyzz_0, g_0_yyy_0_xxxyyzz_1, g_0_yyy_0_xxxyzzz_0, g_0_yyy_0_xxxyzzz_1, g_0_yyy_0_xxyyyyy_0, g_0_yyy_0_xxyyyyy_1, g_0_yyy_0_xxyyyyz_0, g_0_yyy_0_xxyyyyz_1, g_0_yyy_0_xxyyyzz_0, g_0_yyy_0_xxyyyzz_1, g_0_yyy_0_xxyyzzz_0, g_0_yyy_0_xxyyzzz_1, g_0_yyy_0_xxyzzzz_0, g_0_yyy_0_xxyzzzz_1, g_0_yyy_0_xyyyyyy_0, g_0_yyy_0_xyyyyyy_1, g_0_yyy_0_xyyyyyz_0, g_0_yyy_0_xyyyyyz_1, g_0_yyy_0_xyyyyzz_0, g_0_yyy_0_xyyyyzz_1, g_0_yyy_0_xyyyzzz_0, g_0_yyy_0_xyyyzzz_1, g_0_yyy_0_xyyzzzz_0, g_0_yyy_0_xyyzzzz_1, g_0_yyy_0_xyzzzzz_0, g_0_yyy_0_xyzzzzz_1, g_0_yyy_0_yyyyyyy_0, g_0_yyy_0_yyyyyyy_1, g_0_yyy_0_yyyyyyz_0, g_0_yyy_0_yyyyyyz_1, g_0_yyy_0_yyyyyzz_0, g_0_yyy_0_yyyyyzz_1, g_0_yyy_0_yyyyzzz_0, g_0_yyy_0_yyyyzzz_1, g_0_yyy_0_yyyzzzz_0, g_0_yyy_0_yyyzzzz_1, g_0_yyy_0_yyzzzzz_0, g_0_yyy_0_yyzzzzz_1, g_0_yyy_0_yzzzzzz_0, g_0_yyy_0_yzzzzzz_1, g_0_yyy_0_zzzzzzz_0, g_0_yyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyy_0_xxxxxxx_0[i] = 2.0 * g_0_xxy_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxxx_0[i] * pb_y + g_0_xxyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxxxy_0[i] = g_0_yyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxxy_0[i] * pb_x + g_0_xyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxxz_0[i] = 2.0 * g_0_xxy_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxxz_0[i] * pb_y + g_0_xxyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxxyy_0[i] = g_0_yyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxyy_0[i] * pb_x + g_0_xyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxyz_0[i] = g_0_yyy_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxyz_0[i] * pb_x + g_0_xyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxzz_0[i] = 2.0 * g_0_xxy_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxzz_0[i] * pb_y + g_0_xxyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxyyy_0[i] = g_0_yyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxyyy_0[i] * pb_x + g_0_xyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxyyz_0[i] = g_0_yyy_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxyyz_0[i] * pb_x + g_0_xyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxyzz_0[i] = g_0_yyy_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxyzz_0[i] * pb_x + g_0_xyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxzzz_0[i] = 2.0 * g_0_xxy_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxzzz_0[i] * pb_y + g_0_xxyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxyyyy_0[i] = g_0_yyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyyyy_0[i] * pb_x + g_0_xyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyyyz_0[i] = g_0_yyy_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyyyz_0[i] * pb_x + g_0_xyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyyzz_0[i] = g_0_yyy_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyyzz_0[i] * pb_x + g_0_xyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyzzz_0[i] = g_0_yyy_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyzzz_0[i] * pb_x + g_0_xyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxzzzz_0[i] = 2.0 * g_0_xxy_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxzzzz_0[i] * pb_y + g_0_xxyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxyyyyy_0[i] = g_0_yyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyyyy_0[i] * pb_x + g_0_xyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyyyz_0[i] = g_0_yyy_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyyyz_0[i] * pb_x + g_0_xyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyyzz_0[i] = g_0_yyy_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyyzz_0[i] * pb_x + g_0_xyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyzzz_0[i] = g_0_yyy_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyzzz_0[i] * pb_x + g_0_xyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyzzzz_0[i] = g_0_yyy_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyzzzz_0[i] * pb_x + g_0_xyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxzzzzz_0[i] = 2.0 * g_0_xxy_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxzzzzz_0[i] * pb_y + g_0_xxyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xyyyyyy_0[i] = g_0_yyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyyyy_0[i] * pb_x + g_0_xyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyyyz_0[i] = g_0_yyy_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyyyz_0[i] * pb_x + g_0_xyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyyzz_0[i] = g_0_yyy_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyyzz_0[i] * pb_x + g_0_xyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyzzz_0[i] = g_0_yyy_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyzzz_0[i] * pb_x + g_0_xyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyzzzz_0[i] = g_0_yyy_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyzzzz_0[i] * pb_x + g_0_xyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyzzzzz_0[i] = g_0_yyy_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyzzzzz_0[i] * pb_x + g_0_xyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xzzzzzz_0[i] = 2.0 * g_0_xxy_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xzzzzzz_0[i] * pb_y + g_0_xxyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_yyyyyyy_0[i] = g_0_yyy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyyy_0[i] * pb_x + g_0_xyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyyyz_0[i] = g_0_yyy_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyyz_0[i] * pb_x + g_0_xyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyyzz_0[i] = g_0_yyy_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyzz_0[i] * pb_x + g_0_xyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyzzz_0[i] = g_0_yyy_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyzzz_0[i] * pb_x + g_0_xyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyzzzz_0[i] = g_0_yyy_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyzzzz_0[i] * pb_x + g_0_xyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyzzzzz_0[i] = g_0_yyy_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzzzzz_0[i] * pb_x + g_0_xyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yzzzzzz_0[i] = g_0_yyy_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzzzzz_0[i] * pb_x + g_0_xyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_zzzzzzz_0[i] = g_0_yyy_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_zzzzzzz_0[i] * pb_x + g_0_xyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 252-288 components of targeted buffer : SHSK

    auto g_0_xxyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 252);

    auto g_0_xxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 253);

    auto g_0_xxyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 254);

    auto g_0_xxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 255);

    auto g_0_xxyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 256);

    auto g_0_xxyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 257);

    auto g_0_xxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 258);

    auto g_0_xxyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 259);

    auto g_0_xxyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 260);

    auto g_0_xxyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 261);

    auto g_0_xxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 262);

    auto g_0_xxyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 263);

    auto g_0_xxyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 264);

    auto g_0_xxyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 265);

    auto g_0_xxyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 266);

    auto g_0_xxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 267);

    auto g_0_xxyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 268);

    auto g_0_xxyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 269);

    auto g_0_xxyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 270);

    auto g_0_xxyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 271);

    auto g_0_xxyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 272);

    auto g_0_xxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 273);

    auto g_0_xxyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 274);

    auto g_0_xxyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 275);

    auto g_0_xxyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 276);

    auto g_0_xxyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 277);

    auto g_0_xxyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 278);

    auto g_0_xxyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 279);

    auto g_0_xxyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 280);

    auto g_0_xxyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 281);

    auto g_0_xxyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 282);

    auto g_0_xxyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 283);

    auto g_0_xxyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 284);

    auto g_0_xxyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 285);

    auto g_0_xxyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 286);

    auto g_0_xxyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 287);

    #pragma omp simd aligned(g_0_xxyy_0_xxxxxx_1, g_0_xxyy_0_xxxxxxx_0, g_0_xxyy_0_xxxxxxx_1, g_0_xxyy_0_xxxxxxy_0, g_0_xxyy_0_xxxxxxy_1, g_0_xxyy_0_xxxxxxz_0, g_0_xxyy_0_xxxxxxz_1, g_0_xxyy_0_xxxxxy_1, g_0_xxyy_0_xxxxxyy_0, g_0_xxyy_0_xxxxxyy_1, g_0_xxyy_0_xxxxxyz_0, g_0_xxyy_0_xxxxxyz_1, g_0_xxyy_0_xxxxxz_1, g_0_xxyy_0_xxxxxzz_0, g_0_xxyy_0_xxxxxzz_1, g_0_xxyy_0_xxxxyy_1, g_0_xxyy_0_xxxxyyy_0, g_0_xxyy_0_xxxxyyy_1, g_0_xxyy_0_xxxxyyz_0, g_0_xxyy_0_xxxxyyz_1, g_0_xxyy_0_xxxxyz_1, g_0_xxyy_0_xxxxyzz_0, g_0_xxyy_0_xxxxyzz_1, g_0_xxyy_0_xxxxzz_1, g_0_xxyy_0_xxxxzzz_0, g_0_xxyy_0_xxxxzzz_1, g_0_xxyy_0_xxxyyy_1, g_0_xxyy_0_xxxyyyy_0, g_0_xxyy_0_xxxyyyy_1, g_0_xxyy_0_xxxyyyz_0, g_0_xxyy_0_xxxyyyz_1, g_0_xxyy_0_xxxyyz_1, g_0_xxyy_0_xxxyyzz_0, g_0_xxyy_0_xxxyyzz_1, g_0_xxyy_0_xxxyzz_1, g_0_xxyy_0_xxxyzzz_0, g_0_xxyy_0_xxxyzzz_1, g_0_xxyy_0_xxxzzz_1, g_0_xxyy_0_xxxzzzz_0, g_0_xxyy_0_xxxzzzz_1, g_0_xxyy_0_xxyyyy_1, g_0_xxyy_0_xxyyyyy_0, g_0_xxyy_0_xxyyyyy_1, g_0_xxyy_0_xxyyyyz_0, g_0_xxyy_0_xxyyyyz_1, g_0_xxyy_0_xxyyyz_1, g_0_xxyy_0_xxyyyzz_0, g_0_xxyy_0_xxyyyzz_1, g_0_xxyy_0_xxyyzz_1, g_0_xxyy_0_xxyyzzz_0, g_0_xxyy_0_xxyyzzz_1, g_0_xxyy_0_xxyzzz_1, g_0_xxyy_0_xxyzzzz_0, g_0_xxyy_0_xxyzzzz_1, g_0_xxyy_0_xxzzzz_1, g_0_xxyy_0_xxzzzzz_0, g_0_xxyy_0_xxzzzzz_1, g_0_xxyy_0_xyyyyy_1, g_0_xxyy_0_xyyyyyy_0, g_0_xxyy_0_xyyyyyy_1, g_0_xxyy_0_xyyyyyz_0, g_0_xxyy_0_xyyyyyz_1, g_0_xxyy_0_xyyyyz_1, g_0_xxyy_0_xyyyyzz_0, g_0_xxyy_0_xyyyyzz_1, g_0_xxyy_0_xyyyzz_1, g_0_xxyy_0_xyyyzzz_0, g_0_xxyy_0_xyyyzzz_1, g_0_xxyy_0_xyyzzz_1, g_0_xxyy_0_xyyzzzz_0, g_0_xxyy_0_xyyzzzz_1, g_0_xxyy_0_xyzzzz_1, g_0_xxyy_0_xyzzzzz_0, g_0_xxyy_0_xyzzzzz_1, g_0_xxyy_0_xzzzzz_1, g_0_xxyy_0_xzzzzzz_0, g_0_xxyy_0_xzzzzzz_1, g_0_xxyy_0_yyyyyy_1, g_0_xxyy_0_yyyyyyy_0, g_0_xxyy_0_yyyyyyy_1, g_0_xxyy_0_yyyyyyz_0, g_0_xxyy_0_yyyyyyz_1, g_0_xxyy_0_yyyyyz_1, g_0_xxyy_0_yyyyyzz_0, g_0_xxyy_0_yyyyyzz_1, g_0_xxyy_0_yyyyzz_1, g_0_xxyy_0_yyyyzzz_0, g_0_xxyy_0_yyyyzzz_1, g_0_xxyy_0_yyyzzz_1, g_0_xxyy_0_yyyzzzz_0, g_0_xxyy_0_yyyzzzz_1, g_0_xxyy_0_yyzzzz_1, g_0_xxyy_0_yyzzzzz_0, g_0_xxyy_0_yyzzzzz_1, g_0_xxyy_0_yzzzzz_1, g_0_xxyy_0_yzzzzzz_0, g_0_xxyy_0_yzzzzzz_1, g_0_xxyy_0_zzzzzz_1, g_0_xxyy_0_zzzzzzz_0, g_0_xxyy_0_zzzzzzz_1, g_0_xxyyz_0_xxxxxxx_0, g_0_xxyyz_0_xxxxxxy_0, g_0_xxyyz_0_xxxxxxz_0, g_0_xxyyz_0_xxxxxyy_0, g_0_xxyyz_0_xxxxxyz_0, g_0_xxyyz_0_xxxxxzz_0, g_0_xxyyz_0_xxxxyyy_0, g_0_xxyyz_0_xxxxyyz_0, g_0_xxyyz_0_xxxxyzz_0, g_0_xxyyz_0_xxxxzzz_0, g_0_xxyyz_0_xxxyyyy_0, g_0_xxyyz_0_xxxyyyz_0, g_0_xxyyz_0_xxxyyzz_0, g_0_xxyyz_0_xxxyzzz_0, g_0_xxyyz_0_xxxzzzz_0, g_0_xxyyz_0_xxyyyyy_0, g_0_xxyyz_0_xxyyyyz_0, g_0_xxyyz_0_xxyyyzz_0, g_0_xxyyz_0_xxyyzzz_0, g_0_xxyyz_0_xxyzzzz_0, g_0_xxyyz_0_xxzzzzz_0, g_0_xxyyz_0_xyyyyyy_0, g_0_xxyyz_0_xyyyyyz_0, g_0_xxyyz_0_xyyyyzz_0, g_0_xxyyz_0_xyyyzzz_0, g_0_xxyyz_0_xyyzzzz_0, g_0_xxyyz_0_xyzzzzz_0, g_0_xxyyz_0_xzzzzzz_0, g_0_xxyyz_0_yyyyyyy_0, g_0_xxyyz_0_yyyyyyz_0, g_0_xxyyz_0_yyyyyzz_0, g_0_xxyyz_0_yyyyzzz_0, g_0_xxyyz_0_yyyzzzz_0, g_0_xxyyz_0_yyzzzzz_0, g_0_xxyyz_0_yzzzzzz_0, g_0_xxyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyz_0_xxxxxxx_0[i] = g_0_xxyy_0_xxxxxxx_0[i] * pb_z + g_0_xxyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxxy_0[i] = g_0_xxyy_0_xxxxxxy_0[i] * pb_z + g_0_xxyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxxz_0[i] = g_0_xxyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxxz_0[i] * pb_z + g_0_xxyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxyy_0[i] = g_0_xxyy_0_xxxxxyy_0[i] * pb_z + g_0_xxyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxyz_0[i] = g_0_xxyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxyz_0[i] * pb_z + g_0_xxyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxzz_0[i] * pb_z + g_0_xxyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyyy_0[i] = g_0_xxyy_0_xxxxyyy_0[i] * pb_z + g_0_xxyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyyz_0[i] = g_0_xxyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyyz_0[i] * pb_z + g_0_xxyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyzz_0[i] * pb_z + g_0_xxyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxzzz_0[i] * pb_z + g_0_xxyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyyy_0[i] = g_0_xxyy_0_xxxyyyy_0[i] * pb_z + g_0_xxyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyyz_0[i] = g_0_xxyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyyz_0[i] * pb_z + g_0_xxyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyzz_0[i] * pb_z + g_0_xxyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyzzz_0[i] * pb_z + g_0_xxyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxzzzz_0[i] * pb_z + g_0_xxyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyyy_0[i] = g_0_xxyy_0_xxyyyyy_0[i] * pb_z + g_0_xxyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyyz_0[i] = g_0_xxyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyyz_0[i] * pb_z + g_0_xxyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyzz_0[i] * pb_z + g_0_xxyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyzzz_0[i] * pb_z + g_0_xxyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyzzzz_0[i] * pb_z + g_0_xxyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxzzzzz_0[i] * pb_z + g_0_xxyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyyy_0[i] = g_0_xxyy_0_xyyyyyy_0[i] * pb_z + g_0_xxyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyyz_0[i] = g_0_xxyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyyz_0[i] * pb_z + g_0_xxyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyzz_0[i] * pb_z + g_0_xxyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyzzz_0[i] * pb_z + g_0_xxyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyzzzz_0[i] * pb_z + g_0_xxyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyzzzzz_0[i] * pb_z + g_0_xxyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xzzzzzz_0[i] * pb_z + g_0_xxyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyyy_0[i] = g_0_xxyy_0_yyyyyyy_0[i] * pb_z + g_0_xxyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyyz_0[i] = g_0_xxyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyyyz_0[i] * pb_z + g_0_xxyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyyzz_0[i] * pb_z + g_0_xxyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyzzz_0[i] * pb_z + g_0_xxyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyzzzz_0[i] * pb_z + g_0_xxyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyzzzzz_0[i] * pb_z + g_0_xxyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yzzzzzz_0[i] * pb_z + g_0_xxyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_zzzzzzz_0[i] * pb_z + g_0_xxyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 288-324 components of targeted buffer : SHSK

    auto g_0_xxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 288);

    auto g_0_xxyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 289);

    auto g_0_xxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 290);

    auto g_0_xxyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 291);

    auto g_0_xxyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 292);

    auto g_0_xxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 293);

    auto g_0_xxyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 294);

    auto g_0_xxyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 295);

    auto g_0_xxyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 296);

    auto g_0_xxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 297);

    auto g_0_xxyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 298);

    auto g_0_xxyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 299);

    auto g_0_xxyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 300);

    auto g_0_xxyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 301);

    auto g_0_xxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 302);

    auto g_0_xxyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 303);

    auto g_0_xxyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 304);

    auto g_0_xxyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 305);

    auto g_0_xxyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 306);

    auto g_0_xxyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 307);

    auto g_0_xxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 308);

    auto g_0_xxyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 309);

    auto g_0_xxyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 310);

    auto g_0_xxyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 311);

    auto g_0_xxyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 312);

    auto g_0_xxyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 313);

    auto g_0_xxyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 314);

    auto g_0_xxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 315);

    auto g_0_xxyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 316);

    auto g_0_xxyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 317);

    auto g_0_xxyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 318);

    auto g_0_xxyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 319);

    auto g_0_xxyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 320);

    auto g_0_xxyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 321);

    auto g_0_xxyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 322);

    auto g_0_xxyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 323);

    #pragma omp simd aligned(g_0_xxyzz_0_xxxxxxx_0, g_0_xxyzz_0_xxxxxxy_0, g_0_xxyzz_0_xxxxxxz_0, g_0_xxyzz_0_xxxxxyy_0, g_0_xxyzz_0_xxxxxyz_0, g_0_xxyzz_0_xxxxxzz_0, g_0_xxyzz_0_xxxxyyy_0, g_0_xxyzz_0_xxxxyyz_0, g_0_xxyzz_0_xxxxyzz_0, g_0_xxyzz_0_xxxxzzz_0, g_0_xxyzz_0_xxxyyyy_0, g_0_xxyzz_0_xxxyyyz_0, g_0_xxyzz_0_xxxyyzz_0, g_0_xxyzz_0_xxxyzzz_0, g_0_xxyzz_0_xxxzzzz_0, g_0_xxyzz_0_xxyyyyy_0, g_0_xxyzz_0_xxyyyyz_0, g_0_xxyzz_0_xxyyyzz_0, g_0_xxyzz_0_xxyyzzz_0, g_0_xxyzz_0_xxyzzzz_0, g_0_xxyzz_0_xxzzzzz_0, g_0_xxyzz_0_xyyyyyy_0, g_0_xxyzz_0_xyyyyyz_0, g_0_xxyzz_0_xyyyyzz_0, g_0_xxyzz_0_xyyyzzz_0, g_0_xxyzz_0_xyyzzzz_0, g_0_xxyzz_0_xyzzzzz_0, g_0_xxyzz_0_xzzzzzz_0, g_0_xxyzz_0_yyyyyyy_0, g_0_xxyzz_0_yyyyyyz_0, g_0_xxyzz_0_yyyyyzz_0, g_0_xxyzz_0_yyyyzzz_0, g_0_xxyzz_0_yyyzzzz_0, g_0_xxyzz_0_yyzzzzz_0, g_0_xxyzz_0_yzzzzzz_0, g_0_xxyzz_0_zzzzzzz_0, g_0_xxzz_0_xxxxxx_1, g_0_xxzz_0_xxxxxxx_0, g_0_xxzz_0_xxxxxxx_1, g_0_xxzz_0_xxxxxxy_0, g_0_xxzz_0_xxxxxxy_1, g_0_xxzz_0_xxxxxxz_0, g_0_xxzz_0_xxxxxxz_1, g_0_xxzz_0_xxxxxy_1, g_0_xxzz_0_xxxxxyy_0, g_0_xxzz_0_xxxxxyy_1, g_0_xxzz_0_xxxxxyz_0, g_0_xxzz_0_xxxxxyz_1, g_0_xxzz_0_xxxxxz_1, g_0_xxzz_0_xxxxxzz_0, g_0_xxzz_0_xxxxxzz_1, g_0_xxzz_0_xxxxyy_1, g_0_xxzz_0_xxxxyyy_0, g_0_xxzz_0_xxxxyyy_1, g_0_xxzz_0_xxxxyyz_0, g_0_xxzz_0_xxxxyyz_1, g_0_xxzz_0_xxxxyz_1, g_0_xxzz_0_xxxxyzz_0, g_0_xxzz_0_xxxxyzz_1, g_0_xxzz_0_xxxxzz_1, g_0_xxzz_0_xxxxzzz_0, g_0_xxzz_0_xxxxzzz_1, g_0_xxzz_0_xxxyyy_1, g_0_xxzz_0_xxxyyyy_0, g_0_xxzz_0_xxxyyyy_1, g_0_xxzz_0_xxxyyyz_0, g_0_xxzz_0_xxxyyyz_1, g_0_xxzz_0_xxxyyz_1, g_0_xxzz_0_xxxyyzz_0, g_0_xxzz_0_xxxyyzz_1, g_0_xxzz_0_xxxyzz_1, g_0_xxzz_0_xxxyzzz_0, g_0_xxzz_0_xxxyzzz_1, g_0_xxzz_0_xxxzzz_1, g_0_xxzz_0_xxxzzzz_0, g_0_xxzz_0_xxxzzzz_1, g_0_xxzz_0_xxyyyy_1, g_0_xxzz_0_xxyyyyy_0, g_0_xxzz_0_xxyyyyy_1, g_0_xxzz_0_xxyyyyz_0, g_0_xxzz_0_xxyyyyz_1, g_0_xxzz_0_xxyyyz_1, g_0_xxzz_0_xxyyyzz_0, g_0_xxzz_0_xxyyyzz_1, g_0_xxzz_0_xxyyzz_1, g_0_xxzz_0_xxyyzzz_0, g_0_xxzz_0_xxyyzzz_1, g_0_xxzz_0_xxyzzz_1, g_0_xxzz_0_xxyzzzz_0, g_0_xxzz_0_xxyzzzz_1, g_0_xxzz_0_xxzzzz_1, g_0_xxzz_0_xxzzzzz_0, g_0_xxzz_0_xxzzzzz_1, g_0_xxzz_0_xyyyyy_1, g_0_xxzz_0_xyyyyyy_0, g_0_xxzz_0_xyyyyyy_1, g_0_xxzz_0_xyyyyyz_0, g_0_xxzz_0_xyyyyyz_1, g_0_xxzz_0_xyyyyz_1, g_0_xxzz_0_xyyyyzz_0, g_0_xxzz_0_xyyyyzz_1, g_0_xxzz_0_xyyyzz_1, g_0_xxzz_0_xyyyzzz_0, g_0_xxzz_0_xyyyzzz_1, g_0_xxzz_0_xyyzzz_1, g_0_xxzz_0_xyyzzzz_0, g_0_xxzz_0_xyyzzzz_1, g_0_xxzz_0_xyzzzz_1, g_0_xxzz_0_xyzzzzz_0, g_0_xxzz_0_xyzzzzz_1, g_0_xxzz_0_xzzzzz_1, g_0_xxzz_0_xzzzzzz_0, g_0_xxzz_0_xzzzzzz_1, g_0_xxzz_0_yyyyyy_1, g_0_xxzz_0_yyyyyyy_0, g_0_xxzz_0_yyyyyyy_1, g_0_xxzz_0_yyyyyyz_0, g_0_xxzz_0_yyyyyyz_1, g_0_xxzz_0_yyyyyz_1, g_0_xxzz_0_yyyyyzz_0, g_0_xxzz_0_yyyyyzz_1, g_0_xxzz_0_yyyyzz_1, g_0_xxzz_0_yyyyzzz_0, g_0_xxzz_0_yyyyzzz_1, g_0_xxzz_0_yyyzzz_1, g_0_xxzz_0_yyyzzzz_0, g_0_xxzz_0_yyyzzzz_1, g_0_xxzz_0_yyzzzz_1, g_0_xxzz_0_yyzzzzz_0, g_0_xxzz_0_yyzzzzz_1, g_0_xxzz_0_yzzzzz_1, g_0_xxzz_0_yzzzzzz_0, g_0_xxzz_0_yzzzzzz_1, g_0_xxzz_0_zzzzzz_1, g_0_xxzz_0_zzzzzzz_0, g_0_xxzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzz_0_xxxxxxx_0[i] = g_0_xxzz_0_xxxxxxx_0[i] * pb_y + g_0_xxzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxxy_0[i] = g_0_xxzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxxy_0[i] * pb_y + g_0_xxzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxxz_0[i] = g_0_xxzz_0_xxxxxxz_0[i] * pb_y + g_0_xxzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxyy_0[i] * pb_y + g_0_xxzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxyz_0[i] = g_0_xxzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxyz_0[i] * pb_y + g_0_xxzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxzz_0[i] = g_0_xxzz_0_xxxxxzz_0[i] * pb_y + g_0_xxzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyyy_0[i] * pb_y + g_0_xxzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyyz_0[i] * pb_y + g_0_xxzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyzz_0[i] = g_0_xxzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyzz_0[i] * pb_y + g_0_xxzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxzzz_0[i] = g_0_xxzz_0_xxxxzzz_0[i] * pb_y + g_0_xxzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyyy_0[i] * pb_y + g_0_xxzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyyz_0[i] * pb_y + g_0_xxzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyzz_0[i] * pb_y + g_0_xxzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyzzz_0[i] = g_0_xxzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyzzz_0[i] * pb_y + g_0_xxzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxzzzz_0[i] = g_0_xxzz_0_xxxzzzz_0[i] * pb_y + g_0_xxzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyyy_0[i] * pb_y + g_0_xxzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyyz_0[i] * pb_y + g_0_xxzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyzz_0[i] * pb_y + g_0_xxzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyzzz_0[i] * pb_y + g_0_xxzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyzzzz_0[i] = g_0_xxzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyzzzz_0[i] * pb_y + g_0_xxzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxzzzzz_0[i] = g_0_xxzz_0_xxzzzzz_0[i] * pb_y + g_0_xxzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyyy_0[i] * pb_y + g_0_xxzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyyz_0[i] * pb_y + g_0_xxzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyzz_0[i] * pb_y + g_0_xxzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyzzz_0[i] * pb_y + g_0_xxzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyzzzz_0[i] * pb_y + g_0_xxzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyzzzzz_0[i] = g_0_xxzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyzzzzz_0[i] * pb_y + g_0_xxzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xzzzzzz_0[i] = g_0_xxzz_0_xzzzzzz_0[i] * pb_y + g_0_xxzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyyy_0[i] * pb_y + g_0_xxzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyyz_0[i] * pb_y + g_0_xxzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyzz_0[i] * pb_y + g_0_xxzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyzzz_0[i] * pb_y + g_0_xxzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyzzzz_0[i] * pb_y + g_0_xxzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyzzzzz_0[i] * pb_y + g_0_xxzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yzzzzzz_0[i] = g_0_xxzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yzzzzzz_0[i] * pb_y + g_0_xxzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_zzzzzzz_0[i] = g_0_xxzz_0_zzzzzzz_0[i] * pb_y + g_0_xxzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 324-360 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_xxz_0_xxxxxxx_0, g_0_xxz_0_xxxxxxx_1, g_0_xxz_0_xxxxxxy_0, g_0_xxz_0_xxxxxxy_1, g_0_xxz_0_xxxxxyy_0, g_0_xxz_0_xxxxxyy_1, g_0_xxz_0_xxxxyyy_0, g_0_xxz_0_xxxxyyy_1, g_0_xxz_0_xxxyyyy_0, g_0_xxz_0_xxxyyyy_1, g_0_xxz_0_xxyyyyy_0, g_0_xxz_0_xxyyyyy_1, g_0_xxz_0_xyyyyyy_0, g_0_xxz_0_xyyyyyy_1, g_0_xxzz_0_xxxxxxx_0, g_0_xxzz_0_xxxxxxx_1, g_0_xxzz_0_xxxxxxy_0, g_0_xxzz_0_xxxxxxy_1, g_0_xxzz_0_xxxxxyy_0, g_0_xxzz_0_xxxxxyy_1, g_0_xxzz_0_xxxxyyy_0, g_0_xxzz_0_xxxxyyy_1, g_0_xxzz_0_xxxyyyy_0, g_0_xxzz_0_xxxyyyy_1, g_0_xxzz_0_xxyyyyy_0, g_0_xxzz_0_xxyyyyy_1, g_0_xxzz_0_xyyyyyy_0, g_0_xxzz_0_xyyyyyy_1, g_0_xxzzz_0_xxxxxxx_0, g_0_xxzzz_0_xxxxxxy_0, g_0_xxzzz_0_xxxxxxz_0, g_0_xxzzz_0_xxxxxyy_0, g_0_xxzzz_0_xxxxxyz_0, g_0_xxzzz_0_xxxxxzz_0, g_0_xxzzz_0_xxxxyyy_0, g_0_xxzzz_0_xxxxyyz_0, g_0_xxzzz_0_xxxxyzz_0, g_0_xxzzz_0_xxxxzzz_0, g_0_xxzzz_0_xxxyyyy_0, g_0_xxzzz_0_xxxyyyz_0, g_0_xxzzz_0_xxxyyzz_0, g_0_xxzzz_0_xxxyzzz_0, g_0_xxzzz_0_xxxzzzz_0, g_0_xxzzz_0_xxyyyyy_0, g_0_xxzzz_0_xxyyyyz_0, g_0_xxzzz_0_xxyyyzz_0, g_0_xxzzz_0_xxyyzzz_0, g_0_xxzzz_0_xxyzzzz_0, g_0_xxzzz_0_xxzzzzz_0, g_0_xxzzz_0_xyyyyyy_0, g_0_xxzzz_0_xyyyyyz_0, g_0_xxzzz_0_xyyyyzz_0, g_0_xxzzz_0_xyyyzzz_0, g_0_xxzzz_0_xyyzzzz_0, g_0_xxzzz_0_xyzzzzz_0, g_0_xxzzz_0_xzzzzzz_0, g_0_xxzzz_0_yyyyyyy_0, g_0_xxzzz_0_yyyyyyz_0, g_0_xxzzz_0_yyyyyzz_0, g_0_xxzzz_0_yyyyzzz_0, g_0_xxzzz_0_yyyzzzz_0, g_0_xxzzz_0_yyzzzzz_0, g_0_xxzzz_0_yzzzzzz_0, g_0_xxzzz_0_zzzzzzz_0, g_0_xzzz_0_xxxxxxz_0, g_0_xzzz_0_xxxxxxz_1, g_0_xzzz_0_xxxxxyz_0, g_0_xzzz_0_xxxxxyz_1, g_0_xzzz_0_xxxxxz_1, g_0_xzzz_0_xxxxxzz_0, g_0_xzzz_0_xxxxxzz_1, g_0_xzzz_0_xxxxyyz_0, g_0_xzzz_0_xxxxyyz_1, g_0_xzzz_0_xxxxyz_1, g_0_xzzz_0_xxxxyzz_0, g_0_xzzz_0_xxxxyzz_1, g_0_xzzz_0_xxxxzz_1, g_0_xzzz_0_xxxxzzz_0, g_0_xzzz_0_xxxxzzz_1, g_0_xzzz_0_xxxyyyz_0, g_0_xzzz_0_xxxyyyz_1, g_0_xzzz_0_xxxyyz_1, g_0_xzzz_0_xxxyyzz_0, g_0_xzzz_0_xxxyyzz_1, g_0_xzzz_0_xxxyzz_1, g_0_xzzz_0_xxxyzzz_0, g_0_xzzz_0_xxxyzzz_1, g_0_xzzz_0_xxxzzz_1, g_0_xzzz_0_xxxzzzz_0, g_0_xzzz_0_xxxzzzz_1, g_0_xzzz_0_xxyyyyz_0, g_0_xzzz_0_xxyyyyz_1, g_0_xzzz_0_xxyyyz_1, g_0_xzzz_0_xxyyyzz_0, g_0_xzzz_0_xxyyyzz_1, g_0_xzzz_0_xxyyzz_1, g_0_xzzz_0_xxyyzzz_0, g_0_xzzz_0_xxyyzzz_1, g_0_xzzz_0_xxyzzz_1, g_0_xzzz_0_xxyzzzz_0, g_0_xzzz_0_xxyzzzz_1, g_0_xzzz_0_xxzzzz_1, g_0_xzzz_0_xxzzzzz_0, g_0_xzzz_0_xxzzzzz_1, g_0_xzzz_0_xyyyyyz_0, g_0_xzzz_0_xyyyyyz_1, g_0_xzzz_0_xyyyyz_1, g_0_xzzz_0_xyyyyzz_0, g_0_xzzz_0_xyyyyzz_1, g_0_xzzz_0_xyyyzz_1, g_0_xzzz_0_xyyyzzz_0, g_0_xzzz_0_xyyyzzz_1, g_0_xzzz_0_xyyzzz_1, g_0_xzzz_0_xyyzzzz_0, g_0_xzzz_0_xyyzzzz_1, g_0_xzzz_0_xyzzzz_1, g_0_xzzz_0_xyzzzzz_0, g_0_xzzz_0_xyzzzzz_1, g_0_xzzz_0_xzzzzz_1, g_0_xzzz_0_xzzzzzz_0, g_0_xzzz_0_xzzzzzz_1, g_0_xzzz_0_yyyyyyy_0, g_0_xzzz_0_yyyyyyy_1, g_0_xzzz_0_yyyyyyz_0, g_0_xzzz_0_yyyyyyz_1, g_0_xzzz_0_yyyyyz_1, g_0_xzzz_0_yyyyyzz_0, g_0_xzzz_0_yyyyyzz_1, g_0_xzzz_0_yyyyzz_1, g_0_xzzz_0_yyyyzzz_0, g_0_xzzz_0_yyyyzzz_1, g_0_xzzz_0_yyyzzz_1, g_0_xzzz_0_yyyzzzz_0, g_0_xzzz_0_yyyzzzz_1, g_0_xzzz_0_yyzzzz_1, g_0_xzzz_0_yyzzzzz_0, g_0_xzzz_0_yyzzzzz_1, g_0_xzzz_0_yzzzzz_1, g_0_xzzz_0_yzzzzzz_0, g_0_xzzz_0_yzzzzzz_1, g_0_xzzz_0_zzzzzz_1, g_0_xzzz_0_zzzzzzz_0, g_0_xzzz_0_zzzzzzz_1, g_0_zzz_0_xxxxxxz_0, g_0_zzz_0_xxxxxxz_1, g_0_zzz_0_xxxxxyz_0, g_0_zzz_0_xxxxxyz_1, g_0_zzz_0_xxxxxzz_0, g_0_zzz_0_xxxxxzz_1, g_0_zzz_0_xxxxyyz_0, g_0_zzz_0_xxxxyyz_1, g_0_zzz_0_xxxxyzz_0, g_0_zzz_0_xxxxyzz_1, g_0_zzz_0_xxxxzzz_0, g_0_zzz_0_xxxxzzz_1, g_0_zzz_0_xxxyyyz_0, g_0_zzz_0_xxxyyyz_1, g_0_zzz_0_xxxyyzz_0, g_0_zzz_0_xxxyyzz_1, g_0_zzz_0_xxxyzzz_0, g_0_zzz_0_xxxyzzz_1, g_0_zzz_0_xxxzzzz_0, g_0_zzz_0_xxxzzzz_1, g_0_zzz_0_xxyyyyz_0, g_0_zzz_0_xxyyyyz_1, g_0_zzz_0_xxyyyzz_0, g_0_zzz_0_xxyyyzz_1, g_0_zzz_0_xxyyzzz_0, g_0_zzz_0_xxyyzzz_1, g_0_zzz_0_xxyzzzz_0, g_0_zzz_0_xxyzzzz_1, g_0_zzz_0_xxzzzzz_0, g_0_zzz_0_xxzzzzz_1, g_0_zzz_0_xyyyyyz_0, g_0_zzz_0_xyyyyyz_1, g_0_zzz_0_xyyyyzz_0, g_0_zzz_0_xyyyyzz_1, g_0_zzz_0_xyyyzzz_0, g_0_zzz_0_xyyyzzz_1, g_0_zzz_0_xyyzzzz_0, g_0_zzz_0_xyyzzzz_1, g_0_zzz_0_xyzzzzz_0, g_0_zzz_0_xyzzzzz_1, g_0_zzz_0_xzzzzzz_0, g_0_zzz_0_xzzzzzz_1, g_0_zzz_0_yyyyyyy_0, g_0_zzz_0_yyyyyyy_1, g_0_zzz_0_yyyyyyz_0, g_0_zzz_0_yyyyyyz_1, g_0_zzz_0_yyyyyzz_0, g_0_zzz_0_yyyyyzz_1, g_0_zzz_0_yyyyzzz_0, g_0_zzz_0_yyyyzzz_1, g_0_zzz_0_yyyzzzz_0, g_0_zzz_0_yyyzzzz_1, g_0_zzz_0_yyzzzzz_0, g_0_zzz_0_yyzzzzz_1, g_0_zzz_0_yzzzzzz_0, g_0_zzz_0_yzzzzzz_1, g_0_zzz_0_zzzzzzz_0, g_0_zzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzz_0_xxxxxxx_0[i] = 2.0 * g_0_xxz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxxx_0[i] * pb_z + g_0_xxzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxxy_0[i] = 2.0 * g_0_xxz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxxy_0[i] * pb_z + g_0_xxzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxxz_0[i] = g_0_zzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxxz_0[i] * pb_x + g_0_xzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxyy_0[i] * pb_z + g_0_xxzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxyz_0[i] = g_0_zzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxyz_0[i] * pb_x + g_0_xzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxxzz_0[i] = g_0_zzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxzz_0[i] * pb_x + g_0_xzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxyyy_0[i] = 2.0 * g_0_xxz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxyyy_0[i] * pb_z + g_0_xxzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxyyz_0[i] = g_0_zzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxyyz_0[i] * pb_x + g_0_xzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxyzz_0[i] = g_0_zzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxyzz_0[i] * pb_x + g_0_xzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxzzz_0[i] = g_0_zzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxzzz_0[i] * pb_x + g_0_xzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyyyy_0[i] = 2.0 * g_0_xxz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxyyyy_0[i] * pb_z + g_0_xxzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxyyyz_0[i] = g_0_zzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxyyyz_0[i] * pb_x + g_0_xzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyyzz_0[i] = g_0_zzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxyyzz_0[i] * pb_x + g_0_xzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyzzz_0[i] = g_0_zzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxyzzz_0[i] * pb_x + g_0_xzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxzzzz_0[i] = g_0_zzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxzzzz_0[i] * pb_x + g_0_xzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyyyy_0[i] = 2.0 * g_0_xxz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxyyyyy_0[i] * pb_z + g_0_xxzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxyyyyz_0[i] = g_0_zzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyyyyz_0[i] * pb_x + g_0_xzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyyzz_0[i] = g_0_zzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyyyzz_0[i] * pb_x + g_0_xzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyzzz_0[i] = g_0_zzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyyzzz_0[i] * pb_x + g_0_xzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyzzzz_0[i] = g_0_zzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyzzzz_0[i] * pb_x + g_0_xzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxzzzzz_0[i] = g_0_zzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxzzzzz_0[i] * pb_x + g_0_xzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyyyy_0[i] = 2.0 * g_0_xxz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xyyyyyy_0[i] * pb_z + g_0_xxzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xyyyyyz_0[i] = g_0_zzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyyyyz_0[i] * pb_x + g_0_xzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyyzz_0[i] = g_0_zzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyyyzz_0[i] * pb_x + g_0_xzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyzzz_0[i] = g_0_zzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyyzzz_0[i] * pb_x + g_0_xzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyzzzz_0[i] = g_0_zzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyzzzz_0[i] * pb_x + g_0_xzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyzzzzz_0[i] = g_0_zzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyzzzzz_0[i] * pb_x + g_0_xzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xzzzzzz_0[i] = g_0_zzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xzzzzzz_0[i] * pb_x + g_0_xzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyyy_0[i] = g_0_zzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyyy_0[i] * pb_x + g_0_xzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyyz_0[i] = g_0_zzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyyz_0[i] * pb_x + g_0_xzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyzz_0[i] = g_0_zzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyzz_0[i] * pb_x + g_0_xzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyzzz_0[i] = g_0_zzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyzzz_0[i] * pb_x + g_0_xzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyzzzz_0[i] = g_0_zzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyzzzz_0[i] * pb_x + g_0_xzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyzzzzz_0[i] = g_0_zzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzzzzz_0[i] * pb_x + g_0_xzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yzzzzzz_0[i] = g_0_zzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzzzzz_0[i] * pb_x + g_0_xzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_zzzzzzz_0[i] = g_0_zzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzzzzz_0[i] * pb_x + g_0_xzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 360-396 components of targeted buffer : SHSK

    auto g_0_xyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 360);

    auto g_0_xyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 361);

    auto g_0_xyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 362);

    auto g_0_xyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 363);

    auto g_0_xyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 364);

    auto g_0_xyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 365);

    auto g_0_xyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 366);

    auto g_0_xyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 367);

    auto g_0_xyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 368);

    auto g_0_xyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 369);

    auto g_0_xyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 370);

    auto g_0_xyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 371);

    auto g_0_xyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 372);

    auto g_0_xyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 373);

    auto g_0_xyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 374);

    auto g_0_xyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 375);

    auto g_0_xyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 376);

    auto g_0_xyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 377);

    auto g_0_xyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 378);

    auto g_0_xyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 379);

    auto g_0_xyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 380);

    auto g_0_xyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 381);

    auto g_0_xyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 382);

    auto g_0_xyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 383);

    auto g_0_xyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 384);

    auto g_0_xyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 385);

    auto g_0_xyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 386);

    auto g_0_xyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 387);

    auto g_0_xyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 388);

    auto g_0_xyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 389);

    auto g_0_xyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 390);

    auto g_0_xyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 391);

    auto g_0_xyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 392);

    auto g_0_xyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 393);

    auto g_0_xyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 394);

    auto g_0_xyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 395);

    #pragma omp simd aligned(g_0_xyyyy_0_xxxxxxx_0, g_0_xyyyy_0_xxxxxxy_0, g_0_xyyyy_0_xxxxxxz_0, g_0_xyyyy_0_xxxxxyy_0, g_0_xyyyy_0_xxxxxyz_0, g_0_xyyyy_0_xxxxxzz_0, g_0_xyyyy_0_xxxxyyy_0, g_0_xyyyy_0_xxxxyyz_0, g_0_xyyyy_0_xxxxyzz_0, g_0_xyyyy_0_xxxxzzz_0, g_0_xyyyy_0_xxxyyyy_0, g_0_xyyyy_0_xxxyyyz_0, g_0_xyyyy_0_xxxyyzz_0, g_0_xyyyy_0_xxxyzzz_0, g_0_xyyyy_0_xxxzzzz_0, g_0_xyyyy_0_xxyyyyy_0, g_0_xyyyy_0_xxyyyyz_0, g_0_xyyyy_0_xxyyyzz_0, g_0_xyyyy_0_xxyyzzz_0, g_0_xyyyy_0_xxyzzzz_0, g_0_xyyyy_0_xxzzzzz_0, g_0_xyyyy_0_xyyyyyy_0, g_0_xyyyy_0_xyyyyyz_0, g_0_xyyyy_0_xyyyyzz_0, g_0_xyyyy_0_xyyyzzz_0, g_0_xyyyy_0_xyyzzzz_0, g_0_xyyyy_0_xyzzzzz_0, g_0_xyyyy_0_xzzzzzz_0, g_0_xyyyy_0_yyyyyyy_0, g_0_xyyyy_0_yyyyyyz_0, g_0_xyyyy_0_yyyyyzz_0, g_0_xyyyy_0_yyyyzzz_0, g_0_xyyyy_0_yyyzzzz_0, g_0_xyyyy_0_yyzzzzz_0, g_0_xyyyy_0_yzzzzzz_0, g_0_xyyyy_0_zzzzzzz_0, g_0_yyyy_0_xxxxxx_1, g_0_yyyy_0_xxxxxxx_0, g_0_yyyy_0_xxxxxxx_1, g_0_yyyy_0_xxxxxxy_0, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxxz_0, g_0_yyyy_0_xxxxxxz_1, g_0_yyyy_0_xxxxxy_1, g_0_yyyy_0_xxxxxyy_0, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxxyz_0, g_0_yyyy_0_xxxxxyz_1, g_0_yyyy_0_xxxxxz_1, g_0_yyyy_0_xxxxxzz_0, g_0_yyyy_0_xxxxxzz_1, g_0_yyyy_0_xxxxyy_1, g_0_yyyy_0_xxxxyyy_0, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxxyyz_0, g_0_yyyy_0_xxxxyyz_1, g_0_yyyy_0_xxxxyz_1, g_0_yyyy_0_xxxxyzz_0, g_0_yyyy_0_xxxxyzz_1, g_0_yyyy_0_xxxxzz_1, g_0_yyyy_0_xxxxzzz_0, g_0_yyyy_0_xxxxzzz_1, g_0_yyyy_0_xxxyyy_1, g_0_yyyy_0_xxxyyyy_0, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxxyyyz_0, g_0_yyyy_0_xxxyyyz_1, g_0_yyyy_0_xxxyyz_1, g_0_yyyy_0_xxxyyzz_0, g_0_yyyy_0_xxxyyzz_1, g_0_yyyy_0_xxxyzz_1, g_0_yyyy_0_xxxyzzz_0, g_0_yyyy_0_xxxyzzz_1, g_0_yyyy_0_xxxzzz_1, g_0_yyyy_0_xxxzzzz_0, g_0_yyyy_0_xxxzzzz_1, g_0_yyyy_0_xxyyyy_1, g_0_yyyy_0_xxyyyyy_0, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xxyyyyz_0, g_0_yyyy_0_xxyyyyz_1, g_0_yyyy_0_xxyyyz_1, g_0_yyyy_0_xxyyyzz_0, g_0_yyyy_0_xxyyyzz_1, g_0_yyyy_0_xxyyzz_1, g_0_yyyy_0_xxyyzzz_0, g_0_yyyy_0_xxyyzzz_1, g_0_yyyy_0_xxyzzz_1, g_0_yyyy_0_xxyzzzz_0, g_0_yyyy_0_xxyzzzz_1, g_0_yyyy_0_xxzzzz_1, g_0_yyyy_0_xxzzzzz_0, g_0_yyyy_0_xxzzzzz_1, g_0_yyyy_0_xyyyyy_1, g_0_yyyy_0_xyyyyyy_0, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_xyyyyyz_0, g_0_yyyy_0_xyyyyyz_1, g_0_yyyy_0_xyyyyz_1, g_0_yyyy_0_xyyyyzz_0, g_0_yyyy_0_xyyyyzz_1, g_0_yyyy_0_xyyyzz_1, g_0_yyyy_0_xyyyzzz_0, g_0_yyyy_0_xyyyzzz_1, g_0_yyyy_0_xyyzzz_1, g_0_yyyy_0_xyyzzzz_0, g_0_yyyy_0_xyyzzzz_1, g_0_yyyy_0_xyzzzz_1, g_0_yyyy_0_xyzzzzz_0, g_0_yyyy_0_xyzzzzz_1, g_0_yyyy_0_xzzzzz_1, g_0_yyyy_0_xzzzzzz_0, g_0_yyyy_0_xzzzzzz_1, g_0_yyyy_0_yyyyyy_1, g_0_yyyy_0_yyyyyyy_0, g_0_yyyy_0_yyyyyyy_1, g_0_yyyy_0_yyyyyyz_0, g_0_yyyy_0_yyyyyyz_1, g_0_yyyy_0_yyyyyz_1, g_0_yyyy_0_yyyyyzz_0, g_0_yyyy_0_yyyyyzz_1, g_0_yyyy_0_yyyyzz_1, g_0_yyyy_0_yyyyzzz_0, g_0_yyyy_0_yyyyzzz_1, g_0_yyyy_0_yyyzzz_1, g_0_yyyy_0_yyyzzzz_0, g_0_yyyy_0_yyyzzzz_1, g_0_yyyy_0_yyzzzz_1, g_0_yyyy_0_yyzzzzz_0, g_0_yyyy_0_yyzzzzz_1, g_0_yyyy_0_yzzzzz_1, g_0_yyyy_0_yzzzzzz_0, g_0_yyyy_0_yzzzzzz_1, g_0_yyyy_0_zzzzzz_1, g_0_yyyy_0_zzzzzzz_0, g_0_yyyy_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyy_0_xxxxxxx_0[i] = 7.0 * g_0_yyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxx_0[i] * pb_x + g_0_yyyy_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxxy_0[i] = 6.0 * g_0_yyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxy_0[i] * pb_x + g_0_yyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxxz_0[i] = 6.0 * g_0_yyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxz_0[i] * pb_x + g_0_yyyy_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxyy_0[i] = 5.0 * g_0_yyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyy_0[i] * pb_x + g_0_yyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxyz_0[i] = 5.0 * g_0_yyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyz_0[i] * pb_x + g_0_yyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxzz_0[i] = 5.0 * g_0_yyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxzz_0[i] * pb_x + g_0_yyyy_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyyy_0[i] = 4.0 * g_0_yyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyy_0[i] * pb_x + g_0_yyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyyz_0[i] = 4.0 * g_0_yyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyz_0[i] * pb_x + g_0_yyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyzz_0[i] = 4.0 * g_0_yyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyzz_0[i] * pb_x + g_0_yyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxzzz_0[i] = 4.0 * g_0_yyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxzzz_0[i] * pb_x + g_0_yyyy_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyyy_0[i] = 3.0 * g_0_yyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyy_0[i] * pb_x + g_0_yyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyyz_0[i] = 3.0 * g_0_yyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyz_0[i] * pb_x + g_0_yyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyzz_0[i] = 3.0 * g_0_yyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyzz_0[i] * pb_x + g_0_yyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyzzz_0[i] = 3.0 * g_0_yyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyzzz_0[i] * pb_x + g_0_yyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxzzzz_0[i] = 3.0 * g_0_yyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxzzzz_0[i] * pb_x + g_0_yyyy_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyyy_0[i] = 2.0 * g_0_yyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyy_0[i] * pb_x + g_0_yyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyyz_0[i] = 2.0 * g_0_yyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyz_0[i] * pb_x + g_0_yyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyzz_0[i] = 2.0 * g_0_yyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyzz_0[i] * pb_x + g_0_yyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyzzz_0[i] = 2.0 * g_0_yyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzzz_0[i] * pb_x + g_0_yyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyzzzz_0[i] = 2.0 * g_0_yyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzzzz_0[i] * pb_x + g_0_yyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxzzzzz_0[i] = 2.0 * g_0_yyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzzzzz_0[i] * pb_x + g_0_yyyy_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyyy_0[i] = g_0_yyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyy_0[i] * pb_x + g_0_yyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyyz_0[i] = g_0_yyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyz_0[i] * pb_x + g_0_yyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyzz_0[i] = g_0_yyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyzz_0[i] * pb_x + g_0_yyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyzzz_0[i] = g_0_yyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzzz_0[i] * pb_x + g_0_yyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyzzzz_0[i] = g_0_yyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzzz_0[i] * pb_x + g_0_yyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyzzzzz_0[i] = g_0_yyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzzzz_0[i] * pb_x + g_0_yyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xzzzzzz_0[i] = g_0_yyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzzzzz_0[i] * pb_x + g_0_yyyy_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyyy_0[i] = g_0_yyyy_0_yyyyyyy_0[i] * pb_x + g_0_yyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyyz_0[i] = g_0_yyyy_0_yyyyyyz_0[i] * pb_x + g_0_yyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyzz_0[i] = g_0_yyyy_0_yyyyyzz_0[i] * pb_x + g_0_yyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyzzz_0[i] = g_0_yyyy_0_yyyyzzz_0[i] * pb_x + g_0_yyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyzzzz_0[i] = g_0_yyyy_0_yyyzzzz_0[i] * pb_x + g_0_yyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyzzzzz_0[i] = g_0_yyyy_0_yyzzzzz_0[i] * pb_x + g_0_yyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yzzzzzz_0[i] = g_0_yyyy_0_yzzzzzz_0[i] * pb_x + g_0_yyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_zzzzzzz_0[i] = g_0_yyyy_0_zzzzzzz_0[i] * pb_x + g_0_yyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 396-432 components of targeted buffer : SHSK

    auto g_0_xyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 396);

    auto g_0_xyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 397);

    auto g_0_xyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 398);

    auto g_0_xyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 399);

    auto g_0_xyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 400);

    auto g_0_xyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 401);

    auto g_0_xyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 402);

    auto g_0_xyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 403);

    auto g_0_xyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 404);

    auto g_0_xyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 405);

    auto g_0_xyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 406);

    auto g_0_xyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 407);

    auto g_0_xyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 408);

    auto g_0_xyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 409);

    auto g_0_xyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 410);

    auto g_0_xyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 411);

    auto g_0_xyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 412);

    auto g_0_xyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 413);

    auto g_0_xyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 414);

    auto g_0_xyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 415);

    auto g_0_xyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 416);

    auto g_0_xyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 417);

    auto g_0_xyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 418);

    auto g_0_xyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 419);

    auto g_0_xyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 420);

    auto g_0_xyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 421);

    auto g_0_xyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 422);

    auto g_0_xyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 423);

    auto g_0_xyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 424);

    auto g_0_xyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 425);

    auto g_0_xyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 426);

    auto g_0_xyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 427);

    auto g_0_xyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 428);

    auto g_0_xyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 429);

    auto g_0_xyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 430);

    auto g_0_xyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 431);

    #pragma omp simd aligned(g_0_xyyy_0_xxxxxxx_0, g_0_xyyy_0_xxxxxxx_1, g_0_xyyy_0_xxxxxxy_0, g_0_xyyy_0_xxxxxxy_1, g_0_xyyy_0_xxxxxyy_0, g_0_xyyy_0_xxxxxyy_1, g_0_xyyy_0_xxxxyyy_0, g_0_xyyy_0_xxxxyyy_1, g_0_xyyy_0_xxxyyyy_0, g_0_xyyy_0_xxxyyyy_1, g_0_xyyy_0_xxyyyyy_0, g_0_xyyy_0_xxyyyyy_1, g_0_xyyy_0_xyyyyyy_0, g_0_xyyy_0_xyyyyyy_1, g_0_xyyyz_0_xxxxxxx_0, g_0_xyyyz_0_xxxxxxy_0, g_0_xyyyz_0_xxxxxxz_0, g_0_xyyyz_0_xxxxxyy_0, g_0_xyyyz_0_xxxxxyz_0, g_0_xyyyz_0_xxxxxzz_0, g_0_xyyyz_0_xxxxyyy_0, g_0_xyyyz_0_xxxxyyz_0, g_0_xyyyz_0_xxxxyzz_0, g_0_xyyyz_0_xxxxzzz_0, g_0_xyyyz_0_xxxyyyy_0, g_0_xyyyz_0_xxxyyyz_0, g_0_xyyyz_0_xxxyyzz_0, g_0_xyyyz_0_xxxyzzz_0, g_0_xyyyz_0_xxxzzzz_0, g_0_xyyyz_0_xxyyyyy_0, g_0_xyyyz_0_xxyyyyz_0, g_0_xyyyz_0_xxyyyzz_0, g_0_xyyyz_0_xxyyzzz_0, g_0_xyyyz_0_xxyzzzz_0, g_0_xyyyz_0_xxzzzzz_0, g_0_xyyyz_0_xyyyyyy_0, g_0_xyyyz_0_xyyyyyz_0, g_0_xyyyz_0_xyyyyzz_0, g_0_xyyyz_0_xyyyzzz_0, g_0_xyyyz_0_xyyzzzz_0, g_0_xyyyz_0_xyzzzzz_0, g_0_xyyyz_0_xzzzzzz_0, g_0_xyyyz_0_yyyyyyy_0, g_0_xyyyz_0_yyyyyyz_0, g_0_xyyyz_0_yyyyyzz_0, g_0_xyyyz_0_yyyyzzz_0, g_0_xyyyz_0_yyyzzzz_0, g_0_xyyyz_0_yyzzzzz_0, g_0_xyyyz_0_yzzzzzz_0, g_0_xyyyz_0_zzzzzzz_0, g_0_yyyz_0_xxxxxxz_0, g_0_yyyz_0_xxxxxxz_1, g_0_yyyz_0_xxxxxyz_0, g_0_yyyz_0_xxxxxyz_1, g_0_yyyz_0_xxxxxz_1, g_0_yyyz_0_xxxxxzz_0, g_0_yyyz_0_xxxxxzz_1, g_0_yyyz_0_xxxxyyz_0, g_0_yyyz_0_xxxxyyz_1, g_0_yyyz_0_xxxxyz_1, g_0_yyyz_0_xxxxyzz_0, g_0_yyyz_0_xxxxyzz_1, g_0_yyyz_0_xxxxzz_1, g_0_yyyz_0_xxxxzzz_0, g_0_yyyz_0_xxxxzzz_1, g_0_yyyz_0_xxxyyyz_0, g_0_yyyz_0_xxxyyyz_1, g_0_yyyz_0_xxxyyz_1, g_0_yyyz_0_xxxyyzz_0, g_0_yyyz_0_xxxyyzz_1, g_0_yyyz_0_xxxyzz_1, g_0_yyyz_0_xxxyzzz_0, g_0_yyyz_0_xxxyzzz_1, g_0_yyyz_0_xxxzzz_1, g_0_yyyz_0_xxxzzzz_0, g_0_yyyz_0_xxxzzzz_1, g_0_yyyz_0_xxyyyyz_0, g_0_yyyz_0_xxyyyyz_1, g_0_yyyz_0_xxyyyz_1, g_0_yyyz_0_xxyyyzz_0, g_0_yyyz_0_xxyyyzz_1, g_0_yyyz_0_xxyyzz_1, g_0_yyyz_0_xxyyzzz_0, g_0_yyyz_0_xxyyzzz_1, g_0_yyyz_0_xxyzzz_1, g_0_yyyz_0_xxyzzzz_0, g_0_yyyz_0_xxyzzzz_1, g_0_yyyz_0_xxzzzz_1, g_0_yyyz_0_xxzzzzz_0, g_0_yyyz_0_xxzzzzz_1, g_0_yyyz_0_xyyyyyz_0, g_0_yyyz_0_xyyyyyz_1, g_0_yyyz_0_xyyyyz_1, g_0_yyyz_0_xyyyyzz_0, g_0_yyyz_0_xyyyyzz_1, g_0_yyyz_0_xyyyzz_1, g_0_yyyz_0_xyyyzzz_0, g_0_yyyz_0_xyyyzzz_1, g_0_yyyz_0_xyyzzz_1, g_0_yyyz_0_xyyzzzz_0, g_0_yyyz_0_xyyzzzz_1, g_0_yyyz_0_xyzzzz_1, g_0_yyyz_0_xyzzzzz_0, g_0_yyyz_0_xyzzzzz_1, g_0_yyyz_0_xzzzzz_1, g_0_yyyz_0_xzzzzzz_0, g_0_yyyz_0_xzzzzzz_1, g_0_yyyz_0_yyyyyyy_0, g_0_yyyz_0_yyyyyyy_1, g_0_yyyz_0_yyyyyyz_0, g_0_yyyz_0_yyyyyyz_1, g_0_yyyz_0_yyyyyz_1, g_0_yyyz_0_yyyyyzz_0, g_0_yyyz_0_yyyyyzz_1, g_0_yyyz_0_yyyyzz_1, g_0_yyyz_0_yyyyzzz_0, g_0_yyyz_0_yyyyzzz_1, g_0_yyyz_0_yyyzzz_1, g_0_yyyz_0_yyyzzzz_0, g_0_yyyz_0_yyyzzzz_1, g_0_yyyz_0_yyzzzz_1, g_0_yyyz_0_yyzzzzz_0, g_0_yyyz_0_yyzzzzz_1, g_0_yyyz_0_yzzzzz_1, g_0_yyyz_0_yzzzzzz_0, g_0_yyyz_0_yzzzzzz_1, g_0_yyyz_0_zzzzzz_1, g_0_yyyz_0_zzzzzzz_0, g_0_yyyz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyz_0_xxxxxxx_0[i] = g_0_xyyy_0_xxxxxxx_0[i] * pb_z + g_0_xyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxxy_0[i] = g_0_xyyy_0_xxxxxxy_0[i] * pb_z + g_0_xyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxxz_0[i] * pb_x + g_0_yyyz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxxyy_0[i] = g_0_xyyy_0_xxxxxyy_0[i] * pb_z + g_0_xyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxyz_0[i] * pb_x + g_0_yyyz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxzz_0[i] * pb_x + g_0_yyyz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxyyy_0[i] = g_0_xyyy_0_xxxxyyy_0[i] * pb_z + g_0_xyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxyyz_0[i] * pb_x + g_0_yyyz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxyzz_0[i] * pb_x + g_0_yyyz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxzzz_0[i] * pb_x + g_0_yyyz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyyyy_0[i] = g_0_xyyy_0_xxxyyyy_0[i] * pb_z + g_0_xyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyyyz_0[i] * pb_x + g_0_yyyz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyyzz_0[i] * pb_x + g_0_yyyz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyzzz_0[i] * pb_x + g_0_yyyz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxzzzz_0[i] * pb_x + g_0_yyyz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyyyy_0[i] = g_0_xyyy_0_xxyyyyy_0[i] * pb_z + g_0_xyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyyyz_0[i] * pb_x + g_0_yyyz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyyzz_0[i] * pb_x + g_0_yyyz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyzzz_0[i] * pb_x + g_0_yyyz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyzzzz_0[i] * pb_x + g_0_yyyz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxzzzzz_0[i] * pb_x + g_0_yyyz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyyyy_0[i] = g_0_xyyy_0_xyyyyyy_0[i] * pb_z + g_0_xyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xyyyyyz_0[i] = g_0_yyyz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyyyz_0[i] * pb_x + g_0_yyyz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyyzz_0[i] = g_0_yyyz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyyzz_0[i] * pb_x + g_0_yyyz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyzzz_0[i] = g_0_yyyz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyzzz_0[i] * pb_x + g_0_yyyz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyzzzz_0[i] = g_0_yyyz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyzzzz_0[i] * pb_x + g_0_yyyz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyzzzzz_0[i] = g_0_yyyz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyzzzzz_0[i] * pb_x + g_0_yyyz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xzzzzzz_0[i] = g_0_yyyz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xzzzzzz_0[i] * pb_x + g_0_yyyz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyyy_0[i] = g_0_yyyz_0_yyyyyyy_0[i] * pb_x + g_0_yyyz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyyz_0[i] = g_0_yyyz_0_yyyyyyz_0[i] * pb_x + g_0_yyyz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyzz_0[i] = g_0_yyyz_0_yyyyyzz_0[i] * pb_x + g_0_yyyz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyzzz_0[i] = g_0_yyyz_0_yyyyzzz_0[i] * pb_x + g_0_yyyz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyzzzz_0[i] = g_0_yyyz_0_yyyzzzz_0[i] * pb_x + g_0_yyyz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyzzzzz_0[i] = g_0_yyyz_0_yyzzzzz_0[i] * pb_x + g_0_yyyz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yzzzzzz_0[i] = g_0_yyyz_0_yzzzzzz_0[i] * pb_x + g_0_yyyz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_zzzzzzz_0[i] = g_0_yyyz_0_zzzzzzz_0[i] * pb_x + g_0_yyyz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 432-468 components of targeted buffer : SHSK

    auto g_0_xyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 432);

    auto g_0_xyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 433);

    auto g_0_xyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 434);

    auto g_0_xyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 435);

    auto g_0_xyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 436);

    auto g_0_xyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 437);

    auto g_0_xyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 438);

    auto g_0_xyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 439);

    auto g_0_xyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 440);

    auto g_0_xyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 441);

    auto g_0_xyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 442);

    auto g_0_xyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 443);

    auto g_0_xyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 444);

    auto g_0_xyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 445);

    auto g_0_xyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 446);

    auto g_0_xyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 447);

    auto g_0_xyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 448);

    auto g_0_xyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 449);

    auto g_0_xyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 450);

    auto g_0_xyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 451);

    auto g_0_xyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 452);

    auto g_0_xyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 453);

    auto g_0_xyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 454);

    auto g_0_xyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 455);

    auto g_0_xyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 456);

    auto g_0_xyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 457);

    auto g_0_xyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 458);

    auto g_0_xyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 459);

    auto g_0_xyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 460);

    auto g_0_xyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 461);

    auto g_0_xyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 462);

    auto g_0_xyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 463);

    auto g_0_xyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 464);

    auto g_0_xyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 465);

    auto g_0_xyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 466);

    auto g_0_xyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 467);

    #pragma omp simd aligned(g_0_xyyzz_0_xxxxxxx_0, g_0_xyyzz_0_xxxxxxy_0, g_0_xyyzz_0_xxxxxxz_0, g_0_xyyzz_0_xxxxxyy_0, g_0_xyyzz_0_xxxxxyz_0, g_0_xyyzz_0_xxxxxzz_0, g_0_xyyzz_0_xxxxyyy_0, g_0_xyyzz_0_xxxxyyz_0, g_0_xyyzz_0_xxxxyzz_0, g_0_xyyzz_0_xxxxzzz_0, g_0_xyyzz_0_xxxyyyy_0, g_0_xyyzz_0_xxxyyyz_0, g_0_xyyzz_0_xxxyyzz_0, g_0_xyyzz_0_xxxyzzz_0, g_0_xyyzz_0_xxxzzzz_0, g_0_xyyzz_0_xxyyyyy_0, g_0_xyyzz_0_xxyyyyz_0, g_0_xyyzz_0_xxyyyzz_0, g_0_xyyzz_0_xxyyzzz_0, g_0_xyyzz_0_xxyzzzz_0, g_0_xyyzz_0_xxzzzzz_0, g_0_xyyzz_0_xyyyyyy_0, g_0_xyyzz_0_xyyyyyz_0, g_0_xyyzz_0_xyyyyzz_0, g_0_xyyzz_0_xyyyzzz_0, g_0_xyyzz_0_xyyzzzz_0, g_0_xyyzz_0_xyzzzzz_0, g_0_xyyzz_0_xzzzzzz_0, g_0_xyyzz_0_yyyyyyy_0, g_0_xyyzz_0_yyyyyyz_0, g_0_xyyzz_0_yyyyyzz_0, g_0_xyyzz_0_yyyyzzz_0, g_0_xyyzz_0_yyyzzzz_0, g_0_xyyzz_0_yyzzzzz_0, g_0_xyyzz_0_yzzzzzz_0, g_0_xyyzz_0_zzzzzzz_0, g_0_yyzz_0_xxxxxx_1, g_0_yyzz_0_xxxxxxx_0, g_0_yyzz_0_xxxxxxx_1, g_0_yyzz_0_xxxxxxy_0, g_0_yyzz_0_xxxxxxy_1, g_0_yyzz_0_xxxxxxz_0, g_0_yyzz_0_xxxxxxz_1, g_0_yyzz_0_xxxxxy_1, g_0_yyzz_0_xxxxxyy_0, g_0_yyzz_0_xxxxxyy_1, g_0_yyzz_0_xxxxxyz_0, g_0_yyzz_0_xxxxxyz_1, g_0_yyzz_0_xxxxxz_1, g_0_yyzz_0_xxxxxzz_0, g_0_yyzz_0_xxxxxzz_1, g_0_yyzz_0_xxxxyy_1, g_0_yyzz_0_xxxxyyy_0, g_0_yyzz_0_xxxxyyy_1, g_0_yyzz_0_xxxxyyz_0, g_0_yyzz_0_xxxxyyz_1, g_0_yyzz_0_xxxxyz_1, g_0_yyzz_0_xxxxyzz_0, g_0_yyzz_0_xxxxyzz_1, g_0_yyzz_0_xxxxzz_1, g_0_yyzz_0_xxxxzzz_0, g_0_yyzz_0_xxxxzzz_1, g_0_yyzz_0_xxxyyy_1, g_0_yyzz_0_xxxyyyy_0, g_0_yyzz_0_xxxyyyy_1, g_0_yyzz_0_xxxyyyz_0, g_0_yyzz_0_xxxyyyz_1, g_0_yyzz_0_xxxyyz_1, g_0_yyzz_0_xxxyyzz_0, g_0_yyzz_0_xxxyyzz_1, g_0_yyzz_0_xxxyzz_1, g_0_yyzz_0_xxxyzzz_0, g_0_yyzz_0_xxxyzzz_1, g_0_yyzz_0_xxxzzz_1, g_0_yyzz_0_xxxzzzz_0, g_0_yyzz_0_xxxzzzz_1, g_0_yyzz_0_xxyyyy_1, g_0_yyzz_0_xxyyyyy_0, g_0_yyzz_0_xxyyyyy_1, g_0_yyzz_0_xxyyyyz_0, g_0_yyzz_0_xxyyyyz_1, g_0_yyzz_0_xxyyyz_1, g_0_yyzz_0_xxyyyzz_0, g_0_yyzz_0_xxyyyzz_1, g_0_yyzz_0_xxyyzz_1, g_0_yyzz_0_xxyyzzz_0, g_0_yyzz_0_xxyyzzz_1, g_0_yyzz_0_xxyzzz_1, g_0_yyzz_0_xxyzzzz_0, g_0_yyzz_0_xxyzzzz_1, g_0_yyzz_0_xxzzzz_1, g_0_yyzz_0_xxzzzzz_0, g_0_yyzz_0_xxzzzzz_1, g_0_yyzz_0_xyyyyy_1, g_0_yyzz_0_xyyyyyy_0, g_0_yyzz_0_xyyyyyy_1, g_0_yyzz_0_xyyyyyz_0, g_0_yyzz_0_xyyyyyz_1, g_0_yyzz_0_xyyyyz_1, g_0_yyzz_0_xyyyyzz_0, g_0_yyzz_0_xyyyyzz_1, g_0_yyzz_0_xyyyzz_1, g_0_yyzz_0_xyyyzzz_0, g_0_yyzz_0_xyyyzzz_1, g_0_yyzz_0_xyyzzz_1, g_0_yyzz_0_xyyzzzz_0, g_0_yyzz_0_xyyzzzz_1, g_0_yyzz_0_xyzzzz_1, g_0_yyzz_0_xyzzzzz_0, g_0_yyzz_0_xyzzzzz_1, g_0_yyzz_0_xzzzzz_1, g_0_yyzz_0_xzzzzzz_0, g_0_yyzz_0_xzzzzzz_1, g_0_yyzz_0_yyyyyy_1, g_0_yyzz_0_yyyyyyy_0, g_0_yyzz_0_yyyyyyy_1, g_0_yyzz_0_yyyyyyz_0, g_0_yyzz_0_yyyyyyz_1, g_0_yyzz_0_yyyyyz_1, g_0_yyzz_0_yyyyyzz_0, g_0_yyzz_0_yyyyyzz_1, g_0_yyzz_0_yyyyzz_1, g_0_yyzz_0_yyyyzzz_0, g_0_yyzz_0_yyyyzzz_1, g_0_yyzz_0_yyyzzz_1, g_0_yyzz_0_yyyzzzz_0, g_0_yyzz_0_yyyzzzz_1, g_0_yyzz_0_yyzzzz_1, g_0_yyzz_0_yyzzzzz_0, g_0_yyzz_0_yyzzzzz_1, g_0_yyzz_0_yzzzzz_1, g_0_yyzz_0_yzzzzzz_0, g_0_yyzz_0_yzzzzzz_1, g_0_yyzz_0_zzzzzz_1, g_0_yyzz_0_zzzzzzz_0, g_0_yyzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxx_0[i] * pb_x + g_0_yyzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxy_0[i] * pb_x + g_0_yyzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxz_0[i] * pb_x + g_0_yyzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxyy_0[i] * pb_x + g_0_yyzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxyz_0[i] * pb_x + g_0_yyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxzz_0[i] * pb_x + g_0_yyzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyyy_0[i] * pb_x + g_0_yyzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyyz_0[i] * pb_x + g_0_yyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyzz_0[i] * pb_x + g_0_yyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxzzz_0[i] * pb_x + g_0_yyzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyyy_0[i] * pb_x + g_0_yyzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyyz_0[i] * pb_x + g_0_yyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyzz_0[i] * pb_x + g_0_yyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyzzz_0[i] * pb_x + g_0_yyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxzzzz_0[i] * pb_x + g_0_yyzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyyy_0[i] * pb_x + g_0_yyzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyyz_0[i] * pb_x + g_0_yyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyzz_0[i] * pb_x + g_0_yyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyzzz_0[i] * pb_x + g_0_yyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyzzzz_0[i] * pb_x + g_0_yyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxzzzzz_0[i] * pb_x + g_0_yyzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyyy_0[i] = g_0_yyzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyyy_0[i] * pb_x + g_0_yyzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyyz_0[i] = g_0_yyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyyz_0[i] * pb_x + g_0_yyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyzz_0[i] = g_0_yyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyzz_0[i] * pb_x + g_0_yyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyzzz_0[i] = g_0_yyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyzzz_0[i] * pb_x + g_0_yyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyzzzz_0[i] = g_0_yyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyzzzz_0[i] * pb_x + g_0_yyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyzzzzz_0[i] = g_0_yyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyzzzzz_0[i] * pb_x + g_0_yyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xzzzzzz_0[i] = g_0_yyzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xzzzzzz_0[i] * pb_x + g_0_yyzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyyy_0[i] = g_0_yyzz_0_yyyyyyy_0[i] * pb_x + g_0_yyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyyz_0[i] = g_0_yyzz_0_yyyyyyz_0[i] * pb_x + g_0_yyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyzz_0[i] = g_0_yyzz_0_yyyyyzz_0[i] * pb_x + g_0_yyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyzzz_0[i] = g_0_yyzz_0_yyyyzzz_0[i] * pb_x + g_0_yyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyzzzz_0[i] = g_0_yyzz_0_yyyzzzz_0[i] * pb_x + g_0_yyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyzzzzz_0[i] = g_0_yyzz_0_yyzzzzz_0[i] * pb_x + g_0_yyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yzzzzzz_0[i] = g_0_yyzz_0_yzzzzzz_0[i] * pb_x + g_0_yyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_zzzzzzz_0[i] = g_0_yyzz_0_zzzzzzz_0[i] * pb_x + g_0_yyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 468-504 components of targeted buffer : SHSK

    auto g_0_xyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 468);

    auto g_0_xyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 469);

    auto g_0_xyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 470);

    auto g_0_xyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 471);

    auto g_0_xyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 472);

    auto g_0_xyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 473);

    auto g_0_xyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 474);

    auto g_0_xyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 475);

    auto g_0_xyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 476);

    auto g_0_xyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 477);

    auto g_0_xyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 478);

    auto g_0_xyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 479);

    auto g_0_xyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 480);

    auto g_0_xyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 481);

    auto g_0_xyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 482);

    auto g_0_xyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 483);

    auto g_0_xyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 484);

    auto g_0_xyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 485);

    auto g_0_xyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 486);

    auto g_0_xyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 487);

    auto g_0_xyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 488);

    auto g_0_xyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 489);

    auto g_0_xyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 490);

    auto g_0_xyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 491);

    auto g_0_xyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 492);

    auto g_0_xyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 493);

    auto g_0_xyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 494);

    auto g_0_xyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 495);

    auto g_0_xyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 496);

    auto g_0_xyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 497);

    auto g_0_xyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 498);

    auto g_0_xyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 499);

    auto g_0_xyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 500);

    auto g_0_xyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 501);

    auto g_0_xyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 502);

    auto g_0_xyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 503);

    #pragma omp simd aligned(g_0_xyzzz_0_xxxxxxx_0, g_0_xyzzz_0_xxxxxxy_0, g_0_xyzzz_0_xxxxxxz_0, g_0_xyzzz_0_xxxxxyy_0, g_0_xyzzz_0_xxxxxyz_0, g_0_xyzzz_0_xxxxxzz_0, g_0_xyzzz_0_xxxxyyy_0, g_0_xyzzz_0_xxxxyyz_0, g_0_xyzzz_0_xxxxyzz_0, g_0_xyzzz_0_xxxxzzz_0, g_0_xyzzz_0_xxxyyyy_0, g_0_xyzzz_0_xxxyyyz_0, g_0_xyzzz_0_xxxyyzz_0, g_0_xyzzz_0_xxxyzzz_0, g_0_xyzzz_0_xxxzzzz_0, g_0_xyzzz_0_xxyyyyy_0, g_0_xyzzz_0_xxyyyyz_0, g_0_xyzzz_0_xxyyyzz_0, g_0_xyzzz_0_xxyyzzz_0, g_0_xyzzz_0_xxyzzzz_0, g_0_xyzzz_0_xxzzzzz_0, g_0_xyzzz_0_xyyyyyy_0, g_0_xyzzz_0_xyyyyyz_0, g_0_xyzzz_0_xyyyyzz_0, g_0_xyzzz_0_xyyyzzz_0, g_0_xyzzz_0_xyyzzzz_0, g_0_xyzzz_0_xyzzzzz_0, g_0_xyzzz_0_xzzzzzz_0, g_0_xyzzz_0_yyyyyyy_0, g_0_xyzzz_0_yyyyyyz_0, g_0_xyzzz_0_yyyyyzz_0, g_0_xyzzz_0_yyyyzzz_0, g_0_xyzzz_0_yyyzzzz_0, g_0_xyzzz_0_yyzzzzz_0, g_0_xyzzz_0_yzzzzzz_0, g_0_xyzzz_0_zzzzzzz_0, g_0_xzzz_0_xxxxxxx_0, g_0_xzzz_0_xxxxxxx_1, g_0_xzzz_0_xxxxxxz_0, g_0_xzzz_0_xxxxxxz_1, g_0_xzzz_0_xxxxxzz_0, g_0_xzzz_0_xxxxxzz_1, g_0_xzzz_0_xxxxzzz_0, g_0_xzzz_0_xxxxzzz_1, g_0_xzzz_0_xxxzzzz_0, g_0_xzzz_0_xxxzzzz_1, g_0_xzzz_0_xxzzzzz_0, g_0_xzzz_0_xxzzzzz_1, g_0_xzzz_0_xzzzzzz_0, g_0_xzzz_0_xzzzzzz_1, g_0_yzzz_0_xxxxxxy_0, g_0_yzzz_0_xxxxxxy_1, g_0_yzzz_0_xxxxxy_1, g_0_yzzz_0_xxxxxyy_0, g_0_yzzz_0_xxxxxyy_1, g_0_yzzz_0_xxxxxyz_0, g_0_yzzz_0_xxxxxyz_1, g_0_yzzz_0_xxxxyy_1, g_0_yzzz_0_xxxxyyy_0, g_0_yzzz_0_xxxxyyy_1, g_0_yzzz_0_xxxxyyz_0, g_0_yzzz_0_xxxxyyz_1, g_0_yzzz_0_xxxxyz_1, g_0_yzzz_0_xxxxyzz_0, g_0_yzzz_0_xxxxyzz_1, g_0_yzzz_0_xxxyyy_1, g_0_yzzz_0_xxxyyyy_0, g_0_yzzz_0_xxxyyyy_1, g_0_yzzz_0_xxxyyyz_0, g_0_yzzz_0_xxxyyyz_1, g_0_yzzz_0_xxxyyz_1, g_0_yzzz_0_xxxyyzz_0, g_0_yzzz_0_xxxyyzz_1, g_0_yzzz_0_xxxyzz_1, g_0_yzzz_0_xxxyzzz_0, g_0_yzzz_0_xxxyzzz_1, g_0_yzzz_0_xxyyyy_1, g_0_yzzz_0_xxyyyyy_0, g_0_yzzz_0_xxyyyyy_1, g_0_yzzz_0_xxyyyyz_0, g_0_yzzz_0_xxyyyyz_1, g_0_yzzz_0_xxyyyz_1, g_0_yzzz_0_xxyyyzz_0, g_0_yzzz_0_xxyyyzz_1, g_0_yzzz_0_xxyyzz_1, g_0_yzzz_0_xxyyzzz_0, g_0_yzzz_0_xxyyzzz_1, g_0_yzzz_0_xxyzzz_1, g_0_yzzz_0_xxyzzzz_0, g_0_yzzz_0_xxyzzzz_1, g_0_yzzz_0_xyyyyy_1, g_0_yzzz_0_xyyyyyy_0, g_0_yzzz_0_xyyyyyy_1, g_0_yzzz_0_xyyyyyz_0, g_0_yzzz_0_xyyyyyz_1, g_0_yzzz_0_xyyyyz_1, g_0_yzzz_0_xyyyyzz_0, g_0_yzzz_0_xyyyyzz_1, g_0_yzzz_0_xyyyzz_1, g_0_yzzz_0_xyyyzzz_0, g_0_yzzz_0_xyyyzzz_1, g_0_yzzz_0_xyyzzz_1, g_0_yzzz_0_xyyzzzz_0, g_0_yzzz_0_xyyzzzz_1, g_0_yzzz_0_xyzzzz_1, g_0_yzzz_0_xyzzzzz_0, g_0_yzzz_0_xyzzzzz_1, g_0_yzzz_0_yyyyyy_1, g_0_yzzz_0_yyyyyyy_0, g_0_yzzz_0_yyyyyyy_1, g_0_yzzz_0_yyyyyyz_0, g_0_yzzz_0_yyyyyyz_1, g_0_yzzz_0_yyyyyz_1, g_0_yzzz_0_yyyyyzz_0, g_0_yzzz_0_yyyyyzz_1, g_0_yzzz_0_yyyyzz_1, g_0_yzzz_0_yyyyzzz_0, g_0_yzzz_0_yyyyzzz_1, g_0_yzzz_0_yyyzzz_1, g_0_yzzz_0_yyyzzzz_0, g_0_yzzz_0_yyyzzzz_1, g_0_yzzz_0_yyzzzz_1, g_0_yzzz_0_yyzzzzz_0, g_0_yzzz_0_yyzzzzz_1, g_0_yzzz_0_yzzzzz_1, g_0_yzzz_0_yzzzzzz_0, g_0_yzzz_0_yzzzzzz_1, g_0_yzzz_0_zzzzzzz_0, g_0_yzzz_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzz_0_xxxxxxx_0[i] = g_0_xzzz_0_xxxxxxx_0[i] * pb_y + g_0_xzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxxy_0[i] * pb_x + g_0_yzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxxz_0[i] = g_0_xzzz_0_xxxxxxz_0[i] * pb_y + g_0_xzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxyy_0[i] * pb_x + g_0_yzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxyz_0[i] * pb_x + g_0_yzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxzz_0[i] = g_0_xzzz_0_xxxxxzz_0[i] * pb_y + g_0_xzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyyy_0[i] * pb_x + g_0_yzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyyz_0[i] * pb_x + g_0_yzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyzz_0[i] * pb_x + g_0_yzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxzzz_0[i] = g_0_xzzz_0_xxxxzzz_0[i] * pb_y + g_0_xzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyyy_0[i] * pb_x + g_0_yzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyyz_0[i] * pb_x + g_0_yzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyzz_0[i] * pb_x + g_0_yzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyzzz_0[i] * pb_x + g_0_yzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxzzzz_0[i] = g_0_xzzz_0_xxxzzzz_0[i] * pb_y + g_0_xzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyyy_0[i] * pb_x + g_0_yzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyyz_0[i] * pb_x + g_0_yzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyzz_0[i] * pb_x + g_0_yzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyzzz_0[i] * pb_x + g_0_yzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyzzzz_0[i] * pb_x + g_0_yzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxzzzzz_0[i] = g_0_xzzz_0_xxzzzzz_0[i] * pb_y + g_0_xzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xyyyyyy_0[i] = g_0_yzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyyy_0[i] * pb_x + g_0_yzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyyyz_0[i] = g_0_yzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyyz_0[i] * pb_x + g_0_yzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyyzz_0[i] = g_0_yzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyzz_0[i] * pb_x + g_0_yzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyzzz_0[i] = g_0_yzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyzzz_0[i] * pb_x + g_0_yzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyzzzz_0[i] = g_0_yzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyzzzz_0[i] * pb_x + g_0_yzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyzzzzz_0[i] = g_0_yzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyzzzzz_0[i] * pb_x + g_0_yzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xzzzzzz_0[i] = g_0_xzzz_0_xzzzzzz_0[i] * pb_y + g_0_xzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_yyyyyyy_0[i] = g_0_yzzz_0_yyyyyyy_0[i] * pb_x + g_0_yzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyyyz_0[i] = g_0_yzzz_0_yyyyyyz_0[i] * pb_x + g_0_yzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyyzz_0[i] = g_0_yzzz_0_yyyyyzz_0[i] * pb_x + g_0_yzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyzzz_0[i] = g_0_yzzz_0_yyyyzzz_0[i] * pb_x + g_0_yzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyzzzz_0[i] = g_0_yzzz_0_yyyzzzz_0[i] * pb_x + g_0_yzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyzzzzz_0[i] = g_0_yzzz_0_yyzzzzz_0[i] * pb_x + g_0_yzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yzzzzzz_0[i] = g_0_yzzz_0_yzzzzzz_0[i] * pb_x + g_0_yzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_zzzzzzz_0[i] = g_0_yzzz_0_zzzzzzz_0[i] * pb_x + g_0_yzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 504-540 components of targeted buffer : SHSK

    auto g_0_xzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 504);

    auto g_0_xzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 505);

    auto g_0_xzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 506);

    auto g_0_xzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 507);

    auto g_0_xzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 508);

    auto g_0_xzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 509);

    auto g_0_xzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 510);

    auto g_0_xzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 511);

    auto g_0_xzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 512);

    auto g_0_xzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 513);

    auto g_0_xzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 514);

    auto g_0_xzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 515);

    auto g_0_xzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 516);

    auto g_0_xzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 517);

    auto g_0_xzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 518);

    auto g_0_xzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 519);

    auto g_0_xzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 520);

    auto g_0_xzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 521);

    auto g_0_xzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 522);

    auto g_0_xzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 523);

    auto g_0_xzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 524);

    auto g_0_xzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 525);

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

    #pragma omp simd aligned(g_0_xzzzz_0_xxxxxxx_0, g_0_xzzzz_0_xxxxxxy_0, g_0_xzzzz_0_xxxxxxz_0, g_0_xzzzz_0_xxxxxyy_0, g_0_xzzzz_0_xxxxxyz_0, g_0_xzzzz_0_xxxxxzz_0, g_0_xzzzz_0_xxxxyyy_0, g_0_xzzzz_0_xxxxyyz_0, g_0_xzzzz_0_xxxxyzz_0, g_0_xzzzz_0_xxxxzzz_0, g_0_xzzzz_0_xxxyyyy_0, g_0_xzzzz_0_xxxyyyz_0, g_0_xzzzz_0_xxxyyzz_0, g_0_xzzzz_0_xxxyzzz_0, g_0_xzzzz_0_xxxzzzz_0, g_0_xzzzz_0_xxyyyyy_0, g_0_xzzzz_0_xxyyyyz_0, g_0_xzzzz_0_xxyyyzz_0, g_0_xzzzz_0_xxyyzzz_0, g_0_xzzzz_0_xxyzzzz_0, g_0_xzzzz_0_xxzzzzz_0, g_0_xzzzz_0_xyyyyyy_0, g_0_xzzzz_0_xyyyyyz_0, g_0_xzzzz_0_xyyyyzz_0, g_0_xzzzz_0_xyyyzzz_0, g_0_xzzzz_0_xyyzzzz_0, g_0_xzzzz_0_xyzzzzz_0, g_0_xzzzz_0_xzzzzzz_0, g_0_xzzzz_0_yyyyyyy_0, g_0_xzzzz_0_yyyyyyz_0, g_0_xzzzz_0_yyyyyzz_0, g_0_xzzzz_0_yyyyzzz_0, g_0_xzzzz_0_yyyzzzz_0, g_0_xzzzz_0_yyzzzzz_0, g_0_xzzzz_0_yzzzzzz_0, g_0_xzzzz_0_zzzzzzz_0, g_0_zzzz_0_xxxxxx_1, g_0_zzzz_0_xxxxxxx_0, g_0_zzzz_0_xxxxxxx_1, g_0_zzzz_0_xxxxxxy_0, g_0_zzzz_0_xxxxxxy_1, g_0_zzzz_0_xxxxxxz_0, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxy_1, g_0_zzzz_0_xxxxxyy_0, g_0_zzzz_0_xxxxxyy_1, g_0_zzzz_0_xxxxxyz_0, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxz_1, g_0_zzzz_0_xxxxxzz_0, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxyy_1, g_0_zzzz_0_xxxxyyy_0, g_0_zzzz_0_xxxxyyy_1, g_0_zzzz_0_xxxxyyz_0, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyz_1, g_0_zzzz_0_xxxxyzz_0, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxzz_1, g_0_zzzz_0_xxxxzzz_0, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxyyy_1, g_0_zzzz_0_xxxyyyy_0, g_0_zzzz_0_xxxyyyy_1, g_0_zzzz_0_xxxyyyz_0, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyz_1, g_0_zzzz_0_xxxyyzz_0, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyzz_1, g_0_zzzz_0_xxxyzzz_0, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxzzz_1, g_0_zzzz_0_xxxzzzz_0, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxyyyy_1, g_0_zzzz_0_xxyyyyy_0, g_0_zzzz_0_xxyyyyy_1, g_0_zzzz_0_xxyyyyz_0, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyz_1, g_0_zzzz_0_xxyyyzz_0, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyzz_1, g_0_zzzz_0_xxyyzzz_0, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyzzz_1, g_0_zzzz_0_xxyzzzz_0, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxzzzz_1, g_0_zzzz_0_xxzzzzz_0, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xyyyyy_1, g_0_zzzz_0_xyyyyyy_0, g_0_zzzz_0_xyyyyyy_1, g_0_zzzz_0_xyyyyyz_0, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyz_1, g_0_zzzz_0_xyyyyzz_0, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyzz_1, g_0_zzzz_0_xyyyzzz_0, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyzzz_1, g_0_zzzz_0_xyyzzzz_0, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyzzzz_1, g_0_zzzz_0_xyzzzzz_0, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xzzzzz_1, g_0_zzzz_0_xzzzzzz_0, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_yyyyyy_1, g_0_zzzz_0_yyyyyyy_0, g_0_zzzz_0_yyyyyyy_1, g_0_zzzz_0_yyyyyyz_0, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyz_1, g_0_zzzz_0_yyyyyzz_0, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyzz_1, g_0_zzzz_0_yyyyzzz_0, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyzzz_1, g_0_zzzz_0_yyyzzzz_0, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyzzzz_1, g_0_zzzz_0_yyzzzzz_0, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yzzzzz_1, g_0_zzzz_0_yzzzzzz_0, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_zzzzzz_1, g_0_zzzz_0_zzzzzzz_0, g_0_zzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzz_0_xxxxxxx_0[i] = 7.0 * g_0_zzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxx_0[i] * pb_x + g_0_zzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_zzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxy_0[i] * pb_x + g_0_zzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxxz_0[i] = 6.0 * g_0_zzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxz_0[i] * pb_x + g_0_zzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_zzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyy_0[i] * pb_x + g_0_zzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_zzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyz_0[i] * pb_x + g_0_zzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxzz_0[i] = 5.0 * g_0_zzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxzz_0[i] * pb_x + g_0_zzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_zzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyy_0[i] * pb_x + g_0_zzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_zzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyz_0[i] * pb_x + g_0_zzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_zzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyzz_0[i] * pb_x + g_0_zzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxzzz_0[i] = 4.0 * g_0_zzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxzzz_0[i] * pb_x + g_0_zzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_zzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyy_0[i] * pb_x + g_0_zzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyz_0[i] * pb_x + g_0_zzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_zzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyzz_0[i] * pb_x + g_0_zzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_zzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzzz_0[i] * pb_x + g_0_zzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_zzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxzzzz_0[i] * pb_x + g_0_zzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_zzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyy_0[i] * pb_x + g_0_zzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_zzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyz_0[i] * pb_x + g_0_zzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_zzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyzz_0[i] * pb_x + g_0_zzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzzz_0[i] * pb_x + g_0_zzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_zzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzzz_0[i] * pb_x + g_0_zzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_zzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzzzzz_0[i] * pb_x + g_0_zzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyyy_0[i] = g_0_zzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyy_0[i] * pb_x + g_0_zzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyyz_0[i] = g_0_zzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyz_0[i] * pb_x + g_0_zzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyzz_0[i] = g_0_zzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyzz_0[i] * pb_x + g_0_zzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyzzz_0[i] = g_0_zzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzzz_0[i] * pb_x + g_0_zzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyzzzz_0[i] = g_0_zzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzzz_0[i] * pb_x + g_0_zzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyzzzzz_0[i] = g_0_zzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzzz_0[i] * pb_x + g_0_zzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xzzzzzz_0[i] = g_0_zzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzzzzz_0[i] * pb_x + g_0_zzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyyy_0[i] = g_0_zzzz_0_yyyyyyy_0[i] * pb_x + g_0_zzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyyz_0[i] = g_0_zzzz_0_yyyyyyz_0[i] * pb_x + g_0_zzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyzz_0[i] = g_0_zzzz_0_yyyyyzz_0[i] * pb_x + g_0_zzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyzzz_0[i] = g_0_zzzz_0_yyyyzzz_0[i] * pb_x + g_0_zzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyzzzz_0[i] = g_0_zzzz_0_yyyzzzz_0[i] * pb_x + g_0_zzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyzzzzz_0[i] = g_0_zzzz_0_yyzzzzz_0[i] * pb_x + g_0_zzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yzzzzzz_0[i] = g_0_zzzz_0_yzzzzzz_0[i] * pb_x + g_0_zzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_zzzzzzz_0[i] = g_0_zzzz_0_zzzzzzz_0[i] * pb_x + g_0_zzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 540-576 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_yyy_0_xxxxxxx_0, g_0_yyy_0_xxxxxxx_1, g_0_yyy_0_xxxxxxy_0, g_0_yyy_0_xxxxxxy_1, g_0_yyy_0_xxxxxxz_0, g_0_yyy_0_xxxxxxz_1, g_0_yyy_0_xxxxxyy_0, g_0_yyy_0_xxxxxyy_1, g_0_yyy_0_xxxxxyz_0, g_0_yyy_0_xxxxxyz_1, g_0_yyy_0_xxxxxzz_0, g_0_yyy_0_xxxxxzz_1, g_0_yyy_0_xxxxyyy_0, g_0_yyy_0_xxxxyyy_1, g_0_yyy_0_xxxxyyz_0, g_0_yyy_0_xxxxyyz_1, g_0_yyy_0_xxxxyzz_0, g_0_yyy_0_xxxxyzz_1, g_0_yyy_0_xxxxzzz_0, g_0_yyy_0_xxxxzzz_1, g_0_yyy_0_xxxyyyy_0, g_0_yyy_0_xxxyyyy_1, g_0_yyy_0_xxxyyyz_0, g_0_yyy_0_xxxyyyz_1, g_0_yyy_0_xxxyyzz_0, g_0_yyy_0_xxxyyzz_1, g_0_yyy_0_xxxyzzz_0, g_0_yyy_0_xxxyzzz_1, g_0_yyy_0_xxxzzzz_0, g_0_yyy_0_xxxzzzz_1, g_0_yyy_0_xxyyyyy_0, g_0_yyy_0_xxyyyyy_1, g_0_yyy_0_xxyyyyz_0, g_0_yyy_0_xxyyyyz_1, g_0_yyy_0_xxyyyzz_0, g_0_yyy_0_xxyyyzz_1, g_0_yyy_0_xxyyzzz_0, g_0_yyy_0_xxyyzzz_1, g_0_yyy_0_xxyzzzz_0, g_0_yyy_0_xxyzzzz_1, g_0_yyy_0_xxzzzzz_0, g_0_yyy_0_xxzzzzz_1, g_0_yyy_0_xyyyyyy_0, g_0_yyy_0_xyyyyyy_1, g_0_yyy_0_xyyyyyz_0, g_0_yyy_0_xyyyyyz_1, g_0_yyy_0_xyyyyzz_0, g_0_yyy_0_xyyyyzz_1, g_0_yyy_0_xyyyzzz_0, g_0_yyy_0_xyyyzzz_1, g_0_yyy_0_xyyzzzz_0, g_0_yyy_0_xyyzzzz_1, g_0_yyy_0_xyzzzzz_0, g_0_yyy_0_xyzzzzz_1, g_0_yyy_0_xzzzzzz_0, g_0_yyy_0_xzzzzzz_1, g_0_yyy_0_yyyyyyy_0, g_0_yyy_0_yyyyyyy_1, g_0_yyy_0_yyyyyyz_0, g_0_yyy_0_yyyyyyz_1, g_0_yyy_0_yyyyyzz_0, g_0_yyy_0_yyyyyzz_1, g_0_yyy_0_yyyyzzz_0, g_0_yyy_0_yyyyzzz_1, g_0_yyy_0_yyyzzzz_0, g_0_yyy_0_yyyzzzz_1, g_0_yyy_0_yyzzzzz_0, g_0_yyy_0_yyzzzzz_1, g_0_yyy_0_yzzzzzz_0, g_0_yyy_0_yzzzzzz_1, g_0_yyy_0_zzzzzzz_0, g_0_yyy_0_zzzzzzz_1, g_0_yyyy_0_xxxxxx_1, g_0_yyyy_0_xxxxxxx_0, g_0_yyyy_0_xxxxxxx_1, g_0_yyyy_0_xxxxxxy_0, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxxz_0, g_0_yyyy_0_xxxxxxz_1, g_0_yyyy_0_xxxxxy_1, g_0_yyyy_0_xxxxxyy_0, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxxyz_0, g_0_yyyy_0_xxxxxyz_1, g_0_yyyy_0_xxxxxz_1, g_0_yyyy_0_xxxxxzz_0, g_0_yyyy_0_xxxxxzz_1, g_0_yyyy_0_xxxxyy_1, g_0_yyyy_0_xxxxyyy_0, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxxyyz_0, g_0_yyyy_0_xxxxyyz_1, g_0_yyyy_0_xxxxyz_1, g_0_yyyy_0_xxxxyzz_0, g_0_yyyy_0_xxxxyzz_1, g_0_yyyy_0_xxxxzz_1, g_0_yyyy_0_xxxxzzz_0, g_0_yyyy_0_xxxxzzz_1, g_0_yyyy_0_xxxyyy_1, g_0_yyyy_0_xxxyyyy_0, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxxyyyz_0, g_0_yyyy_0_xxxyyyz_1, g_0_yyyy_0_xxxyyz_1, g_0_yyyy_0_xxxyyzz_0, g_0_yyyy_0_xxxyyzz_1, g_0_yyyy_0_xxxyzz_1, g_0_yyyy_0_xxxyzzz_0, g_0_yyyy_0_xxxyzzz_1, g_0_yyyy_0_xxxzzz_1, g_0_yyyy_0_xxxzzzz_0, g_0_yyyy_0_xxxzzzz_1, g_0_yyyy_0_xxyyyy_1, g_0_yyyy_0_xxyyyyy_0, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xxyyyyz_0, g_0_yyyy_0_xxyyyyz_1, g_0_yyyy_0_xxyyyz_1, g_0_yyyy_0_xxyyyzz_0, g_0_yyyy_0_xxyyyzz_1, g_0_yyyy_0_xxyyzz_1, g_0_yyyy_0_xxyyzzz_0, g_0_yyyy_0_xxyyzzz_1, g_0_yyyy_0_xxyzzz_1, g_0_yyyy_0_xxyzzzz_0, g_0_yyyy_0_xxyzzzz_1, g_0_yyyy_0_xxzzzz_1, g_0_yyyy_0_xxzzzzz_0, g_0_yyyy_0_xxzzzzz_1, g_0_yyyy_0_xyyyyy_1, g_0_yyyy_0_xyyyyyy_0, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_xyyyyyz_0, g_0_yyyy_0_xyyyyyz_1, g_0_yyyy_0_xyyyyz_1, g_0_yyyy_0_xyyyyzz_0, g_0_yyyy_0_xyyyyzz_1, g_0_yyyy_0_xyyyzz_1, g_0_yyyy_0_xyyyzzz_0, g_0_yyyy_0_xyyyzzz_1, g_0_yyyy_0_xyyzzz_1, g_0_yyyy_0_xyyzzzz_0, g_0_yyyy_0_xyyzzzz_1, g_0_yyyy_0_xyzzzz_1, g_0_yyyy_0_xyzzzzz_0, g_0_yyyy_0_xyzzzzz_1, g_0_yyyy_0_xzzzzz_1, g_0_yyyy_0_xzzzzzz_0, g_0_yyyy_0_xzzzzzz_1, g_0_yyyy_0_yyyyyy_1, g_0_yyyy_0_yyyyyyy_0, g_0_yyyy_0_yyyyyyy_1, g_0_yyyy_0_yyyyyyz_0, g_0_yyyy_0_yyyyyyz_1, g_0_yyyy_0_yyyyyz_1, g_0_yyyy_0_yyyyyzz_0, g_0_yyyy_0_yyyyyzz_1, g_0_yyyy_0_yyyyzz_1, g_0_yyyy_0_yyyyzzz_0, g_0_yyyy_0_yyyyzzz_1, g_0_yyyy_0_yyyzzz_1, g_0_yyyy_0_yyyzzzz_0, g_0_yyyy_0_yyyzzzz_1, g_0_yyyy_0_yyzzzz_1, g_0_yyyy_0_yyzzzzz_0, g_0_yyyy_0_yyzzzzz_1, g_0_yyyy_0_yzzzzz_1, g_0_yyyy_0_yzzzzzz_0, g_0_yyyy_0_yzzzzzz_1, g_0_yyyy_0_zzzzzz_1, g_0_yyyy_0_zzzzzzz_0, g_0_yyyy_0_zzzzzzz_1, g_0_yyyyy_0_xxxxxxx_0, g_0_yyyyy_0_xxxxxxy_0, g_0_yyyyy_0_xxxxxxz_0, g_0_yyyyy_0_xxxxxyy_0, g_0_yyyyy_0_xxxxxyz_0, g_0_yyyyy_0_xxxxxzz_0, g_0_yyyyy_0_xxxxyyy_0, g_0_yyyyy_0_xxxxyyz_0, g_0_yyyyy_0_xxxxyzz_0, g_0_yyyyy_0_xxxxzzz_0, g_0_yyyyy_0_xxxyyyy_0, g_0_yyyyy_0_xxxyyyz_0, g_0_yyyyy_0_xxxyyzz_0, g_0_yyyyy_0_xxxyzzz_0, g_0_yyyyy_0_xxxzzzz_0, g_0_yyyyy_0_xxyyyyy_0, g_0_yyyyy_0_xxyyyyz_0, g_0_yyyyy_0_xxyyyzz_0, g_0_yyyyy_0_xxyyzzz_0, g_0_yyyyy_0_xxyzzzz_0, g_0_yyyyy_0_xxzzzzz_0, g_0_yyyyy_0_xyyyyyy_0, g_0_yyyyy_0_xyyyyyz_0, g_0_yyyyy_0_xyyyyzz_0, g_0_yyyyy_0_xyyyzzz_0, g_0_yyyyy_0_xyyzzzz_0, g_0_yyyyy_0_xyzzzzz_0, g_0_yyyyy_0_xzzzzzz_0, g_0_yyyyy_0_yyyyyyy_0, g_0_yyyyy_0_yyyyyyz_0, g_0_yyyyy_0_yyyyyzz_0, g_0_yyyyy_0_yyyyzzz_0, g_0_yyyyy_0_yyyzzzz_0, g_0_yyyyy_0_yyzzzzz_0, g_0_yyyyy_0_yzzzzzz_0, g_0_yyyyy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyy_0_xxxxxxx_0[i] = 4.0 * g_0_yyy_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxxx_0[i] * pb_y + g_0_yyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxxy_0[i] = 4.0 * g_0_yyy_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxy_0[i] * pb_y + g_0_yyyy_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxxz_0[i] = 4.0 * g_0_yyy_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxxz_0[i] * pb_y + g_0_yyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxyy_0[i] = 4.0 * g_0_yyy_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyy_0[i] * pb_y + g_0_yyyy_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxyz_0[i] = 4.0 * g_0_yyy_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyz_0[i] * pb_y + g_0_yyyy_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxzz_0[i] = 4.0 * g_0_yyy_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxzz_0[i] * pb_y + g_0_yyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyyy_0[i] = 4.0 * g_0_yyy_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyy_0[i] * pb_y + g_0_yyyy_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyyz_0[i] = 4.0 * g_0_yyy_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyz_0[i] * pb_y + g_0_yyyy_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyzz_0[i] = 4.0 * g_0_yyy_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyzz_0[i] * pb_y + g_0_yyyy_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxzzz_0[i] = 4.0 * g_0_yyy_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxzzz_0[i] * pb_y + g_0_yyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyyy_0[i] = 4.0 * g_0_yyy_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyy_0[i] * pb_y + g_0_yyyy_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyyz_0[i] = 4.0 * g_0_yyy_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyz_0[i] * pb_y + g_0_yyyy_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyzz_0[i] = 4.0 * g_0_yyy_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyzz_0[i] * pb_y + g_0_yyyy_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyzzz_0[i] = 4.0 * g_0_yyy_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyzzz_0[i] * pb_y + g_0_yyyy_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxzzzz_0[i] = 4.0 * g_0_yyy_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxzzzz_0[i] * pb_y + g_0_yyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyyy_0[i] = 4.0 * g_0_yyy_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyy_0[i] * pb_y + g_0_yyyy_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyyz_0[i] = 4.0 * g_0_yyy_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyz_0[i] * pb_y + g_0_yyyy_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyzz_0[i] = 4.0 * g_0_yyy_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyzz_0[i] * pb_y + g_0_yyyy_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyzzz_0[i] = 4.0 * g_0_yyy_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzzz_0[i] * pb_y + g_0_yyyy_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyzzzz_0[i] = 4.0 * g_0_yyy_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzzzz_0[i] * pb_y + g_0_yyyy_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxzzzzz_0[i] = 4.0 * g_0_yyy_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzzzzz_0[i] * pb_y + g_0_yyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyyy_0[i] = 4.0 * g_0_yyy_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_yyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyy_0[i] * pb_y + g_0_yyyy_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyyz_0[i] = 4.0 * g_0_yyy_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyz_0[i] * pb_y + g_0_yyyy_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyzz_0[i] = 4.0 * g_0_yyy_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyzz_0[i] * pb_y + g_0_yyyy_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyzzz_0[i] = 4.0 * g_0_yyy_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzzz_0[i] * pb_y + g_0_yyyy_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyzzzz_0[i] = 4.0 * g_0_yyy_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzzz_0[i] * pb_y + g_0_yyyy_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyzzzzz_0[i] = 4.0 * g_0_yyy_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzzzz_0[i] * pb_y + g_0_yyyy_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xzzzzzz_0[i] = 4.0 * g_0_yyy_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzzzzz_0[i] * pb_y + g_0_yyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyyy_0[i] = 4.0 * g_0_yyy_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyyy_1[i] * fti_ab_0 + 7.0 * g_0_yyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyyy_0[i] * pb_y + g_0_yyyy_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyyz_0[i] = 4.0 * g_0_yyy_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyyz_0[i] * pb_y + g_0_yyyy_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyzz_0[i] = 4.0 * g_0_yyy_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyzz_0[i] * pb_y + g_0_yyyy_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyzzz_0[i] = 4.0 * g_0_yyy_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyzzz_0[i] * pb_y + g_0_yyyy_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyzzzz_0[i] = 4.0 * g_0_yyy_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyzzzz_0[i] * pb_y + g_0_yyyy_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyzzzzz_0[i] = 4.0 * g_0_yyy_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzzzzz_0[i] * pb_y + g_0_yyyy_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yzzzzzz_0[i] = 4.0 * g_0_yyy_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzzzzzz_0[i] * pb_y + g_0_yyyy_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_zzzzzzz_0[i] = 4.0 * g_0_yyy_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzzzzz_0[i] * pb_y + g_0_yyyy_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 576-612 components of targeted buffer : SHSK

    auto g_0_yyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 576);

    auto g_0_yyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 577);

    auto g_0_yyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 578);

    auto g_0_yyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 579);

    auto g_0_yyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 580);

    auto g_0_yyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 581);

    auto g_0_yyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 582);

    auto g_0_yyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 583);

    auto g_0_yyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 584);

    auto g_0_yyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 585);

    auto g_0_yyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 586);

    auto g_0_yyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 587);

    auto g_0_yyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 588);

    auto g_0_yyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 589);

    auto g_0_yyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 590);

    auto g_0_yyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 591);

    auto g_0_yyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 592);

    auto g_0_yyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 593);

    auto g_0_yyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 594);

    auto g_0_yyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 595);

    auto g_0_yyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 596);

    auto g_0_yyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 597);

    auto g_0_yyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 598);

    auto g_0_yyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 599);

    auto g_0_yyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 600);

    auto g_0_yyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 601);

    auto g_0_yyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 602);

    auto g_0_yyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 603);

    auto g_0_yyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 604);

    auto g_0_yyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 605);

    auto g_0_yyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 606);

    auto g_0_yyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 607);

    auto g_0_yyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 608);

    auto g_0_yyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 609);

    auto g_0_yyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 610);

    auto g_0_yyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 611);

    #pragma omp simd aligned(g_0_yyyy_0_xxxxxx_1, g_0_yyyy_0_xxxxxxx_0, g_0_yyyy_0_xxxxxxx_1, g_0_yyyy_0_xxxxxxy_0, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxxz_0, g_0_yyyy_0_xxxxxxz_1, g_0_yyyy_0_xxxxxy_1, g_0_yyyy_0_xxxxxyy_0, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxxyz_0, g_0_yyyy_0_xxxxxyz_1, g_0_yyyy_0_xxxxxz_1, g_0_yyyy_0_xxxxxzz_0, g_0_yyyy_0_xxxxxzz_1, g_0_yyyy_0_xxxxyy_1, g_0_yyyy_0_xxxxyyy_0, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxxyyz_0, g_0_yyyy_0_xxxxyyz_1, g_0_yyyy_0_xxxxyz_1, g_0_yyyy_0_xxxxyzz_0, g_0_yyyy_0_xxxxyzz_1, g_0_yyyy_0_xxxxzz_1, g_0_yyyy_0_xxxxzzz_0, g_0_yyyy_0_xxxxzzz_1, g_0_yyyy_0_xxxyyy_1, g_0_yyyy_0_xxxyyyy_0, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxxyyyz_0, g_0_yyyy_0_xxxyyyz_1, g_0_yyyy_0_xxxyyz_1, g_0_yyyy_0_xxxyyzz_0, g_0_yyyy_0_xxxyyzz_1, g_0_yyyy_0_xxxyzz_1, g_0_yyyy_0_xxxyzzz_0, g_0_yyyy_0_xxxyzzz_1, g_0_yyyy_0_xxxzzz_1, g_0_yyyy_0_xxxzzzz_0, g_0_yyyy_0_xxxzzzz_1, g_0_yyyy_0_xxyyyy_1, g_0_yyyy_0_xxyyyyy_0, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xxyyyyz_0, g_0_yyyy_0_xxyyyyz_1, g_0_yyyy_0_xxyyyz_1, g_0_yyyy_0_xxyyyzz_0, g_0_yyyy_0_xxyyyzz_1, g_0_yyyy_0_xxyyzz_1, g_0_yyyy_0_xxyyzzz_0, g_0_yyyy_0_xxyyzzz_1, g_0_yyyy_0_xxyzzz_1, g_0_yyyy_0_xxyzzzz_0, g_0_yyyy_0_xxyzzzz_1, g_0_yyyy_0_xxzzzz_1, g_0_yyyy_0_xxzzzzz_0, g_0_yyyy_0_xxzzzzz_1, g_0_yyyy_0_xyyyyy_1, g_0_yyyy_0_xyyyyyy_0, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_xyyyyyz_0, g_0_yyyy_0_xyyyyyz_1, g_0_yyyy_0_xyyyyz_1, g_0_yyyy_0_xyyyyzz_0, g_0_yyyy_0_xyyyyzz_1, g_0_yyyy_0_xyyyzz_1, g_0_yyyy_0_xyyyzzz_0, g_0_yyyy_0_xyyyzzz_1, g_0_yyyy_0_xyyzzz_1, g_0_yyyy_0_xyyzzzz_0, g_0_yyyy_0_xyyzzzz_1, g_0_yyyy_0_xyzzzz_1, g_0_yyyy_0_xyzzzzz_0, g_0_yyyy_0_xyzzzzz_1, g_0_yyyy_0_xzzzzz_1, g_0_yyyy_0_xzzzzzz_0, g_0_yyyy_0_xzzzzzz_1, g_0_yyyy_0_yyyyyy_1, g_0_yyyy_0_yyyyyyy_0, g_0_yyyy_0_yyyyyyy_1, g_0_yyyy_0_yyyyyyz_0, g_0_yyyy_0_yyyyyyz_1, g_0_yyyy_0_yyyyyz_1, g_0_yyyy_0_yyyyyzz_0, g_0_yyyy_0_yyyyyzz_1, g_0_yyyy_0_yyyyzz_1, g_0_yyyy_0_yyyyzzz_0, g_0_yyyy_0_yyyyzzz_1, g_0_yyyy_0_yyyzzz_1, g_0_yyyy_0_yyyzzzz_0, g_0_yyyy_0_yyyzzzz_1, g_0_yyyy_0_yyzzzz_1, g_0_yyyy_0_yyzzzzz_0, g_0_yyyy_0_yyzzzzz_1, g_0_yyyy_0_yzzzzz_1, g_0_yyyy_0_yzzzzzz_0, g_0_yyyy_0_yzzzzzz_1, g_0_yyyy_0_zzzzzz_1, g_0_yyyy_0_zzzzzzz_0, g_0_yyyy_0_zzzzzzz_1, g_0_yyyyz_0_xxxxxxx_0, g_0_yyyyz_0_xxxxxxy_0, g_0_yyyyz_0_xxxxxxz_0, g_0_yyyyz_0_xxxxxyy_0, g_0_yyyyz_0_xxxxxyz_0, g_0_yyyyz_0_xxxxxzz_0, g_0_yyyyz_0_xxxxyyy_0, g_0_yyyyz_0_xxxxyyz_0, g_0_yyyyz_0_xxxxyzz_0, g_0_yyyyz_0_xxxxzzz_0, g_0_yyyyz_0_xxxyyyy_0, g_0_yyyyz_0_xxxyyyz_0, g_0_yyyyz_0_xxxyyzz_0, g_0_yyyyz_0_xxxyzzz_0, g_0_yyyyz_0_xxxzzzz_0, g_0_yyyyz_0_xxyyyyy_0, g_0_yyyyz_0_xxyyyyz_0, g_0_yyyyz_0_xxyyyzz_0, g_0_yyyyz_0_xxyyzzz_0, g_0_yyyyz_0_xxyzzzz_0, g_0_yyyyz_0_xxzzzzz_0, g_0_yyyyz_0_xyyyyyy_0, g_0_yyyyz_0_xyyyyyz_0, g_0_yyyyz_0_xyyyyzz_0, g_0_yyyyz_0_xyyyzzz_0, g_0_yyyyz_0_xyyzzzz_0, g_0_yyyyz_0_xyzzzzz_0, g_0_yyyyz_0_xzzzzzz_0, g_0_yyyyz_0_yyyyyyy_0, g_0_yyyyz_0_yyyyyyz_0, g_0_yyyyz_0_yyyyyzz_0, g_0_yyyyz_0_yyyyzzz_0, g_0_yyyyz_0_yyyzzzz_0, g_0_yyyyz_0_yyzzzzz_0, g_0_yyyyz_0_yzzzzzz_0, g_0_yyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyz_0_xxxxxxx_0[i] = g_0_yyyy_0_xxxxxxx_0[i] * pb_z + g_0_yyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxxy_0[i] = g_0_yyyy_0_xxxxxxy_0[i] * pb_z + g_0_yyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxxz_0[i] = g_0_yyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxz_0[i] * pb_z + g_0_yyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxyy_0[i] = g_0_yyyy_0_xxxxxyy_0[i] * pb_z + g_0_yyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxyz_0[i] = g_0_yyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyz_0[i] * pb_z + g_0_yyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_yyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxzz_0[i] * pb_z + g_0_yyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyyy_0[i] = g_0_yyyy_0_xxxxyyy_0[i] * pb_z + g_0_yyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyyz_0[i] = g_0_yyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyz_0[i] * pb_z + g_0_yyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_yyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyzz_0[i] * pb_z + g_0_yyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_yyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxzzz_0[i] * pb_z + g_0_yyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyyy_0[i] = g_0_yyyy_0_xxxyyyy_0[i] * pb_z + g_0_yyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyyz_0[i] = g_0_yyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyz_0[i] * pb_z + g_0_yyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_yyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyzz_0[i] * pb_z + g_0_yyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyzzz_0[i] * pb_z + g_0_yyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_yyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxzzzz_0[i] * pb_z + g_0_yyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyyy_0[i] = g_0_yyyy_0_xxyyyyy_0[i] * pb_z + g_0_yyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyyz_0[i] = g_0_yyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyz_0[i] * pb_z + g_0_yyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyzz_0[i] * pb_z + g_0_yyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_yyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzzz_0[i] * pb_z + g_0_yyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_yyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzzzz_0[i] * pb_z + g_0_yyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_yyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzzzzz_0[i] * pb_z + g_0_yyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyyy_0[i] = g_0_yyyy_0_xyyyyyy_0[i] * pb_z + g_0_yyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyyz_0[i] = g_0_yyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyz_0[i] * pb_z + g_0_yyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_yyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyzz_0[i] * pb_z + g_0_yyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_yyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzzz_0[i] * pb_z + g_0_yyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_yyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzzz_0[i] * pb_z + g_0_yyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_yyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzzzz_0[i] * pb_z + g_0_yyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_yyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzzzzz_0[i] * pb_z + g_0_yyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyyy_0[i] = g_0_yyyy_0_yyyyyyy_0[i] * pb_z + g_0_yyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyyz_0[i] = g_0_yyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyyz_0[i] * pb_z + g_0_yyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_yyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyzz_0[i] * pb_z + g_0_yyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_yyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyzzz_0[i] * pb_z + g_0_yyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_yyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyzzzz_0[i] * pb_z + g_0_yyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_yyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzzzzz_0[i] * pb_z + g_0_yyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_yyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzzzzzz_0[i] * pb_z + g_0_yyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_yyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_zzzzzzz_0[i] * pb_z + g_0_yyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 612-648 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_yyy_0_xxxxxxy_0, g_0_yyy_0_xxxxxxy_1, g_0_yyy_0_xxxxxyy_0, g_0_yyy_0_xxxxxyy_1, g_0_yyy_0_xxxxyyy_0, g_0_yyy_0_xxxxyyy_1, g_0_yyy_0_xxxyyyy_0, g_0_yyy_0_xxxyyyy_1, g_0_yyy_0_xxyyyyy_0, g_0_yyy_0_xxyyyyy_1, g_0_yyy_0_xyyyyyy_0, g_0_yyy_0_xyyyyyy_1, g_0_yyy_0_yyyyyyy_0, g_0_yyy_0_yyyyyyy_1, g_0_yyyz_0_xxxxxxy_0, g_0_yyyz_0_xxxxxxy_1, g_0_yyyz_0_xxxxxyy_0, g_0_yyyz_0_xxxxxyy_1, g_0_yyyz_0_xxxxyyy_0, g_0_yyyz_0_xxxxyyy_1, g_0_yyyz_0_xxxyyyy_0, g_0_yyyz_0_xxxyyyy_1, g_0_yyyz_0_xxyyyyy_0, g_0_yyyz_0_xxyyyyy_1, g_0_yyyz_0_xyyyyyy_0, g_0_yyyz_0_xyyyyyy_1, g_0_yyyz_0_yyyyyyy_0, g_0_yyyz_0_yyyyyyy_1, g_0_yyyzz_0_xxxxxxx_0, g_0_yyyzz_0_xxxxxxy_0, g_0_yyyzz_0_xxxxxxz_0, g_0_yyyzz_0_xxxxxyy_0, g_0_yyyzz_0_xxxxxyz_0, g_0_yyyzz_0_xxxxxzz_0, g_0_yyyzz_0_xxxxyyy_0, g_0_yyyzz_0_xxxxyyz_0, g_0_yyyzz_0_xxxxyzz_0, g_0_yyyzz_0_xxxxzzz_0, g_0_yyyzz_0_xxxyyyy_0, g_0_yyyzz_0_xxxyyyz_0, g_0_yyyzz_0_xxxyyzz_0, g_0_yyyzz_0_xxxyzzz_0, g_0_yyyzz_0_xxxzzzz_0, g_0_yyyzz_0_xxyyyyy_0, g_0_yyyzz_0_xxyyyyz_0, g_0_yyyzz_0_xxyyyzz_0, g_0_yyyzz_0_xxyyzzz_0, g_0_yyyzz_0_xxyzzzz_0, g_0_yyyzz_0_xxzzzzz_0, g_0_yyyzz_0_xyyyyyy_0, g_0_yyyzz_0_xyyyyyz_0, g_0_yyyzz_0_xyyyyzz_0, g_0_yyyzz_0_xyyyzzz_0, g_0_yyyzz_0_xyyzzzz_0, g_0_yyyzz_0_xyzzzzz_0, g_0_yyyzz_0_xzzzzzz_0, g_0_yyyzz_0_yyyyyyy_0, g_0_yyyzz_0_yyyyyyz_0, g_0_yyyzz_0_yyyyyzz_0, g_0_yyyzz_0_yyyyzzz_0, g_0_yyyzz_0_yyyzzzz_0, g_0_yyyzz_0_yyzzzzz_0, g_0_yyyzz_0_yzzzzzz_0, g_0_yyyzz_0_zzzzzzz_0, g_0_yyzz_0_xxxxxxx_0, g_0_yyzz_0_xxxxxxx_1, g_0_yyzz_0_xxxxxxz_0, g_0_yyzz_0_xxxxxxz_1, g_0_yyzz_0_xxxxxyz_0, g_0_yyzz_0_xxxxxyz_1, g_0_yyzz_0_xxxxxz_1, g_0_yyzz_0_xxxxxzz_0, g_0_yyzz_0_xxxxxzz_1, g_0_yyzz_0_xxxxyyz_0, g_0_yyzz_0_xxxxyyz_1, g_0_yyzz_0_xxxxyz_1, g_0_yyzz_0_xxxxyzz_0, g_0_yyzz_0_xxxxyzz_1, g_0_yyzz_0_xxxxzz_1, g_0_yyzz_0_xxxxzzz_0, g_0_yyzz_0_xxxxzzz_1, g_0_yyzz_0_xxxyyyz_0, g_0_yyzz_0_xxxyyyz_1, g_0_yyzz_0_xxxyyz_1, g_0_yyzz_0_xxxyyzz_0, g_0_yyzz_0_xxxyyzz_1, g_0_yyzz_0_xxxyzz_1, g_0_yyzz_0_xxxyzzz_0, g_0_yyzz_0_xxxyzzz_1, g_0_yyzz_0_xxxzzz_1, g_0_yyzz_0_xxxzzzz_0, g_0_yyzz_0_xxxzzzz_1, g_0_yyzz_0_xxyyyyz_0, g_0_yyzz_0_xxyyyyz_1, g_0_yyzz_0_xxyyyz_1, g_0_yyzz_0_xxyyyzz_0, g_0_yyzz_0_xxyyyzz_1, g_0_yyzz_0_xxyyzz_1, g_0_yyzz_0_xxyyzzz_0, g_0_yyzz_0_xxyyzzz_1, g_0_yyzz_0_xxyzzz_1, g_0_yyzz_0_xxyzzzz_0, g_0_yyzz_0_xxyzzzz_1, g_0_yyzz_0_xxzzzz_1, g_0_yyzz_0_xxzzzzz_0, g_0_yyzz_0_xxzzzzz_1, g_0_yyzz_0_xyyyyyz_0, g_0_yyzz_0_xyyyyyz_1, g_0_yyzz_0_xyyyyz_1, g_0_yyzz_0_xyyyyzz_0, g_0_yyzz_0_xyyyyzz_1, g_0_yyzz_0_xyyyzz_1, g_0_yyzz_0_xyyyzzz_0, g_0_yyzz_0_xyyyzzz_1, g_0_yyzz_0_xyyzzz_1, g_0_yyzz_0_xyyzzzz_0, g_0_yyzz_0_xyyzzzz_1, g_0_yyzz_0_xyzzzz_1, g_0_yyzz_0_xyzzzzz_0, g_0_yyzz_0_xyzzzzz_1, g_0_yyzz_0_xzzzzz_1, g_0_yyzz_0_xzzzzzz_0, g_0_yyzz_0_xzzzzzz_1, g_0_yyzz_0_yyyyyyz_0, g_0_yyzz_0_yyyyyyz_1, g_0_yyzz_0_yyyyyz_1, g_0_yyzz_0_yyyyyzz_0, g_0_yyzz_0_yyyyyzz_1, g_0_yyzz_0_yyyyzz_1, g_0_yyzz_0_yyyyzzz_0, g_0_yyzz_0_yyyyzzz_1, g_0_yyzz_0_yyyzzz_1, g_0_yyzz_0_yyyzzzz_0, g_0_yyzz_0_yyyzzzz_1, g_0_yyzz_0_yyzzzz_1, g_0_yyzz_0_yyzzzzz_0, g_0_yyzz_0_yyzzzzz_1, g_0_yyzz_0_yzzzzz_1, g_0_yyzz_0_yzzzzzz_0, g_0_yyzz_0_yzzzzzz_1, g_0_yyzz_0_zzzzzz_1, g_0_yyzz_0_zzzzzzz_0, g_0_yyzz_0_zzzzzzz_1, g_0_yzz_0_xxxxxxx_0, g_0_yzz_0_xxxxxxx_1, g_0_yzz_0_xxxxxxz_0, g_0_yzz_0_xxxxxxz_1, g_0_yzz_0_xxxxxyz_0, g_0_yzz_0_xxxxxyz_1, g_0_yzz_0_xxxxxzz_0, g_0_yzz_0_xxxxxzz_1, g_0_yzz_0_xxxxyyz_0, g_0_yzz_0_xxxxyyz_1, g_0_yzz_0_xxxxyzz_0, g_0_yzz_0_xxxxyzz_1, g_0_yzz_0_xxxxzzz_0, g_0_yzz_0_xxxxzzz_1, g_0_yzz_0_xxxyyyz_0, g_0_yzz_0_xxxyyyz_1, g_0_yzz_0_xxxyyzz_0, g_0_yzz_0_xxxyyzz_1, g_0_yzz_0_xxxyzzz_0, g_0_yzz_0_xxxyzzz_1, g_0_yzz_0_xxxzzzz_0, g_0_yzz_0_xxxzzzz_1, g_0_yzz_0_xxyyyyz_0, g_0_yzz_0_xxyyyyz_1, g_0_yzz_0_xxyyyzz_0, g_0_yzz_0_xxyyyzz_1, g_0_yzz_0_xxyyzzz_0, g_0_yzz_0_xxyyzzz_1, g_0_yzz_0_xxyzzzz_0, g_0_yzz_0_xxyzzzz_1, g_0_yzz_0_xxzzzzz_0, g_0_yzz_0_xxzzzzz_1, g_0_yzz_0_xyyyyyz_0, g_0_yzz_0_xyyyyyz_1, g_0_yzz_0_xyyyyzz_0, g_0_yzz_0_xyyyyzz_1, g_0_yzz_0_xyyyzzz_0, g_0_yzz_0_xyyyzzz_1, g_0_yzz_0_xyyzzzz_0, g_0_yzz_0_xyyzzzz_1, g_0_yzz_0_xyzzzzz_0, g_0_yzz_0_xyzzzzz_1, g_0_yzz_0_xzzzzzz_0, g_0_yzz_0_xzzzzzz_1, g_0_yzz_0_yyyyyyz_0, g_0_yzz_0_yyyyyyz_1, g_0_yzz_0_yyyyyzz_0, g_0_yzz_0_yyyyyzz_1, g_0_yzz_0_yyyyzzz_0, g_0_yzz_0_yyyyzzz_1, g_0_yzz_0_yyyzzzz_0, g_0_yzz_0_yyyzzzz_1, g_0_yzz_0_yyzzzzz_0, g_0_yzz_0_yyzzzzz_1, g_0_yzz_0_yzzzzzz_0, g_0_yzz_0_yzzzzzz_1, g_0_yzz_0_zzzzzzz_0, g_0_yzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzz_0_xxxxxxx_0[i] = 2.0 * g_0_yzz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxx_0[i] * pb_y + g_0_yyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxxy_0[i] = g_0_yyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxxxy_0[i] * pb_z + g_0_yyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxxxz_0[i] = 2.0 * g_0_yzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxz_0[i] * pb_y + g_0_yyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxyy_0[i] = g_0_yyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxxyy_0[i] * pb_z + g_0_yyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxxyz_0[i] = 2.0 * g_0_yzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxyz_0[i] * pb_y + g_0_yyzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxzz_0[i] = 2.0 * g_0_yzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxzz_0[i] * pb_y + g_0_yyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxyyy_0[i] = g_0_yyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxyyy_0[i] * pb_z + g_0_yyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxyyz_0[i] = 2.0 * g_0_yzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyyz_0[i] * pb_y + g_0_yyzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxyzz_0[i] = 2.0 * g_0_yzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyzz_0[i] * pb_y + g_0_yyzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxzzz_0[i] = 2.0 * g_0_yzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxzzz_0[i] * pb_y + g_0_yyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyyyy_0[i] = g_0_yyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxyyyy_0[i] * pb_z + g_0_yyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxyyyz_0[i] = 2.0 * g_0_yzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyyz_0[i] * pb_y + g_0_yyzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyyzz_0[i] = 2.0 * g_0_yzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyzz_0[i] * pb_y + g_0_yyzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyzzz_0[i] = 2.0 * g_0_yzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyzzz_0[i] * pb_y + g_0_yyzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxzzzz_0[i] = 2.0 * g_0_yzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxzzzz_0[i] * pb_y + g_0_yyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyyyy_0[i] = g_0_yyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxyyyyy_0[i] * pb_z + g_0_yyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyyz_0[i] * pb_y + g_0_yyzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyzz_0[i] * pb_y + g_0_yyzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyzzz_0[i] * pb_y + g_0_yyzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyzzzz_0[i] * pb_y + g_0_yyzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxzzzzz_0[i] = 2.0 * g_0_yzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzzzzz_0[i] * pb_y + g_0_yyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyyyy_0[i] = g_0_yyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xyyyyyy_0[i] * pb_z + g_0_yyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xyyyyyz_0[i] = 2.0 * g_0_yzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyyz_0[i] * pb_y + g_0_yyzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyyzz_0[i] = 2.0 * g_0_yzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyzz_0[i] * pb_y + g_0_yyzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyzzz_0[i] = 2.0 * g_0_yzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyzzz_0[i] * pb_y + g_0_yyzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyzzzz_0[i] = 2.0 * g_0_yzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyzzzz_0[i] * pb_y + g_0_yyzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyzzzzz_0[i] = 2.0 * g_0_yzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyzzzzz_0[i] * pb_y + g_0_yyzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xzzzzzz_0[i] = 2.0 * g_0_yzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzzzzz_0[i] * pb_y + g_0_yyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyyyy_0[i] = g_0_yyy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_yyyyyyy_0[i] * pb_z + g_0_yyyz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_yyyyyyz_0[i] = 2.0 * g_0_yzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyyyz_0[i] * pb_y + g_0_yyzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyyzz_0[i] = 2.0 * g_0_yzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyyzz_0[i] * pb_y + g_0_yyzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyzzz_0[i] = 2.0 * g_0_yzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyzzz_0[i] * pb_y + g_0_yyzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyzzzz_0[i] = 2.0 * g_0_yzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyzzzz_0[i] * pb_y + g_0_yyzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyzzzzz_0[i] = 2.0 * g_0_yzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyzzzzz_0[i] * pb_y + g_0_yyzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yzzzzzz_0[i] = 2.0 * g_0_yzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yzzzzzz_0[i] * pb_y + g_0_yyzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_zzzzzzz_0[i] = 2.0 * g_0_yzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzzzzz_0[i] * pb_y + g_0_yyzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 648-684 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_yyz_0_xxxxxxy_0, g_0_yyz_0_xxxxxxy_1, g_0_yyz_0_xxxxxyy_0, g_0_yyz_0_xxxxxyy_1, g_0_yyz_0_xxxxyyy_0, g_0_yyz_0_xxxxyyy_1, g_0_yyz_0_xxxyyyy_0, g_0_yyz_0_xxxyyyy_1, g_0_yyz_0_xxyyyyy_0, g_0_yyz_0_xxyyyyy_1, g_0_yyz_0_xyyyyyy_0, g_0_yyz_0_xyyyyyy_1, g_0_yyz_0_yyyyyyy_0, g_0_yyz_0_yyyyyyy_1, g_0_yyzz_0_xxxxxxy_0, g_0_yyzz_0_xxxxxxy_1, g_0_yyzz_0_xxxxxyy_0, g_0_yyzz_0_xxxxxyy_1, g_0_yyzz_0_xxxxyyy_0, g_0_yyzz_0_xxxxyyy_1, g_0_yyzz_0_xxxyyyy_0, g_0_yyzz_0_xxxyyyy_1, g_0_yyzz_0_xxyyyyy_0, g_0_yyzz_0_xxyyyyy_1, g_0_yyzz_0_xyyyyyy_0, g_0_yyzz_0_xyyyyyy_1, g_0_yyzz_0_yyyyyyy_0, g_0_yyzz_0_yyyyyyy_1, g_0_yyzzz_0_xxxxxxx_0, g_0_yyzzz_0_xxxxxxy_0, g_0_yyzzz_0_xxxxxxz_0, g_0_yyzzz_0_xxxxxyy_0, g_0_yyzzz_0_xxxxxyz_0, g_0_yyzzz_0_xxxxxzz_0, g_0_yyzzz_0_xxxxyyy_0, g_0_yyzzz_0_xxxxyyz_0, g_0_yyzzz_0_xxxxyzz_0, g_0_yyzzz_0_xxxxzzz_0, g_0_yyzzz_0_xxxyyyy_0, g_0_yyzzz_0_xxxyyyz_0, g_0_yyzzz_0_xxxyyzz_0, g_0_yyzzz_0_xxxyzzz_0, g_0_yyzzz_0_xxxzzzz_0, g_0_yyzzz_0_xxyyyyy_0, g_0_yyzzz_0_xxyyyyz_0, g_0_yyzzz_0_xxyyyzz_0, g_0_yyzzz_0_xxyyzzz_0, g_0_yyzzz_0_xxyzzzz_0, g_0_yyzzz_0_xxzzzzz_0, g_0_yyzzz_0_xyyyyyy_0, g_0_yyzzz_0_xyyyyyz_0, g_0_yyzzz_0_xyyyyzz_0, g_0_yyzzz_0_xyyyzzz_0, g_0_yyzzz_0_xyyzzzz_0, g_0_yyzzz_0_xyzzzzz_0, g_0_yyzzz_0_xzzzzzz_0, g_0_yyzzz_0_yyyyyyy_0, g_0_yyzzz_0_yyyyyyz_0, g_0_yyzzz_0_yyyyyzz_0, g_0_yyzzz_0_yyyyzzz_0, g_0_yyzzz_0_yyyzzzz_0, g_0_yyzzz_0_yyzzzzz_0, g_0_yyzzz_0_yzzzzzz_0, g_0_yyzzz_0_zzzzzzz_0, g_0_yzzz_0_xxxxxxx_0, g_0_yzzz_0_xxxxxxx_1, g_0_yzzz_0_xxxxxxz_0, g_0_yzzz_0_xxxxxxz_1, g_0_yzzz_0_xxxxxyz_0, g_0_yzzz_0_xxxxxyz_1, g_0_yzzz_0_xxxxxz_1, g_0_yzzz_0_xxxxxzz_0, g_0_yzzz_0_xxxxxzz_1, g_0_yzzz_0_xxxxyyz_0, g_0_yzzz_0_xxxxyyz_1, g_0_yzzz_0_xxxxyz_1, g_0_yzzz_0_xxxxyzz_0, g_0_yzzz_0_xxxxyzz_1, g_0_yzzz_0_xxxxzz_1, g_0_yzzz_0_xxxxzzz_0, g_0_yzzz_0_xxxxzzz_1, g_0_yzzz_0_xxxyyyz_0, g_0_yzzz_0_xxxyyyz_1, g_0_yzzz_0_xxxyyz_1, g_0_yzzz_0_xxxyyzz_0, g_0_yzzz_0_xxxyyzz_1, g_0_yzzz_0_xxxyzz_1, g_0_yzzz_0_xxxyzzz_0, g_0_yzzz_0_xxxyzzz_1, g_0_yzzz_0_xxxzzz_1, g_0_yzzz_0_xxxzzzz_0, g_0_yzzz_0_xxxzzzz_1, g_0_yzzz_0_xxyyyyz_0, g_0_yzzz_0_xxyyyyz_1, g_0_yzzz_0_xxyyyz_1, g_0_yzzz_0_xxyyyzz_0, g_0_yzzz_0_xxyyyzz_1, g_0_yzzz_0_xxyyzz_1, g_0_yzzz_0_xxyyzzz_0, g_0_yzzz_0_xxyyzzz_1, g_0_yzzz_0_xxyzzz_1, g_0_yzzz_0_xxyzzzz_0, g_0_yzzz_0_xxyzzzz_1, g_0_yzzz_0_xxzzzz_1, g_0_yzzz_0_xxzzzzz_0, g_0_yzzz_0_xxzzzzz_1, g_0_yzzz_0_xyyyyyz_0, g_0_yzzz_0_xyyyyyz_1, g_0_yzzz_0_xyyyyz_1, g_0_yzzz_0_xyyyyzz_0, g_0_yzzz_0_xyyyyzz_1, g_0_yzzz_0_xyyyzz_1, g_0_yzzz_0_xyyyzzz_0, g_0_yzzz_0_xyyyzzz_1, g_0_yzzz_0_xyyzzz_1, g_0_yzzz_0_xyyzzzz_0, g_0_yzzz_0_xyyzzzz_1, g_0_yzzz_0_xyzzzz_1, g_0_yzzz_0_xyzzzzz_0, g_0_yzzz_0_xyzzzzz_1, g_0_yzzz_0_xzzzzz_1, g_0_yzzz_0_xzzzzzz_0, g_0_yzzz_0_xzzzzzz_1, g_0_yzzz_0_yyyyyyz_0, g_0_yzzz_0_yyyyyyz_1, g_0_yzzz_0_yyyyyz_1, g_0_yzzz_0_yyyyyzz_0, g_0_yzzz_0_yyyyyzz_1, g_0_yzzz_0_yyyyzz_1, g_0_yzzz_0_yyyyzzz_0, g_0_yzzz_0_yyyyzzz_1, g_0_yzzz_0_yyyzzz_1, g_0_yzzz_0_yyyzzzz_0, g_0_yzzz_0_yyyzzzz_1, g_0_yzzz_0_yyzzzz_1, g_0_yzzz_0_yyzzzzz_0, g_0_yzzz_0_yyzzzzz_1, g_0_yzzz_0_yzzzzz_1, g_0_yzzz_0_yzzzzzz_0, g_0_yzzz_0_yzzzzzz_1, g_0_yzzz_0_zzzzzz_1, g_0_yzzz_0_zzzzzzz_0, g_0_yzzz_0_zzzzzzz_1, g_0_zzz_0_xxxxxxx_0, g_0_zzz_0_xxxxxxx_1, g_0_zzz_0_xxxxxxz_0, g_0_zzz_0_xxxxxxz_1, g_0_zzz_0_xxxxxyz_0, g_0_zzz_0_xxxxxyz_1, g_0_zzz_0_xxxxxzz_0, g_0_zzz_0_xxxxxzz_1, g_0_zzz_0_xxxxyyz_0, g_0_zzz_0_xxxxyyz_1, g_0_zzz_0_xxxxyzz_0, g_0_zzz_0_xxxxyzz_1, g_0_zzz_0_xxxxzzz_0, g_0_zzz_0_xxxxzzz_1, g_0_zzz_0_xxxyyyz_0, g_0_zzz_0_xxxyyyz_1, g_0_zzz_0_xxxyyzz_0, g_0_zzz_0_xxxyyzz_1, g_0_zzz_0_xxxyzzz_0, g_0_zzz_0_xxxyzzz_1, g_0_zzz_0_xxxzzzz_0, g_0_zzz_0_xxxzzzz_1, g_0_zzz_0_xxyyyyz_0, g_0_zzz_0_xxyyyyz_1, g_0_zzz_0_xxyyyzz_0, g_0_zzz_0_xxyyyzz_1, g_0_zzz_0_xxyyzzz_0, g_0_zzz_0_xxyyzzz_1, g_0_zzz_0_xxyzzzz_0, g_0_zzz_0_xxyzzzz_1, g_0_zzz_0_xxzzzzz_0, g_0_zzz_0_xxzzzzz_1, g_0_zzz_0_xyyyyyz_0, g_0_zzz_0_xyyyyyz_1, g_0_zzz_0_xyyyyzz_0, g_0_zzz_0_xyyyyzz_1, g_0_zzz_0_xyyyzzz_0, g_0_zzz_0_xyyyzzz_1, g_0_zzz_0_xyyzzzz_0, g_0_zzz_0_xyyzzzz_1, g_0_zzz_0_xyzzzzz_0, g_0_zzz_0_xyzzzzz_1, g_0_zzz_0_xzzzzzz_0, g_0_zzz_0_xzzzzzz_1, g_0_zzz_0_yyyyyyz_0, g_0_zzz_0_yyyyyyz_1, g_0_zzz_0_yyyyyzz_0, g_0_zzz_0_yyyyyzz_1, g_0_zzz_0_yyyyzzz_0, g_0_zzz_0_yyyyzzz_1, g_0_zzz_0_yyyzzzz_0, g_0_zzz_0_yyyzzzz_1, g_0_zzz_0_yyzzzzz_0, g_0_zzz_0_yyzzzzz_1, g_0_zzz_0_yzzzzzz_0, g_0_zzz_0_yzzzzzz_1, g_0_zzz_0_zzzzzzz_0, g_0_zzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzz_0_xxxxxxx_0[i] = g_0_zzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxxx_0[i] * pb_y + g_0_yzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxxy_0[i] = 2.0 * g_0_yyz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxy_0[i] * pb_z + g_0_yyzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxxxz_0[i] = g_0_zzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxxz_0[i] * pb_y + g_0_yzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_yyz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxyy_0[i] * pb_z + g_0_yyzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxxyz_0[i] = g_0_zzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxyz_0[i] * pb_y + g_0_yzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxzz_0[i] = g_0_zzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxzz_0[i] * pb_y + g_0_yzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxyyy_0[i] = 2.0 * g_0_yyz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxyyy_0[i] * pb_z + g_0_yyzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxyyz_0[i] = g_0_zzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyyz_0[i] * pb_y + g_0_yzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxyzz_0[i] = g_0_zzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyzz_0[i] * pb_y + g_0_yzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxzzz_0[i] = g_0_zzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxzzz_0[i] * pb_y + g_0_yzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyyyy_0[i] = 2.0 * g_0_yyz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxyyyy_0[i] * pb_z + g_0_yyzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxyyyz_0[i] = g_0_zzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyyz_0[i] * pb_y + g_0_yzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyyzz_0[i] = g_0_zzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyzz_0[i] * pb_y + g_0_yzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyzzz_0[i] = g_0_zzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyzzz_0[i] * pb_y + g_0_yzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxzzzz_0[i] = g_0_zzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxzzzz_0[i] * pb_y + g_0_yzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxyyyyy_0[i] * pb_z + g_0_yyzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxyyyyz_0[i] = g_0_zzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyyz_0[i] * pb_y + g_0_yzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyyzz_0[i] = g_0_zzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyzz_0[i] * pb_y + g_0_yzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyzzz_0[i] = g_0_zzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyzzz_0[i] * pb_y + g_0_yzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyzzzz_0[i] = g_0_zzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyzzzz_0[i] * pb_y + g_0_yzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxzzzzz_0[i] = g_0_zzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzzzzz_0[i] * pb_y + g_0_yzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyyyy_0[i] = 2.0 * g_0_yyz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xyyyyyy_0[i] * pb_z + g_0_yyzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xyyyyyz_0[i] = g_0_zzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyyz_0[i] * pb_y + g_0_yzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyyzz_0[i] = g_0_zzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyzz_0[i] * pb_y + g_0_yzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyzzz_0[i] = g_0_zzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyzzz_0[i] * pb_y + g_0_yzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyzzzz_0[i] = g_0_zzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyzzzz_0[i] * pb_y + g_0_yzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyzzzzz_0[i] = g_0_zzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyzzzzz_0[i] * pb_y + g_0_yzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xzzzzzz_0[i] = g_0_zzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzzzzz_0[i] * pb_y + g_0_yzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyyyy_0[i] = 2.0 * g_0_yyz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_yyyyyyy_0[i] * pb_z + g_0_yyzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_yyyyyyz_0[i] = g_0_zzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyyyyz_0[i] * pb_y + g_0_yzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyyzz_0[i] = g_0_zzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyyyzz_0[i] * pb_y + g_0_yzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyzzz_0[i] = g_0_zzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyyzzz_0[i] * pb_y + g_0_yzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyzzzz_0[i] = g_0_zzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyzzzz_0[i] * pb_y + g_0_yzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyzzzzz_0[i] = g_0_zzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyzzzzz_0[i] * pb_y + g_0_yzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yzzzzzz_0[i] = g_0_zzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yzzzzzz_0[i] * pb_y + g_0_yzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_zzzzzzz_0[i] = g_0_zzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzzzzz_0[i] * pb_y + g_0_yzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 684-720 components of targeted buffer : SHSK

    auto g_0_yzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_shsk + 684);

    auto g_0_yzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_shsk + 685);

    auto g_0_yzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_shsk + 686);

    auto g_0_yzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_shsk + 687);

    auto g_0_yzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_shsk + 688);

    auto g_0_yzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_shsk + 689);

    auto g_0_yzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_shsk + 690);

    auto g_0_yzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_shsk + 691);

    auto g_0_yzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_shsk + 692);

    auto g_0_yzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_shsk + 693);

    auto g_0_yzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_shsk + 694);

    auto g_0_yzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_shsk + 695);

    auto g_0_yzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_shsk + 696);

    auto g_0_yzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_shsk + 697);

    auto g_0_yzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_shsk + 698);

    auto g_0_yzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 699);

    auto g_0_yzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 700);

    auto g_0_yzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 701);

    auto g_0_yzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 702);

    auto g_0_yzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 703);

    auto g_0_yzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 704);

    auto g_0_yzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 705);

    auto g_0_yzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 706);

    auto g_0_yzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 707);

    auto g_0_yzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 708);

    auto g_0_yzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 709);

    auto g_0_yzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 710);

    auto g_0_yzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 711);

    auto g_0_yzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_shsk + 712);

    auto g_0_yzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_shsk + 713);

    auto g_0_yzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_shsk + 714);

    auto g_0_yzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_shsk + 715);

    auto g_0_yzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_shsk + 716);

    auto g_0_yzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 717);

    auto g_0_yzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 718);

    auto g_0_yzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_shsk + 719);

    #pragma omp simd aligned(g_0_yzzzz_0_xxxxxxx_0, g_0_yzzzz_0_xxxxxxy_0, g_0_yzzzz_0_xxxxxxz_0, g_0_yzzzz_0_xxxxxyy_0, g_0_yzzzz_0_xxxxxyz_0, g_0_yzzzz_0_xxxxxzz_0, g_0_yzzzz_0_xxxxyyy_0, g_0_yzzzz_0_xxxxyyz_0, g_0_yzzzz_0_xxxxyzz_0, g_0_yzzzz_0_xxxxzzz_0, g_0_yzzzz_0_xxxyyyy_0, g_0_yzzzz_0_xxxyyyz_0, g_0_yzzzz_0_xxxyyzz_0, g_0_yzzzz_0_xxxyzzz_0, g_0_yzzzz_0_xxxzzzz_0, g_0_yzzzz_0_xxyyyyy_0, g_0_yzzzz_0_xxyyyyz_0, g_0_yzzzz_0_xxyyyzz_0, g_0_yzzzz_0_xxyyzzz_0, g_0_yzzzz_0_xxyzzzz_0, g_0_yzzzz_0_xxzzzzz_0, g_0_yzzzz_0_xyyyyyy_0, g_0_yzzzz_0_xyyyyyz_0, g_0_yzzzz_0_xyyyyzz_0, g_0_yzzzz_0_xyyyzzz_0, g_0_yzzzz_0_xyyzzzz_0, g_0_yzzzz_0_xyzzzzz_0, g_0_yzzzz_0_xzzzzzz_0, g_0_yzzzz_0_yyyyyyy_0, g_0_yzzzz_0_yyyyyyz_0, g_0_yzzzz_0_yyyyyzz_0, g_0_yzzzz_0_yyyyzzz_0, g_0_yzzzz_0_yyyzzzz_0, g_0_yzzzz_0_yyzzzzz_0, g_0_yzzzz_0_yzzzzzz_0, g_0_yzzzz_0_zzzzzzz_0, g_0_zzzz_0_xxxxxx_1, g_0_zzzz_0_xxxxxxx_0, g_0_zzzz_0_xxxxxxx_1, g_0_zzzz_0_xxxxxxy_0, g_0_zzzz_0_xxxxxxy_1, g_0_zzzz_0_xxxxxxz_0, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxy_1, g_0_zzzz_0_xxxxxyy_0, g_0_zzzz_0_xxxxxyy_1, g_0_zzzz_0_xxxxxyz_0, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxz_1, g_0_zzzz_0_xxxxxzz_0, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxyy_1, g_0_zzzz_0_xxxxyyy_0, g_0_zzzz_0_xxxxyyy_1, g_0_zzzz_0_xxxxyyz_0, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyz_1, g_0_zzzz_0_xxxxyzz_0, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxzz_1, g_0_zzzz_0_xxxxzzz_0, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxyyy_1, g_0_zzzz_0_xxxyyyy_0, g_0_zzzz_0_xxxyyyy_1, g_0_zzzz_0_xxxyyyz_0, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyz_1, g_0_zzzz_0_xxxyyzz_0, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyzz_1, g_0_zzzz_0_xxxyzzz_0, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxzzz_1, g_0_zzzz_0_xxxzzzz_0, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxyyyy_1, g_0_zzzz_0_xxyyyyy_0, g_0_zzzz_0_xxyyyyy_1, g_0_zzzz_0_xxyyyyz_0, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyz_1, g_0_zzzz_0_xxyyyzz_0, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyzz_1, g_0_zzzz_0_xxyyzzz_0, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyzzz_1, g_0_zzzz_0_xxyzzzz_0, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxzzzz_1, g_0_zzzz_0_xxzzzzz_0, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xyyyyy_1, g_0_zzzz_0_xyyyyyy_0, g_0_zzzz_0_xyyyyyy_1, g_0_zzzz_0_xyyyyyz_0, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyz_1, g_0_zzzz_0_xyyyyzz_0, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyzz_1, g_0_zzzz_0_xyyyzzz_0, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyzzz_1, g_0_zzzz_0_xyyzzzz_0, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyzzzz_1, g_0_zzzz_0_xyzzzzz_0, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xzzzzz_1, g_0_zzzz_0_xzzzzzz_0, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_yyyyyy_1, g_0_zzzz_0_yyyyyyy_0, g_0_zzzz_0_yyyyyyy_1, g_0_zzzz_0_yyyyyyz_0, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyz_1, g_0_zzzz_0_yyyyyzz_0, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyzz_1, g_0_zzzz_0_yyyyzzz_0, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyzzz_1, g_0_zzzz_0_yyyzzzz_0, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyzzzz_1, g_0_zzzz_0_yyzzzzz_0, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yzzzzz_1, g_0_zzzz_0_yzzzzzz_0, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_zzzzzz_1, g_0_zzzz_0_zzzzzzz_0, g_0_zzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzz_0_xxxxxxx_0[i] = g_0_zzzz_0_xxxxxxx_0[i] * pb_y + g_0_zzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxxy_0[i] = g_0_zzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxy_0[i] * pb_y + g_0_zzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxxz_0[i] = g_0_zzzz_0_xxxxxxz_0[i] * pb_y + g_0_zzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxyy_0[i] = 2.0 * g_0_zzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyy_0[i] * pb_y + g_0_zzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxyz_0[i] = g_0_zzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyz_0[i] * pb_y + g_0_zzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxzz_0[i] = g_0_zzzz_0_xxxxxzz_0[i] * pb_y + g_0_zzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_zzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyy_0[i] * pb_y + g_0_zzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_zzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyz_0[i] * pb_y + g_0_zzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyzz_0[i] = g_0_zzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyzz_0[i] * pb_y + g_0_zzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxzzz_0[i] = g_0_zzzz_0_xxxxzzz_0[i] * pb_y + g_0_zzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_zzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyy_0[i] * pb_y + g_0_zzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyz_0[i] * pb_y + g_0_zzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_zzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyzz_0[i] * pb_y + g_0_zzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyzzz_0[i] = g_0_zzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzzz_0[i] * pb_y + g_0_zzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxzzzz_0[i] = g_0_zzzz_0_xxxzzzz_0[i] * pb_y + g_0_zzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_zzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyy_0[i] * pb_y + g_0_zzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyyz_0[i] = 4.0 * g_0_zzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyz_0[i] * pb_y + g_0_zzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_zzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyzz_0[i] * pb_y + g_0_zzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzzz_0[i] * pb_y + g_0_zzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyzzzz_0[i] = g_0_zzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzzz_0[i] * pb_y + g_0_zzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxzzzzz_0[i] = g_0_zzzz_0_xxzzzzz_0[i] * pb_y + g_0_zzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyyy_0[i] = 6.0 * g_0_zzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyy_0[i] * pb_y + g_0_zzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyyz_0[i] = 5.0 * g_0_zzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyz_0[i] * pb_y + g_0_zzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyzz_0[i] = 4.0 * g_0_zzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyzz_0[i] * pb_y + g_0_zzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_zzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzzz_0[i] * pb_y + g_0_zzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_zzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzzz_0[i] * pb_y + g_0_zzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyzzzzz_0[i] = g_0_zzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzzz_0[i] * pb_y + g_0_zzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xzzzzzz_0[i] = g_0_zzzz_0_xzzzzzz_0[i] * pb_y + g_0_zzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyyy_0[i] = 7.0 * g_0_zzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyyy_0[i] * pb_y + g_0_zzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyyz_0[i] = 6.0 * g_0_zzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyyz_0[i] * pb_y + g_0_zzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyzz_0[i] = 5.0 * g_0_zzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyzz_0[i] * pb_y + g_0_zzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyzzz_0[i] = 4.0 * g_0_zzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyzzz_0[i] * pb_y + g_0_zzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_zzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyzzzz_0[i] * pb_y + g_0_zzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_zzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzzzzz_0[i] * pb_y + g_0_zzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yzzzzzz_0[i] = g_0_zzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzzzzz_0[i] * pb_y + g_0_zzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_zzzzzzz_0[i] = g_0_zzzz_0_zzzzzzz_0[i] * pb_y + g_0_zzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 720-756 components of targeted buffer : SHSK

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

    #pragma omp simd aligned(g_0_zzz_0_xxxxxxx_0, g_0_zzz_0_xxxxxxx_1, g_0_zzz_0_xxxxxxy_0, g_0_zzz_0_xxxxxxy_1, g_0_zzz_0_xxxxxxz_0, g_0_zzz_0_xxxxxxz_1, g_0_zzz_0_xxxxxyy_0, g_0_zzz_0_xxxxxyy_1, g_0_zzz_0_xxxxxyz_0, g_0_zzz_0_xxxxxyz_1, g_0_zzz_0_xxxxxzz_0, g_0_zzz_0_xxxxxzz_1, g_0_zzz_0_xxxxyyy_0, g_0_zzz_0_xxxxyyy_1, g_0_zzz_0_xxxxyyz_0, g_0_zzz_0_xxxxyyz_1, g_0_zzz_0_xxxxyzz_0, g_0_zzz_0_xxxxyzz_1, g_0_zzz_0_xxxxzzz_0, g_0_zzz_0_xxxxzzz_1, g_0_zzz_0_xxxyyyy_0, g_0_zzz_0_xxxyyyy_1, g_0_zzz_0_xxxyyyz_0, g_0_zzz_0_xxxyyyz_1, g_0_zzz_0_xxxyyzz_0, g_0_zzz_0_xxxyyzz_1, g_0_zzz_0_xxxyzzz_0, g_0_zzz_0_xxxyzzz_1, g_0_zzz_0_xxxzzzz_0, g_0_zzz_0_xxxzzzz_1, g_0_zzz_0_xxyyyyy_0, g_0_zzz_0_xxyyyyy_1, g_0_zzz_0_xxyyyyz_0, g_0_zzz_0_xxyyyyz_1, g_0_zzz_0_xxyyyzz_0, g_0_zzz_0_xxyyyzz_1, g_0_zzz_0_xxyyzzz_0, g_0_zzz_0_xxyyzzz_1, g_0_zzz_0_xxyzzzz_0, g_0_zzz_0_xxyzzzz_1, g_0_zzz_0_xxzzzzz_0, g_0_zzz_0_xxzzzzz_1, g_0_zzz_0_xyyyyyy_0, g_0_zzz_0_xyyyyyy_1, g_0_zzz_0_xyyyyyz_0, g_0_zzz_0_xyyyyyz_1, g_0_zzz_0_xyyyyzz_0, g_0_zzz_0_xyyyyzz_1, g_0_zzz_0_xyyyzzz_0, g_0_zzz_0_xyyyzzz_1, g_0_zzz_0_xyyzzzz_0, g_0_zzz_0_xyyzzzz_1, g_0_zzz_0_xyzzzzz_0, g_0_zzz_0_xyzzzzz_1, g_0_zzz_0_xzzzzzz_0, g_0_zzz_0_xzzzzzz_1, g_0_zzz_0_yyyyyyy_0, g_0_zzz_0_yyyyyyy_1, g_0_zzz_0_yyyyyyz_0, g_0_zzz_0_yyyyyyz_1, g_0_zzz_0_yyyyyzz_0, g_0_zzz_0_yyyyyzz_1, g_0_zzz_0_yyyyzzz_0, g_0_zzz_0_yyyyzzz_1, g_0_zzz_0_yyyzzzz_0, g_0_zzz_0_yyyzzzz_1, g_0_zzz_0_yyzzzzz_0, g_0_zzz_0_yyzzzzz_1, g_0_zzz_0_yzzzzzz_0, g_0_zzz_0_yzzzzzz_1, g_0_zzz_0_zzzzzzz_0, g_0_zzz_0_zzzzzzz_1, g_0_zzzz_0_xxxxxx_1, g_0_zzzz_0_xxxxxxx_0, g_0_zzzz_0_xxxxxxx_1, g_0_zzzz_0_xxxxxxy_0, g_0_zzzz_0_xxxxxxy_1, g_0_zzzz_0_xxxxxxz_0, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxy_1, g_0_zzzz_0_xxxxxyy_0, g_0_zzzz_0_xxxxxyy_1, g_0_zzzz_0_xxxxxyz_0, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxz_1, g_0_zzzz_0_xxxxxzz_0, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxyy_1, g_0_zzzz_0_xxxxyyy_0, g_0_zzzz_0_xxxxyyy_1, g_0_zzzz_0_xxxxyyz_0, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyz_1, g_0_zzzz_0_xxxxyzz_0, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxzz_1, g_0_zzzz_0_xxxxzzz_0, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxyyy_1, g_0_zzzz_0_xxxyyyy_0, g_0_zzzz_0_xxxyyyy_1, g_0_zzzz_0_xxxyyyz_0, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyz_1, g_0_zzzz_0_xxxyyzz_0, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyzz_1, g_0_zzzz_0_xxxyzzz_0, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxzzz_1, g_0_zzzz_0_xxxzzzz_0, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxyyyy_1, g_0_zzzz_0_xxyyyyy_0, g_0_zzzz_0_xxyyyyy_1, g_0_zzzz_0_xxyyyyz_0, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyz_1, g_0_zzzz_0_xxyyyzz_0, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyzz_1, g_0_zzzz_0_xxyyzzz_0, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyzzz_1, g_0_zzzz_0_xxyzzzz_0, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxzzzz_1, g_0_zzzz_0_xxzzzzz_0, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xyyyyy_1, g_0_zzzz_0_xyyyyyy_0, g_0_zzzz_0_xyyyyyy_1, g_0_zzzz_0_xyyyyyz_0, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyz_1, g_0_zzzz_0_xyyyyzz_0, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyzz_1, g_0_zzzz_0_xyyyzzz_0, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyzzz_1, g_0_zzzz_0_xyyzzzz_0, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyzzzz_1, g_0_zzzz_0_xyzzzzz_0, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xzzzzz_1, g_0_zzzz_0_xzzzzzz_0, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_yyyyyy_1, g_0_zzzz_0_yyyyyyy_0, g_0_zzzz_0_yyyyyyy_1, g_0_zzzz_0_yyyyyyz_0, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyz_1, g_0_zzzz_0_yyyyyzz_0, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyzz_1, g_0_zzzz_0_yyyyzzz_0, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyzzz_1, g_0_zzzz_0_yyyzzzz_0, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyzzzz_1, g_0_zzzz_0_yyzzzzz_0, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yzzzzz_1, g_0_zzzz_0_yzzzzzz_0, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_zzzzzz_1, g_0_zzzz_0_zzzzzzz_0, g_0_zzzz_0_zzzzzzz_1, g_0_zzzzz_0_xxxxxxx_0, g_0_zzzzz_0_xxxxxxy_0, g_0_zzzzz_0_xxxxxxz_0, g_0_zzzzz_0_xxxxxyy_0, g_0_zzzzz_0_xxxxxyz_0, g_0_zzzzz_0_xxxxxzz_0, g_0_zzzzz_0_xxxxyyy_0, g_0_zzzzz_0_xxxxyyz_0, g_0_zzzzz_0_xxxxyzz_0, g_0_zzzzz_0_xxxxzzz_0, g_0_zzzzz_0_xxxyyyy_0, g_0_zzzzz_0_xxxyyyz_0, g_0_zzzzz_0_xxxyyzz_0, g_0_zzzzz_0_xxxyzzz_0, g_0_zzzzz_0_xxxzzzz_0, g_0_zzzzz_0_xxyyyyy_0, g_0_zzzzz_0_xxyyyyz_0, g_0_zzzzz_0_xxyyyzz_0, g_0_zzzzz_0_xxyyzzz_0, g_0_zzzzz_0_xxyzzzz_0, g_0_zzzzz_0_xxzzzzz_0, g_0_zzzzz_0_xyyyyyy_0, g_0_zzzzz_0_xyyyyyz_0, g_0_zzzzz_0_xyyyyzz_0, g_0_zzzzz_0_xyyyzzz_0, g_0_zzzzz_0_xyyzzzz_0, g_0_zzzzz_0_xyzzzzz_0, g_0_zzzzz_0_xzzzzzz_0, g_0_zzzzz_0_yyyyyyy_0, g_0_zzzzz_0_yyyyyyz_0, g_0_zzzzz_0_yyyyyzz_0, g_0_zzzzz_0_yyyyzzz_0, g_0_zzzzz_0_yyyzzzz_0, g_0_zzzzz_0_yyzzzzz_0, g_0_zzzzz_0_yzzzzzz_0, g_0_zzzzz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzz_0_xxxxxxx_0[i] = 4.0 * g_0_zzz_0_xxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxxx_0[i] * pb_z + g_0_zzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxxy_0[i] = 4.0 * g_0_zzz_0_xxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxxy_0[i] * pb_z + g_0_zzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxxz_0[i] = 4.0 * g_0_zzz_0_xxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxz_0[i] * pb_z + g_0_zzzz_0_xxxxxxz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxyy_0[i] = 4.0 * g_0_zzz_0_xxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxyy_0[i] * pb_z + g_0_zzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxyz_0[i] = 4.0 * g_0_zzz_0_xxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyz_0[i] * pb_z + g_0_zzzz_0_xxxxxyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxzz_0[i] = 4.0 * g_0_zzz_0_xxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxzz_0[i] * pb_z + g_0_zzzz_0_xxxxxzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_zzz_0_xxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxyyy_0[i] * pb_z + g_0_zzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_zzz_0_xxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyz_0[i] * pb_z + g_0_zzzz_0_xxxxyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_zzz_0_xxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyzz_0[i] * pb_z + g_0_zzzz_0_xxxxyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxzzz_0[i] = 4.0 * g_0_zzz_0_xxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxzzz_0[i] * pb_z + g_0_zzzz_0_xxxxzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_zzz_0_xxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxyyyy_0[i] * pb_z + g_0_zzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyyz_0[i] = 4.0 * g_0_zzz_0_xxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyz_0[i] * pb_z + g_0_zzzz_0_xxxyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyzz_0[i] = 4.0 * g_0_zzz_0_xxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyzz_0[i] * pb_z + g_0_zzzz_0_xxxyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyzzz_0[i] = 4.0 * g_0_zzz_0_xxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzzz_0[i] * pb_z + g_0_zzzz_0_xxxyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxzzzz_0[i] = 4.0 * g_0_zzz_0_xxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxzzzz_0[i] * pb_z + g_0_zzzz_0_xxxzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyyy_0[i] = 4.0 * g_0_zzz_0_xxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxyyyyy_0[i] * pb_z + g_0_zzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyyz_0[i] = 4.0 * g_0_zzz_0_xxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyz_0[i] * pb_z + g_0_zzzz_0_xxyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyzz_0[i] = 4.0 * g_0_zzz_0_xxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyzz_0[i] * pb_z + g_0_zzzz_0_xxyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyzzz_0[i] = 4.0 * g_0_zzz_0_xxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzzz_0[i] * pb_z + g_0_zzzz_0_xxyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyzzzz_0[i] = 4.0 * g_0_zzz_0_xxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzzz_0[i] * pb_z + g_0_zzzz_0_xxyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxzzzzz_0[i] = 4.0 * g_0_zzz_0_xxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzzzzz_0[i] * pb_z + g_0_zzzz_0_xxzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyyy_0[i] = 4.0 * g_0_zzz_0_xyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xyyyyyy_0[i] * pb_z + g_0_zzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyyz_0[i] = 4.0 * g_0_zzz_0_xyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyz_0[i] * pb_z + g_0_zzzz_0_xyyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyzz_0[i] = 4.0 * g_0_zzz_0_xyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyzz_0[i] * pb_z + g_0_zzzz_0_xyyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyzzz_0[i] = 4.0 * g_0_zzz_0_xyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzzz_0[i] * pb_z + g_0_zzzz_0_xyyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyzzzz_0[i] = 4.0 * g_0_zzz_0_xyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzzz_0[i] * pb_z + g_0_zzzz_0_xyyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyzzzzz_0[i] = 4.0 * g_0_zzz_0_xyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzzz_0[i] * pb_z + g_0_zzzz_0_xyzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xzzzzzz_0[i] = 4.0 * g_0_zzz_0_xzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzzzzz_0[i] * pb_z + g_0_zzzz_0_xzzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyyy_0[i] = 4.0 * g_0_zzz_0_yyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_yyyyyyy_0[i] * pb_z + g_0_zzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyyz_0[i] = 4.0 * g_0_zzz_0_yyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyyz_0[i] * pb_z + g_0_zzzz_0_yyyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyzz_0[i] = 4.0 * g_0_zzz_0_yyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyzz_0[i] * pb_z + g_0_zzzz_0_yyyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyzzz_0[i] = 4.0 * g_0_zzz_0_yyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyzzz_0[i] * pb_z + g_0_zzzz_0_yyyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyzzzz_0[i] = 4.0 * g_0_zzz_0_yyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyzzzz_0[i] * pb_z + g_0_zzzz_0_yyyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyzzzzz_0[i] = 4.0 * g_0_zzz_0_yyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzzzzz_0[i] * pb_z + g_0_zzzz_0_yyzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yzzzzzz_0[i] = 4.0 * g_0_zzz_0_yzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzzzzz_0[i] * pb_z + g_0_zzzz_0_yzzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_zzzzzzz_0[i] = 4.0 * g_0_zzz_0_zzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_zzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_zzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_zzzzzzz_0[i] * pb_z + g_0_zzzz_0_zzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

