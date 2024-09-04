#include "ElectronRepulsionPrimRecSHSH.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_shsh(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_shsh,
                                  size_t idx_eri_0_sfsh,
                                  size_t idx_eri_1_sfsh,
                                  size_t idx_eri_1_sgsg,
                                  size_t idx_eri_0_sgsh,
                                  size_t idx_eri_1_sgsh,
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

    /// Set up components of auxilary buffer : SFSH

    auto g_0_xxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh);

    auto g_0_xxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 1);

    auto g_0_xxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 2);

    auto g_0_xxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 3);

    auto g_0_xxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 4);

    auto g_0_xxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 5);

    auto g_0_xxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 6);

    auto g_0_xxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 7);

    auto g_0_xxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 8);

    auto g_0_xxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 9);

    auto g_0_xxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 10);

    auto g_0_xxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 11);

    auto g_0_xxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 12);

    auto g_0_xxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 13);

    auto g_0_xxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 14);

    auto g_0_xxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 15);

    auto g_0_xxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 16);

    auto g_0_xxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 17);

    auto g_0_xxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 18);

    auto g_0_xxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 19);

    auto g_0_xxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 20);

    auto g_0_xxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 21);

    auto g_0_xxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 23);

    auto g_0_xxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 26);

    auto g_0_xxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 30);

    auto g_0_xxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 35);

    auto g_0_xxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 42);

    auto g_0_xxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 43);

    auto g_0_xxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 45);

    auto g_0_xxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 48);

    auto g_0_xxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 52);

    auto g_0_xyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 64);

    auto g_0_xyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 66);

    auto g_0_xyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 67);

    auto g_0_xyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 69);

    auto g_0_xyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 70);

    auto g_0_xyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 71);

    auto g_0_xyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 73);

    auto g_0_xyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 74);

    auto g_0_xyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 75);

    auto g_0_xyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 76);

    auto g_0_xyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 78);

    auto g_0_xyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 79);

    auto g_0_xyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 80);

    auto g_0_xyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 81);

    auto g_0_xyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 82);

    auto g_0_xyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 83);

    auto g_0_xzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 107);

    auto g_0_xzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 109);

    auto g_0_xzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 110);

    auto g_0_xzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 112);

    auto g_0_xzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 113);

    auto g_0_xzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 114);

    auto g_0_xzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 116);

    auto g_0_xzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 117);

    auto g_0_xzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 118);

    auto g_0_xzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 119);

    auto g_0_xzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 120);

    auto g_0_xzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 121);

    auto g_0_xzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 122);

    auto g_0_xzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 123);

    auto g_0_xzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 124);

    auto g_0_xzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 125);

    auto g_0_yyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 126);

    auto g_0_yyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 127);

    auto g_0_yyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 128);

    auto g_0_yyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 129);

    auto g_0_yyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 130);

    auto g_0_yyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 131);

    auto g_0_yyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 132);

    auto g_0_yyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 133);

    auto g_0_yyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 134);

    auto g_0_yyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 135);

    auto g_0_yyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 136);

    auto g_0_yyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 137);

    auto g_0_yyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 138);

    auto g_0_yyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 139);

    auto g_0_yyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 140);

    auto g_0_yyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 141);

    auto g_0_yyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 142);

    auto g_0_yyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 143);

    auto g_0_yyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 144);

    auto g_0_yyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 145);

    auto g_0_yyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 146);

    auto g_0_yyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 148);

    auto g_0_yyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 150);

    auto g_0_yyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 153);

    auto g_0_yyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 157);

    auto g_0_yyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 162);

    auto g_0_yzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 168);

    auto g_0_yzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 170);

    auto g_0_yzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 172);

    auto g_0_yzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 173);

    auto g_0_yzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 175);

    auto g_0_yzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 176);

    auto g_0_yzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 177);

    auto g_0_yzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 179);

    auto g_0_yzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 180);

    auto g_0_yzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 181);

    auto g_0_yzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 182);

    auto g_0_yzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 184);

    auto g_0_yzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 185);

    auto g_0_yzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 186);

    auto g_0_yzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 187);

    auto g_0_yzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 188);

    auto g_0_zzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 189);

    auto g_0_zzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 190);

    auto g_0_zzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 191);

    auto g_0_zzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 192);

    auto g_0_zzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 193);

    auto g_0_zzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 194);

    auto g_0_zzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 195);

    auto g_0_zzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 196);

    auto g_0_zzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 197);

    auto g_0_zzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 198);

    auto g_0_zzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 199);

    auto g_0_zzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 200);

    auto g_0_zzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 201);

    auto g_0_zzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 202);

    auto g_0_zzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 203);

    auto g_0_zzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 204);

    auto g_0_zzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 205);

    auto g_0_zzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 206);

    auto g_0_zzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 207);

    auto g_0_zzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 208);

    auto g_0_zzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 209);

    /// Set up components of auxilary buffer : SFSH

    auto g_0_xxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh);

    auto g_0_xxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 1);

    auto g_0_xxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 2);

    auto g_0_xxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 3);

    auto g_0_xxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 4);

    auto g_0_xxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 5);

    auto g_0_xxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 6);

    auto g_0_xxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 7);

    auto g_0_xxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 8);

    auto g_0_xxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 9);

    auto g_0_xxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 10);

    auto g_0_xxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 11);

    auto g_0_xxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 12);

    auto g_0_xxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 13);

    auto g_0_xxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 14);

    auto g_0_xxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 15);

    auto g_0_xxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 16);

    auto g_0_xxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 17);

    auto g_0_xxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 18);

    auto g_0_xxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 19);

    auto g_0_xxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 20);

    auto g_0_xxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 21);

    auto g_0_xxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 23);

    auto g_0_xxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 26);

    auto g_0_xxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 30);

    auto g_0_xxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 35);

    auto g_0_xxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 42);

    auto g_0_xxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 43);

    auto g_0_xxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 45);

    auto g_0_xxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 48);

    auto g_0_xxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 52);

    auto g_0_xyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 64);

    auto g_0_xyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 66);

    auto g_0_xyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 67);

    auto g_0_xyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 69);

    auto g_0_xyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 70);

    auto g_0_xyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 71);

    auto g_0_xyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 73);

    auto g_0_xyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 74);

    auto g_0_xyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 75);

    auto g_0_xyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 76);

    auto g_0_xyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 78);

    auto g_0_xyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 79);

    auto g_0_xyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 80);

    auto g_0_xyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 81);

    auto g_0_xyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 82);

    auto g_0_xyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 83);

    auto g_0_xzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 107);

    auto g_0_xzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 109);

    auto g_0_xzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 110);

    auto g_0_xzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 112);

    auto g_0_xzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 113);

    auto g_0_xzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 114);

    auto g_0_xzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 116);

    auto g_0_xzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 117);

    auto g_0_xzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 118);

    auto g_0_xzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 119);

    auto g_0_xzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 120);

    auto g_0_xzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 121);

    auto g_0_xzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 122);

    auto g_0_xzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 123);

    auto g_0_xzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 124);

    auto g_0_xzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 125);

    auto g_0_yyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 126);

    auto g_0_yyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 127);

    auto g_0_yyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 128);

    auto g_0_yyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 129);

    auto g_0_yyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 130);

    auto g_0_yyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 131);

    auto g_0_yyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 132);

    auto g_0_yyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 133);

    auto g_0_yyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 134);

    auto g_0_yyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 135);

    auto g_0_yyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 136);

    auto g_0_yyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 137);

    auto g_0_yyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 138);

    auto g_0_yyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 139);

    auto g_0_yyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 140);

    auto g_0_yyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 141);

    auto g_0_yyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 142);

    auto g_0_yyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 143);

    auto g_0_yyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 144);

    auto g_0_yyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 145);

    auto g_0_yyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 146);

    auto g_0_yyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 148);

    auto g_0_yyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 150);

    auto g_0_yyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 153);

    auto g_0_yyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 157);

    auto g_0_yyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 162);

    auto g_0_yzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 168);

    auto g_0_yzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 170);

    auto g_0_yzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 172);

    auto g_0_yzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 173);

    auto g_0_yzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 175);

    auto g_0_yzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 176);

    auto g_0_yzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 177);

    auto g_0_yzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 179);

    auto g_0_yzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 180);

    auto g_0_yzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 181);

    auto g_0_yzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 182);

    auto g_0_yzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 184);

    auto g_0_yzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 185);

    auto g_0_yzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 186);

    auto g_0_yzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 187);

    auto g_0_yzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 188);

    auto g_0_zzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 189);

    auto g_0_zzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 190);

    auto g_0_zzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 191);

    auto g_0_zzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 192);

    auto g_0_zzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 193);

    auto g_0_zzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 194);

    auto g_0_zzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 195);

    auto g_0_zzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 196);

    auto g_0_zzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 197);

    auto g_0_zzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 198);

    auto g_0_zzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 199);

    auto g_0_zzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 200);

    auto g_0_zzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 201);

    auto g_0_zzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 202);

    auto g_0_zzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 203);

    auto g_0_zzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 204);

    auto g_0_zzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 205);

    auto g_0_zzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 206);

    auto g_0_zzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 207);

    auto g_0_zzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 208);

    auto g_0_zzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 209);

    /// Set up components of auxilary buffer : SGSG

    auto g_0_xxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg);

    auto g_0_xxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 1);

    auto g_0_xxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 2);

    auto g_0_xxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 3);

    auto g_0_xxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 4);

    auto g_0_xxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 5);

    auto g_0_xxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 6);

    auto g_0_xxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 7);

    auto g_0_xxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 8);

    auto g_0_xxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 9);

    auto g_0_xxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 10);

    auto g_0_xxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 11);

    auto g_0_xxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 12);

    auto g_0_xxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 13);

    auto g_0_xxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 14);

    auto g_0_xxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 32);

    auto g_0_xxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 34);

    auto g_0_xxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 35);

    auto g_0_xxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 37);

    auto g_0_xxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 38);

    auto g_0_xxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 39);

    auto g_0_xxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 41);

    auto g_0_xxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 42);

    auto g_0_xxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 43);

    auto g_0_xxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 44);

    auto g_0_xxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 45);

    auto g_0_xxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 46);

    auto g_0_xxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 47);

    auto g_0_xxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 48);

    auto g_0_xxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 49);

    auto g_0_xxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 50);

    auto g_0_xxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 51);

    auto g_0_xxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 52);

    auto g_0_xxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 53);

    auto g_0_xxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 54);

    auto g_0_xxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 55);

    auto g_0_xxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 56);

    auto g_0_xxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 57);

    auto g_0_xxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 58);

    auto g_0_xxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 59);

    auto g_0_xxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 75);

    auto g_0_xxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 76);

    auto g_0_xxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 77);

    auto g_0_xxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 78);

    auto g_0_xxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 79);

    auto g_0_xxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 80);

    auto g_0_xxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 81);

    auto g_0_xxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 82);

    auto g_0_xxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 83);

    auto g_0_xxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 84);

    auto g_0_xxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 85);

    auto g_0_xxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 86);

    auto g_0_xxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 87);

    auto g_0_xxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 88);

    auto g_0_xxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 89);

    auto g_0_xyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 91);

    auto g_0_xyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 93);

    auto g_0_xyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 94);

    auto g_0_xyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 96);

    auto g_0_xyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 97);

    auto g_0_xyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 98);

    auto g_0_xyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 100);

    auto g_0_xyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 101);

    auto g_0_xyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 102);

    auto g_0_xyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 103);

    auto g_0_xzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 137);

    auto g_0_xzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 139);

    auto g_0_xzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 140);

    auto g_0_xzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 142);

    auto g_0_xzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 143);

    auto g_0_xzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 144);

    auto g_0_xzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 146);

    auto g_0_xzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 147);

    auto g_0_xzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 148);

    auto g_0_xzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 149);

    auto g_0_yyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 150);

    auto g_0_yyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 151);

    auto g_0_yyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 152);

    auto g_0_yyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 153);

    auto g_0_yyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 154);

    auto g_0_yyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 155);

    auto g_0_yyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 156);

    auto g_0_yyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 157);

    auto g_0_yyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 158);

    auto g_0_yyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 159);

    auto g_0_yyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 160);

    auto g_0_yyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 161);

    auto g_0_yyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 162);

    auto g_0_yyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 163);

    auto g_0_yyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 164);

    auto g_0_yyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 167);

    auto g_0_yyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 169);

    auto g_0_yyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 170);

    auto g_0_yyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 172);

    auto g_0_yyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 173);

    auto g_0_yyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 174);

    auto g_0_yyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 176);

    auto g_0_yyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 177);

    auto g_0_yyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 178);

    auto g_0_yyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 179);

    auto g_0_yyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 180);

    auto g_0_yyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 181);

    auto g_0_yyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 182);

    auto g_0_yyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 183);

    auto g_0_yyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 184);

    auto g_0_yyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 185);

    auto g_0_yyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 186);

    auto g_0_yyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 187);

    auto g_0_yyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 188);

    auto g_0_yyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 189);

    auto g_0_yyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 190);

    auto g_0_yyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 191);

    auto g_0_yyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 192);

    auto g_0_yyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 193);

    auto g_0_yyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 194);

    auto g_0_yzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 196);

    auto g_0_yzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 197);

    auto g_0_yzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 198);

    auto g_0_yzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 199);

    auto g_0_yzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 200);

    auto g_0_yzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 201);

    auto g_0_yzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 202);

    auto g_0_yzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 203);

    auto g_0_yzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 204);

    auto g_0_yzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 205);

    auto g_0_yzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 206);

    auto g_0_yzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 207);

    auto g_0_yzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 208);

    auto g_0_yzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 209);

    auto g_0_zzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 210);

    auto g_0_zzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 211);

    auto g_0_zzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 212);

    auto g_0_zzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 213);

    auto g_0_zzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 214);

    auto g_0_zzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 215);

    auto g_0_zzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 216);

    auto g_0_zzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 217);

    auto g_0_zzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 218);

    auto g_0_zzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 219);

    auto g_0_zzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 220);

    auto g_0_zzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 221);

    auto g_0_zzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 222);

    auto g_0_zzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 223);

    auto g_0_zzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 224);

    /// Set up components of auxilary buffer : SGSH

    auto g_0_xxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh);

    auto g_0_xxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 1);

    auto g_0_xxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 2);

    auto g_0_xxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 3);

    auto g_0_xxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 4);

    auto g_0_xxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 5);

    auto g_0_xxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 6);

    auto g_0_xxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 7);

    auto g_0_xxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 8);

    auto g_0_xxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 9);

    auto g_0_xxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 10);

    auto g_0_xxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 11);

    auto g_0_xxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 12);

    auto g_0_xxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 13);

    auto g_0_xxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 14);

    auto g_0_xxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 15);

    auto g_0_xxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 16);

    auto g_0_xxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 17);

    auto g_0_xxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 18);

    auto g_0_xxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 19);

    auto g_0_xxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 20);

    auto g_0_xxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 21);

    auto g_0_xxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 22);

    auto g_0_xxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 23);

    auto g_0_xxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 24);

    auto g_0_xxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 26);

    auto g_0_xxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 27);

    auto g_0_xxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 30);

    auto g_0_xxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 31);

    auto g_0_xxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 35);

    auto g_0_xxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 36);

    auto g_0_xxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 42);

    auto g_0_xxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 43);

    auto g_0_xxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 44);

    auto g_0_xxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 45);

    auto g_0_xxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 46);

    auto g_0_xxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 47);

    auto g_0_xxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 48);

    auto g_0_xxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 49);

    auto g_0_xxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 50);

    auto g_0_xxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 51);

    auto g_0_xxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 52);

    auto g_0_xxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 53);

    auto g_0_xxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 54);

    auto g_0_xxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 55);

    auto g_0_xxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 56);

    auto g_0_xxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 58);

    auto g_0_xxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 59);

    auto g_0_xxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 60);

    auto g_0_xxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 61);

    auto g_0_xxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 62);

    auto g_0_xxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 63);

    auto g_0_xxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 64);

    auto g_0_xxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 65);

    auto g_0_xxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 66);

    auto g_0_xxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 67);

    auto g_0_xxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 68);

    auto g_0_xxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 69);

    auto g_0_xxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 70);

    auto g_0_xxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 71);

    auto g_0_xxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 72);

    auto g_0_xxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 73);

    auto g_0_xxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 74);

    auto g_0_xxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 75);

    auto g_0_xxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 76);

    auto g_0_xxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 77);

    auto g_0_xxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 78);

    auto g_0_xxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 79);

    auto g_0_xxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 80);

    auto g_0_xxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 81);

    auto g_0_xxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 82);

    auto g_0_xxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 83);

    auto g_0_xxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 105);

    auto g_0_xxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 106);

    auto g_0_xxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 107);

    auto g_0_xxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 108);

    auto g_0_xxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 109);

    auto g_0_xxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 110);

    auto g_0_xxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 111);

    auto g_0_xxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 112);

    auto g_0_xxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 113);

    auto g_0_xxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 114);

    auto g_0_xxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 115);

    auto g_0_xxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 116);

    auto g_0_xxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 117);

    auto g_0_xxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 118);

    auto g_0_xxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 119);

    auto g_0_xxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 120);

    auto g_0_xxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 121);

    auto g_0_xxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 122);

    auto g_0_xxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 123);

    auto g_0_xxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 124);

    auto g_0_xxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 125);

    auto g_0_xyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 126);

    auto g_0_xyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 127);

    auto g_0_xyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 129);

    auto g_0_xyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 130);

    auto g_0_xyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 132);

    auto g_0_xyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 133);

    auto g_0_xyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 134);

    auto g_0_xyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 136);

    auto g_0_xyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 137);

    auto g_0_xyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 138);

    auto g_0_xyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 139);

    auto g_0_xyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 141);

    auto g_0_xyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 142);

    auto g_0_xyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 143);

    auto g_0_xyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 144);

    auto g_0_xyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 145);

    auto g_0_xyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 146);

    auto g_0_xzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 189);

    auto g_0_xzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 191);

    auto g_0_xzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 193);

    auto g_0_xzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 194);

    auto g_0_xzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 196);

    auto g_0_xzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 197);

    auto g_0_xzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 198);

    auto g_0_xzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 200);

    auto g_0_xzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 201);

    auto g_0_xzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 202);

    auto g_0_xzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 203);

    auto g_0_xzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 204);

    auto g_0_xzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 205);

    auto g_0_xzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 206);

    auto g_0_xzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 207);

    auto g_0_xzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 208);

    auto g_0_xzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 209);

    auto g_0_yyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 210);

    auto g_0_yyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 211);

    auto g_0_yyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 212);

    auto g_0_yyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 213);

    auto g_0_yyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 214);

    auto g_0_yyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 215);

    auto g_0_yyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 216);

    auto g_0_yyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 217);

    auto g_0_yyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 218);

    auto g_0_yyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 219);

    auto g_0_yyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 220);

    auto g_0_yyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 221);

    auto g_0_yyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 222);

    auto g_0_yyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 223);

    auto g_0_yyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 224);

    auto g_0_yyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 225);

    auto g_0_yyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 226);

    auto g_0_yyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 227);

    auto g_0_yyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 228);

    auto g_0_yyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 229);

    auto g_0_yyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 230);

    auto g_0_yyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 232);

    auto g_0_yyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 233);

    auto g_0_yyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 234);

    auto g_0_yyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 235);

    auto g_0_yyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 236);

    auto g_0_yyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 237);

    auto g_0_yyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 238);

    auto g_0_yyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 239);

    auto g_0_yyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 240);

    auto g_0_yyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 241);

    auto g_0_yyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 242);

    auto g_0_yyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 243);

    auto g_0_yyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 244);

    auto g_0_yyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 245);

    auto g_0_yyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 246);

    auto g_0_yyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 247);

    auto g_0_yyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 248);

    auto g_0_yyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 249);

    auto g_0_yyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 250);

    auto g_0_yyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 251);

    auto g_0_yyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 252);

    auto g_0_yyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 253);

    auto g_0_yyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 254);

    auto g_0_yyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 255);

    auto g_0_yyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 256);

    auto g_0_yyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 257);

    auto g_0_yyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 258);

    auto g_0_yyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 259);

    auto g_0_yyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 260);

    auto g_0_yyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 261);

    auto g_0_yyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 262);

    auto g_0_yyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 263);

    auto g_0_yyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 264);

    auto g_0_yyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 265);

    auto g_0_yyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 266);

    auto g_0_yyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 267);

    auto g_0_yyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 268);

    auto g_0_yyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 269);

    auto g_0_yyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 270);

    auto g_0_yyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 271);

    auto g_0_yyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 272);

    auto g_0_yzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 273);

    auto g_0_yzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 274);

    auto g_0_yzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 275);

    auto g_0_yzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 276);

    auto g_0_yzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 277);

    auto g_0_yzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 278);

    auto g_0_yzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 279);

    auto g_0_yzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 280);

    auto g_0_yzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 281);

    auto g_0_yzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 282);

    auto g_0_yzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 283);

    auto g_0_yzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 284);

    auto g_0_yzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 285);

    auto g_0_yzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 286);

    auto g_0_yzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 287);

    auto g_0_yzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 288);

    auto g_0_yzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 289);

    auto g_0_yzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 290);

    auto g_0_yzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 291);

    auto g_0_yzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 292);

    auto g_0_yzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 293);

    auto g_0_zzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 294);

    auto g_0_zzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 295);

    auto g_0_zzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 296);

    auto g_0_zzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 297);

    auto g_0_zzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 298);

    auto g_0_zzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 299);

    auto g_0_zzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 300);

    auto g_0_zzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 301);

    auto g_0_zzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 302);

    auto g_0_zzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 303);

    auto g_0_zzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 304);

    auto g_0_zzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 305);

    auto g_0_zzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 306);

    auto g_0_zzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 307);

    auto g_0_zzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 308);

    auto g_0_zzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 309);

    auto g_0_zzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 310);

    auto g_0_zzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 311);

    auto g_0_zzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 312);

    auto g_0_zzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 313);

    auto g_0_zzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 314);

    /// Set up components of auxilary buffer : SGSH

    auto g_0_xxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh);

    auto g_0_xxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 1);

    auto g_0_xxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 2);

    auto g_0_xxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 3);

    auto g_0_xxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 4);

    auto g_0_xxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 5);

    auto g_0_xxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 6);

    auto g_0_xxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 7);

    auto g_0_xxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 8);

    auto g_0_xxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 9);

    auto g_0_xxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 10);

    auto g_0_xxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 11);

    auto g_0_xxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 12);

    auto g_0_xxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 13);

    auto g_0_xxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 14);

    auto g_0_xxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 15);

    auto g_0_xxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 16);

    auto g_0_xxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 17);

    auto g_0_xxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 18);

    auto g_0_xxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 19);

    auto g_0_xxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 20);

    auto g_0_xxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 21);

    auto g_0_xxxy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 22);

    auto g_0_xxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 23);

    auto g_0_xxxy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 24);

    auto g_0_xxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 26);

    auto g_0_xxxy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 27);

    auto g_0_xxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 30);

    auto g_0_xxxy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 31);

    auto g_0_xxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 35);

    auto g_0_xxxy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 36);

    auto g_0_xxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 42);

    auto g_0_xxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 43);

    auto g_0_xxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 44);

    auto g_0_xxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 45);

    auto g_0_xxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 46);

    auto g_0_xxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 47);

    auto g_0_xxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 48);

    auto g_0_xxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 49);

    auto g_0_xxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 50);

    auto g_0_xxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 51);

    auto g_0_xxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 52);

    auto g_0_xxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 53);

    auto g_0_xxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 54);

    auto g_0_xxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 55);

    auto g_0_xxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 56);

    auto g_0_xxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 58);

    auto g_0_xxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 59);

    auto g_0_xxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 60);

    auto g_0_xxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 61);

    auto g_0_xxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 62);

    auto g_0_xxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 63);

    auto g_0_xxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 64);

    auto g_0_xxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 65);

    auto g_0_xxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 66);

    auto g_0_xxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 67);

    auto g_0_xxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 68);

    auto g_0_xxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 69);

    auto g_0_xxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 70);

    auto g_0_xxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 71);

    auto g_0_xxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 72);

    auto g_0_xxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 73);

    auto g_0_xxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 74);

    auto g_0_xxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 75);

    auto g_0_xxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 76);

    auto g_0_xxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 77);

    auto g_0_xxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 78);

    auto g_0_xxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 79);

    auto g_0_xxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 80);

    auto g_0_xxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 81);

    auto g_0_xxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 82);

    auto g_0_xxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 83);

    auto g_0_xxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 105);

    auto g_0_xxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 106);

    auto g_0_xxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 107);

    auto g_0_xxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 108);

    auto g_0_xxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 109);

    auto g_0_xxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 110);

    auto g_0_xxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 111);

    auto g_0_xxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 112);

    auto g_0_xxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 113);

    auto g_0_xxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 114);

    auto g_0_xxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 115);

    auto g_0_xxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 116);

    auto g_0_xxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 117);

    auto g_0_xxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 118);

    auto g_0_xxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 119);

    auto g_0_xxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 120);

    auto g_0_xxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 121);

    auto g_0_xxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 122);

    auto g_0_xxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 123);

    auto g_0_xxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 124);

    auto g_0_xxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 125);

    auto g_0_xyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 126);

    auto g_0_xyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 127);

    auto g_0_xyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 129);

    auto g_0_xyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 130);

    auto g_0_xyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 132);

    auto g_0_xyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 133);

    auto g_0_xyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 134);

    auto g_0_xyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 136);

    auto g_0_xyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 137);

    auto g_0_xyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 138);

    auto g_0_xyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 139);

    auto g_0_xyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 141);

    auto g_0_xyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 142);

    auto g_0_xyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 143);

    auto g_0_xyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 144);

    auto g_0_xyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 145);

    auto g_0_xyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 146);

    auto g_0_xzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 189);

    auto g_0_xzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 191);

    auto g_0_xzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 193);

    auto g_0_xzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 194);

    auto g_0_xzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 196);

    auto g_0_xzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 197);

    auto g_0_xzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 198);

    auto g_0_xzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 200);

    auto g_0_xzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 201);

    auto g_0_xzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 202);

    auto g_0_xzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 203);

    auto g_0_xzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 204);

    auto g_0_xzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 205);

    auto g_0_xzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 206);

    auto g_0_xzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 207);

    auto g_0_xzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 208);

    auto g_0_xzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 209);

    auto g_0_yyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 210);

    auto g_0_yyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 211);

    auto g_0_yyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 212);

    auto g_0_yyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 213);

    auto g_0_yyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 214);

    auto g_0_yyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 215);

    auto g_0_yyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 216);

    auto g_0_yyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 217);

    auto g_0_yyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 218);

    auto g_0_yyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 219);

    auto g_0_yyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 220);

    auto g_0_yyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 221);

    auto g_0_yyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 222);

    auto g_0_yyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 223);

    auto g_0_yyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 224);

    auto g_0_yyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 225);

    auto g_0_yyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 226);

    auto g_0_yyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 227);

    auto g_0_yyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 228);

    auto g_0_yyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 229);

    auto g_0_yyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 230);

    auto g_0_yyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 232);

    auto g_0_yyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 233);

    auto g_0_yyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 234);

    auto g_0_yyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 235);

    auto g_0_yyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 236);

    auto g_0_yyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 237);

    auto g_0_yyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 238);

    auto g_0_yyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 239);

    auto g_0_yyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 240);

    auto g_0_yyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 241);

    auto g_0_yyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 242);

    auto g_0_yyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 243);

    auto g_0_yyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 244);

    auto g_0_yyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 245);

    auto g_0_yyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 246);

    auto g_0_yyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 247);

    auto g_0_yyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 248);

    auto g_0_yyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 249);

    auto g_0_yyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 250);

    auto g_0_yyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 251);

    auto g_0_yyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 252);

    auto g_0_yyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 253);

    auto g_0_yyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 254);

    auto g_0_yyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 255);

    auto g_0_yyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 256);

    auto g_0_yyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 257);

    auto g_0_yyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 258);

    auto g_0_yyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 259);

    auto g_0_yyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 260);

    auto g_0_yyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 261);

    auto g_0_yyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 262);

    auto g_0_yyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 263);

    auto g_0_yyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 264);

    auto g_0_yyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 265);

    auto g_0_yyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 266);

    auto g_0_yyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 267);

    auto g_0_yyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 268);

    auto g_0_yyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 269);

    auto g_0_yyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 270);

    auto g_0_yyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 271);

    auto g_0_yyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 272);

    auto g_0_yzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 273);

    auto g_0_yzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 274);

    auto g_0_yzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 275);

    auto g_0_yzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 276);

    auto g_0_yzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 277);

    auto g_0_yzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 278);

    auto g_0_yzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 279);

    auto g_0_yzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 280);

    auto g_0_yzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 281);

    auto g_0_yzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 282);

    auto g_0_yzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 283);

    auto g_0_yzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 284);

    auto g_0_yzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 285);

    auto g_0_yzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 286);

    auto g_0_yzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 287);

    auto g_0_yzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 288);

    auto g_0_yzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 289);

    auto g_0_yzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 290);

    auto g_0_yzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 291);

    auto g_0_yzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 292);

    auto g_0_yzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 293);

    auto g_0_zzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sgsh + 294);

    auto g_0_zzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sgsh + 295);

    auto g_0_zzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sgsh + 296);

    auto g_0_zzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sgsh + 297);

    auto g_0_zzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sgsh + 298);

    auto g_0_zzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sgsh + 299);

    auto g_0_zzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sgsh + 300);

    auto g_0_zzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sgsh + 301);

    auto g_0_zzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sgsh + 302);

    auto g_0_zzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sgsh + 303);

    auto g_0_zzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 304);

    auto g_0_zzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 305);

    auto g_0_zzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 306);

    auto g_0_zzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 307);

    auto g_0_zzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 308);

    auto g_0_zzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sgsh + 309);

    auto g_0_zzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sgsh + 310);

    auto g_0_zzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sgsh + 311);

    auto g_0_zzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sgsh + 312);

    auto g_0_zzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 313);

    auto g_0_zzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sgsh + 314);

    /// Set up 0-21 components of targeted buffer : SHSH

    auto g_0_xxxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh);

    auto g_0_xxxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 1);

    auto g_0_xxxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 2);

    auto g_0_xxxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 3);

    auto g_0_xxxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 4);

    auto g_0_xxxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 5);

    auto g_0_xxxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 6);

    auto g_0_xxxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 7);

    auto g_0_xxxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 8);

    auto g_0_xxxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 9);

    auto g_0_xxxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 10);

    auto g_0_xxxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 11);

    auto g_0_xxxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 12);

    auto g_0_xxxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 13);

    auto g_0_xxxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 14);

    auto g_0_xxxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 15);

    auto g_0_xxxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 16);

    auto g_0_xxxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 17);

    auto g_0_xxxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 18);

    auto g_0_xxxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 19);

    auto g_0_xxxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 20);

    #pragma omp simd aligned(g_0_xxx_0_xxxxx_0, g_0_xxx_0_xxxxx_1, g_0_xxx_0_xxxxy_0, g_0_xxx_0_xxxxy_1, g_0_xxx_0_xxxxz_0, g_0_xxx_0_xxxxz_1, g_0_xxx_0_xxxyy_0, g_0_xxx_0_xxxyy_1, g_0_xxx_0_xxxyz_0, g_0_xxx_0_xxxyz_1, g_0_xxx_0_xxxzz_0, g_0_xxx_0_xxxzz_1, g_0_xxx_0_xxyyy_0, g_0_xxx_0_xxyyy_1, g_0_xxx_0_xxyyz_0, g_0_xxx_0_xxyyz_1, g_0_xxx_0_xxyzz_0, g_0_xxx_0_xxyzz_1, g_0_xxx_0_xxzzz_0, g_0_xxx_0_xxzzz_1, g_0_xxx_0_xyyyy_0, g_0_xxx_0_xyyyy_1, g_0_xxx_0_xyyyz_0, g_0_xxx_0_xyyyz_1, g_0_xxx_0_xyyzz_0, g_0_xxx_0_xyyzz_1, g_0_xxx_0_xyzzz_0, g_0_xxx_0_xyzzz_1, g_0_xxx_0_xzzzz_0, g_0_xxx_0_xzzzz_1, g_0_xxx_0_yyyyy_0, g_0_xxx_0_yyyyy_1, g_0_xxx_0_yyyyz_0, g_0_xxx_0_yyyyz_1, g_0_xxx_0_yyyzz_0, g_0_xxx_0_yyyzz_1, g_0_xxx_0_yyzzz_0, g_0_xxx_0_yyzzz_1, g_0_xxx_0_yzzzz_0, g_0_xxx_0_yzzzz_1, g_0_xxx_0_zzzzz_0, g_0_xxx_0_zzzzz_1, g_0_xxxx_0_xxxx_1, g_0_xxxx_0_xxxxx_0, g_0_xxxx_0_xxxxx_1, g_0_xxxx_0_xxxxy_0, g_0_xxxx_0_xxxxy_1, g_0_xxxx_0_xxxxz_0, g_0_xxxx_0_xxxxz_1, g_0_xxxx_0_xxxy_1, g_0_xxxx_0_xxxyy_0, g_0_xxxx_0_xxxyy_1, g_0_xxxx_0_xxxyz_0, g_0_xxxx_0_xxxyz_1, g_0_xxxx_0_xxxz_1, g_0_xxxx_0_xxxzz_0, g_0_xxxx_0_xxxzz_1, g_0_xxxx_0_xxyy_1, g_0_xxxx_0_xxyyy_0, g_0_xxxx_0_xxyyy_1, g_0_xxxx_0_xxyyz_0, g_0_xxxx_0_xxyyz_1, g_0_xxxx_0_xxyz_1, g_0_xxxx_0_xxyzz_0, g_0_xxxx_0_xxyzz_1, g_0_xxxx_0_xxzz_1, g_0_xxxx_0_xxzzz_0, g_0_xxxx_0_xxzzz_1, g_0_xxxx_0_xyyy_1, g_0_xxxx_0_xyyyy_0, g_0_xxxx_0_xyyyy_1, g_0_xxxx_0_xyyyz_0, g_0_xxxx_0_xyyyz_1, g_0_xxxx_0_xyyz_1, g_0_xxxx_0_xyyzz_0, g_0_xxxx_0_xyyzz_1, g_0_xxxx_0_xyzz_1, g_0_xxxx_0_xyzzz_0, g_0_xxxx_0_xyzzz_1, g_0_xxxx_0_xzzz_1, g_0_xxxx_0_xzzzz_0, g_0_xxxx_0_xzzzz_1, g_0_xxxx_0_yyyy_1, g_0_xxxx_0_yyyyy_0, g_0_xxxx_0_yyyyy_1, g_0_xxxx_0_yyyyz_0, g_0_xxxx_0_yyyyz_1, g_0_xxxx_0_yyyz_1, g_0_xxxx_0_yyyzz_0, g_0_xxxx_0_yyyzz_1, g_0_xxxx_0_yyzz_1, g_0_xxxx_0_yyzzz_0, g_0_xxxx_0_yyzzz_1, g_0_xxxx_0_yzzz_1, g_0_xxxx_0_yzzzz_0, g_0_xxxx_0_yzzzz_1, g_0_xxxx_0_zzzz_1, g_0_xxxx_0_zzzzz_0, g_0_xxxx_0_zzzzz_1, g_0_xxxxx_0_xxxxx_0, g_0_xxxxx_0_xxxxy_0, g_0_xxxxx_0_xxxxz_0, g_0_xxxxx_0_xxxyy_0, g_0_xxxxx_0_xxxyz_0, g_0_xxxxx_0_xxxzz_0, g_0_xxxxx_0_xxyyy_0, g_0_xxxxx_0_xxyyz_0, g_0_xxxxx_0_xxyzz_0, g_0_xxxxx_0_xxzzz_0, g_0_xxxxx_0_xyyyy_0, g_0_xxxxx_0_xyyyz_0, g_0_xxxxx_0_xyyzz_0, g_0_xxxxx_0_xyzzz_0, g_0_xxxxx_0_xzzzz_0, g_0_xxxxx_0_yyyyy_0, g_0_xxxxx_0_yyyyz_0, g_0_xxxxx_0_yyyzz_0, g_0_xxxxx_0_yyzzz_0, g_0_xxxxx_0_yzzzz_0, g_0_xxxxx_0_zzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_xxxxx_0[i] = 4.0 * g_0_xxx_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxx_1[i] * fti_ab_0 + 5.0 * g_0_xxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxx_0[i] * pb_x + g_0_xxxx_0_xxxxx_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxy_0[i] = 4.0 * g_0_xxx_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxy_0[i] * pb_x + g_0_xxxx_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxz_0[i] = 4.0 * g_0_xxx_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxz_0[i] * pb_x + g_0_xxxx_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyy_0[i] = 4.0 * g_0_xxx_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyy_0[i] * pb_x + g_0_xxxx_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyz_0[i] = 4.0 * g_0_xxx_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyz_0[i] * pb_x + g_0_xxxx_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxzz_0[i] = 4.0 * g_0_xxx_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxzz_0[i] * pb_x + g_0_xxxx_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyy_0[i] = 4.0 * g_0_xxx_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyy_0[i] * pb_x + g_0_xxxx_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyz_0[i] = 4.0 * g_0_xxx_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyz_0[i] * pb_x + g_0_xxxx_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyzz_0[i] = 4.0 * g_0_xxx_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzz_0[i] * pb_x + g_0_xxxx_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxzzz_0[i] = 4.0 * g_0_xxx_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzzz_0[i] * pb_x + g_0_xxxx_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyy_0[i] = 4.0 * g_0_xxx_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyy_0[i] * pb_x + g_0_xxxx_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyz_0[i] = 4.0 * g_0_xxx_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyz_0[i] * pb_x + g_0_xxxx_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyzz_0[i] = 4.0 * g_0_xxx_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzz_0[i] * pb_x + g_0_xxxx_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyzzz_0[i] = 4.0 * g_0_xxx_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzz_0[i] * pb_x + g_0_xxxx_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xzzzz_0[i] = 4.0 * g_0_xxx_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzzzz_0[i] * pb_x + g_0_xxxx_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyy_0[i] = 4.0 * g_0_xxx_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyy_0[i] * pb_x + g_0_xxxx_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyz_0[i] = 4.0 * g_0_xxx_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyz_0[i] * pb_x + g_0_xxxx_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyzz_0[i] = 4.0 * g_0_xxx_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyzz_0[i] * pb_x + g_0_xxxx_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyzzz_0[i] = 4.0 * g_0_xxx_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzzz_0[i] * pb_x + g_0_xxxx_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yzzzz_0[i] = 4.0 * g_0_xxx_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzzz_0[i] * pb_x + g_0_xxxx_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_zzzzz_0[i] = 4.0 * g_0_xxx_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzzz_0[i] * pb_x + g_0_xxxx_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 21-42 components of targeted buffer : SHSH

    auto g_0_xxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 21);

    auto g_0_xxxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 22);

    auto g_0_xxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 23);

    auto g_0_xxxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 24);

    auto g_0_xxxxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 25);

    auto g_0_xxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 26);

    auto g_0_xxxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 27);

    auto g_0_xxxxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 28);

    auto g_0_xxxxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 29);

    auto g_0_xxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 30);

    auto g_0_xxxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 31);

    auto g_0_xxxxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 32);

    auto g_0_xxxxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 33);

    auto g_0_xxxxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 34);

    auto g_0_xxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 35);

    auto g_0_xxxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 36);

    auto g_0_xxxxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 37);

    auto g_0_xxxxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 38);

    auto g_0_xxxxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 39);

    auto g_0_xxxxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 40);

    auto g_0_xxxxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 41);

    #pragma omp simd aligned(g_0_xxxx_0_xxxx_1, g_0_xxxx_0_xxxxx_0, g_0_xxxx_0_xxxxx_1, g_0_xxxx_0_xxxxy_0, g_0_xxxx_0_xxxxy_1, g_0_xxxx_0_xxxxz_0, g_0_xxxx_0_xxxxz_1, g_0_xxxx_0_xxxy_1, g_0_xxxx_0_xxxyy_0, g_0_xxxx_0_xxxyy_1, g_0_xxxx_0_xxxyz_0, g_0_xxxx_0_xxxyz_1, g_0_xxxx_0_xxxz_1, g_0_xxxx_0_xxxzz_0, g_0_xxxx_0_xxxzz_1, g_0_xxxx_0_xxyy_1, g_0_xxxx_0_xxyyy_0, g_0_xxxx_0_xxyyy_1, g_0_xxxx_0_xxyyz_0, g_0_xxxx_0_xxyyz_1, g_0_xxxx_0_xxyz_1, g_0_xxxx_0_xxyzz_0, g_0_xxxx_0_xxyzz_1, g_0_xxxx_0_xxzz_1, g_0_xxxx_0_xxzzz_0, g_0_xxxx_0_xxzzz_1, g_0_xxxx_0_xyyy_1, g_0_xxxx_0_xyyyy_0, g_0_xxxx_0_xyyyy_1, g_0_xxxx_0_xyyyz_0, g_0_xxxx_0_xyyyz_1, g_0_xxxx_0_xyyz_1, g_0_xxxx_0_xyyzz_0, g_0_xxxx_0_xyyzz_1, g_0_xxxx_0_xyzz_1, g_0_xxxx_0_xyzzz_0, g_0_xxxx_0_xyzzz_1, g_0_xxxx_0_xzzz_1, g_0_xxxx_0_xzzzz_0, g_0_xxxx_0_xzzzz_1, g_0_xxxx_0_yyyy_1, g_0_xxxx_0_yyyyy_0, g_0_xxxx_0_yyyyy_1, g_0_xxxx_0_yyyyz_0, g_0_xxxx_0_yyyyz_1, g_0_xxxx_0_yyyz_1, g_0_xxxx_0_yyyzz_0, g_0_xxxx_0_yyyzz_1, g_0_xxxx_0_yyzz_1, g_0_xxxx_0_yyzzz_0, g_0_xxxx_0_yyzzz_1, g_0_xxxx_0_yzzz_1, g_0_xxxx_0_yzzzz_0, g_0_xxxx_0_yzzzz_1, g_0_xxxx_0_zzzz_1, g_0_xxxx_0_zzzzz_0, g_0_xxxx_0_zzzzz_1, g_0_xxxxy_0_xxxxx_0, g_0_xxxxy_0_xxxxy_0, g_0_xxxxy_0_xxxxz_0, g_0_xxxxy_0_xxxyy_0, g_0_xxxxy_0_xxxyz_0, g_0_xxxxy_0_xxxzz_0, g_0_xxxxy_0_xxyyy_0, g_0_xxxxy_0_xxyyz_0, g_0_xxxxy_0_xxyzz_0, g_0_xxxxy_0_xxzzz_0, g_0_xxxxy_0_xyyyy_0, g_0_xxxxy_0_xyyyz_0, g_0_xxxxy_0_xyyzz_0, g_0_xxxxy_0_xyzzz_0, g_0_xxxxy_0_xzzzz_0, g_0_xxxxy_0_yyyyy_0, g_0_xxxxy_0_yyyyz_0, g_0_xxxxy_0_yyyzz_0, g_0_xxxxy_0_yyzzz_0, g_0_xxxxy_0_yzzzz_0, g_0_xxxxy_0_zzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxy_0_xxxxx_0[i] = g_0_xxxx_0_xxxxx_0[i] * pb_y + g_0_xxxx_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxy_0[i] = g_0_xxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxy_0[i] * pb_y + g_0_xxxx_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxz_0[i] = g_0_xxxx_0_xxxxz_0[i] * pb_y + g_0_xxxx_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyy_0[i] = 2.0 * g_0_xxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyy_0[i] * pb_y + g_0_xxxx_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyz_0[i] = g_0_xxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyz_0[i] * pb_y + g_0_xxxx_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxzz_0[i] = g_0_xxxx_0_xxxzz_0[i] * pb_y + g_0_xxxx_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyy_0[i] = 3.0 * g_0_xxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyy_0[i] * pb_y + g_0_xxxx_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyz_0[i] = 2.0 * g_0_xxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyz_0[i] * pb_y + g_0_xxxx_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyzz_0[i] = g_0_xxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzz_0[i] * pb_y + g_0_xxxx_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxzzz_0[i] = g_0_xxxx_0_xxzzz_0[i] * pb_y + g_0_xxxx_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyy_0[i] = 4.0 * g_0_xxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyy_0[i] * pb_y + g_0_xxxx_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyz_0[i] = 3.0 * g_0_xxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyz_0[i] * pb_y + g_0_xxxx_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyzz_0[i] = 2.0 * g_0_xxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzz_0[i] * pb_y + g_0_xxxx_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyzzz_0[i] = g_0_xxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzz_0[i] * pb_y + g_0_xxxx_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xzzzz_0[i] = g_0_xxxx_0_xzzzz_0[i] * pb_y + g_0_xxxx_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyy_0[i] = 5.0 * g_0_xxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyy_0[i] * pb_y + g_0_xxxx_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyz_0[i] = 4.0 * g_0_xxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyz_0[i] * pb_y + g_0_xxxx_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyzz_0[i] = 3.0 * g_0_xxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyzz_0[i] * pb_y + g_0_xxxx_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyzzz_0[i] = 2.0 * g_0_xxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzzz_0[i] * pb_y + g_0_xxxx_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yzzzz_0[i] = g_0_xxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzzz_0[i] * pb_y + g_0_xxxx_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_zzzzz_0[i] = g_0_xxxx_0_zzzzz_0[i] * pb_y + g_0_xxxx_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 42-63 components of targeted buffer : SHSH

    auto g_0_xxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 42);

    auto g_0_xxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 43);

    auto g_0_xxxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 44);

    auto g_0_xxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 45);

    auto g_0_xxxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 46);

    auto g_0_xxxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 47);

    auto g_0_xxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 48);

    auto g_0_xxxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 49);

    auto g_0_xxxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 50);

    auto g_0_xxxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 51);

    auto g_0_xxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 52);

    auto g_0_xxxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 53);

    auto g_0_xxxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 54);

    auto g_0_xxxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 55);

    auto g_0_xxxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 56);

    auto g_0_xxxxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 57);

    auto g_0_xxxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 58);

    auto g_0_xxxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 59);

    auto g_0_xxxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 60);

    auto g_0_xxxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 61);

    auto g_0_xxxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 62);

    #pragma omp simd aligned(g_0_xxxx_0_xxxx_1, g_0_xxxx_0_xxxxx_0, g_0_xxxx_0_xxxxx_1, g_0_xxxx_0_xxxxy_0, g_0_xxxx_0_xxxxy_1, g_0_xxxx_0_xxxxz_0, g_0_xxxx_0_xxxxz_1, g_0_xxxx_0_xxxy_1, g_0_xxxx_0_xxxyy_0, g_0_xxxx_0_xxxyy_1, g_0_xxxx_0_xxxyz_0, g_0_xxxx_0_xxxyz_1, g_0_xxxx_0_xxxz_1, g_0_xxxx_0_xxxzz_0, g_0_xxxx_0_xxxzz_1, g_0_xxxx_0_xxyy_1, g_0_xxxx_0_xxyyy_0, g_0_xxxx_0_xxyyy_1, g_0_xxxx_0_xxyyz_0, g_0_xxxx_0_xxyyz_1, g_0_xxxx_0_xxyz_1, g_0_xxxx_0_xxyzz_0, g_0_xxxx_0_xxyzz_1, g_0_xxxx_0_xxzz_1, g_0_xxxx_0_xxzzz_0, g_0_xxxx_0_xxzzz_1, g_0_xxxx_0_xyyy_1, g_0_xxxx_0_xyyyy_0, g_0_xxxx_0_xyyyy_1, g_0_xxxx_0_xyyyz_0, g_0_xxxx_0_xyyyz_1, g_0_xxxx_0_xyyz_1, g_0_xxxx_0_xyyzz_0, g_0_xxxx_0_xyyzz_1, g_0_xxxx_0_xyzz_1, g_0_xxxx_0_xyzzz_0, g_0_xxxx_0_xyzzz_1, g_0_xxxx_0_xzzz_1, g_0_xxxx_0_xzzzz_0, g_0_xxxx_0_xzzzz_1, g_0_xxxx_0_yyyy_1, g_0_xxxx_0_yyyyy_0, g_0_xxxx_0_yyyyy_1, g_0_xxxx_0_yyyyz_0, g_0_xxxx_0_yyyyz_1, g_0_xxxx_0_yyyz_1, g_0_xxxx_0_yyyzz_0, g_0_xxxx_0_yyyzz_1, g_0_xxxx_0_yyzz_1, g_0_xxxx_0_yyzzz_0, g_0_xxxx_0_yyzzz_1, g_0_xxxx_0_yzzz_1, g_0_xxxx_0_yzzzz_0, g_0_xxxx_0_yzzzz_1, g_0_xxxx_0_zzzz_1, g_0_xxxx_0_zzzzz_0, g_0_xxxx_0_zzzzz_1, g_0_xxxxz_0_xxxxx_0, g_0_xxxxz_0_xxxxy_0, g_0_xxxxz_0_xxxxz_0, g_0_xxxxz_0_xxxyy_0, g_0_xxxxz_0_xxxyz_0, g_0_xxxxz_0_xxxzz_0, g_0_xxxxz_0_xxyyy_0, g_0_xxxxz_0_xxyyz_0, g_0_xxxxz_0_xxyzz_0, g_0_xxxxz_0_xxzzz_0, g_0_xxxxz_0_xyyyy_0, g_0_xxxxz_0_xyyyz_0, g_0_xxxxz_0_xyyzz_0, g_0_xxxxz_0_xyzzz_0, g_0_xxxxz_0_xzzzz_0, g_0_xxxxz_0_yyyyy_0, g_0_xxxxz_0_yyyyz_0, g_0_xxxxz_0_yyyzz_0, g_0_xxxxz_0_yyzzz_0, g_0_xxxxz_0_yzzzz_0, g_0_xxxxz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxz_0_xxxxx_0[i] = g_0_xxxx_0_xxxxx_0[i] * pb_z + g_0_xxxx_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxy_0[i] = g_0_xxxx_0_xxxxy_0[i] * pb_z + g_0_xxxx_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxz_0[i] = g_0_xxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxz_0[i] * pb_z + g_0_xxxx_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyy_0[i] = g_0_xxxx_0_xxxyy_0[i] * pb_z + g_0_xxxx_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyz_0[i] = g_0_xxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyz_0[i] * pb_z + g_0_xxxx_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxzz_0[i] = 2.0 * g_0_xxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxzz_0[i] * pb_z + g_0_xxxx_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyy_0[i] = g_0_xxxx_0_xxyyy_0[i] * pb_z + g_0_xxxx_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyz_0[i] = g_0_xxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyz_0[i] * pb_z + g_0_xxxx_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyzz_0[i] = 2.0 * g_0_xxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzz_0[i] * pb_z + g_0_xxxx_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxzzz_0[i] = 3.0 * g_0_xxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzzz_0[i] * pb_z + g_0_xxxx_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyy_0[i] = g_0_xxxx_0_xyyyy_0[i] * pb_z + g_0_xxxx_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyz_0[i] = g_0_xxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyz_0[i] * pb_z + g_0_xxxx_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyzz_0[i] = 2.0 * g_0_xxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzz_0[i] * pb_z + g_0_xxxx_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyzzz_0[i] = 3.0 * g_0_xxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzz_0[i] * pb_z + g_0_xxxx_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xzzzz_0[i] = 4.0 * g_0_xxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzzzz_0[i] * pb_z + g_0_xxxx_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyy_0[i] = g_0_xxxx_0_yyyyy_0[i] * pb_z + g_0_xxxx_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyz_0[i] = g_0_xxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyz_0[i] * pb_z + g_0_xxxx_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyzz_0[i] = 2.0 * g_0_xxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyzz_0[i] * pb_z + g_0_xxxx_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyzzz_0[i] = 3.0 * g_0_xxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzzz_0[i] * pb_z + g_0_xxxx_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yzzzz_0[i] = 4.0 * g_0_xxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzzz_0[i] * pb_z + g_0_xxxx_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_zzzzz_0[i] = 5.0 * g_0_xxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_zzzzz_0[i] * pb_z + g_0_xxxx_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 63-84 components of targeted buffer : SHSH

    auto g_0_xxxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 63);

    auto g_0_xxxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 64);

    auto g_0_xxxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 65);

    auto g_0_xxxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 66);

    auto g_0_xxxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 67);

    auto g_0_xxxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 68);

    auto g_0_xxxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 69);

    auto g_0_xxxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 70);

    auto g_0_xxxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 71);

    auto g_0_xxxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 72);

    auto g_0_xxxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 73);

    auto g_0_xxxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 74);

    auto g_0_xxxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 75);

    auto g_0_xxxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 76);

    auto g_0_xxxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 77);

    auto g_0_xxxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 78);

    auto g_0_xxxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 79);

    auto g_0_xxxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 80);

    auto g_0_xxxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 81);

    auto g_0_xxxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 82);

    auto g_0_xxxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 83);

    #pragma omp simd aligned(g_0_xxx_0_xxxxx_0, g_0_xxx_0_xxxxx_1, g_0_xxx_0_xxxxz_0, g_0_xxx_0_xxxxz_1, g_0_xxx_0_xxxzz_0, g_0_xxx_0_xxxzz_1, g_0_xxx_0_xxzzz_0, g_0_xxx_0_xxzzz_1, g_0_xxx_0_xzzzz_0, g_0_xxx_0_xzzzz_1, g_0_xxxy_0_xxxxx_0, g_0_xxxy_0_xxxxx_1, g_0_xxxy_0_xxxxz_0, g_0_xxxy_0_xxxxz_1, g_0_xxxy_0_xxxzz_0, g_0_xxxy_0_xxxzz_1, g_0_xxxy_0_xxzzz_0, g_0_xxxy_0_xxzzz_1, g_0_xxxy_0_xzzzz_0, g_0_xxxy_0_xzzzz_1, g_0_xxxyy_0_xxxxx_0, g_0_xxxyy_0_xxxxy_0, g_0_xxxyy_0_xxxxz_0, g_0_xxxyy_0_xxxyy_0, g_0_xxxyy_0_xxxyz_0, g_0_xxxyy_0_xxxzz_0, g_0_xxxyy_0_xxyyy_0, g_0_xxxyy_0_xxyyz_0, g_0_xxxyy_0_xxyzz_0, g_0_xxxyy_0_xxzzz_0, g_0_xxxyy_0_xyyyy_0, g_0_xxxyy_0_xyyyz_0, g_0_xxxyy_0_xyyzz_0, g_0_xxxyy_0_xyzzz_0, g_0_xxxyy_0_xzzzz_0, g_0_xxxyy_0_yyyyy_0, g_0_xxxyy_0_yyyyz_0, g_0_xxxyy_0_yyyzz_0, g_0_xxxyy_0_yyzzz_0, g_0_xxxyy_0_yzzzz_0, g_0_xxxyy_0_zzzzz_0, g_0_xxyy_0_xxxxy_0, g_0_xxyy_0_xxxxy_1, g_0_xxyy_0_xxxy_1, g_0_xxyy_0_xxxyy_0, g_0_xxyy_0_xxxyy_1, g_0_xxyy_0_xxxyz_0, g_0_xxyy_0_xxxyz_1, g_0_xxyy_0_xxyy_1, g_0_xxyy_0_xxyyy_0, g_0_xxyy_0_xxyyy_1, g_0_xxyy_0_xxyyz_0, g_0_xxyy_0_xxyyz_1, g_0_xxyy_0_xxyz_1, g_0_xxyy_0_xxyzz_0, g_0_xxyy_0_xxyzz_1, g_0_xxyy_0_xyyy_1, g_0_xxyy_0_xyyyy_0, g_0_xxyy_0_xyyyy_1, g_0_xxyy_0_xyyyz_0, g_0_xxyy_0_xyyyz_1, g_0_xxyy_0_xyyz_1, g_0_xxyy_0_xyyzz_0, g_0_xxyy_0_xyyzz_1, g_0_xxyy_0_xyzz_1, g_0_xxyy_0_xyzzz_0, g_0_xxyy_0_xyzzz_1, g_0_xxyy_0_yyyy_1, g_0_xxyy_0_yyyyy_0, g_0_xxyy_0_yyyyy_1, g_0_xxyy_0_yyyyz_0, g_0_xxyy_0_yyyyz_1, g_0_xxyy_0_yyyz_1, g_0_xxyy_0_yyyzz_0, g_0_xxyy_0_yyyzz_1, g_0_xxyy_0_yyzz_1, g_0_xxyy_0_yyzzz_0, g_0_xxyy_0_yyzzz_1, g_0_xxyy_0_yzzz_1, g_0_xxyy_0_yzzzz_0, g_0_xxyy_0_yzzzz_1, g_0_xxyy_0_zzzzz_0, g_0_xxyy_0_zzzzz_1, g_0_xyy_0_xxxxy_0, g_0_xyy_0_xxxxy_1, g_0_xyy_0_xxxyy_0, g_0_xyy_0_xxxyy_1, g_0_xyy_0_xxxyz_0, g_0_xyy_0_xxxyz_1, g_0_xyy_0_xxyyy_0, g_0_xyy_0_xxyyy_1, g_0_xyy_0_xxyyz_0, g_0_xyy_0_xxyyz_1, g_0_xyy_0_xxyzz_0, g_0_xyy_0_xxyzz_1, g_0_xyy_0_xyyyy_0, g_0_xyy_0_xyyyy_1, g_0_xyy_0_xyyyz_0, g_0_xyy_0_xyyyz_1, g_0_xyy_0_xyyzz_0, g_0_xyy_0_xyyzz_1, g_0_xyy_0_xyzzz_0, g_0_xyy_0_xyzzz_1, g_0_xyy_0_yyyyy_0, g_0_xyy_0_yyyyy_1, g_0_xyy_0_yyyyz_0, g_0_xyy_0_yyyyz_1, g_0_xyy_0_yyyzz_0, g_0_xyy_0_yyyzz_1, g_0_xyy_0_yyzzz_0, g_0_xyy_0_yyzzz_1, g_0_xyy_0_yzzzz_0, g_0_xyy_0_yzzzz_1, g_0_xyy_0_zzzzz_0, g_0_xyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyy_0_xxxxx_0[i] = g_0_xxx_0_xxxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxx_0[i] * pb_y + g_0_xxxy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxy_0[i] = 2.0 * g_0_xyy_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxy_0[i] * pb_x + g_0_xxyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxz_0[i] = g_0_xxx_0_xxxxz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxz_0[i] * pb_y + g_0_xxxy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxyy_0[i] = 2.0 * g_0_xyy_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyy_0[i] * pb_x + g_0_xxyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyz_0[i] = 2.0 * g_0_xyy_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyz_0[i] * pb_x + g_0_xxyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxzz_0[i] = g_0_xxx_0_xxxzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxzz_0[i] * pb_y + g_0_xxxy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxyyy_0[i] = 2.0 * g_0_xyy_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyy_0[i] * pb_x + g_0_xxyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyz_0[i] = 2.0 * g_0_xyy_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyz_0[i] * pb_x + g_0_xxyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyzz_0[i] = 2.0 * g_0_xyy_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyzz_0[i] * pb_x + g_0_xxyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxzzz_0[i] = g_0_xxx_0_xxzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxzzz_0[i] * pb_y + g_0_xxxy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xyyyy_0[i] = 2.0 * g_0_xyy_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyy_0[i] * pb_x + g_0_xxyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyz_0[i] = 2.0 * g_0_xyy_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyz_0[i] * pb_x + g_0_xxyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyzz_0[i] = 2.0 * g_0_xyy_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyzz_0[i] * pb_x + g_0_xxyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyzzz_0[i] = 2.0 * g_0_xyy_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyzzz_0[i] * pb_x + g_0_xxyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xzzzz_0[i] = g_0_xxx_0_xzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xzzzz_0[i] * pb_y + g_0_xxxy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_yyyyy_0[i] = 2.0 * g_0_xyy_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyy_0[i] * pb_x + g_0_xxyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyz_0[i] = 2.0 * g_0_xyy_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyz_0[i] * pb_x + g_0_xxyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyzz_0[i] = 2.0 * g_0_xyy_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyzz_0[i] * pb_x + g_0_xxyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyzzz_0[i] = 2.0 * g_0_xyy_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzzz_0[i] * pb_x + g_0_xxyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yzzzz_0[i] = 2.0 * g_0_xyy_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzzz_0[i] * pb_x + g_0_xxyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_zzzzz_0[i] = 2.0 * g_0_xyy_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_zzzzz_0[i] * pb_x + g_0_xxyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 84-105 components of targeted buffer : SHSH

    auto g_0_xxxyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 84);

    auto g_0_xxxyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 85);

    auto g_0_xxxyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 86);

    auto g_0_xxxyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 87);

    auto g_0_xxxyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 88);

    auto g_0_xxxyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 89);

    auto g_0_xxxyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 90);

    auto g_0_xxxyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 91);

    auto g_0_xxxyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 92);

    auto g_0_xxxyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 93);

    auto g_0_xxxyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 94);

    auto g_0_xxxyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 95);

    auto g_0_xxxyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 96);

    auto g_0_xxxyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 97);

    auto g_0_xxxyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 98);

    auto g_0_xxxyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 99);

    auto g_0_xxxyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 100);

    auto g_0_xxxyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 101);

    auto g_0_xxxyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 102);

    auto g_0_xxxyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 103);

    auto g_0_xxxyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 104);

    #pragma omp simd aligned(g_0_xxxy_0_xxxxy_0, g_0_xxxy_0_xxxxy_1, g_0_xxxy_0_xxxyy_0, g_0_xxxy_0_xxxyy_1, g_0_xxxy_0_xxyyy_0, g_0_xxxy_0_xxyyy_1, g_0_xxxy_0_xyyyy_0, g_0_xxxy_0_xyyyy_1, g_0_xxxy_0_yyyyy_0, g_0_xxxy_0_yyyyy_1, g_0_xxxyz_0_xxxxx_0, g_0_xxxyz_0_xxxxy_0, g_0_xxxyz_0_xxxxz_0, g_0_xxxyz_0_xxxyy_0, g_0_xxxyz_0_xxxyz_0, g_0_xxxyz_0_xxxzz_0, g_0_xxxyz_0_xxyyy_0, g_0_xxxyz_0_xxyyz_0, g_0_xxxyz_0_xxyzz_0, g_0_xxxyz_0_xxzzz_0, g_0_xxxyz_0_xyyyy_0, g_0_xxxyz_0_xyyyz_0, g_0_xxxyz_0_xyyzz_0, g_0_xxxyz_0_xyzzz_0, g_0_xxxyz_0_xzzzz_0, g_0_xxxyz_0_yyyyy_0, g_0_xxxyz_0_yyyyz_0, g_0_xxxyz_0_yyyzz_0, g_0_xxxyz_0_yyzzz_0, g_0_xxxyz_0_yzzzz_0, g_0_xxxyz_0_zzzzz_0, g_0_xxxz_0_xxxxx_0, g_0_xxxz_0_xxxxx_1, g_0_xxxz_0_xxxxz_0, g_0_xxxz_0_xxxxz_1, g_0_xxxz_0_xxxyz_0, g_0_xxxz_0_xxxyz_1, g_0_xxxz_0_xxxz_1, g_0_xxxz_0_xxxzz_0, g_0_xxxz_0_xxxzz_1, g_0_xxxz_0_xxyyz_0, g_0_xxxz_0_xxyyz_1, g_0_xxxz_0_xxyz_1, g_0_xxxz_0_xxyzz_0, g_0_xxxz_0_xxyzz_1, g_0_xxxz_0_xxzz_1, g_0_xxxz_0_xxzzz_0, g_0_xxxz_0_xxzzz_1, g_0_xxxz_0_xyyyz_0, g_0_xxxz_0_xyyyz_1, g_0_xxxz_0_xyyz_1, g_0_xxxz_0_xyyzz_0, g_0_xxxz_0_xyyzz_1, g_0_xxxz_0_xyzz_1, g_0_xxxz_0_xyzzz_0, g_0_xxxz_0_xyzzz_1, g_0_xxxz_0_xzzz_1, g_0_xxxz_0_xzzzz_0, g_0_xxxz_0_xzzzz_1, g_0_xxxz_0_yyyyz_0, g_0_xxxz_0_yyyyz_1, g_0_xxxz_0_yyyz_1, g_0_xxxz_0_yyyzz_0, g_0_xxxz_0_yyyzz_1, g_0_xxxz_0_yyzz_1, g_0_xxxz_0_yyzzz_0, g_0_xxxz_0_yyzzz_1, g_0_xxxz_0_yzzz_1, g_0_xxxz_0_yzzzz_0, g_0_xxxz_0_yzzzz_1, g_0_xxxz_0_zzzz_1, g_0_xxxz_0_zzzzz_0, g_0_xxxz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyz_0_xxxxx_0[i] = g_0_xxxz_0_xxxxx_0[i] * pb_y + g_0_xxxz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxy_0[i] = g_0_xxxy_0_xxxxy_0[i] * pb_z + g_0_xxxy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxz_0[i] = g_0_xxxz_0_xxxxz_0[i] * pb_y + g_0_xxxz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyy_0[i] = g_0_xxxy_0_xxxyy_0[i] * pb_z + g_0_xxxy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxyz_0[i] = g_0_xxxz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyz_0[i] * pb_y + g_0_xxxz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxzz_0[i] = g_0_xxxz_0_xxxzz_0[i] * pb_y + g_0_xxxz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyy_0[i] = g_0_xxxy_0_xxyyy_0[i] * pb_z + g_0_xxxy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxyyz_0[i] = 2.0 * g_0_xxxz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyz_0[i] * pb_y + g_0_xxxz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyzz_0[i] = g_0_xxxz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyzz_0[i] * pb_y + g_0_xxxz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxzzz_0[i] = g_0_xxxz_0_xxzzz_0[i] * pb_y + g_0_xxxz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyy_0[i] = g_0_xxxy_0_xyyyy_0[i] * pb_z + g_0_xxxy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xyyyz_0[i] = 3.0 * g_0_xxxz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyz_0[i] * pb_y + g_0_xxxz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyzz_0[i] = 2.0 * g_0_xxxz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyzz_0[i] * pb_y + g_0_xxxz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyzzz_0[i] = g_0_xxxz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyzzz_0[i] * pb_y + g_0_xxxz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xzzzz_0[i] = g_0_xxxz_0_xzzzz_0[i] * pb_y + g_0_xxxz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyy_0[i] = g_0_xxxy_0_yyyyy_0[i] * pb_z + g_0_xxxy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_yyyyz_0[i] = 4.0 * g_0_xxxz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyz_0[i] * pb_y + g_0_xxxz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyzz_0[i] = 3.0 * g_0_xxxz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyzz_0[i] * pb_y + g_0_xxxz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyzzz_0[i] = 2.0 * g_0_xxxz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyzzz_0[i] * pb_y + g_0_xxxz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yzzzz_0[i] = g_0_xxxz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yzzzz_0[i] * pb_y + g_0_xxxz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_zzzzz_0[i] = g_0_xxxz_0_zzzzz_0[i] * pb_y + g_0_xxxz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 105-126 components of targeted buffer : SHSH

    auto g_0_xxxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 105);

    auto g_0_xxxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 106);

    auto g_0_xxxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 107);

    auto g_0_xxxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 108);

    auto g_0_xxxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 109);

    auto g_0_xxxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 110);

    auto g_0_xxxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 111);

    auto g_0_xxxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 112);

    auto g_0_xxxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 113);

    auto g_0_xxxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 114);

    auto g_0_xxxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 115);

    auto g_0_xxxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 116);

    auto g_0_xxxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 117);

    auto g_0_xxxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 118);

    auto g_0_xxxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 119);

    auto g_0_xxxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 120);

    auto g_0_xxxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 121);

    auto g_0_xxxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 122);

    auto g_0_xxxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 123);

    auto g_0_xxxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 124);

    auto g_0_xxxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 125);

    #pragma omp simd aligned(g_0_xxx_0_xxxxx_0, g_0_xxx_0_xxxxx_1, g_0_xxx_0_xxxxy_0, g_0_xxx_0_xxxxy_1, g_0_xxx_0_xxxyy_0, g_0_xxx_0_xxxyy_1, g_0_xxx_0_xxyyy_0, g_0_xxx_0_xxyyy_1, g_0_xxx_0_xyyyy_0, g_0_xxx_0_xyyyy_1, g_0_xxxz_0_xxxxx_0, g_0_xxxz_0_xxxxx_1, g_0_xxxz_0_xxxxy_0, g_0_xxxz_0_xxxxy_1, g_0_xxxz_0_xxxyy_0, g_0_xxxz_0_xxxyy_1, g_0_xxxz_0_xxyyy_0, g_0_xxxz_0_xxyyy_1, g_0_xxxz_0_xyyyy_0, g_0_xxxz_0_xyyyy_1, g_0_xxxzz_0_xxxxx_0, g_0_xxxzz_0_xxxxy_0, g_0_xxxzz_0_xxxxz_0, g_0_xxxzz_0_xxxyy_0, g_0_xxxzz_0_xxxyz_0, g_0_xxxzz_0_xxxzz_0, g_0_xxxzz_0_xxyyy_0, g_0_xxxzz_0_xxyyz_0, g_0_xxxzz_0_xxyzz_0, g_0_xxxzz_0_xxzzz_0, g_0_xxxzz_0_xyyyy_0, g_0_xxxzz_0_xyyyz_0, g_0_xxxzz_0_xyyzz_0, g_0_xxxzz_0_xyzzz_0, g_0_xxxzz_0_xzzzz_0, g_0_xxxzz_0_yyyyy_0, g_0_xxxzz_0_yyyyz_0, g_0_xxxzz_0_yyyzz_0, g_0_xxxzz_0_yyzzz_0, g_0_xxxzz_0_yzzzz_0, g_0_xxxzz_0_zzzzz_0, g_0_xxzz_0_xxxxz_0, g_0_xxzz_0_xxxxz_1, g_0_xxzz_0_xxxyz_0, g_0_xxzz_0_xxxyz_1, g_0_xxzz_0_xxxz_1, g_0_xxzz_0_xxxzz_0, g_0_xxzz_0_xxxzz_1, g_0_xxzz_0_xxyyz_0, g_0_xxzz_0_xxyyz_1, g_0_xxzz_0_xxyz_1, g_0_xxzz_0_xxyzz_0, g_0_xxzz_0_xxyzz_1, g_0_xxzz_0_xxzz_1, g_0_xxzz_0_xxzzz_0, g_0_xxzz_0_xxzzz_1, g_0_xxzz_0_xyyyz_0, g_0_xxzz_0_xyyyz_1, g_0_xxzz_0_xyyz_1, g_0_xxzz_0_xyyzz_0, g_0_xxzz_0_xyyzz_1, g_0_xxzz_0_xyzz_1, g_0_xxzz_0_xyzzz_0, g_0_xxzz_0_xyzzz_1, g_0_xxzz_0_xzzz_1, g_0_xxzz_0_xzzzz_0, g_0_xxzz_0_xzzzz_1, g_0_xxzz_0_yyyyy_0, g_0_xxzz_0_yyyyy_1, g_0_xxzz_0_yyyyz_0, g_0_xxzz_0_yyyyz_1, g_0_xxzz_0_yyyz_1, g_0_xxzz_0_yyyzz_0, g_0_xxzz_0_yyyzz_1, g_0_xxzz_0_yyzz_1, g_0_xxzz_0_yyzzz_0, g_0_xxzz_0_yyzzz_1, g_0_xxzz_0_yzzz_1, g_0_xxzz_0_yzzzz_0, g_0_xxzz_0_yzzzz_1, g_0_xxzz_0_zzzz_1, g_0_xxzz_0_zzzzz_0, g_0_xxzz_0_zzzzz_1, g_0_xzz_0_xxxxz_0, g_0_xzz_0_xxxxz_1, g_0_xzz_0_xxxyz_0, g_0_xzz_0_xxxyz_1, g_0_xzz_0_xxxzz_0, g_0_xzz_0_xxxzz_1, g_0_xzz_0_xxyyz_0, g_0_xzz_0_xxyyz_1, g_0_xzz_0_xxyzz_0, g_0_xzz_0_xxyzz_1, g_0_xzz_0_xxzzz_0, g_0_xzz_0_xxzzz_1, g_0_xzz_0_xyyyz_0, g_0_xzz_0_xyyyz_1, g_0_xzz_0_xyyzz_0, g_0_xzz_0_xyyzz_1, g_0_xzz_0_xyzzz_0, g_0_xzz_0_xyzzz_1, g_0_xzz_0_xzzzz_0, g_0_xzz_0_xzzzz_1, g_0_xzz_0_yyyyy_0, g_0_xzz_0_yyyyy_1, g_0_xzz_0_yyyyz_0, g_0_xzz_0_yyyyz_1, g_0_xzz_0_yyyzz_0, g_0_xzz_0_yyyzz_1, g_0_xzz_0_yyzzz_0, g_0_xzz_0_yyzzz_1, g_0_xzz_0_yzzzz_0, g_0_xzz_0_yzzzz_1, g_0_xzz_0_zzzzz_0, g_0_xzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzz_0_xxxxx_0[i] = g_0_xxx_0_xxxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxx_0[i] * pb_z + g_0_xxxz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxy_0[i] = g_0_xxx_0_xxxxy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxy_0[i] * pb_z + g_0_xxxz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxz_0[i] = 2.0 * g_0_xzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxz_0[i] * pb_x + g_0_xxzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyy_0[i] = g_0_xxx_0_xxxyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxyy_0[i] * pb_z + g_0_xxxz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxyz_0[i] = 2.0 * g_0_xzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyz_0[i] * pb_x + g_0_xxzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxzz_0[i] = 2.0 * g_0_xzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxzz_0[i] * pb_x + g_0_xxzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyy_0[i] = g_0_xxx_0_xxyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxyyy_0[i] * pb_z + g_0_xxxz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxyyz_0[i] = 2.0 * g_0_xzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyz_0[i] * pb_x + g_0_xxzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyzz_0[i] = 2.0 * g_0_xzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyzz_0[i] * pb_x + g_0_xxzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxzzz_0[i] = 2.0 * g_0_xzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxzzz_0[i] * pb_x + g_0_xxzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyy_0[i] = g_0_xxx_0_xyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xyyyy_0[i] * pb_z + g_0_xxxz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xyyyz_0[i] = 2.0 * g_0_xzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyz_0[i] * pb_x + g_0_xxzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyzz_0[i] = 2.0 * g_0_xzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyzz_0[i] * pb_x + g_0_xxzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyzzz_0[i] = 2.0 * g_0_xzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyzzz_0[i] * pb_x + g_0_xxzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xzzzz_0[i] = 2.0 * g_0_xzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xzzzz_0[i] * pb_x + g_0_xxzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyy_0[i] = 2.0 * g_0_xzz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyy_0[i] * pb_x + g_0_xxzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyz_0[i] = 2.0 * g_0_xzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyz_0[i] * pb_x + g_0_xxzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyzz_0[i] = 2.0 * g_0_xzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyzz_0[i] * pb_x + g_0_xxzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyzzz_0[i] = 2.0 * g_0_xzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzzz_0[i] * pb_x + g_0_xxzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yzzzz_0[i] = 2.0 * g_0_xzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzzz_0[i] * pb_x + g_0_xxzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_zzzzz_0[i] = 2.0 * g_0_xzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzzz_0[i] * pb_x + g_0_xxzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 126-147 components of targeted buffer : SHSH

    auto g_0_xxyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 126);

    auto g_0_xxyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 127);

    auto g_0_xxyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 128);

    auto g_0_xxyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 129);

    auto g_0_xxyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 130);

    auto g_0_xxyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 131);

    auto g_0_xxyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 132);

    auto g_0_xxyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 133);

    auto g_0_xxyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 134);

    auto g_0_xxyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 135);

    auto g_0_xxyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 136);

    auto g_0_xxyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 137);

    auto g_0_xxyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 138);

    auto g_0_xxyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 139);

    auto g_0_xxyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 140);

    auto g_0_xxyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 141);

    auto g_0_xxyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 142);

    auto g_0_xxyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 143);

    auto g_0_xxyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 144);

    auto g_0_xxyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 145);

    auto g_0_xxyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 146);

    #pragma omp simd aligned(g_0_xxy_0_xxxxx_0, g_0_xxy_0_xxxxx_1, g_0_xxy_0_xxxxz_0, g_0_xxy_0_xxxxz_1, g_0_xxy_0_xxxzz_0, g_0_xxy_0_xxxzz_1, g_0_xxy_0_xxzzz_0, g_0_xxy_0_xxzzz_1, g_0_xxy_0_xzzzz_0, g_0_xxy_0_xzzzz_1, g_0_xxyy_0_xxxxx_0, g_0_xxyy_0_xxxxx_1, g_0_xxyy_0_xxxxz_0, g_0_xxyy_0_xxxxz_1, g_0_xxyy_0_xxxzz_0, g_0_xxyy_0_xxxzz_1, g_0_xxyy_0_xxzzz_0, g_0_xxyy_0_xxzzz_1, g_0_xxyy_0_xzzzz_0, g_0_xxyy_0_xzzzz_1, g_0_xxyyy_0_xxxxx_0, g_0_xxyyy_0_xxxxy_0, g_0_xxyyy_0_xxxxz_0, g_0_xxyyy_0_xxxyy_0, g_0_xxyyy_0_xxxyz_0, g_0_xxyyy_0_xxxzz_0, g_0_xxyyy_0_xxyyy_0, g_0_xxyyy_0_xxyyz_0, g_0_xxyyy_0_xxyzz_0, g_0_xxyyy_0_xxzzz_0, g_0_xxyyy_0_xyyyy_0, g_0_xxyyy_0_xyyyz_0, g_0_xxyyy_0_xyyzz_0, g_0_xxyyy_0_xyzzz_0, g_0_xxyyy_0_xzzzz_0, g_0_xxyyy_0_yyyyy_0, g_0_xxyyy_0_yyyyz_0, g_0_xxyyy_0_yyyzz_0, g_0_xxyyy_0_yyzzz_0, g_0_xxyyy_0_yzzzz_0, g_0_xxyyy_0_zzzzz_0, g_0_xyyy_0_xxxxy_0, g_0_xyyy_0_xxxxy_1, g_0_xyyy_0_xxxy_1, g_0_xyyy_0_xxxyy_0, g_0_xyyy_0_xxxyy_1, g_0_xyyy_0_xxxyz_0, g_0_xyyy_0_xxxyz_1, g_0_xyyy_0_xxyy_1, g_0_xyyy_0_xxyyy_0, g_0_xyyy_0_xxyyy_1, g_0_xyyy_0_xxyyz_0, g_0_xyyy_0_xxyyz_1, g_0_xyyy_0_xxyz_1, g_0_xyyy_0_xxyzz_0, g_0_xyyy_0_xxyzz_1, g_0_xyyy_0_xyyy_1, g_0_xyyy_0_xyyyy_0, g_0_xyyy_0_xyyyy_1, g_0_xyyy_0_xyyyz_0, g_0_xyyy_0_xyyyz_1, g_0_xyyy_0_xyyz_1, g_0_xyyy_0_xyyzz_0, g_0_xyyy_0_xyyzz_1, g_0_xyyy_0_xyzz_1, g_0_xyyy_0_xyzzz_0, g_0_xyyy_0_xyzzz_1, g_0_xyyy_0_yyyy_1, g_0_xyyy_0_yyyyy_0, g_0_xyyy_0_yyyyy_1, g_0_xyyy_0_yyyyz_0, g_0_xyyy_0_yyyyz_1, g_0_xyyy_0_yyyz_1, g_0_xyyy_0_yyyzz_0, g_0_xyyy_0_yyyzz_1, g_0_xyyy_0_yyzz_1, g_0_xyyy_0_yyzzz_0, g_0_xyyy_0_yyzzz_1, g_0_xyyy_0_yzzz_1, g_0_xyyy_0_yzzzz_0, g_0_xyyy_0_yzzzz_1, g_0_xyyy_0_zzzzz_0, g_0_xyyy_0_zzzzz_1, g_0_yyy_0_xxxxy_0, g_0_yyy_0_xxxxy_1, g_0_yyy_0_xxxyy_0, g_0_yyy_0_xxxyy_1, g_0_yyy_0_xxxyz_0, g_0_yyy_0_xxxyz_1, g_0_yyy_0_xxyyy_0, g_0_yyy_0_xxyyy_1, g_0_yyy_0_xxyyz_0, g_0_yyy_0_xxyyz_1, g_0_yyy_0_xxyzz_0, g_0_yyy_0_xxyzz_1, g_0_yyy_0_xyyyy_0, g_0_yyy_0_xyyyy_1, g_0_yyy_0_xyyyz_0, g_0_yyy_0_xyyyz_1, g_0_yyy_0_xyyzz_0, g_0_yyy_0_xyyzz_1, g_0_yyy_0_xyzzz_0, g_0_yyy_0_xyzzz_1, g_0_yyy_0_yyyyy_0, g_0_yyy_0_yyyyy_1, g_0_yyy_0_yyyyz_0, g_0_yyy_0_yyyyz_1, g_0_yyy_0_yyyzz_0, g_0_yyy_0_yyyzz_1, g_0_yyy_0_yyzzz_0, g_0_yyy_0_yyzzz_1, g_0_yyy_0_yzzzz_0, g_0_yyy_0_yzzzz_1, g_0_yyy_0_zzzzz_0, g_0_yyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyy_0_xxxxx_0[i] = 2.0 * g_0_xxy_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxx_0[i] * pb_y + g_0_xxyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxy_0[i] = g_0_yyy_0_xxxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxy_0[i] * pb_x + g_0_xyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxz_0[i] = 2.0 * g_0_xxy_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxz_0[i] * pb_y + g_0_xxyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxyy_0[i] = g_0_yyy_0_xxxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyy_0[i] * pb_x + g_0_xyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyz_0[i] = g_0_yyy_0_xxxyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyz_0[i] * pb_x + g_0_xyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxzz_0[i] = 2.0 * g_0_xxy_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxzz_0[i] * pb_y + g_0_xxyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxyyy_0[i] = g_0_yyy_0_xxyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyy_0[i] * pb_x + g_0_xyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyz_0[i] = g_0_yyy_0_xxyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyz_0[i] * pb_x + g_0_xyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyzz_0[i] = g_0_yyy_0_xxyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyzz_0[i] * pb_x + g_0_xyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxzzz_0[i] = 2.0 * g_0_xxy_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxzzz_0[i] * pb_y + g_0_xxyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xyyyy_0[i] = g_0_yyy_0_xyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyy_0[i] * pb_x + g_0_xyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyz_0[i] = g_0_yyy_0_xyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyz_0[i] * pb_x + g_0_xyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyzz_0[i] = g_0_yyy_0_xyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyzz_0[i] * pb_x + g_0_xyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyzzz_0[i] = g_0_yyy_0_xyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyzzz_0[i] * pb_x + g_0_xyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xzzzz_0[i] = 2.0 * g_0_xxy_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xzzzz_0[i] * pb_y + g_0_xxyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_yyyyy_0[i] = g_0_yyy_0_yyyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyy_0[i] * pb_x + g_0_xyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyz_0[i] = g_0_yyy_0_yyyyz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyz_0[i] * pb_x + g_0_xyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyzz_0[i] = g_0_yyy_0_yyyzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyzz_0[i] * pb_x + g_0_xyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyzzz_0[i] = g_0_yyy_0_yyzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzzz_0[i] * pb_x + g_0_xyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yzzzz_0[i] = g_0_yyy_0_yzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzzz_0[i] * pb_x + g_0_xyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_zzzzz_0[i] = g_0_yyy_0_zzzzz_0[i] * fi_ab_0 - g_0_yyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_zzzzz_0[i] * pb_x + g_0_xyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 147-168 components of targeted buffer : SHSH

    auto g_0_xxyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 147);

    auto g_0_xxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 148);

    auto g_0_xxyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 149);

    auto g_0_xxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 150);

    auto g_0_xxyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 151);

    auto g_0_xxyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 152);

    auto g_0_xxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 153);

    auto g_0_xxyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 154);

    auto g_0_xxyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 155);

    auto g_0_xxyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 156);

    auto g_0_xxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 157);

    auto g_0_xxyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 158);

    auto g_0_xxyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 159);

    auto g_0_xxyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 160);

    auto g_0_xxyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 161);

    auto g_0_xxyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 162);

    auto g_0_xxyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 163);

    auto g_0_xxyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 164);

    auto g_0_xxyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 165);

    auto g_0_xxyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 166);

    auto g_0_xxyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 167);

    #pragma omp simd aligned(g_0_xxyy_0_xxxx_1, g_0_xxyy_0_xxxxx_0, g_0_xxyy_0_xxxxx_1, g_0_xxyy_0_xxxxy_0, g_0_xxyy_0_xxxxy_1, g_0_xxyy_0_xxxxz_0, g_0_xxyy_0_xxxxz_1, g_0_xxyy_0_xxxy_1, g_0_xxyy_0_xxxyy_0, g_0_xxyy_0_xxxyy_1, g_0_xxyy_0_xxxyz_0, g_0_xxyy_0_xxxyz_1, g_0_xxyy_0_xxxz_1, g_0_xxyy_0_xxxzz_0, g_0_xxyy_0_xxxzz_1, g_0_xxyy_0_xxyy_1, g_0_xxyy_0_xxyyy_0, g_0_xxyy_0_xxyyy_1, g_0_xxyy_0_xxyyz_0, g_0_xxyy_0_xxyyz_1, g_0_xxyy_0_xxyz_1, g_0_xxyy_0_xxyzz_0, g_0_xxyy_0_xxyzz_1, g_0_xxyy_0_xxzz_1, g_0_xxyy_0_xxzzz_0, g_0_xxyy_0_xxzzz_1, g_0_xxyy_0_xyyy_1, g_0_xxyy_0_xyyyy_0, g_0_xxyy_0_xyyyy_1, g_0_xxyy_0_xyyyz_0, g_0_xxyy_0_xyyyz_1, g_0_xxyy_0_xyyz_1, g_0_xxyy_0_xyyzz_0, g_0_xxyy_0_xyyzz_1, g_0_xxyy_0_xyzz_1, g_0_xxyy_0_xyzzz_0, g_0_xxyy_0_xyzzz_1, g_0_xxyy_0_xzzz_1, g_0_xxyy_0_xzzzz_0, g_0_xxyy_0_xzzzz_1, g_0_xxyy_0_yyyy_1, g_0_xxyy_0_yyyyy_0, g_0_xxyy_0_yyyyy_1, g_0_xxyy_0_yyyyz_0, g_0_xxyy_0_yyyyz_1, g_0_xxyy_0_yyyz_1, g_0_xxyy_0_yyyzz_0, g_0_xxyy_0_yyyzz_1, g_0_xxyy_0_yyzz_1, g_0_xxyy_0_yyzzz_0, g_0_xxyy_0_yyzzz_1, g_0_xxyy_0_yzzz_1, g_0_xxyy_0_yzzzz_0, g_0_xxyy_0_yzzzz_1, g_0_xxyy_0_zzzz_1, g_0_xxyy_0_zzzzz_0, g_0_xxyy_0_zzzzz_1, g_0_xxyyz_0_xxxxx_0, g_0_xxyyz_0_xxxxy_0, g_0_xxyyz_0_xxxxz_0, g_0_xxyyz_0_xxxyy_0, g_0_xxyyz_0_xxxyz_0, g_0_xxyyz_0_xxxzz_0, g_0_xxyyz_0_xxyyy_0, g_0_xxyyz_0_xxyyz_0, g_0_xxyyz_0_xxyzz_0, g_0_xxyyz_0_xxzzz_0, g_0_xxyyz_0_xyyyy_0, g_0_xxyyz_0_xyyyz_0, g_0_xxyyz_0_xyyzz_0, g_0_xxyyz_0_xyzzz_0, g_0_xxyyz_0_xzzzz_0, g_0_xxyyz_0_yyyyy_0, g_0_xxyyz_0_yyyyz_0, g_0_xxyyz_0_yyyzz_0, g_0_xxyyz_0_yyzzz_0, g_0_xxyyz_0_yzzzz_0, g_0_xxyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyz_0_xxxxx_0[i] = g_0_xxyy_0_xxxxx_0[i] * pb_z + g_0_xxyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxy_0[i] = g_0_xxyy_0_xxxxy_0[i] * pb_z + g_0_xxyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxz_0[i] = g_0_xxyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxz_0[i] * pb_z + g_0_xxyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyy_0[i] = g_0_xxyy_0_xxxyy_0[i] * pb_z + g_0_xxyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyz_0[i] = g_0_xxyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyz_0[i] * pb_z + g_0_xxyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxzz_0[i] = 2.0 * g_0_xxyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxzz_0[i] * pb_z + g_0_xxyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyy_0[i] = g_0_xxyy_0_xxyyy_0[i] * pb_z + g_0_xxyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyz_0[i] = g_0_xxyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyz_0[i] * pb_z + g_0_xxyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyzz_0[i] = 2.0 * g_0_xxyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyzz_0[i] * pb_z + g_0_xxyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxzzz_0[i] = 3.0 * g_0_xxyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxzzz_0[i] * pb_z + g_0_xxyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyy_0[i] = g_0_xxyy_0_xyyyy_0[i] * pb_z + g_0_xxyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyz_0[i] = g_0_xxyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyz_0[i] * pb_z + g_0_xxyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyzz_0[i] = 2.0 * g_0_xxyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyzz_0[i] * pb_z + g_0_xxyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyzzz_0[i] = 3.0 * g_0_xxyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyzzz_0[i] * pb_z + g_0_xxyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xzzzz_0[i] = 4.0 * g_0_xxyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xzzzz_0[i] * pb_z + g_0_xxyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyy_0[i] = g_0_xxyy_0_yyyyy_0[i] * pb_z + g_0_xxyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyz_0[i] = g_0_xxyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyz_0[i] * pb_z + g_0_xxyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyzz_0[i] = 2.0 * g_0_xxyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyzz_0[i] * pb_z + g_0_xxyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyzzz_0[i] = 3.0 * g_0_xxyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyzzz_0[i] * pb_z + g_0_xxyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yzzzz_0[i] = 4.0 * g_0_xxyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yzzzz_0[i] * pb_z + g_0_xxyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_zzzzz_0[i] = 5.0 * g_0_xxyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_zzzzz_0[i] * pb_z + g_0_xxyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 168-189 components of targeted buffer : SHSH

    auto g_0_xxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 168);

    auto g_0_xxyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 169);

    auto g_0_xxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 170);

    auto g_0_xxyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 171);

    auto g_0_xxyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 172);

    auto g_0_xxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 173);

    auto g_0_xxyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 174);

    auto g_0_xxyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 175);

    auto g_0_xxyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 176);

    auto g_0_xxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 177);

    auto g_0_xxyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 178);

    auto g_0_xxyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 179);

    auto g_0_xxyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 180);

    auto g_0_xxyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 181);

    auto g_0_xxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 182);

    auto g_0_xxyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 183);

    auto g_0_xxyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 184);

    auto g_0_xxyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 185);

    auto g_0_xxyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 186);

    auto g_0_xxyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 187);

    auto g_0_xxyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 188);

    #pragma omp simd aligned(g_0_xxyzz_0_xxxxx_0, g_0_xxyzz_0_xxxxy_0, g_0_xxyzz_0_xxxxz_0, g_0_xxyzz_0_xxxyy_0, g_0_xxyzz_0_xxxyz_0, g_0_xxyzz_0_xxxzz_0, g_0_xxyzz_0_xxyyy_0, g_0_xxyzz_0_xxyyz_0, g_0_xxyzz_0_xxyzz_0, g_0_xxyzz_0_xxzzz_0, g_0_xxyzz_0_xyyyy_0, g_0_xxyzz_0_xyyyz_0, g_0_xxyzz_0_xyyzz_0, g_0_xxyzz_0_xyzzz_0, g_0_xxyzz_0_xzzzz_0, g_0_xxyzz_0_yyyyy_0, g_0_xxyzz_0_yyyyz_0, g_0_xxyzz_0_yyyzz_0, g_0_xxyzz_0_yyzzz_0, g_0_xxyzz_0_yzzzz_0, g_0_xxyzz_0_zzzzz_0, g_0_xxzz_0_xxxx_1, g_0_xxzz_0_xxxxx_0, g_0_xxzz_0_xxxxx_1, g_0_xxzz_0_xxxxy_0, g_0_xxzz_0_xxxxy_1, g_0_xxzz_0_xxxxz_0, g_0_xxzz_0_xxxxz_1, g_0_xxzz_0_xxxy_1, g_0_xxzz_0_xxxyy_0, g_0_xxzz_0_xxxyy_1, g_0_xxzz_0_xxxyz_0, g_0_xxzz_0_xxxyz_1, g_0_xxzz_0_xxxz_1, g_0_xxzz_0_xxxzz_0, g_0_xxzz_0_xxxzz_1, g_0_xxzz_0_xxyy_1, g_0_xxzz_0_xxyyy_0, g_0_xxzz_0_xxyyy_1, g_0_xxzz_0_xxyyz_0, g_0_xxzz_0_xxyyz_1, g_0_xxzz_0_xxyz_1, g_0_xxzz_0_xxyzz_0, g_0_xxzz_0_xxyzz_1, g_0_xxzz_0_xxzz_1, g_0_xxzz_0_xxzzz_0, g_0_xxzz_0_xxzzz_1, g_0_xxzz_0_xyyy_1, g_0_xxzz_0_xyyyy_0, g_0_xxzz_0_xyyyy_1, g_0_xxzz_0_xyyyz_0, g_0_xxzz_0_xyyyz_1, g_0_xxzz_0_xyyz_1, g_0_xxzz_0_xyyzz_0, g_0_xxzz_0_xyyzz_1, g_0_xxzz_0_xyzz_1, g_0_xxzz_0_xyzzz_0, g_0_xxzz_0_xyzzz_1, g_0_xxzz_0_xzzz_1, g_0_xxzz_0_xzzzz_0, g_0_xxzz_0_xzzzz_1, g_0_xxzz_0_yyyy_1, g_0_xxzz_0_yyyyy_0, g_0_xxzz_0_yyyyy_1, g_0_xxzz_0_yyyyz_0, g_0_xxzz_0_yyyyz_1, g_0_xxzz_0_yyyz_1, g_0_xxzz_0_yyyzz_0, g_0_xxzz_0_yyyzz_1, g_0_xxzz_0_yyzz_1, g_0_xxzz_0_yyzzz_0, g_0_xxzz_0_yyzzz_1, g_0_xxzz_0_yzzz_1, g_0_xxzz_0_yzzzz_0, g_0_xxzz_0_yzzzz_1, g_0_xxzz_0_zzzz_1, g_0_xxzz_0_zzzzz_0, g_0_xxzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzz_0_xxxxx_0[i] = g_0_xxzz_0_xxxxx_0[i] * pb_y + g_0_xxzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxy_0[i] = g_0_xxzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxy_0[i] * pb_y + g_0_xxzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxz_0[i] = g_0_xxzz_0_xxxxz_0[i] * pb_y + g_0_xxzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyy_0[i] = 2.0 * g_0_xxzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyy_0[i] * pb_y + g_0_xxzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyz_0[i] = g_0_xxzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyz_0[i] * pb_y + g_0_xxzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxzz_0[i] = g_0_xxzz_0_xxxzz_0[i] * pb_y + g_0_xxzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyy_0[i] = 3.0 * g_0_xxzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyy_0[i] * pb_y + g_0_xxzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyz_0[i] = 2.0 * g_0_xxzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyz_0[i] * pb_y + g_0_xxzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyzz_0[i] = g_0_xxzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyzz_0[i] * pb_y + g_0_xxzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxzzz_0[i] = g_0_xxzz_0_xxzzz_0[i] * pb_y + g_0_xxzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyy_0[i] = 4.0 * g_0_xxzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyy_0[i] * pb_y + g_0_xxzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyz_0[i] = 3.0 * g_0_xxzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyz_0[i] * pb_y + g_0_xxzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyzz_0[i] = 2.0 * g_0_xxzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyzz_0[i] * pb_y + g_0_xxzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyzzz_0[i] = g_0_xxzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyzzz_0[i] * pb_y + g_0_xxzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xzzzz_0[i] = g_0_xxzz_0_xzzzz_0[i] * pb_y + g_0_xxzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyy_0[i] = 5.0 * g_0_xxzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyy_0[i] * pb_y + g_0_xxzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyz_0[i] = 4.0 * g_0_xxzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyz_0[i] * pb_y + g_0_xxzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyzz_0[i] = 3.0 * g_0_xxzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyzz_0[i] * pb_y + g_0_xxzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyzzz_0[i] = 2.0 * g_0_xxzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyzzz_0[i] * pb_y + g_0_xxzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yzzzz_0[i] = g_0_xxzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yzzzz_0[i] * pb_y + g_0_xxzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_zzzzz_0[i] = g_0_xxzz_0_zzzzz_0[i] * pb_y + g_0_xxzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 189-210 components of targeted buffer : SHSH

    auto g_0_xxzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 189);

    auto g_0_xxzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 190);

    auto g_0_xxzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 191);

    auto g_0_xxzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 192);

    auto g_0_xxzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 193);

    auto g_0_xxzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 194);

    auto g_0_xxzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 195);

    auto g_0_xxzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 196);

    auto g_0_xxzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 197);

    auto g_0_xxzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 198);

    auto g_0_xxzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 199);

    auto g_0_xxzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 200);

    auto g_0_xxzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 201);

    auto g_0_xxzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 202);

    auto g_0_xxzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 203);

    auto g_0_xxzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 204);

    auto g_0_xxzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 205);

    auto g_0_xxzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 206);

    auto g_0_xxzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 207);

    auto g_0_xxzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 208);

    auto g_0_xxzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 209);

    #pragma omp simd aligned(g_0_xxz_0_xxxxx_0, g_0_xxz_0_xxxxx_1, g_0_xxz_0_xxxxy_0, g_0_xxz_0_xxxxy_1, g_0_xxz_0_xxxyy_0, g_0_xxz_0_xxxyy_1, g_0_xxz_0_xxyyy_0, g_0_xxz_0_xxyyy_1, g_0_xxz_0_xyyyy_0, g_0_xxz_0_xyyyy_1, g_0_xxzz_0_xxxxx_0, g_0_xxzz_0_xxxxx_1, g_0_xxzz_0_xxxxy_0, g_0_xxzz_0_xxxxy_1, g_0_xxzz_0_xxxyy_0, g_0_xxzz_0_xxxyy_1, g_0_xxzz_0_xxyyy_0, g_0_xxzz_0_xxyyy_1, g_0_xxzz_0_xyyyy_0, g_0_xxzz_0_xyyyy_1, g_0_xxzzz_0_xxxxx_0, g_0_xxzzz_0_xxxxy_0, g_0_xxzzz_0_xxxxz_0, g_0_xxzzz_0_xxxyy_0, g_0_xxzzz_0_xxxyz_0, g_0_xxzzz_0_xxxzz_0, g_0_xxzzz_0_xxyyy_0, g_0_xxzzz_0_xxyyz_0, g_0_xxzzz_0_xxyzz_0, g_0_xxzzz_0_xxzzz_0, g_0_xxzzz_0_xyyyy_0, g_0_xxzzz_0_xyyyz_0, g_0_xxzzz_0_xyyzz_0, g_0_xxzzz_0_xyzzz_0, g_0_xxzzz_0_xzzzz_0, g_0_xxzzz_0_yyyyy_0, g_0_xxzzz_0_yyyyz_0, g_0_xxzzz_0_yyyzz_0, g_0_xxzzz_0_yyzzz_0, g_0_xxzzz_0_yzzzz_0, g_0_xxzzz_0_zzzzz_0, g_0_xzzz_0_xxxxz_0, g_0_xzzz_0_xxxxz_1, g_0_xzzz_0_xxxyz_0, g_0_xzzz_0_xxxyz_1, g_0_xzzz_0_xxxz_1, g_0_xzzz_0_xxxzz_0, g_0_xzzz_0_xxxzz_1, g_0_xzzz_0_xxyyz_0, g_0_xzzz_0_xxyyz_1, g_0_xzzz_0_xxyz_1, g_0_xzzz_0_xxyzz_0, g_0_xzzz_0_xxyzz_1, g_0_xzzz_0_xxzz_1, g_0_xzzz_0_xxzzz_0, g_0_xzzz_0_xxzzz_1, g_0_xzzz_0_xyyyz_0, g_0_xzzz_0_xyyyz_1, g_0_xzzz_0_xyyz_1, g_0_xzzz_0_xyyzz_0, g_0_xzzz_0_xyyzz_1, g_0_xzzz_0_xyzz_1, g_0_xzzz_0_xyzzz_0, g_0_xzzz_0_xyzzz_1, g_0_xzzz_0_xzzz_1, g_0_xzzz_0_xzzzz_0, g_0_xzzz_0_xzzzz_1, g_0_xzzz_0_yyyyy_0, g_0_xzzz_0_yyyyy_1, g_0_xzzz_0_yyyyz_0, g_0_xzzz_0_yyyyz_1, g_0_xzzz_0_yyyz_1, g_0_xzzz_0_yyyzz_0, g_0_xzzz_0_yyyzz_1, g_0_xzzz_0_yyzz_1, g_0_xzzz_0_yyzzz_0, g_0_xzzz_0_yyzzz_1, g_0_xzzz_0_yzzz_1, g_0_xzzz_0_yzzzz_0, g_0_xzzz_0_yzzzz_1, g_0_xzzz_0_zzzz_1, g_0_xzzz_0_zzzzz_0, g_0_xzzz_0_zzzzz_1, g_0_zzz_0_xxxxz_0, g_0_zzz_0_xxxxz_1, g_0_zzz_0_xxxyz_0, g_0_zzz_0_xxxyz_1, g_0_zzz_0_xxxzz_0, g_0_zzz_0_xxxzz_1, g_0_zzz_0_xxyyz_0, g_0_zzz_0_xxyyz_1, g_0_zzz_0_xxyzz_0, g_0_zzz_0_xxyzz_1, g_0_zzz_0_xxzzz_0, g_0_zzz_0_xxzzz_1, g_0_zzz_0_xyyyz_0, g_0_zzz_0_xyyyz_1, g_0_zzz_0_xyyzz_0, g_0_zzz_0_xyyzz_1, g_0_zzz_0_xyzzz_0, g_0_zzz_0_xyzzz_1, g_0_zzz_0_xzzzz_0, g_0_zzz_0_xzzzz_1, g_0_zzz_0_yyyyy_0, g_0_zzz_0_yyyyy_1, g_0_zzz_0_yyyyz_0, g_0_zzz_0_yyyyz_1, g_0_zzz_0_yyyzz_0, g_0_zzz_0_yyyzz_1, g_0_zzz_0_yyzzz_0, g_0_zzz_0_yyzzz_1, g_0_zzz_0_yzzzz_0, g_0_zzz_0_yzzzz_1, g_0_zzz_0_zzzzz_0, g_0_zzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzz_0_xxxxx_0[i] = 2.0 * g_0_xxz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxx_0[i] * pb_z + g_0_xxzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxy_0[i] = 2.0 * g_0_xxz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxy_0[i] * pb_z + g_0_xxzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxz_0[i] = g_0_zzz_0_xxxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxz_0[i] * pb_x + g_0_xzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyy_0[i] = 2.0 * g_0_xxz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxyy_0[i] * pb_z + g_0_xxzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxyz_0[i] = g_0_zzz_0_xxxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxyz_0[i] * pb_x + g_0_xzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxzz_0[i] = g_0_zzz_0_xxxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxzz_0[i] * pb_x + g_0_xzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyy_0[i] = 2.0 * g_0_xxz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxyyy_0[i] * pb_z + g_0_xxzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxyyz_0[i] = g_0_zzz_0_xxyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyyz_0[i] * pb_x + g_0_xzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyzz_0[i] = g_0_zzz_0_xxyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyzz_0[i] * pb_x + g_0_xzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxzzz_0[i] = g_0_zzz_0_xxzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxzzz_0[i] * pb_x + g_0_xzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyy_0[i] = 2.0 * g_0_xxz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xyyyy_0[i] * pb_z + g_0_xxzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xyyyz_0[i] = g_0_zzz_0_xyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyyz_0[i] * pb_x + g_0_xzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyzz_0[i] = g_0_zzz_0_xyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyzz_0[i] * pb_x + g_0_xzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyzzz_0[i] = g_0_zzz_0_xyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyzzz_0[i] * pb_x + g_0_xzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xzzzz_0[i] = g_0_zzz_0_xzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xzzzz_0[i] * pb_x + g_0_xzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyy_0[i] = g_0_zzz_0_yyyyy_0[i] * fi_ab_0 - g_0_zzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyy_0[i] * pb_x + g_0_xzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyz_0[i] = g_0_zzz_0_yyyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyz_0[i] * pb_x + g_0_xzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyzz_0[i] = g_0_zzz_0_yyyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyzz_0[i] * pb_x + g_0_xzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyzzz_0[i] = g_0_zzz_0_yyzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzzz_0[i] * pb_x + g_0_xzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yzzzz_0[i] = g_0_zzz_0_yzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzzz_0[i] * pb_x + g_0_xzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_zzzzz_0[i] = g_0_zzz_0_zzzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzzz_0[i] * pb_x + g_0_xzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 210-231 components of targeted buffer : SHSH

    auto g_0_xyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 210);

    auto g_0_xyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 211);

    auto g_0_xyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 212);

    auto g_0_xyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 213);

    auto g_0_xyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 214);

    auto g_0_xyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 215);

    auto g_0_xyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 216);

    auto g_0_xyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 217);

    auto g_0_xyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 218);

    auto g_0_xyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 219);

    auto g_0_xyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 220);

    auto g_0_xyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 221);

    auto g_0_xyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 222);

    auto g_0_xyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 223);

    auto g_0_xyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 224);

    auto g_0_xyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 225);

    auto g_0_xyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 226);

    auto g_0_xyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 227);

    auto g_0_xyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 228);

    auto g_0_xyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 229);

    auto g_0_xyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 230);

    #pragma omp simd aligned(g_0_xyyyy_0_xxxxx_0, g_0_xyyyy_0_xxxxy_0, g_0_xyyyy_0_xxxxz_0, g_0_xyyyy_0_xxxyy_0, g_0_xyyyy_0_xxxyz_0, g_0_xyyyy_0_xxxzz_0, g_0_xyyyy_0_xxyyy_0, g_0_xyyyy_0_xxyyz_0, g_0_xyyyy_0_xxyzz_0, g_0_xyyyy_0_xxzzz_0, g_0_xyyyy_0_xyyyy_0, g_0_xyyyy_0_xyyyz_0, g_0_xyyyy_0_xyyzz_0, g_0_xyyyy_0_xyzzz_0, g_0_xyyyy_0_xzzzz_0, g_0_xyyyy_0_yyyyy_0, g_0_xyyyy_0_yyyyz_0, g_0_xyyyy_0_yyyzz_0, g_0_xyyyy_0_yyzzz_0, g_0_xyyyy_0_yzzzz_0, g_0_xyyyy_0_zzzzz_0, g_0_yyyy_0_xxxx_1, g_0_yyyy_0_xxxxx_0, g_0_yyyy_0_xxxxx_1, g_0_yyyy_0_xxxxy_0, g_0_yyyy_0_xxxxy_1, g_0_yyyy_0_xxxxz_0, g_0_yyyy_0_xxxxz_1, g_0_yyyy_0_xxxy_1, g_0_yyyy_0_xxxyy_0, g_0_yyyy_0_xxxyy_1, g_0_yyyy_0_xxxyz_0, g_0_yyyy_0_xxxyz_1, g_0_yyyy_0_xxxz_1, g_0_yyyy_0_xxxzz_0, g_0_yyyy_0_xxxzz_1, g_0_yyyy_0_xxyy_1, g_0_yyyy_0_xxyyy_0, g_0_yyyy_0_xxyyy_1, g_0_yyyy_0_xxyyz_0, g_0_yyyy_0_xxyyz_1, g_0_yyyy_0_xxyz_1, g_0_yyyy_0_xxyzz_0, g_0_yyyy_0_xxyzz_1, g_0_yyyy_0_xxzz_1, g_0_yyyy_0_xxzzz_0, g_0_yyyy_0_xxzzz_1, g_0_yyyy_0_xyyy_1, g_0_yyyy_0_xyyyy_0, g_0_yyyy_0_xyyyy_1, g_0_yyyy_0_xyyyz_0, g_0_yyyy_0_xyyyz_1, g_0_yyyy_0_xyyz_1, g_0_yyyy_0_xyyzz_0, g_0_yyyy_0_xyyzz_1, g_0_yyyy_0_xyzz_1, g_0_yyyy_0_xyzzz_0, g_0_yyyy_0_xyzzz_1, g_0_yyyy_0_xzzz_1, g_0_yyyy_0_xzzzz_0, g_0_yyyy_0_xzzzz_1, g_0_yyyy_0_yyyy_1, g_0_yyyy_0_yyyyy_0, g_0_yyyy_0_yyyyy_1, g_0_yyyy_0_yyyyz_0, g_0_yyyy_0_yyyyz_1, g_0_yyyy_0_yyyz_1, g_0_yyyy_0_yyyzz_0, g_0_yyyy_0_yyyzz_1, g_0_yyyy_0_yyzz_1, g_0_yyyy_0_yyzzz_0, g_0_yyyy_0_yyzzz_1, g_0_yyyy_0_yzzz_1, g_0_yyyy_0_yzzzz_0, g_0_yyyy_0_yzzzz_1, g_0_yyyy_0_zzzz_1, g_0_yyyy_0_zzzzz_0, g_0_yyyy_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyy_0_xxxxx_0[i] = 5.0 * g_0_yyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxx_0[i] * pb_x + g_0_yyyy_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxy_0[i] = 4.0 * g_0_yyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxy_0[i] * pb_x + g_0_yyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxz_0[i] = 4.0 * g_0_yyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxz_0[i] * pb_x + g_0_yyyy_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyy_0[i] = 3.0 * g_0_yyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyy_0[i] * pb_x + g_0_yyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyz_0[i] = 3.0 * g_0_yyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyz_0[i] * pb_x + g_0_yyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxzz_0[i] = 3.0 * g_0_yyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxzz_0[i] * pb_x + g_0_yyyy_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyy_0[i] = 2.0 * g_0_yyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyy_0[i] * pb_x + g_0_yyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyz_0[i] = 2.0 * g_0_yyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyz_0[i] * pb_x + g_0_yyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyzz_0[i] = 2.0 * g_0_yyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzz_0[i] * pb_x + g_0_yyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxzzz_0[i] = 2.0 * g_0_yyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzzz_0[i] * pb_x + g_0_yyyy_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyy_0[i] = g_0_yyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyy_0[i] * pb_x + g_0_yyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyz_0[i] = g_0_yyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyz_0[i] * pb_x + g_0_yyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyzz_0[i] = g_0_yyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzz_0[i] * pb_x + g_0_yyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyzzz_0[i] = g_0_yyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzz_0[i] * pb_x + g_0_yyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xzzzz_0[i] = g_0_yyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzzz_0[i] * pb_x + g_0_yyyy_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyy_0[i] = g_0_yyyy_0_yyyyy_0[i] * pb_x + g_0_yyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyz_0[i] = g_0_yyyy_0_yyyyz_0[i] * pb_x + g_0_yyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyzz_0[i] = g_0_yyyy_0_yyyzz_0[i] * pb_x + g_0_yyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyzzz_0[i] = g_0_yyyy_0_yyzzz_0[i] * pb_x + g_0_yyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yzzzz_0[i] = g_0_yyyy_0_yzzzz_0[i] * pb_x + g_0_yyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_zzzzz_0[i] = g_0_yyyy_0_zzzzz_0[i] * pb_x + g_0_yyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 231-252 components of targeted buffer : SHSH

    auto g_0_xyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 231);

    auto g_0_xyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 232);

    auto g_0_xyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 233);

    auto g_0_xyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 234);

    auto g_0_xyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 235);

    auto g_0_xyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 236);

    auto g_0_xyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 237);

    auto g_0_xyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 238);

    auto g_0_xyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 239);

    auto g_0_xyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 240);

    auto g_0_xyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 241);

    auto g_0_xyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 242);

    auto g_0_xyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 243);

    auto g_0_xyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 244);

    auto g_0_xyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 245);

    auto g_0_xyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 246);

    auto g_0_xyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 247);

    auto g_0_xyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 248);

    auto g_0_xyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 249);

    auto g_0_xyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 250);

    auto g_0_xyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 251);

    #pragma omp simd aligned(g_0_xyyy_0_xxxxx_0, g_0_xyyy_0_xxxxx_1, g_0_xyyy_0_xxxxy_0, g_0_xyyy_0_xxxxy_1, g_0_xyyy_0_xxxyy_0, g_0_xyyy_0_xxxyy_1, g_0_xyyy_0_xxyyy_0, g_0_xyyy_0_xxyyy_1, g_0_xyyy_0_xyyyy_0, g_0_xyyy_0_xyyyy_1, g_0_xyyyz_0_xxxxx_0, g_0_xyyyz_0_xxxxy_0, g_0_xyyyz_0_xxxxz_0, g_0_xyyyz_0_xxxyy_0, g_0_xyyyz_0_xxxyz_0, g_0_xyyyz_0_xxxzz_0, g_0_xyyyz_0_xxyyy_0, g_0_xyyyz_0_xxyyz_0, g_0_xyyyz_0_xxyzz_0, g_0_xyyyz_0_xxzzz_0, g_0_xyyyz_0_xyyyy_0, g_0_xyyyz_0_xyyyz_0, g_0_xyyyz_0_xyyzz_0, g_0_xyyyz_0_xyzzz_0, g_0_xyyyz_0_xzzzz_0, g_0_xyyyz_0_yyyyy_0, g_0_xyyyz_0_yyyyz_0, g_0_xyyyz_0_yyyzz_0, g_0_xyyyz_0_yyzzz_0, g_0_xyyyz_0_yzzzz_0, g_0_xyyyz_0_zzzzz_0, g_0_yyyz_0_xxxxz_0, g_0_yyyz_0_xxxxz_1, g_0_yyyz_0_xxxyz_0, g_0_yyyz_0_xxxyz_1, g_0_yyyz_0_xxxz_1, g_0_yyyz_0_xxxzz_0, g_0_yyyz_0_xxxzz_1, g_0_yyyz_0_xxyyz_0, g_0_yyyz_0_xxyyz_1, g_0_yyyz_0_xxyz_1, g_0_yyyz_0_xxyzz_0, g_0_yyyz_0_xxyzz_1, g_0_yyyz_0_xxzz_1, g_0_yyyz_0_xxzzz_0, g_0_yyyz_0_xxzzz_1, g_0_yyyz_0_xyyyz_0, g_0_yyyz_0_xyyyz_1, g_0_yyyz_0_xyyz_1, g_0_yyyz_0_xyyzz_0, g_0_yyyz_0_xyyzz_1, g_0_yyyz_0_xyzz_1, g_0_yyyz_0_xyzzz_0, g_0_yyyz_0_xyzzz_1, g_0_yyyz_0_xzzz_1, g_0_yyyz_0_xzzzz_0, g_0_yyyz_0_xzzzz_1, g_0_yyyz_0_yyyyy_0, g_0_yyyz_0_yyyyy_1, g_0_yyyz_0_yyyyz_0, g_0_yyyz_0_yyyyz_1, g_0_yyyz_0_yyyz_1, g_0_yyyz_0_yyyzz_0, g_0_yyyz_0_yyyzz_1, g_0_yyyz_0_yyzz_1, g_0_yyyz_0_yyzzz_0, g_0_yyyz_0_yyzzz_1, g_0_yyyz_0_yzzz_1, g_0_yyyz_0_yzzzz_0, g_0_yyyz_0_yzzzz_1, g_0_yyyz_0_zzzz_1, g_0_yyyz_0_zzzzz_0, g_0_yyyz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyz_0_xxxxx_0[i] = g_0_xyyy_0_xxxxx_0[i] * pb_z + g_0_xyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxy_0[i] = g_0_xyyy_0_xxxxy_0[i] * pb_z + g_0_xyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxz_0[i] = 4.0 * g_0_yyyz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxz_0[i] * pb_x + g_0_yyyz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyy_0[i] = g_0_xyyy_0_xxxyy_0[i] * pb_z + g_0_xyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxyz_0[i] = 3.0 * g_0_yyyz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyz_0[i] * pb_x + g_0_yyyz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxzz_0[i] = 3.0 * g_0_yyyz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxzz_0[i] * pb_x + g_0_yyyz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyy_0[i] = g_0_xyyy_0_xxyyy_0[i] * pb_z + g_0_xyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxyyz_0[i] = 2.0 * g_0_yyyz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyz_0[i] * pb_x + g_0_yyyz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyzz_0[i] = 2.0 * g_0_yyyz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyzz_0[i] * pb_x + g_0_yyyz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxzzz_0[i] = 2.0 * g_0_yyyz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxzzz_0[i] * pb_x + g_0_yyyz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyy_0[i] = g_0_xyyy_0_xyyyy_0[i] * pb_z + g_0_xyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xyyyz_0[i] = g_0_yyyz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyz_0[i] * pb_x + g_0_yyyz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyzz_0[i] = g_0_yyyz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyzz_0[i] * pb_x + g_0_yyyz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyzzz_0[i] = g_0_yyyz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyzzz_0[i] * pb_x + g_0_yyyz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xzzzz_0[i] = g_0_yyyz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xzzzz_0[i] * pb_x + g_0_yyyz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyy_0[i] = g_0_yyyz_0_yyyyy_0[i] * pb_x + g_0_yyyz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyz_0[i] = g_0_yyyz_0_yyyyz_0[i] * pb_x + g_0_yyyz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyzz_0[i] = g_0_yyyz_0_yyyzz_0[i] * pb_x + g_0_yyyz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyzzz_0[i] = g_0_yyyz_0_yyzzz_0[i] * pb_x + g_0_yyyz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yzzzz_0[i] = g_0_yyyz_0_yzzzz_0[i] * pb_x + g_0_yyyz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_zzzzz_0[i] = g_0_yyyz_0_zzzzz_0[i] * pb_x + g_0_yyyz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 252-273 components of targeted buffer : SHSH

    auto g_0_xyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 252);

    auto g_0_xyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 253);

    auto g_0_xyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 254);

    auto g_0_xyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 255);

    auto g_0_xyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 256);

    auto g_0_xyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 257);

    auto g_0_xyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 258);

    auto g_0_xyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 259);

    auto g_0_xyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 260);

    auto g_0_xyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 261);

    auto g_0_xyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 262);

    auto g_0_xyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 263);

    auto g_0_xyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 264);

    auto g_0_xyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 265);

    auto g_0_xyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 266);

    auto g_0_xyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 267);

    auto g_0_xyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 268);

    auto g_0_xyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 269);

    auto g_0_xyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 270);

    auto g_0_xyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 271);

    auto g_0_xyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 272);

    #pragma omp simd aligned(g_0_xyyzz_0_xxxxx_0, g_0_xyyzz_0_xxxxy_0, g_0_xyyzz_0_xxxxz_0, g_0_xyyzz_0_xxxyy_0, g_0_xyyzz_0_xxxyz_0, g_0_xyyzz_0_xxxzz_0, g_0_xyyzz_0_xxyyy_0, g_0_xyyzz_0_xxyyz_0, g_0_xyyzz_0_xxyzz_0, g_0_xyyzz_0_xxzzz_0, g_0_xyyzz_0_xyyyy_0, g_0_xyyzz_0_xyyyz_0, g_0_xyyzz_0_xyyzz_0, g_0_xyyzz_0_xyzzz_0, g_0_xyyzz_0_xzzzz_0, g_0_xyyzz_0_yyyyy_0, g_0_xyyzz_0_yyyyz_0, g_0_xyyzz_0_yyyzz_0, g_0_xyyzz_0_yyzzz_0, g_0_xyyzz_0_yzzzz_0, g_0_xyyzz_0_zzzzz_0, g_0_yyzz_0_xxxx_1, g_0_yyzz_0_xxxxx_0, g_0_yyzz_0_xxxxx_1, g_0_yyzz_0_xxxxy_0, g_0_yyzz_0_xxxxy_1, g_0_yyzz_0_xxxxz_0, g_0_yyzz_0_xxxxz_1, g_0_yyzz_0_xxxy_1, g_0_yyzz_0_xxxyy_0, g_0_yyzz_0_xxxyy_1, g_0_yyzz_0_xxxyz_0, g_0_yyzz_0_xxxyz_1, g_0_yyzz_0_xxxz_1, g_0_yyzz_0_xxxzz_0, g_0_yyzz_0_xxxzz_1, g_0_yyzz_0_xxyy_1, g_0_yyzz_0_xxyyy_0, g_0_yyzz_0_xxyyy_1, g_0_yyzz_0_xxyyz_0, g_0_yyzz_0_xxyyz_1, g_0_yyzz_0_xxyz_1, g_0_yyzz_0_xxyzz_0, g_0_yyzz_0_xxyzz_1, g_0_yyzz_0_xxzz_1, g_0_yyzz_0_xxzzz_0, g_0_yyzz_0_xxzzz_1, g_0_yyzz_0_xyyy_1, g_0_yyzz_0_xyyyy_0, g_0_yyzz_0_xyyyy_1, g_0_yyzz_0_xyyyz_0, g_0_yyzz_0_xyyyz_1, g_0_yyzz_0_xyyz_1, g_0_yyzz_0_xyyzz_0, g_0_yyzz_0_xyyzz_1, g_0_yyzz_0_xyzz_1, g_0_yyzz_0_xyzzz_0, g_0_yyzz_0_xyzzz_1, g_0_yyzz_0_xzzz_1, g_0_yyzz_0_xzzzz_0, g_0_yyzz_0_xzzzz_1, g_0_yyzz_0_yyyy_1, g_0_yyzz_0_yyyyy_0, g_0_yyzz_0_yyyyy_1, g_0_yyzz_0_yyyyz_0, g_0_yyzz_0_yyyyz_1, g_0_yyzz_0_yyyz_1, g_0_yyzz_0_yyyzz_0, g_0_yyzz_0_yyyzz_1, g_0_yyzz_0_yyzz_1, g_0_yyzz_0_yyzzz_0, g_0_yyzz_0_yyzzz_1, g_0_yyzz_0_yzzz_1, g_0_yyzz_0_yzzzz_0, g_0_yyzz_0_yzzzz_1, g_0_yyzz_0_zzzz_1, g_0_yyzz_0_zzzzz_0, g_0_yyzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzz_0_xxxxx_0[i] = 5.0 * g_0_yyzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxx_0[i] * pb_x + g_0_yyzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxy_0[i] = 4.0 * g_0_yyzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxy_0[i] * pb_x + g_0_yyzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxz_0[i] = 4.0 * g_0_yyzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxz_0[i] * pb_x + g_0_yyzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyy_0[i] = 3.0 * g_0_yyzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyy_0[i] * pb_x + g_0_yyzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyz_0[i] = 3.0 * g_0_yyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyz_0[i] * pb_x + g_0_yyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxzz_0[i] = 3.0 * g_0_yyzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxzz_0[i] * pb_x + g_0_yyzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyy_0[i] = 2.0 * g_0_yyzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyy_0[i] * pb_x + g_0_yyzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyz_0[i] = 2.0 * g_0_yyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyz_0[i] * pb_x + g_0_yyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyzz_0[i] = 2.0 * g_0_yyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyzz_0[i] * pb_x + g_0_yyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxzzz_0[i] = 2.0 * g_0_yyzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxzzz_0[i] * pb_x + g_0_yyzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyy_0[i] = g_0_yyzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyy_0[i] * pb_x + g_0_yyzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyz_0[i] = g_0_yyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyz_0[i] * pb_x + g_0_yyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyzz_0[i] = g_0_yyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyzz_0[i] * pb_x + g_0_yyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyzzz_0[i] = g_0_yyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyzzz_0[i] * pb_x + g_0_yyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xzzzz_0[i] = g_0_yyzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xzzzz_0[i] * pb_x + g_0_yyzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyy_0[i] = g_0_yyzz_0_yyyyy_0[i] * pb_x + g_0_yyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyz_0[i] = g_0_yyzz_0_yyyyz_0[i] * pb_x + g_0_yyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyzz_0[i] = g_0_yyzz_0_yyyzz_0[i] * pb_x + g_0_yyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyzzz_0[i] = g_0_yyzz_0_yyzzz_0[i] * pb_x + g_0_yyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yzzzz_0[i] = g_0_yyzz_0_yzzzz_0[i] * pb_x + g_0_yyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_zzzzz_0[i] = g_0_yyzz_0_zzzzz_0[i] * pb_x + g_0_yyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 273-294 components of targeted buffer : SHSH

    auto g_0_xyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 273);

    auto g_0_xyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 274);

    auto g_0_xyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 275);

    auto g_0_xyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 276);

    auto g_0_xyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 277);

    auto g_0_xyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 278);

    auto g_0_xyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 279);

    auto g_0_xyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 280);

    auto g_0_xyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 281);

    auto g_0_xyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 282);

    auto g_0_xyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 283);

    auto g_0_xyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 284);

    auto g_0_xyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 285);

    auto g_0_xyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 286);

    auto g_0_xyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 287);

    auto g_0_xyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 288);

    auto g_0_xyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 289);

    auto g_0_xyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 290);

    auto g_0_xyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 291);

    auto g_0_xyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 292);

    auto g_0_xyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 293);

    #pragma omp simd aligned(g_0_xyzzz_0_xxxxx_0, g_0_xyzzz_0_xxxxy_0, g_0_xyzzz_0_xxxxz_0, g_0_xyzzz_0_xxxyy_0, g_0_xyzzz_0_xxxyz_0, g_0_xyzzz_0_xxxzz_0, g_0_xyzzz_0_xxyyy_0, g_0_xyzzz_0_xxyyz_0, g_0_xyzzz_0_xxyzz_0, g_0_xyzzz_0_xxzzz_0, g_0_xyzzz_0_xyyyy_0, g_0_xyzzz_0_xyyyz_0, g_0_xyzzz_0_xyyzz_0, g_0_xyzzz_0_xyzzz_0, g_0_xyzzz_0_xzzzz_0, g_0_xyzzz_0_yyyyy_0, g_0_xyzzz_0_yyyyz_0, g_0_xyzzz_0_yyyzz_0, g_0_xyzzz_0_yyzzz_0, g_0_xyzzz_0_yzzzz_0, g_0_xyzzz_0_zzzzz_0, g_0_xzzz_0_xxxxx_0, g_0_xzzz_0_xxxxx_1, g_0_xzzz_0_xxxxz_0, g_0_xzzz_0_xxxxz_1, g_0_xzzz_0_xxxzz_0, g_0_xzzz_0_xxxzz_1, g_0_xzzz_0_xxzzz_0, g_0_xzzz_0_xxzzz_1, g_0_xzzz_0_xzzzz_0, g_0_xzzz_0_xzzzz_1, g_0_yzzz_0_xxxxy_0, g_0_yzzz_0_xxxxy_1, g_0_yzzz_0_xxxy_1, g_0_yzzz_0_xxxyy_0, g_0_yzzz_0_xxxyy_1, g_0_yzzz_0_xxxyz_0, g_0_yzzz_0_xxxyz_1, g_0_yzzz_0_xxyy_1, g_0_yzzz_0_xxyyy_0, g_0_yzzz_0_xxyyy_1, g_0_yzzz_0_xxyyz_0, g_0_yzzz_0_xxyyz_1, g_0_yzzz_0_xxyz_1, g_0_yzzz_0_xxyzz_0, g_0_yzzz_0_xxyzz_1, g_0_yzzz_0_xyyy_1, g_0_yzzz_0_xyyyy_0, g_0_yzzz_0_xyyyy_1, g_0_yzzz_0_xyyyz_0, g_0_yzzz_0_xyyyz_1, g_0_yzzz_0_xyyz_1, g_0_yzzz_0_xyyzz_0, g_0_yzzz_0_xyyzz_1, g_0_yzzz_0_xyzz_1, g_0_yzzz_0_xyzzz_0, g_0_yzzz_0_xyzzz_1, g_0_yzzz_0_yyyy_1, g_0_yzzz_0_yyyyy_0, g_0_yzzz_0_yyyyy_1, g_0_yzzz_0_yyyyz_0, g_0_yzzz_0_yyyyz_1, g_0_yzzz_0_yyyz_1, g_0_yzzz_0_yyyzz_0, g_0_yzzz_0_yyyzz_1, g_0_yzzz_0_yyzz_1, g_0_yzzz_0_yyzzz_0, g_0_yzzz_0_yyzzz_1, g_0_yzzz_0_yzzz_1, g_0_yzzz_0_yzzzz_0, g_0_yzzz_0_yzzzz_1, g_0_yzzz_0_zzzzz_0, g_0_yzzz_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzz_0_xxxxx_0[i] = g_0_xzzz_0_xxxxx_0[i] * pb_y + g_0_xzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxy_0[i] = 4.0 * g_0_yzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxy_0[i] * pb_x + g_0_yzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxz_0[i] = g_0_xzzz_0_xxxxz_0[i] * pb_y + g_0_xzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxyy_0[i] = 3.0 * g_0_yzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyy_0[i] * pb_x + g_0_yzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyz_0[i] = 3.0 * g_0_yzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyz_0[i] * pb_x + g_0_yzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxzz_0[i] = g_0_xzzz_0_xxxzz_0[i] * pb_y + g_0_xzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxyyy_0[i] = 2.0 * g_0_yzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyy_0[i] * pb_x + g_0_yzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyz_0[i] = 2.0 * g_0_yzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyz_0[i] * pb_x + g_0_yzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyzz_0[i] = 2.0 * g_0_yzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyzz_0[i] * pb_x + g_0_yzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxzzz_0[i] = g_0_xzzz_0_xxzzz_0[i] * pb_y + g_0_xzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xyyyy_0[i] = g_0_yzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyy_0[i] * pb_x + g_0_yzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyz_0[i] = g_0_yzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyz_0[i] * pb_x + g_0_yzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyzz_0[i] = g_0_yzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyzz_0[i] * pb_x + g_0_yzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyzzz_0[i] = g_0_yzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyzzz_0[i] * pb_x + g_0_yzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xzzzz_0[i] = g_0_xzzz_0_xzzzz_0[i] * pb_y + g_0_xzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_yyyyy_0[i] = g_0_yzzz_0_yyyyy_0[i] * pb_x + g_0_yzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyz_0[i] = g_0_yzzz_0_yyyyz_0[i] * pb_x + g_0_yzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyzz_0[i] = g_0_yzzz_0_yyyzz_0[i] * pb_x + g_0_yzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyzzz_0[i] = g_0_yzzz_0_yyzzz_0[i] * pb_x + g_0_yzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yzzzz_0[i] = g_0_yzzz_0_yzzzz_0[i] * pb_x + g_0_yzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_zzzzz_0[i] = g_0_yzzz_0_zzzzz_0[i] * pb_x + g_0_yzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 294-315 components of targeted buffer : SHSH

    auto g_0_xzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 294);

    auto g_0_xzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 295);

    auto g_0_xzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 296);

    auto g_0_xzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 297);

    auto g_0_xzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 298);

    auto g_0_xzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 299);

    auto g_0_xzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 300);

    auto g_0_xzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 301);

    auto g_0_xzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 302);

    auto g_0_xzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 303);

    auto g_0_xzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 304);

    auto g_0_xzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 305);

    auto g_0_xzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 306);

    auto g_0_xzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 307);

    auto g_0_xzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 308);

    auto g_0_xzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 309);

    auto g_0_xzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 310);

    auto g_0_xzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 311);

    auto g_0_xzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 312);

    auto g_0_xzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 313);

    auto g_0_xzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 314);

    #pragma omp simd aligned(g_0_xzzzz_0_xxxxx_0, g_0_xzzzz_0_xxxxy_0, g_0_xzzzz_0_xxxxz_0, g_0_xzzzz_0_xxxyy_0, g_0_xzzzz_0_xxxyz_0, g_0_xzzzz_0_xxxzz_0, g_0_xzzzz_0_xxyyy_0, g_0_xzzzz_0_xxyyz_0, g_0_xzzzz_0_xxyzz_0, g_0_xzzzz_0_xxzzz_0, g_0_xzzzz_0_xyyyy_0, g_0_xzzzz_0_xyyyz_0, g_0_xzzzz_0_xyyzz_0, g_0_xzzzz_0_xyzzz_0, g_0_xzzzz_0_xzzzz_0, g_0_xzzzz_0_yyyyy_0, g_0_xzzzz_0_yyyyz_0, g_0_xzzzz_0_yyyzz_0, g_0_xzzzz_0_yyzzz_0, g_0_xzzzz_0_yzzzz_0, g_0_xzzzz_0_zzzzz_0, g_0_zzzz_0_xxxx_1, g_0_zzzz_0_xxxxx_0, g_0_zzzz_0_xxxxx_1, g_0_zzzz_0_xxxxy_0, g_0_zzzz_0_xxxxy_1, g_0_zzzz_0_xxxxz_0, g_0_zzzz_0_xxxxz_1, g_0_zzzz_0_xxxy_1, g_0_zzzz_0_xxxyy_0, g_0_zzzz_0_xxxyy_1, g_0_zzzz_0_xxxyz_0, g_0_zzzz_0_xxxyz_1, g_0_zzzz_0_xxxz_1, g_0_zzzz_0_xxxzz_0, g_0_zzzz_0_xxxzz_1, g_0_zzzz_0_xxyy_1, g_0_zzzz_0_xxyyy_0, g_0_zzzz_0_xxyyy_1, g_0_zzzz_0_xxyyz_0, g_0_zzzz_0_xxyyz_1, g_0_zzzz_0_xxyz_1, g_0_zzzz_0_xxyzz_0, g_0_zzzz_0_xxyzz_1, g_0_zzzz_0_xxzz_1, g_0_zzzz_0_xxzzz_0, g_0_zzzz_0_xxzzz_1, g_0_zzzz_0_xyyy_1, g_0_zzzz_0_xyyyy_0, g_0_zzzz_0_xyyyy_1, g_0_zzzz_0_xyyyz_0, g_0_zzzz_0_xyyyz_1, g_0_zzzz_0_xyyz_1, g_0_zzzz_0_xyyzz_0, g_0_zzzz_0_xyyzz_1, g_0_zzzz_0_xyzz_1, g_0_zzzz_0_xyzzz_0, g_0_zzzz_0_xyzzz_1, g_0_zzzz_0_xzzz_1, g_0_zzzz_0_xzzzz_0, g_0_zzzz_0_xzzzz_1, g_0_zzzz_0_yyyy_1, g_0_zzzz_0_yyyyy_0, g_0_zzzz_0_yyyyy_1, g_0_zzzz_0_yyyyz_0, g_0_zzzz_0_yyyyz_1, g_0_zzzz_0_yyyz_1, g_0_zzzz_0_yyyzz_0, g_0_zzzz_0_yyyzz_1, g_0_zzzz_0_yyzz_1, g_0_zzzz_0_yyzzz_0, g_0_zzzz_0_yyzzz_1, g_0_zzzz_0_yzzz_1, g_0_zzzz_0_yzzzz_0, g_0_zzzz_0_yzzzz_1, g_0_zzzz_0_zzzz_1, g_0_zzzz_0_zzzzz_0, g_0_zzzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzz_0_xxxxx_0[i] = 5.0 * g_0_zzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxx_0[i] * pb_x + g_0_zzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxy_0[i] = 4.0 * g_0_zzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxy_0[i] * pb_x + g_0_zzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxz_0[i] = 4.0 * g_0_zzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxz_0[i] * pb_x + g_0_zzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyy_0[i] = 3.0 * g_0_zzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyy_0[i] * pb_x + g_0_zzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyz_0[i] = 3.0 * g_0_zzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyz_0[i] * pb_x + g_0_zzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxzz_0[i] = 3.0 * g_0_zzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxzz_0[i] * pb_x + g_0_zzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyy_0[i] = 2.0 * g_0_zzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyy_0[i] * pb_x + g_0_zzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyz_0[i] = 2.0 * g_0_zzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyz_0[i] * pb_x + g_0_zzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyzz_0[i] = 2.0 * g_0_zzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzz_0[i] * pb_x + g_0_zzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxzzz_0[i] = 2.0 * g_0_zzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzzz_0[i] * pb_x + g_0_zzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyy_0[i] = g_0_zzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyy_0[i] * pb_x + g_0_zzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyz_0[i] = g_0_zzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyz_0[i] * pb_x + g_0_zzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyzz_0[i] = g_0_zzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzz_0[i] * pb_x + g_0_zzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyzzz_0[i] = g_0_zzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzz_0[i] * pb_x + g_0_zzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xzzzz_0[i] = g_0_zzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzzz_0[i] * pb_x + g_0_zzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyy_0[i] = g_0_zzzz_0_yyyyy_0[i] * pb_x + g_0_zzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyz_0[i] = g_0_zzzz_0_yyyyz_0[i] * pb_x + g_0_zzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyzz_0[i] = g_0_zzzz_0_yyyzz_0[i] * pb_x + g_0_zzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyzzz_0[i] = g_0_zzzz_0_yyzzz_0[i] * pb_x + g_0_zzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yzzzz_0[i] = g_0_zzzz_0_yzzzz_0[i] * pb_x + g_0_zzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_zzzzz_0[i] = g_0_zzzz_0_zzzzz_0[i] * pb_x + g_0_zzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 315-336 components of targeted buffer : SHSH

    auto g_0_yyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 315);

    auto g_0_yyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 316);

    auto g_0_yyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 317);

    auto g_0_yyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 318);

    auto g_0_yyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 319);

    auto g_0_yyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 320);

    auto g_0_yyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 321);

    auto g_0_yyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 322);

    auto g_0_yyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 323);

    auto g_0_yyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 324);

    auto g_0_yyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 325);

    auto g_0_yyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 326);

    auto g_0_yyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 327);

    auto g_0_yyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 328);

    auto g_0_yyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 329);

    auto g_0_yyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 330);

    auto g_0_yyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 331);

    auto g_0_yyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 332);

    auto g_0_yyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 333);

    auto g_0_yyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 334);

    auto g_0_yyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 335);

    #pragma omp simd aligned(g_0_yyy_0_xxxxx_0, g_0_yyy_0_xxxxx_1, g_0_yyy_0_xxxxy_0, g_0_yyy_0_xxxxy_1, g_0_yyy_0_xxxxz_0, g_0_yyy_0_xxxxz_1, g_0_yyy_0_xxxyy_0, g_0_yyy_0_xxxyy_1, g_0_yyy_0_xxxyz_0, g_0_yyy_0_xxxyz_1, g_0_yyy_0_xxxzz_0, g_0_yyy_0_xxxzz_1, g_0_yyy_0_xxyyy_0, g_0_yyy_0_xxyyy_1, g_0_yyy_0_xxyyz_0, g_0_yyy_0_xxyyz_1, g_0_yyy_0_xxyzz_0, g_0_yyy_0_xxyzz_1, g_0_yyy_0_xxzzz_0, g_0_yyy_0_xxzzz_1, g_0_yyy_0_xyyyy_0, g_0_yyy_0_xyyyy_1, g_0_yyy_0_xyyyz_0, g_0_yyy_0_xyyyz_1, g_0_yyy_0_xyyzz_0, g_0_yyy_0_xyyzz_1, g_0_yyy_0_xyzzz_0, g_0_yyy_0_xyzzz_1, g_0_yyy_0_xzzzz_0, g_0_yyy_0_xzzzz_1, g_0_yyy_0_yyyyy_0, g_0_yyy_0_yyyyy_1, g_0_yyy_0_yyyyz_0, g_0_yyy_0_yyyyz_1, g_0_yyy_0_yyyzz_0, g_0_yyy_0_yyyzz_1, g_0_yyy_0_yyzzz_0, g_0_yyy_0_yyzzz_1, g_0_yyy_0_yzzzz_0, g_0_yyy_0_yzzzz_1, g_0_yyy_0_zzzzz_0, g_0_yyy_0_zzzzz_1, g_0_yyyy_0_xxxx_1, g_0_yyyy_0_xxxxx_0, g_0_yyyy_0_xxxxx_1, g_0_yyyy_0_xxxxy_0, g_0_yyyy_0_xxxxy_1, g_0_yyyy_0_xxxxz_0, g_0_yyyy_0_xxxxz_1, g_0_yyyy_0_xxxy_1, g_0_yyyy_0_xxxyy_0, g_0_yyyy_0_xxxyy_1, g_0_yyyy_0_xxxyz_0, g_0_yyyy_0_xxxyz_1, g_0_yyyy_0_xxxz_1, g_0_yyyy_0_xxxzz_0, g_0_yyyy_0_xxxzz_1, g_0_yyyy_0_xxyy_1, g_0_yyyy_0_xxyyy_0, g_0_yyyy_0_xxyyy_1, g_0_yyyy_0_xxyyz_0, g_0_yyyy_0_xxyyz_1, g_0_yyyy_0_xxyz_1, g_0_yyyy_0_xxyzz_0, g_0_yyyy_0_xxyzz_1, g_0_yyyy_0_xxzz_1, g_0_yyyy_0_xxzzz_0, g_0_yyyy_0_xxzzz_1, g_0_yyyy_0_xyyy_1, g_0_yyyy_0_xyyyy_0, g_0_yyyy_0_xyyyy_1, g_0_yyyy_0_xyyyz_0, g_0_yyyy_0_xyyyz_1, g_0_yyyy_0_xyyz_1, g_0_yyyy_0_xyyzz_0, g_0_yyyy_0_xyyzz_1, g_0_yyyy_0_xyzz_1, g_0_yyyy_0_xyzzz_0, g_0_yyyy_0_xyzzz_1, g_0_yyyy_0_xzzz_1, g_0_yyyy_0_xzzzz_0, g_0_yyyy_0_xzzzz_1, g_0_yyyy_0_yyyy_1, g_0_yyyy_0_yyyyy_0, g_0_yyyy_0_yyyyy_1, g_0_yyyy_0_yyyyz_0, g_0_yyyy_0_yyyyz_1, g_0_yyyy_0_yyyz_1, g_0_yyyy_0_yyyzz_0, g_0_yyyy_0_yyyzz_1, g_0_yyyy_0_yyzz_1, g_0_yyyy_0_yyzzz_0, g_0_yyyy_0_yyzzz_1, g_0_yyyy_0_yzzz_1, g_0_yyyy_0_yzzzz_0, g_0_yyyy_0_yzzzz_1, g_0_yyyy_0_zzzz_1, g_0_yyyy_0_zzzzz_0, g_0_yyyy_0_zzzzz_1, g_0_yyyyy_0_xxxxx_0, g_0_yyyyy_0_xxxxy_0, g_0_yyyyy_0_xxxxz_0, g_0_yyyyy_0_xxxyy_0, g_0_yyyyy_0_xxxyz_0, g_0_yyyyy_0_xxxzz_0, g_0_yyyyy_0_xxyyy_0, g_0_yyyyy_0_xxyyz_0, g_0_yyyyy_0_xxyzz_0, g_0_yyyyy_0_xxzzz_0, g_0_yyyyy_0_xyyyy_0, g_0_yyyyy_0_xyyyz_0, g_0_yyyyy_0_xyyzz_0, g_0_yyyyy_0_xyzzz_0, g_0_yyyyy_0_xzzzz_0, g_0_yyyyy_0_yyyyy_0, g_0_yyyyy_0_yyyyz_0, g_0_yyyyy_0_yyyzz_0, g_0_yyyyy_0_yyzzz_0, g_0_yyyyy_0_yzzzz_0, g_0_yyyyy_0_zzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyy_0_xxxxx_0[i] = 4.0 * g_0_yyy_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxx_0[i] * pb_y + g_0_yyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxy_0[i] = 4.0 * g_0_yyy_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxy_0[i] * pb_y + g_0_yyyy_0_xxxxy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxz_0[i] = 4.0 * g_0_yyy_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxz_0[i] * pb_y + g_0_yyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyy_0[i] = 4.0 * g_0_yyy_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyy_0[i] * pb_y + g_0_yyyy_0_xxxyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyz_0[i] = 4.0 * g_0_yyy_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyz_0[i] * pb_y + g_0_yyyy_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxzz_0[i] = 4.0 * g_0_yyy_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxzz_0[i] * pb_y + g_0_yyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyy_0[i] = 4.0 * g_0_yyy_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyy_0[i] * pb_y + g_0_yyyy_0_xxyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyz_0[i] = 4.0 * g_0_yyy_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyz_0[i] * pb_y + g_0_yyyy_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyzz_0[i] = 4.0 * g_0_yyy_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzz_0[i] * pb_y + g_0_yyyy_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxzzz_0[i] = 4.0 * g_0_yyy_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzzz_0[i] * pb_y + g_0_yyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyy_0[i] = 4.0 * g_0_yyy_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyy_0[i] * pb_y + g_0_yyyy_0_xyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyz_0[i] = 4.0 * g_0_yyy_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyz_0[i] * pb_y + g_0_yyyy_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyzz_0[i] = 4.0 * g_0_yyy_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzz_0[i] * pb_y + g_0_yyyy_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyzzz_0[i] = 4.0 * g_0_yyy_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzz_0[i] * pb_y + g_0_yyyy_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xzzzz_0[i] = 4.0 * g_0_yyy_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzzz_0[i] * pb_y + g_0_yyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyy_0[i] = 4.0 * g_0_yyy_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyy_0[i] * pb_y + g_0_yyyy_0_yyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyz_0[i] = 4.0 * g_0_yyy_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyz_0[i] * pb_y + g_0_yyyy_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyzz_0[i] = 4.0 * g_0_yyy_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyzz_0[i] * pb_y + g_0_yyyy_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyzzz_0[i] = 4.0 * g_0_yyy_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzzz_0[i] * pb_y + g_0_yyyy_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yzzzz_0[i] = 4.0 * g_0_yyy_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzzzz_0[i] * pb_y + g_0_yyyy_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_zzzzz_0[i] = 4.0 * g_0_yyy_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzzz_0[i] * pb_y + g_0_yyyy_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 336-357 components of targeted buffer : SHSH

    auto g_0_yyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 336);

    auto g_0_yyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 337);

    auto g_0_yyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 338);

    auto g_0_yyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 339);

    auto g_0_yyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 340);

    auto g_0_yyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 341);

    auto g_0_yyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 342);

    auto g_0_yyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 343);

    auto g_0_yyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 344);

    auto g_0_yyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 345);

    auto g_0_yyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 346);

    auto g_0_yyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 347);

    auto g_0_yyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 348);

    auto g_0_yyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 349);

    auto g_0_yyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 350);

    auto g_0_yyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 351);

    auto g_0_yyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 352);

    auto g_0_yyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 353);

    auto g_0_yyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 354);

    auto g_0_yyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 355);

    auto g_0_yyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 356);

    #pragma omp simd aligned(g_0_yyyy_0_xxxx_1, g_0_yyyy_0_xxxxx_0, g_0_yyyy_0_xxxxx_1, g_0_yyyy_0_xxxxy_0, g_0_yyyy_0_xxxxy_1, g_0_yyyy_0_xxxxz_0, g_0_yyyy_0_xxxxz_1, g_0_yyyy_0_xxxy_1, g_0_yyyy_0_xxxyy_0, g_0_yyyy_0_xxxyy_1, g_0_yyyy_0_xxxyz_0, g_0_yyyy_0_xxxyz_1, g_0_yyyy_0_xxxz_1, g_0_yyyy_0_xxxzz_0, g_0_yyyy_0_xxxzz_1, g_0_yyyy_0_xxyy_1, g_0_yyyy_0_xxyyy_0, g_0_yyyy_0_xxyyy_1, g_0_yyyy_0_xxyyz_0, g_0_yyyy_0_xxyyz_1, g_0_yyyy_0_xxyz_1, g_0_yyyy_0_xxyzz_0, g_0_yyyy_0_xxyzz_1, g_0_yyyy_0_xxzz_1, g_0_yyyy_0_xxzzz_0, g_0_yyyy_0_xxzzz_1, g_0_yyyy_0_xyyy_1, g_0_yyyy_0_xyyyy_0, g_0_yyyy_0_xyyyy_1, g_0_yyyy_0_xyyyz_0, g_0_yyyy_0_xyyyz_1, g_0_yyyy_0_xyyz_1, g_0_yyyy_0_xyyzz_0, g_0_yyyy_0_xyyzz_1, g_0_yyyy_0_xyzz_1, g_0_yyyy_0_xyzzz_0, g_0_yyyy_0_xyzzz_1, g_0_yyyy_0_xzzz_1, g_0_yyyy_0_xzzzz_0, g_0_yyyy_0_xzzzz_1, g_0_yyyy_0_yyyy_1, g_0_yyyy_0_yyyyy_0, g_0_yyyy_0_yyyyy_1, g_0_yyyy_0_yyyyz_0, g_0_yyyy_0_yyyyz_1, g_0_yyyy_0_yyyz_1, g_0_yyyy_0_yyyzz_0, g_0_yyyy_0_yyyzz_1, g_0_yyyy_0_yyzz_1, g_0_yyyy_0_yyzzz_0, g_0_yyyy_0_yyzzz_1, g_0_yyyy_0_yzzz_1, g_0_yyyy_0_yzzzz_0, g_0_yyyy_0_yzzzz_1, g_0_yyyy_0_zzzz_1, g_0_yyyy_0_zzzzz_0, g_0_yyyy_0_zzzzz_1, g_0_yyyyz_0_xxxxx_0, g_0_yyyyz_0_xxxxy_0, g_0_yyyyz_0_xxxxz_0, g_0_yyyyz_0_xxxyy_0, g_0_yyyyz_0_xxxyz_0, g_0_yyyyz_0_xxxzz_0, g_0_yyyyz_0_xxyyy_0, g_0_yyyyz_0_xxyyz_0, g_0_yyyyz_0_xxyzz_0, g_0_yyyyz_0_xxzzz_0, g_0_yyyyz_0_xyyyy_0, g_0_yyyyz_0_xyyyz_0, g_0_yyyyz_0_xyyzz_0, g_0_yyyyz_0_xyzzz_0, g_0_yyyyz_0_xzzzz_0, g_0_yyyyz_0_yyyyy_0, g_0_yyyyz_0_yyyyz_0, g_0_yyyyz_0_yyyzz_0, g_0_yyyyz_0_yyzzz_0, g_0_yyyyz_0_yzzzz_0, g_0_yyyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyz_0_xxxxx_0[i] = g_0_yyyy_0_xxxxx_0[i] * pb_z + g_0_yyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxy_0[i] = g_0_yyyy_0_xxxxy_0[i] * pb_z + g_0_yyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxz_0[i] = g_0_yyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxz_0[i] * pb_z + g_0_yyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyy_0[i] = g_0_yyyy_0_xxxyy_0[i] * pb_z + g_0_yyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyz_0[i] = g_0_yyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyz_0[i] * pb_z + g_0_yyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxzz_0[i] = 2.0 * g_0_yyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxzz_0[i] * pb_z + g_0_yyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyy_0[i] = g_0_yyyy_0_xxyyy_0[i] * pb_z + g_0_yyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyz_0[i] = g_0_yyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyz_0[i] * pb_z + g_0_yyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyzz_0[i] = 2.0 * g_0_yyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzz_0[i] * pb_z + g_0_yyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxzzz_0[i] = 3.0 * g_0_yyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzzz_0[i] * pb_z + g_0_yyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyy_0[i] = g_0_yyyy_0_xyyyy_0[i] * pb_z + g_0_yyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyz_0[i] = g_0_yyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyz_0[i] * pb_z + g_0_yyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyzz_0[i] = 2.0 * g_0_yyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzz_0[i] * pb_z + g_0_yyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyzzz_0[i] = 3.0 * g_0_yyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzz_0[i] * pb_z + g_0_yyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xzzzz_0[i] = 4.0 * g_0_yyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzzz_0[i] * pb_z + g_0_yyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyy_0[i] = g_0_yyyy_0_yyyyy_0[i] * pb_z + g_0_yyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyz_0[i] = g_0_yyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyz_0[i] * pb_z + g_0_yyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyzz_0[i] = 2.0 * g_0_yyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyzz_0[i] * pb_z + g_0_yyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyzzz_0[i] = 3.0 * g_0_yyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzzz_0[i] * pb_z + g_0_yyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yzzzz_0[i] = 4.0 * g_0_yyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzzzz_0[i] * pb_z + g_0_yyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_zzzzz_0[i] = 5.0 * g_0_yyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_zzzzz_0[i] * pb_z + g_0_yyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 357-378 components of targeted buffer : SHSH

    auto g_0_yyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 357);

    auto g_0_yyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 358);

    auto g_0_yyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 359);

    auto g_0_yyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 360);

    auto g_0_yyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 361);

    auto g_0_yyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 362);

    auto g_0_yyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 363);

    auto g_0_yyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 364);

    auto g_0_yyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 365);

    auto g_0_yyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 366);

    auto g_0_yyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 367);

    auto g_0_yyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 368);

    auto g_0_yyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 369);

    auto g_0_yyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 370);

    auto g_0_yyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 371);

    auto g_0_yyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 372);

    auto g_0_yyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 373);

    auto g_0_yyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 374);

    auto g_0_yyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 375);

    auto g_0_yyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 376);

    auto g_0_yyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 377);

    #pragma omp simd aligned(g_0_yyy_0_xxxxy_0, g_0_yyy_0_xxxxy_1, g_0_yyy_0_xxxyy_0, g_0_yyy_0_xxxyy_1, g_0_yyy_0_xxyyy_0, g_0_yyy_0_xxyyy_1, g_0_yyy_0_xyyyy_0, g_0_yyy_0_xyyyy_1, g_0_yyy_0_yyyyy_0, g_0_yyy_0_yyyyy_1, g_0_yyyz_0_xxxxy_0, g_0_yyyz_0_xxxxy_1, g_0_yyyz_0_xxxyy_0, g_0_yyyz_0_xxxyy_1, g_0_yyyz_0_xxyyy_0, g_0_yyyz_0_xxyyy_1, g_0_yyyz_0_xyyyy_0, g_0_yyyz_0_xyyyy_1, g_0_yyyz_0_yyyyy_0, g_0_yyyz_0_yyyyy_1, g_0_yyyzz_0_xxxxx_0, g_0_yyyzz_0_xxxxy_0, g_0_yyyzz_0_xxxxz_0, g_0_yyyzz_0_xxxyy_0, g_0_yyyzz_0_xxxyz_0, g_0_yyyzz_0_xxxzz_0, g_0_yyyzz_0_xxyyy_0, g_0_yyyzz_0_xxyyz_0, g_0_yyyzz_0_xxyzz_0, g_0_yyyzz_0_xxzzz_0, g_0_yyyzz_0_xyyyy_0, g_0_yyyzz_0_xyyyz_0, g_0_yyyzz_0_xyyzz_0, g_0_yyyzz_0_xyzzz_0, g_0_yyyzz_0_xzzzz_0, g_0_yyyzz_0_yyyyy_0, g_0_yyyzz_0_yyyyz_0, g_0_yyyzz_0_yyyzz_0, g_0_yyyzz_0_yyzzz_0, g_0_yyyzz_0_yzzzz_0, g_0_yyyzz_0_zzzzz_0, g_0_yyzz_0_xxxxx_0, g_0_yyzz_0_xxxxx_1, g_0_yyzz_0_xxxxz_0, g_0_yyzz_0_xxxxz_1, g_0_yyzz_0_xxxyz_0, g_0_yyzz_0_xxxyz_1, g_0_yyzz_0_xxxz_1, g_0_yyzz_0_xxxzz_0, g_0_yyzz_0_xxxzz_1, g_0_yyzz_0_xxyyz_0, g_0_yyzz_0_xxyyz_1, g_0_yyzz_0_xxyz_1, g_0_yyzz_0_xxyzz_0, g_0_yyzz_0_xxyzz_1, g_0_yyzz_0_xxzz_1, g_0_yyzz_0_xxzzz_0, g_0_yyzz_0_xxzzz_1, g_0_yyzz_0_xyyyz_0, g_0_yyzz_0_xyyyz_1, g_0_yyzz_0_xyyz_1, g_0_yyzz_0_xyyzz_0, g_0_yyzz_0_xyyzz_1, g_0_yyzz_0_xyzz_1, g_0_yyzz_0_xyzzz_0, g_0_yyzz_0_xyzzz_1, g_0_yyzz_0_xzzz_1, g_0_yyzz_0_xzzzz_0, g_0_yyzz_0_xzzzz_1, g_0_yyzz_0_yyyyz_0, g_0_yyzz_0_yyyyz_1, g_0_yyzz_0_yyyz_1, g_0_yyzz_0_yyyzz_0, g_0_yyzz_0_yyyzz_1, g_0_yyzz_0_yyzz_1, g_0_yyzz_0_yyzzz_0, g_0_yyzz_0_yyzzz_1, g_0_yyzz_0_yzzz_1, g_0_yyzz_0_yzzzz_0, g_0_yyzz_0_yzzzz_1, g_0_yyzz_0_zzzz_1, g_0_yyzz_0_zzzzz_0, g_0_yyzz_0_zzzzz_1, g_0_yzz_0_xxxxx_0, g_0_yzz_0_xxxxx_1, g_0_yzz_0_xxxxz_0, g_0_yzz_0_xxxxz_1, g_0_yzz_0_xxxyz_0, g_0_yzz_0_xxxyz_1, g_0_yzz_0_xxxzz_0, g_0_yzz_0_xxxzz_1, g_0_yzz_0_xxyyz_0, g_0_yzz_0_xxyyz_1, g_0_yzz_0_xxyzz_0, g_0_yzz_0_xxyzz_1, g_0_yzz_0_xxzzz_0, g_0_yzz_0_xxzzz_1, g_0_yzz_0_xyyyz_0, g_0_yzz_0_xyyyz_1, g_0_yzz_0_xyyzz_0, g_0_yzz_0_xyyzz_1, g_0_yzz_0_xyzzz_0, g_0_yzz_0_xyzzz_1, g_0_yzz_0_xzzzz_0, g_0_yzz_0_xzzzz_1, g_0_yzz_0_yyyyz_0, g_0_yzz_0_yyyyz_1, g_0_yzz_0_yyyzz_0, g_0_yzz_0_yyyzz_1, g_0_yzz_0_yyzzz_0, g_0_yzz_0_yyzzz_1, g_0_yzz_0_yzzzz_0, g_0_yzz_0_yzzzz_1, g_0_yzz_0_zzzzz_0, g_0_yzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzz_0_xxxxx_0[i] = 2.0 * g_0_yzz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxx_0[i] * pb_y + g_0_yyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxy_0[i] = g_0_yyy_0_xxxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxy_0[i] * pb_z + g_0_yyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxz_0[i] = 2.0 * g_0_yzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxz_0[i] * pb_y + g_0_yyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyy_0[i] = g_0_yyy_0_xxxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxyy_0[i] * pb_z + g_0_yyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxyz_0[i] = 2.0 * g_0_yzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyz_0[i] * pb_y + g_0_yyzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxzz_0[i] = 2.0 * g_0_yzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxzz_0[i] * pb_y + g_0_yyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyy_0[i] = g_0_yyy_0_xxyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxyyy_0[i] * pb_z + g_0_yyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxyyz_0[i] = 2.0 * g_0_yzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyz_0[i] * pb_y + g_0_yyzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyzz_0[i] = 2.0 * g_0_yzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyzz_0[i] * pb_y + g_0_yyzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxzzz_0[i] = 2.0 * g_0_yzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzzz_0[i] * pb_y + g_0_yyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyy_0[i] = g_0_yyy_0_xyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xyyyy_0[i] * pb_z + g_0_yyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xyyyz_0[i] = 2.0 * g_0_yzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyz_0[i] * pb_y + g_0_yyzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyzz_0[i] = 2.0 * g_0_yzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyzz_0[i] * pb_y + g_0_yyzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyzzz_0[i] = 2.0 * g_0_yzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyzzz_0[i] * pb_y + g_0_yyzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xzzzz_0[i] = 2.0 * g_0_yzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzzz_0[i] * pb_y + g_0_yyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyy_0[i] = g_0_yyy_0_yyyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_yyyyy_0[i] * pb_z + g_0_yyyz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_yyyyz_0[i] = 2.0 * g_0_yzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyz_0[i] * pb_y + g_0_yyzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyzz_0[i] = 2.0 * g_0_yzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyzz_0[i] * pb_y + g_0_yyzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyzzz_0[i] = 2.0 * g_0_yzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyzzz_0[i] * pb_y + g_0_yyzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yzzzz_0[i] = 2.0 * g_0_yzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yzzzz_0[i] * pb_y + g_0_yyzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_zzzzz_0[i] = 2.0 * g_0_yzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzzz_0[i] * pb_y + g_0_yyzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 378-399 components of targeted buffer : SHSH

    auto g_0_yyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 378);

    auto g_0_yyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 379);

    auto g_0_yyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 380);

    auto g_0_yyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 381);

    auto g_0_yyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 382);

    auto g_0_yyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 383);

    auto g_0_yyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 384);

    auto g_0_yyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 385);

    auto g_0_yyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 386);

    auto g_0_yyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 387);

    auto g_0_yyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 388);

    auto g_0_yyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 389);

    auto g_0_yyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 390);

    auto g_0_yyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 391);

    auto g_0_yyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 392);

    auto g_0_yyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 393);

    auto g_0_yyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 394);

    auto g_0_yyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 395);

    auto g_0_yyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 396);

    auto g_0_yyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 397);

    auto g_0_yyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 398);

    #pragma omp simd aligned(g_0_yyz_0_xxxxy_0, g_0_yyz_0_xxxxy_1, g_0_yyz_0_xxxyy_0, g_0_yyz_0_xxxyy_1, g_0_yyz_0_xxyyy_0, g_0_yyz_0_xxyyy_1, g_0_yyz_0_xyyyy_0, g_0_yyz_0_xyyyy_1, g_0_yyz_0_yyyyy_0, g_0_yyz_0_yyyyy_1, g_0_yyzz_0_xxxxy_0, g_0_yyzz_0_xxxxy_1, g_0_yyzz_0_xxxyy_0, g_0_yyzz_0_xxxyy_1, g_0_yyzz_0_xxyyy_0, g_0_yyzz_0_xxyyy_1, g_0_yyzz_0_xyyyy_0, g_0_yyzz_0_xyyyy_1, g_0_yyzz_0_yyyyy_0, g_0_yyzz_0_yyyyy_1, g_0_yyzzz_0_xxxxx_0, g_0_yyzzz_0_xxxxy_0, g_0_yyzzz_0_xxxxz_0, g_0_yyzzz_0_xxxyy_0, g_0_yyzzz_0_xxxyz_0, g_0_yyzzz_0_xxxzz_0, g_0_yyzzz_0_xxyyy_0, g_0_yyzzz_0_xxyyz_0, g_0_yyzzz_0_xxyzz_0, g_0_yyzzz_0_xxzzz_0, g_0_yyzzz_0_xyyyy_0, g_0_yyzzz_0_xyyyz_0, g_0_yyzzz_0_xyyzz_0, g_0_yyzzz_0_xyzzz_0, g_0_yyzzz_0_xzzzz_0, g_0_yyzzz_0_yyyyy_0, g_0_yyzzz_0_yyyyz_0, g_0_yyzzz_0_yyyzz_0, g_0_yyzzz_0_yyzzz_0, g_0_yyzzz_0_yzzzz_0, g_0_yyzzz_0_zzzzz_0, g_0_yzzz_0_xxxxx_0, g_0_yzzz_0_xxxxx_1, g_0_yzzz_0_xxxxz_0, g_0_yzzz_0_xxxxz_1, g_0_yzzz_0_xxxyz_0, g_0_yzzz_0_xxxyz_1, g_0_yzzz_0_xxxz_1, g_0_yzzz_0_xxxzz_0, g_0_yzzz_0_xxxzz_1, g_0_yzzz_0_xxyyz_0, g_0_yzzz_0_xxyyz_1, g_0_yzzz_0_xxyz_1, g_0_yzzz_0_xxyzz_0, g_0_yzzz_0_xxyzz_1, g_0_yzzz_0_xxzz_1, g_0_yzzz_0_xxzzz_0, g_0_yzzz_0_xxzzz_1, g_0_yzzz_0_xyyyz_0, g_0_yzzz_0_xyyyz_1, g_0_yzzz_0_xyyz_1, g_0_yzzz_0_xyyzz_0, g_0_yzzz_0_xyyzz_1, g_0_yzzz_0_xyzz_1, g_0_yzzz_0_xyzzz_0, g_0_yzzz_0_xyzzz_1, g_0_yzzz_0_xzzz_1, g_0_yzzz_0_xzzzz_0, g_0_yzzz_0_xzzzz_1, g_0_yzzz_0_yyyyz_0, g_0_yzzz_0_yyyyz_1, g_0_yzzz_0_yyyz_1, g_0_yzzz_0_yyyzz_0, g_0_yzzz_0_yyyzz_1, g_0_yzzz_0_yyzz_1, g_0_yzzz_0_yyzzz_0, g_0_yzzz_0_yyzzz_1, g_0_yzzz_0_yzzz_1, g_0_yzzz_0_yzzzz_0, g_0_yzzz_0_yzzzz_1, g_0_yzzz_0_zzzz_1, g_0_yzzz_0_zzzzz_0, g_0_yzzz_0_zzzzz_1, g_0_zzz_0_xxxxx_0, g_0_zzz_0_xxxxx_1, g_0_zzz_0_xxxxz_0, g_0_zzz_0_xxxxz_1, g_0_zzz_0_xxxyz_0, g_0_zzz_0_xxxyz_1, g_0_zzz_0_xxxzz_0, g_0_zzz_0_xxxzz_1, g_0_zzz_0_xxyyz_0, g_0_zzz_0_xxyyz_1, g_0_zzz_0_xxyzz_0, g_0_zzz_0_xxyzz_1, g_0_zzz_0_xxzzz_0, g_0_zzz_0_xxzzz_1, g_0_zzz_0_xyyyz_0, g_0_zzz_0_xyyyz_1, g_0_zzz_0_xyyzz_0, g_0_zzz_0_xyyzz_1, g_0_zzz_0_xyzzz_0, g_0_zzz_0_xyzzz_1, g_0_zzz_0_xzzzz_0, g_0_zzz_0_xzzzz_1, g_0_zzz_0_yyyyz_0, g_0_zzz_0_yyyyz_1, g_0_zzz_0_yyyzz_0, g_0_zzz_0_yyyzz_1, g_0_zzz_0_yyzzz_0, g_0_zzz_0_yyzzz_1, g_0_zzz_0_yzzzz_0, g_0_zzz_0_yzzzz_1, g_0_zzz_0_zzzzz_0, g_0_zzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzz_0_xxxxx_0[i] = g_0_zzz_0_xxxxx_0[i] * fi_ab_0 - g_0_zzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxx_0[i] * pb_y + g_0_yzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxy_0[i] = 2.0 * g_0_yyz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxy_0[i] * pb_z + g_0_yyzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxz_0[i] = g_0_zzz_0_xxxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxz_0[i] * pb_y + g_0_yzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyy_0[i] = 2.0 * g_0_yyz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxyy_0[i] * pb_z + g_0_yyzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxyz_0[i] = g_0_zzz_0_xxxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyz_0[i] * pb_y + g_0_yzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxzz_0[i] = g_0_zzz_0_xxxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxzz_0[i] * pb_y + g_0_yzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyy_0[i] = 2.0 * g_0_yyz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxyyy_0[i] * pb_z + g_0_yyzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxyyz_0[i] = g_0_zzz_0_xxyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyz_0[i] * pb_y + g_0_yzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyzz_0[i] = g_0_zzz_0_xxyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyzz_0[i] * pb_y + g_0_yzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxzzz_0[i] = g_0_zzz_0_xxzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzzz_0[i] * pb_y + g_0_yzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyy_0[i] = 2.0 * g_0_yyz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xyyyy_0[i] * pb_z + g_0_yyzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xyyyz_0[i] = g_0_zzz_0_xyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyz_0[i] * pb_y + g_0_yzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyzz_0[i] = g_0_zzz_0_xyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyzz_0[i] * pb_y + g_0_yzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyzzz_0[i] = g_0_zzz_0_xyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyzzz_0[i] * pb_y + g_0_yzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xzzzz_0[i] = g_0_zzz_0_xzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzzz_0[i] * pb_y + g_0_yzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyy_0[i] = 2.0 * g_0_yyz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_yyyyy_0[i] * pb_z + g_0_yyzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_yyyyz_0[i] = g_0_zzz_0_yyyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyyz_0[i] * pb_y + g_0_yzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyzz_0[i] = g_0_zzz_0_yyyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyzz_0[i] * pb_y + g_0_yzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyzzz_0[i] = g_0_zzz_0_yyzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyzzz_0[i] * pb_y + g_0_yzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yzzzz_0[i] = g_0_zzz_0_yzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yzzzz_0[i] * pb_y + g_0_yzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_zzzzz_0[i] = g_0_zzz_0_zzzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzzz_0[i] * pb_y + g_0_yzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 399-420 components of targeted buffer : SHSH

    auto g_0_yzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 399);

    auto g_0_yzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 400);

    auto g_0_yzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 401);

    auto g_0_yzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 402);

    auto g_0_yzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 403);

    auto g_0_yzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 404);

    auto g_0_yzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 405);

    auto g_0_yzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 406);

    auto g_0_yzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 407);

    auto g_0_yzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 408);

    auto g_0_yzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 409);

    auto g_0_yzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 410);

    auto g_0_yzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 411);

    auto g_0_yzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 412);

    auto g_0_yzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 413);

    auto g_0_yzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 414);

    auto g_0_yzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 415);

    auto g_0_yzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 416);

    auto g_0_yzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 417);

    auto g_0_yzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 418);

    auto g_0_yzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 419);

    #pragma omp simd aligned(g_0_yzzzz_0_xxxxx_0, g_0_yzzzz_0_xxxxy_0, g_0_yzzzz_0_xxxxz_0, g_0_yzzzz_0_xxxyy_0, g_0_yzzzz_0_xxxyz_0, g_0_yzzzz_0_xxxzz_0, g_0_yzzzz_0_xxyyy_0, g_0_yzzzz_0_xxyyz_0, g_0_yzzzz_0_xxyzz_0, g_0_yzzzz_0_xxzzz_0, g_0_yzzzz_0_xyyyy_0, g_0_yzzzz_0_xyyyz_0, g_0_yzzzz_0_xyyzz_0, g_0_yzzzz_0_xyzzz_0, g_0_yzzzz_0_xzzzz_0, g_0_yzzzz_0_yyyyy_0, g_0_yzzzz_0_yyyyz_0, g_0_yzzzz_0_yyyzz_0, g_0_yzzzz_0_yyzzz_0, g_0_yzzzz_0_yzzzz_0, g_0_yzzzz_0_zzzzz_0, g_0_zzzz_0_xxxx_1, g_0_zzzz_0_xxxxx_0, g_0_zzzz_0_xxxxx_1, g_0_zzzz_0_xxxxy_0, g_0_zzzz_0_xxxxy_1, g_0_zzzz_0_xxxxz_0, g_0_zzzz_0_xxxxz_1, g_0_zzzz_0_xxxy_1, g_0_zzzz_0_xxxyy_0, g_0_zzzz_0_xxxyy_1, g_0_zzzz_0_xxxyz_0, g_0_zzzz_0_xxxyz_1, g_0_zzzz_0_xxxz_1, g_0_zzzz_0_xxxzz_0, g_0_zzzz_0_xxxzz_1, g_0_zzzz_0_xxyy_1, g_0_zzzz_0_xxyyy_0, g_0_zzzz_0_xxyyy_1, g_0_zzzz_0_xxyyz_0, g_0_zzzz_0_xxyyz_1, g_0_zzzz_0_xxyz_1, g_0_zzzz_0_xxyzz_0, g_0_zzzz_0_xxyzz_1, g_0_zzzz_0_xxzz_1, g_0_zzzz_0_xxzzz_0, g_0_zzzz_0_xxzzz_1, g_0_zzzz_0_xyyy_1, g_0_zzzz_0_xyyyy_0, g_0_zzzz_0_xyyyy_1, g_0_zzzz_0_xyyyz_0, g_0_zzzz_0_xyyyz_1, g_0_zzzz_0_xyyz_1, g_0_zzzz_0_xyyzz_0, g_0_zzzz_0_xyyzz_1, g_0_zzzz_0_xyzz_1, g_0_zzzz_0_xyzzz_0, g_0_zzzz_0_xyzzz_1, g_0_zzzz_0_xzzz_1, g_0_zzzz_0_xzzzz_0, g_0_zzzz_0_xzzzz_1, g_0_zzzz_0_yyyy_1, g_0_zzzz_0_yyyyy_0, g_0_zzzz_0_yyyyy_1, g_0_zzzz_0_yyyyz_0, g_0_zzzz_0_yyyyz_1, g_0_zzzz_0_yyyz_1, g_0_zzzz_0_yyyzz_0, g_0_zzzz_0_yyyzz_1, g_0_zzzz_0_yyzz_1, g_0_zzzz_0_yyzzz_0, g_0_zzzz_0_yyzzz_1, g_0_zzzz_0_yzzz_1, g_0_zzzz_0_yzzzz_0, g_0_zzzz_0_yzzzz_1, g_0_zzzz_0_zzzz_1, g_0_zzzz_0_zzzzz_0, g_0_zzzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzz_0_xxxxx_0[i] = g_0_zzzz_0_xxxxx_0[i] * pb_y + g_0_zzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxy_0[i] = g_0_zzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxy_0[i] * pb_y + g_0_zzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxz_0[i] = g_0_zzzz_0_xxxxz_0[i] * pb_y + g_0_zzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyy_0[i] = 2.0 * g_0_zzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyy_0[i] * pb_y + g_0_zzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyz_0[i] = g_0_zzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyz_0[i] * pb_y + g_0_zzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxzz_0[i] = g_0_zzzz_0_xxxzz_0[i] * pb_y + g_0_zzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyy_0[i] = 3.0 * g_0_zzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyy_0[i] * pb_y + g_0_zzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyz_0[i] = 2.0 * g_0_zzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyz_0[i] * pb_y + g_0_zzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyzz_0[i] = g_0_zzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzz_0[i] * pb_y + g_0_zzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxzzz_0[i] = g_0_zzzz_0_xxzzz_0[i] * pb_y + g_0_zzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyy_0[i] = 4.0 * g_0_zzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyy_0[i] * pb_y + g_0_zzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyz_0[i] = 3.0 * g_0_zzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyz_0[i] * pb_y + g_0_zzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyzz_0[i] = 2.0 * g_0_zzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzz_0[i] * pb_y + g_0_zzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyzzz_0[i] = g_0_zzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzz_0[i] * pb_y + g_0_zzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xzzzz_0[i] = g_0_zzzz_0_xzzzz_0[i] * pb_y + g_0_zzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyy_0[i] = 5.0 * g_0_zzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyy_0[i] * pb_y + g_0_zzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyz_0[i] = 4.0 * g_0_zzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyz_0[i] * pb_y + g_0_zzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyzz_0[i] = 3.0 * g_0_zzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyzz_0[i] * pb_y + g_0_zzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyzzz_0[i] = 2.0 * g_0_zzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzzz_0[i] * pb_y + g_0_zzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yzzzz_0[i] = g_0_zzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzzz_0[i] * pb_y + g_0_zzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_zzzzz_0[i] = g_0_zzzz_0_zzzzz_0[i] * pb_y + g_0_zzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 420-441 components of targeted buffer : SHSH

    auto g_0_zzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_shsh + 420);

    auto g_0_zzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_shsh + 421);

    auto g_0_zzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_shsh + 422);

    auto g_0_zzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_shsh + 423);

    auto g_0_zzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_shsh + 424);

    auto g_0_zzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_shsh + 425);

    auto g_0_zzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_shsh + 426);

    auto g_0_zzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_shsh + 427);

    auto g_0_zzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_shsh + 428);

    auto g_0_zzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_shsh + 429);

    auto g_0_zzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_shsh + 430);

    auto g_0_zzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_shsh + 431);

    auto g_0_zzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_shsh + 432);

    auto g_0_zzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_shsh + 433);

    auto g_0_zzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_shsh + 434);

    auto g_0_zzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_shsh + 435);

    auto g_0_zzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_shsh + 436);

    auto g_0_zzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_shsh + 437);

    auto g_0_zzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_shsh + 438);

    auto g_0_zzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_shsh + 439);

    auto g_0_zzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_shsh + 440);

    #pragma omp simd aligned(g_0_zzz_0_xxxxx_0, g_0_zzz_0_xxxxx_1, g_0_zzz_0_xxxxy_0, g_0_zzz_0_xxxxy_1, g_0_zzz_0_xxxxz_0, g_0_zzz_0_xxxxz_1, g_0_zzz_0_xxxyy_0, g_0_zzz_0_xxxyy_1, g_0_zzz_0_xxxyz_0, g_0_zzz_0_xxxyz_1, g_0_zzz_0_xxxzz_0, g_0_zzz_0_xxxzz_1, g_0_zzz_0_xxyyy_0, g_0_zzz_0_xxyyy_1, g_0_zzz_0_xxyyz_0, g_0_zzz_0_xxyyz_1, g_0_zzz_0_xxyzz_0, g_0_zzz_0_xxyzz_1, g_0_zzz_0_xxzzz_0, g_0_zzz_0_xxzzz_1, g_0_zzz_0_xyyyy_0, g_0_zzz_0_xyyyy_1, g_0_zzz_0_xyyyz_0, g_0_zzz_0_xyyyz_1, g_0_zzz_0_xyyzz_0, g_0_zzz_0_xyyzz_1, g_0_zzz_0_xyzzz_0, g_0_zzz_0_xyzzz_1, g_0_zzz_0_xzzzz_0, g_0_zzz_0_xzzzz_1, g_0_zzz_0_yyyyy_0, g_0_zzz_0_yyyyy_1, g_0_zzz_0_yyyyz_0, g_0_zzz_0_yyyyz_1, g_0_zzz_0_yyyzz_0, g_0_zzz_0_yyyzz_1, g_0_zzz_0_yyzzz_0, g_0_zzz_0_yyzzz_1, g_0_zzz_0_yzzzz_0, g_0_zzz_0_yzzzz_1, g_0_zzz_0_zzzzz_0, g_0_zzz_0_zzzzz_1, g_0_zzzz_0_xxxx_1, g_0_zzzz_0_xxxxx_0, g_0_zzzz_0_xxxxx_1, g_0_zzzz_0_xxxxy_0, g_0_zzzz_0_xxxxy_1, g_0_zzzz_0_xxxxz_0, g_0_zzzz_0_xxxxz_1, g_0_zzzz_0_xxxy_1, g_0_zzzz_0_xxxyy_0, g_0_zzzz_0_xxxyy_1, g_0_zzzz_0_xxxyz_0, g_0_zzzz_0_xxxyz_1, g_0_zzzz_0_xxxz_1, g_0_zzzz_0_xxxzz_0, g_0_zzzz_0_xxxzz_1, g_0_zzzz_0_xxyy_1, g_0_zzzz_0_xxyyy_0, g_0_zzzz_0_xxyyy_1, g_0_zzzz_0_xxyyz_0, g_0_zzzz_0_xxyyz_1, g_0_zzzz_0_xxyz_1, g_0_zzzz_0_xxyzz_0, g_0_zzzz_0_xxyzz_1, g_0_zzzz_0_xxzz_1, g_0_zzzz_0_xxzzz_0, g_0_zzzz_0_xxzzz_1, g_0_zzzz_0_xyyy_1, g_0_zzzz_0_xyyyy_0, g_0_zzzz_0_xyyyy_1, g_0_zzzz_0_xyyyz_0, g_0_zzzz_0_xyyyz_1, g_0_zzzz_0_xyyz_1, g_0_zzzz_0_xyyzz_0, g_0_zzzz_0_xyyzz_1, g_0_zzzz_0_xyzz_1, g_0_zzzz_0_xyzzz_0, g_0_zzzz_0_xyzzz_1, g_0_zzzz_0_xzzz_1, g_0_zzzz_0_xzzzz_0, g_0_zzzz_0_xzzzz_1, g_0_zzzz_0_yyyy_1, g_0_zzzz_0_yyyyy_0, g_0_zzzz_0_yyyyy_1, g_0_zzzz_0_yyyyz_0, g_0_zzzz_0_yyyyz_1, g_0_zzzz_0_yyyz_1, g_0_zzzz_0_yyyzz_0, g_0_zzzz_0_yyyzz_1, g_0_zzzz_0_yyzz_1, g_0_zzzz_0_yyzzz_0, g_0_zzzz_0_yyzzz_1, g_0_zzzz_0_yzzz_1, g_0_zzzz_0_yzzzz_0, g_0_zzzz_0_yzzzz_1, g_0_zzzz_0_zzzz_1, g_0_zzzz_0_zzzzz_0, g_0_zzzz_0_zzzzz_1, g_0_zzzzz_0_xxxxx_0, g_0_zzzzz_0_xxxxy_0, g_0_zzzzz_0_xxxxz_0, g_0_zzzzz_0_xxxyy_0, g_0_zzzzz_0_xxxyz_0, g_0_zzzzz_0_xxxzz_0, g_0_zzzzz_0_xxyyy_0, g_0_zzzzz_0_xxyyz_0, g_0_zzzzz_0_xxyzz_0, g_0_zzzzz_0_xxzzz_0, g_0_zzzzz_0_xyyyy_0, g_0_zzzzz_0_xyyyz_0, g_0_zzzzz_0_xyyzz_0, g_0_zzzzz_0_xyzzz_0, g_0_zzzzz_0_xzzzz_0, g_0_zzzzz_0_yyyyy_0, g_0_zzzzz_0_yyyyz_0, g_0_zzzzz_0_yyyzz_0, g_0_zzzzz_0_yyzzz_0, g_0_zzzzz_0_yzzzz_0, g_0_zzzzz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzz_0_xxxxx_0[i] = 4.0 * g_0_zzz_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxx_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxx_0[i] * pb_z + g_0_zzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxy_0[i] = 4.0 * g_0_zzz_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxy_0[i] * pb_z + g_0_zzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxz_0[i] = 4.0 * g_0_zzz_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxz_0[i] * pb_z + g_0_zzzz_0_xxxxz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyy_0[i] = 4.0 * g_0_zzz_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxyy_0[i] * pb_z + g_0_zzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyz_0[i] = 4.0 * g_0_zzz_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyz_0[i] * pb_z + g_0_zzzz_0_xxxyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxzz_0[i] = 4.0 * g_0_zzz_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxzz_0[i] * pb_z + g_0_zzzz_0_xxxzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyy_0[i] = 4.0 * g_0_zzz_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxyyy_0[i] * pb_z + g_0_zzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyz_0[i] = 4.0 * g_0_zzz_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyz_0[i] * pb_z + g_0_zzzz_0_xxyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyzz_0[i] = 4.0 * g_0_zzz_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzz_0[i] * pb_z + g_0_zzzz_0_xxyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxzzz_0[i] = 4.0 * g_0_zzz_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzzz_0[i] * pb_z + g_0_zzzz_0_xxzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyy_0[i] = 4.0 * g_0_zzz_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xyyyy_0[i] * pb_z + g_0_zzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyz_0[i] = 4.0 * g_0_zzz_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyz_0[i] * pb_z + g_0_zzzz_0_xyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyzz_0[i] = 4.0 * g_0_zzz_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzz_0[i] * pb_z + g_0_zzzz_0_xyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyzzz_0[i] = 4.0 * g_0_zzz_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzz_0[i] * pb_z + g_0_zzzz_0_xyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xzzzz_0[i] = 4.0 * g_0_zzz_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzzz_0[i] * pb_z + g_0_zzzz_0_xzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyy_0[i] = 4.0 * g_0_zzz_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_yyyyy_0[i] * pb_z + g_0_zzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyz_0[i] = 4.0 * g_0_zzz_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyz_0[i] * pb_z + g_0_zzzz_0_yyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyzz_0[i] = 4.0 * g_0_zzz_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyzz_0[i] * pb_z + g_0_zzzz_0_yyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyzzz_0[i] = 4.0 * g_0_zzz_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzzz_0[i] * pb_z + g_0_zzzz_0_yyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yzzzz_0[i] = 4.0 * g_0_zzz_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzzz_0[i] * pb_z + g_0_zzzz_0_yzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_zzzzz_0[i] = 4.0 * g_0_zzz_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_zzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_zzzzz_0[i] * pb_z + g_0_zzzz_0_zzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

