#include "ThreeCenterElectronRepulsionPrimRecHSH.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hsh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsh,
                                 size_t idx_eri_0_fsh,
                                 size_t idx_eri_1_fsh,
                                 size_t idx_eri_1_gsg,
                                 size_t idx_eri_1_gsh,
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

    /// Set up components of auxilary buffer : FSH

    auto g_xxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh);

    auto g_xxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 1);

    auto g_xxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 2);

    auto g_xxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 3);

    auto g_xxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 4);

    auto g_xxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 5);

    auto g_xxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 6);

    auto g_xxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 7);

    auto g_xxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 8);

    auto g_xxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 9);

    auto g_xxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 10);

    auto g_xxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 11);

    auto g_xxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 12);

    auto g_xxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 13);

    auto g_xxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 14);

    auto g_xxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 15);

    auto g_xxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 16);

    auto g_xxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 17);

    auto g_xxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 18);

    auto g_xxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 19);

    auto g_xxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 20);

    auto g_xxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 21);

    auto g_xxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 23);

    auto g_xxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 26);

    auto g_xxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 30);

    auto g_xxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 35);

    auto g_xxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 42);

    auto g_xxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 43);

    auto g_xxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 45);

    auto g_xxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 48);

    auto g_xxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 52);

    auto g_xyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 64);

    auto g_xyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 66);

    auto g_xyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 67);

    auto g_xyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 69);

    auto g_xyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 70);

    auto g_xyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 71);

    auto g_xyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 73);

    auto g_xyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 74);

    auto g_xyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 75);

    auto g_xyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 76);

    auto g_xyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 78);

    auto g_xyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 79);

    auto g_xyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 80);

    auto g_xyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 81);

    auto g_xyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 82);

    auto g_xyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 83);

    auto g_xzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 107);

    auto g_xzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 109);

    auto g_xzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 110);

    auto g_xzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 112);

    auto g_xzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 113);

    auto g_xzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 114);

    auto g_xzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 116);

    auto g_xzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 117);

    auto g_xzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 118);

    auto g_xzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 119);

    auto g_xzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 120);

    auto g_xzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 121);

    auto g_xzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 122);

    auto g_xzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 123);

    auto g_xzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 124);

    auto g_xzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 125);

    auto g_yyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 126);

    auto g_yyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 127);

    auto g_yyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 128);

    auto g_yyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 129);

    auto g_yyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 130);

    auto g_yyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 131);

    auto g_yyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 132);

    auto g_yyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 133);

    auto g_yyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 134);

    auto g_yyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 135);

    auto g_yyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 136);

    auto g_yyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 137);

    auto g_yyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 138);

    auto g_yyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 139);

    auto g_yyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 140);

    auto g_yyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 141);

    auto g_yyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 142);

    auto g_yyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 143);

    auto g_yyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 144);

    auto g_yyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 145);

    auto g_yyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 146);

    auto g_yyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 148);

    auto g_yyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 150);

    auto g_yyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 153);

    auto g_yyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 157);

    auto g_yyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 162);

    auto g_yzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 168);

    auto g_yzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 170);

    auto g_yzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 172);

    auto g_yzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 173);

    auto g_yzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 175);

    auto g_yzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 176);

    auto g_yzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 177);

    auto g_yzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 179);

    auto g_yzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 180);

    auto g_yzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 181);

    auto g_yzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 182);

    auto g_yzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 184);

    auto g_yzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 185);

    auto g_yzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 186);

    auto g_yzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 187);

    auto g_yzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 188);

    auto g_zzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 189);

    auto g_zzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 190);

    auto g_zzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 191);

    auto g_zzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 192);

    auto g_zzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 193);

    auto g_zzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 194);

    auto g_zzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 195);

    auto g_zzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 196);

    auto g_zzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 197);

    auto g_zzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 198);

    auto g_zzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 199);

    auto g_zzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 200);

    auto g_zzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 201);

    auto g_zzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 202);

    auto g_zzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 203);

    auto g_zzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 204);

    auto g_zzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 205);

    auto g_zzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 206);

    auto g_zzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 207);

    auto g_zzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 208);

    auto g_zzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 209);

    /// Set up components of auxilary buffer : FSH

    auto g_xxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh);

    auto g_xxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 1);

    auto g_xxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 2);

    auto g_xxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 3);

    auto g_xxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 4);

    auto g_xxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 5);

    auto g_xxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 6);

    auto g_xxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 7);

    auto g_xxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 8);

    auto g_xxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 9);

    auto g_xxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 10);

    auto g_xxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 11);

    auto g_xxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 12);

    auto g_xxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 13);

    auto g_xxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 14);

    auto g_xxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 15);

    auto g_xxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 16);

    auto g_xxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 17);

    auto g_xxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 18);

    auto g_xxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 19);

    auto g_xxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 20);

    auto g_xxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 21);

    auto g_xxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 23);

    auto g_xxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 26);

    auto g_xxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 30);

    auto g_xxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 35);

    auto g_xxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 42);

    auto g_xxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 43);

    auto g_xxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 45);

    auto g_xxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 48);

    auto g_xxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 52);

    auto g_xyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 64);

    auto g_xyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 66);

    auto g_xyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 67);

    auto g_xyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 69);

    auto g_xyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 70);

    auto g_xyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 71);

    auto g_xyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 73);

    auto g_xyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 74);

    auto g_xyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 75);

    auto g_xyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 76);

    auto g_xyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 78);

    auto g_xyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 79);

    auto g_xyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 80);

    auto g_xyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 81);

    auto g_xyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 82);

    auto g_xyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 83);

    auto g_xzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 107);

    auto g_xzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 109);

    auto g_xzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 110);

    auto g_xzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 112);

    auto g_xzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 113);

    auto g_xzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 114);

    auto g_xzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 116);

    auto g_xzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 117);

    auto g_xzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 118);

    auto g_xzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 119);

    auto g_xzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 120);

    auto g_xzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 121);

    auto g_xzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 122);

    auto g_xzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 123);

    auto g_xzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 124);

    auto g_xzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 125);

    auto g_yyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 126);

    auto g_yyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 127);

    auto g_yyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 128);

    auto g_yyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 129);

    auto g_yyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 130);

    auto g_yyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 131);

    auto g_yyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 132);

    auto g_yyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 133);

    auto g_yyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 134);

    auto g_yyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 135);

    auto g_yyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 136);

    auto g_yyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 137);

    auto g_yyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 138);

    auto g_yyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 139);

    auto g_yyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 140);

    auto g_yyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 141);

    auto g_yyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 142);

    auto g_yyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 143);

    auto g_yyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 144);

    auto g_yyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 145);

    auto g_yyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 146);

    auto g_yyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 148);

    auto g_yyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 150);

    auto g_yyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 153);

    auto g_yyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 157);

    auto g_yyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 162);

    auto g_yzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 168);

    auto g_yzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 170);

    auto g_yzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 172);

    auto g_yzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 173);

    auto g_yzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 175);

    auto g_yzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 176);

    auto g_yzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 177);

    auto g_yzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 179);

    auto g_yzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 180);

    auto g_yzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 181);

    auto g_yzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 182);

    auto g_yzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 184);

    auto g_yzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 185);

    auto g_yzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 186);

    auto g_yzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 187);

    auto g_yzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 188);

    auto g_zzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 189);

    auto g_zzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 190);

    auto g_zzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 191);

    auto g_zzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 192);

    auto g_zzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 193);

    auto g_zzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 194);

    auto g_zzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 195);

    auto g_zzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 196);

    auto g_zzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 197);

    auto g_zzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 198);

    auto g_zzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 199);

    auto g_zzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 200);

    auto g_zzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 201);

    auto g_zzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 202);

    auto g_zzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 203);

    auto g_zzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 204);

    auto g_zzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 205);

    auto g_zzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 206);

    auto g_zzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 207);

    auto g_zzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 208);

    auto g_zzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 209);

    /// Set up components of auxilary buffer : GSG

    auto g_xxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg);

    auto g_xxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 1);

    auto g_xxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 2);

    auto g_xxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 3);

    auto g_xxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 4);

    auto g_xxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 5);

    auto g_xxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 6);

    auto g_xxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 7);

    auto g_xxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 8);

    auto g_xxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 9);

    auto g_xxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 10);

    auto g_xxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 11);

    auto g_xxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 12);

    auto g_xxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 13);

    auto g_xxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 14);

    auto g_xxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 32);

    auto g_xxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 34);

    auto g_xxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 35);

    auto g_xxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 37);

    auto g_xxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 38);

    auto g_xxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 39);

    auto g_xxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 41);

    auto g_xxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 42);

    auto g_xxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 43);

    auto g_xxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 44);

    auto g_xxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 45);

    auto g_xxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 46);

    auto g_xxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 47);

    auto g_xxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 48);

    auto g_xxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 49);

    auto g_xxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 50);

    auto g_xxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 51);

    auto g_xxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 52);

    auto g_xxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 53);

    auto g_xxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 54);

    auto g_xxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 55);

    auto g_xxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 56);

    auto g_xxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 57);

    auto g_xxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 58);

    auto g_xxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 59);

    auto g_xxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 75);

    auto g_xxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 76);

    auto g_xxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 77);

    auto g_xxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 78);

    auto g_xxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 79);

    auto g_xxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 80);

    auto g_xxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 81);

    auto g_xxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 82);

    auto g_xxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 83);

    auto g_xxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 84);

    auto g_xxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 85);

    auto g_xxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 86);

    auto g_xxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 87);

    auto g_xxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 88);

    auto g_xxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 89);

    auto g_xyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 91);

    auto g_xyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 93);

    auto g_xyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 94);

    auto g_xyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 96);

    auto g_xyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 97);

    auto g_xyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 98);

    auto g_xyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 100);

    auto g_xyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 101);

    auto g_xyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 102);

    auto g_xyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 103);

    auto g_xzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 137);

    auto g_xzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 139);

    auto g_xzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 140);

    auto g_xzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 142);

    auto g_xzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 143);

    auto g_xzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 144);

    auto g_xzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 146);

    auto g_xzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 147);

    auto g_xzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 148);

    auto g_xzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 149);

    auto g_yyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 150);

    auto g_yyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 151);

    auto g_yyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 152);

    auto g_yyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 153);

    auto g_yyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 154);

    auto g_yyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 155);

    auto g_yyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 156);

    auto g_yyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 157);

    auto g_yyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 158);

    auto g_yyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 159);

    auto g_yyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 160);

    auto g_yyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 161);

    auto g_yyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 162);

    auto g_yyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 163);

    auto g_yyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 164);

    auto g_yyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 167);

    auto g_yyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 169);

    auto g_yyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 170);

    auto g_yyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 172);

    auto g_yyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 173);

    auto g_yyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 174);

    auto g_yyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 176);

    auto g_yyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 177);

    auto g_yyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 178);

    auto g_yyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 179);

    auto g_yyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 180);

    auto g_yyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 181);

    auto g_yyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 182);

    auto g_yyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 183);

    auto g_yyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 184);

    auto g_yyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 185);

    auto g_yyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 186);

    auto g_yyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 187);

    auto g_yyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 188);

    auto g_yyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 189);

    auto g_yyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 190);

    auto g_yyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 191);

    auto g_yyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 192);

    auto g_yyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 193);

    auto g_yyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 194);

    auto g_yzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 196);

    auto g_yzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 197);

    auto g_yzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 198);

    auto g_yzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 199);

    auto g_yzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 200);

    auto g_yzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 201);

    auto g_yzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 202);

    auto g_yzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 203);

    auto g_yzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 204);

    auto g_yzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 205);

    auto g_yzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 206);

    auto g_yzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 207);

    auto g_yzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 208);

    auto g_yzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 209);

    auto g_zzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 210);

    auto g_zzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 211);

    auto g_zzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 212);

    auto g_zzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 213);

    auto g_zzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 214);

    auto g_zzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 215);

    auto g_zzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 216);

    auto g_zzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 217);

    auto g_zzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 218);

    auto g_zzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 219);

    auto g_zzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 220);

    auto g_zzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 221);

    auto g_zzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 222);

    auto g_zzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 223);

    auto g_zzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 224);

    /// Set up components of auxilary buffer : GSH

    auto g_xxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh);

    auto g_xxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 1);

    auto g_xxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 2);

    auto g_xxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 3);

    auto g_xxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 4);

    auto g_xxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 5);

    auto g_xxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 6);

    auto g_xxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 7);

    auto g_xxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 8);

    auto g_xxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 9);

    auto g_xxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 10);

    auto g_xxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 11);

    auto g_xxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 12);

    auto g_xxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 13);

    auto g_xxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 14);

    auto g_xxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 15);

    auto g_xxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 16);

    auto g_xxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 17);

    auto g_xxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 18);

    auto g_xxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 19);

    auto g_xxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 20);

    auto g_xxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 21);

    auto g_xxxy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 22);

    auto g_xxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 23);

    auto g_xxxy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 24);

    auto g_xxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 26);

    auto g_xxxy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 27);

    auto g_xxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 30);

    auto g_xxxy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 31);

    auto g_xxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 35);

    auto g_xxxy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 36);

    auto g_xxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 42);

    auto g_xxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 43);

    auto g_xxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 44);

    auto g_xxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 45);

    auto g_xxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 46);

    auto g_xxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 47);

    auto g_xxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 48);

    auto g_xxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 49);

    auto g_xxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 50);

    auto g_xxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 51);

    auto g_xxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 52);

    auto g_xxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 53);

    auto g_xxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 54);

    auto g_xxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 55);

    auto g_xxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 56);

    auto g_xxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 58);

    auto g_xxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 59);

    auto g_xxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 60);

    auto g_xxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 61);

    auto g_xxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 62);

    auto g_xxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 63);

    auto g_xxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 64);

    auto g_xxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 65);

    auto g_xxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 66);

    auto g_xxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 67);

    auto g_xxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 68);

    auto g_xxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 69);

    auto g_xxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 70);

    auto g_xxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 71);

    auto g_xxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 72);

    auto g_xxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 73);

    auto g_xxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 74);

    auto g_xxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 75);

    auto g_xxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 76);

    auto g_xxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 77);

    auto g_xxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 78);

    auto g_xxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 79);

    auto g_xxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 80);

    auto g_xxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 81);

    auto g_xxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 82);

    auto g_xxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 83);

    auto g_xxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 105);

    auto g_xxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 106);

    auto g_xxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 107);

    auto g_xxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 108);

    auto g_xxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 109);

    auto g_xxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 110);

    auto g_xxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 111);

    auto g_xxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 112);

    auto g_xxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 113);

    auto g_xxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 114);

    auto g_xxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 115);

    auto g_xxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 116);

    auto g_xxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 117);

    auto g_xxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 118);

    auto g_xxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 119);

    auto g_xxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 120);

    auto g_xxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 121);

    auto g_xxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 122);

    auto g_xxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 123);

    auto g_xxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 124);

    auto g_xxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 125);

    auto g_xyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 126);

    auto g_xyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 127);

    auto g_xyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 129);

    auto g_xyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 130);

    auto g_xyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 132);

    auto g_xyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 133);

    auto g_xyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 134);

    auto g_xyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 136);

    auto g_xyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 137);

    auto g_xyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 138);

    auto g_xyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 139);

    auto g_xyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 141);

    auto g_xyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 142);

    auto g_xyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 143);

    auto g_xyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 144);

    auto g_xyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 145);

    auto g_xyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 146);

    auto g_xzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 189);

    auto g_xzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 191);

    auto g_xzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 193);

    auto g_xzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 194);

    auto g_xzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 196);

    auto g_xzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 197);

    auto g_xzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 198);

    auto g_xzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 200);

    auto g_xzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 201);

    auto g_xzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 202);

    auto g_xzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 203);

    auto g_xzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 204);

    auto g_xzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 205);

    auto g_xzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 206);

    auto g_xzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 207);

    auto g_xzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 208);

    auto g_xzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 209);

    auto g_yyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 210);

    auto g_yyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 211);

    auto g_yyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 212);

    auto g_yyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 213);

    auto g_yyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 214);

    auto g_yyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 215);

    auto g_yyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 216);

    auto g_yyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 217);

    auto g_yyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 218);

    auto g_yyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 219);

    auto g_yyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 220);

    auto g_yyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 221);

    auto g_yyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 222);

    auto g_yyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 223);

    auto g_yyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 224);

    auto g_yyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 225);

    auto g_yyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 226);

    auto g_yyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 227);

    auto g_yyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 228);

    auto g_yyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 229);

    auto g_yyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 230);

    auto g_yyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 232);

    auto g_yyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 233);

    auto g_yyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 234);

    auto g_yyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 235);

    auto g_yyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 236);

    auto g_yyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 237);

    auto g_yyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 238);

    auto g_yyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 239);

    auto g_yyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 240);

    auto g_yyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 241);

    auto g_yyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 242);

    auto g_yyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 243);

    auto g_yyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 244);

    auto g_yyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 245);

    auto g_yyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 246);

    auto g_yyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 247);

    auto g_yyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 248);

    auto g_yyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 249);

    auto g_yyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 250);

    auto g_yyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 251);

    auto g_yyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 252);

    auto g_yyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 253);

    auto g_yyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 254);

    auto g_yyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 255);

    auto g_yyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 256);

    auto g_yyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 257);

    auto g_yyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 258);

    auto g_yyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 259);

    auto g_yyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 260);

    auto g_yyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 261);

    auto g_yyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 262);

    auto g_yyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 263);

    auto g_yyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 264);

    auto g_yyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 265);

    auto g_yyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 266);

    auto g_yyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 267);

    auto g_yyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 268);

    auto g_yyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 269);

    auto g_yyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 270);

    auto g_yyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 271);

    auto g_yyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 272);

    auto g_yzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 273);

    auto g_yzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 274);

    auto g_yzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 275);

    auto g_yzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 276);

    auto g_yzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 277);

    auto g_yzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 278);

    auto g_yzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 279);

    auto g_yzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 280);

    auto g_yzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 281);

    auto g_yzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 282);

    auto g_yzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 283);

    auto g_yzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 284);

    auto g_yzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 285);

    auto g_yzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 286);

    auto g_yzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 287);

    auto g_yzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 288);

    auto g_yzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 289);

    auto g_yzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 290);

    auto g_yzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 291);

    auto g_yzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 292);

    auto g_yzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 293);

    auto g_zzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_gsh + 294);

    auto g_zzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_gsh + 295);

    auto g_zzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_gsh + 296);

    auto g_zzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_gsh + 297);

    auto g_zzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_gsh + 298);

    auto g_zzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_gsh + 299);

    auto g_zzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_gsh + 300);

    auto g_zzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_gsh + 301);

    auto g_zzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_gsh + 302);

    auto g_zzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_gsh + 303);

    auto g_zzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_gsh + 304);

    auto g_zzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_gsh + 305);

    auto g_zzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_gsh + 306);

    auto g_zzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_gsh + 307);

    auto g_zzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_gsh + 308);

    auto g_zzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_gsh + 309);

    auto g_zzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_gsh + 310);

    auto g_zzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_gsh + 311);

    auto g_zzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_gsh + 312);

    auto g_zzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_gsh + 313);

    auto g_zzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_gsh + 314);

    /// Set up 0-21 components of targeted buffer : HSH

    auto g_xxxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh);

    auto g_xxxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 1);

    auto g_xxxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 2);

    auto g_xxxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 3);

    auto g_xxxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 4);

    auto g_xxxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 5);

    auto g_xxxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 6);

    auto g_xxxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 7);

    auto g_xxxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 8);

    auto g_xxxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 9);

    auto g_xxxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 10);

    auto g_xxxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 11);

    auto g_xxxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 12);

    auto g_xxxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 13);

    auto g_xxxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 14);

    auto g_xxxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 15);

    auto g_xxxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 16);

    auto g_xxxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 17);

    auto g_xxxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 18);

    auto g_xxxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 19);

    auto g_xxxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 20);

    #pragma omp simd aligned(g_xxx_0_xxxxx_0, g_xxx_0_xxxxx_1, g_xxx_0_xxxxy_0, g_xxx_0_xxxxy_1, g_xxx_0_xxxxz_0, g_xxx_0_xxxxz_1, g_xxx_0_xxxyy_0, g_xxx_0_xxxyy_1, g_xxx_0_xxxyz_0, g_xxx_0_xxxyz_1, g_xxx_0_xxxzz_0, g_xxx_0_xxxzz_1, g_xxx_0_xxyyy_0, g_xxx_0_xxyyy_1, g_xxx_0_xxyyz_0, g_xxx_0_xxyyz_1, g_xxx_0_xxyzz_0, g_xxx_0_xxyzz_1, g_xxx_0_xxzzz_0, g_xxx_0_xxzzz_1, g_xxx_0_xyyyy_0, g_xxx_0_xyyyy_1, g_xxx_0_xyyyz_0, g_xxx_0_xyyyz_1, g_xxx_0_xyyzz_0, g_xxx_0_xyyzz_1, g_xxx_0_xyzzz_0, g_xxx_0_xyzzz_1, g_xxx_0_xzzzz_0, g_xxx_0_xzzzz_1, g_xxx_0_yyyyy_0, g_xxx_0_yyyyy_1, g_xxx_0_yyyyz_0, g_xxx_0_yyyyz_1, g_xxx_0_yyyzz_0, g_xxx_0_yyyzz_1, g_xxx_0_yyzzz_0, g_xxx_0_yyzzz_1, g_xxx_0_yzzzz_0, g_xxx_0_yzzzz_1, g_xxx_0_zzzzz_0, g_xxx_0_zzzzz_1, g_xxxx_0_xxxx_1, g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxy_1, g_xxxx_0_xxxxz_1, g_xxxx_0_xxxy_1, g_xxxx_0_xxxyy_1, g_xxxx_0_xxxyz_1, g_xxxx_0_xxxz_1, g_xxxx_0_xxxzz_1, g_xxxx_0_xxyy_1, g_xxxx_0_xxyyy_1, g_xxxx_0_xxyyz_1, g_xxxx_0_xxyz_1, g_xxxx_0_xxyzz_1, g_xxxx_0_xxzz_1, g_xxxx_0_xxzzz_1, g_xxxx_0_xyyy_1, g_xxxx_0_xyyyy_1, g_xxxx_0_xyyyz_1, g_xxxx_0_xyyz_1, g_xxxx_0_xyyzz_1, g_xxxx_0_xyzz_1, g_xxxx_0_xyzzz_1, g_xxxx_0_xzzz_1, g_xxxx_0_xzzzz_1, g_xxxx_0_yyyy_1, g_xxxx_0_yyyyy_1, g_xxxx_0_yyyyz_1, g_xxxx_0_yyyz_1, g_xxxx_0_yyyzz_1, g_xxxx_0_yyzz_1, g_xxxx_0_yyzzz_1, g_xxxx_0_yzzz_1, g_xxxx_0_yzzzz_1, g_xxxx_0_zzzz_1, g_xxxx_0_zzzzz_1, g_xxxxx_0_xxxxx_0, g_xxxxx_0_xxxxy_0, g_xxxxx_0_xxxxz_0, g_xxxxx_0_xxxyy_0, g_xxxxx_0_xxxyz_0, g_xxxxx_0_xxxzz_0, g_xxxxx_0_xxyyy_0, g_xxxxx_0_xxyyz_0, g_xxxxx_0_xxyzz_0, g_xxxxx_0_xxzzz_0, g_xxxxx_0_xyyyy_0, g_xxxxx_0_xyyyz_0, g_xxxxx_0_xyyzz_0, g_xxxxx_0_xyzzz_0, g_xxxxx_0_xzzzz_0, g_xxxxx_0_yyyyy_0, g_xxxxx_0_yyyyz_0, g_xxxxx_0_yyyzz_0, g_xxxxx_0_yyzzz_0, g_xxxxx_0_yzzzz_0, g_xxxxx_0_zzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxx_0_xxxxx_0[i] = 4.0 * g_xxx_0_xxxxx_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxx_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxx_1[i] * wa_x[i];

        g_xxxxx_0_xxxxy_0[i] = 4.0 * g_xxx_0_xxxxy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxz_0[i] = 4.0 * g_xxx_0_xxxxz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyy_0[i] = 4.0 * g_xxx_0_xxxyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxyz_0[i] = 4.0 * g_xxx_0_xxxyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxzz_0[i] = 4.0 * g_xxx_0_xxxzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyy_0[i] = 4.0 * g_xxx_0_xxyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxyyz_0[i] = 4.0 * g_xxx_0_xxyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxyzz_0[i] = 4.0 * g_xxx_0_xxyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxzzz_0[i] = 4.0 * g_xxx_0_xxzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxx_0_xxzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyy_0[i] = 4.0 * g_xxx_0_xyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyy_1[i] * fz_be_0 + g_xxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xyyyz_0[i] = 4.0 * g_xxx_0_xyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyz_1[i] * fz_be_0 + g_xxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xyyzz_0[i] = 4.0 * g_xxx_0_xyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyzz_1[i] * fz_be_0 + g_xxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xyzzz_0[i] = 4.0 * g_xxx_0_xyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyzzz_1[i] * fz_be_0 + g_xxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xzzzz_0[i] = 4.0 * g_xxx_0_xzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xzzzz_1[i] * fz_be_0 + g_xxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxx_0_xzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyy_0[i] = 4.0 * g_xxx_0_yyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyy_1[i] * fz_be_0 + g_xxxx_0_yyyyy_1[i] * wa_x[i];

        g_xxxxx_0_yyyyz_0[i] = 4.0 * g_xxx_0_yyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyz_1[i] * fz_be_0 + g_xxxx_0_yyyyz_1[i] * wa_x[i];

        g_xxxxx_0_yyyzz_0[i] = 4.0 * g_xxx_0_yyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyzz_1[i] * fz_be_0 + g_xxxx_0_yyyzz_1[i] * wa_x[i];

        g_xxxxx_0_yyzzz_0[i] = 4.0 * g_xxx_0_yyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyzzz_1[i] * fz_be_0 + g_xxxx_0_yyzzz_1[i] * wa_x[i];

        g_xxxxx_0_yzzzz_0[i] = 4.0 * g_xxx_0_yzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yzzzz_1[i] * fz_be_0 + g_xxxx_0_yzzzz_1[i] * wa_x[i];

        g_xxxxx_0_zzzzz_0[i] = 4.0 * g_xxx_0_zzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_zzzzz_1[i] * fz_be_0 + g_xxxx_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 21-42 components of targeted buffer : HSH

    auto g_xxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 21);

    auto g_xxxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 22);

    auto g_xxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 23);

    auto g_xxxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 24);

    auto g_xxxxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 25);

    auto g_xxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 26);

    auto g_xxxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 27);

    auto g_xxxxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 28);

    auto g_xxxxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 29);

    auto g_xxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 30);

    auto g_xxxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 31);

    auto g_xxxxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 32);

    auto g_xxxxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 33);

    auto g_xxxxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 34);

    auto g_xxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 35);

    auto g_xxxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 36);

    auto g_xxxxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 37);

    auto g_xxxxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 38);

    auto g_xxxxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 39);

    auto g_xxxxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 40);

    auto g_xxxxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 41);

    #pragma omp simd aligned(g_xxxx_0_xxxx_1, g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxy_1, g_xxxx_0_xxxxz_1, g_xxxx_0_xxxy_1, g_xxxx_0_xxxyy_1, g_xxxx_0_xxxyz_1, g_xxxx_0_xxxz_1, g_xxxx_0_xxxzz_1, g_xxxx_0_xxyy_1, g_xxxx_0_xxyyy_1, g_xxxx_0_xxyyz_1, g_xxxx_0_xxyz_1, g_xxxx_0_xxyzz_1, g_xxxx_0_xxzz_1, g_xxxx_0_xxzzz_1, g_xxxx_0_xyyy_1, g_xxxx_0_xyyyy_1, g_xxxx_0_xyyyz_1, g_xxxx_0_xyyz_1, g_xxxx_0_xyyzz_1, g_xxxx_0_xyzz_1, g_xxxx_0_xyzzz_1, g_xxxx_0_xzzz_1, g_xxxx_0_xzzzz_1, g_xxxx_0_yyyy_1, g_xxxx_0_yyyyy_1, g_xxxx_0_yyyyz_1, g_xxxx_0_yyyz_1, g_xxxx_0_yyyzz_1, g_xxxx_0_yyzz_1, g_xxxx_0_yyzzz_1, g_xxxx_0_yzzz_1, g_xxxx_0_yzzzz_1, g_xxxx_0_zzzz_1, g_xxxx_0_zzzzz_1, g_xxxxy_0_xxxxx_0, g_xxxxy_0_xxxxy_0, g_xxxxy_0_xxxxz_0, g_xxxxy_0_xxxyy_0, g_xxxxy_0_xxxyz_0, g_xxxxy_0_xxxzz_0, g_xxxxy_0_xxyyy_0, g_xxxxy_0_xxyyz_0, g_xxxxy_0_xxyzz_0, g_xxxxy_0_xxzzz_0, g_xxxxy_0_xyyyy_0, g_xxxxy_0_xyyyz_0, g_xxxxy_0_xyyzz_0, g_xxxxy_0_xyzzz_0, g_xxxxy_0_xzzzz_0, g_xxxxy_0_yyyyy_0, g_xxxxy_0_yyyyz_0, g_xxxxy_0_yyyzz_0, g_xxxxy_0_yyzzz_0, g_xxxxy_0_yzzzz_0, g_xxxxy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxy_0_xxxxx_0[i] = g_xxxx_0_xxxxx_1[i] * wa_y[i];

        g_xxxxy_0_xxxxy_0[i] = g_xxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxz_0[i] = g_xxxx_0_xxxxz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyy_0[i] = 2.0 * g_xxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxyz_0[i] = g_xxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxzz_0[i] = g_xxxx_0_xxxzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyy_0[i] = 3.0 * g_xxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxyyz_0[i] = 2.0 * g_xxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxyzz_0[i] = g_xxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxzzz_0[i] = g_xxxx_0_xxzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyy_0[i] = 4.0 * g_xxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xyyyz_0[i] = 3.0 * g_xxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xyyzz_0[i] = 2.0 * g_xxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xyzzz_0[i] = g_xxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xzzzz_0[i] = g_xxxx_0_xzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyy_0[i] = 5.0 * g_xxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxx_0_yyyyy_1[i] * wa_y[i];

        g_xxxxy_0_yyyyz_0[i] = 4.0 * g_xxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxx_0_yyyyz_1[i] * wa_y[i];

        g_xxxxy_0_yyyzz_0[i] = 3.0 * g_xxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxx_0_yyyzz_1[i] * wa_y[i];

        g_xxxxy_0_yyzzz_0[i] = 2.0 * g_xxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxx_0_yyzzz_1[i] * wa_y[i];

        g_xxxxy_0_yzzzz_0[i] = g_xxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxx_0_yzzzz_1[i] * wa_y[i];

        g_xxxxy_0_zzzzz_0[i] = g_xxxx_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 42-63 components of targeted buffer : HSH

    auto g_xxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 42);

    auto g_xxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 43);

    auto g_xxxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 44);

    auto g_xxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 45);

    auto g_xxxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 46);

    auto g_xxxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 47);

    auto g_xxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 48);

    auto g_xxxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 49);

    auto g_xxxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 50);

    auto g_xxxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 51);

    auto g_xxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 52);

    auto g_xxxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 53);

    auto g_xxxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 54);

    auto g_xxxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 55);

    auto g_xxxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 56);

    auto g_xxxxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 57);

    auto g_xxxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 58);

    auto g_xxxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 59);

    auto g_xxxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 60);

    auto g_xxxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 61);

    auto g_xxxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 62);

    #pragma omp simd aligned(g_xxxx_0_xxxx_1, g_xxxx_0_xxxxx_1, g_xxxx_0_xxxxy_1, g_xxxx_0_xxxxz_1, g_xxxx_0_xxxy_1, g_xxxx_0_xxxyy_1, g_xxxx_0_xxxyz_1, g_xxxx_0_xxxz_1, g_xxxx_0_xxxzz_1, g_xxxx_0_xxyy_1, g_xxxx_0_xxyyy_1, g_xxxx_0_xxyyz_1, g_xxxx_0_xxyz_1, g_xxxx_0_xxyzz_1, g_xxxx_0_xxzz_1, g_xxxx_0_xxzzz_1, g_xxxx_0_xyyy_1, g_xxxx_0_xyyyy_1, g_xxxx_0_xyyyz_1, g_xxxx_0_xyyz_1, g_xxxx_0_xyyzz_1, g_xxxx_0_xyzz_1, g_xxxx_0_xyzzz_1, g_xxxx_0_xzzz_1, g_xxxx_0_xzzzz_1, g_xxxx_0_yyyy_1, g_xxxx_0_yyyyy_1, g_xxxx_0_yyyyz_1, g_xxxx_0_yyyz_1, g_xxxx_0_yyyzz_1, g_xxxx_0_yyzz_1, g_xxxx_0_yyzzz_1, g_xxxx_0_yzzz_1, g_xxxx_0_yzzzz_1, g_xxxx_0_zzzz_1, g_xxxx_0_zzzzz_1, g_xxxxz_0_xxxxx_0, g_xxxxz_0_xxxxy_0, g_xxxxz_0_xxxxz_0, g_xxxxz_0_xxxyy_0, g_xxxxz_0_xxxyz_0, g_xxxxz_0_xxxzz_0, g_xxxxz_0_xxyyy_0, g_xxxxz_0_xxyyz_0, g_xxxxz_0_xxyzz_0, g_xxxxz_0_xxzzz_0, g_xxxxz_0_xyyyy_0, g_xxxxz_0_xyyyz_0, g_xxxxz_0_xyyzz_0, g_xxxxz_0_xyzzz_0, g_xxxxz_0_xzzzz_0, g_xxxxz_0_yyyyy_0, g_xxxxz_0_yyyyz_0, g_xxxxz_0_yyyzz_0, g_xxxxz_0_yyzzz_0, g_xxxxz_0_yzzzz_0, g_xxxxz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxz_0_xxxxx_0[i] = g_xxxx_0_xxxxx_1[i] * wa_z[i];

        g_xxxxz_0_xxxxy_0[i] = g_xxxx_0_xxxxy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxz_0[i] = g_xxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyy_0[i] = g_xxxx_0_xxxyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxyz_0[i] = g_xxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxzz_0[i] = 2.0 * g_xxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyy_0[i] = g_xxxx_0_xxyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxyyz_0[i] = g_xxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxyzz_0[i] = 2.0 * g_xxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxx_0_xxyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxzzz_0[i] = 3.0 * g_xxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxx_0_xxzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyy_0[i] = g_xxxx_0_xyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xyyyz_0[i] = g_xxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xyyzz_0[i] = 2.0 * g_xxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xyzzz_0[i] = 3.0 * g_xxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xzzzz_0[i] = 4.0 * g_xxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxx_0_xzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyy_0[i] = g_xxxx_0_yyyyy_1[i] * wa_z[i];

        g_xxxxz_0_yyyyz_0[i] = g_xxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxx_0_yyyyz_1[i] * wa_z[i];

        g_xxxxz_0_yyyzz_0[i] = 2.0 * g_xxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxx_0_yyyzz_1[i] * wa_z[i];

        g_xxxxz_0_yyzzz_0[i] = 3.0 * g_xxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxx_0_yyzzz_1[i] * wa_z[i];

        g_xxxxz_0_yzzzz_0[i] = 4.0 * g_xxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxx_0_yzzzz_1[i] * wa_z[i];

        g_xxxxz_0_zzzzz_0[i] = 5.0 * g_xxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxx_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 63-84 components of targeted buffer : HSH

    auto g_xxxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 63);

    auto g_xxxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 64);

    auto g_xxxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 65);

    auto g_xxxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 66);

    auto g_xxxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 67);

    auto g_xxxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 68);

    auto g_xxxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 69);

    auto g_xxxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 70);

    auto g_xxxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 71);

    auto g_xxxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 72);

    auto g_xxxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 73);

    auto g_xxxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 74);

    auto g_xxxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 75);

    auto g_xxxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 76);

    auto g_xxxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 77);

    auto g_xxxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 78);

    auto g_xxxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 79);

    auto g_xxxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 80);

    auto g_xxxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 81);

    auto g_xxxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 82);

    auto g_xxxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 83);

    #pragma omp simd aligned(g_xxx_0_xxxxx_0, g_xxx_0_xxxxx_1, g_xxx_0_xxxxz_0, g_xxx_0_xxxxz_1, g_xxx_0_xxxzz_0, g_xxx_0_xxxzz_1, g_xxx_0_xxzzz_0, g_xxx_0_xxzzz_1, g_xxx_0_xzzzz_0, g_xxx_0_xzzzz_1, g_xxxy_0_xxxxx_1, g_xxxy_0_xxxxz_1, g_xxxy_0_xxxzz_1, g_xxxy_0_xxzzz_1, g_xxxy_0_xzzzz_1, g_xxxyy_0_xxxxx_0, g_xxxyy_0_xxxxy_0, g_xxxyy_0_xxxxz_0, g_xxxyy_0_xxxyy_0, g_xxxyy_0_xxxyz_0, g_xxxyy_0_xxxzz_0, g_xxxyy_0_xxyyy_0, g_xxxyy_0_xxyyz_0, g_xxxyy_0_xxyzz_0, g_xxxyy_0_xxzzz_0, g_xxxyy_0_xyyyy_0, g_xxxyy_0_xyyyz_0, g_xxxyy_0_xyyzz_0, g_xxxyy_0_xyzzz_0, g_xxxyy_0_xzzzz_0, g_xxxyy_0_yyyyy_0, g_xxxyy_0_yyyyz_0, g_xxxyy_0_yyyzz_0, g_xxxyy_0_yyzzz_0, g_xxxyy_0_yzzzz_0, g_xxxyy_0_zzzzz_0, g_xxyy_0_xxxxy_1, g_xxyy_0_xxxy_1, g_xxyy_0_xxxyy_1, g_xxyy_0_xxxyz_1, g_xxyy_0_xxyy_1, g_xxyy_0_xxyyy_1, g_xxyy_0_xxyyz_1, g_xxyy_0_xxyz_1, g_xxyy_0_xxyzz_1, g_xxyy_0_xyyy_1, g_xxyy_0_xyyyy_1, g_xxyy_0_xyyyz_1, g_xxyy_0_xyyz_1, g_xxyy_0_xyyzz_1, g_xxyy_0_xyzz_1, g_xxyy_0_xyzzz_1, g_xxyy_0_yyyy_1, g_xxyy_0_yyyyy_1, g_xxyy_0_yyyyz_1, g_xxyy_0_yyyz_1, g_xxyy_0_yyyzz_1, g_xxyy_0_yyzz_1, g_xxyy_0_yyzzz_1, g_xxyy_0_yzzz_1, g_xxyy_0_yzzzz_1, g_xxyy_0_zzzzz_1, g_xyy_0_xxxxy_0, g_xyy_0_xxxxy_1, g_xyy_0_xxxyy_0, g_xyy_0_xxxyy_1, g_xyy_0_xxxyz_0, g_xyy_0_xxxyz_1, g_xyy_0_xxyyy_0, g_xyy_0_xxyyy_1, g_xyy_0_xxyyz_0, g_xyy_0_xxyyz_1, g_xyy_0_xxyzz_0, g_xyy_0_xxyzz_1, g_xyy_0_xyyyy_0, g_xyy_0_xyyyy_1, g_xyy_0_xyyyz_0, g_xyy_0_xyyyz_1, g_xyy_0_xyyzz_0, g_xyy_0_xyyzz_1, g_xyy_0_xyzzz_0, g_xyy_0_xyzzz_1, g_xyy_0_yyyyy_0, g_xyy_0_yyyyy_1, g_xyy_0_yyyyz_0, g_xyy_0_yyyyz_1, g_xyy_0_yyyzz_0, g_xyy_0_yyyzz_1, g_xyy_0_yyzzz_0, g_xyy_0_yyzzz_1, g_xyy_0_yzzzz_0, g_xyy_0_yzzzz_1, g_xyy_0_zzzzz_0, g_xyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyy_0_xxxxx_0[i] = g_xxx_0_xxxxx_0[i] * fbe_0 - g_xxx_0_xxxxx_1[i] * fz_be_0 + g_xxxy_0_xxxxx_1[i] * wa_y[i];

        g_xxxyy_0_xxxxy_0[i] = 2.0 * g_xyy_0_xxxxy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxy_1[i] * fi_acd_0 + g_xxyy_0_xxxxy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxz_0[i] = g_xxx_0_xxxxz_0[i] * fbe_0 - g_xxx_0_xxxxz_1[i] * fz_be_0 + g_xxxy_0_xxxxz_1[i] * wa_y[i];

        g_xxxyy_0_xxxyy_0[i] = 2.0 * g_xyy_0_xxxyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyy_1[i] * fi_acd_0 + g_xxyy_0_xxxyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxyz_0[i] = 2.0 * g_xyy_0_xxxyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyz_1[i] * fi_acd_0 + g_xxyy_0_xxxyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxzz_0[i] = g_xxx_0_xxxzz_0[i] * fbe_0 - g_xxx_0_xxxzz_1[i] * fz_be_0 + g_xxxy_0_xxxzz_1[i] * wa_y[i];

        g_xxxyy_0_xxyyy_0[i] = 2.0 * g_xyy_0_xxyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyy_1[i] * fi_acd_0 + g_xxyy_0_xxyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxyyz_0[i] = 2.0 * g_xyy_0_xxyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyz_1[i] * fi_acd_0 + g_xxyy_0_xxyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxyzz_0[i] = 2.0 * g_xyy_0_xxyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyzz_1[i] * fi_acd_0 + g_xxyy_0_xxyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxzzz_0[i] = g_xxx_0_xxzzz_0[i] * fbe_0 - g_xxx_0_xxzzz_1[i] * fz_be_0 + g_xxxy_0_xxzzz_1[i] * wa_y[i];

        g_xxxyy_0_xyyyy_0[i] = 2.0 * g_xyy_0_xyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyy_1[i] * fz_be_0 + g_xxyy_0_yyyy_1[i] * fi_acd_0 + g_xxyy_0_xyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xyyyz_0[i] = 2.0 * g_xyy_0_xyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyz_1[i] * fz_be_0 + g_xxyy_0_yyyz_1[i] * fi_acd_0 + g_xxyy_0_xyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xyyzz_0[i] = 2.0 * g_xyy_0_xyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyzz_1[i] * fz_be_0 + g_xxyy_0_yyzz_1[i] * fi_acd_0 + g_xxyy_0_xyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xyzzz_0[i] = 2.0 * g_xyy_0_xyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyzzz_1[i] * fz_be_0 + g_xxyy_0_yzzz_1[i] * fi_acd_0 + g_xxyy_0_xyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xzzzz_0[i] = g_xxx_0_xzzzz_0[i] * fbe_0 - g_xxx_0_xzzzz_1[i] * fz_be_0 + g_xxxy_0_xzzzz_1[i] * wa_y[i];

        g_xxxyy_0_yyyyy_0[i] = 2.0 * g_xyy_0_yyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyy_1[i] * fz_be_0 + g_xxyy_0_yyyyy_1[i] * wa_x[i];

        g_xxxyy_0_yyyyz_0[i] = 2.0 * g_xyy_0_yyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyz_1[i] * fz_be_0 + g_xxyy_0_yyyyz_1[i] * wa_x[i];

        g_xxxyy_0_yyyzz_0[i] = 2.0 * g_xyy_0_yyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyzz_1[i] * fz_be_0 + g_xxyy_0_yyyzz_1[i] * wa_x[i];

        g_xxxyy_0_yyzzz_0[i] = 2.0 * g_xyy_0_yyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyzzz_1[i] * fz_be_0 + g_xxyy_0_yyzzz_1[i] * wa_x[i];

        g_xxxyy_0_yzzzz_0[i] = 2.0 * g_xyy_0_yzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yzzzz_1[i] * fz_be_0 + g_xxyy_0_yzzzz_1[i] * wa_x[i];

        g_xxxyy_0_zzzzz_0[i] = 2.0 * g_xyy_0_zzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_zzzzz_1[i] * fz_be_0 + g_xxyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 84-105 components of targeted buffer : HSH

    auto g_xxxyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 84);

    auto g_xxxyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 85);

    auto g_xxxyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 86);

    auto g_xxxyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 87);

    auto g_xxxyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 88);

    auto g_xxxyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 89);

    auto g_xxxyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 90);

    auto g_xxxyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 91);

    auto g_xxxyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 92);

    auto g_xxxyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 93);

    auto g_xxxyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 94);

    auto g_xxxyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 95);

    auto g_xxxyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 96);

    auto g_xxxyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 97);

    auto g_xxxyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 98);

    auto g_xxxyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 99);

    auto g_xxxyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 100);

    auto g_xxxyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 101);

    auto g_xxxyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 102);

    auto g_xxxyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 103);

    auto g_xxxyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 104);

    #pragma omp simd aligned(g_xxxy_0_xxxxy_1, g_xxxy_0_xxxyy_1, g_xxxy_0_xxyyy_1, g_xxxy_0_xyyyy_1, g_xxxy_0_yyyyy_1, g_xxxyz_0_xxxxx_0, g_xxxyz_0_xxxxy_0, g_xxxyz_0_xxxxz_0, g_xxxyz_0_xxxyy_0, g_xxxyz_0_xxxyz_0, g_xxxyz_0_xxxzz_0, g_xxxyz_0_xxyyy_0, g_xxxyz_0_xxyyz_0, g_xxxyz_0_xxyzz_0, g_xxxyz_0_xxzzz_0, g_xxxyz_0_xyyyy_0, g_xxxyz_0_xyyyz_0, g_xxxyz_0_xyyzz_0, g_xxxyz_0_xyzzz_0, g_xxxyz_0_xzzzz_0, g_xxxyz_0_yyyyy_0, g_xxxyz_0_yyyyz_0, g_xxxyz_0_yyyzz_0, g_xxxyz_0_yyzzz_0, g_xxxyz_0_yzzzz_0, g_xxxyz_0_zzzzz_0, g_xxxz_0_xxxxx_1, g_xxxz_0_xxxxz_1, g_xxxz_0_xxxyz_1, g_xxxz_0_xxxz_1, g_xxxz_0_xxxzz_1, g_xxxz_0_xxyyz_1, g_xxxz_0_xxyz_1, g_xxxz_0_xxyzz_1, g_xxxz_0_xxzz_1, g_xxxz_0_xxzzz_1, g_xxxz_0_xyyyz_1, g_xxxz_0_xyyz_1, g_xxxz_0_xyyzz_1, g_xxxz_0_xyzz_1, g_xxxz_0_xyzzz_1, g_xxxz_0_xzzz_1, g_xxxz_0_xzzzz_1, g_xxxz_0_yyyyz_1, g_xxxz_0_yyyz_1, g_xxxz_0_yyyzz_1, g_xxxz_0_yyzz_1, g_xxxz_0_yyzzz_1, g_xxxz_0_yzzz_1, g_xxxz_0_yzzzz_1, g_xxxz_0_zzzz_1, g_xxxz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyz_0_xxxxx_0[i] = g_xxxz_0_xxxxx_1[i] * wa_y[i];

        g_xxxyz_0_xxxxy_0[i] = g_xxxy_0_xxxxy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxz_0[i] = g_xxxz_0_xxxxz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyy_0[i] = g_xxxy_0_xxxyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxyz_0[i] = g_xxxz_0_xxxz_1[i] * fi_acd_0 + g_xxxz_0_xxxyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxzz_0[i] = g_xxxz_0_xxxzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyy_0[i] = g_xxxy_0_xxyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxyyz_0[i] = 2.0 * g_xxxz_0_xxyz_1[i] * fi_acd_0 + g_xxxz_0_xxyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxyzz_0[i] = g_xxxz_0_xxzz_1[i] * fi_acd_0 + g_xxxz_0_xxyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxzzz_0[i] = g_xxxz_0_xxzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyy_0[i] = g_xxxy_0_xyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xyyyz_0[i] = 3.0 * g_xxxz_0_xyyz_1[i] * fi_acd_0 + g_xxxz_0_xyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xyyzz_0[i] = 2.0 * g_xxxz_0_xyzz_1[i] * fi_acd_0 + g_xxxz_0_xyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xyzzz_0[i] = g_xxxz_0_xzzz_1[i] * fi_acd_0 + g_xxxz_0_xyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xzzzz_0[i] = g_xxxz_0_xzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyy_0[i] = g_xxxy_0_yyyyy_1[i] * wa_z[i];

        g_xxxyz_0_yyyyz_0[i] = 4.0 * g_xxxz_0_yyyz_1[i] * fi_acd_0 + g_xxxz_0_yyyyz_1[i] * wa_y[i];

        g_xxxyz_0_yyyzz_0[i] = 3.0 * g_xxxz_0_yyzz_1[i] * fi_acd_0 + g_xxxz_0_yyyzz_1[i] * wa_y[i];

        g_xxxyz_0_yyzzz_0[i] = 2.0 * g_xxxz_0_yzzz_1[i] * fi_acd_0 + g_xxxz_0_yyzzz_1[i] * wa_y[i];

        g_xxxyz_0_yzzzz_0[i] = g_xxxz_0_zzzz_1[i] * fi_acd_0 + g_xxxz_0_yzzzz_1[i] * wa_y[i];

        g_xxxyz_0_zzzzz_0[i] = g_xxxz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 105-126 components of targeted buffer : HSH

    auto g_xxxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 105);

    auto g_xxxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 106);

    auto g_xxxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 107);

    auto g_xxxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 108);

    auto g_xxxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 109);

    auto g_xxxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 110);

    auto g_xxxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 111);

    auto g_xxxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 112);

    auto g_xxxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 113);

    auto g_xxxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 114);

    auto g_xxxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 115);

    auto g_xxxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 116);

    auto g_xxxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 117);

    auto g_xxxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 118);

    auto g_xxxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 119);

    auto g_xxxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 120);

    auto g_xxxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 121);

    auto g_xxxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 122);

    auto g_xxxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 123);

    auto g_xxxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 124);

    auto g_xxxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 125);

    #pragma omp simd aligned(g_xxx_0_xxxxx_0, g_xxx_0_xxxxx_1, g_xxx_0_xxxxy_0, g_xxx_0_xxxxy_1, g_xxx_0_xxxyy_0, g_xxx_0_xxxyy_1, g_xxx_0_xxyyy_0, g_xxx_0_xxyyy_1, g_xxx_0_xyyyy_0, g_xxx_0_xyyyy_1, g_xxxz_0_xxxxx_1, g_xxxz_0_xxxxy_1, g_xxxz_0_xxxyy_1, g_xxxz_0_xxyyy_1, g_xxxz_0_xyyyy_1, g_xxxzz_0_xxxxx_0, g_xxxzz_0_xxxxy_0, g_xxxzz_0_xxxxz_0, g_xxxzz_0_xxxyy_0, g_xxxzz_0_xxxyz_0, g_xxxzz_0_xxxzz_0, g_xxxzz_0_xxyyy_0, g_xxxzz_0_xxyyz_0, g_xxxzz_0_xxyzz_0, g_xxxzz_0_xxzzz_0, g_xxxzz_0_xyyyy_0, g_xxxzz_0_xyyyz_0, g_xxxzz_0_xyyzz_0, g_xxxzz_0_xyzzz_0, g_xxxzz_0_xzzzz_0, g_xxxzz_0_yyyyy_0, g_xxxzz_0_yyyyz_0, g_xxxzz_0_yyyzz_0, g_xxxzz_0_yyzzz_0, g_xxxzz_0_yzzzz_0, g_xxxzz_0_zzzzz_0, g_xxzz_0_xxxxz_1, g_xxzz_0_xxxyz_1, g_xxzz_0_xxxz_1, g_xxzz_0_xxxzz_1, g_xxzz_0_xxyyz_1, g_xxzz_0_xxyz_1, g_xxzz_0_xxyzz_1, g_xxzz_0_xxzz_1, g_xxzz_0_xxzzz_1, g_xxzz_0_xyyyz_1, g_xxzz_0_xyyz_1, g_xxzz_0_xyyzz_1, g_xxzz_0_xyzz_1, g_xxzz_0_xyzzz_1, g_xxzz_0_xzzz_1, g_xxzz_0_xzzzz_1, g_xxzz_0_yyyyy_1, g_xxzz_0_yyyyz_1, g_xxzz_0_yyyz_1, g_xxzz_0_yyyzz_1, g_xxzz_0_yyzz_1, g_xxzz_0_yyzzz_1, g_xxzz_0_yzzz_1, g_xxzz_0_yzzzz_1, g_xxzz_0_zzzz_1, g_xxzz_0_zzzzz_1, g_xzz_0_xxxxz_0, g_xzz_0_xxxxz_1, g_xzz_0_xxxyz_0, g_xzz_0_xxxyz_1, g_xzz_0_xxxzz_0, g_xzz_0_xxxzz_1, g_xzz_0_xxyyz_0, g_xzz_0_xxyyz_1, g_xzz_0_xxyzz_0, g_xzz_0_xxyzz_1, g_xzz_0_xxzzz_0, g_xzz_0_xxzzz_1, g_xzz_0_xyyyz_0, g_xzz_0_xyyyz_1, g_xzz_0_xyyzz_0, g_xzz_0_xyyzz_1, g_xzz_0_xyzzz_0, g_xzz_0_xyzzz_1, g_xzz_0_xzzzz_0, g_xzz_0_xzzzz_1, g_xzz_0_yyyyy_0, g_xzz_0_yyyyy_1, g_xzz_0_yyyyz_0, g_xzz_0_yyyyz_1, g_xzz_0_yyyzz_0, g_xzz_0_yyyzz_1, g_xzz_0_yyzzz_0, g_xzz_0_yyzzz_1, g_xzz_0_yzzzz_0, g_xzz_0_yzzzz_1, g_xzz_0_zzzzz_0, g_xzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzz_0_xxxxx_0[i] = g_xxx_0_xxxxx_0[i] * fbe_0 - g_xxx_0_xxxxx_1[i] * fz_be_0 + g_xxxz_0_xxxxx_1[i] * wa_z[i];

        g_xxxzz_0_xxxxy_0[i] = g_xxx_0_xxxxy_0[i] * fbe_0 - g_xxx_0_xxxxy_1[i] * fz_be_0 + g_xxxz_0_xxxxy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxz_0[i] = 2.0 * g_xzz_0_xxxxz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxz_1[i] * fi_acd_0 + g_xxzz_0_xxxxz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyy_0[i] = g_xxx_0_xxxyy_0[i] * fbe_0 - g_xxx_0_xxxyy_1[i] * fz_be_0 + g_xxxz_0_xxxyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxyz_0[i] = 2.0 * g_xzz_0_xxxyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyz_1[i] * fi_acd_0 + g_xxzz_0_xxxyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxzz_0[i] = 2.0 * g_xzz_0_xxxzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxzz_1[i] * fi_acd_0 + g_xxzz_0_xxxzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyy_0[i] = g_xxx_0_xxyyy_0[i] * fbe_0 - g_xxx_0_xxyyy_1[i] * fz_be_0 + g_xxxz_0_xxyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxyyz_0[i] = 2.0 * g_xzz_0_xxyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyz_1[i] * fi_acd_0 + g_xxzz_0_xxyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxyzz_0[i] = 2.0 * g_xzz_0_xxyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyzz_1[i] * fi_acd_0 + g_xxzz_0_xxyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxzzz_0[i] = 2.0 * g_xzz_0_xxzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xzzz_1[i] * fi_acd_0 + g_xxzz_0_xxzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyy_0[i] = g_xxx_0_xyyyy_0[i] * fbe_0 - g_xxx_0_xyyyy_1[i] * fz_be_0 + g_xxxz_0_xyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xyyyz_0[i] = 2.0 * g_xzz_0_xyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyz_1[i] * fz_be_0 + g_xxzz_0_yyyz_1[i] * fi_acd_0 + g_xxzz_0_xyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xyyzz_0[i] = 2.0 * g_xzz_0_xyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyzz_1[i] * fz_be_0 + g_xxzz_0_yyzz_1[i] * fi_acd_0 + g_xxzz_0_xyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xyzzz_0[i] = 2.0 * g_xzz_0_xyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyzzz_1[i] * fz_be_0 + g_xxzz_0_yzzz_1[i] * fi_acd_0 + g_xxzz_0_xyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xzzzz_0[i] = 2.0 * g_xzz_0_xzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xzzzz_1[i] * fz_be_0 + g_xxzz_0_zzzz_1[i] * fi_acd_0 + g_xxzz_0_xzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyy_0[i] = 2.0 * g_xzz_0_yyyyy_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyy_1[i] * fz_be_0 + g_xxzz_0_yyyyy_1[i] * wa_x[i];

        g_xxxzz_0_yyyyz_0[i] = 2.0 * g_xzz_0_yyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyz_1[i] * fz_be_0 + g_xxzz_0_yyyyz_1[i] * wa_x[i];

        g_xxxzz_0_yyyzz_0[i] = 2.0 * g_xzz_0_yyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyzz_1[i] * fz_be_0 + g_xxzz_0_yyyzz_1[i] * wa_x[i];

        g_xxxzz_0_yyzzz_0[i] = 2.0 * g_xzz_0_yyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyzzz_1[i] * fz_be_0 + g_xxzz_0_yyzzz_1[i] * wa_x[i];

        g_xxxzz_0_yzzzz_0[i] = 2.0 * g_xzz_0_yzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yzzzz_1[i] * fz_be_0 + g_xxzz_0_yzzzz_1[i] * wa_x[i];

        g_xxxzz_0_zzzzz_0[i] = 2.0 * g_xzz_0_zzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_zzzzz_1[i] * fz_be_0 + g_xxzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 126-147 components of targeted buffer : HSH

    auto g_xxyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 126);

    auto g_xxyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 127);

    auto g_xxyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 128);

    auto g_xxyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 129);

    auto g_xxyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 130);

    auto g_xxyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 131);

    auto g_xxyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 132);

    auto g_xxyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 133);

    auto g_xxyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 134);

    auto g_xxyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 135);

    auto g_xxyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 136);

    auto g_xxyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 137);

    auto g_xxyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 138);

    auto g_xxyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 139);

    auto g_xxyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 140);

    auto g_xxyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 141);

    auto g_xxyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 142);

    auto g_xxyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 143);

    auto g_xxyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 144);

    auto g_xxyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 145);

    auto g_xxyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 146);

    #pragma omp simd aligned(g_xxy_0_xxxxx_0, g_xxy_0_xxxxx_1, g_xxy_0_xxxxz_0, g_xxy_0_xxxxz_1, g_xxy_0_xxxzz_0, g_xxy_0_xxxzz_1, g_xxy_0_xxzzz_0, g_xxy_0_xxzzz_1, g_xxy_0_xzzzz_0, g_xxy_0_xzzzz_1, g_xxyy_0_xxxxx_1, g_xxyy_0_xxxxz_1, g_xxyy_0_xxxzz_1, g_xxyy_0_xxzzz_1, g_xxyy_0_xzzzz_1, g_xxyyy_0_xxxxx_0, g_xxyyy_0_xxxxy_0, g_xxyyy_0_xxxxz_0, g_xxyyy_0_xxxyy_0, g_xxyyy_0_xxxyz_0, g_xxyyy_0_xxxzz_0, g_xxyyy_0_xxyyy_0, g_xxyyy_0_xxyyz_0, g_xxyyy_0_xxyzz_0, g_xxyyy_0_xxzzz_0, g_xxyyy_0_xyyyy_0, g_xxyyy_0_xyyyz_0, g_xxyyy_0_xyyzz_0, g_xxyyy_0_xyzzz_0, g_xxyyy_0_xzzzz_0, g_xxyyy_0_yyyyy_0, g_xxyyy_0_yyyyz_0, g_xxyyy_0_yyyzz_0, g_xxyyy_0_yyzzz_0, g_xxyyy_0_yzzzz_0, g_xxyyy_0_zzzzz_0, g_xyyy_0_xxxxy_1, g_xyyy_0_xxxy_1, g_xyyy_0_xxxyy_1, g_xyyy_0_xxxyz_1, g_xyyy_0_xxyy_1, g_xyyy_0_xxyyy_1, g_xyyy_0_xxyyz_1, g_xyyy_0_xxyz_1, g_xyyy_0_xxyzz_1, g_xyyy_0_xyyy_1, g_xyyy_0_xyyyy_1, g_xyyy_0_xyyyz_1, g_xyyy_0_xyyz_1, g_xyyy_0_xyyzz_1, g_xyyy_0_xyzz_1, g_xyyy_0_xyzzz_1, g_xyyy_0_yyyy_1, g_xyyy_0_yyyyy_1, g_xyyy_0_yyyyz_1, g_xyyy_0_yyyz_1, g_xyyy_0_yyyzz_1, g_xyyy_0_yyzz_1, g_xyyy_0_yyzzz_1, g_xyyy_0_yzzz_1, g_xyyy_0_yzzzz_1, g_xyyy_0_zzzzz_1, g_yyy_0_xxxxy_0, g_yyy_0_xxxxy_1, g_yyy_0_xxxyy_0, g_yyy_0_xxxyy_1, g_yyy_0_xxxyz_0, g_yyy_0_xxxyz_1, g_yyy_0_xxyyy_0, g_yyy_0_xxyyy_1, g_yyy_0_xxyyz_0, g_yyy_0_xxyyz_1, g_yyy_0_xxyzz_0, g_yyy_0_xxyzz_1, g_yyy_0_xyyyy_0, g_yyy_0_xyyyy_1, g_yyy_0_xyyyz_0, g_yyy_0_xyyyz_1, g_yyy_0_xyyzz_0, g_yyy_0_xyyzz_1, g_yyy_0_xyzzz_0, g_yyy_0_xyzzz_1, g_yyy_0_yyyyy_0, g_yyy_0_yyyyy_1, g_yyy_0_yyyyz_0, g_yyy_0_yyyyz_1, g_yyy_0_yyyzz_0, g_yyy_0_yyyzz_1, g_yyy_0_yyzzz_0, g_yyy_0_yyzzz_1, g_yyy_0_yzzzz_0, g_yyy_0_yzzzz_1, g_yyy_0_zzzzz_0, g_yyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyy_0_xxxxx_0[i] = 2.0 * g_xxy_0_xxxxx_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxx_1[i] * fz_be_0 + g_xxyy_0_xxxxx_1[i] * wa_y[i];

        g_xxyyy_0_xxxxy_0[i] = g_yyy_0_xxxxy_0[i] * fbe_0 - g_yyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxy_1[i] * fi_acd_0 + g_xyyy_0_xxxxy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxz_0[i] = 2.0 * g_xxy_0_xxxxz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxz_1[i] * fz_be_0 + g_xxyy_0_xxxxz_1[i] * wa_y[i];

        g_xxyyy_0_xxxyy_0[i] = g_yyy_0_xxxyy_0[i] * fbe_0 - g_yyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyy_1[i] * fi_acd_0 + g_xyyy_0_xxxyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxyz_0[i] = g_yyy_0_xxxyz_0[i] * fbe_0 - g_yyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyz_1[i] * fi_acd_0 + g_xyyy_0_xxxyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxzz_0[i] = 2.0 * g_xxy_0_xxxzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxzz_1[i] * fz_be_0 + g_xxyy_0_xxxzz_1[i] * wa_y[i];

        g_xxyyy_0_xxyyy_0[i] = g_yyy_0_xxyyy_0[i] * fbe_0 - g_yyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyy_1[i] * fi_acd_0 + g_xyyy_0_xxyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxyyz_0[i] = g_yyy_0_xxyyz_0[i] * fbe_0 - g_yyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyz_1[i] * fi_acd_0 + g_xyyy_0_xxyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxyzz_0[i] = g_yyy_0_xxyzz_0[i] * fbe_0 - g_yyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyzz_1[i] * fi_acd_0 + g_xyyy_0_xxyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxzzz_0[i] = 2.0 * g_xxy_0_xxzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxzzz_1[i] * fz_be_0 + g_xxyy_0_xxzzz_1[i] * wa_y[i];

        g_xxyyy_0_xyyyy_0[i] = g_yyy_0_xyyyy_0[i] * fbe_0 - g_yyy_0_xyyyy_1[i] * fz_be_0 + g_xyyy_0_yyyy_1[i] * fi_acd_0 + g_xyyy_0_xyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xyyyz_0[i] = g_yyy_0_xyyyz_0[i] * fbe_0 - g_yyy_0_xyyyz_1[i] * fz_be_0 + g_xyyy_0_yyyz_1[i] * fi_acd_0 + g_xyyy_0_xyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xyyzz_0[i] = g_yyy_0_xyyzz_0[i] * fbe_0 - g_yyy_0_xyyzz_1[i] * fz_be_0 + g_xyyy_0_yyzz_1[i] * fi_acd_0 + g_xyyy_0_xyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xyzzz_0[i] = g_yyy_0_xyzzz_0[i] * fbe_0 - g_yyy_0_xyzzz_1[i] * fz_be_0 + g_xyyy_0_yzzz_1[i] * fi_acd_0 + g_xyyy_0_xyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xzzzz_0[i] = 2.0 * g_xxy_0_xzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xzzzz_1[i] * fz_be_0 + g_xxyy_0_xzzzz_1[i] * wa_y[i];

        g_xxyyy_0_yyyyy_0[i] = g_yyy_0_yyyyy_0[i] * fbe_0 - g_yyy_0_yyyyy_1[i] * fz_be_0 + g_xyyy_0_yyyyy_1[i] * wa_x[i];

        g_xxyyy_0_yyyyz_0[i] = g_yyy_0_yyyyz_0[i] * fbe_0 - g_yyy_0_yyyyz_1[i] * fz_be_0 + g_xyyy_0_yyyyz_1[i] * wa_x[i];

        g_xxyyy_0_yyyzz_0[i] = g_yyy_0_yyyzz_0[i] * fbe_0 - g_yyy_0_yyyzz_1[i] * fz_be_0 + g_xyyy_0_yyyzz_1[i] * wa_x[i];

        g_xxyyy_0_yyzzz_0[i] = g_yyy_0_yyzzz_0[i] * fbe_0 - g_yyy_0_yyzzz_1[i] * fz_be_0 + g_xyyy_0_yyzzz_1[i] * wa_x[i];

        g_xxyyy_0_yzzzz_0[i] = g_yyy_0_yzzzz_0[i] * fbe_0 - g_yyy_0_yzzzz_1[i] * fz_be_0 + g_xyyy_0_yzzzz_1[i] * wa_x[i];

        g_xxyyy_0_zzzzz_0[i] = g_yyy_0_zzzzz_0[i] * fbe_0 - g_yyy_0_zzzzz_1[i] * fz_be_0 + g_xyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 147-168 components of targeted buffer : HSH

    auto g_xxyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 147);

    auto g_xxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 148);

    auto g_xxyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 149);

    auto g_xxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 150);

    auto g_xxyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 151);

    auto g_xxyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 152);

    auto g_xxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 153);

    auto g_xxyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 154);

    auto g_xxyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 155);

    auto g_xxyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 156);

    auto g_xxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 157);

    auto g_xxyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 158);

    auto g_xxyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 159);

    auto g_xxyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 160);

    auto g_xxyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 161);

    auto g_xxyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 162);

    auto g_xxyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 163);

    auto g_xxyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 164);

    auto g_xxyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 165);

    auto g_xxyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 166);

    auto g_xxyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 167);

    #pragma omp simd aligned(g_xxyy_0_xxxx_1, g_xxyy_0_xxxxx_1, g_xxyy_0_xxxxy_1, g_xxyy_0_xxxxz_1, g_xxyy_0_xxxy_1, g_xxyy_0_xxxyy_1, g_xxyy_0_xxxyz_1, g_xxyy_0_xxxz_1, g_xxyy_0_xxxzz_1, g_xxyy_0_xxyy_1, g_xxyy_0_xxyyy_1, g_xxyy_0_xxyyz_1, g_xxyy_0_xxyz_1, g_xxyy_0_xxyzz_1, g_xxyy_0_xxzz_1, g_xxyy_0_xxzzz_1, g_xxyy_0_xyyy_1, g_xxyy_0_xyyyy_1, g_xxyy_0_xyyyz_1, g_xxyy_0_xyyz_1, g_xxyy_0_xyyzz_1, g_xxyy_0_xyzz_1, g_xxyy_0_xyzzz_1, g_xxyy_0_xzzz_1, g_xxyy_0_xzzzz_1, g_xxyy_0_yyyy_1, g_xxyy_0_yyyyy_1, g_xxyy_0_yyyyz_1, g_xxyy_0_yyyz_1, g_xxyy_0_yyyzz_1, g_xxyy_0_yyzz_1, g_xxyy_0_yyzzz_1, g_xxyy_0_yzzz_1, g_xxyy_0_yzzzz_1, g_xxyy_0_zzzz_1, g_xxyy_0_zzzzz_1, g_xxyyz_0_xxxxx_0, g_xxyyz_0_xxxxy_0, g_xxyyz_0_xxxxz_0, g_xxyyz_0_xxxyy_0, g_xxyyz_0_xxxyz_0, g_xxyyz_0_xxxzz_0, g_xxyyz_0_xxyyy_0, g_xxyyz_0_xxyyz_0, g_xxyyz_0_xxyzz_0, g_xxyyz_0_xxzzz_0, g_xxyyz_0_xyyyy_0, g_xxyyz_0_xyyyz_0, g_xxyyz_0_xyyzz_0, g_xxyyz_0_xyzzz_0, g_xxyyz_0_xzzzz_0, g_xxyyz_0_yyyyy_0, g_xxyyz_0_yyyyz_0, g_xxyyz_0_yyyzz_0, g_xxyyz_0_yyzzz_0, g_xxyyz_0_yzzzz_0, g_xxyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyz_0_xxxxx_0[i] = g_xxyy_0_xxxxx_1[i] * wa_z[i];

        g_xxyyz_0_xxxxy_0[i] = g_xxyy_0_xxxxy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxz_0[i] = g_xxyy_0_xxxx_1[i] * fi_acd_0 + g_xxyy_0_xxxxz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyy_0[i] = g_xxyy_0_xxxyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxyz_0[i] = g_xxyy_0_xxxy_1[i] * fi_acd_0 + g_xxyy_0_xxxyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxzz_0[i] = 2.0 * g_xxyy_0_xxxz_1[i] * fi_acd_0 + g_xxyy_0_xxxzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyy_0[i] = g_xxyy_0_xxyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxyyz_0[i] = g_xxyy_0_xxyy_1[i] * fi_acd_0 + g_xxyy_0_xxyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxyzz_0[i] = 2.0 * g_xxyy_0_xxyz_1[i] * fi_acd_0 + g_xxyy_0_xxyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxzzz_0[i] = 3.0 * g_xxyy_0_xxzz_1[i] * fi_acd_0 + g_xxyy_0_xxzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyy_0[i] = g_xxyy_0_xyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xyyyz_0[i] = g_xxyy_0_xyyy_1[i] * fi_acd_0 + g_xxyy_0_xyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xyyzz_0[i] = 2.0 * g_xxyy_0_xyyz_1[i] * fi_acd_0 + g_xxyy_0_xyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xyzzz_0[i] = 3.0 * g_xxyy_0_xyzz_1[i] * fi_acd_0 + g_xxyy_0_xyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xzzzz_0[i] = 4.0 * g_xxyy_0_xzzz_1[i] * fi_acd_0 + g_xxyy_0_xzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyy_0[i] = g_xxyy_0_yyyyy_1[i] * wa_z[i];

        g_xxyyz_0_yyyyz_0[i] = g_xxyy_0_yyyy_1[i] * fi_acd_0 + g_xxyy_0_yyyyz_1[i] * wa_z[i];

        g_xxyyz_0_yyyzz_0[i] = 2.0 * g_xxyy_0_yyyz_1[i] * fi_acd_0 + g_xxyy_0_yyyzz_1[i] * wa_z[i];

        g_xxyyz_0_yyzzz_0[i] = 3.0 * g_xxyy_0_yyzz_1[i] * fi_acd_0 + g_xxyy_0_yyzzz_1[i] * wa_z[i];

        g_xxyyz_0_yzzzz_0[i] = 4.0 * g_xxyy_0_yzzz_1[i] * fi_acd_0 + g_xxyy_0_yzzzz_1[i] * wa_z[i];

        g_xxyyz_0_zzzzz_0[i] = 5.0 * g_xxyy_0_zzzz_1[i] * fi_acd_0 + g_xxyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 168-189 components of targeted buffer : HSH

    auto g_xxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 168);

    auto g_xxyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 169);

    auto g_xxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 170);

    auto g_xxyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 171);

    auto g_xxyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 172);

    auto g_xxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 173);

    auto g_xxyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 174);

    auto g_xxyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 175);

    auto g_xxyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 176);

    auto g_xxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 177);

    auto g_xxyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 178);

    auto g_xxyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 179);

    auto g_xxyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 180);

    auto g_xxyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 181);

    auto g_xxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 182);

    auto g_xxyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 183);

    auto g_xxyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 184);

    auto g_xxyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 185);

    auto g_xxyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 186);

    auto g_xxyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 187);

    auto g_xxyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 188);

    #pragma omp simd aligned(g_xxyzz_0_xxxxx_0, g_xxyzz_0_xxxxy_0, g_xxyzz_0_xxxxz_0, g_xxyzz_0_xxxyy_0, g_xxyzz_0_xxxyz_0, g_xxyzz_0_xxxzz_0, g_xxyzz_0_xxyyy_0, g_xxyzz_0_xxyyz_0, g_xxyzz_0_xxyzz_0, g_xxyzz_0_xxzzz_0, g_xxyzz_0_xyyyy_0, g_xxyzz_0_xyyyz_0, g_xxyzz_0_xyyzz_0, g_xxyzz_0_xyzzz_0, g_xxyzz_0_xzzzz_0, g_xxyzz_0_yyyyy_0, g_xxyzz_0_yyyyz_0, g_xxyzz_0_yyyzz_0, g_xxyzz_0_yyzzz_0, g_xxyzz_0_yzzzz_0, g_xxyzz_0_zzzzz_0, g_xxzz_0_xxxx_1, g_xxzz_0_xxxxx_1, g_xxzz_0_xxxxy_1, g_xxzz_0_xxxxz_1, g_xxzz_0_xxxy_1, g_xxzz_0_xxxyy_1, g_xxzz_0_xxxyz_1, g_xxzz_0_xxxz_1, g_xxzz_0_xxxzz_1, g_xxzz_0_xxyy_1, g_xxzz_0_xxyyy_1, g_xxzz_0_xxyyz_1, g_xxzz_0_xxyz_1, g_xxzz_0_xxyzz_1, g_xxzz_0_xxzz_1, g_xxzz_0_xxzzz_1, g_xxzz_0_xyyy_1, g_xxzz_0_xyyyy_1, g_xxzz_0_xyyyz_1, g_xxzz_0_xyyz_1, g_xxzz_0_xyyzz_1, g_xxzz_0_xyzz_1, g_xxzz_0_xyzzz_1, g_xxzz_0_xzzz_1, g_xxzz_0_xzzzz_1, g_xxzz_0_yyyy_1, g_xxzz_0_yyyyy_1, g_xxzz_0_yyyyz_1, g_xxzz_0_yyyz_1, g_xxzz_0_yyyzz_1, g_xxzz_0_yyzz_1, g_xxzz_0_yyzzz_1, g_xxzz_0_yzzz_1, g_xxzz_0_yzzzz_1, g_xxzz_0_zzzz_1, g_xxzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzz_0_xxxxx_0[i] = g_xxzz_0_xxxxx_1[i] * wa_y[i];

        g_xxyzz_0_xxxxy_0[i] = g_xxzz_0_xxxx_1[i] * fi_acd_0 + g_xxzz_0_xxxxy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxz_0[i] = g_xxzz_0_xxxxz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyy_0[i] = 2.0 * g_xxzz_0_xxxy_1[i] * fi_acd_0 + g_xxzz_0_xxxyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxyz_0[i] = g_xxzz_0_xxxz_1[i] * fi_acd_0 + g_xxzz_0_xxxyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxzz_0[i] = g_xxzz_0_xxxzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyy_0[i] = 3.0 * g_xxzz_0_xxyy_1[i] * fi_acd_0 + g_xxzz_0_xxyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxyyz_0[i] = 2.0 * g_xxzz_0_xxyz_1[i] * fi_acd_0 + g_xxzz_0_xxyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxyzz_0[i] = g_xxzz_0_xxzz_1[i] * fi_acd_0 + g_xxzz_0_xxyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxzzz_0[i] = g_xxzz_0_xxzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyy_0[i] = 4.0 * g_xxzz_0_xyyy_1[i] * fi_acd_0 + g_xxzz_0_xyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xyyyz_0[i] = 3.0 * g_xxzz_0_xyyz_1[i] * fi_acd_0 + g_xxzz_0_xyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xyyzz_0[i] = 2.0 * g_xxzz_0_xyzz_1[i] * fi_acd_0 + g_xxzz_0_xyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xyzzz_0[i] = g_xxzz_0_xzzz_1[i] * fi_acd_0 + g_xxzz_0_xyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xzzzz_0[i] = g_xxzz_0_xzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyy_0[i] = 5.0 * g_xxzz_0_yyyy_1[i] * fi_acd_0 + g_xxzz_0_yyyyy_1[i] * wa_y[i];

        g_xxyzz_0_yyyyz_0[i] = 4.0 * g_xxzz_0_yyyz_1[i] * fi_acd_0 + g_xxzz_0_yyyyz_1[i] * wa_y[i];

        g_xxyzz_0_yyyzz_0[i] = 3.0 * g_xxzz_0_yyzz_1[i] * fi_acd_0 + g_xxzz_0_yyyzz_1[i] * wa_y[i];

        g_xxyzz_0_yyzzz_0[i] = 2.0 * g_xxzz_0_yzzz_1[i] * fi_acd_0 + g_xxzz_0_yyzzz_1[i] * wa_y[i];

        g_xxyzz_0_yzzzz_0[i] = g_xxzz_0_zzzz_1[i] * fi_acd_0 + g_xxzz_0_yzzzz_1[i] * wa_y[i];

        g_xxyzz_0_zzzzz_0[i] = g_xxzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 189-210 components of targeted buffer : HSH

    auto g_xxzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 189);

    auto g_xxzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 190);

    auto g_xxzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 191);

    auto g_xxzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 192);

    auto g_xxzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 193);

    auto g_xxzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 194);

    auto g_xxzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 195);

    auto g_xxzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 196);

    auto g_xxzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 197);

    auto g_xxzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 198);

    auto g_xxzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 199);

    auto g_xxzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 200);

    auto g_xxzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 201);

    auto g_xxzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 202);

    auto g_xxzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 203);

    auto g_xxzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 204);

    auto g_xxzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 205);

    auto g_xxzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 206);

    auto g_xxzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 207);

    auto g_xxzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 208);

    auto g_xxzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 209);

    #pragma omp simd aligned(g_xxz_0_xxxxx_0, g_xxz_0_xxxxx_1, g_xxz_0_xxxxy_0, g_xxz_0_xxxxy_1, g_xxz_0_xxxyy_0, g_xxz_0_xxxyy_1, g_xxz_0_xxyyy_0, g_xxz_0_xxyyy_1, g_xxz_0_xyyyy_0, g_xxz_0_xyyyy_1, g_xxzz_0_xxxxx_1, g_xxzz_0_xxxxy_1, g_xxzz_0_xxxyy_1, g_xxzz_0_xxyyy_1, g_xxzz_0_xyyyy_1, g_xxzzz_0_xxxxx_0, g_xxzzz_0_xxxxy_0, g_xxzzz_0_xxxxz_0, g_xxzzz_0_xxxyy_0, g_xxzzz_0_xxxyz_0, g_xxzzz_0_xxxzz_0, g_xxzzz_0_xxyyy_0, g_xxzzz_0_xxyyz_0, g_xxzzz_0_xxyzz_0, g_xxzzz_0_xxzzz_0, g_xxzzz_0_xyyyy_0, g_xxzzz_0_xyyyz_0, g_xxzzz_0_xyyzz_0, g_xxzzz_0_xyzzz_0, g_xxzzz_0_xzzzz_0, g_xxzzz_0_yyyyy_0, g_xxzzz_0_yyyyz_0, g_xxzzz_0_yyyzz_0, g_xxzzz_0_yyzzz_0, g_xxzzz_0_yzzzz_0, g_xxzzz_0_zzzzz_0, g_xzzz_0_xxxxz_1, g_xzzz_0_xxxyz_1, g_xzzz_0_xxxz_1, g_xzzz_0_xxxzz_1, g_xzzz_0_xxyyz_1, g_xzzz_0_xxyz_1, g_xzzz_0_xxyzz_1, g_xzzz_0_xxzz_1, g_xzzz_0_xxzzz_1, g_xzzz_0_xyyyz_1, g_xzzz_0_xyyz_1, g_xzzz_0_xyyzz_1, g_xzzz_0_xyzz_1, g_xzzz_0_xyzzz_1, g_xzzz_0_xzzz_1, g_xzzz_0_xzzzz_1, g_xzzz_0_yyyyy_1, g_xzzz_0_yyyyz_1, g_xzzz_0_yyyz_1, g_xzzz_0_yyyzz_1, g_xzzz_0_yyzz_1, g_xzzz_0_yyzzz_1, g_xzzz_0_yzzz_1, g_xzzz_0_yzzzz_1, g_xzzz_0_zzzz_1, g_xzzz_0_zzzzz_1, g_zzz_0_xxxxz_0, g_zzz_0_xxxxz_1, g_zzz_0_xxxyz_0, g_zzz_0_xxxyz_1, g_zzz_0_xxxzz_0, g_zzz_0_xxxzz_1, g_zzz_0_xxyyz_0, g_zzz_0_xxyyz_1, g_zzz_0_xxyzz_0, g_zzz_0_xxyzz_1, g_zzz_0_xxzzz_0, g_zzz_0_xxzzz_1, g_zzz_0_xyyyz_0, g_zzz_0_xyyyz_1, g_zzz_0_xyyzz_0, g_zzz_0_xyyzz_1, g_zzz_0_xyzzz_0, g_zzz_0_xyzzz_1, g_zzz_0_xzzzz_0, g_zzz_0_xzzzz_1, g_zzz_0_yyyyy_0, g_zzz_0_yyyyy_1, g_zzz_0_yyyyz_0, g_zzz_0_yyyyz_1, g_zzz_0_yyyzz_0, g_zzz_0_yyyzz_1, g_zzz_0_yyzzz_0, g_zzz_0_yyzzz_1, g_zzz_0_yzzzz_0, g_zzz_0_yzzzz_1, g_zzz_0_zzzzz_0, g_zzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzz_0_xxxxx_0[i] = 2.0 * g_xxz_0_xxxxx_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxx_1[i] * fz_be_0 + g_xxzz_0_xxxxx_1[i] * wa_z[i];

        g_xxzzz_0_xxxxy_0[i] = 2.0 * g_xxz_0_xxxxy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxy_1[i] * fz_be_0 + g_xxzz_0_xxxxy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxz_0[i] = g_zzz_0_xxxxz_0[i] * fbe_0 - g_zzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxz_1[i] * fi_acd_0 + g_xzzz_0_xxxxz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyy_0[i] = 2.0 * g_xxz_0_xxxyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxyy_1[i] * fz_be_0 + g_xxzz_0_xxxyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxyz_0[i] = g_zzz_0_xxxyz_0[i] * fbe_0 - g_zzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyz_1[i] * fi_acd_0 + g_xzzz_0_xxxyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxzz_0[i] = g_zzz_0_xxxzz_0[i] * fbe_0 - g_zzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxzz_1[i] * fi_acd_0 + g_xzzz_0_xxxzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyy_0[i] = 2.0 * g_xxz_0_xxyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxyyy_1[i] * fz_be_0 + g_xxzz_0_xxyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxyyz_0[i] = g_zzz_0_xxyyz_0[i] * fbe_0 - g_zzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyz_1[i] * fi_acd_0 + g_xzzz_0_xxyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxyzz_0[i] = g_zzz_0_xxyzz_0[i] * fbe_0 - g_zzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyzz_1[i] * fi_acd_0 + g_xzzz_0_xxyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxzzz_0[i] = g_zzz_0_xxzzz_0[i] * fbe_0 - g_zzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xzzz_1[i] * fi_acd_0 + g_xzzz_0_xxzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyy_0[i] = 2.0 * g_xxz_0_xyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xyyyy_1[i] * fz_be_0 + g_xxzz_0_xyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xyyyz_0[i] = g_zzz_0_xyyyz_0[i] * fbe_0 - g_zzz_0_xyyyz_1[i] * fz_be_0 + g_xzzz_0_yyyz_1[i] * fi_acd_0 + g_xzzz_0_xyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xyyzz_0[i] = g_zzz_0_xyyzz_0[i] * fbe_0 - g_zzz_0_xyyzz_1[i] * fz_be_0 + g_xzzz_0_yyzz_1[i] * fi_acd_0 + g_xzzz_0_xyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xyzzz_0[i] = g_zzz_0_xyzzz_0[i] * fbe_0 - g_zzz_0_xyzzz_1[i] * fz_be_0 + g_xzzz_0_yzzz_1[i] * fi_acd_0 + g_xzzz_0_xyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xzzzz_0[i] = g_zzz_0_xzzzz_0[i] * fbe_0 - g_zzz_0_xzzzz_1[i] * fz_be_0 + g_xzzz_0_zzzz_1[i] * fi_acd_0 + g_xzzz_0_xzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyy_0[i] = g_zzz_0_yyyyy_0[i] * fbe_0 - g_zzz_0_yyyyy_1[i] * fz_be_0 + g_xzzz_0_yyyyy_1[i] * wa_x[i];

        g_xxzzz_0_yyyyz_0[i] = g_zzz_0_yyyyz_0[i] * fbe_0 - g_zzz_0_yyyyz_1[i] * fz_be_0 + g_xzzz_0_yyyyz_1[i] * wa_x[i];

        g_xxzzz_0_yyyzz_0[i] = g_zzz_0_yyyzz_0[i] * fbe_0 - g_zzz_0_yyyzz_1[i] * fz_be_0 + g_xzzz_0_yyyzz_1[i] * wa_x[i];

        g_xxzzz_0_yyzzz_0[i] = g_zzz_0_yyzzz_0[i] * fbe_0 - g_zzz_0_yyzzz_1[i] * fz_be_0 + g_xzzz_0_yyzzz_1[i] * wa_x[i];

        g_xxzzz_0_yzzzz_0[i] = g_zzz_0_yzzzz_0[i] * fbe_0 - g_zzz_0_yzzzz_1[i] * fz_be_0 + g_xzzz_0_yzzzz_1[i] * wa_x[i];

        g_xxzzz_0_zzzzz_0[i] = g_zzz_0_zzzzz_0[i] * fbe_0 - g_zzz_0_zzzzz_1[i] * fz_be_0 + g_xzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 210-231 components of targeted buffer : HSH

    auto g_xyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 210);

    auto g_xyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 211);

    auto g_xyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 212);

    auto g_xyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 213);

    auto g_xyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 214);

    auto g_xyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 215);

    auto g_xyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 216);

    auto g_xyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 217);

    auto g_xyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 218);

    auto g_xyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 219);

    auto g_xyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 220);

    auto g_xyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 221);

    auto g_xyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 222);

    auto g_xyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 223);

    auto g_xyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 224);

    auto g_xyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 225);

    auto g_xyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 226);

    auto g_xyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 227);

    auto g_xyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 228);

    auto g_xyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 229);

    auto g_xyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 230);

    #pragma omp simd aligned(g_xyyyy_0_xxxxx_0, g_xyyyy_0_xxxxy_0, g_xyyyy_0_xxxxz_0, g_xyyyy_0_xxxyy_0, g_xyyyy_0_xxxyz_0, g_xyyyy_0_xxxzz_0, g_xyyyy_0_xxyyy_0, g_xyyyy_0_xxyyz_0, g_xyyyy_0_xxyzz_0, g_xyyyy_0_xxzzz_0, g_xyyyy_0_xyyyy_0, g_xyyyy_0_xyyyz_0, g_xyyyy_0_xyyzz_0, g_xyyyy_0_xyzzz_0, g_xyyyy_0_xzzzz_0, g_xyyyy_0_yyyyy_0, g_xyyyy_0_yyyyz_0, g_xyyyy_0_yyyzz_0, g_xyyyy_0_yyzzz_0, g_xyyyy_0_yzzzz_0, g_xyyyy_0_zzzzz_0, g_yyyy_0_xxxx_1, g_yyyy_0_xxxxx_1, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxxz_1, g_yyyy_0_xxxy_1, g_yyyy_0_xxxyy_1, g_yyyy_0_xxxyz_1, g_yyyy_0_xxxz_1, g_yyyy_0_xxxzz_1, g_yyyy_0_xxyy_1, g_yyyy_0_xxyyy_1, g_yyyy_0_xxyyz_1, g_yyyy_0_xxyz_1, g_yyyy_0_xxyzz_1, g_yyyy_0_xxzz_1, g_yyyy_0_xxzzz_1, g_yyyy_0_xyyy_1, g_yyyy_0_xyyyy_1, g_yyyy_0_xyyyz_1, g_yyyy_0_xyyz_1, g_yyyy_0_xyyzz_1, g_yyyy_0_xyzz_1, g_yyyy_0_xyzzz_1, g_yyyy_0_xzzz_1, g_yyyy_0_xzzzz_1, g_yyyy_0_yyyy_1, g_yyyy_0_yyyyy_1, g_yyyy_0_yyyyz_1, g_yyyy_0_yyyz_1, g_yyyy_0_yyyzz_1, g_yyyy_0_yyzz_1, g_yyyy_0_yyzzz_1, g_yyyy_0_yzzz_1, g_yyyy_0_yzzzz_1, g_yyyy_0_zzzz_1, g_yyyy_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyy_0_xxxxx_0[i] = 5.0 * g_yyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxx_1[i] * wa_x[i];

        g_xyyyy_0_xxxxy_0[i] = 4.0 * g_yyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxz_0[i] = 4.0 * g_yyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyy_0[i] = 3.0 * g_yyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxyz_0[i] = 3.0 * g_yyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxzz_0[i] = 3.0 * g_yyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyy_0[i] = 2.0 * g_yyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxyyz_0[i] = 2.0 * g_yyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxyzz_0[i] = 2.0 * g_yyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxzzz_0[i] = 2.0 * g_yyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyy_0_xxzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyy_0[i] = g_yyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xyyyz_0[i] = g_yyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xyyzz_0[i] = g_yyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xyzzz_0[i] = g_yyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xzzzz_0[i] = g_yyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyy_0_xzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyy_0[i] = g_yyyy_0_yyyyy_1[i] * wa_x[i];

        g_xyyyy_0_yyyyz_0[i] = g_yyyy_0_yyyyz_1[i] * wa_x[i];

        g_xyyyy_0_yyyzz_0[i] = g_yyyy_0_yyyzz_1[i] * wa_x[i];

        g_xyyyy_0_yyzzz_0[i] = g_yyyy_0_yyzzz_1[i] * wa_x[i];

        g_xyyyy_0_yzzzz_0[i] = g_yyyy_0_yzzzz_1[i] * wa_x[i];

        g_xyyyy_0_zzzzz_0[i] = g_yyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 231-252 components of targeted buffer : HSH

    auto g_xyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 231);

    auto g_xyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 232);

    auto g_xyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 233);

    auto g_xyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 234);

    auto g_xyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 235);

    auto g_xyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 236);

    auto g_xyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 237);

    auto g_xyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 238);

    auto g_xyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 239);

    auto g_xyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 240);

    auto g_xyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 241);

    auto g_xyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 242);

    auto g_xyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 243);

    auto g_xyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 244);

    auto g_xyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 245);

    auto g_xyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 246);

    auto g_xyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 247);

    auto g_xyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 248);

    auto g_xyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 249);

    auto g_xyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 250);

    auto g_xyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 251);

    #pragma omp simd aligned(g_xyyy_0_xxxxx_1, g_xyyy_0_xxxxy_1, g_xyyy_0_xxxyy_1, g_xyyy_0_xxyyy_1, g_xyyy_0_xyyyy_1, g_xyyyz_0_xxxxx_0, g_xyyyz_0_xxxxy_0, g_xyyyz_0_xxxxz_0, g_xyyyz_0_xxxyy_0, g_xyyyz_0_xxxyz_0, g_xyyyz_0_xxxzz_0, g_xyyyz_0_xxyyy_0, g_xyyyz_0_xxyyz_0, g_xyyyz_0_xxyzz_0, g_xyyyz_0_xxzzz_0, g_xyyyz_0_xyyyy_0, g_xyyyz_0_xyyyz_0, g_xyyyz_0_xyyzz_0, g_xyyyz_0_xyzzz_0, g_xyyyz_0_xzzzz_0, g_xyyyz_0_yyyyy_0, g_xyyyz_0_yyyyz_0, g_xyyyz_0_yyyzz_0, g_xyyyz_0_yyzzz_0, g_xyyyz_0_yzzzz_0, g_xyyyz_0_zzzzz_0, g_yyyz_0_xxxxz_1, g_yyyz_0_xxxyz_1, g_yyyz_0_xxxz_1, g_yyyz_0_xxxzz_1, g_yyyz_0_xxyyz_1, g_yyyz_0_xxyz_1, g_yyyz_0_xxyzz_1, g_yyyz_0_xxzz_1, g_yyyz_0_xxzzz_1, g_yyyz_0_xyyyz_1, g_yyyz_0_xyyz_1, g_yyyz_0_xyyzz_1, g_yyyz_0_xyzz_1, g_yyyz_0_xyzzz_1, g_yyyz_0_xzzz_1, g_yyyz_0_xzzzz_1, g_yyyz_0_yyyyy_1, g_yyyz_0_yyyyz_1, g_yyyz_0_yyyz_1, g_yyyz_0_yyyzz_1, g_yyyz_0_yyzz_1, g_yyyz_0_yyzzz_1, g_yyyz_0_yzzz_1, g_yyyz_0_yzzzz_1, g_yyyz_0_zzzz_1, g_yyyz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyz_0_xxxxx_0[i] = g_xyyy_0_xxxxx_1[i] * wa_z[i];

        g_xyyyz_0_xxxxy_0[i] = g_xyyy_0_xxxxy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxz_0[i] = 4.0 * g_yyyz_0_xxxz_1[i] * fi_acd_0 + g_yyyz_0_xxxxz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyy_0[i] = g_xyyy_0_xxxyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxyz_0[i] = 3.0 * g_yyyz_0_xxyz_1[i] * fi_acd_0 + g_yyyz_0_xxxyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxzz_0[i] = 3.0 * g_yyyz_0_xxzz_1[i] * fi_acd_0 + g_yyyz_0_xxxzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyy_0[i] = g_xyyy_0_xxyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxyyz_0[i] = 2.0 * g_yyyz_0_xyyz_1[i] * fi_acd_0 + g_yyyz_0_xxyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxyzz_0[i] = 2.0 * g_yyyz_0_xyzz_1[i] * fi_acd_0 + g_yyyz_0_xxyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxzzz_0[i] = 2.0 * g_yyyz_0_xzzz_1[i] * fi_acd_0 + g_yyyz_0_xxzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyy_0[i] = g_xyyy_0_xyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xyyyz_0[i] = g_yyyz_0_yyyz_1[i] * fi_acd_0 + g_yyyz_0_xyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xyyzz_0[i] = g_yyyz_0_yyzz_1[i] * fi_acd_0 + g_yyyz_0_xyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xyzzz_0[i] = g_yyyz_0_yzzz_1[i] * fi_acd_0 + g_yyyz_0_xyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xzzzz_0[i] = g_yyyz_0_zzzz_1[i] * fi_acd_0 + g_yyyz_0_xzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyy_0[i] = g_yyyz_0_yyyyy_1[i] * wa_x[i];

        g_xyyyz_0_yyyyz_0[i] = g_yyyz_0_yyyyz_1[i] * wa_x[i];

        g_xyyyz_0_yyyzz_0[i] = g_yyyz_0_yyyzz_1[i] * wa_x[i];

        g_xyyyz_0_yyzzz_0[i] = g_yyyz_0_yyzzz_1[i] * wa_x[i];

        g_xyyyz_0_yzzzz_0[i] = g_yyyz_0_yzzzz_1[i] * wa_x[i];

        g_xyyyz_0_zzzzz_0[i] = g_yyyz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 252-273 components of targeted buffer : HSH

    auto g_xyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 252);

    auto g_xyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 253);

    auto g_xyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 254);

    auto g_xyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 255);

    auto g_xyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 256);

    auto g_xyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 257);

    auto g_xyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 258);

    auto g_xyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 259);

    auto g_xyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 260);

    auto g_xyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 261);

    auto g_xyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 262);

    auto g_xyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 263);

    auto g_xyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 264);

    auto g_xyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 265);

    auto g_xyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 266);

    auto g_xyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 267);

    auto g_xyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 268);

    auto g_xyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 269);

    auto g_xyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 270);

    auto g_xyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 271);

    auto g_xyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 272);

    #pragma omp simd aligned(g_xyyzz_0_xxxxx_0, g_xyyzz_0_xxxxy_0, g_xyyzz_0_xxxxz_0, g_xyyzz_0_xxxyy_0, g_xyyzz_0_xxxyz_0, g_xyyzz_0_xxxzz_0, g_xyyzz_0_xxyyy_0, g_xyyzz_0_xxyyz_0, g_xyyzz_0_xxyzz_0, g_xyyzz_0_xxzzz_0, g_xyyzz_0_xyyyy_0, g_xyyzz_0_xyyyz_0, g_xyyzz_0_xyyzz_0, g_xyyzz_0_xyzzz_0, g_xyyzz_0_xzzzz_0, g_xyyzz_0_yyyyy_0, g_xyyzz_0_yyyyz_0, g_xyyzz_0_yyyzz_0, g_xyyzz_0_yyzzz_0, g_xyyzz_0_yzzzz_0, g_xyyzz_0_zzzzz_0, g_yyzz_0_xxxx_1, g_yyzz_0_xxxxx_1, g_yyzz_0_xxxxy_1, g_yyzz_0_xxxxz_1, g_yyzz_0_xxxy_1, g_yyzz_0_xxxyy_1, g_yyzz_0_xxxyz_1, g_yyzz_0_xxxz_1, g_yyzz_0_xxxzz_1, g_yyzz_0_xxyy_1, g_yyzz_0_xxyyy_1, g_yyzz_0_xxyyz_1, g_yyzz_0_xxyz_1, g_yyzz_0_xxyzz_1, g_yyzz_0_xxzz_1, g_yyzz_0_xxzzz_1, g_yyzz_0_xyyy_1, g_yyzz_0_xyyyy_1, g_yyzz_0_xyyyz_1, g_yyzz_0_xyyz_1, g_yyzz_0_xyyzz_1, g_yyzz_0_xyzz_1, g_yyzz_0_xyzzz_1, g_yyzz_0_xzzz_1, g_yyzz_0_xzzzz_1, g_yyzz_0_yyyy_1, g_yyzz_0_yyyyy_1, g_yyzz_0_yyyyz_1, g_yyzz_0_yyyz_1, g_yyzz_0_yyyzz_1, g_yyzz_0_yyzz_1, g_yyzz_0_yyzzz_1, g_yyzz_0_yzzz_1, g_yyzz_0_yzzzz_1, g_yyzz_0_zzzz_1, g_yyzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzz_0_xxxxx_0[i] = 5.0 * g_yyzz_0_xxxx_1[i] * fi_acd_0 + g_yyzz_0_xxxxx_1[i] * wa_x[i];

        g_xyyzz_0_xxxxy_0[i] = 4.0 * g_yyzz_0_xxxy_1[i] * fi_acd_0 + g_yyzz_0_xxxxy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxz_0[i] = 4.0 * g_yyzz_0_xxxz_1[i] * fi_acd_0 + g_yyzz_0_xxxxz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyy_0[i] = 3.0 * g_yyzz_0_xxyy_1[i] * fi_acd_0 + g_yyzz_0_xxxyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxyz_0[i] = 3.0 * g_yyzz_0_xxyz_1[i] * fi_acd_0 + g_yyzz_0_xxxyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxzz_0[i] = 3.0 * g_yyzz_0_xxzz_1[i] * fi_acd_0 + g_yyzz_0_xxxzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyy_0[i] = 2.0 * g_yyzz_0_xyyy_1[i] * fi_acd_0 + g_yyzz_0_xxyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxyyz_0[i] = 2.0 * g_yyzz_0_xyyz_1[i] * fi_acd_0 + g_yyzz_0_xxyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxyzz_0[i] = 2.0 * g_yyzz_0_xyzz_1[i] * fi_acd_0 + g_yyzz_0_xxyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxzzz_0[i] = 2.0 * g_yyzz_0_xzzz_1[i] * fi_acd_0 + g_yyzz_0_xxzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyy_0[i] = g_yyzz_0_yyyy_1[i] * fi_acd_0 + g_yyzz_0_xyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xyyyz_0[i] = g_yyzz_0_yyyz_1[i] * fi_acd_0 + g_yyzz_0_xyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xyyzz_0[i] = g_yyzz_0_yyzz_1[i] * fi_acd_0 + g_yyzz_0_xyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xyzzz_0[i] = g_yyzz_0_yzzz_1[i] * fi_acd_0 + g_yyzz_0_xyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xzzzz_0[i] = g_yyzz_0_zzzz_1[i] * fi_acd_0 + g_yyzz_0_xzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyy_0[i] = g_yyzz_0_yyyyy_1[i] * wa_x[i];

        g_xyyzz_0_yyyyz_0[i] = g_yyzz_0_yyyyz_1[i] * wa_x[i];

        g_xyyzz_0_yyyzz_0[i] = g_yyzz_0_yyyzz_1[i] * wa_x[i];

        g_xyyzz_0_yyzzz_0[i] = g_yyzz_0_yyzzz_1[i] * wa_x[i];

        g_xyyzz_0_yzzzz_0[i] = g_yyzz_0_yzzzz_1[i] * wa_x[i];

        g_xyyzz_0_zzzzz_0[i] = g_yyzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 273-294 components of targeted buffer : HSH

    auto g_xyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 273);

    auto g_xyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 274);

    auto g_xyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 275);

    auto g_xyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 276);

    auto g_xyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 277);

    auto g_xyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 278);

    auto g_xyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 279);

    auto g_xyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 280);

    auto g_xyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 281);

    auto g_xyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 282);

    auto g_xyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 283);

    auto g_xyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 284);

    auto g_xyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 285);

    auto g_xyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 286);

    auto g_xyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 287);

    auto g_xyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 288);

    auto g_xyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 289);

    auto g_xyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 290);

    auto g_xyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 291);

    auto g_xyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 292);

    auto g_xyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 293);

    #pragma omp simd aligned(g_xyzzz_0_xxxxx_0, g_xyzzz_0_xxxxy_0, g_xyzzz_0_xxxxz_0, g_xyzzz_0_xxxyy_0, g_xyzzz_0_xxxyz_0, g_xyzzz_0_xxxzz_0, g_xyzzz_0_xxyyy_0, g_xyzzz_0_xxyyz_0, g_xyzzz_0_xxyzz_0, g_xyzzz_0_xxzzz_0, g_xyzzz_0_xyyyy_0, g_xyzzz_0_xyyyz_0, g_xyzzz_0_xyyzz_0, g_xyzzz_0_xyzzz_0, g_xyzzz_0_xzzzz_0, g_xyzzz_0_yyyyy_0, g_xyzzz_0_yyyyz_0, g_xyzzz_0_yyyzz_0, g_xyzzz_0_yyzzz_0, g_xyzzz_0_yzzzz_0, g_xyzzz_0_zzzzz_0, g_xzzz_0_xxxxx_1, g_xzzz_0_xxxxz_1, g_xzzz_0_xxxzz_1, g_xzzz_0_xxzzz_1, g_xzzz_0_xzzzz_1, g_yzzz_0_xxxxy_1, g_yzzz_0_xxxy_1, g_yzzz_0_xxxyy_1, g_yzzz_0_xxxyz_1, g_yzzz_0_xxyy_1, g_yzzz_0_xxyyy_1, g_yzzz_0_xxyyz_1, g_yzzz_0_xxyz_1, g_yzzz_0_xxyzz_1, g_yzzz_0_xyyy_1, g_yzzz_0_xyyyy_1, g_yzzz_0_xyyyz_1, g_yzzz_0_xyyz_1, g_yzzz_0_xyyzz_1, g_yzzz_0_xyzz_1, g_yzzz_0_xyzzz_1, g_yzzz_0_yyyy_1, g_yzzz_0_yyyyy_1, g_yzzz_0_yyyyz_1, g_yzzz_0_yyyz_1, g_yzzz_0_yyyzz_1, g_yzzz_0_yyzz_1, g_yzzz_0_yyzzz_1, g_yzzz_0_yzzz_1, g_yzzz_0_yzzzz_1, g_yzzz_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzz_0_xxxxx_0[i] = g_xzzz_0_xxxxx_1[i] * wa_y[i];

        g_xyzzz_0_xxxxy_0[i] = 4.0 * g_yzzz_0_xxxy_1[i] * fi_acd_0 + g_yzzz_0_xxxxy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxz_0[i] = g_xzzz_0_xxxxz_1[i] * wa_y[i];

        g_xyzzz_0_xxxyy_0[i] = 3.0 * g_yzzz_0_xxyy_1[i] * fi_acd_0 + g_yzzz_0_xxxyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxyz_0[i] = 3.0 * g_yzzz_0_xxyz_1[i] * fi_acd_0 + g_yzzz_0_xxxyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxzz_0[i] = g_xzzz_0_xxxzz_1[i] * wa_y[i];

        g_xyzzz_0_xxyyy_0[i] = 2.0 * g_yzzz_0_xyyy_1[i] * fi_acd_0 + g_yzzz_0_xxyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxyyz_0[i] = 2.0 * g_yzzz_0_xyyz_1[i] * fi_acd_0 + g_yzzz_0_xxyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxyzz_0[i] = 2.0 * g_yzzz_0_xyzz_1[i] * fi_acd_0 + g_yzzz_0_xxyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxzzz_0[i] = g_xzzz_0_xxzzz_1[i] * wa_y[i];

        g_xyzzz_0_xyyyy_0[i] = g_yzzz_0_yyyy_1[i] * fi_acd_0 + g_yzzz_0_xyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xyyyz_0[i] = g_yzzz_0_yyyz_1[i] * fi_acd_0 + g_yzzz_0_xyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xyyzz_0[i] = g_yzzz_0_yyzz_1[i] * fi_acd_0 + g_yzzz_0_xyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xyzzz_0[i] = g_yzzz_0_yzzz_1[i] * fi_acd_0 + g_yzzz_0_xyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xzzzz_0[i] = g_xzzz_0_xzzzz_1[i] * wa_y[i];

        g_xyzzz_0_yyyyy_0[i] = g_yzzz_0_yyyyy_1[i] * wa_x[i];

        g_xyzzz_0_yyyyz_0[i] = g_yzzz_0_yyyyz_1[i] * wa_x[i];

        g_xyzzz_0_yyyzz_0[i] = g_yzzz_0_yyyzz_1[i] * wa_x[i];

        g_xyzzz_0_yyzzz_0[i] = g_yzzz_0_yyzzz_1[i] * wa_x[i];

        g_xyzzz_0_yzzzz_0[i] = g_yzzz_0_yzzzz_1[i] * wa_x[i];

        g_xyzzz_0_zzzzz_0[i] = g_yzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 294-315 components of targeted buffer : HSH

    auto g_xzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 294);

    auto g_xzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 295);

    auto g_xzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 296);

    auto g_xzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 297);

    auto g_xzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 298);

    auto g_xzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 299);

    auto g_xzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 300);

    auto g_xzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 301);

    auto g_xzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 302);

    auto g_xzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 303);

    auto g_xzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 304);

    auto g_xzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 305);

    auto g_xzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 306);

    auto g_xzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 307);

    auto g_xzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 308);

    auto g_xzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 309);

    auto g_xzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 310);

    auto g_xzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 311);

    auto g_xzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 312);

    auto g_xzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 313);

    auto g_xzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 314);

    #pragma omp simd aligned(g_xzzzz_0_xxxxx_0, g_xzzzz_0_xxxxy_0, g_xzzzz_0_xxxxz_0, g_xzzzz_0_xxxyy_0, g_xzzzz_0_xxxyz_0, g_xzzzz_0_xxxzz_0, g_xzzzz_0_xxyyy_0, g_xzzzz_0_xxyyz_0, g_xzzzz_0_xxyzz_0, g_xzzzz_0_xxzzz_0, g_xzzzz_0_xyyyy_0, g_xzzzz_0_xyyyz_0, g_xzzzz_0_xyyzz_0, g_xzzzz_0_xyzzz_0, g_xzzzz_0_xzzzz_0, g_xzzzz_0_yyyyy_0, g_xzzzz_0_yyyyz_0, g_xzzzz_0_yyyzz_0, g_xzzzz_0_yyzzz_0, g_xzzzz_0_yzzzz_0, g_xzzzz_0_zzzzz_0, g_zzzz_0_xxxx_1, g_zzzz_0_xxxxx_1, g_zzzz_0_xxxxy_1, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxy_1, g_zzzz_0_xxxyy_1, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxz_1, g_zzzz_0_xxxzz_1, g_zzzz_0_xxyy_1, g_zzzz_0_xxyyy_1, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyz_1, g_zzzz_0_xxyzz_1, g_zzzz_0_xxzz_1, g_zzzz_0_xxzzz_1, g_zzzz_0_xyyy_1, g_zzzz_0_xyyyy_1, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyz_1, g_zzzz_0_xyyzz_1, g_zzzz_0_xyzz_1, g_zzzz_0_xyzzz_1, g_zzzz_0_xzzz_1, g_zzzz_0_xzzzz_1, g_zzzz_0_yyyy_1, g_zzzz_0_yyyyy_1, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyz_1, g_zzzz_0_yyyzz_1, g_zzzz_0_yyzz_1, g_zzzz_0_yyzzz_1, g_zzzz_0_yzzz_1, g_zzzz_0_yzzzz_1, g_zzzz_0_zzzz_1, g_zzzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzz_0_xxxxx_0[i] = 5.0 * g_zzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxx_1[i] * wa_x[i];

        g_xzzzz_0_xxxxy_0[i] = 4.0 * g_zzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxz_0[i] = 4.0 * g_zzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyy_0[i] = 3.0 * g_zzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxyz_0[i] = 3.0 * g_zzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxzz_0[i] = 3.0 * g_zzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyy_0[i] = 2.0 * g_zzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxyyz_0[i] = 2.0 * g_zzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxyzz_0[i] = 2.0 * g_zzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxzzz_0[i] = 2.0 * g_zzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzz_0_xxzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyy_0[i] = g_zzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xyyyz_0[i] = g_zzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xyyzz_0[i] = g_zzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xyzzz_0[i] = g_zzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xzzzz_0[i] = g_zzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzz_0_xzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyy_0[i] = g_zzzz_0_yyyyy_1[i] * wa_x[i];

        g_xzzzz_0_yyyyz_0[i] = g_zzzz_0_yyyyz_1[i] * wa_x[i];

        g_xzzzz_0_yyyzz_0[i] = g_zzzz_0_yyyzz_1[i] * wa_x[i];

        g_xzzzz_0_yyzzz_0[i] = g_zzzz_0_yyzzz_1[i] * wa_x[i];

        g_xzzzz_0_yzzzz_0[i] = g_zzzz_0_yzzzz_1[i] * wa_x[i];

        g_xzzzz_0_zzzzz_0[i] = g_zzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 315-336 components of targeted buffer : HSH

    auto g_yyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 315);

    auto g_yyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 316);

    auto g_yyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 317);

    auto g_yyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 318);

    auto g_yyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 319);

    auto g_yyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 320);

    auto g_yyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 321);

    auto g_yyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 322);

    auto g_yyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 323);

    auto g_yyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 324);

    auto g_yyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 325);

    auto g_yyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 326);

    auto g_yyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 327);

    auto g_yyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 328);

    auto g_yyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 329);

    auto g_yyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 330);

    auto g_yyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 331);

    auto g_yyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 332);

    auto g_yyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 333);

    auto g_yyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 334);

    auto g_yyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 335);

    #pragma omp simd aligned(g_yyy_0_xxxxx_0, g_yyy_0_xxxxx_1, g_yyy_0_xxxxy_0, g_yyy_0_xxxxy_1, g_yyy_0_xxxxz_0, g_yyy_0_xxxxz_1, g_yyy_0_xxxyy_0, g_yyy_0_xxxyy_1, g_yyy_0_xxxyz_0, g_yyy_0_xxxyz_1, g_yyy_0_xxxzz_0, g_yyy_0_xxxzz_1, g_yyy_0_xxyyy_0, g_yyy_0_xxyyy_1, g_yyy_0_xxyyz_0, g_yyy_0_xxyyz_1, g_yyy_0_xxyzz_0, g_yyy_0_xxyzz_1, g_yyy_0_xxzzz_0, g_yyy_0_xxzzz_1, g_yyy_0_xyyyy_0, g_yyy_0_xyyyy_1, g_yyy_0_xyyyz_0, g_yyy_0_xyyyz_1, g_yyy_0_xyyzz_0, g_yyy_0_xyyzz_1, g_yyy_0_xyzzz_0, g_yyy_0_xyzzz_1, g_yyy_0_xzzzz_0, g_yyy_0_xzzzz_1, g_yyy_0_yyyyy_0, g_yyy_0_yyyyy_1, g_yyy_0_yyyyz_0, g_yyy_0_yyyyz_1, g_yyy_0_yyyzz_0, g_yyy_0_yyyzz_1, g_yyy_0_yyzzz_0, g_yyy_0_yyzzz_1, g_yyy_0_yzzzz_0, g_yyy_0_yzzzz_1, g_yyy_0_zzzzz_0, g_yyy_0_zzzzz_1, g_yyyy_0_xxxx_1, g_yyyy_0_xxxxx_1, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxxz_1, g_yyyy_0_xxxy_1, g_yyyy_0_xxxyy_1, g_yyyy_0_xxxyz_1, g_yyyy_0_xxxz_1, g_yyyy_0_xxxzz_1, g_yyyy_0_xxyy_1, g_yyyy_0_xxyyy_1, g_yyyy_0_xxyyz_1, g_yyyy_0_xxyz_1, g_yyyy_0_xxyzz_1, g_yyyy_0_xxzz_1, g_yyyy_0_xxzzz_1, g_yyyy_0_xyyy_1, g_yyyy_0_xyyyy_1, g_yyyy_0_xyyyz_1, g_yyyy_0_xyyz_1, g_yyyy_0_xyyzz_1, g_yyyy_0_xyzz_1, g_yyyy_0_xyzzz_1, g_yyyy_0_xzzz_1, g_yyyy_0_xzzzz_1, g_yyyy_0_yyyy_1, g_yyyy_0_yyyyy_1, g_yyyy_0_yyyyz_1, g_yyyy_0_yyyz_1, g_yyyy_0_yyyzz_1, g_yyyy_0_yyzz_1, g_yyyy_0_yyzzz_1, g_yyyy_0_yzzz_1, g_yyyy_0_yzzzz_1, g_yyyy_0_zzzz_1, g_yyyy_0_zzzzz_1, g_yyyyy_0_xxxxx_0, g_yyyyy_0_xxxxy_0, g_yyyyy_0_xxxxz_0, g_yyyyy_0_xxxyy_0, g_yyyyy_0_xxxyz_0, g_yyyyy_0_xxxzz_0, g_yyyyy_0_xxyyy_0, g_yyyyy_0_xxyyz_0, g_yyyyy_0_xxyzz_0, g_yyyyy_0_xxzzz_0, g_yyyyy_0_xyyyy_0, g_yyyyy_0_xyyyz_0, g_yyyyy_0_xyyzz_0, g_yyyyy_0_xyzzz_0, g_yyyyy_0_xzzzz_0, g_yyyyy_0_yyyyy_0, g_yyyyy_0_yyyyz_0, g_yyyyy_0_yyyzz_0, g_yyyyy_0_yyzzz_0, g_yyyyy_0_yzzzz_0, g_yyyyy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyy_0_xxxxx_0[i] = 4.0 * g_yyy_0_xxxxx_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxx_1[i] * fz_be_0 + g_yyyy_0_xxxxx_1[i] * wa_y[i];

        g_yyyyy_0_xxxxy_0[i] = 4.0 * g_yyy_0_xxxxy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxy_1[i] * fz_be_0 + g_yyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxz_0[i] = 4.0 * g_yyy_0_xxxxz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxz_1[i] * fz_be_0 + g_yyyy_0_xxxxz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyy_0[i] = 4.0 * g_yyy_0_xxxyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyy_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxyz_0[i] = 4.0 * g_yyy_0_xxxyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyz_1[i] * fz_be_0 + g_yyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxzz_0[i] = 4.0 * g_yyy_0_xxxzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxzz_1[i] * fz_be_0 + g_yyyy_0_xxxzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyy_0[i] = 4.0 * g_yyy_0_xxyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyy_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxyyz_0[i] = 4.0 * g_yyy_0_xxyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxyzz_0[i] = 4.0 * g_yyy_0_xxyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyzz_1[i] * fz_be_0 + g_yyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxzzz_0[i] = 4.0 * g_yyy_0_xxzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxzzz_1[i] * fz_be_0 + g_yyyy_0_xxzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyy_0[i] = 4.0 * g_yyy_0_xyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyy_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xyyyz_0[i] = 4.0 * g_yyy_0_xyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xyyzz_0[i] = 4.0 * g_yyy_0_xyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xyzzz_0[i] = 4.0 * g_yyy_0_xyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyzzz_1[i] * fz_be_0 + g_yyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xzzzz_0[i] = 4.0 * g_yyy_0_xzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xzzzz_1[i] * fz_be_0 + g_yyyy_0_xzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyy_0[i] = 4.0 * g_yyy_0_yyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyy_1[i] * fz_be_0 + 5.0 * g_yyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyy_0_yyyyy_1[i] * wa_y[i];

        g_yyyyy_0_yyyyz_0[i] = 4.0 * g_yyy_0_yyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyy_0_yyyyz_1[i] * wa_y[i];

        g_yyyyy_0_yyyzz_0[i] = 4.0 * g_yyy_0_yyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyy_0_yyyzz_1[i] * wa_y[i];

        g_yyyyy_0_yyzzz_0[i] = 4.0 * g_yyy_0_yyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyy_0_yyzzz_1[i] * wa_y[i];

        g_yyyyy_0_yzzzz_0[i] = 4.0 * g_yyy_0_yzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yzzzz_1[i] * fz_be_0 + g_yyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyy_0_yzzzz_1[i] * wa_y[i];

        g_yyyyy_0_zzzzz_0[i] = 4.0 * g_yyy_0_zzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_zzzzz_1[i] * fz_be_0 + g_yyyy_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 336-357 components of targeted buffer : HSH

    auto g_yyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 336);

    auto g_yyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 337);

    auto g_yyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 338);

    auto g_yyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 339);

    auto g_yyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 340);

    auto g_yyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 341);

    auto g_yyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 342);

    auto g_yyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 343);

    auto g_yyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 344);

    auto g_yyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 345);

    auto g_yyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 346);

    auto g_yyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 347);

    auto g_yyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 348);

    auto g_yyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 349);

    auto g_yyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 350);

    auto g_yyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 351);

    auto g_yyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 352);

    auto g_yyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 353);

    auto g_yyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 354);

    auto g_yyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 355);

    auto g_yyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 356);

    #pragma omp simd aligned(g_yyyy_0_xxxx_1, g_yyyy_0_xxxxx_1, g_yyyy_0_xxxxy_1, g_yyyy_0_xxxxz_1, g_yyyy_0_xxxy_1, g_yyyy_0_xxxyy_1, g_yyyy_0_xxxyz_1, g_yyyy_0_xxxz_1, g_yyyy_0_xxxzz_1, g_yyyy_0_xxyy_1, g_yyyy_0_xxyyy_1, g_yyyy_0_xxyyz_1, g_yyyy_0_xxyz_1, g_yyyy_0_xxyzz_1, g_yyyy_0_xxzz_1, g_yyyy_0_xxzzz_1, g_yyyy_0_xyyy_1, g_yyyy_0_xyyyy_1, g_yyyy_0_xyyyz_1, g_yyyy_0_xyyz_1, g_yyyy_0_xyyzz_1, g_yyyy_0_xyzz_1, g_yyyy_0_xyzzz_1, g_yyyy_0_xzzz_1, g_yyyy_0_xzzzz_1, g_yyyy_0_yyyy_1, g_yyyy_0_yyyyy_1, g_yyyy_0_yyyyz_1, g_yyyy_0_yyyz_1, g_yyyy_0_yyyzz_1, g_yyyy_0_yyzz_1, g_yyyy_0_yyzzz_1, g_yyyy_0_yzzz_1, g_yyyy_0_yzzzz_1, g_yyyy_0_zzzz_1, g_yyyy_0_zzzzz_1, g_yyyyz_0_xxxxx_0, g_yyyyz_0_xxxxy_0, g_yyyyz_0_xxxxz_0, g_yyyyz_0_xxxyy_0, g_yyyyz_0_xxxyz_0, g_yyyyz_0_xxxzz_0, g_yyyyz_0_xxyyy_0, g_yyyyz_0_xxyyz_0, g_yyyyz_0_xxyzz_0, g_yyyyz_0_xxzzz_0, g_yyyyz_0_xyyyy_0, g_yyyyz_0_xyyyz_0, g_yyyyz_0_xyyzz_0, g_yyyyz_0_xyzzz_0, g_yyyyz_0_xzzzz_0, g_yyyyz_0_yyyyy_0, g_yyyyz_0_yyyyz_0, g_yyyyz_0_yyyzz_0, g_yyyyz_0_yyzzz_0, g_yyyyz_0_yzzzz_0, g_yyyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyz_0_xxxxx_0[i] = g_yyyy_0_xxxxx_1[i] * wa_z[i];

        g_yyyyz_0_xxxxy_0[i] = g_yyyy_0_xxxxy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxz_0[i] = g_yyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyy_0[i] = g_yyyy_0_xxxyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxyz_0[i] = g_yyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxzz_0[i] = 2.0 * g_yyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyy_0[i] = g_yyyy_0_xxyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxyyz_0[i] = g_yyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxyzz_0[i] = 2.0 * g_yyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyy_0_xxyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxzzz_0[i] = 3.0 * g_yyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyy_0_xxzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyy_0[i] = g_yyyy_0_xyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xyyyz_0[i] = g_yyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xyyzz_0[i] = 2.0 * g_yyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xyzzz_0[i] = 3.0 * g_yyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xzzzz_0[i] = 4.0 * g_yyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyy_0_xzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyy_0[i] = g_yyyy_0_yyyyy_1[i] * wa_z[i];

        g_yyyyz_0_yyyyz_0[i] = g_yyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyy_0_yyyyz_1[i] * wa_z[i];

        g_yyyyz_0_yyyzz_0[i] = 2.0 * g_yyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyy_0_yyyzz_1[i] * wa_z[i];

        g_yyyyz_0_yyzzz_0[i] = 3.0 * g_yyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyy_0_yyzzz_1[i] * wa_z[i];

        g_yyyyz_0_yzzzz_0[i] = 4.0 * g_yyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyy_0_yzzzz_1[i] * wa_z[i];

        g_yyyyz_0_zzzzz_0[i] = 5.0 * g_yyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 357-378 components of targeted buffer : HSH

    auto g_yyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 357);

    auto g_yyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 358);

    auto g_yyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 359);

    auto g_yyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 360);

    auto g_yyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 361);

    auto g_yyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 362);

    auto g_yyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 363);

    auto g_yyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 364);

    auto g_yyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 365);

    auto g_yyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 366);

    auto g_yyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 367);

    auto g_yyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 368);

    auto g_yyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 369);

    auto g_yyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 370);

    auto g_yyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 371);

    auto g_yyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 372);

    auto g_yyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 373);

    auto g_yyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 374);

    auto g_yyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 375);

    auto g_yyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 376);

    auto g_yyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 377);

    #pragma omp simd aligned(g_yyy_0_xxxxy_0, g_yyy_0_xxxxy_1, g_yyy_0_xxxyy_0, g_yyy_0_xxxyy_1, g_yyy_0_xxyyy_0, g_yyy_0_xxyyy_1, g_yyy_0_xyyyy_0, g_yyy_0_xyyyy_1, g_yyy_0_yyyyy_0, g_yyy_0_yyyyy_1, g_yyyz_0_xxxxy_1, g_yyyz_0_xxxyy_1, g_yyyz_0_xxyyy_1, g_yyyz_0_xyyyy_1, g_yyyz_0_yyyyy_1, g_yyyzz_0_xxxxx_0, g_yyyzz_0_xxxxy_0, g_yyyzz_0_xxxxz_0, g_yyyzz_0_xxxyy_0, g_yyyzz_0_xxxyz_0, g_yyyzz_0_xxxzz_0, g_yyyzz_0_xxyyy_0, g_yyyzz_0_xxyyz_0, g_yyyzz_0_xxyzz_0, g_yyyzz_0_xxzzz_0, g_yyyzz_0_xyyyy_0, g_yyyzz_0_xyyyz_0, g_yyyzz_0_xyyzz_0, g_yyyzz_0_xyzzz_0, g_yyyzz_0_xzzzz_0, g_yyyzz_0_yyyyy_0, g_yyyzz_0_yyyyz_0, g_yyyzz_0_yyyzz_0, g_yyyzz_0_yyzzz_0, g_yyyzz_0_yzzzz_0, g_yyyzz_0_zzzzz_0, g_yyzz_0_xxxxx_1, g_yyzz_0_xxxxz_1, g_yyzz_0_xxxyz_1, g_yyzz_0_xxxz_1, g_yyzz_0_xxxzz_1, g_yyzz_0_xxyyz_1, g_yyzz_0_xxyz_1, g_yyzz_0_xxyzz_1, g_yyzz_0_xxzz_1, g_yyzz_0_xxzzz_1, g_yyzz_0_xyyyz_1, g_yyzz_0_xyyz_1, g_yyzz_0_xyyzz_1, g_yyzz_0_xyzz_1, g_yyzz_0_xyzzz_1, g_yyzz_0_xzzz_1, g_yyzz_0_xzzzz_1, g_yyzz_0_yyyyz_1, g_yyzz_0_yyyz_1, g_yyzz_0_yyyzz_1, g_yyzz_0_yyzz_1, g_yyzz_0_yyzzz_1, g_yyzz_0_yzzz_1, g_yyzz_0_yzzzz_1, g_yyzz_0_zzzz_1, g_yyzz_0_zzzzz_1, g_yzz_0_xxxxx_0, g_yzz_0_xxxxx_1, g_yzz_0_xxxxz_0, g_yzz_0_xxxxz_1, g_yzz_0_xxxyz_0, g_yzz_0_xxxyz_1, g_yzz_0_xxxzz_0, g_yzz_0_xxxzz_1, g_yzz_0_xxyyz_0, g_yzz_0_xxyyz_1, g_yzz_0_xxyzz_0, g_yzz_0_xxyzz_1, g_yzz_0_xxzzz_0, g_yzz_0_xxzzz_1, g_yzz_0_xyyyz_0, g_yzz_0_xyyyz_1, g_yzz_0_xyyzz_0, g_yzz_0_xyyzz_1, g_yzz_0_xyzzz_0, g_yzz_0_xyzzz_1, g_yzz_0_xzzzz_0, g_yzz_0_xzzzz_1, g_yzz_0_yyyyz_0, g_yzz_0_yyyyz_1, g_yzz_0_yyyzz_0, g_yzz_0_yyyzz_1, g_yzz_0_yyzzz_0, g_yzz_0_yyzzz_1, g_yzz_0_yzzzz_0, g_yzz_0_yzzzz_1, g_yzz_0_zzzzz_0, g_yzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzz_0_xxxxx_0[i] = 2.0 * g_yzz_0_xxxxx_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxx_1[i] * fz_be_0 + g_yyzz_0_xxxxx_1[i] * wa_y[i];

        g_yyyzz_0_xxxxy_0[i] = g_yyy_0_xxxxy_0[i] * fbe_0 - g_yyy_0_xxxxy_1[i] * fz_be_0 + g_yyyz_0_xxxxy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxz_0[i] = 2.0 * g_yzz_0_xxxxz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxz_1[i] * fz_be_0 + g_yyzz_0_xxxxz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyy_0[i] = g_yyy_0_xxxyy_0[i] * fbe_0 - g_yyy_0_xxxyy_1[i] * fz_be_0 + g_yyyz_0_xxxyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxyz_0[i] = 2.0 * g_yzz_0_xxxyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyz_1[i] * fz_be_0 + g_yyzz_0_xxxz_1[i] * fi_acd_0 + g_yyzz_0_xxxyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxzz_0[i] = 2.0 * g_yzz_0_xxxzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxzz_1[i] * fz_be_0 + g_yyzz_0_xxxzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyy_0[i] = g_yyy_0_xxyyy_0[i] * fbe_0 - g_yyy_0_xxyyy_1[i] * fz_be_0 + g_yyyz_0_xxyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxyyz_0[i] = 2.0 * g_yzz_0_xxyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxyz_1[i] * fi_acd_0 + g_yyzz_0_xxyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxyzz_0[i] = 2.0 * g_yzz_0_xxyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyzz_1[i] * fz_be_0 + g_yyzz_0_xxzz_1[i] * fi_acd_0 + g_yyzz_0_xxyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxzzz_0[i] = 2.0 * g_yzz_0_xxzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxzzz_1[i] * fz_be_0 + g_yyzz_0_xxzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyy_0[i] = g_yyy_0_xyyyy_0[i] * fbe_0 - g_yyy_0_xyyyy_1[i] * fz_be_0 + g_yyyz_0_xyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xyyyz_0[i] = 2.0 * g_yzz_0_xyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xyyz_1[i] * fi_acd_0 + g_yyzz_0_xyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xyyzz_0[i] = 2.0 * g_yzz_0_xyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xyzz_1[i] * fi_acd_0 + g_yyzz_0_xyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xyzzz_0[i] = 2.0 * g_yzz_0_xyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyzzz_1[i] * fz_be_0 + g_yyzz_0_xzzz_1[i] * fi_acd_0 + g_yyzz_0_xyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xzzzz_0[i] = 2.0 * g_yzz_0_xzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xzzzz_1[i] * fz_be_0 + g_yyzz_0_xzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyy_0[i] = g_yyy_0_yyyyy_0[i] * fbe_0 - g_yyy_0_yyyyy_1[i] * fz_be_0 + g_yyyz_0_yyyyy_1[i] * wa_z[i];

        g_yyyzz_0_yyyyz_0[i] = 2.0 * g_yzz_0_yyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_yyyz_1[i] * fi_acd_0 + g_yyzz_0_yyyyz_1[i] * wa_y[i];

        g_yyyzz_0_yyyzz_0[i] = 2.0 * g_yzz_0_yyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_yyzz_1[i] * fi_acd_0 + g_yyzz_0_yyyzz_1[i] * wa_y[i];

        g_yyyzz_0_yyzzz_0[i] = 2.0 * g_yzz_0_yyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_yzzz_1[i] * fi_acd_0 + g_yyzz_0_yyzzz_1[i] * wa_y[i];

        g_yyyzz_0_yzzzz_0[i] = 2.0 * g_yzz_0_yzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yzzzz_1[i] * fz_be_0 + g_yyzz_0_zzzz_1[i] * fi_acd_0 + g_yyzz_0_yzzzz_1[i] * wa_y[i];

        g_yyyzz_0_zzzzz_0[i] = 2.0 * g_yzz_0_zzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_zzzzz_1[i] * fz_be_0 + g_yyzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 378-399 components of targeted buffer : HSH

    auto g_yyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 378);

    auto g_yyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 379);

    auto g_yyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 380);

    auto g_yyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 381);

    auto g_yyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 382);

    auto g_yyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 383);

    auto g_yyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 384);

    auto g_yyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 385);

    auto g_yyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 386);

    auto g_yyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 387);

    auto g_yyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 388);

    auto g_yyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 389);

    auto g_yyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 390);

    auto g_yyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 391);

    auto g_yyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 392);

    auto g_yyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 393);

    auto g_yyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 394);

    auto g_yyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 395);

    auto g_yyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 396);

    auto g_yyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 397);

    auto g_yyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 398);

    #pragma omp simd aligned(g_yyz_0_xxxxy_0, g_yyz_0_xxxxy_1, g_yyz_0_xxxyy_0, g_yyz_0_xxxyy_1, g_yyz_0_xxyyy_0, g_yyz_0_xxyyy_1, g_yyz_0_xyyyy_0, g_yyz_0_xyyyy_1, g_yyz_0_yyyyy_0, g_yyz_0_yyyyy_1, g_yyzz_0_xxxxy_1, g_yyzz_0_xxxyy_1, g_yyzz_0_xxyyy_1, g_yyzz_0_xyyyy_1, g_yyzz_0_yyyyy_1, g_yyzzz_0_xxxxx_0, g_yyzzz_0_xxxxy_0, g_yyzzz_0_xxxxz_0, g_yyzzz_0_xxxyy_0, g_yyzzz_0_xxxyz_0, g_yyzzz_0_xxxzz_0, g_yyzzz_0_xxyyy_0, g_yyzzz_0_xxyyz_0, g_yyzzz_0_xxyzz_0, g_yyzzz_0_xxzzz_0, g_yyzzz_0_xyyyy_0, g_yyzzz_0_xyyyz_0, g_yyzzz_0_xyyzz_0, g_yyzzz_0_xyzzz_0, g_yyzzz_0_xzzzz_0, g_yyzzz_0_yyyyy_0, g_yyzzz_0_yyyyz_0, g_yyzzz_0_yyyzz_0, g_yyzzz_0_yyzzz_0, g_yyzzz_0_yzzzz_0, g_yyzzz_0_zzzzz_0, g_yzzz_0_xxxxx_1, g_yzzz_0_xxxxz_1, g_yzzz_0_xxxyz_1, g_yzzz_0_xxxz_1, g_yzzz_0_xxxzz_1, g_yzzz_0_xxyyz_1, g_yzzz_0_xxyz_1, g_yzzz_0_xxyzz_1, g_yzzz_0_xxzz_1, g_yzzz_0_xxzzz_1, g_yzzz_0_xyyyz_1, g_yzzz_0_xyyz_1, g_yzzz_0_xyyzz_1, g_yzzz_0_xyzz_1, g_yzzz_0_xyzzz_1, g_yzzz_0_xzzz_1, g_yzzz_0_xzzzz_1, g_yzzz_0_yyyyz_1, g_yzzz_0_yyyz_1, g_yzzz_0_yyyzz_1, g_yzzz_0_yyzz_1, g_yzzz_0_yyzzz_1, g_yzzz_0_yzzz_1, g_yzzz_0_yzzzz_1, g_yzzz_0_zzzz_1, g_yzzz_0_zzzzz_1, g_zzz_0_xxxxx_0, g_zzz_0_xxxxx_1, g_zzz_0_xxxxz_0, g_zzz_0_xxxxz_1, g_zzz_0_xxxyz_0, g_zzz_0_xxxyz_1, g_zzz_0_xxxzz_0, g_zzz_0_xxxzz_1, g_zzz_0_xxyyz_0, g_zzz_0_xxyyz_1, g_zzz_0_xxyzz_0, g_zzz_0_xxyzz_1, g_zzz_0_xxzzz_0, g_zzz_0_xxzzz_1, g_zzz_0_xyyyz_0, g_zzz_0_xyyyz_1, g_zzz_0_xyyzz_0, g_zzz_0_xyyzz_1, g_zzz_0_xyzzz_0, g_zzz_0_xyzzz_1, g_zzz_0_xzzzz_0, g_zzz_0_xzzzz_1, g_zzz_0_yyyyz_0, g_zzz_0_yyyyz_1, g_zzz_0_yyyzz_0, g_zzz_0_yyyzz_1, g_zzz_0_yyzzz_0, g_zzz_0_yyzzz_1, g_zzz_0_yzzzz_0, g_zzz_0_yzzzz_1, g_zzz_0_zzzzz_0, g_zzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzz_0_xxxxx_0[i] = g_zzz_0_xxxxx_0[i] * fbe_0 - g_zzz_0_xxxxx_1[i] * fz_be_0 + g_yzzz_0_xxxxx_1[i] * wa_y[i];

        g_yyzzz_0_xxxxy_0[i] = 2.0 * g_yyz_0_xxxxy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxy_1[i] * fz_be_0 + g_yyzz_0_xxxxy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxz_0[i] = g_zzz_0_xxxxz_0[i] * fbe_0 - g_zzz_0_xxxxz_1[i] * fz_be_0 + g_yzzz_0_xxxxz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyy_0[i] = 2.0 * g_yyz_0_xxxyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxyy_1[i] * fz_be_0 + g_yyzz_0_xxxyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxyz_0[i] = g_zzz_0_xxxyz_0[i] * fbe_0 - g_zzz_0_xxxyz_1[i] * fz_be_0 + g_yzzz_0_xxxz_1[i] * fi_acd_0 + g_yzzz_0_xxxyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxzz_0[i] = g_zzz_0_xxxzz_0[i] * fbe_0 - g_zzz_0_xxxzz_1[i] * fz_be_0 + g_yzzz_0_xxxzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyy_0[i] = 2.0 * g_yyz_0_xxyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxyyy_1[i] * fz_be_0 + g_yyzz_0_xxyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxyyz_0[i] = g_zzz_0_xxyyz_0[i] * fbe_0 - g_zzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxyz_1[i] * fi_acd_0 + g_yzzz_0_xxyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxyzz_0[i] = g_zzz_0_xxyzz_0[i] * fbe_0 - g_zzz_0_xxyzz_1[i] * fz_be_0 + g_yzzz_0_xxzz_1[i] * fi_acd_0 + g_yzzz_0_xxyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxzzz_0[i] = g_zzz_0_xxzzz_0[i] * fbe_0 - g_zzz_0_xxzzz_1[i] * fz_be_0 + g_yzzz_0_xxzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyy_0[i] = 2.0 * g_yyz_0_xyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xyyyy_1[i] * fz_be_0 + g_yyzz_0_xyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xyyyz_0[i] = g_zzz_0_xyyyz_0[i] * fbe_0 - g_zzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xyyz_1[i] * fi_acd_0 + g_yzzz_0_xyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xyyzz_0[i] = g_zzz_0_xyyzz_0[i] * fbe_0 - g_zzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xyzz_1[i] * fi_acd_0 + g_yzzz_0_xyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xyzzz_0[i] = g_zzz_0_xyzzz_0[i] * fbe_0 - g_zzz_0_xyzzz_1[i] * fz_be_0 + g_yzzz_0_xzzz_1[i] * fi_acd_0 + g_yzzz_0_xyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xzzzz_0[i] = g_zzz_0_xzzzz_0[i] * fbe_0 - g_zzz_0_xzzzz_1[i] * fz_be_0 + g_yzzz_0_xzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyy_0[i] = 2.0 * g_yyz_0_yyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_yyyyy_1[i] * fz_be_0 + g_yyzz_0_yyyyy_1[i] * wa_z[i];

        g_yyzzz_0_yyyyz_0[i] = g_zzz_0_yyyyz_0[i] * fbe_0 - g_zzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_yyyz_1[i] * fi_acd_0 + g_yzzz_0_yyyyz_1[i] * wa_y[i];

        g_yyzzz_0_yyyzz_0[i] = g_zzz_0_yyyzz_0[i] * fbe_0 - g_zzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_yyzz_1[i] * fi_acd_0 + g_yzzz_0_yyyzz_1[i] * wa_y[i];

        g_yyzzz_0_yyzzz_0[i] = g_zzz_0_yyzzz_0[i] * fbe_0 - g_zzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_yzzz_1[i] * fi_acd_0 + g_yzzz_0_yyzzz_1[i] * wa_y[i];

        g_yyzzz_0_yzzzz_0[i] = g_zzz_0_yzzzz_0[i] * fbe_0 - g_zzz_0_yzzzz_1[i] * fz_be_0 + g_yzzz_0_zzzz_1[i] * fi_acd_0 + g_yzzz_0_yzzzz_1[i] * wa_y[i];

        g_yyzzz_0_zzzzz_0[i] = g_zzz_0_zzzzz_0[i] * fbe_0 - g_zzz_0_zzzzz_1[i] * fz_be_0 + g_yzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 399-420 components of targeted buffer : HSH

    auto g_yzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 399);

    auto g_yzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 400);

    auto g_yzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 401);

    auto g_yzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 402);

    auto g_yzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 403);

    auto g_yzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 404);

    auto g_yzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 405);

    auto g_yzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 406);

    auto g_yzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 407);

    auto g_yzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 408);

    auto g_yzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 409);

    auto g_yzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 410);

    auto g_yzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 411);

    auto g_yzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 412);

    auto g_yzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 413);

    auto g_yzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 414);

    auto g_yzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 415);

    auto g_yzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 416);

    auto g_yzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 417);

    auto g_yzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 418);

    auto g_yzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 419);

    #pragma omp simd aligned(g_yzzzz_0_xxxxx_0, g_yzzzz_0_xxxxy_0, g_yzzzz_0_xxxxz_0, g_yzzzz_0_xxxyy_0, g_yzzzz_0_xxxyz_0, g_yzzzz_0_xxxzz_0, g_yzzzz_0_xxyyy_0, g_yzzzz_0_xxyyz_0, g_yzzzz_0_xxyzz_0, g_yzzzz_0_xxzzz_0, g_yzzzz_0_xyyyy_0, g_yzzzz_0_xyyyz_0, g_yzzzz_0_xyyzz_0, g_yzzzz_0_xyzzz_0, g_yzzzz_0_xzzzz_0, g_yzzzz_0_yyyyy_0, g_yzzzz_0_yyyyz_0, g_yzzzz_0_yyyzz_0, g_yzzzz_0_yyzzz_0, g_yzzzz_0_yzzzz_0, g_yzzzz_0_zzzzz_0, g_zzzz_0_xxxx_1, g_zzzz_0_xxxxx_1, g_zzzz_0_xxxxy_1, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxy_1, g_zzzz_0_xxxyy_1, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxz_1, g_zzzz_0_xxxzz_1, g_zzzz_0_xxyy_1, g_zzzz_0_xxyyy_1, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyz_1, g_zzzz_0_xxyzz_1, g_zzzz_0_xxzz_1, g_zzzz_0_xxzzz_1, g_zzzz_0_xyyy_1, g_zzzz_0_xyyyy_1, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyz_1, g_zzzz_0_xyyzz_1, g_zzzz_0_xyzz_1, g_zzzz_0_xyzzz_1, g_zzzz_0_xzzz_1, g_zzzz_0_xzzzz_1, g_zzzz_0_yyyy_1, g_zzzz_0_yyyyy_1, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyz_1, g_zzzz_0_yyyzz_1, g_zzzz_0_yyzz_1, g_zzzz_0_yyzzz_1, g_zzzz_0_yzzz_1, g_zzzz_0_yzzzz_1, g_zzzz_0_zzzz_1, g_zzzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzz_0_xxxxx_0[i] = g_zzzz_0_xxxxx_1[i] * wa_y[i];

        g_yzzzz_0_xxxxy_0[i] = g_zzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxz_0[i] = g_zzzz_0_xxxxz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyy_0[i] = 2.0 * g_zzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxyz_0[i] = g_zzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxzz_0[i] = g_zzzz_0_xxxzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyy_0[i] = 3.0 * g_zzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxyyz_0[i] = 2.0 * g_zzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxyzz_0[i] = g_zzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxzzz_0[i] = g_zzzz_0_xxzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyy_0[i] = 4.0 * g_zzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xyyyz_0[i] = 3.0 * g_zzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xyyzz_0[i] = 2.0 * g_zzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xyzzz_0[i] = g_zzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xzzzz_0[i] = g_zzzz_0_xzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyy_0[i] = 5.0 * g_zzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzz_0_yyyyy_1[i] * wa_y[i];

        g_yzzzz_0_yyyyz_0[i] = 4.0 * g_zzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzz_0_yyyyz_1[i] * wa_y[i];

        g_yzzzz_0_yyyzz_0[i] = 3.0 * g_zzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzz_0_yyyzz_1[i] * wa_y[i];

        g_yzzzz_0_yyzzz_0[i] = 2.0 * g_zzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzz_0_yyzzz_1[i] * wa_y[i];

        g_yzzzz_0_yzzzz_0[i] = g_zzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzz_0_yzzzz_1[i] * wa_y[i];

        g_yzzzz_0_zzzzz_0[i] = g_zzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 420-441 components of targeted buffer : HSH

    auto g_zzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 420);

    auto g_zzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 421);

    auto g_zzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 422);

    auto g_zzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 423);

    auto g_zzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 424);

    auto g_zzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 425);

    auto g_zzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 426);

    auto g_zzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 427);

    auto g_zzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 428);

    auto g_zzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 429);

    auto g_zzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 430);

    auto g_zzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 431);

    auto g_zzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 432);

    auto g_zzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 433);

    auto g_zzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 434);

    auto g_zzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 435);

    auto g_zzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 436);

    auto g_zzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 437);

    auto g_zzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 438);

    auto g_zzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 439);

    auto g_zzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 440);

    #pragma omp simd aligned(g_zzz_0_xxxxx_0, g_zzz_0_xxxxx_1, g_zzz_0_xxxxy_0, g_zzz_0_xxxxy_1, g_zzz_0_xxxxz_0, g_zzz_0_xxxxz_1, g_zzz_0_xxxyy_0, g_zzz_0_xxxyy_1, g_zzz_0_xxxyz_0, g_zzz_0_xxxyz_1, g_zzz_0_xxxzz_0, g_zzz_0_xxxzz_1, g_zzz_0_xxyyy_0, g_zzz_0_xxyyy_1, g_zzz_0_xxyyz_0, g_zzz_0_xxyyz_1, g_zzz_0_xxyzz_0, g_zzz_0_xxyzz_1, g_zzz_0_xxzzz_0, g_zzz_0_xxzzz_1, g_zzz_0_xyyyy_0, g_zzz_0_xyyyy_1, g_zzz_0_xyyyz_0, g_zzz_0_xyyyz_1, g_zzz_0_xyyzz_0, g_zzz_0_xyyzz_1, g_zzz_0_xyzzz_0, g_zzz_0_xyzzz_1, g_zzz_0_xzzzz_0, g_zzz_0_xzzzz_1, g_zzz_0_yyyyy_0, g_zzz_0_yyyyy_1, g_zzz_0_yyyyz_0, g_zzz_0_yyyyz_1, g_zzz_0_yyyzz_0, g_zzz_0_yyyzz_1, g_zzz_0_yyzzz_0, g_zzz_0_yyzzz_1, g_zzz_0_yzzzz_0, g_zzz_0_yzzzz_1, g_zzz_0_zzzzz_0, g_zzz_0_zzzzz_1, g_zzzz_0_xxxx_1, g_zzzz_0_xxxxx_1, g_zzzz_0_xxxxy_1, g_zzzz_0_xxxxz_1, g_zzzz_0_xxxy_1, g_zzzz_0_xxxyy_1, g_zzzz_0_xxxyz_1, g_zzzz_0_xxxz_1, g_zzzz_0_xxxzz_1, g_zzzz_0_xxyy_1, g_zzzz_0_xxyyy_1, g_zzzz_0_xxyyz_1, g_zzzz_0_xxyz_1, g_zzzz_0_xxyzz_1, g_zzzz_0_xxzz_1, g_zzzz_0_xxzzz_1, g_zzzz_0_xyyy_1, g_zzzz_0_xyyyy_1, g_zzzz_0_xyyyz_1, g_zzzz_0_xyyz_1, g_zzzz_0_xyyzz_1, g_zzzz_0_xyzz_1, g_zzzz_0_xyzzz_1, g_zzzz_0_xzzz_1, g_zzzz_0_xzzzz_1, g_zzzz_0_yyyy_1, g_zzzz_0_yyyyy_1, g_zzzz_0_yyyyz_1, g_zzzz_0_yyyz_1, g_zzzz_0_yyyzz_1, g_zzzz_0_yyzz_1, g_zzzz_0_yyzzz_1, g_zzzz_0_yzzz_1, g_zzzz_0_yzzzz_1, g_zzzz_0_zzzz_1, g_zzzz_0_zzzzz_1, g_zzzzz_0_xxxxx_0, g_zzzzz_0_xxxxy_0, g_zzzzz_0_xxxxz_0, g_zzzzz_0_xxxyy_0, g_zzzzz_0_xxxyz_0, g_zzzzz_0_xxxzz_0, g_zzzzz_0_xxyyy_0, g_zzzzz_0_xxyyz_0, g_zzzzz_0_xxyzz_0, g_zzzzz_0_xxzzz_0, g_zzzzz_0_xyyyy_0, g_zzzzz_0_xyyyz_0, g_zzzzz_0_xyyzz_0, g_zzzzz_0_xyzzz_0, g_zzzzz_0_xzzzz_0, g_zzzzz_0_yyyyy_0, g_zzzzz_0_yyyyz_0, g_zzzzz_0_yyyzz_0, g_zzzzz_0_yyzzz_0, g_zzzzz_0_yzzzz_0, g_zzzzz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzz_0_xxxxx_0[i] = 4.0 * g_zzz_0_xxxxx_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxx_1[i] * fz_be_0 + g_zzzz_0_xxxxx_1[i] * wa_z[i];

        g_zzzzz_0_xxxxy_0[i] = 4.0 * g_zzz_0_xxxxy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxy_1[i] * fz_be_0 + g_zzzz_0_xxxxy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxz_0[i] = 4.0 * g_zzz_0_xxxxz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxz_1[i] * fz_be_0 + g_zzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyy_0[i] = 4.0 * g_zzz_0_xxxyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyy_1[i] * fz_be_0 + g_zzzz_0_xxxyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxyz_0[i] = 4.0 * g_zzz_0_xxxyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyz_1[i] * fz_be_0 + g_zzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxzz_0[i] = 4.0 * g_zzz_0_xxxzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyy_0[i] = 4.0 * g_zzz_0_xxyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyy_1[i] * fz_be_0 + g_zzzz_0_xxyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxyyz_0[i] = 4.0 * g_zzz_0_xxyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyz_1[i] * fz_be_0 + g_zzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxyzz_0[i] = 4.0 * g_zzz_0_xxyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzz_0_xxyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxzzz_0[i] = 4.0 * g_zzz_0_xxzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzz_0_xxzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyy_0[i] = 4.0 * g_zzz_0_xyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyy_1[i] * fz_be_0 + g_zzzz_0_xyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xyyyz_0[i] = 4.0 * g_zzz_0_xyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyz_1[i] * fz_be_0 + g_zzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xyyzz_0[i] = 4.0 * g_zzz_0_xyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xyzzz_0[i] = 4.0 * g_zzz_0_xyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xzzzz_0[i] = 4.0 * g_zzz_0_xzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzz_0_xzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyy_0[i] = 4.0 * g_zzz_0_yyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyy_1[i] * fz_be_0 + g_zzzz_0_yyyyy_1[i] * wa_z[i];

        g_zzzzz_0_yyyyz_0[i] = 4.0 * g_zzz_0_yyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyz_1[i] * fz_be_0 + g_zzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzz_0_yyyyz_1[i] * wa_z[i];

        g_zzzzz_0_yyyzz_0[i] = 4.0 * g_zzz_0_yyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzz_0_yyyzz_1[i] * wa_z[i];

        g_zzzzz_0_yyzzz_0[i] = 4.0 * g_zzz_0_yyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzz_0_yyzzz_1[i] * wa_z[i];

        g_zzzzz_0_yzzzz_0[i] = 4.0 * g_zzz_0_yzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzz_0_yzzzz_1[i] * wa_z[i];

        g_zzzzz_0_zzzzz_0[i] = 4.0 * g_zzz_0_zzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_zzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzz_0_zzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

