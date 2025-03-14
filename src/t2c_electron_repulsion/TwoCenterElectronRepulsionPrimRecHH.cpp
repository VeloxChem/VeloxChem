#include "TwoCenterElectronRepulsionPrimRecHH.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_hh(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_hh,
                                const size_t idx_eri_0_fh,
                                const size_t idx_eri_1_fh,
                                const size_t idx_eri_1_gg,
                                const size_t idx_eri_1_gh,
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

    // Set up components of auxiliary buffer : FH

    auto g_xxx_xxxxx_0 = pbuffer.data(idx_eri_0_fh);

    auto g_xxx_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 1);

    auto g_xxx_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 2);

    auto g_xxx_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 3);

    auto g_xxx_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 4);

    auto g_xxx_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 5);

    auto g_xxx_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 6);

    auto g_xxx_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 7);

    auto g_xxx_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 8);

    auto g_xxx_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 9);

    auto g_xxx_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 10);

    auto g_xxx_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 11);

    auto g_xxx_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 12);

    auto g_xxx_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 13);

    auto g_xxx_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 14);

    auto g_xxx_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 15);

    auto g_xxx_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 16);

    auto g_xxx_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 17);

    auto g_xxx_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 18);

    auto g_xxx_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 19);

    auto g_xxx_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 20);

    auto g_xxy_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 21);

    auto g_xxy_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 23);

    auto g_xxy_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 26);

    auto g_xxy_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 30);

    auto g_xxy_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 35);

    auto g_xxz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 42);

    auto g_xxz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 43);

    auto g_xxz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 45);

    auto g_xxz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 48);

    auto g_xxz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 52);

    auto g_xyy_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 64);

    auto g_xyy_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 66);

    auto g_xyy_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 67);

    auto g_xyy_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 69);

    auto g_xyy_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 70);

    auto g_xyy_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 71);

    auto g_xyy_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 73);

    auto g_xyy_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 74);

    auto g_xyy_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 75);

    auto g_xyy_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 76);

    auto g_xyy_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 78);

    auto g_xyy_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 79);

    auto g_xyy_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 80);

    auto g_xyy_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 81);

    auto g_xyy_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 82);

    auto g_xyy_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 83);

    auto g_xzz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 107);

    auto g_xzz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 109);

    auto g_xzz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 110);

    auto g_xzz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 112);

    auto g_xzz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 113);

    auto g_xzz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 114);

    auto g_xzz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 116);

    auto g_xzz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 117);

    auto g_xzz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 118);

    auto g_xzz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 119);

    auto g_xzz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 120);

    auto g_xzz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 121);

    auto g_xzz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 122);

    auto g_xzz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 123);

    auto g_xzz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 124);

    auto g_xzz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 125);

    auto g_yyy_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 126);

    auto g_yyy_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 127);

    auto g_yyy_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 128);

    auto g_yyy_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 129);

    auto g_yyy_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 130);

    auto g_yyy_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 131);

    auto g_yyy_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 132);

    auto g_yyy_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 133);

    auto g_yyy_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 134);

    auto g_yyy_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 135);

    auto g_yyy_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 136);

    auto g_yyy_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 137);

    auto g_yyy_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 138);

    auto g_yyy_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 139);

    auto g_yyy_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 140);

    auto g_yyy_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 141);

    auto g_yyy_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 142);

    auto g_yyy_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 143);

    auto g_yyy_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 144);

    auto g_yyy_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 145);

    auto g_yyy_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 146);

    auto g_yyz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 148);

    auto g_yyz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 150);

    auto g_yyz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 153);

    auto g_yyz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 157);

    auto g_yyz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 162);

    auto g_yzz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 168);

    auto g_yzz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 170);

    auto g_yzz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 172);

    auto g_yzz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 173);

    auto g_yzz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 175);

    auto g_yzz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 176);

    auto g_yzz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 177);

    auto g_yzz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 179);

    auto g_yzz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 180);

    auto g_yzz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 181);

    auto g_yzz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 182);

    auto g_yzz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 184);

    auto g_yzz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 185);

    auto g_yzz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 186);

    auto g_yzz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 187);

    auto g_yzz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 188);

    auto g_zzz_xxxxx_0 = pbuffer.data(idx_eri_0_fh + 189);

    auto g_zzz_xxxxy_0 = pbuffer.data(idx_eri_0_fh + 190);

    auto g_zzz_xxxxz_0 = pbuffer.data(idx_eri_0_fh + 191);

    auto g_zzz_xxxyy_0 = pbuffer.data(idx_eri_0_fh + 192);

    auto g_zzz_xxxyz_0 = pbuffer.data(idx_eri_0_fh + 193);

    auto g_zzz_xxxzz_0 = pbuffer.data(idx_eri_0_fh + 194);

    auto g_zzz_xxyyy_0 = pbuffer.data(idx_eri_0_fh + 195);

    auto g_zzz_xxyyz_0 = pbuffer.data(idx_eri_0_fh + 196);

    auto g_zzz_xxyzz_0 = pbuffer.data(idx_eri_0_fh + 197);

    auto g_zzz_xxzzz_0 = pbuffer.data(idx_eri_0_fh + 198);

    auto g_zzz_xyyyy_0 = pbuffer.data(idx_eri_0_fh + 199);

    auto g_zzz_xyyyz_0 = pbuffer.data(idx_eri_0_fh + 200);

    auto g_zzz_xyyzz_0 = pbuffer.data(idx_eri_0_fh + 201);

    auto g_zzz_xyzzz_0 = pbuffer.data(idx_eri_0_fh + 202);

    auto g_zzz_xzzzz_0 = pbuffer.data(idx_eri_0_fh + 203);

    auto g_zzz_yyyyy_0 = pbuffer.data(idx_eri_0_fh + 204);

    auto g_zzz_yyyyz_0 = pbuffer.data(idx_eri_0_fh + 205);

    auto g_zzz_yyyzz_0 = pbuffer.data(idx_eri_0_fh + 206);

    auto g_zzz_yyzzz_0 = pbuffer.data(idx_eri_0_fh + 207);

    auto g_zzz_yzzzz_0 = pbuffer.data(idx_eri_0_fh + 208);

    auto g_zzz_zzzzz_0 = pbuffer.data(idx_eri_0_fh + 209);

    // Set up components of auxiliary buffer : FH

    auto g_xxx_xxxxx_1 = pbuffer.data(idx_eri_1_fh);

    auto g_xxx_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 1);

    auto g_xxx_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 2);

    auto g_xxx_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 3);

    auto g_xxx_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 4);

    auto g_xxx_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 5);

    auto g_xxx_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 6);

    auto g_xxx_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 7);

    auto g_xxx_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 8);

    auto g_xxx_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 9);

    auto g_xxx_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 10);

    auto g_xxx_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 11);

    auto g_xxx_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 12);

    auto g_xxx_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 13);

    auto g_xxx_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 14);

    auto g_xxx_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 15);

    auto g_xxx_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 16);

    auto g_xxx_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 17);

    auto g_xxx_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 18);

    auto g_xxx_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 19);

    auto g_xxx_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 20);

    auto g_xxy_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 21);

    auto g_xxy_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 23);

    auto g_xxy_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 26);

    auto g_xxy_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 30);

    auto g_xxy_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 35);

    auto g_xxz_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 42);

    auto g_xxz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 43);

    auto g_xxz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 45);

    auto g_xxz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 48);

    auto g_xxz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 52);

    auto g_xyy_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 64);

    auto g_xyy_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 66);

    auto g_xyy_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 67);

    auto g_xyy_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 69);

    auto g_xyy_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 70);

    auto g_xyy_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 71);

    auto g_xyy_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 73);

    auto g_xyy_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 74);

    auto g_xyy_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 75);

    auto g_xyy_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 76);

    auto g_xyy_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 78);

    auto g_xyy_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 79);

    auto g_xyy_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 80);

    auto g_xyy_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 81);

    auto g_xyy_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 82);

    auto g_xyy_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 83);

    auto g_xzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 107);

    auto g_xzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 109);

    auto g_xzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 110);

    auto g_xzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 112);

    auto g_xzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 113);

    auto g_xzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 114);

    auto g_xzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 116);

    auto g_xzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 117);

    auto g_xzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 118);

    auto g_xzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 119);

    auto g_xzz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 120);

    auto g_xzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 121);

    auto g_xzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 122);

    auto g_xzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 123);

    auto g_xzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 124);

    auto g_xzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 125);

    auto g_yyy_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 126);

    auto g_yyy_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 127);

    auto g_yyy_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 128);

    auto g_yyy_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 129);

    auto g_yyy_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 130);

    auto g_yyy_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 131);

    auto g_yyy_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 132);

    auto g_yyy_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 133);

    auto g_yyy_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 134);

    auto g_yyy_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 135);

    auto g_yyy_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 136);

    auto g_yyy_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 137);

    auto g_yyy_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 138);

    auto g_yyy_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 139);

    auto g_yyy_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 140);

    auto g_yyy_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 141);

    auto g_yyy_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 142);

    auto g_yyy_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 143);

    auto g_yyy_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 144);

    auto g_yyy_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 145);

    auto g_yyy_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 146);

    auto g_yyz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 148);

    auto g_yyz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 150);

    auto g_yyz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 153);

    auto g_yyz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 157);

    auto g_yyz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 162);

    auto g_yzz_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 168);

    auto g_yzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 170);

    auto g_yzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 172);

    auto g_yzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 173);

    auto g_yzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 175);

    auto g_yzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 176);

    auto g_yzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 177);

    auto g_yzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 179);

    auto g_yzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 180);

    auto g_yzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 181);

    auto g_yzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 182);

    auto g_yzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 184);

    auto g_yzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 185);

    auto g_yzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 186);

    auto g_yzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 187);

    auto g_yzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 188);

    auto g_zzz_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 189);

    auto g_zzz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 190);

    auto g_zzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 191);

    auto g_zzz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 192);

    auto g_zzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 193);

    auto g_zzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 194);

    auto g_zzz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 195);

    auto g_zzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 196);

    auto g_zzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 197);

    auto g_zzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 198);

    auto g_zzz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 199);

    auto g_zzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 200);

    auto g_zzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 201);

    auto g_zzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 202);

    auto g_zzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 203);

    auto g_zzz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 204);

    auto g_zzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 205);

    auto g_zzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 206);

    auto g_zzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 207);

    auto g_zzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 208);

    auto g_zzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 209);

    // Set up components of auxiliary buffer : GG

    auto g_xxxx_xxxx_1 = pbuffer.data(idx_eri_1_gg);

    auto g_xxxx_xxxy_1 = pbuffer.data(idx_eri_1_gg + 1);

    auto g_xxxx_xxxz_1 = pbuffer.data(idx_eri_1_gg + 2);

    auto g_xxxx_xxyy_1 = pbuffer.data(idx_eri_1_gg + 3);

    auto g_xxxx_xxyz_1 = pbuffer.data(idx_eri_1_gg + 4);

    auto g_xxxx_xxzz_1 = pbuffer.data(idx_eri_1_gg + 5);

    auto g_xxxx_xyyy_1 = pbuffer.data(idx_eri_1_gg + 6);

    auto g_xxxx_xyyz_1 = pbuffer.data(idx_eri_1_gg + 7);

    auto g_xxxx_xyzz_1 = pbuffer.data(idx_eri_1_gg + 8);

    auto g_xxxx_xzzz_1 = pbuffer.data(idx_eri_1_gg + 9);

    auto g_xxxx_yyyy_1 = pbuffer.data(idx_eri_1_gg + 10);

    auto g_xxxx_yyyz_1 = pbuffer.data(idx_eri_1_gg + 11);

    auto g_xxxx_yyzz_1 = pbuffer.data(idx_eri_1_gg + 12);

    auto g_xxxx_yzzz_1 = pbuffer.data(idx_eri_1_gg + 13);

    auto g_xxxx_zzzz_1 = pbuffer.data(idx_eri_1_gg + 14);

    auto g_xxxz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 32);

    auto g_xxxz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 34);

    auto g_xxxz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 35);

    auto g_xxxz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 37);

    auto g_xxxz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 38);

    auto g_xxxz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 39);

    auto g_xxxz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 41);

    auto g_xxxz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 42);

    auto g_xxxz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 43);

    auto g_xxxz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 44);

    auto g_xxyy_xxxx_1 = pbuffer.data(idx_eri_1_gg + 45);

    auto g_xxyy_xxxy_1 = pbuffer.data(idx_eri_1_gg + 46);

    auto g_xxyy_xxxz_1 = pbuffer.data(idx_eri_1_gg + 47);

    auto g_xxyy_xxyy_1 = pbuffer.data(idx_eri_1_gg + 48);

    auto g_xxyy_xxyz_1 = pbuffer.data(idx_eri_1_gg + 49);

    auto g_xxyy_xxzz_1 = pbuffer.data(idx_eri_1_gg + 50);

    auto g_xxyy_xyyy_1 = pbuffer.data(idx_eri_1_gg + 51);

    auto g_xxyy_xyyz_1 = pbuffer.data(idx_eri_1_gg + 52);

    auto g_xxyy_xyzz_1 = pbuffer.data(idx_eri_1_gg + 53);

    auto g_xxyy_xzzz_1 = pbuffer.data(idx_eri_1_gg + 54);

    auto g_xxyy_yyyy_1 = pbuffer.data(idx_eri_1_gg + 55);

    auto g_xxyy_yyyz_1 = pbuffer.data(idx_eri_1_gg + 56);

    auto g_xxyy_yyzz_1 = pbuffer.data(idx_eri_1_gg + 57);

    auto g_xxyy_yzzz_1 = pbuffer.data(idx_eri_1_gg + 58);

    auto g_xxyy_zzzz_1 = pbuffer.data(idx_eri_1_gg + 59);

    auto g_xxzz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 75);

    auto g_xxzz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 76);

    auto g_xxzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 77);

    auto g_xxzz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 78);

    auto g_xxzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 79);

    auto g_xxzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 80);

    auto g_xxzz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 81);

    auto g_xxzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 82);

    auto g_xxzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 83);

    auto g_xxzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 84);

    auto g_xxzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 85);

    auto g_xxzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 86);

    auto g_xxzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 87);

    auto g_xxzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 88);

    auto g_xxzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 89);

    auto g_xyyy_xxxy_1 = pbuffer.data(idx_eri_1_gg + 91);

    auto g_xyyy_xxyy_1 = pbuffer.data(idx_eri_1_gg + 93);

    auto g_xyyy_xxyz_1 = pbuffer.data(idx_eri_1_gg + 94);

    auto g_xyyy_xyyy_1 = pbuffer.data(idx_eri_1_gg + 96);

    auto g_xyyy_xyyz_1 = pbuffer.data(idx_eri_1_gg + 97);

    auto g_xyyy_xyzz_1 = pbuffer.data(idx_eri_1_gg + 98);

    auto g_xyyy_yyyy_1 = pbuffer.data(idx_eri_1_gg + 100);

    auto g_xyyy_yyyz_1 = pbuffer.data(idx_eri_1_gg + 101);

    auto g_xyyy_yyzz_1 = pbuffer.data(idx_eri_1_gg + 102);

    auto g_xyyy_yzzz_1 = pbuffer.data(idx_eri_1_gg + 103);

    auto g_xzzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 137);

    auto g_xzzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 139);

    auto g_xzzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 140);

    auto g_xzzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 142);

    auto g_xzzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 143);

    auto g_xzzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 144);

    auto g_xzzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 146);

    auto g_xzzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 147);

    auto g_xzzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 148);

    auto g_xzzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 149);

    auto g_yyyy_xxxx_1 = pbuffer.data(idx_eri_1_gg + 150);

    auto g_yyyy_xxxy_1 = pbuffer.data(idx_eri_1_gg + 151);

    auto g_yyyy_xxxz_1 = pbuffer.data(idx_eri_1_gg + 152);

    auto g_yyyy_xxyy_1 = pbuffer.data(idx_eri_1_gg + 153);

    auto g_yyyy_xxyz_1 = pbuffer.data(idx_eri_1_gg + 154);

    auto g_yyyy_xxzz_1 = pbuffer.data(idx_eri_1_gg + 155);

    auto g_yyyy_xyyy_1 = pbuffer.data(idx_eri_1_gg + 156);

    auto g_yyyy_xyyz_1 = pbuffer.data(idx_eri_1_gg + 157);

    auto g_yyyy_xyzz_1 = pbuffer.data(idx_eri_1_gg + 158);

    auto g_yyyy_xzzz_1 = pbuffer.data(idx_eri_1_gg + 159);

    auto g_yyyy_yyyy_1 = pbuffer.data(idx_eri_1_gg + 160);

    auto g_yyyy_yyyz_1 = pbuffer.data(idx_eri_1_gg + 161);

    auto g_yyyy_yyzz_1 = pbuffer.data(idx_eri_1_gg + 162);

    auto g_yyyy_yzzz_1 = pbuffer.data(idx_eri_1_gg + 163);

    auto g_yyyy_zzzz_1 = pbuffer.data(idx_eri_1_gg + 164);

    auto g_yyyz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 167);

    auto g_yyyz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 169);

    auto g_yyyz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 170);

    auto g_yyyz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 172);

    auto g_yyyz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 173);

    auto g_yyyz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 174);

    auto g_yyyz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 176);

    auto g_yyyz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 177);

    auto g_yyyz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 178);

    auto g_yyyz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 179);

    auto g_yyzz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 180);

    auto g_yyzz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 181);

    auto g_yyzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 182);

    auto g_yyzz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 183);

    auto g_yyzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 184);

    auto g_yyzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 185);

    auto g_yyzz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 186);

    auto g_yyzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 187);

    auto g_yyzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 188);

    auto g_yyzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 189);

    auto g_yyzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 190);

    auto g_yyzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 191);

    auto g_yyzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 192);

    auto g_yyzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 193);

    auto g_yyzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 194);

    auto g_yzzz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 196);

    auto g_yzzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 197);

    auto g_yzzz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 198);

    auto g_yzzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 199);

    auto g_yzzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 200);

    auto g_yzzz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 201);

    auto g_yzzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 202);

    auto g_yzzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 203);

    auto g_yzzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 204);

    auto g_yzzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 205);

    auto g_yzzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 206);

    auto g_yzzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 207);

    auto g_yzzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 208);

    auto g_yzzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 209);

    auto g_zzzz_xxxx_1 = pbuffer.data(idx_eri_1_gg + 210);

    auto g_zzzz_xxxy_1 = pbuffer.data(idx_eri_1_gg + 211);

    auto g_zzzz_xxxz_1 = pbuffer.data(idx_eri_1_gg + 212);

    auto g_zzzz_xxyy_1 = pbuffer.data(idx_eri_1_gg + 213);

    auto g_zzzz_xxyz_1 = pbuffer.data(idx_eri_1_gg + 214);

    auto g_zzzz_xxzz_1 = pbuffer.data(idx_eri_1_gg + 215);

    auto g_zzzz_xyyy_1 = pbuffer.data(idx_eri_1_gg + 216);

    auto g_zzzz_xyyz_1 = pbuffer.data(idx_eri_1_gg + 217);

    auto g_zzzz_xyzz_1 = pbuffer.data(idx_eri_1_gg + 218);

    auto g_zzzz_xzzz_1 = pbuffer.data(idx_eri_1_gg + 219);

    auto g_zzzz_yyyy_1 = pbuffer.data(idx_eri_1_gg + 220);

    auto g_zzzz_yyyz_1 = pbuffer.data(idx_eri_1_gg + 221);

    auto g_zzzz_yyzz_1 = pbuffer.data(idx_eri_1_gg + 222);

    auto g_zzzz_yzzz_1 = pbuffer.data(idx_eri_1_gg + 223);

    auto g_zzzz_zzzz_1 = pbuffer.data(idx_eri_1_gg + 224);

    // Set up components of auxiliary buffer : GH

    auto g_xxxx_xxxxx_1 = pbuffer.data(idx_eri_1_gh);

    auto g_xxxx_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 1);

    auto g_xxxx_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 2);

    auto g_xxxx_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 3);

    auto g_xxxx_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 4);

    auto g_xxxx_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 5);

    auto g_xxxx_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 6);

    auto g_xxxx_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 7);

    auto g_xxxx_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 8);

    auto g_xxxx_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 9);

    auto g_xxxx_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 10);

    auto g_xxxx_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 11);

    auto g_xxxx_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 12);

    auto g_xxxx_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 13);

    auto g_xxxx_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 14);

    auto g_xxxx_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 15);

    auto g_xxxx_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 16);

    auto g_xxxx_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 17);

    auto g_xxxx_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 18);

    auto g_xxxx_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 19);

    auto g_xxxx_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 20);

    auto g_xxxy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 21);

    auto g_xxxy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 22);

    auto g_xxxy_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 23);

    auto g_xxxy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 24);

    auto g_xxxy_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 26);

    auto g_xxxy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 27);

    auto g_xxxy_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 30);

    auto g_xxxy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 31);

    auto g_xxxy_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 35);

    auto g_xxxy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 36);

    auto g_xxxz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 42);

    auto g_xxxz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 43);

    auto g_xxxz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 44);

    auto g_xxxz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 45);

    auto g_xxxz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 46);

    auto g_xxxz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 47);

    auto g_xxxz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 48);

    auto g_xxxz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 49);

    auto g_xxxz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 50);

    auto g_xxxz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 51);

    auto g_xxxz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 52);

    auto g_xxxz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 53);

    auto g_xxxz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 54);

    auto g_xxxz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 55);

    auto g_xxxz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 56);

    auto g_xxxz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 58);

    auto g_xxxz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 59);

    auto g_xxxz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 60);

    auto g_xxxz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 61);

    auto g_xxxz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 62);

    auto g_xxyy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 63);

    auto g_xxyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 64);

    auto g_xxyy_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 65);

    auto g_xxyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 66);

    auto g_xxyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 67);

    auto g_xxyy_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 68);

    auto g_xxyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 69);

    auto g_xxyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 70);

    auto g_xxyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 71);

    auto g_xxyy_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 72);

    auto g_xxyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 73);

    auto g_xxyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 74);

    auto g_xxyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 75);

    auto g_xxyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 76);

    auto g_xxyy_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 77);

    auto g_xxyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 78);

    auto g_xxyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 79);

    auto g_xxyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 80);

    auto g_xxyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 81);

    auto g_xxyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 82);

    auto g_xxyy_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 83);

    auto g_xxzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 105);

    auto g_xxzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 106);

    auto g_xxzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 107);

    auto g_xxzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 108);

    auto g_xxzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 109);

    auto g_xxzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 110);

    auto g_xxzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 111);

    auto g_xxzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 112);

    auto g_xxzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 113);

    auto g_xxzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 114);

    auto g_xxzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 115);

    auto g_xxzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 116);

    auto g_xxzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 117);

    auto g_xxzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 118);

    auto g_xxzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 119);

    auto g_xxzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 120);

    auto g_xxzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 121);

    auto g_xxzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 122);

    auto g_xxzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 123);

    auto g_xxzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 124);

    auto g_xxzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 125);

    auto g_xyyy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 126);

    auto g_xyyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 127);

    auto g_xyyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 129);

    auto g_xyyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 130);

    auto g_xyyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 132);

    auto g_xyyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 133);

    auto g_xyyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 134);

    auto g_xyyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 136);

    auto g_xyyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 137);

    auto g_xyyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 138);

    auto g_xyyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 139);

    auto g_xyyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 141);

    auto g_xyyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 142);

    auto g_xyyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 143);

    auto g_xyyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 144);

    auto g_xyyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 145);

    auto g_xyyy_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 146);

    auto g_xzzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 189);

    auto g_xzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 191);

    auto g_xzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 193);

    auto g_xzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 194);

    auto g_xzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 196);

    auto g_xzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 197);

    auto g_xzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 198);

    auto g_xzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 200);

    auto g_xzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 201);

    auto g_xzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 202);

    auto g_xzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 203);

    auto g_xzzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 204);

    auto g_xzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 205);

    auto g_xzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 206);

    auto g_xzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 207);

    auto g_xzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 208);

    auto g_xzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 209);

    auto g_yyyy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 210);

    auto g_yyyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 211);

    auto g_yyyy_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 212);

    auto g_yyyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 213);

    auto g_yyyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 214);

    auto g_yyyy_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 215);

    auto g_yyyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 216);

    auto g_yyyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 217);

    auto g_yyyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 218);

    auto g_yyyy_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 219);

    auto g_yyyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 220);

    auto g_yyyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 221);

    auto g_yyyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 222);

    auto g_yyyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 223);

    auto g_yyyy_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 224);

    auto g_yyyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 225);

    auto g_yyyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 226);

    auto g_yyyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 227);

    auto g_yyyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 228);

    auto g_yyyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 229);

    auto g_yyyy_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 230);

    auto g_yyyz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 232);

    auto g_yyyz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 233);

    auto g_yyyz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 234);

    auto g_yyyz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 235);

    auto g_yyyz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 236);

    auto g_yyyz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 237);

    auto g_yyyz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 238);

    auto g_yyyz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 239);

    auto g_yyyz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 240);

    auto g_yyyz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 241);

    auto g_yyyz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 242);

    auto g_yyyz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 243);

    auto g_yyyz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 244);

    auto g_yyyz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 245);

    auto g_yyyz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 246);

    auto g_yyyz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 247);

    auto g_yyyz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 248);

    auto g_yyyz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 249);

    auto g_yyyz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 250);

    auto g_yyyz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 251);

    auto g_yyzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 252);

    auto g_yyzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 253);

    auto g_yyzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 254);

    auto g_yyzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 255);

    auto g_yyzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 256);

    auto g_yyzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 257);

    auto g_yyzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 258);

    auto g_yyzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 259);

    auto g_yyzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 260);

    auto g_yyzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 261);

    auto g_yyzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 262);

    auto g_yyzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 263);

    auto g_yyzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 264);

    auto g_yyzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 265);

    auto g_yyzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 266);

    auto g_yyzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 267);

    auto g_yyzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 268);

    auto g_yyzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 269);

    auto g_yyzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 270);

    auto g_yyzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 271);

    auto g_yyzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 272);

    auto g_yzzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 273);

    auto g_yzzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 274);

    auto g_yzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 275);

    auto g_yzzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 276);

    auto g_yzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 277);

    auto g_yzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 278);

    auto g_yzzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 279);

    auto g_yzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 280);

    auto g_yzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 281);

    auto g_yzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 282);

    auto g_yzzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 283);

    auto g_yzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 284);

    auto g_yzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 285);

    auto g_yzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 286);

    auto g_yzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 287);

    auto g_yzzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 288);

    auto g_yzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 289);

    auto g_yzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 290);

    auto g_yzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 291);

    auto g_yzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 292);

    auto g_yzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 293);

    auto g_zzzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 294);

    auto g_zzzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 295);

    auto g_zzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 296);

    auto g_zzzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 297);

    auto g_zzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 298);

    auto g_zzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 299);

    auto g_zzzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 300);

    auto g_zzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 301);

    auto g_zzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 302);

    auto g_zzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 303);

    auto g_zzzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 304);

    auto g_zzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 305);

    auto g_zzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 306);

    auto g_zzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 307);

    auto g_zzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 308);

    auto g_zzzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 309);

    auto g_zzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 310);

    auto g_zzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 311);

    auto g_zzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 312);

    auto g_zzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 313);

    auto g_zzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 314);

    // Set up 0-21 components of targeted buffer : HH

    auto g_xxxxx_xxxxx_0 = pbuffer.data(idx_eri_0_hh);

    auto g_xxxxx_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 1);

    auto g_xxxxx_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 2);

    auto g_xxxxx_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 3);

    auto g_xxxxx_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 4);

    auto g_xxxxx_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 5);

    auto g_xxxxx_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 6);

    auto g_xxxxx_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 7);

    auto g_xxxxx_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 8);

    auto g_xxxxx_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 9);

    auto g_xxxxx_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 10);

    auto g_xxxxx_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 11);

    auto g_xxxxx_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 12);

    auto g_xxxxx_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 13);

    auto g_xxxxx_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 14);

    auto g_xxxxx_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 15);

    auto g_xxxxx_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 16);

    auto g_xxxxx_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 17);

    auto g_xxxxx_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 18);

    auto g_xxxxx_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 19);

    auto g_xxxxx_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 20);

    #pragma omp simd aligned(g_xxx_xxxxx_0, g_xxx_xxxxx_1, g_xxx_xxxxy_0, g_xxx_xxxxy_1, g_xxx_xxxxz_0, g_xxx_xxxxz_1, g_xxx_xxxyy_0, g_xxx_xxxyy_1, g_xxx_xxxyz_0, g_xxx_xxxyz_1, g_xxx_xxxzz_0, g_xxx_xxxzz_1, g_xxx_xxyyy_0, g_xxx_xxyyy_1, g_xxx_xxyyz_0, g_xxx_xxyyz_1, g_xxx_xxyzz_0, g_xxx_xxyzz_1, g_xxx_xxzzz_0, g_xxx_xxzzz_1, g_xxx_xyyyy_0, g_xxx_xyyyy_1, g_xxx_xyyyz_0, g_xxx_xyyyz_1, g_xxx_xyyzz_0, g_xxx_xyyzz_1, g_xxx_xyzzz_0, g_xxx_xyzzz_1, g_xxx_xzzzz_0, g_xxx_xzzzz_1, g_xxx_yyyyy_0, g_xxx_yyyyy_1, g_xxx_yyyyz_0, g_xxx_yyyyz_1, g_xxx_yyyzz_0, g_xxx_yyyzz_1, g_xxx_yyzzz_0, g_xxx_yyzzz_1, g_xxx_yzzzz_0, g_xxx_yzzzz_1, g_xxx_zzzzz_0, g_xxx_zzzzz_1, g_xxxx_xxxx_1, g_xxxx_xxxxx_1, g_xxxx_xxxxy_1, g_xxxx_xxxxz_1, g_xxxx_xxxy_1, g_xxxx_xxxyy_1, g_xxxx_xxxyz_1, g_xxxx_xxxz_1, g_xxxx_xxxzz_1, g_xxxx_xxyy_1, g_xxxx_xxyyy_1, g_xxxx_xxyyz_1, g_xxxx_xxyz_1, g_xxxx_xxyzz_1, g_xxxx_xxzz_1, g_xxxx_xxzzz_1, g_xxxx_xyyy_1, g_xxxx_xyyyy_1, g_xxxx_xyyyz_1, g_xxxx_xyyz_1, g_xxxx_xyyzz_1, g_xxxx_xyzz_1, g_xxxx_xyzzz_1, g_xxxx_xzzz_1, g_xxxx_xzzzz_1, g_xxxx_yyyy_1, g_xxxx_yyyyy_1, g_xxxx_yyyyz_1, g_xxxx_yyyz_1, g_xxxx_yyyzz_1, g_xxxx_yyzz_1, g_xxxx_yyzzz_1, g_xxxx_yzzz_1, g_xxxx_yzzzz_1, g_xxxx_zzzz_1, g_xxxx_zzzzz_1, g_xxxxx_xxxxx_0, g_xxxxx_xxxxy_0, g_xxxxx_xxxxz_0, g_xxxxx_xxxyy_0, g_xxxxx_xxxyz_0, g_xxxxx_xxxzz_0, g_xxxxx_xxyyy_0, g_xxxxx_xxyyz_0, g_xxxxx_xxyzz_0, g_xxxxx_xxzzz_0, g_xxxxx_xyyyy_0, g_xxxxx_xyyyz_0, g_xxxxx_xyyzz_0, g_xxxxx_xyzzz_0, g_xxxxx_xzzzz_0, g_xxxxx_yyyyy_0, g_xxxxx_yyyyz_0, g_xxxxx_yyyzz_0, g_xxxxx_yyzzz_0, g_xxxxx_yzzzz_0, g_xxxxx_zzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxx_xxxxx_0[i] = 4.0 * g_xxx_xxxxx_0[i] * fbe_0 - 4.0 * g_xxx_xxxxx_1[i] * fz_be_0 + 5.0 * g_xxxx_xxxx_1[i] * fe_0 + g_xxxx_xxxxx_1[i] * pa_x[i];

        g_xxxxx_xxxxy_0[i] = 4.0 * g_xxx_xxxxy_0[i] * fbe_0 - 4.0 * g_xxx_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxx_xxxy_1[i] * fe_0 + g_xxxx_xxxxy_1[i] * pa_x[i];

        g_xxxxx_xxxxz_0[i] = 4.0 * g_xxx_xxxxz_0[i] * fbe_0 - 4.0 * g_xxx_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxx_xxxz_1[i] * fe_0 + g_xxxx_xxxxz_1[i] * pa_x[i];

        g_xxxxx_xxxyy_0[i] = 4.0 * g_xxx_xxxyy_0[i] * fbe_0 - 4.0 * g_xxx_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxx_xxyy_1[i] * fe_0 + g_xxxx_xxxyy_1[i] * pa_x[i];

        g_xxxxx_xxxyz_0[i] = 4.0 * g_xxx_xxxyz_0[i] * fbe_0 - 4.0 * g_xxx_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxx_xxyz_1[i] * fe_0 + g_xxxx_xxxyz_1[i] * pa_x[i];

        g_xxxxx_xxxzz_0[i] = 4.0 * g_xxx_xxxzz_0[i] * fbe_0 - 4.0 * g_xxx_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxx_xxzz_1[i] * fe_0 + g_xxxx_xxxzz_1[i] * pa_x[i];

        g_xxxxx_xxyyy_0[i] = 4.0 * g_xxx_xxyyy_0[i] * fbe_0 - 4.0 * g_xxx_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxx_xyyy_1[i] * fe_0 + g_xxxx_xxyyy_1[i] * pa_x[i];

        g_xxxxx_xxyyz_0[i] = 4.0 * g_xxx_xxyyz_0[i] * fbe_0 - 4.0 * g_xxx_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxx_xyyz_1[i] * fe_0 + g_xxxx_xxyyz_1[i] * pa_x[i];

        g_xxxxx_xxyzz_0[i] = 4.0 * g_xxx_xxyzz_0[i] * fbe_0 - 4.0 * g_xxx_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxx_xyzz_1[i] * fe_0 + g_xxxx_xxyzz_1[i] * pa_x[i];

        g_xxxxx_xxzzz_0[i] = 4.0 * g_xxx_xxzzz_0[i] * fbe_0 - 4.0 * g_xxx_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_xzzz_1[i] * fe_0 + g_xxxx_xxzzz_1[i] * pa_x[i];

        g_xxxxx_xyyyy_0[i] = 4.0 * g_xxx_xyyyy_0[i] * fbe_0 - 4.0 * g_xxx_xyyyy_1[i] * fz_be_0 + g_xxxx_yyyy_1[i] * fe_0 + g_xxxx_xyyyy_1[i] * pa_x[i];

        g_xxxxx_xyyyz_0[i] = 4.0 * g_xxx_xyyyz_0[i] * fbe_0 - 4.0 * g_xxx_xyyyz_1[i] * fz_be_0 + g_xxxx_yyyz_1[i] * fe_0 + g_xxxx_xyyyz_1[i] * pa_x[i];

        g_xxxxx_xyyzz_0[i] = 4.0 * g_xxx_xyyzz_0[i] * fbe_0 - 4.0 * g_xxx_xyyzz_1[i] * fz_be_0 + g_xxxx_yyzz_1[i] * fe_0 + g_xxxx_xyyzz_1[i] * pa_x[i];

        g_xxxxx_xyzzz_0[i] = 4.0 * g_xxx_xyzzz_0[i] * fbe_0 - 4.0 * g_xxx_xyzzz_1[i] * fz_be_0 + g_xxxx_yzzz_1[i] * fe_0 + g_xxxx_xyzzz_1[i] * pa_x[i];

        g_xxxxx_xzzzz_0[i] = 4.0 * g_xxx_xzzzz_0[i] * fbe_0 - 4.0 * g_xxx_xzzzz_1[i] * fz_be_0 + g_xxxx_zzzz_1[i] * fe_0 + g_xxxx_xzzzz_1[i] * pa_x[i];

        g_xxxxx_yyyyy_0[i] = 4.0 * g_xxx_yyyyy_0[i] * fbe_0 - 4.0 * g_xxx_yyyyy_1[i] * fz_be_0 + g_xxxx_yyyyy_1[i] * pa_x[i];

        g_xxxxx_yyyyz_0[i] = 4.0 * g_xxx_yyyyz_0[i] * fbe_0 - 4.0 * g_xxx_yyyyz_1[i] * fz_be_0 + g_xxxx_yyyyz_1[i] * pa_x[i];

        g_xxxxx_yyyzz_0[i] = 4.0 * g_xxx_yyyzz_0[i] * fbe_0 - 4.0 * g_xxx_yyyzz_1[i] * fz_be_0 + g_xxxx_yyyzz_1[i] * pa_x[i];

        g_xxxxx_yyzzz_0[i] = 4.0 * g_xxx_yyzzz_0[i] * fbe_0 - 4.0 * g_xxx_yyzzz_1[i] * fz_be_0 + g_xxxx_yyzzz_1[i] * pa_x[i];

        g_xxxxx_yzzzz_0[i] = 4.0 * g_xxx_yzzzz_0[i] * fbe_0 - 4.0 * g_xxx_yzzzz_1[i] * fz_be_0 + g_xxxx_yzzzz_1[i] * pa_x[i];

        g_xxxxx_zzzzz_0[i] = 4.0 * g_xxx_zzzzz_0[i] * fbe_0 - 4.0 * g_xxx_zzzzz_1[i] * fz_be_0 + g_xxxx_zzzzz_1[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : HH

    auto g_xxxxy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 21);

    auto g_xxxxy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 22);

    auto g_xxxxy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 23);

    auto g_xxxxy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 24);

    auto g_xxxxy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 25);

    auto g_xxxxy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 26);

    auto g_xxxxy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 27);

    auto g_xxxxy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 28);

    auto g_xxxxy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 29);

    auto g_xxxxy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 30);

    auto g_xxxxy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 31);

    auto g_xxxxy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 32);

    auto g_xxxxy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 33);

    auto g_xxxxy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 34);

    auto g_xxxxy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 35);

    auto g_xxxxy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 36);

    auto g_xxxxy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 37);

    auto g_xxxxy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 38);

    auto g_xxxxy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 39);

    auto g_xxxxy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 40);

    auto g_xxxxy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 41);

    #pragma omp simd aligned(g_xxxx_xxxx_1, g_xxxx_xxxxx_1, g_xxxx_xxxxy_1, g_xxxx_xxxxz_1, g_xxxx_xxxy_1, g_xxxx_xxxyy_1, g_xxxx_xxxyz_1, g_xxxx_xxxz_1, g_xxxx_xxxzz_1, g_xxxx_xxyy_1, g_xxxx_xxyyy_1, g_xxxx_xxyyz_1, g_xxxx_xxyz_1, g_xxxx_xxyzz_1, g_xxxx_xxzz_1, g_xxxx_xxzzz_1, g_xxxx_xyyy_1, g_xxxx_xyyyy_1, g_xxxx_xyyyz_1, g_xxxx_xyyz_1, g_xxxx_xyyzz_1, g_xxxx_xyzz_1, g_xxxx_xyzzz_1, g_xxxx_xzzz_1, g_xxxx_xzzzz_1, g_xxxx_yyyy_1, g_xxxx_yyyyy_1, g_xxxx_yyyyz_1, g_xxxx_yyyz_1, g_xxxx_yyyzz_1, g_xxxx_yyzz_1, g_xxxx_yyzzz_1, g_xxxx_yzzz_1, g_xxxx_yzzzz_1, g_xxxx_zzzz_1, g_xxxx_zzzzz_1, g_xxxxy_xxxxx_0, g_xxxxy_xxxxy_0, g_xxxxy_xxxxz_0, g_xxxxy_xxxyy_0, g_xxxxy_xxxyz_0, g_xxxxy_xxxzz_0, g_xxxxy_xxyyy_0, g_xxxxy_xxyyz_0, g_xxxxy_xxyzz_0, g_xxxxy_xxzzz_0, g_xxxxy_xyyyy_0, g_xxxxy_xyyyz_0, g_xxxxy_xyyzz_0, g_xxxxy_xyzzz_0, g_xxxxy_xzzzz_0, g_xxxxy_yyyyy_0, g_xxxxy_yyyyz_0, g_xxxxy_yyyzz_0, g_xxxxy_yyzzz_0, g_xxxxy_yzzzz_0, g_xxxxy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxy_xxxxx_0[i] = g_xxxx_xxxxx_1[i] * pa_y[i];

        g_xxxxy_xxxxy_0[i] = g_xxxx_xxxx_1[i] * fe_0 + g_xxxx_xxxxy_1[i] * pa_y[i];

        g_xxxxy_xxxxz_0[i] = g_xxxx_xxxxz_1[i] * pa_y[i];

        g_xxxxy_xxxyy_0[i] = 2.0 * g_xxxx_xxxy_1[i] * fe_0 + g_xxxx_xxxyy_1[i] * pa_y[i];

        g_xxxxy_xxxyz_0[i] = g_xxxx_xxxz_1[i] * fe_0 + g_xxxx_xxxyz_1[i] * pa_y[i];

        g_xxxxy_xxxzz_0[i] = g_xxxx_xxxzz_1[i] * pa_y[i];

        g_xxxxy_xxyyy_0[i] = 3.0 * g_xxxx_xxyy_1[i] * fe_0 + g_xxxx_xxyyy_1[i] * pa_y[i];

        g_xxxxy_xxyyz_0[i] = 2.0 * g_xxxx_xxyz_1[i] * fe_0 + g_xxxx_xxyyz_1[i] * pa_y[i];

        g_xxxxy_xxyzz_0[i] = g_xxxx_xxzz_1[i] * fe_0 + g_xxxx_xxyzz_1[i] * pa_y[i];

        g_xxxxy_xxzzz_0[i] = g_xxxx_xxzzz_1[i] * pa_y[i];

        g_xxxxy_xyyyy_0[i] = 4.0 * g_xxxx_xyyy_1[i] * fe_0 + g_xxxx_xyyyy_1[i] * pa_y[i];

        g_xxxxy_xyyyz_0[i] = 3.0 * g_xxxx_xyyz_1[i] * fe_0 + g_xxxx_xyyyz_1[i] * pa_y[i];

        g_xxxxy_xyyzz_0[i] = 2.0 * g_xxxx_xyzz_1[i] * fe_0 + g_xxxx_xyyzz_1[i] * pa_y[i];

        g_xxxxy_xyzzz_0[i] = g_xxxx_xzzz_1[i] * fe_0 + g_xxxx_xyzzz_1[i] * pa_y[i];

        g_xxxxy_xzzzz_0[i] = g_xxxx_xzzzz_1[i] * pa_y[i];

        g_xxxxy_yyyyy_0[i] = 5.0 * g_xxxx_yyyy_1[i] * fe_0 + g_xxxx_yyyyy_1[i] * pa_y[i];

        g_xxxxy_yyyyz_0[i] = 4.0 * g_xxxx_yyyz_1[i] * fe_0 + g_xxxx_yyyyz_1[i] * pa_y[i];

        g_xxxxy_yyyzz_0[i] = 3.0 * g_xxxx_yyzz_1[i] * fe_0 + g_xxxx_yyyzz_1[i] * pa_y[i];

        g_xxxxy_yyzzz_0[i] = 2.0 * g_xxxx_yzzz_1[i] * fe_0 + g_xxxx_yyzzz_1[i] * pa_y[i];

        g_xxxxy_yzzzz_0[i] = g_xxxx_zzzz_1[i] * fe_0 + g_xxxx_yzzzz_1[i] * pa_y[i];

        g_xxxxy_zzzzz_0[i] = g_xxxx_zzzzz_1[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : HH

    auto g_xxxxz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 42);

    auto g_xxxxz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 43);

    auto g_xxxxz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 44);

    auto g_xxxxz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 45);

    auto g_xxxxz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 46);

    auto g_xxxxz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 47);

    auto g_xxxxz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 48);

    auto g_xxxxz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 49);

    auto g_xxxxz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 50);

    auto g_xxxxz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 51);

    auto g_xxxxz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 52);

    auto g_xxxxz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 53);

    auto g_xxxxz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 54);

    auto g_xxxxz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 55);

    auto g_xxxxz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 56);

    auto g_xxxxz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 57);

    auto g_xxxxz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 58);

    auto g_xxxxz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 59);

    auto g_xxxxz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 60);

    auto g_xxxxz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 61);

    auto g_xxxxz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 62);

    #pragma omp simd aligned(g_xxxx_xxxx_1, g_xxxx_xxxxx_1, g_xxxx_xxxxy_1, g_xxxx_xxxxz_1, g_xxxx_xxxy_1, g_xxxx_xxxyy_1, g_xxxx_xxxyz_1, g_xxxx_xxxz_1, g_xxxx_xxxzz_1, g_xxxx_xxyy_1, g_xxxx_xxyyy_1, g_xxxx_xxyyz_1, g_xxxx_xxyz_1, g_xxxx_xxyzz_1, g_xxxx_xxzz_1, g_xxxx_xxzzz_1, g_xxxx_xyyy_1, g_xxxx_xyyyy_1, g_xxxx_xyyyz_1, g_xxxx_xyyz_1, g_xxxx_xyyzz_1, g_xxxx_xyzz_1, g_xxxx_xyzzz_1, g_xxxx_xzzz_1, g_xxxx_xzzzz_1, g_xxxx_yyyy_1, g_xxxx_yyyyy_1, g_xxxx_yyyyz_1, g_xxxx_yyyz_1, g_xxxx_yyyzz_1, g_xxxx_yyzz_1, g_xxxx_yyzzz_1, g_xxxx_yzzz_1, g_xxxx_yzzzz_1, g_xxxx_zzzz_1, g_xxxx_zzzzz_1, g_xxxxz_xxxxx_0, g_xxxxz_xxxxy_0, g_xxxxz_xxxxz_0, g_xxxxz_xxxyy_0, g_xxxxz_xxxyz_0, g_xxxxz_xxxzz_0, g_xxxxz_xxyyy_0, g_xxxxz_xxyyz_0, g_xxxxz_xxyzz_0, g_xxxxz_xxzzz_0, g_xxxxz_xyyyy_0, g_xxxxz_xyyyz_0, g_xxxxz_xyyzz_0, g_xxxxz_xyzzz_0, g_xxxxz_xzzzz_0, g_xxxxz_yyyyy_0, g_xxxxz_yyyyz_0, g_xxxxz_yyyzz_0, g_xxxxz_yyzzz_0, g_xxxxz_yzzzz_0, g_xxxxz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxz_xxxxx_0[i] = g_xxxx_xxxxx_1[i] * pa_z[i];

        g_xxxxz_xxxxy_0[i] = g_xxxx_xxxxy_1[i] * pa_z[i];

        g_xxxxz_xxxxz_0[i] = g_xxxx_xxxx_1[i] * fe_0 + g_xxxx_xxxxz_1[i] * pa_z[i];

        g_xxxxz_xxxyy_0[i] = g_xxxx_xxxyy_1[i] * pa_z[i];

        g_xxxxz_xxxyz_0[i] = g_xxxx_xxxy_1[i] * fe_0 + g_xxxx_xxxyz_1[i] * pa_z[i];

        g_xxxxz_xxxzz_0[i] = 2.0 * g_xxxx_xxxz_1[i] * fe_0 + g_xxxx_xxxzz_1[i] * pa_z[i];

        g_xxxxz_xxyyy_0[i] = g_xxxx_xxyyy_1[i] * pa_z[i];

        g_xxxxz_xxyyz_0[i] = g_xxxx_xxyy_1[i] * fe_0 + g_xxxx_xxyyz_1[i] * pa_z[i];

        g_xxxxz_xxyzz_0[i] = 2.0 * g_xxxx_xxyz_1[i] * fe_0 + g_xxxx_xxyzz_1[i] * pa_z[i];

        g_xxxxz_xxzzz_0[i] = 3.0 * g_xxxx_xxzz_1[i] * fe_0 + g_xxxx_xxzzz_1[i] * pa_z[i];

        g_xxxxz_xyyyy_0[i] = g_xxxx_xyyyy_1[i] * pa_z[i];

        g_xxxxz_xyyyz_0[i] = g_xxxx_xyyy_1[i] * fe_0 + g_xxxx_xyyyz_1[i] * pa_z[i];

        g_xxxxz_xyyzz_0[i] = 2.0 * g_xxxx_xyyz_1[i] * fe_0 + g_xxxx_xyyzz_1[i] * pa_z[i];

        g_xxxxz_xyzzz_0[i] = 3.0 * g_xxxx_xyzz_1[i] * fe_0 + g_xxxx_xyzzz_1[i] * pa_z[i];

        g_xxxxz_xzzzz_0[i] = 4.0 * g_xxxx_xzzz_1[i] * fe_0 + g_xxxx_xzzzz_1[i] * pa_z[i];

        g_xxxxz_yyyyy_0[i] = g_xxxx_yyyyy_1[i] * pa_z[i];

        g_xxxxz_yyyyz_0[i] = g_xxxx_yyyy_1[i] * fe_0 + g_xxxx_yyyyz_1[i] * pa_z[i];

        g_xxxxz_yyyzz_0[i] = 2.0 * g_xxxx_yyyz_1[i] * fe_0 + g_xxxx_yyyzz_1[i] * pa_z[i];

        g_xxxxz_yyzzz_0[i] = 3.0 * g_xxxx_yyzz_1[i] * fe_0 + g_xxxx_yyzzz_1[i] * pa_z[i];

        g_xxxxz_yzzzz_0[i] = 4.0 * g_xxxx_yzzz_1[i] * fe_0 + g_xxxx_yzzzz_1[i] * pa_z[i];

        g_xxxxz_zzzzz_0[i] = 5.0 * g_xxxx_zzzz_1[i] * fe_0 + g_xxxx_zzzzz_1[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : HH

    auto g_xxxyy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 63);

    auto g_xxxyy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 64);

    auto g_xxxyy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 65);

    auto g_xxxyy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 66);

    auto g_xxxyy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 67);

    auto g_xxxyy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 68);

    auto g_xxxyy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 69);

    auto g_xxxyy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 70);

    auto g_xxxyy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 71);

    auto g_xxxyy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 72);

    auto g_xxxyy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 73);

    auto g_xxxyy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 74);

    auto g_xxxyy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 75);

    auto g_xxxyy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 76);

    auto g_xxxyy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 77);

    auto g_xxxyy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 78);

    auto g_xxxyy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 79);

    auto g_xxxyy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 80);

    auto g_xxxyy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 81);

    auto g_xxxyy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 82);

    auto g_xxxyy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 83);

    #pragma omp simd aligned(g_xxx_xxxxx_0, g_xxx_xxxxx_1, g_xxx_xxxxz_0, g_xxx_xxxxz_1, g_xxx_xxxzz_0, g_xxx_xxxzz_1, g_xxx_xxzzz_0, g_xxx_xxzzz_1, g_xxx_xzzzz_0, g_xxx_xzzzz_1, g_xxxy_xxxxx_1, g_xxxy_xxxxz_1, g_xxxy_xxxzz_1, g_xxxy_xxzzz_1, g_xxxy_xzzzz_1, g_xxxyy_xxxxx_0, g_xxxyy_xxxxy_0, g_xxxyy_xxxxz_0, g_xxxyy_xxxyy_0, g_xxxyy_xxxyz_0, g_xxxyy_xxxzz_0, g_xxxyy_xxyyy_0, g_xxxyy_xxyyz_0, g_xxxyy_xxyzz_0, g_xxxyy_xxzzz_0, g_xxxyy_xyyyy_0, g_xxxyy_xyyyz_0, g_xxxyy_xyyzz_0, g_xxxyy_xyzzz_0, g_xxxyy_xzzzz_0, g_xxxyy_yyyyy_0, g_xxxyy_yyyyz_0, g_xxxyy_yyyzz_0, g_xxxyy_yyzzz_0, g_xxxyy_yzzzz_0, g_xxxyy_zzzzz_0, g_xxyy_xxxxy_1, g_xxyy_xxxy_1, g_xxyy_xxxyy_1, g_xxyy_xxxyz_1, g_xxyy_xxyy_1, g_xxyy_xxyyy_1, g_xxyy_xxyyz_1, g_xxyy_xxyz_1, g_xxyy_xxyzz_1, g_xxyy_xyyy_1, g_xxyy_xyyyy_1, g_xxyy_xyyyz_1, g_xxyy_xyyz_1, g_xxyy_xyyzz_1, g_xxyy_xyzz_1, g_xxyy_xyzzz_1, g_xxyy_yyyy_1, g_xxyy_yyyyy_1, g_xxyy_yyyyz_1, g_xxyy_yyyz_1, g_xxyy_yyyzz_1, g_xxyy_yyzz_1, g_xxyy_yyzzz_1, g_xxyy_yzzz_1, g_xxyy_yzzzz_1, g_xxyy_zzzzz_1, g_xyy_xxxxy_0, g_xyy_xxxxy_1, g_xyy_xxxyy_0, g_xyy_xxxyy_1, g_xyy_xxxyz_0, g_xyy_xxxyz_1, g_xyy_xxyyy_0, g_xyy_xxyyy_1, g_xyy_xxyyz_0, g_xyy_xxyyz_1, g_xyy_xxyzz_0, g_xyy_xxyzz_1, g_xyy_xyyyy_0, g_xyy_xyyyy_1, g_xyy_xyyyz_0, g_xyy_xyyyz_1, g_xyy_xyyzz_0, g_xyy_xyyzz_1, g_xyy_xyzzz_0, g_xyy_xyzzz_1, g_xyy_yyyyy_0, g_xyy_yyyyy_1, g_xyy_yyyyz_0, g_xyy_yyyyz_1, g_xyy_yyyzz_0, g_xyy_yyyzz_1, g_xyy_yyzzz_0, g_xyy_yyzzz_1, g_xyy_yzzzz_0, g_xyy_yzzzz_1, g_xyy_zzzzz_0, g_xyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyy_xxxxx_0[i] = g_xxx_xxxxx_0[i] * fbe_0 - g_xxx_xxxxx_1[i] * fz_be_0 + g_xxxy_xxxxx_1[i] * pa_y[i];

        g_xxxyy_xxxxy_0[i] = 2.0 * g_xyy_xxxxy_0[i] * fbe_0 - 2.0 * g_xyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxyy_xxxy_1[i] * fe_0 + g_xxyy_xxxxy_1[i] * pa_x[i];

        g_xxxyy_xxxxz_0[i] = g_xxx_xxxxz_0[i] * fbe_0 - g_xxx_xxxxz_1[i] * fz_be_0 + g_xxxy_xxxxz_1[i] * pa_y[i];

        g_xxxyy_xxxyy_0[i] = 2.0 * g_xyy_xxxyy_0[i] * fbe_0 - 2.0 * g_xyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxyy_xxyy_1[i] * fe_0 + g_xxyy_xxxyy_1[i] * pa_x[i];

        g_xxxyy_xxxyz_0[i] = 2.0 * g_xyy_xxxyz_0[i] * fbe_0 - 2.0 * g_xyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxyy_xxyz_1[i] * fe_0 + g_xxyy_xxxyz_1[i] * pa_x[i];

        g_xxxyy_xxxzz_0[i] = g_xxx_xxxzz_0[i] * fbe_0 - g_xxx_xxxzz_1[i] * fz_be_0 + g_xxxy_xxxzz_1[i] * pa_y[i];

        g_xxxyy_xxyyy_0[i] = 2.0 * g_xyy_xxyyy_0[i] * fbe_0 - 2.0 * g_xyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxyy_xyyy_1[i] * fe_0 + g_xxyy_xxyyy_1[i] * pa_x[i];

        g_xxxyy_xxyyz_0[i] = 2.0 * g_xyy_xxyyz_0[i] * fbe_0 - 2.0 * g_xyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxyy_xyyz_1[i] * fe_0 + g_xxyy_xxyyz_1[i] * pa_x[i];

        g_xxxyy_xxyzz_0[i] = 2.0 * g_xyy_xxyzz_0[i] * fbe_0 - 2.0 * g_xyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxyy_xyzz_1[i] * fe_0 + g_xxyy_xxyzz_1[i] * pa_x[i];

        g_xxxyy_xxzzz_0[i] = g_xxx_xxzzz_0[i] * fbe_0 - g_xxx_xxzzz_1[i] * fz_be_0 + g_xxxy_xxzzz_1[i] * pa_y[i];

        g_xxxyy_xyyyy_0[i] = 2.0 * g_xyy_xyyyy_0[i] * fbe_0 - 2.0 * g_xyy_xyyyy_1[i] * fz_be_0 + g_xxyy_yyyy_1[i] * fe_0 + g_xxyy_xyyyy_1[i] * pa_x[i];

        g_xxxyy_xyyyz_0[i] = 2.0 * g_xyy_xyyyz_0[i] * fbe_0 - 2.0 * g_xyy_xyyyz_1[i] * fz_be_0 + g_xxyy_yyyz_1[i] * fe_0 + g_xxyy_xyyyz_1[i] * pa_x[i];

        g_xxxyy_xyyzz_0[i] = 2.0 * g_xyy_xyyzz_0[i] * fbe_0 - 2.0 * g_xyy_xyyzz_1[i] * fz_be_0 + g_xxyy_yyzz_1[i] * fe_0 + g_xxyy_xyyzz_1[i] * pa_x[i];

        g_xxxyy_xyzzz_0[i] = 2.0 * g_xyy_xyzzz_0[i] * fbe_0 - 2.0 * g_xyy_xyzzz_1[i] * fz_be_0 + g_xxyy_yzzz_1[i] * fe_0 + g_xxyy_xyzzz_1[i] * pa_x[i];

        g_xxxyy_xzzzz_0[i] = g_xxx_xzzzz_0[i] * fbe_0 - g_xxx_xzzzz_1[i] * fz_be_0 + g_xxxy_xzzzz_1[i] * pa_y[i];

        g_xxxyy_yyyyy_0[i] = 2.0 * g_xyy_yyyyy_0[i] * fbe_0 - 2.0 * g_xyy_yyyyy_1[i] * fz_be_0 + g_xxyy_yyyyy_1[i] * pa_x[i];

        g_xxxyy_yyyyz_0[i] = 2.0 * g_xyy_yyyyz_0[i] * fbe_0 - 2.0 * g_xyy_yyyyz_1[i] * fz_be_0 + g_xxyy_yyyyz_1[i] * pa_x[i];

        g_xxxyy_yyyzz_0[i] = 2.0 * g_xyy_yyyzz_0[i] * fbe_0 - 2.0 * g_xyy_yyyzz_1[i] * fz_be_0 + g_xxyy_yyyzz_1[i] * pa_x[i];

        g_xxxyy_yyzzz_0[i] = 2.0 * g_xyy_yyzzz_0[i] * fbe_0 - 2.0 * g_xyy_yyzzz_1[i] * fz_be_0 + g_xxyy_yyzzz_1[i] * pa_x[i];

        g_xxxyy_yzzzz_0[i] = 2.0 * g_xyy_yzzzz_0[i] * fbe_0 - 2.0 * g_xyy_yzzzz_1[i] * fz_be_0 + g_xxyy_yzzzz_1[i] * pa_x[i];

        g_xxxyy_zzzzz_0[i] = 2.0 * g_xyy_zzzzz_0[i] * fbe_0 - 2.0 * g_xyy_zzzzz_1[i] * fz_be_0 + g_xxyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : HH

    auto g_xxxyz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 84);

    auto g_xxxyz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 85);

    auto g_xxxyz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 86);

    auto g_xxxyz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 87);

    auto g_xxxyz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 88);

    auto g_xxxyz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 89);

    auto g_xxxyz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 90);

    auto g_xxxyz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 91);

    auto g_xxxyz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 92);

    auto g_xxxyz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 93);

    auto g_xxxyz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 94);

    auto g_xxxyz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 95);

    auto g_xxxyz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 96);

    auto g_xxxyz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 97);

    auto g_xxxyz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 98);

    auto g_xxxyz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 99);

    auto g_xxxyz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 100);

    auto g_xxxyz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 101);

    auto g_xxxyz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 102);

    auto g_xxxyz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 103);

    auto g_xxxyz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 104);

    #pragma omp simd aligned(g_xxxy_xxxxy_1, g_xxxy_xxxyy_1, g_xxxy_xxyyy_1, g_xxxy_xyyyy_1, g_xxxy_yyyyy_1, g_xxxyz_xxxxx_0, g_xxxyz_xxxxy_0, g_xxxyz_xxxxz_0, g_xxxyz_xxxyy_0, g_xxxyz_xxxyz_0, g_xxxyz_xxxzz_0, g_xxxyz_xxyyy_0, g_xxxyz_xxyyz_0, g_xxxyz_xxyzz_0, g_xxxyz_xxzzz_0, g_xxxyz_xyyyy_0, g_xxxyz_xyyyz_0, g_xxxyz_xyyzz_0, g_xxxyz_xyzzz_0, g_xxxyz_xzzzz_0, g_xxxyz_yyyyy_0, g_xxxyz_yyyyz_0, g_xxxyz_yyyzz_0, g_xxxyz_yyzzz_0, g_xxxyz_yzzzz_0, g_xxxyz_zzzzz_0, g_xxxz_xxxxx_1, g_xxxz_xxxxz_1, g_xxxz_xxxyz_1, g_xxxz_xxxz_1, g_xxxz_xxxzz_1, g_xxxz_xxyyz_1, g_xxxz_xxyz_1, g_xxxz_xxyzz_1, g_xxxz_xxzz_1, g_xxxz_xxzzz_1, g_xxxz_xyyyz_1, g_xxxz_xyyz_1, g_xxxz_xyyzz_1, g_xxxz_xyzz_1, g_xxxz_xyzzz_1, g_xxxz_xzzz_1, g_xxxz_xzzzz_1, g_xxxz_yyyyz_1, g_xxxz_yyyz_1, g_xxxz_yyyzz_1, g_xxxz_yyzz_1, g_xxxz_yyzzz_1, g_xxxz_yzzz_1, g_xxxz_yzzzz_1, g_xxxz_zzzz_1, g_xxxz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyz_xxxxx_0[i] = g_xxxz_xxxxx_1[i] * pa_y[i];

        g_xxxyz_xxxxy_0[i] = g_xxxy_xxxxy_1[i] * pa_z[i];

        g_xxxyz_xxxxz_0[i] = g_xxxz_xxxxz_1[i] * pa_y[i];

        g_xxxyz_xxxyy_0[i] = g_xxxy_xxxyy_1[i] * pa_z[i];

        g_xxxyz_xxxyz_0[i] = g_xxxz_xxxz_1[i] * fe_0 + g_xxxz_xxxyz_1[i] * pa_y[i];

        g_xxxyz_xxxzz_0[i] = g_xxxz_xxxzz_1[i] * pa_y[i];

        g_xxxyz_xxyyy_0[i] = g_xxxy_xxyyy_1[i] * pa_z[i];

        g_xxxyz_xxyyz_0[i] = 2.0 * g_xxxz_xxyz_1[i] * fe_0 + g_xxxz_xxyyz_1[i] * pa_y[i];

        g_xxxyz_xxyzz_0[i] = g_xxxz_xxzz_1[i] * fe_0 + g_xxxz_xxyzz_1[i] * pa_y[i];

        g_xxxyz_xxzzz_0[i] = g_xxxz_xxzzz_1[i] * pa_y[i];

        g_xxxyz_xyyyy_0[i] = g_xxxy_xyyyy_1[i] * pa_z[i];

        g_xxxyz_xyyyz_0[i] = 3.0 * g_xxxz_xyyz_1[i] * fe_0 + g_xxxz_xyyyz_1[i] * pa_y[i];

        g_xxxyz_xyyzz_0[i] = 2.0 * g_xxxz_xyzz_1[i] * fe_0 + g_xxxz_xyyzz_1[i] * pa_y[i];

        g_xxxyz_xyzzz_0[i] = g_xxxz_xzzz_1[i] * fe_0 + g_xxxz_xyzzz_1[i] * pa_y[i];

        g_xxxyz_xzzzz_0[i] = g_xxxz_xzzzz_1[i] * pa_y[i];

        g_xxxyz_yyyyy_0[i] = g_xxxy_yyyyy_1[i] * pa_z[i];

        g_xxxyz_yyyyz_0[i] = 4.0 * g_xxxz_yyyz_1[i] * fe_0 + g_xxxz_yyyyz_1[i] * pa_y[i];

        g_xxxyz_yyyzz_0[i] = 3.0 * g_xxxz_yyzz_1[i] * fe_0 + g_xxxz_yyyzz_1[i] * pa_y[i];

        g_xxxyz_yyzzz_0[i] = 2.0 * g_xxxz_yzzz_1[i] * fe_0 + g_xxxz_yyzzz_1[i] * pa_y[i];

        g_xxxyz_yzzzz_0[i] = g_xxxz_zzzz_1[i] * fe_0 + g_xxxz_yzzzz_1[i] * pa_y[i];

        g_xxxyz_zzzzz_0[i] = g_xxxz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : HH

    auto g_xxxzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 105);

    auto g_xxxzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 106);

    auto g_xxxzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 107);

    auto g_xxxzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 108);

    auto g_xxxzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 109);

    auto g_xxxzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 110);

    auto g_xxxzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 111);

    auto g_xxxzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 112);

    auto g_xxxzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 113);

    auto g_xxxzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 114);

    auto g_xxxzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 115);

    auto g_xxxzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 116);

    auto g_xxxzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 117);

    auto g_xxxzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 118);

    auto g_xxxzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 119);

    auto g_xxxzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 120);

    auto g_xxxzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 121);

    auto g_xxxzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 122);

    auto g_xxxzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 123);

    auto g_xxxzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 124);

    auto g_xxxzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 125);

    #pragma omp simd aligned(g_xxx_xxxxx_0, g_xxx_xxxxx_1, g_xxx_xxxxy_0, g_xxx_xxxxy_1, g_xxx_xxxyy_0, g_xxx_xxxyy_1, g_xxx_xxyyy_0, g_xxx_xxyyy_1, g_xxx_xyyyy_0, g_xxx_xyyyy_1, g_xxxz_xxxxx_1, g_xxxz_xxxxy_1, g_xxxz_xxxyy_1, g_xxxz_xxyyy_1, g_xxxz_xyyyy_1, g_xxxzz_xxxxx_0, g_xxxzz_xxxxy_0, g_xxxzz_xxxxz_0, g_xxxzz_xxxyy_0, g_xxxzz_xxxyz_0, g_xxxzz_xxxzz_0, g_xxxzz_xxyyy_0, g_xxxzz_xxyyz_0, g_xxxzz_xxyzz_0, g_xxxzz_xxzzz_0, g_xxxzz_xyyyy_0, g_xxxzz_xyyyz_0, g_xxxzz_xyyzz_0, g_xxxzz_xyzzz_0, g_xxxzz_xzzzz_0, g_xxxzz_yyyyy_0, g_xxxzz_yyyyz_0, g_xxxzz_yyyzz_0, g_xxxzz_yyzzz_0, g_xxxzz_yzzzz_0, g_xxxzz_zzzzz_0, g_xxzz_xxxxz_1, g_xxzz_xxxyz_1, g_xxzz_xxxz_1, g_xxzz_xxxzz_1, g_xxzz_xxyyz_1, g_xxzz_xxyz_1, g_xxzz_xxyzz_1, g_xxzz_xxzz_1, g_xxzz_xxzzz_1, g_xxzz_xyyyz_1, g_xxzz_xyyz_1, g_xxzz_xyyzz_1, g_xxzz_xyzz_1, g_xxzz_xyzzz_1, g_xxzz_xzzz_1, g_xxzz_xzzzz_1, g_xxzz_yyyyy_1, g_xxzz_yyyyz_1, g_xxzz_yyyz_1, g_xxzz_yyyzz_1, g_xxzz_yyzz_1, g_xxzz_yyzzz_1, g_xxzz_yzzz_1, g_xxzz_yzzzz_1, g_xxzz_zzzz_1, g_xxzz_zzzzz_1, g_xzz_xxxxz_0, g_xzz_xxxxz_1, g_xzz_xxxyz_0, g_xzz_xxxyz_1, g_xzz_xxxzz_0, g_xzz_xxxzz_1, g_xzz_xxyyz_0, g_xzz_xxyyz_1, g_xzz_xxyzz_0, g_xzz_xxyzz_1, g_xzz_xxzzz_0, g_xzz_xxzzz_1, g_xzz_xyyyz_0, g_xzz_xyyyz_1, g_xzz_xyyzz_0, g_xzz_xyyzz_1, g_xzz_xyzzz_0, g_xzz_xyzzz_1, g_xzz_xzzzz_0, g_xzz_xzzzz_1, g_xzz_yyyyy_0, g_xzz_yyyyy_1, g_xzz_yyyyz_0, g_xzz_yyyyz_1, g_xzz_yyyzz_0, g_xzz_yyyzz_1, g_xzz_yyzzz_0, g_xzz_yyzzz_1, g_xzz_yzzzz_0, g_xzz_yzzzz_1, g_xzz_zzzzz_0, g_xzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzz_xxxxx_0[i] = g_xxx_xxxxx_0[i] * fbe_0 - g_xxx_xxxxx_1[i] * fz_be_0 + g_xxxz_xxxxx_1[i] * pa_z[i];

        g_xxxzz_xxxxy_0[i] = g_xxx_xxxxy_0[i] * fbe_0 - g_xxx_xxxxy_1[i] * fz_be_0 + g_xxxz_xxxxy_1[i] * pa_z[i];

        g_xxxzz_xxxxz_0[i] = 2.0 * g_xzz_xxxxz_0[i] * fbe_0 - 2.0 * g_xzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxzz_xxxz_1[i] * fe_0 + g_xxzz_xxxxz_1[i] * pa_x[i];

        g_xxxzz_xxxyy_0[i] = g_xxx_xxxyy_0[i] * fbe_0 - g_xxx_xxxyy_1[i] * fz_be_0 + g_xxxz_xxxyy_1[i] * pa_z[i];

        g_xxxzz_xxxyz_0[i] = 2.0 * g_xzz_xxxyz_0[i] * fbe_0 - 2.0 * g_xzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxzz_xxyz_1[i] * fe_0 + g_xxzz_xxxyz_1[i] * pa_x[i];

        g_xxxzz_xxxzz_0[i] = 2.0 * g_xzz_xxxzz_0[i] * fbe_0 - 2.0 * g_xzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxzz_xxzz_1[i] * fe_0 + g_xxzz_xxxzz_1[i] * pa_x[i];

        g_xxxzz_xxyyy_0[i] = g_xxx_xxyyy_0[i] * fbe_0 - g_xxx_xxyyy_1[i] * fz_be_0 + g_xxxz_xxyyy_1[i] * pa_z[i];

        g_xxxzz_xxyyz_0[i] = 2.0 * g_xzz_xxyyz_0[i] * fbe_0 - 2.0 * g_xzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxzz_xyyz_1[i] * fe_0 + g_xxzz_xxyyz_1[i] * pa_x[i];

        g_xxxzz_xxyzz_0[i] = 2.0 * g_xzz_xxyzz_0[i] * fbe_0 - 2.0 * g_xzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxzz_xyzz_1[i] * fe_0 + g_xxzz_xxyzz_1[i] * pa_x[i];

        g_xxxzz_xxzzz_0[i] = 2.0 * g_xzz_xxzzz_0[i] * fbe_0 - 2.0 * g_xzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_xzzz_1[i] * fe_0 + g_xxzz_xxzzz_1[i] * pa_x[i];

        g_xxxzz_xyyyy_0[i] = g_xxx_xyyyy_0[i] * fbe_0 - g_xxx_xyyyy_1[i] * fz_be_0 + g_xxxz_xyyyy_1[i] * pa_z[i];

        g_xxxzz_xyyyz_0[i] = 2.0 * g_xzz_xyyyz_0[i] * fbe_0 - 2.0 * g_xzz_xyyyz_1[i] * fz_be_0 + g_xxzz_yyyz_1[i] * fe_0 + g_xxzz_xyyyz_1[i] * pa_x[i];

        g_xxxzz_xyyzz_0[i] = 2.0 * g_xzz_xyyzz_0[i] * fbe_0 - 2.0 * g_xzz_xyyzz_1[i] * fz_be_0 + g_xxzz_yyzz_1[i] * fe_0 + g_xxzz_xyyzz_1[i] * pa_x[i];

        g_xxxzz_xyzzz_0[i] = 2.0 * g_xzz_xyzzz_0[i] * fbe_0 - 2.0 * g_xzz_xyzzz_1[i] * fz_be_0 + g_xxzz_yzzz_1[i] * fe_0 + g_xxzz_xyzzz_1[i] * pa_x[i];

        g_xxxzz_xzzzz_0[i] = 2.0 * g_xzz_xzzzz_0[i] * fbe_0 - 2.0 * g_xzz_xzzzz_1[i] * fz_be_0 + g_xxzz_zzzz_1[i] * fe_0 + g_xxzz_xzzzz_1[i] * pa_x[i];

        g_xxxzz_yyyyy_0[i] = 2.0 * g_xzz_yyyyy_0[i] * fbe_0 - 2.0 * g_xzz_yyyyy_1[i] * fz_be_0 + g_xxzz_yyyyy_1[i] * pa_x[i];

        g_xxxzz_yyyyz_0[i] = 2.0 * g_xzz_yyyyz_0[i] * fbe_0 - 2.0 * g_xzz_yyyyz_1[i] * fz_be_0 + g_xxzz_yyyyz_1[i] * pa_x[i];

        g_xxxzz_yyyzz_0[i] = 2.0 * g_xzz_yyyzz_0[i] * fbe_0 - 2.0 * g_xzz_yyyzz_1[i] * fz_be_0 + g_xxzz_yyyzz_1[i] * pa_x[i];

        g_xxxzz_yyzzz_0[i] = 2.0 * g_xzz_yyzzz_0[i] * fbe_0 - 2.0 * g_xzz_yyzzz_1[i] * fz_be_0 + g_xxzz_yyzzz_1[i] * pa_x[i];

        g_xxxzz_yzzzz_0[i] = 2.0 * g_xzz_yzzzz_0[i] * fbe_0 - 2.0 * g_xzz_yzzzz_1[i] * fz_be_0 + g_xxzz_yzzzz_1[i] * pa_x[i];

        g_xxxzz_zzzzz_0[i] = 2.0 * g_xzz_zzzzz_0[i] * fbe_0 - 2.0 * g_xzz_zzzzz_1[i] * fz_be_0 + g_xxzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : HH

    auto g_xxyyy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 126);

    auto g_xxyyy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 127);

    auto g_xxyyy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 128);

    auto g_xxyyy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 129);

    auto g_xxyyy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 130);

    auto g_xxyyy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 131);

    auto g_xxyyy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 132);

    auto g_xxyyy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 133);

    auto g_xxyyy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 134);

    auto g_xxyyy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 135);

    auto g_xxyyy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 136);

    auto g_xxyyy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 137);

    auto g_xxyyy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 138);

    auto g_xxyyy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 139);

    auto g_xxyyy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 140);

    auto g_xxyyy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 141);

    auto g_xxyyy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 142);

    auto g_xxyyy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 143);

    auto g_xxyyy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 144);

    auto g_xxyyy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 145);

    auto g_xxyyy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 146);

    #pragma omp simd aligned(g_xxy_xxxxx_0, g_xxy_xxxxx_1, g_xxy_xxxxz_0, g_xxy_xxxxz_1, g_xxy_xxxzz_0, g_xxy_xxxzz_1, g_xxy_xxzzz_0, g_xxy_xxzzz_1, g_xxy_xzzzz_0, g_xxy_xzzzz_1, g_xxyy_xxxxx_1, g_xxyy_xxxxz_1, g_xxyy_xxxzz_1, g_xxyy_xxzzz_1, g_xxyy_xzzzz_1, g_xxyyy_xxxxx_0, g_xxyyy_xxxxy_0, g_xxyyy_xxxxz_0, g_xxyyy_xxxyy_0, g_xxyyy_xxxyz_0, g_xxyyy_xxxzz_0, g_xxyyy_xxyyy_0, g_xxyyy_xxyyz_0, g_xxyyy_xxyzz_0, g_xxyyy_xxzzz_0, g_xxyyy_xyyyy_0, g_xxyyy_xyyyz_0, g_xxyyy_xyyzz_0, g_xxyyy_xyzzz_0, g_xxyyy_xzzzz_0, g_xxyyy_yyyyy_0, g_xxyyy_yyyyz_0, g_xxyyy_yyyzz_0, g_xxyyy_yyzzz_0, g_xxyyy_yzzzz_0, g_xxyyy_zzzzz_0, g_xyyy_xxxxy_1, g_xyyy_xxxy_1, g_xyyy_xxxyy_1, g_xyyy_xxxyz_1, g_xyyy_xxyy_1, g_xyyy_xxyyy_1, g_xyyy_xxyyz_1, g_xyyy_xxyz_1, g_xyyy_xxyzz_1, g_xyyy_xyyy_1, g_xyyy_xyyyy_1, g_xyyy_xyyyz_1, g_xyyy_xyyz_1, g_xyyy_xyyzz_1, g_xyyy_xyzz_1, g_xyyy_xyzzz_1, g_xyyy_yyyy_1, g_xyyy_yyyyy_1, g_xyyy_yyyyz_1, g_xyyy_yyyz_1, g_xyyy_yyyzz_1, g_xyyy_yyzz_1, g_xyyy_yyzzz_1, g_xyyy_yzzz_1, g_xyyy_yzzzz_1, g_xyyy_zzzzz_1, g_yyy_xxxxy_0, g_yyy_xxxxy_1, g_yyy_xxxyy_0, g_yyy_xxxyy_1, g_yyy_xxxyz_0, g_yyy_xxxyz_1, g_yyy_xxyyy_0, g_yyy_xxyyy_1, g_yyy_xxyyz_0, g_yyy_xxyyz_1, g_yyy_xxyzz_0, g_yyy_xxyzz_1, g_yyy_xyyyy_0, g_yyy_xyyyy_1, g_yyy_xyyyz_0, g_yyy_xyyyz_1, g_yyy_xyyzz_0, g_yyy_xyyzz_1, g_yyy_xyzzz_0, g_yyy_xyzzz_1, g_yyy_yyyyy_0, g_yyy_yyyyy_1, g_yyy_yyyyz_0, g_yyy_yyyyz_1, g_yyy_yyyzz_0, g_yyy_yyyzz_1, g_yyy_yyzzz_0, g_yyy_yyzzz_1, g_yyy_yzzzz_0, g_yyy_yzzzz_1, g_yyy_zzzzz_0, g_yyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyy_xxxxx_0[i] = 2.0 * g_xxy_xxxxx_0[i] * fbe_0 - 2.0 * g_xxy_xxxxx_1[i] * fz_be_0 + g_xxyy_xxxxx_1[i] * pa_y[i];

        g_xxyyy_xxxxy_0[i] = g_yyy_xxxxy_0[i] * fbe_0 - g_yyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xyyy_xxxy_1[i] * fe_0 + g_xyyy_xxxxy_1[i] * pa_x[i];

        g_xxyyy_xxxxz_0[i] = 2.0 * g_xxy_xxxxz_0[i] * fbe_0 - 2.0 * g_xxy_xxxxz_1[i] * fz_be_0 + g_xxyy_xxxxz_1[i] * pa_y[i];

        g_xxyyy_xxxyy_0[i] = g_yyy_xxxyy_0[i] * fbe_0 - g_yyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xyyy_xxyy_1[i] * fe_0 + g_xyyy_xxxyy_1[i] * pa_x[i];

        g_xxyyy_xxxyz_0[i] = g_yyy_xxxyz_0[i] * fbe_0 - g_yyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyy_xxyz_1[i] * fe_0 + g_xyyy_xxxyz_1[i] * pa_x[i];

        g_xxyyy_xxxzz_0[i] = 2.0 * g_xxy_xxxzz_0[i] * fbe_0 - 2.0 * g_xxy_xxxzz_1[i] * fz_be_0 + g_xxyy_xxxzz_1[i] * pa_y[i];

        g_xxyyy_xxyyy_0[i] = g_yyy_xxyyy_0[i] * fbe_0 - g_yyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xyyy_xyyy_1[i] * fe_0 + g_xyyy_xxyyy_1[i] * pa_x[i];

        g_xxyyy_xxyyz_0[i] = g_yyy_xxyyz_0[i] * fbe_0 - g_yyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyy_xyyz_1[i] * fe_0 + g_xyyy_xxyyz_1[i] * pa_x[i];

        g_xxyyy_xxyzz_0[i] = g_yyy_xxyzz_0[i] * fbe_0 - g_yyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyy_xyzz_1[i] * fe_0 + g_xyyy_xxyzz_1[i] * pa_x[i];

        g_xxyyy_xxzzz_0[i] = 2.0 * g_xxy_xxzzz_0[i] * fbe_0 - 2.0 * g_xxy_xxzzz_1[i] * fz_be_0 + g_xxyy_xxzzz_1[i] * pa_y[i];

        g_xxyyy_xyyyy_0[i] = g_yyy_xyyyy_0[i] * fbe_0 - g_yyy_xyyyy_1[i] * fz_be_0 + g_xyyy_yyyy_1[i] * fe_0 + g_xyyy_xyyyy_1[i] * pa_x[i];

        g_xxyyy_xyyyz_0[i] = g_yyy_xyyyz_0[i] * fbe_0 - g_yyy_xyyyz_1[i] * fz_be_0 + g_xyyy_yyyz_1[i] * fe_0 + g_xyyy_xyyyz_1[i] * pa_x[i];

        g_xxyyy_xyyzz_0[i] = g_yyy_xyyzz_0[i] * fbe_0 - g_yyy_xyyzz_1[i] * fz_be_0 + g_xyyy_yyzz_1[i] * fe_0 + g_xyyy_xyyzz_1[i] * pa_x[i];

        g_xxyyy_xyzzz_0[i] = g_yyy_xyzzz_0[i] * fbe_0 - g_yyy_xyzzz_1[i] * fz_be_0 + g_xyyy_yzzz_1[i] * fe_0 + g_xyyy_xyzzz_1[i] * pa_x[i];

        g_xxyyy_xzzzz_0[i] = 2.0 * g_xxy_xzzzz_0[i] * fbe_0 - 2.0 * g_xxy_xzzzz_1[i] * fz_be_0 + g_xxyy_xzzzz_1[i] * pa_y[i];

        g_xxyyy_yyyyy_0[i] = g_yyy_yyyyy_0[i] * fbe_0 - g_yyy_yyyyy_1[i] * fz_be_0 + g_xyyy_yyyyy_1[i] * pa_x[i];

        g_xxyyy_yyyyz_0[i] = g_yyy_yyyyz_0[i] * fbe_0 - g_yyy_yyyyz_1[i] * fz_be_0 + g_xyyy_yyyyz_1[i] * pa_x[i];

        g_xxyyy_yyyzz_0[i] = g_yyy_yyyzz_0[i] * fbe_0 - g_yyy_yyyzz_1[i] * fz_be_0 + g_xyyy_yyyzz_1[i] * pa_x[i];

        g_xxyyy_yyzzz_0[i] = g_yyy_yyzzz_0[i] * fbe_0 - g_yyy_yyzzz_1[i] * fz_be_0 + g_xyyy_yyzzz_1[i] * pa_x[i];

        g_xxyyy_yzzzz_0[i] = g_yyy_yzzzz_0[i] * fbe_0 - g_yyy_yzzzz_1[i] * fz_be_0 + g_xyyy_yzzzz_1[i] * pa_x[i];

        g_xxyyy_zzzzz_0[i] = g_yyy_zzzzz_0[i] * fbe_0 - g_yyy_zzzzz_1[i] * fz_be_0 + g_xyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : HH

    auto g_xxyyz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 147);

    auto g_xxyyz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 148);

    auto g_xxyyz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 149);

    auto g_xxyyz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 150);

    auto g_xxyyz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 151);

    auto g_xxyyz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 152);

    auto g_xxyyz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 153);

    auto g_xxyyz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 154);

    auto g_xxyyz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 155);

    auto g_xxyyz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 156);

    auto g_xxyyz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 157);

    auto g_xxyyz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 158);

    auto g_xxyyz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 159);

    auto g_xxyyz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 160);

    auto g_xxyyz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 161);

    auto g_xxyyz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 162);

    auto g_xxyyz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 163);

    auto g_xxyyz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 164);

    auto g_xxyyz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 165);

    auto g_xxyyz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 166);

    auto g_xxyyz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 167);

    #pragma omp simd aligned(g_xxyy_xxxx_1, g_xxyy_xxxxx_1, g_xxyy_xxxxy_1, g_xxyy_xxxxz_1, g_xxyy_xxxy_1, g_xxyy_xxxyy_1, g_xxyy_xxxyz_1, g_xxyy_xxxz_1, g_xxyy_xxxzz_1, g_xxyy_xxyy_1, g_xxyy_xxyyy_1, g_xxyy_xxyyz_1, g_xxyy_xxyz_1, g_xxyy_xxyzz_1, g_xxyy_xxzz_1, g_xxyy_xxzzz_1, g_xxyy_xyyy_1, g_xxyy_xyyyy_1, g_xxyy_xyyyz_1, g_xxyy_xyyz_1, g_xxyy_xyyzz_1, g_xxyy_xyzz_1, g_xxyy_xyzzz_1, g_xxyy_xzzz_1, g_xxyy_xzzzz_1, g_xxyy_yyyy_1, g_xxyy_yyyyy_1, g_xxyy_yyyyz_1, g_xxyy_yyyz_1, g_xxyy_yyyzz_1, g_xxyy_yyzz_1, g_xxyy_yyzzz_1, g_xxyy_yzzz_1, g_xxyy_yzzzz_1, g_xxyy_zzzz_1, g_xxyy_zzzzz_1, g_xxyyz_xxxxx_0, g_xxyyz_xxxxy_0, g_xxyyz_xxxxz_0, g_xxyyz_xxxyy_0, g_xxyyz_xxxyz_0, g_xxyyz_xxxzz_0, g_xxyyz_xxyyy_0, g_xxyyz_xxyyz_0, g_xxyyz_xxyzz_0, g_xxyyz_xxzzz_0, g_xxyyz_xyyyy_0, g_xxyyz_xyyyz_0, g_xxyyz_xyyzz_0, g_xxyyz_xyzzz_0, g_xxyyz_xzzzz_0, g_xxyyz_yyyyy_0, g_xxyyz_yyyyz_0, g_xxyyz_yyyzz_0, g_xxyyz_yyzzz_0, g_xxyyz_yzzzz_0, g_xxyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyz_xxxxx_0[i] = g_xxyy_xxxxx_1[i] * pa_z[i];

        g_xxyyz_xxxxy_0[i] = g_xxyy_xxxxy_1[i] * pa_z[i];

        g_xxyyz_xxxxz_0[i] = g_xxyy_xxxx_1[i] * fe_0 + g_xxyy_xxxxz_1[i] * pa_z[i];

        g_xxyyz_xxxyy_0[i] = g_xxyy_xxxyy_1[i] * pa_z[i];

        g_xxyyz_xxxyz_0[i] = g_xxyy_xxxy_1[i] * fe_0 + g_xxyy_xxxyz_1[i] * pa_z[i];

        g_xxyyz_xxxzz_0[i] = 2.0 * g_xxyy_xxxz_1[i] * fe_0 + g_xxyy_xxxzz_1[i] * pa_z[i];

        g_xxyyz_xxyyy_0[i] = g_xxyy_xxyyy_1[i] * pa_z[i];

        g_xxyyz_xxyyz_0[i] = g_xxyy_xxyy_1[i] * fe_0 + g_xxyy_xxyyz_1[i] * pa_z[i];

        g_xxyyz_xxyzz_0[i] = 2.0 * g_xxyy_xxyz_1[i] * fe_0 + g_xxyy_xxyzz_1[i] * pa_z[i];

        g_xxyyz_xxzzz_0[i] = 3.0 * g_xxyy_xxzz_1[i] * fe_0 + g_xxyy_xxzzz_1[i] * pa_z[i];

        g_xxyyz_xyyyy_0[i] = g_xxyy_xyyyy_1[i] * pa_z[i];

        g_xxyyz_xyyyz_0[i] = g_xxyy_xyyy_1[i] * fe_0 + g_xxyy_xyyyz_1[i] * pa_z[i];

        g_xxyyz_xyyzz_0[i] = 2.0 * g_xxyy_xyyz_1[i] * fe_0 + g_xxyy_xyyzz_1[i] * pa_z[i];

        g_xxyyz_xyzzz_0[i] = 3.0 * g_xxyy_xyzz_1[i] * fe_0 + g_xxyy_xyzzz_1[i] * pa_z[i];

        g_xxyyz_xzzzz_0[i] = 4.0 * g_xxyy_xzzz_1[i] * fe_0 + g_xxyy_xzzzz_1[i] * pa_z[i];

        g_xxyyz_yyyyy_0[i] = g_xxyy_yyyyy_1[i] * pa_z[i];

        g_xxyyz_yyyyz_0[i] = g_xxyy_yyyy_1[i] * fe_0 + g_xxyy_yyyyz_1[i] * pa_z[i];

        g_xxyyz_yyyzz_0[i] = 2.0 * g_xxyy_yyyz_1[i] * fe_0 + g_xxyy_yyyzz_1[i] * pa_z[i];

        g_xxyyz_yyzzz_0[i] = 3.0 * g_xxyy_yyzz_1[i] * fe_0 + g_xxyy_yyzzz_1[i] * pa_z[i];

        g_xxyyz_yzzzz_0[i] = 4.0 * g_xxyy_yzzz_1[i] * fe_0 + g_xxyy_yzzzz_1[i] * pa_z[i];

        g_xxyyz_zzzzz_0[i] = 5.0 * g_xxyy_zzzz_1[i] * fe_0 + g_xxyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 168-189 components of targeted buffer : HH

    auto g_xxyzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 168);

    auto g_xxyzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 169);

    auto g_xxyzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 170);

    auto g_xxyzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 171);

    auto g_xxyzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 172);

    auto g_xxyzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 173);

    auto g_xxyzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 174);

    auto g_xxyzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 175);

    auto g_xxyzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 176);

    auto g_xxyzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 177);

    auto g_xxyzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 178);

    auto g_xxyzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 179);

    auto g_xxyzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 180);

    auto g_xxyzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 181);

    auto g_xxyzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 182);

    auto g_xxyzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 183);

    auto g_xxyzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 184);

    auto g_xxyzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 185);

    auto g_xxyzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 186);

    auto g_xxyzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 187);

    auto g_xxyzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 188);

    #pragma omp simd aligned(g_xxyzz_xxxxx_0, g_xxyzz_xxxxy_0, g_xxyzz_xxxxz_0, g_xxyzz_xxxyy_0, g_xxyzz_xxxyz_0, g_xxyzz_xxxzz_0, g_xxyzz_xxyyy_0, g_xxyzz_xxyyz_0, g_xxyzz_xxyzz_0, g_xxyzz_xxzzz_0, g_xxyzz_xyyyy_0, g_xxyzz_xyyyz_0, g_xxyzz_xyyzz_0, g_xxyzz_xyzzz_0, g_xxyzz_xzzzz_0, g_xxyzz_yyyyy_0, g_xxyzz_yyyyz_0, g_xxyzz_yyyzz_0, g_xxyzz_yyzzz_0, g_xxyzz_yzzzz_0, g_xxyzz_zzzzz_0, g_xxzz_xxxx_1, g_xxzz_xxxxx_1, g_xxzz_xxxxy_1, g_xxzz_xxxxz_1, g_xxzz_xxxy_1, g_xxzz_xxxyy_1, g_xxzz_xxxyz_1, g_xxzz_xxxz_1, g_xxzz_xxxzz_1, g_xxzz_xxyy_1, g_xxzz_xxyyy_1, g_xxzz_xxyyz_1, g_xxzz_xxyz_1, g_xxzz_xxyzz_1, g_xxzz_xxzz_1, g_xxzz_xxzzz_1, g_xxzz_xyyy_1, g_xxzz_xyyyy_1, g_xxzz_xyyyz_1, g_xxzz_xyyz_1, g_xxzz_xyyzz_1, g_xxzz_xyzz_1, g_xxzz_xyzzz_1, g_xxzz_xzzz_1, g_xxzz_xzzzz_1, g_xxzz_yyyy_1, g_xxzz_yyyyy_1, g_xxzz_yyyyz_1, g_xxzz_yyyz_1, g_xxzz_yyyzz_1, g_xxzz_yyzz_1, g_xxzz_yyzzz_1, g_xxzz_yzzz_1, g_xxzz_yzzzz_1, g_xxzz_zzzz_1, g_xxzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzz_xxxxx_0[i] = g_xxzz_xxxxx_1[i] * pa_y[i];

        g_xxyzz_xxxxy_0[i] = g_xxzz_xxxx_1[i] * fe_0 + g_xxzz_xxxxy_1[i] * pa_y[i];

        g_xxyzz_xxxxz_0[i] = g_xxzz_xxxxz_1[i] * pa_y[i];

        g_xxyzz_xxxyy_0[i] = 2.0 * g_xxzz_xxxy_1[i] * fe_0 + g_xxzz_xxxyy_1[i] * pa_y[i];

        g_xxyzz_xxxyz_0[i] = g_xxzz_xxxz_1[i] * fe_0 + g_xxzz_xxxyz_1[i] * pa_y[i];

        g_xxyzz_xxxzz_0[i] = g_xxzz_xxxzz_1[i] * pa_y[i];

        g_xxyzz_xxyyy_0[i] = 3.0 * g_xxzz_xxyy_1[i] * fe_0 + g_xxzz_xxyyy_1[i] * pa_y[i];

        g_xxyzz_xxyyz_0[i] = 2.0 * g_xxzz_xxyz_1[i] * fe_0 + g_xxzz_xxyyz_1[i] * pa_y[i];

        g_xxyzz_xxyzz_0[i] = g_xxzz_xxzz_1[i] * fe_0 + g_xxzz_xxyzz_1[i] * pa_y[i];

        g_xxyzz_xxzzz_0[i] = g_xxzz_xxzzz_1[i] * pa_y[i];

        g_xxyzz_xyyyy_0[i] = 4.0 * g_xxzz_xyyy_1[i] * fe_0 + g_xxzz_xyyyy_1[i] * pa_y[i];

        g_xxyzz_xyyyz_0[i] = 3.0 * g_xxzz_xyyz_1[i] * fe_0 + g_xxzz_xyyyz_1[i] * pa_y[i];

        g_xxyzz_xyyzz_0[i] = 2.0 * g_xxzz_xyzz_1[i] * fe_0 + g_xxzz_xyyzz_1[i] * pa_y[i];

        g_xxyzz_xyzzz_0[i] = g_xxzz_xzzz_1[i] * fe_0 + g_xxzz_xyzzz_1[i] * pa_y[i];

        g_xxyzz_xzzzz_0[i] = g_xxzz_xzzzz_1[i] * pa_y[i];

        g_xxyzz_yyyyy_0[i] = 5.0 * g_xxzz_yyyy_1[i] * fe_0 + g_xxzz_yyyyy_1[i] * pa_y[i];

        g_xxyzz_yyyyz_0[i] = 4.0 * g_xxzz_yyyz_1[i] * fe_0 + g_xxzz_yyyyz_1[i] * pa_y[i];

        g_xxyzz_yyyzz_0[i] = 3.0 * g_xxzz_yyzz_1[i] * fe_0 + g_xxzz_yyyzz_1[i] * pa_y[i];

        g_xxyzz_yyzzz_0[i] = 2.0 * g_xxzz_yzzz_1[i] * fe_0 + g_xxzz_yyzzz_1[i] * pa_y[i];

        g_xxyzz_yzzzz_0[i] = g_xxzz_zzzz_1[i] * fe_0 + g_xxzz_yzzzz_1[i] * pa_y[i];

        g_xxyzz_zzzzz_0[i] = g_xxzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : HH

    auto g_xxzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 189);

    auto g_xxzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 190);

    auto g_xxzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 191);

    auto g_xxzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 192);

    auto g_xxzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 193);

    auto g_xxzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 194);

    auto g_xxzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 195);

    auto g_xxzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 196);

    auto g_xxzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 197);

    auto g_xxzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 198);

    auto g_xxzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 199);

    auto g_xxzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 200);

    auto g_xxzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 201);

    auto g_xxzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 202);

    auto g_xxzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 203);

    auto g_xxzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 204);

    auto g_xxzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 205);

    auto g_xxzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 206);

    auto g_xxzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 207);

    auto g_xxzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 208);

    auto g_xxzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 209);

    #pragma omp simd aligned(g_xxz_xxxxx_0, g_xxz_xxxxx_1, g_xxz_xxxxy_0, g_xxz_xxxxy_1, g_xxz_xxxyy_0, g_xxz_xxxyy_1, g_xxz_xxyyy_0, g_xxz_xxyyy_1, g_xxz_xyyyy_0, g_xxz_xyyyy_1, g_xxzz_xxxxx_1, g_xxzz_xxxxy_1, g_xxzz_xxxyy_1, g_xxzz_xxyyy_1, g_xxzz_xyyyy_1, g_xxzzz_xxxxx_0, g_xxzzz_xxxxy_0, g_xxzzz_xxxxz_0, g_xxzzz_xxxyy_0, g_xxzzz_xxxyz_0, g_xxzzz_xxxzz_0, g_xxzzz_xxyyy_0, g_xxzzz_xxyyz_0, g_xxzzz_xxyzz_0, g_xxzzz_xxzzz_0, g_xxzzz_xyyyy_0, g_xxzzz_xyyyz_0, g_xxzzz_xyyzz_0, g_xxzzz_xyzzz_0, g_xxzzz_xzzzz_0, g_xxzzz_yyyyy_0, g_xxzzz_yyyyz_0, g_xxzzz_yyyzz_0, g_xxzzz_yyzzz_0, g_xxzzz_yzzzz_0, g_xxzzz_zzzzz_0, g_xzzz_xxxxz_1, g_xzzz_xxxyz_1, g_xzzz_xxxz_1, g_xzzz_xxxzz_1, g_xzzz_xxyyz_1, g_xzzz_xxyz_1, g_xzzz_xxyzz_1, g_xzzz_xxzz_1, g_xzzz_xxzzz_1, g_xzzz_xyyyz_1, g_xzzz_xyyz_1, g_xzzz_xyyzz_1, g_xzzz_xyzz_1, g_xzzz_xyzzz_1, g_xzzz_xzzz_1, g_xzzz_xzzzz_1, g_xzzz_yyyyy_1, g_xzzz_yyyyz_1, g_xzzz_yyyz_1, g_xzzz_yyyzz_1, g_xzzz_yyzz_1, g_xzzz_yyzzz_1, g_xzzz_yzzz_1, g_xzzz_yzzzz_1, g_xzzz_zzzz_1, g_xzzz_zzzzz_1, g_zzz_xxxxz_0, g_zzz_xxxxz_1, g_zzz_xxxyz_0, g_zzz_xxxyz_1, g_zzz_xxxzz_0, g_zzz_xxxzz_1, g_zzz_xxyyz_0, g_zzz_xxyyz_1, g_zzz_xxyzz_0, g_zzz_xxyzz_1, g_zzz_xxzzz_0, g_zzz_xxzzz_1, g_zzz_xyyyz_0, g_zzz_xyyyz_1, g_zzz_xyyzz_0, g_zzz_xyyzz_1, g_zzz_xyzzz_0, g_zzz_xyzzz_1, g_zzz_xzzzz_0, g_zzz_xzzzz_1, g_zzz_yyyyy_0, g_zzz_yyyyy_1, g_zzz_yyyyz_0, g_zzz_yyyyz_1, g_zzz_yyyzz_0, g_zzz_yyyzz_1, g_zzz_yyzzz_0, g_zzz_yyzzz_1, g_zzz_yzzzz_0, g_zzz_yzzzz_1, g_zzz_zzzzz_0, g_zzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzz_xxxxx_0[i] = 2.0 * g_xxz_xxxxx_0[i] * fbe_0 - 2.0 * g_xxz_xxxxx_1[i] * fz_be_0 + g_xxzz_xxxxx_1[i] * pa_z[i];

        g_xxzzz_xxxxy_0[i] = 2.0 * g_xxz_xxxxy_0[i] * fbe_0 - 2.0 * g_xxz_xxxxy_1[i] * fz_be_0 + g_xxzz_xxxxy_1[i] * pa_z[i];

        g_xxzzz_xxxxz_0[i] = g_zzz_xxxxz_0[i] * fbe_0 - g_zzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xzzz_xxxz_1[i] * fe_0 + g_xzzz_xxxxz_1[i] * pa_x[i];

        g_xxzzz_xxxyy_0[i] = 2.0 * g_xxz_xxxyy_0[i] * fbe_0 - 2.0 * g_xxz_xxxyy_1[i] * fz_be_0 + g_xxzz_xxxyy_1[i] * pa_z[i];

        g_xxzzz_xxxyz_0[i] = g_zzz_xxxyz_0[i] * fbe_0 - g_zzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xzzz_xxyz_1[i] * fe_0 + g_xzzz_xxxyz_1[i] * pa_x[i];

        g_xxzzz_xxxzz_0[i] = g_zzz_xxxzz_0[i] * fbe_0 - g_zzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xzzz_xxzz_1[i] * fe_0 + g_xzzz_xxxzz_1[i] * pa_x[i];

        g_xxzzz_xxyyy_0[i] = 2.0 * g_xxz_xxyyy_0[i] * fbe_0 - 2.0 * g_xxz_xxyyy_1[i] * fz_be_0 + g_xxzz_xxyyy_1[i] * pa_z[i];

        g_xxzzz_xxyyz_0[i] = g_zzz_xxyyz_0[i] * fbe_0 - g_zzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xzzz_xyyz_1[i] * fe_0 + g_xzzz_xxyyz_1[i] * pa_x[i];

        g_xxzzz_xxyzz_0[i] = g_zzz_xxyzz_0[i] * fbe_0 - g_zzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xzzz_xyzz_1[i] * fe_0 + g_xzzz_xxyzz_1[i] * pa_x[i];

        g_xxzzz_xxzzz_0[i] = g_zzz_xxzzz_0[i] * fbe_0 - g_zzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_xzzz_1[i] * fe_0 + g_xzzz_xxzzz_1[i] * pa_x[i];

        g_xxzzz_xyyyy_0[i] = 2.0 * g_xxz_xyyyy_0[i] * fbe_0 - 2.0 * g_xxz_xyyyy_1[i] * fz_be_0 + g_xxzz_xyyyy_1[i] * pa_z[i];

        g_xxzzz_xyyyz_0[i] = g_zzz_xyyyz_0[i] * fbe_0 - g_zzz_xyyyz_1[i] * fz_be_0 + g_xzzz_yyyz_1[i] * fe_0 + g_xzzz_xyyyz_1[i] * pa_x[i];

        g_xxzzz_xyyzz_0[i] = g_zzz_xyyzz_0[i] * fbe_0 - g_zzz_xyyzz_1[i] * fz_be_0 + g_xzzz_yyzz_1[i] * fe_0 + g_xzzz_xyyzz_1[i] * pa_x[i];

        g_xxzzz_xyzzz_0[i] = g_zzz_xyzzz_0[i] * fbe_0 - g_zzz_xyzzz_1[i] * fz_be_0 + g_xzzz_yzzz_1[i] * fe_0 + g_xzzz_xyzzz_1[i] * pa_x[i];

        g_xxzzz_xzzzz_0[i] = g_zzz_xzzzz_0[i] * fbe_0 - g_zzz_xzzzz_1[i] * fz_be_0 + g_xzzz_zzzz_1[i] * fe_0 + g_xzzz_xzzzz_1[i] * pa_x[i];

        g_xxzzz_yyyyy_0[i] = g_zzz_yyyyy_0[i] * fbe_0 - g_zzz_yyyyy_1[i] * fz_be_0 + g_xzzz_yyyyy_1[i] * pa_x[i];

        g_xxzzz_yyyyz_0[i] = g_zzz_yyyyz_0[i] * fbe_0 - g_zzz_yyyyz_1[i] * fz_be_0 + g_xzzz_yyyyz_1[i] * pa_x[i];

        g_xxzzz_yyyzz_0[i] = g_zzz_yyyzz_0[i] * fbe_0 - g_zzz_yyyzz_1[i] * fz_be_0 + g_xzzz_yyyzz_1[i] * pa_x[i];

        g_xxzzz_yyzzz_0[i] = g_zzz_yyzzz_0[i] * fbe_0 - g_zzz_yyzzz_1[i] * fz_be_0 + g_xzzz_yyzzz_1[i] * pa_x[i];

        g_xxzzz_yzzzz_0[i] = g_zzz_yzzzz_0[i] * fbe_0 - g_zzz_yzzzz_1[i] * fz_be_0 + g_xzzz_yzzzz_1[i] * pa_x[i];

        g_xxzzz_zzzzz_0[i] = g_zzz_zzzzz_0[i] * fbe_0 - g_zzz_zzzzz_1[i] * fz_be_0 + g_xzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : HH

    auto g_xyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 210);

    auto g_xyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 211);

    auto g_xyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 212);

    auto g_xyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 213);

    auto g_xyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 214);

    auto g_xyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 215);

    auto g_xyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 216);

    auto g_xyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 217);

    auto g_xyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 218);

    auto g_xyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 219);

    auto g_xyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 220);

    auto g_xyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 221);

    auto g_xyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 222);

    auto g_xyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 223);

    auto g_xyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 224);

    auto g_xyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 225);

    auto g_xyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 226);

    auto g_xyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 227);

    auto g_xyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 228);

    auto g_xyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 229);

    auto g_xyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 230);

    #pragma omp simd aligned(g_xyyyy_xxxxx_0, g_xyyyy_xxxxy_0, g_xyyyy_xxxxz_0, g_xyyyy_xxxyy_0, g_xyyyy_xxxyz_0, g_xyyyy_xxxzz_0, g_xyyyy_xxyyy_0, g_xyyyy_xxyyz_0, g_xyyyy_xxyzz_0, g_xyyyy_xxzzz_0, g_xyyyy_xyyyy_0, g_xyyyy_xyyyz_0, g_xyyyy_xyyzz_0, g_xyyyy_xyzzz_0, g_xyyyy_xzzzz_0, g_xyyyy_yyyyy_0, g_xyyyy_yyyyz_0, g_xyyyy_yyyzz_0, g_xyyyy_yyzzz_0, g_xyyyy_yzzzz_0, g_xyyyy_zzzzz_0, g_yyyy_xxxx_1, g_yyyy_xxxxx_1, g_yyyy_xxxxy_1, g_yyyy_xxxxz_1, g_yyyy_xxxy_1, g_yyyy_xxxyy_1, g_yyyy_xxxyz_1, g_yyyy_xxxz_1, g_yyyy_xxxzz_1, g_yyyy_xxyy_1, g_yyyy_xxyyy_1, g_yyyy_xxyyz_1, g_yyyy_xxyz_1, g_yyyy_xxyzz_1, g_yyyy_xxzz_1, g_yyyy_xxzzz_1, g_yyyy_xyyy_1, g_yyyy_xyyyy_1, g_yyyy_xyyyz_1, g_yyyy_xyyz_1, g_yyyy_xyyzz_1, g_yyyy_xyzz_1, g_yyyy_xyzzz_1, g_yyyy_xzzz_1, g_yyyy_xzzzz_1, g_yyyy_yyyy_1, g_yyyy_yyyyy_1, g_yyyy_yyyyz_1, g_yyyy_yyyz_1, g_yyyy_yyyzz_1, g_yyyy_yyzz_1, g_yyyy_yyzzz_1, g_yyyy_yzzz_1, g_yyyy_yzzzz_1, g_yyyy_zzzz_1, g_yyyy_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyy_xxxxx_0[i] = 5.0 * g_yyyy_xxxx_1[i] * fe_0 + g_yyyy_xxxxx_1[i] * pa_x[i];

        g_xyyyy_xxxxy_0[i] = 4.0 * g_yyyy_xxxy_1[i] * fe_0 + g_yyyy_xxxxy_1[i] * pa_x[i];

        g_xyyyy_xxxxz_0[i] = 4.0 * g_yyyy_xxxz_1[i] * fe_0 + g_yyyy_xxxxz_1[i] * pa_x[i];

        g_xyyyy_xxxyy_0[i] = 3.0 * g_yyyy_xxyy_1[i] * fe_0 + g_yyyy_xxxyy_1[i] * pa_x[i];

        g_xyyyy_xxxyz_0[i] = 3.0 * g_yyyy_xxyz_1[i] * fe_0 + g_yyyy_xxxyz_1[i] * pa_x[i];

        g_xyyyy_xxxzz_0[i] = 3.0 * g_yyyy_xxzz_1[i] * fe_0 + g_yyyy_xxxzz_1[i] * pa_x[i];

        g_xyyyy_xxyyy_0[i] = 2.0 * g_yyyy_xyyy_1[i] * fe_0 + g_yyyy_xxyyy_1[i] * pa_x[i];

        g_xyyyy_xxyyz_0[i] = 2.0 * g_yyyy_xyyz_1[i] * fe_0 + g_yyyy_xxyyz_1[i] * pa_x[i];

        g_xyyyy_xxyzz_0[i] = 2.0 * g_yyyy_xyzz_1[i] * fe_0 + g_yyyy_xxyzz_1[i] * pa_x[i];

        g_xyyyy_xxzzz_0[i] = 2.0 * g_yyyy_xzzz_1[i] * fe_0 + g_yyyy_xxzzz_1[i] * pa_x[i];

        g_xyyyy_xyyyy_0[i] = g_yyyy_yyyy_1[i] * fe_0 + g_yyyy_xyyyy_1[i] * pa_x[i];

        g_xyyyy_xyyyz_0[i] = g_yyyy_yyyz_1[i] * fe_0 + g_yyyy_xyyyz_1[i] * pa_x[i];

        g_xyyyy_xyyzz_0[i] = g_yyyy_yyzz_1[i] * fe_0 + g_yyyy_xyyzz_1[i] * pa_x[i];

        g_xyyyy_xyzzz_0[i] = g_yyyy_yzzz_1[i] * fe_0 + g_yyyy_xyzzz_1[i] * pa_x[i];

        g_xyyyy_xzzzz_0[i] = g_yyyy_zzzz_1[i] * fe_0 + g_yyyy_xzzzz_1[i] * pa_x[i];

        g_xyyyy_yyyyy_0[i] = g_yyyy_yyyyy_1[i] * pa_x[i];

        g_xyyyy_yyyyz_0[i] = g_yyyy_yyyyz_1[i] * pa_x[i];

        g_xyyyy_yyyzz_0[i] = g_yyyy_yyyzz_1[i] * pa_x[i];

        g_xyyyy_yyzzz_0[i] = g_yyyy_yyzzz_1[i] * pa_x[i];

        g_xyyyy_yzzzz_0[i] = g_yyyy_yzzzz_1[i] * pa_x[i];

        g_xyyyy_zzzzz_0[i] = g_yyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 231-252 components of targeted buffer : HH

    auto g_xyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 231);

    auto g_xyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 232);

    auto g_xyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 233);

    auto g_xyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 234);

    auto g_xyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 235);

    auto g_xyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 236);

    auto g_xyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 237);

    auto g_xyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 238);

    auto g_xyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 239);

    auto g_xyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 240);

    auto g_xyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 241);

    auto g_xyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 242);

    auto g_xyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 243);

    auto g_xyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 244);

    auto g_xyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 245);

    auto g_xyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 246);

    auto g_xyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 247);

    auto g_xyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 248);

    auto g_xyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 249);

    auto g_xyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 250);

    auto g_xyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 251);

    #pragma omp simd aligned(g_xyyy_xxxxx_1, g_xyyy_xxxxy_1, g_xyyy_xxxyy_1, g_xyyy_xxyyy_1, g_xyyy_xyyyy_1, g_xyyyz_xxxxx_0, g_xyyyz_xxxxy_0, g_xyyyz_xxxxz_0, g_xyyyz_xxxyy_0, g_xyyyz_xxxyz_0, g_xyyyz_xxxzz_0, g_xyyyz_xxyyy_0, g_xyyyz_xxyyz_0, g_xyyyz_xxyzz_0, g_xyyyz_xxzzz_0, g_xyyyz_xyyyy_0, g_xyyyz_xyyyz_0, g_xyyyz_xyyzz_0, g_xyyyz_xyzzz_0, g_xyyyz_xzzzz_0, g_xyyyz_yyyyy_0, g_xyyyz_yyyyz_0, g_xyyyz_yyyzz_0, g_xyyyz_yyzzz_0, g_xyyyz_yzzzz_0, g_xyyyz_zzzzz_0, g_yyyz_xxxxz_1, g_yyyz_xxxyz_1, g_yyyz_xxxz_1, g_yyyz_xxxzz_1, g_yyyz_xxyyz_1, g_yyyz_xxyz_1, g_yyyz_xxyzz_1, g_yyyz_xxzz_1, g_yyyz_xxzzz_1, g_yyyz_xyyyz_1, g_yyyz_xyyz_1, g_yyyz_xyyzz_1, g_yyyz_xyzz_1, g_yyyz_xyzzz_1, g_yyyz_xzzz_1, g_yyyz_xzzzz_1, g_yyyz_yyyyy_1, g_yyyz_yyyyz_1, g_yyyz_yyyz_1, g_yyyz_yyyzz_1, g_yyyz_yyzz_1, g_yyyz_yyzzz_1, g_yyyz_yzzz_1, g_yyyz_yzzzz_1, g_yyyz_zzzz_1, g_yyyz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyz_xxxxx_0[i] = g_xyyy_xxxxx_1[i] * pa_z[i];

        g_xyyyz_xxxxy_0[i] = g_xyyy_xxxxy_1[i] * pa_z[i];

        g_xyyyz_xxxxz_0[i] = 4.0 * g_yyyz_xxxz_1[i] * fe_0 + g_yyyz_xxxxz_1[i] * pa_x[i];

        g_xyyyz_xxxyy_0[i] = g_xyyy_xxxyy_1[i] * pa_z[i];

        g_xyyyz_xxxyz_0[i] = 3.0 * g_yyyz_xxyz_1[i] * fe_0 + g_yyyz_xxxyz_1[i] * pa_x[i];

        g_xyyyz_xxxzz_0[i] = 3.0 * g_yyyz_xxzz_1[i] * fe_0 + g_yyyz_xxxzz_1[i] * pa_x[i];

        g_xyyyz_xxyyy_0[i] = g_xyyy_xxyyy_1[i] * pa_z[i];

        g_xyyyz_xxyyz_0[i] = 2.0 * g_yyyz_xyyz_1[i] * fe_0 + g_yyyz_xxyyz_1[i] * pa_x[i];

        g_xyyyz_xxyzz_0[i] = 2.0 * g_yyyz_xyzz_1[i] * fe_0 + g_yyyz_xxyzz_1[i] * pa_x[i];

        g_xyyyz_xxzzz_0[i] = 2.0 * g_yyyz_xzzz_1[i] * fe_0 + g_yyyz_xxzzz_1[i] * pa_x[i];

        g_xyyyz_xyyyy_0[i] = g_xyyy_xyyyy_1[i] * pa_z[i];

        g_xyyyz_xyyyz_0[i] = g_yyyz_yyyz_1[i] * fe_0 + g_yyyz_xyyyz_1[i] * pa_x[i];

        g_xyyyz_xyyzz_0[i] = g_yyyz_yyzz_1[i] * fe_0 + g_yyyz_xyyzz_1[i] * pa_x[i];

        g_xyyyz_xyzzz_0[i] = g_yyyz_yzzz_1[i] * fe_0 + g_yyyz_xyzzz_1[i] * pa_x[i];

        g_xyyyz_xzzzz_0[i] = g_yyyz_zzzz_1[i] * fe_0 + g_yyyz_xzzzz_1[i] * pa_x[i];

        g_xyyyz_yyyyy_0[i] = g_yyyz_yyyyy_1[i] * pa_x[i];

        g_xyyyz_yyyyz_0[i] = g_yyyz_yyyyz_1[i] * pa_x[i];

        g_xyyyz_yyyzz_0[i] = g_yyyz_yyyzz_1[i] * pa_x[i];

        g_xyyyz_yyzzz_0[i] = g_yyyz_yyzzz_1[i] * pa_x[i];

        g_xyyyz_yzzzz_0[i] = g_yyyz_yzzzz_1[i] * pa_x[i];

        g_xyyyz_zzzzz_0[i] = g_yyyz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 252-273 components of targeted buffer : HH

    auto g_xyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 252);

    auto g_xyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 253);

    auto g_xyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 254);

    auto g_xyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 255);

    auto g_xyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 256);

    auto g_xyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 257);

    auto g_xyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 258);

    auto g_xyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 259);

    auto g_xyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 260);

    auto g_xyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 261);

    auto g_xyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 262);

    auto g_xyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 263);

    auto g_xyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 264);

    auto g_xyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 265);

    auto g_xyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 266);

    auto g_xyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 267);

    auto g_xyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 268);

    auto g_xyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 269);

    auto g_xyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 270);

    auto g_xyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 271);

    auto g_xyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 272);

    #pragma omp simd aligned(g_xyyzz_xxxxx_0, g_xyyzz_xxxxy_0, g_xyyzz_xxxxz_0, g_xyyzz_xxxyy_0, g_xyyzz_xxxyz_0, g_xyyzz_xxxzz_0, g_xyyzz_xxyyy_0, g_xyyzz_xxyyz_0, g_xyyzz_xxyzz_0, g_xyyzz_xxzzz_0, g_xyyzz_xyyyy_0, g_xyyzz_xyyyz_0, g_xyyzz_xyyzz_0, g_xyyzz_xyzzz_0, g_xyyzz_xzzzz_0, g_xyyzz_yyyyy_0, g_xyyzz_yyyyz_0, g_xyyzz_yyyzz_0, g_xyyzz_yyzzz_0, g_xyyzz_yzzzz_0, g_xyyzz_zzzzz_0, g_yyzz_xxxx_1, g_yyzz_xxxxx_1, g_yyzz_xxxxy_1, g_yyzz_xxxxz_1, g_yyzz_xxxy_1, g_yyzz_xxxyy_1, g_yyzz_xxxyz_1, g_yyzz_xxxz_1, g_yyzz_xxxzz_1, g_yyzz_xxyy_1, g_yyzz_xxyyy_1, g_yyzz_xxyyz_1, g_yyzz_xxyz_1, g_yyzz_xxyzz_1, g_yyzz_xxzz_1, g_yyzz_xxzzz_1, g_yyzz_xyyy_1, g_yyzz_xyyyy_1, g_yyzz_xyyyz_1, g_yyzz_xyyz_1, g_yyzz_xyyzz_1, g_yyzz_xyzz_1, g_yyzz_xyzzz_1, g_yyzz_xzzz_1, g_yyzz_xzzzz_1, g_yyzz_yyyy_1, g_yyzz_yyyyy_1, g_yyzz_yyyyz_1, g_yyzz_yyyz_1, g_yyzz_yyyzz_1, g_yyzz_yyzz_1, g_yyzz_yyzzz_1, g_yyzz_yzzz_1, g_yyzz_yzzzz_1, g_yyzz_zzzz_1, g_yyzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzz_xxxxx_0[i] = 5.0 * g_yyzz_xxxx_1[i] * fe_0 + g_yyzz_xxxxx_1[i] * pa_x[i];

        g_xyyzz_xxxxy_0[i] = 4.0 * g_yyzz_xxxy_1[i] * fe_0 + g_yyzz_xxxxy_1[i] * pa_x[i];

        g_xyyzz_xxxxz_0[i] = 4.0 * g_yyzz_xxxz_1[i] * fe_0 + g_yyzz_xxxxz_1[i] * pa_x[i];

        g_xyyzz_xxxyy_0[i] = 3.0 * g_yyzz_xxyy_1[i] * fe_0 + g_yyzz_xxxyy_1[i] * pa_x[i];

        g_xyyzz_xxxyz_0[i] = 3.0 * g_yyzz_xxyz_1[i] * fe_0 + g_yyzz_xxxyz_1[i] * pa_x[i];

        g_xyyzz_xxxzz_0[i] = 3.0 * g_yyzz_xxzz_1[i] * fe_0 + g_yyzz_xxxzz_1[i] * pa_x[i];

        g_xyyzz_xxyyy_0[i] = 2.0 * g_yyzz_xyyy_1[i] * fe_0 + g_yyzz_xxyyy_1[i] * pa_x[i];

        g_xyyzz_xxyyz_0[i] = 2.0 * g_yyzz_xyyz_1[i] * fe_0 + g_yyzz_xxyyz_1[i] * pa_x[i];

        g_xyyzz_xxyzz_0[i] = 2.0 * g_yyzz_xyzz_1[i] * fe_0 + g_yyzz_xxyzz_1[i] * pa_x[i];

        g_xyyzz_xxzzz_0[i] = 2.0 * g_yyzz_xzzz_1[i] * fe_0 + g_yyzz_xxzzz_1[i] * pa_x[i];

        g_xyyzz_xyyyy_0[i] = g_yyzz_yyyy_1[i] * fe_0 + g_yyzz_xyyyy_1[i] * pa_x[i];

        g_xyyzz_xyyyz_0[i] = g_yyzz_yyyz_1[i] * fe_0 + g_yyzz_xyyyz_1[i] * pa_x[i];

        g_xyyzz_xyyzz_0[i] = g_yyzz_yyzz_1[i] * fe_0 + g_yyzz_xyyzz_1[i] * pa_x[i];

        g_xyyzz_xyzzz_0[i] = g_yyzz_yzzz_1[i] * fe_0 + g_yyzz_xyzzz_1[i] * pa_x[i];

        g_xyyzz_xzzzz_0[i] = g_yyzz_zzzz_1[i] * fe_0 + g_yyzz_xzzzz_1[i] * pa_x[i];

        g_xyyzz_yyyyy_0[i] = g_yyzz_yyyyy_1[i] * pa_x[i];

        g_xyyzz_yyyyz_0[i] = g_yyzz_yyyyz_1[i] * pa_x[i];

        g_xyyzz_yyyzz_0[i] = g_yyzz_yyyzz_1[i] * pa_x[i];

        g_xyyzz_yyzzz_0[i] = g_yyzz_yyzzz_1[i] * pa_x[i];

        g_xyyzz_yzzzz_0[i] = g_yyzz_yzzzz_1[i] * pa_x[i];

        g_xyyzz_zzzzz_0[i] = g_yyzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 273-294 components of targeted buffer : HH

    auto g_xyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 273);

    auto g_xyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 274);

    auto g_xyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 275);

    auto g_xyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 276);

    auto g_xyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 277);

    auto g_xyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 278);

    auto g_xyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 279);

    auto g_xyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 280);

    auto g_xyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 281);

    auto g_xyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 282);

    auto g_xyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 283);

    auto g_xyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 284);

    auto g_xyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 285);

    auto g_xyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 286);

    auto g_xyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 287);

    auto g_xyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 288);

    auto g_xyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 289);

    auto g_xyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 290);

    auto g_xyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 291);

    auto g_xyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 292);

    auto g_xyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 293);

    #pragma omp simd aligned(g_xyzzz_xxxxx_0, g_xyzzz_xxxxy_0, g_xyzzz_xxxxz_0, g_xyzzz_xxxyy_0, g_xyzzz_xxxyz_0, g_xyzzz_xxxzz_0, g_xyzzz_xxyyy_0, g_xyzzz_xxyyz_0, g_xyzzz_xxyzz_0, g_xyzzz_xxzzz_0, g_xyzzz_xyyyy_0, g_xyzzz_xyyyz_0, g_xyzzz_xyyzz_0, g_xyzzz_xyzzz_0, g_xyzzz_xzzzz_0, g_xyzzz_yyyyy_0, g_xyzzz_yyyyz_0, g_xyzzz_yyyzz_0, g_xyzzz_yyzzz_0, g_xyzzz_yzzzz_0, g_xyzzz_zzzzz_0, g_xzzz_xxxxx_1, g_xzzz_xxxxz_1, g_xzzz_xxxzz_1, g_xzzz_xxzzz_1, g_xzzz_xzzzz_1, g_yzzz_xxxxy_1, g_yzzz_xxxy_1, g_yzzz_xxxyy_1, g_yzzz_xxxyz_1, g_yzzz_xxyy_1, g_yzzz_xxyyy_1, g_yzzz_xxyyz_1, g_yzzz_xxyz_1, g_yzzz_xxyzz_1, g_yzzz_xyyy_1, g_yzzz_xyyyy_1, g_yzzz_xyyyz_1, g_yzzz_xyyz_1, g_yzzz_xyyzz_1, g_yzzz_xyzz_1, g_yzzz_xyzzz_1, g_yzzz_yyyy_1, g_yzzz_yyyyy_1, g_yzzz_yyyyz_1, g_yzzz_yyyz_1, g_yzzz_yyyzz_1, g_yzzz_yyzz_1, g_yzzz_yyzzz_1, g_yzzz_yzzz_1, g_yzzz_yzzzz_1, g_yzzz_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzz_xxxxx_0[i] = g_xzzz_xxxxx_1[i] * pa_y[i];

        g_xyzzz_xxxxy_0[i] = 4.0 * g_yzzz_xxxy_1[i] * fe_0 + g_yzzz_xxxxy_1[i] * pa_x[i];

        g_xyzzz_xxxxz_0[i] = g_xzzz_xxxxz_1[i] * pa_y[i];

        g_xyzzz_xxxyy_0[i] = 3.0 * g_yzzz_xxyy_1[i] * fe_0 + g_yzzz_xxxyy_1[i] * pa_x[i];

        g_xyzzz_xxxyz_0[i] = 3.0 * g_yzzz_xxyz_1[i] * fe_0 + g_yzzz_xxxyz_1[i] * pa_x[i];

        g_xyzzz_xxxzz_0[i] = g_xzzz_xxxzz_1[i] * pa_y[i];

        g_xyzzz_xxyyy_0[i] = 2.0 * g_yzzz_xyyy_1[i] * fe_0 + g_yzzz_xxyyy_1[i] * pa_x[i];

        g_xyzzz_xxyyz_0[i] = 2.0 * g_yzzz_xyyz_1[i] * fe_0 + g_yzzz_xxyyz_1[i] * pa_x[i];

        g_xyzzz_xxyzz_0[i] = 2.0 * g_yzzz_xyzz_1[i] * fe_0 + g_yzzz_xxyzz_1[i] * pa_x[i];

        g_xyzzz_xxzzz_0[i] = g_xzzz_xxzzz_1[i] * pa_y[i];

        g_xyzzz_xyyyy_0[i] = g_yzzz_yyyy_1[i] * fe_0 + g_yzzz_xyyyy_1[i] * pa_x[i];

        g_xyzzz_xyyyz_0[i] = g_yzzz_yyyz_1[i] * fe_0 + g_yzzz_xyyyz_1[i] * pa_x[i];

        g_xyzzz_xyyzz_0[i] = g_yzzz_yyzz_1[i] * fe_0 + g_yzzz_xyyzz_1[i] * pa_x[i];

        g_xyzzz_xyzzz_0[i] = g_yzzz_yzzz_1[i] * fe_0 + g_yzzz_xyzzz_1[i] * pa_x[i];

        g_xyzzz_xzzzz_0[i] = g_xzzz_xzzzz_1[i] * pa_y[i];

        g_xyzzz_yyyyy_0[i] = g_yzzz_yyyyy_1[i] * pa_x[i];

        g_xyzzz_yyyyz_0[i] = g_yzzz_yyyyz_1[i] * pa_x[i];

        g_xyzzz_yyyzz_0[i] = g_yzzz_yyyzz_1[i] * pa_x[i];

        g_xyzzz_yyzzz_0[i] = g_yzzz_yyzzz_1[i] * pa_x[i];

        g_xyzzz_yzzzz_0[i] = g_yzzz_yzzzz_1[i] * pa_x[i];

        g_xyzzz_zzzzz_0[i] = g_yzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 294-315 components of targeted buffer : HH

    auto g_xzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 294);

    auto g_xzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 295);

    auto g_xzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 296);

    auto g_xzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 297);

    auto g_xzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 298);

    auto g_xzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 299);

    auto g_xzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 300);

    auto g_xzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 301);

    auto g_xzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 302);

    auto g_xzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 303);

    auto g_xzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 304);

    auto g_xzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 305);

    auto g_xzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 306);

    auto g_xzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 307);

    auto g_xzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 308);

    auto g_xzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 309);

    auto g_xzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 310);

    auto g_xzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 311);

    auto g_xzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 312);

    auto g_xzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 313);

    auto g_xzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 314);

    #pragma omp simd aligned(g_xzzzz_xxxxx_0, g_xzzzz_xxxxy_0, g_xzzzz_xxxxz_0, g_xzzzz_xxxyy_0, g_xzzzz_xxxyz_0, g_xzzzz_xxxzz_0, g_xzzzz_xxyyy_0, g_xzzzz_xxyyz_0, g_xzzzz_xxyzz_0, g_xzzzz_xxzzz_0, g_xzzzz_xyyyy_0, g_xzzzz_xyyyz_0, g_xzzzz_xyyzz_0, g_xzzzz_xyzzz_0, g_xzzzz_xzzzz_0, g_xzzzz_yyyyy_0, g_xzzzz_yyyyz_0, g_xzzzz_yyyzz_0, g_xzzzz_yyzzz_0, g_xzzzz_yzzzz_0, g_xzzzz_zzzzz_0, g_zzzz_xxxx_1, g_zzzz_xxxxx_1, g_zzzz_xxxxy_1, g_zzzz_xxxxz_1, g_zzzz_xxxy_1, g_zzzz_xxxyy_1, g_zzzz_xxxyz_1, g_zzzz_xxxz_1, g_zzzz_xxxzz_1, g_zzzz_xxyy_1, g_zzzz_xxyyy_1, g_zzzz_xxyyz_1, g_zzzz_xxyz_1, g_zzzz_xxyzz_1, g_zzzz_xxzz_1, g_zzzz_xxzzz_1, g_zzzz_xyyy_1, g_zzzz_xyyyy_1, g_zzzz_xyyyz_1, g_zzzz_xyyz_1, g_zzzz_xyyzz_1, g_zzzz_xyzz_1, g_zzzz_xyzzz_1, g_zzzz_xzzz_1, g_zzzz_xzzzz_1, g_zzzz_yyyy_1, g_zzzz_yyyyy_1, g_zzzz_yyyyz_1, g_zzzz_yyyz_1, g_zzzz_yyyzz_1, g_zzzz_yyzz_1, g_zzzz_yyzzz_1, g_zzzz_yzzz_1, g_zzzz_yzzzz_1, g_zzzz_zzzz_1, g_zzzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzz_xxxxx_0[i] = 5.0 * g_zzzz_xxxx_1[i] * fe_0 + g_zzzz_xxxxx_1[i] * pa_x[i];

        g_xzzzz_xxxxy_0[i] = 4.0 * g_zzzz_xxxy_1[i] * fe_0 + g_zzzz_xxxxy_1[i] * pa_x[i];

        g_xzzzz_xxxxz_0[i] = 4.0 * g_zzzz_xxxz_1[i] * fe_0 + g_zzzz_xxxxz_1[i] * pa_x[i];

        g_xzzzz_xxxyy_0[i] = 3.0 * g_zzzz_xxyy_1[i] * fe_0 + g_zzzz_xxxyy_1[i] * pa_x[i];

        g_xzzzz_xxxyz_0[i] = 3.0 * g_zzzz_xxyz_1[i] * fe_0 + g_zzzz_xxxyz_1[i] * pa_x[i];

        g_xzzzz_xxxzz_0[i] = 3.0 * g_zzzz_xxzz_1[i] * fe_0 + g_zzzz_xxxzz_1[i] * pa_x[i];

        g_xzzzz_xxyyy_0[i] = 2.0 * g_zzzz_xyyy_1[i] * fe_0 + g_zzzz_xxyyy_1[i] * pa_x[i];

        g_xzzzz_xxyyz_0[i] = 2.0 * g_zzzz_xyyz_1[i] * fe_0 + g_zzzz_xxyyz_1[i] * pa_x[i];

        g_xzzzz_xxyzz_0[i] = 2.0 * g_zzzz_xyzz_1[i] * fe_0 + g_zzzz_xxyzz_1[i] * pa_x[i];

        g_xzzzz_xxzzz_0[i] = 2.0 * g_zzzz_xzzz_1[i] * fe_0 + g_zzzz_xxzzz_1[i] * pa_x[i];

        g_xzzzz_xyyyy_0[i] = g_zzzz_yyyy_1[i] * fe_0 + g_zzzz_xyyyy_1[i] * pa_x[i];

        g_xzzzz_xyyyz_0[i] = g_zzzz_yyyz_1[i] * fe_0 + g_zzzz_xyyyz_1[i] * pa_x[i];

        g_xzzzz_xyyzz_0[i] = g_zzzz_yyzz_1[i] * fe_0 + g_zzzz_xyyzz_1[i] * pa_x[i];

        g_xzzzz_xyzzz_0[i] = g_zzzz_yzzz_1[i] * fe_0 + g_zzzz_xyzzz_1[i] * pa_x[i];

        g_xzzzz_xzzzz_0[i] = g_zzzz_zzzz_1[i] * fe_0 + g_zzzz_xzzzz_1[i] * pa_x[i];

        g_xzzzz_yyyyy_0[i] = g_zzzz_yyyyy_1[i] * pa_x[i];

        g_xzzzz_yyyyz_0[i] = g_zzzz_yyyyz_1[i] * pa_x[i];

        g_xzzzz_yyyzz_0[i] = g_zzzz_yyyzz_1[i] * pa_x[i];

        g_xzzzz_yyzzz_0[i] = g_zzzz_yyzzz_1[i] * pa_x[i];

        g_xzzzz_yzzzz_0[i] = g_zzzz_yzzzz_1[i] * pa_x[i];

        g_xzzzz_zzzzz_0[i] = g_zzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : HH

    auto g_yyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 315);

    auto g_yyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 316);

    auto g_yyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 317);

    auto g_yyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 318);

    auto g_yyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 319);

    auto g_yyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 320);

    auto g_yyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 321);

    auto g_yyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 322);

    auto g_yyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 323);

    auto g_yyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 324);

    auto g_yyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 325);

    auto g_yyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 326);

    auto g_yyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 327);

    auto g_yyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 328);

    auto g_yyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 329);

    auto g_yyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 330);

    auto g_yyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 331);

    auto g_yyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 332);

    auto g_yyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 333);

    auto g_yyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 334);

    auto g_yyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 335);

    #pragma omp simd aligned(g_yyy_xxxxx_0, g_yyy_xxxxx_1, g_yyy_xxxxy_0, g_yyy_xxxxy_1, g_yyy_xxxxz_0, g_yyy_xxxxz_1, g_yyy_xxxyy_0, g_yyy_xxxyy_1, g_yyy_xxxyz_0, g_yyy_xxxyz_1, g_yyy_xxxzz_0, g_yyy_xxxzz_1, g_yyy_xxyyy_0, g_yyy_xxyyy_1, g_yyy_xxyyz_0, g_yyy_xxyyz_1, g_yyy_xxyzz_0, g_yyy_xxyzz_1, g_yyy_xxzzz_0, g_yyy_xxzzz_1, g_yyy_xyyyy_0, g_yyy_xyyyy_1, g_yyy_xyyyz_0, g_yyy_xyyyz_1, g_yyy_xyyzz_0, g_yyy_xyyzz_1, g_yyy_xyzzz_0, g_yyy_xyzzz_1, g_yyy_xzzzz_0, g_yyy_xzzzz_1, g_yyy_yyyyy_0, g_yyy_yyyyy_1, g_yyy_yyyyz_0, g_yyy_yyyyz_1, g_yyy_yyyzz_0, g_yyy_yyyzz_1, g_yyy_yyzzz_0, g_yyy_yyzzz_1, g_yyy_yzzzz_0, g_yyy_yzzzz_1, g_yyy_zzzzz_0, g_yyy_zzzzz_1, g_yyyy_xxxx_1, g_yyyy_xxxxx_1, g_yyyy_xxxxy_1, g_yyyy_xxxxz_1, g_yyyy_xxxy_1, g_yyyy_xxxyy_1, g_yyyy_xxxyz_1, g_yyyy_xxxz_1, g_yyyy_xxxzz_1, g_yyyy_xxyy_1, g_yyyy_xxyyy_1, g_yyyy_xxyyz_1, g_yyyy_xxyz_1, g_yyyy_xxyzz_1, g_yyyy_xxzz_1, g_yyyy_xxzzz_1, g_yyyy_xyyy_1, g_yyyy_xyyyy_1, g_yyyy_xyyyz_1, g_yyyy_xyyz_1, g_yyyy_xyyzz_1, g_yyyy_xyzz_1, g_yyyy_xyzzz_1, g_yyyy_xzzz_1, g_yyyy_xzzzz_1, g_yyyy_yyyy_1, g_yyyy_yyyyy_1, g_yyyy_yyyyz_1, g_yyyy_yyyz_1, g_yyyy_yyyzz_1, g_yyyy_yyzz_1, g_yyyy_yyzzz_1, g_yyyy_yzzz_1, g_yyyy_yzzzz_1, g_yyyy_zzzz_1, g_yyyy_zzzzz_1, g_yyyyy_xxxxx_0, g_yyyyy_xxxxy_0, g_yyyyy_xxxxz_0, g_yyyyy_xxxyy_0, g_yyyyy_xxxyz_0, g_yyyyy_xxxzz_0, g_yyyyy_xxyyy_0, g_yyyyy_xxyyz_0, g_yyyyy_xxyzz_0, g_yyyyy_xxzzz_0, g_yyyyy_xyyyy_0, g_yyyyy_xyyyz_0, g_yyyyy_xyyzz_0, g_yyyyy_xyzzz_0, g_yyyyy_xzzzz_0, g_yyyyy_yyyyy_0, g_yyyyy_yyyyz_0, g_yyyyy_yyyzz_0, g_yyyyy_yyzzz_0, g_yyyyy_yzzzz_0, g_yyyyy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyy_xxxxx_0[i] = 4.0 * g_yyy_xxxxx_0[i] * fbe_0 - 4.0 * g_yyy_xxxxx_1[i] * fz_be_0 + g_yyyy_xxxxx_1[i] * pa_y[i];

        g_yyyyy_xxxxy_0[i] = 4.0 * g_yyy_xxxxy_0[i] * fbe_0 - 4.0 * g_yyy_xxxxy_1[i] * fz_be_0 + g_yyyy_xxxx_1[i] * fe_0 + g_yyyy_xxxxy_1[i] * pa_y[i];

        g_yyyyy_xxxxz_0[i] = 4.0 * g_yyy_xxxxz_0[i] * fbe_0 - 4.0 * g_yyy_xxxxz_1[i] * fz_be_0 + g_yyyy_xxxxz_1[i] * pa_y[i];

        g_yyyyy_xxxyy_0[i] = 4.0 * g_yyy_xxxyy_0[i] * fbe_0 - 4.0 * g_yyy_xxxyy_1[i] * fz_be_0 + 2.0 * g_yyyy_xxxy_1[i] * fe_0 + g_yyyy_xxxyy_1[i] * pa_y[i];

        g_yyyyy_xxxyz_0[i] = 4.0 * g_yyy_xxxyz_0[i] * fbe_0 - 4.0 * g_yyy_xxxyz_1[i] * fz_be_0 + g_yyyy_xxxz_1[i] * fe_0 + g_yyyy_xxxyz_1[i] * pa_y[i];

        g_yyyyy_xxxzz_0[i] = 4.0 * g_yyy_xxxzz_0[i] * fbe_0 - 4.0 * g_yyy_xxxzz_1[i] * fz_be_0 + g_yyyy_xxxzz_1[i] * pa_y[i];

        g_yyyyy_xxyyy_0[i] = 4.0 * g_yyy_xxyyy_0[i] * fbe_0 - 4.0 * g_yyy_xxyyy_1[i] * fz_be_0 + 3.0 * g_yyyy_xxyy_1[i] * fe_0 + g_yyyy_xxyyy_1[i] * pa_y[i];

        g_yyyyy_xxyyz_0[i] = 4.0 * g_yyy_xxyyz_0[i] * fbe_0 - 4.0 * g_yyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyy_xxyz_1[i] * fe_0 + g_yyyy_xxyyz_1[i] * pa_y[i];

        g_yyyyy_xxyzz_0[i] = 4.0 * g_yyy_xxyzz_0[i] * fbe_0 - 4.0 * g_yyy_xxyzz_1[i] * fz_be_0 + g_yyyy_xxzz_1[i] * fe_0 + g_yyyy_xxyzz_1[i] * pa_y[i];

        g_yyyyy_xxzzz_0[i] = 4.0 * g_yyy_xxzzz_0[i] * fbe_0 - 4.0 * g_yyy_xxzzz_1[i] * fz_be_0 + g_yyyy_xxzzz_1[i] * pa_y[i];

        g_yyyyy_xyyyy_0[i] = 4.0 * g_yyy_xyyyy_0[i] * fbe_0 - 4.0 * g_yyy_xyyyy_1[i] * fz_be_0 + 4.0 * g_yyyy_xyyy_1[i] * fe_0 + g_yyyy_xyyyy_1[i] * pa_y[i];

        g_yyyyy_xyyyz_0[i] = 4.0 * g_yyy_xyyyz_0[i] * fbe_0 - 4.0 * g_yyy_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyy_xyyz_1[i] * fe_0 + g_yyyy_xyyyz_1[i] * pa_y[i];

        g_yyyyy_xyyzz_0[i] = 4.0 * g_yyy_xyyzz_0[i] * fbe_0 - 4.0 * g_yyy_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyy_xyzz_1[i] * fe_0 + g_yyyy_xyyzz_1[i] * pa_y[i];

        g_yyyyy_xyzzz_0[i] = 4.0 * g_yyy_xyzzz_0[i] * fbe_0 - 4.0 * g_yyy_xyzzz_1[i] * fz_be_0 + g_yyyy_xzzz_1[i] * fe_0 + g_yyyy_xyzzz_1[i] * pa_y[i];

        g_yyyyy_xzzzz_0[i] = 4.0 * g_yyy_xzzzz_0[i] * fbe_0 - 4.0 * g_yyy_xzzzz_1[i] * fz_be_0 + g_yyyy_xzzzz_1[i] * pa_y[i];

        g_yyyyy_yyyyy_0[i] = 4.0 * g_yyy_yyyyy_0[i] * fbe_0 - 4.0 * g_yyy_yyyyy_1[i] * fz_be_0 + 5.0 * g_yyyy_yyyy_1[i] * fe_0 + g_yyyy_yyyyy_1[i] * pa_y[i];

        g_yyyyy_yyyyz_0[i] = 4.0 * g_yyy_yyyyz_0[i] * fbe_0 - 4.0 * g_yyy_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyy_yyyz_1[i] * fe_0 + g_yyyy_yyyyz_1[i] * pa_y[i];

        g_yyyyy_yyyzz_0[i] = 4.0 * g_yyy_yyyzz_0[i] * fbe_0 - 4.0 * g_yyy_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyy_yyzz_1[i] * fe_0 + g_yyyy_yyyzz_1[i] * pa_y[i];

        g_yyyyy_yyzzz_0[i] = 4.0 * g_yyy_yyzzz_0[i] * fbe_0 - 4.0 * g_yyy_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_yzzz_1[i] * fe_0 + g_yyyy_yyzzz_1[i] * pa_y[i];

        g_yyyyy_yzzzz_0[i] = 4.0 * g_yyy_yzzzz_0[i] * fbe_0 - 4.0 * g_yyy_yzzzz_1[i] * fz_be_0 + g_yyyy_zzzz_1[i] * fe_0 + g_yyyy_yzzzz_1[i] * pa_y[i];

        g_yyyyy_zzzzz_0[i] = 4.0 * g_yyy_zzzzz_0[i] * fbe_0 - 4.0 * g_yyy_zzzzz_1[i] * fz_be_0 + g_yyyy_zzzzz_1[i] * pa_y[i];
    }

    // Set up 336-357 components of targeted buffer : HH

    auto g_yyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 336);

    auto g_yyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 337);

    auto g_yyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 338);

    auto g_yyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 339);

    auto g_yyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 340);

    auto g_yyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 341);

    auto g_yyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 342);

    auto g_yyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 343);

    auto g_yyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 344);

    auto g_yyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 345);

    auto g_yyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 346);

    auto g_yyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 347);

    auto g_yyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 348);

    auto g_yyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 349);

    auto g_yyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 350);

    auto g_yyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 351);

    auto g_yyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 352);

    auto g_yyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 353);

    auto g_yyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 354);

    auto g_yyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 355);

    auto g_yyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 356);

    #pragma omp simd aligned(g_yyyy_xxxx_1, g_yyyy_xxxxx_1, g_yyyy_xxxxy_1, g_yyyy_xxxxz_1, g_yyyy_xxxy_1, g_yyyy_xxxyy_1, g_yyyy_xxxyz_1, g_yyyy_xxxz_1, g_yyyy_xxxzz_1, g_yyyy_xxyy_1, g_yyyy_xxyyy_1, g_yyyy_xxyyz_1, g_yyyy_xxyz_1, g_yyyy_xxyzz_1, g_yyyy_xxzz_1, g_yyyy_xxzzz_1, g_yyyy_xyyy_1, g_yyyy_xyyyy_1, g_yyyy_xyyyz_1, g_yyyy_xyyz_1, g_yyyy_xyyzz_1, g_yyyy_xyzz_1, g_yyyy_xyzzz_1, g_yyyy_xzzz_1, g_yyyy_xzzzz_1, g_yyyy_yyyy_1, g_yyyy_yyyyy_1, g_yyyy_yyyyz_1, g_yyyy_yyyz_1, g_yyyy_yyyzz_1, g_yyyy_yyzz_1, g_yyyy_yyzzz_1, g_yyyy_yzzz_1, g_yyyy_yzzzz_1, g_yyyy_zzzz_1, g_yyyy_zzzzz_1, g_yyyyz_xxxxx_0, g_yyyyz_xxxxy_0, g_yyyyz_xxxxz_0, g_yyyyz_xxxyy_0, g_yyyyz_xxxyz_0, g_yyyyz_xxxzz_0, g_yyyyz_xxyyy_0, g_yyyyz_xxyyz_0, g_yyyyz_xxyzz_0, g_yyyyz_xxzzz_0, g_yyyyz_xyyyy_0, g_yyyyz_xyyyz_0, g_yyyyz_xyyzz_0, g_yyyyz_xyzzz_0, g_yyyyz_xzzzz_0, g_yyyyz_yyyyy_0, g_yyyyz_yyyyz_0, g_yyyyz_yyyzz_0, g_yyyyz_yyzzz_0, g_yyyyz_yzzzz_0, g_yyyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyz_xxxxx_0[i] = g_yyyy_xxxxx_1[i] * pa_z[i];

        g_yyyyz_xxxxy_0[i] = g_yyyy_xxxxy_1[i] * pa_z[i];

        g_yyyyz_xxxxz_0[i] = g_yyyy_xxxx_1[i] * fe_0 + g_yyyy_xxxxz_1[i] * pa_z[i];

        g_yyyyz_xxxyy_0[i] = g_yyyy_xxxyy_1[i] * pa_z[i];

        g_yyyyz_xxxyz_0[i] = g_yyyy_xxxy_1[i] * fe_0 + g_yyyy_xxxyz_1[i] * pa_z[i];

        g_yyyyz_xxxzz_0[i] = 2.0 * g_yyyy_xxxz_1[i] * fe_0 + g_yyyy_xxxzz_1[i] * pa_z[i];

        g_yyyyz_xxyyy_0[i] = g_yyyy_xxyyy_1[i] * pa_z[i];

        g_yyyyz_xxyyz_0[i] = g_yyyy_xxyy_1[i] * fe_0 + g_yyyy_xxyyz_1[i] * pa_z[i];

        g_yyyyz_xxyzz_0[i] = 2.0 * g_yyyy_xxyz_1[i] * fe_0 + g_yyyy_xxyzz_1[i] * pa_z[i];

        g_yyyyz_xxzzz_0[i] = 3.0 * g_yyyy_xxzz_1[i] * fe_0 + g_yyyy_xxzzz_1[i] * pa_z[i];

        g_yyyyz_xyyyy_0[i] = g_yyyy_xyyyy_1[i] * pa_z[i];

        g_yyyyz_xyyyz_0[i] = g_yyyy_xyyy_1[i] * fe_0 + g_yyyy_xyyyz_1[i] * pa_z[i];

        g_yyyyz_xyyzz_0[i] = 2.0 * g_yyyy_xyyz_1[i] * fe_0 + g_yyyy_xyyzz_1[i] * pa_z[i];

        g_yyyyz_xyzzz_0[i] = 3.0 * g_yyyy_xyzz_1[i] * fe_0 + g_yyyy_xyzzz_1[i] * pa_z[i];

        g_yyyyz_xzzzz_0[i] = 4.0 * g_yyyy_xzzz_1[i] * fe_0 + g_yyyy_xzzzz_1[i] * pa_z[i];

        g_yyyyz_yyyyy_0[i] = g_yyyy_yyyyy_1[i] * pa_z[i];

        g_yyyyz_yyyyz_0[i] = g_yyyy_yyyy_1[i] * fe_0 + g_yyyy_yyyyz_1[i] * pa_z[i];

        g_yyyyz_yyyzz_0[i] = 2.0 * g_yyyy_yyyz_1[i] * fe_0 + g_yyyy_yyyzz_1[i] * pa_z[i];

        g_yyyyz_yyzzz_0[i] = 3.0 * g_yyyy_yyzz_1[i] * fe_0 + g_yyyy_yyzzz_1[i] * pa_z[i];

        g_yyyyz_yzzzz_0[i] = 4.0 * g_yyyy_yzzz_1[i] * fe_0 + g_yyyy_yzzzz_1[i] * pa_z[i];

        g_yyyyz_zzzzz_0[i] = 5.0 * g_yyyy_zzzz_1[i] * fe_0 + g_yyyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 357-378 components of targeted buffer : HH

    auto g_yyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 357);

    auto g_yyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 358);

    auto g_yyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 359);

    auto g_yyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 360);

    auto g_yyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 361);

    auto g_yyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 362);

    auto g_yyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 363);

    auto g_yyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 364);

    auto g_yyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 365);

    auto g_yyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 366);

    auto g_yyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 367);

    auto g_yyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 368);

    auto g_yyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 369);

    auto g_yyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 370);

    auto g_yyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 371);

    auto g_yyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 372);

    auto g_yyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 373);

    auto g_yyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 374);

    auto g_yyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 375);

    auto g_yyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 376);

    auto g_yyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 377);

    #pragma omp simd aligned(g_yyy_xxxxy_0, g_yyy_xxxxy_1, g_yyy_xxxyy_0, g_yyy_xxxyy_1, g_yyy_xxyyy_0, g_yyy_xxyyy_1, g_yyy_xyyyy_0, g_yyy_xyyyy_1, g_yyy_yyyyy_0, g_yyy_yyyyy_1, g_yyyz_xxxxy_1, g_yyyz_xxxyy_1, g_yyyz_xxyyy_1, g_yyyz_xyyyy_1, g_yyyz_yyyyy_1, g_yyyzz_xxxxx_0, g_yyyzz_xxxxy_0, g_yyyzz_xxxxz_0, g_yyyzz_xxxyy_0, g_yyyzz_xxxyz_0, g_yyyzz_xxxzz_0, g_yyyzz_xxyyy_0, g_yyyzz_xxyyz_0, g_yyyzz_xxyzz_0, g_yyyzz_xxzzz_0, g_yyyzz_xyyyy_0, g_yyyzz_xyyyz_0, g_yyyzz_xyyzz_0, g_yyyzz_xyzzz_0, g_yyyzz_xzzzz_0, g_yyyzz_yyyyy_0, g_yyyzz_yyyyz_0, g_yyyzz_yyyzz_0, g_yyyzz_yyzzz_0, g_yyyzz_yzzzz_0, g_yyyzz_zzzzz_0, g_yyzz_xxxxx_1, g_yyzz_xxxxz_1, g_yyzz_xxxyz_1, g_yyzz_xxxz_1, g_yyzz_xxxzz_1, g_yyzz_xxyyz_1, g_yyzz_xxyz_1, g_yyzz_xxyzz_1, g_yyzz_xxzz_1, g_yyzz_xxzzz_1, g_yyzz_xyyyz_1, g_yyzz_xyyz_1, g_yyzz_xyyzz_1, g_yyzz_xyzz_1, g_yyzz_xyzzz_1, g_yyzz_xzzz_1, g_yyzz_xzzzz_1, g_yyzz_yyyyz_1, g_yyzz_yyyz_1, g_yyzz_yyyzz_1, g_yyzz_yyzz_1, g_yyzz_yyzzz_1, g_yyzz_yzzz_1, g_yyzz_yzzzz_1, g_yyzz_zzzz_1, g_yyzz_zzzzz_1, g_yzz_xxxxx_0, g_yzz_xxxxx_1, g_yzz_xxxxz_0, g_yzz_xxxxz_1, g_yzz_xxxyz_0, g_yzz_xxxyz_1, g_yzz_xxxzz_0, g_yzz_xxxzz_1, g_yzz_xxyyz_0, g_yzz_xxyyz_1, g_yzz_xxyzz_0, g_yzz_xxyzz_1, g_yzz_xxzzz_0, g_yzz_xxzzz_1, g_yzz_xyyyz_0, g_yzz_xyyyz_1, g_yzz_xyyzz_0, g_yzz_xyyzz_1, g_yzz_xyzzz_0, g_yzz_xyzzz_1, g_yzz_xzzzz_0, g_yzz_xzzzz_1, g_yzz_yyyyz_0, g_yzz_yyyyz_1, g_yzz_yyyzz_0, g_yzz_yyyzz_1, g_yzz_yyzzz_0, g_yzz_yyzzz_1, g_yzz_yzzzz_0, g_yzz_yzzzz_1, g_yzz_zzzzz_0, g_yzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzz_xxxxx_0[i] = 2.0 * g_yzz_xxxxx_0[i] * fbe_0 - 2.0 * g_yzz_xxxxx_1[i] * fz_be_0 + g_yyzz_xxxxx_1[i] * pa_y[i];

        g_yyyzz_xxxxy_0[i] = g_yyy_xxxxy_0[i] * fbe_0 - g_yyy_xxxxy_1[i] * fz_be_0 + g_yyyz_xxxxy_1[i] * pa_z[i];

        g_yyyzz_xxxxz_0[i] = 2.0 * g_yzz_xxxxz_0[i] * fbe_0 - 2.0 * g_yzz_xxxxz_1[i] * fz_be_0 + g_yyzz_xxxxz_1[i] * pa_y[i];

        g_yyyzz_xxxyy_0[i] = g_yyy_xxxyy_0[i] * fbe_0 - g_yyy_xxxyy_1[i] * fz_be_0 + g_yyyz_xxxyy_1[i] * pa_z[i];

        g_yyyzz_xxxyz_0[i] = 2.0 * g_yzz_xxxyz_0[i] * fbe_0 - 2.0 * g_yzz_xxxyz_1[i] * fz_be_0 + g_yyzz_xxxz_1[i] * fe_0 + g_yyzz_xxxyz_1[i] * pa_y[i];

        g_yyyzz_xxxzz_0[i] = 2.0 * g_yzz_xxxzz_0[i] * fbe_0 - 2.0 * g_yzz_xxxzz_1[i] * fz_be_0 + g_yyzz_xxxzz_1[i] * pa_y[i];

        g_yyyzz_xxyyy_0[i] = g_yyy_xxyyy_0[i] * fbe_0 - g_yyy_xxyyy_1[i] * fz_be_0 + g_yyyz_xxyyy_1[i] * pa_z[i];

        g_yyyzz_xxyyz_0[i] = 2.0 * g_yzz_xxyyz_0[i] * fbe_0 - 2.0 * g_yzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyzz_xxyz_1[i] * fe_0 + g_yyzz_xxyyz_1[i] * pa_y[i];

        g_yyyzz_xxyzz_0[i] = 2.0 * g_yzz_xxyzz_0[i] * fbe_0 - 2.0 * g_yzz_xxyzz_1[i] * fz_be_0 + g_yyzz_xxzz_1[i] * fe_0 + g_yyzz_xxyzz_1[i] * pa_y[i];

        g_yyyzz_xxzzz_0[i] = 2.0 * g_yzz_xxzzz_0[i] * fbe_0 - 2.0 * g_yzz_xxzzz_1[i] * fz_be_0 + g_yyzz_xxzzz_1[i] * pa_y[i];

        g_yyyzz_xyyyy_0[i] = g_yyy_xyyyy_0[i] * fbe_0 - g_yyy_xyyyy_1[i] * fz_be_0 + g_yyyz_xyyyy_1[i] * pa_z[i];

        g_yyyzz_xyyyz_0[i] = 2.0 * g_yzz_xyyyz_0[i] * fbe_0 - 2.0 * g_yzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyzz_xyyz_1[i] * fe_0 + g_yyzz_xyyyz_1[i] * pa_y[i];

        g_yyyzz_xyyzz_0[i] = 2.0 * g_yzz_xyyzz_0[i] * fbe_0 - 2.0 * g_yzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyzz_xyzz_1[i] * fe_0 + g_yyzz_xyyzz_1[i] * pa_y[i];

        g_yyyzz_xyzzz_0[i] = 2.0 * g_yzz_xyzzz_0[i] * fbe_0 - 2.0 * g_yzz_xyzzz_1[i] * fz_be_0 + g_yyzz_xzzz_1[i] * fe_0 + g_yyzz_xyzzz_1[i] * pa_y[i];

        g_yyyzz_xzzzz_0[i] = 2.0 * g_yzz_xzzzz_0[i] * fbe_0 - 2.0 * g_yzz_xzzzz_1[i] * fz_be_0 + g_yyzz_xzzzz_1[i] * pa_y[i];

        g_yyyzz_yyyyy_0[i] = g_yyy_yyyyy_0[i] * fbe_0 - g_yyy_yyyyy_1[i] * fz_be_0 + g_yyyz_yyyyy_1[i] * pa_z[i];

        g_yyyzz_yyyyz_0[i] = 2.0 * g_yzz_yyyyz_0[i] * fbe_0 - 2.0 * g_yzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyzz_yyyz_1[i] * fe_0 + g_yyzz_yyyyz_1[i] * pa_y[i];

        g_yyyzz_yyyzz_0[i] = 2.0 * g_yzz_yyyzz_0[i] * fbe_0 - 2.0 * g_yzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyzz_yyzz_1[i] * fe_0 + g_yyzz_yyyzz_1[i] * pa_y[i];

        g_yyyzz_yyzzz_0[i] = 2.0 * g_yzz_yyzzz_0[i] * fbe_0 - 2.0 * g_yzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_yzzz_1[i] * fe_0 + g_yyzz_yyzzz_1[i] * pa_y[i];

        g_yyyzz_yzzzz_0[i] = 2.0 * g_yzz_yzzzz_0[i] * fbe_0 - 2.0 * g_yzz_yzzzz_1[i] * fz_be_0 + g_yyzz_zzzz_1[i] * fe_0 + g_yyzz_yzzzz_1[i] * pa_y[i];

        g_yyyzz_zzzzz_0[i] = 2.0 * g_yzz_zzzzz_0[i] * fbe_0 - 2.0 * g_yzz_zzzzz_1[i] * fz_be_0 + g_yyzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 378-399 components of targeted buffer : HH

    auto g_yyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 378);

    auto g_yyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 379);

    auto g_yyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 380);

    auto g_yyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 381);

    auto g_yyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 382);

    auto g_yyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 383);

    auto g_yyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 384);

    auto g_yyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 385);

    auto g_yyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 386);

    auto g_yyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 387);

    auto g_yyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 388);

    auto g_yyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 389);

    auto g_yyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 390);

    auto g_yyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 391);

    auto g_yyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 392);

    auto g_yyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 393);

    auto g_yyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 394);

    auto g_yyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 395);

    auto g_yyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 396);

    auto g_yyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 397);

    auto g_yyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 398);

    #pragma omp simd aligned(g_yyz_xxxxy_0, g_yyz_xxxxy_1, g_yyz_xxxyy_0, g_yyz_xxxyy_1, g_yyz_xxyyy_0, g_yyz_xxyyy_1, g_yyz_xyyyy_0, g_yyz_xyyyy_1, g_yyz_yyyyy_0, g_yyz_yyyyy_1, g_yyzz_xxxxy_1, g_yyzz_xxxyy_1, g_yyzz_xxyyy_1, g_yyzz_xyyyy_1, g_yyzz_yyyyy_1, g_yyzzz_xxxxx_0, g_yyzzz_xxxxy_0, g_yyzzz_xxxxz_0, g_yyzzz_xxxyy_0, g_yyzzz_xxxyz_0, g_yyzzz_xxxzz_0, g_yyzzz_xxyyy_0, g_yyzzz_xxyyz_0, g_yyzzz_xxyzz_0, g_yyzzz_xxzzz_0, g_yyzzz_xyyyy_0, g_yyzzz_xyyyz_0, g_yyzzz_xyyzz_0, g_yyzzz_xyzzz_0, g_yyzzz_xzzzz_0, g_yyzzz_yyyyy_0, g_yyzzz_yyyyz_0, g_yyzzz_yyyzz_0, g_yyzzz_yyzzz_0, g_yyzzz_yzzzz_0, g_yyzzz_zzzzz_0, g_yzzz_xxxxx_1, g_yzzz_xxxxz_1, g_yzzz_xxxyz_1, g_yzzz_xxxz_1, g_yzzz_xxxzz_1, g_yzzz_xxyyz_1, g_yzzz_xxyz_1, g_yzzz_xxyzz_1, g_yzzz_xxzz_1, g_yzzz_xxzzz_1, g_yzzz_xyyyz_1, g_yzzz_xyyz_1, g_yzzz_xyyzz_1, g_yzzz_xyzz_1, g_yzzz_xyzzz_1, g_yzzz_xzzz_1, g_yzzz_xzzzz_1, g_yzzz_yyyyz_1, g_yzzz_yyyz_1, g_yzzz_yyyzz_1, g_yzzz_yyzz_1, g_yzzz_yyzzz_1, g_yzzz_yzzz_1, g_yzzz_yzzzz_1, g_yzzz_zzzz_1, g_yzzz_zzzzz_1, g_zzz_xxxxx_0, g_zzz_xxxxx_1, g_zzz_xxxxz_0, g_zzz_xxxxz_1, g_zzz_xxxyz_0, g_zzz_xxxyz_1, g_zzz_xxxzz_0, g_zzz_xxxzz_1, g_zzz_xxyyz_0, g_zzz_xxyyz_1, g_zzz_xxyzz_0, g_zzz_xxyzz_1, g_zzz_xxzzz_0, g_zzz_xxzzz_1, g_zzz_xyyyz_0, g_zzz_xyyyz_1, g_zzz_xyyzz_0, g_zzz_xyyzz_1, g_zzz_xyzzz_0, g_zzz_xyzzz_1, g_zzz_xzzzz_0, g_zzz_xzzzz_1, g_zzz_yyyyz_0, g_zzz_yyyyz_1, g_zzz_yyyzz_0, g_zzz_yyyzz_1, g_zzz_yyzzz_0, g_zzz_yyzzz_1, g_zzz_yzzzz_0, g_zzz_yzzzz_1, g_zzz_zzzzz_0, g_zzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzz_xxxxx_0[i] = g_zzz_xxxxx_0[i] * fbe_0 - g_zzz_xxxxx_1[i] * fz_be_0 + g_yzzz_xxxxx_1[i] * pa_y[i];

        g_yyzzz_xxxxy_0[i] = 2.0 * g_yyz_xxxxy_0[i] * fbe_0 - 2.0 * g_yyz_xxxxy_1[i] * fz_be_0 + g_yyzz_xxxxy_1[i] * pa_z[i];

        g_yyzzz_xxxxz_0[i] = g_zzz_xxxxz_0[i] * fbe_0 - g_zzz_xxxxz_1[i] * fz_be_0 + g_yzzz_xxxxz_1[i] * pa_y[i];

        g_yyzzz_xxxyy_0[i] = 2.0 * g_yyz_xxxyy_0[i] * fbe_0 - 2.0 * g_yyz_xxxyy_1[i] * fz_be_0 + g_yyzz_xxxyy_1[i] * pa_z[i];

        g_yyzzz_xxxyz_0[i] = g_zzz_xxxyz_0[i] * fbe_0 - g_zzz_xxxyz_1[i] * fz_be_0 + g_yzzz_xxxz_1[i] * fe_0 + g_yzzz_xxxyz_1[i] * pa_y[i];

        g_yyzzz_xxxzz_0[i] = g_zzz_xxxzz_0[i] * fbe_0 - g_zzz_xxxzz_1[i] * fz_be_0 + g_yzzz_xxxzz_1[i] * pa_y[i];

        g_yyzzz_xxyyy_0[i] = 2.0 * g_yyz_xxyyy_0[i] * fbe_0 - 2.0 * g_yyz_xxyyy_1[i] * fz_be_0 + g_yyzz_xxyyy_1[i] * pa_z[i];

        g_yyzzz_xxyyz_0[i] = g_zzz_xxyyz_0[i] * fbe_0 - g_zzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yzzz_xxyz_1[i] * fe_0 + g_yzzz_xxyyz_1[i] * pa_y[i];

        g_yyzzz_xxyzz_0[i] = g_zzz_xxyzz_0[i] * fbe_0 - g_zzz_xxyzz_1[i] * fz_be_0 + g_yzzz_xxzz_1[i] * fe_0 + g_yzzz_xxyzz_1[i] * pa_y[i];

        g_yyzzz_xxzzz_0[i] = g_zzz_xxzzz_0[i] * fbe_0 - g_zzz_xxzzz_1[i] * fz_be_0 + g_yzzz_xxzzz_1[i] * pa_y[i];

        g_yyzzz_xyyyy_0[i] = 2.0 * g_yyz_xyyyy_0[i] * fbe_0 - 2.0 * g_yyz_xyyyy_1[i] * fz_be_0 + g_yyzz_xyyyy_1[i] * pa_z[i];

        g_yyzzz_xyyyz_0[i] = g_zzz_xyyyz_0[i] * fbe_0 - g_zzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yzzz_xyyz_1[i] * fe_0 + g_yzzz_xyyyz_1[i] * pa_y[i];

        g_yyzzz_xyyzz_0[i] = g_zzz_xyyzz_0[i] * fbe_0 - g_zzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yzzz_xyzz_1[i] * fe_0 + g_yzzz_xyyzz_1[i] * pa_y[i];

        g_yyzzz_xyzzz_0[i] = g_zzz_xyzzz_0[i] * fbe_0 - g_zzz_xyzzz_1[i] * fz_be_0 + g_yzzz_xzzz_1[i] * fe_0 + g_yzzz_xyzzz_1[i] * pa_y[i];

        g_yyzzz_xzzzz_0[i] = g_zzz_xzzzz_0[i] * fbe_0 - g_zzz_xzzzz_1[i] * fz_be_0 + g_yzzz_xzzzz_1[i] * pa_y[i];

        g_yyzzz_yyyyy_0[i] = 2.0 * g_yyz_yyyyy_0[i] * fbe_0 - 2.0 * g_yyz_yyyyy_1[i] * fz_be_0 + g_yyzz_yyyyy_1[i] * pa_z[i];

        g_yyzzz_yyyyz_0[i] = g_zzz_yyyyz_0[i] * fbe_0 - g_zzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yzzz_yyyz_1[i] * fe_0 + g_yzzz_yyyyz_1[i] * pa_y[i];

        g_yyzzz_yyyzz_0[i] = g_zzz_yyyzz_0[i] * fbe_0 - g_zzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yzzz_yyzz_1[i] * fe_0 + g_yzzz_yyyzz_1[i] * pa_y[i];

        g_yyzzz_yyzzz_0[i] = g_zzz_yyzzz_0[i] * fbe_0 - g_zzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_yzzz_1[i] * fe_0 + g_yzzz_yyzzz_1[i] * pa_y[i];

        g_yyzzz_yzzzz_0[i] = g_zzz_yzzzz_0[i] * fbe_0 - g_zzz_yzzzz_1[i] * fz_be_0 + g_yzzz_zzzz_1[i] * fe_0 + g_yzzz_yzzzz_1[i] * pa_y[i];

        g_yyzzz_zzzzz_0[i] = g_zzz_zzzzz_0[i] * fbe_0 - g_zzz_zzzzz_1[i] * fz_be_0 + g_yzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 399-420 components of targeted buffer : HH

    auto g_yzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 399);

    auto g_yzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 400);

    auto g_yzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 401);

    auto g_yzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 402);

    auto g_yzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 403);

    auto g_yzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 404);

    auto g_yzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 405);

    auto g_yzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 406);

    auto g_yzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 407);

    auto g_yzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 408);

    auto g_yzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 409);

    auto g_yzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 410);

    auto g_yzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 411);

    auto g_yzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 412);

    auto g_yzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 413);

    auto g_yzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 414);

    auto g_yzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 415);

    auto g_yzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 416);

    auto g_yzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 417);

    auto g_yzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 418);

    auto g_yzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 419);

    #pragma omp simd aligned(g_yzzzz_xxxxx_0, g_yzzzz_xxxxy_0, g_yzzzz_xxxxz_0, g_yzzzz_xxxyy_0, g_yzzzz_xxxyz_0, g_yzzzz_xxxzz_0, g_yzzzz_xxyyy_0, g_yzzzz_xxyyz_0, g_yzzzz_xxyzz_0, g_yzzzz_xxzzz_0, g_yzzzz_xyyyy_0, g_yzzzz_xyyyz_0, g_yzzzz_xyyzz_0, g_yzzzz_xyzzz_0, g_yzzzz_xzzzz_0, g_yzzzz_yyyyy_0, g_yzzzz_yyyyz_0, g_yzzzz_yyyzz_0, g_yzzzz_yyzzz_0, g_yzzzz_yzzzz_0, g_yzzzz_zzzzz_0, g_zzzz_xxxx_1, g_zzzz_xxxxx_1, g_zzzz_xxxxy_1, g_zzzz_xxxxz_1, g_zzzz_xxxy_1, g_zzzz_xxxyy_1, g_zzzz_xxxyz_1, g_zzzz_xxxz_1, g_zzzz_xxxzz_1, g_zzzz_xxyy_1, g_zzzz_xxyyy_1, g_zzzz_xxyyz_1, g_zzzz_xxyz_1, g_zzzz_xxyzz_1, g_zzzz_xxzz_1, g_zzzz_xxzzz_1, g_zzzz_xyyy_1, g_zzzz_xyyyy_1, g_zzzz_xyyyz_1, g_zzzz_xyyz_1, g_zzzz_xyyzz_1, g_zzzz_xyzz_1, g_zzzz_xyzzz_1, g_zzzz_xzzz_1, g_zzzz_xzzzz_1, g_zzzz_yyyy_1, g_zzzz_yyyyy_1, g_zzzz_yyyyz_1, g_zzzz_yyyz_1, g_zzzz_yyyzz_1, g_zzzz_yyzz_1, g_zzzz_yyzzz_1, g_zzzz_yzzz_1, g_zzzz_yzzzz_1, g_zzzz_zzzz_1, g_zzzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzz_xxxxx_0[i] = g_zzzz_xxxxx_1[i] * pa_y[i];

        g_yzzzz_xxxxy_0[i] = g_zzzz_xxxx_1[i] * fe_0 + g_zzzz_xxxxy_1[i] * pa_y[i];

        g_yzzzz_xxxxz_0[i] = g_zzzz_xxxxz_1[i] * pa_y[i];

        g_yzzzz_xxxyy_0[i] = 2.0 * g_zzzz_xxxy_1[i] * fe_0 + g_zzzz_xxxyy_1[i] * pa_y[i];

        g_yzzzz_xxxyz_0[i] = g_zzzz_xxxz_1[i] * fe_0 + g_zzzz_xxxyz_1[i] * pa_y[i];

        g_yzzzz_xxxzz_0[i] = g_zzzz_xxxzz_1[i] * pa_y[i];

        g_yzzzz_xxyyy_0[i] = 3.0 * g_zzzz_xxyy_1[i] * fe_0 + g_zzzz_xxyyy_1[i] * pa_y[i];

        g_yzzzz_xxyyz_0[i] = 2.0 * g_zzzz_xxyz_1[i] * fe_0 + g_zzzz_xxyyz_1[i] * pa_y[i];

        g_yzzzz_xxyzz_0[i] = g_zzzz_xxzz_1[i] * fe_0 + g_zzzz_xxyzz_1[i] * pa_y[i];

        g_yzzzz_xxzzz_0[i] = g_zzzz_xxzzz_1[i] * pa_y[i];

        g_yzzzz_xyyyy_0[i] = 4.0 * g_zzzz_xyyy_1[i] * fe_0 + g_zzzz_xyyyy_1[i] * pa_y[i];

        g_yzzzz_xyyyz_0[i] = 3.0 * g_zzzz_xyyz_1[i] * fe_0 + g_zzzz_xyyyz_1[i] * pa_y[i];

        g_yzzzz_xyyzz_0[i] = 2.0 * g_zzzz_xyzz_1[i] * fe_0 + g_zzzz_xyyzz_1[i] * pa_y[i];

        g_yzzzz_xyzzz_0[i] = g_zzzz_xzzz_1[i] * fe_0 + g_zzzz_xyzzz_1[i] * pa_y[i];

        g_yzzzz_xzzzz_0[i] = g_zzzz_xzzzz_1[i] * pa_y[i];

        g_yzzzz_yyyyy_0[i] = 5.0 * g_zzzz_yyyy_1[i] * fe_0 + g_zzzz_yyyyy_1[i] * pa_y[i];

        g_yzzzz_yyyyz_0[i] = 4.0 * g_zzzz_yyyz_1[i] * fe_0 + g_zzzz_yyyyz_1[i] * pa_y[i];

        g_yzzzz_yyyzz_0[i] = 3.0 * g_zzzz_yyzz_1[i] * fe_0 + g_zzzz_yyyzz_1[i] * pa_y[i];

        g_yzzzz_yyzzz_0[i] = 2.0 * g_zzzz_yzzz_1[i] * fe_0 + g_zzzz_yyzzz_1[i] * pa_y[i];

        g_yzzzz_yzzzz_0[i] = g_zzzz_zzzz_1[i] * fe_0 + g_zzzz_yzzzz_1[i] * pa_y[i];

        g_yzzzz_zzzzz_0[i] = g_zzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 420-441 components of targeted buffer : HH

    auto g_zzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 420);

    auto g_zzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 421);

    auto g_zzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 422);

    auto g_zzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 423);

    auto g_zzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 424);

    auto g_zzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 425);

    auto g_zzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 426);

    auto g_zzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 427);

    auto g_zzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 428);

    auto g_zzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 429);

    auto g_zzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 430);

    auto g_zzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 431);

    auto g_zzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 432);

    auto g_zzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 433);

    auto g_zzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 434);

    auto g_zzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 435);

    auto g_zzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 436);

    auto g_zzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 437);

    auto g_zzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 438);

    auto g_zzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 439);

    auto g_zzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 440);

    #pragma omp simd aligned(g_zzz_xxxxx_0, g_zzz_xxxxx_1, g_zzz_xxxxy_0, g_zzz_xxxxy_1, g_zzz_xxxxz_0, g_zzz_xxxxz_1, g_zzz_xxxyy_0, g_zzz_xxxyy_1, g_zzz_xxxyz_0, g_zzz_xxxyz_1, g_zzz_xxxzz_0, g_zzz_xxxzz_1, g_zzz_xxyyy_0, g_zzz_xxyyy_1, g_zzz_xxyyz_0, g_zzz_xxyyz_1, g_zzz_xxyzz_0, g_zzz_xxyzz_1, g_zzz_xxzzz_0, g_zzz_xxzzz_1, g_zzz_xyyyy_0, g_zzz_xyyyy_1, g_zzz_xyyyz_0, g_zzz_xyyyz_1, g_zzz_xyyzz_0, g_zzz_xyyzz_1, g_zzz_xyzzz_0, g_zzz_xyzzz_1, g_zzz_xzzzz_0, g_zzz_xzzzz_1, g_zzz_yyyyy_0, g_zzz_yyyyy_1, g_zzz_yyyyz_0, g_zzz_yyyyz_1, g_zzz_yyyzz_0, g_zzz_yyyzz_1, g_zzz_yyzzz_0, g_zzz_yyzzz_1, g_zzz_yzzzz_0, g_zzz_yzzzz_1, g_zzz_zzzzz_0, g_zzz_zzzzz_1, g_zzzz_xxxx_1, g_zzzz_xxxxx_1, g_zzzz_xxxxy_1, g_zzzz_xxxxz_1, g_zzzz_xxxy_1, g_zzzz_xxxyy_1, g_zzzz_xxxyz_1, g_zzzz_xxxz_1, g_zzzz_xxxzz_1, g_zzzz_xxyy_1, g_zzzz_xxyyy_1, g_zzzz_xxyyz_1, g_zzzz_xxyz_1, g_zzzz_xxyzz_1, g_zzzz_xxzz_1, g_zzzz_xxzzz_1, g_zzzz_xyyy_1, g_zzzz_xyyyy_1, g_zzzz_xyyyz_1, g_zzzz_xyyz_1, g_zzzz_xyyzz_1, g_zzzz_xyzz_1, g_zzzz_xyzzz_1, g_zzzz_xzzz_1, g_zzzz_xzzzz_1, g_zzzz_yyyy_1, g_zzzz_yyyyy_1, g_zzzz_yyyyz_1, g_zzzz_yyyz_1, g_zzzz_yyyzz_1, g_zzzz_yyzz_1, g_zzzz_yyzzz_1, g_zzzz_yzzz_1, g_zzzz_yzzzz_1, g_zzzz_zzzz_1, g_zzzz_zzzzz_1, g_zzzzz_xxxxx_0, g_zzzzz_xxxxy_0, g_zzzzz_xxxxz_0, g_zzzzz_xxxyy_0, g_zzzzz_xxxyz_0, g_zzzzz_xxxzz_0, g_zzzzz_xxyyy_0, g_zzzzz_xxyyz_0, g_zzzzz_xxyzz_0, g_zzzzz_xxzzz_0, g_zzzzz_xyyyy_0, g_zzzzz_xyyyz_0, g_zzzzz_xyyzz_0, g_zzzzz_xyzzz_0, g_zzzzz_xzzzz_0, g_zzzzz_yyyyy_0, g_zzzzz_yyyyz_0, g_zzzzz_yyyzz_0, g_zzzzz_yyzzz_0, g_zzzzz_yzzzz_0, g_zzzzz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzz_xxxxx_0[i] = 4.0 * g_zzz_xxxxx_0[i] * fbe_0 - 4.0 * g_zzz_xxxxx_1[i] * fz_be_0 + g_zzzz_xxxxx_1[i] * pa_z[i];

        g_zzzzz_xxxxy_0[i] = 4.0 * g_zzz_xxxxy_0[i] * fbe_0 - 4.0 * g_zzz_xxxxy_1[i] * fz_be_0 + g_zzzz_xxxxy_1[i] * pa_z[i];

        g_zzzzz_xxxxz_0[i] = 4.0 * g_zzz_xxxxz_0[i] * fbe_0 - 4.0 * g_zzz_xxxxz_1[i] * fz_be_0 + g_zzzz_xxxx_1[i] * fe_0 + g_zzzz_xxxxz_1[i] * pa_z[i];

        g_zzzzz_xxxyy_0[i] = 4.0 * g_zzz_xxxyy_0[i] * fbe_0 - 4.0 * g_zzz_xxxyy_1[i] * fz_be_0 + g_zzzz_xxxyy_1[i] * pa_z[i];

        g_zzzzz_xxxyz_0[i] = 4.0 * g_zzz_xxxyz_0[i] * fbe_0 - 4.0 * g_zzz_xxxyz_1[i] * fz_be_0 + g_zzzz_xxxy_1[i] * fe_0 + g_zzzz_xxxyz_1[i] * pa_z[i];

        g_zzzzz_xxxzz_0[i] = 4.0 * g_zzz_xxxzz_0[i] * fbe_0 - 4.0 * g_zzz_xxxzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xxxz_1[i] * fe_0 + g_zzzz_xxxzz_1[i] * pa_z[i];

        g_zzzzz_xxyyy_0[i] = 4.0 * g_zzz_xxyyy_0[i] * fbe_0 - 4.0 * g_zzz_xxyyy_1[i] * fz_be_0 + g_zzzz_xxyyy_1[i] * pa_z[i];

        g_zzzzz_xxyyz_0[i] = 4.0 * g_zzz_xxyyz_0[i] * fbe_0 - 4.0 * g_zzz_xxyyz_1[i] * fz_be_0 + g_zzzz_xxyy_1[i] * fe_0 + g_zzzz_xxyyz_1[i] * pa_z[i];

        g_zzzzz_xxyzz_0[i] = 4.0 * g_zzz_xxyzz_0[i] * fbe_0 - 4.0 * g_zzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xxyz_1[i] * fe_0 + g_zzzz_xxyzz_1[i] * pa_z[i];

        g_zzzzz_xxzzz_0[i] = 4.0 * g_zzz_xxzzz_0[i] * fbe_0 - 4.0 * g_zzz_xxzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_xxzz_1[i] * fe_0 + g_zzzz_xxzzz_1[i] * pa_z[i];

        g_zzzzz_xyyyy_0[i] = 4.0 * g_zzz_xyyyy_0[i] * fbe_0 - 4.0 * g_zzz_xyyyy_1[i] * fz_be_0 + g_zzzz_xyyyy_1[i] * pa_z[i];

        g_zzzzz_xyyyz_0[i] = 4.0 * g_zzz_xyyyz_0[i] * fbe_0 - 4.0 * g_zzz_xyyyz_1[i] * fz_be_0 + g_zzzz_xyyy_1[i] * fe_0 + g_zzzz_xyyyz_1[i] * pa_z[i];

        g_zzzzz_xyyzz_0[i] = 4.0 * g_zzz_xyyzz_0[i] * fbe_0 - 4.0 * g_zzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xyyz_1[i] * fe_0 + g_zzzz_xyyzz_1[i] * pa_z[i];

        g_zzzzz_xyzzz_0[i] = 4.0 * g_zzz_xyzzz_0[i] * fbe_0 - 4.0 * g_zzz_xyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_xyzz_1[i] * fe_0 + g_zzzz_xyzzz_1[i] * pa_z[i];

        g_zzzzz_xzzzz_0[i] = 4.0 * g_zzz_xzzzz_0[i] * fbe_0 - 4.0 * g_zzz_xzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_xzzz_1[i] * fe_0 + g_zzzz_xzzzz_1[i] * pa_z[i];

        g_zzzzz_yyyyy_0[i] = 4.0 * g_zzz_yyyyy_0[i] * fbe_0 - 4.0 * g_zzz_yyyyy_1[i] * fz_be_0 + g_zzzz_yyyyy_1[i] * pa_z[i];

        g_zzzzz_yyyyz_0[i] = 4.0 * g_zzz_yyyyz_0[i] * fbe_0 - 4.0 * g_zzz_yyyyz_1[i] * fz_be_0 + g_zzzz_yyyy_1[i] * fe_0 + g_zzzz_yyyyz_1[i] * pa_z[i];

        g_zzzzz_yyyzz_0[i] = 4.0 * g_zzz_yyyzz_0[i] * fbe_0 - 4.0 * g_zzz_yyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_yyyz_1[i] * fe_0 + g_zzzz_yyyzz_1[i] * pa_z[i];

        g_zzzzz_yyzzz_0[i] = 4.0 * g_zzz_yyzzz_0[i] * fbe_0 - 4.0 * g_zzz_yyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_yyzz_1[i] * fe_0 + g_zzzz_yyzzz_1[i] * pa_z[i];

        g_zzzzz_yzzzz_0[i] = 4.0 * g_zzz_yzzzz_0[i] * fbe_0 - 4.0 * g_zzz_yzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_yzzz_1[i] * fe_0 + g_zzzz_yzzzz_1[i] * pa_z[i];

        g_zzzzz_zzzzz_0[i] = 4.0 * g_zzz_zzzzz_0[i] * fbe_0 - 4.0 * g_zzz_zzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_zzzz_1[i] * fe_0 + g_zzzz_zzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

